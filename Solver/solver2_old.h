#ifndef SOLVER2_H
#define SOLVER2_H

#include <limits>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <fstream>

#include "edge_cv.h"
#include "rp_cv.h"
#include "parametrization.h"


#include <unsupported/Eigen/IterativeSolvers>
#include "LinearOperator.h"
#include "opGMRES.h"

#include "BoxTree.h"

#include "mkl.h"
#include <omp.h>
#include "mpi.h"

#include <fftw3.h>
#include "global_interpolator.h"

#include <filesystem> // Required for std::filesystem
#include <algorithm>  // Required for std::count_if

/**
 * Changes by me :)
 * I have modified the Solver class to take an MPI_Comm as an input
 * world_rank_, world_size_ -> comm_rank_, comm_size_
 * Communication is only done within the communicator
 * 
 * I have change many of the global parameters to now be part of the solver class
 * 
 * I have implemented an automatic level picker for IFGF
 * 
 * I have split declaration and intialization.
 * 
 * Added a function to automatically choose the splits per patch. WARNING This function assumes
 * that each patch is approx the same size, so the number of splits per patch can be done globally.
 * In patches are multiscale, then split them manually, write to files, and then pass them in
 * to avoid having to completely change the rectangular polar code. Also it is assumed that
 * the initial patches resolve the geometry well, (less than one oscillation of the geometry)
 */

 /*
 * Changes for once we get the new code
 * Figure out how to make PS and PANG work using templates, and not the global
 * variables. Will need to change BOX tree fundamentally since interpolator is just a global variable.
 *
 * 
 */

 using namespace Eigen;

int inline G_SEARCH_MAX_ITER = 50;
double inline G_SEARCH_TOL = 1E-12;


struct InterpPatch {

    int imax;
    int jmax;
    int kmax;

    long long zoneID;

    std::vector<double> uNodes; // size imax
    std::vector<double> vNodes; // size jmax

    std::vector<double> uWeights; // size imax
    std::vector<double> vWeights; // size jmax

    std::vector<double> x, y, z; // size imax*jmax*kmax
    std::vector<int> mask; // size imax*jmax*kmax
    std::vector<double> dxdu, dydu, dzdu; // size imax*jmax*kmax
    std::vector<double> dxdv, dydv, dzdv; // size imax*jmax*kmax
    std::vector<double> dS; // size imax*jmax*kmax
    std::vector<double> nuX, nuY, nuZ; // size imax*jmax*kmax

}; 

class Solver 
{

    private:

        int Nu_int_, Nv_int_;
        int Nu_prec_, Nv_prec_;

        int Q_, Qx_, Qy_;

        // This is sneaky. The interp_surface is not broadcase to all ranks,
        // so instead it is stored as an unordered map for ease of use
        std::unordered_map<long long, InterpPatch> interp_surface_;

        std::vector<double> fejer_nodes_u_int_, fejer_weights_u_int_;
        std::vector<double> fejer_nodes_v_int_, fejer_weights_v_int_;
        std::vector<double> fejer_nodes_u_prec_, fejer_weights_u_prec_;
        std::vector<double> fejer_nodes_v_prec_, fejer_weights_v_prec_;

        std::vector<std::complex<double>> Tn_;
        std::vector<std::complex<double>> Tm_;

        int comm_rank_;
        int comm_size_;
        MPI_Comm mpi_comm_;

        std::vector<long long> split_points_;
        std::vector<long long> split_points_2_;

        std::vector<int> recv_counts_;
        std::vector<int> displs_;
        std::vector<int> recv_counts_2_;
        std::vector<int> displs_2_;

        long long patch_low_;
        long long patch_up_;
        long long orig_patch_low_;
        long long orig_patch_up_;
        long long point_low_;
        long long point_up_;

        std::vector<int> flags_domain_u_all_;
        std::vector<int> flags_domain_v_all_;

        std::vector<double> disc_points_x_all_;
        std::vector<double> disc_points_y_all_;
        std::vector<double> disc_points_z_all_;

        std::vector<double> norm_points_x_all_;
        std::vector<double> norm_points_y_all_;
        std::vector<double> norm_points_z_all_;
        
        std::vector<double> dsdtjac_all_;

        std::complex<double> coupling_parameter_ = 1.0;
        bool init_compute_coupling_param_ = true;

        std::vector<long long> sing_and_near_sing_patches_estimate_;
        std::vector<long long> start_sing_and_near_sing_patches_estimate_;
        std::vector<long long> size_sing_and_near_sing_patches_estimate_;

        std::vector<std::complex<double>> precomputations_;
        std::vector<long long> point_precomputations_;
        std::vector<long long> patch_num_precomputations_;

        std::vector<long long> patch_num_coeffs_;

        BoxTree boxes_;

        std::vector<long long> new_order_points_IFGF_;
        std::vector<long long> new_order_points_RP_;

        int P_over_ = P_OVERSAMPLING;

        double *original_patch_;
        double *oversampled_patch_;
        fftw_plan plan_original_patch_;
        fftw_plan plan_oversampled_patch_;

        std::vector<double> argmin_precomputations_u_;
        std::vector<double> argmin_precomputations_v_;

        std::complex<double> wavenumber_;
        std::complex<double> lambda_;

        double proximity_;
        int nlevels_;

        bool split_patch_by_wavenumber_ = false;
        double num_wl_per_patch  = 1.0;

    public:
    ///////////////////////// Problem parameters that can be changed directly /////////////
        // Single Layer = 1
        // Double Layer = 2
        // Combined Layer = 3
        int EQUATION_FORMULATION = 3;

        // -1 interior, 1 exterior problem
        double INT_EXT = 1.0;

        int MAX_GMRES_ITER = 50;
        double TOL_GMRES = 1E-6;

        // If true, negate eta so no false poles of the combined operator lie in the lower
        // half plane
        bool use_flipped_eta = false;

        bool USE_ACCELERATOR = true;

        // Incident field parameters:
        // 0 = PLANE WAVE, 1 = POINT SOURCE
        int PLANE_OR_POINT = 0; 
        // Incident direction of the plane wave
        double PLANE_WAVE_THE = 0.0; // in [0, 2pi)
        double PLANE_WAVE_PHI = M_PI; // in [0, pi]
        
        int NUM_POINT_SOURCES = 0;
        // Specify x,y,z, location of point sources in this array
        std::vector<std::array<double,3>> POINT_SOURCE_CENTER;

        double SPHERE_RADIUS = 1.0;
        std::array<double,3> SPHERE_CENTER = {0.0,0.0,0.0};

        // 0 for sphere, 1 for other
        int GEOMETRY = 0;

        // Currently edge geometries are unsupported
        std::vector<bool> EDGE_FLAG_U_A;
        std::vector<bool> EDGE_FLAG_U_B;
        std::vector<bool> EDGE_FLAG_V_A;
        std::vector<bool> EDGE_FLAG_V_B;

        std::string DIRECTORY = "";
        std::string FILE_NAME = "";

    //////////////////////////////////////////////////////////////////////////////////////////


        void HH(const double p1x, const double p1y, const double p1z,
        const double p2x, const double p2y, const double p2z,
        const double nx, const double ny, const double nz,
        const double dxdsx, const double dxdsy, const double dxdsz,
        const double dxdtx, const double dxdty, const double dxdtz,
        const double dxdsdsx, const double dxdsdsy, const double dxdsdsz,
        const double dxdsdtx, const double dxdsdty, const double dxdsdtz,
        const double dxdtdtx, const double dxdtdty, const double dxdtdtz,
        const std::complex<double> coupling_parameter, std::complex<double> wavenumber,
        std::complex<double>& solution); 

        void HH2(const double p1x, const double p1y, const double p1z,
         const double p2x, const double p2y, const double p2z,
         const double nx, const double ny, const double nz,
         std::complex<double> coupling_parameter, std::complex<double> wavenumber,
         std::complex<double>& solution);

        void HH_far(const double xVers_0, const double xVers_1, const double xVers_2,
            const double y_0, const double y_1, const double y_2,
            const double n_0, const double n_1, const double n_2,
            const std::complex<double> coupling_parameter, const double wavenumber,
            std::complex<double>& solution);

        std::vector<std::complex<double>> compute_spherical_harmonics();
        std::vector<std::complex<double>> compute_spherical_harmonics_solution_exact();
        void solve_spherical_harmonics();

        void compute_parallel_parameters(); 

        void load_interpolated_surface(); 
        
        // Function which is called after the surface is loaded to determine how many
        // splits are necessary.
        void compute_number_patch_splits();

        void compute_fejer_nodes_and_weights();

        void compute_chebyshev_evaluations();
        
        void compute_flags_domain(); // Change

        void compute_discretization_domain(); // Change
  

        void compute_coupling_parameter(); // Change

        void compute_near_singular_patches_estimate(); // Change

        void beta(const double r_0, const double r_1, const double r_2,
                  const long long q, const int flag_u_loc, const int flag_v_loc, 
                  const double u_a_loc, const double u_b_loc, const double v_a_loc, const double v_b_loc,
                  const double ubar_loc, const double vbar_loc,
                  std::vector<std::complex<double>>& prec);
        
        void compute_precomputations(); // Change

        void get_coeffs(const long long q, const std::complex<double>* phi,
                        std::complex<double>* coeffs); // Change
       
        void compute_coeffs(const VectorXcd& phi,
                            std::vector<std::complex<double>>& vec_coeffs); // Change

        void int_near(const std::complex<double>* coeffs,
                      const std::complex<double>* precomputations,
                      std::complex<double>& solution);

        void int_far(const double r_0, const double r_1, const double r_2,
                     const long long npatch, 
                     const std::complex<double>* phi,
                     std::complex<double>& solution);
    
        void compute_integral(const VectorXcd& phi, const std::vector<std::complex<double>>& vec_coeffs,
                              VectorXcd& integral); // Eigen
     

        void iterator_function_unacc(const VectorXcd& phi,
                                     VectorXcd& rhs); // Eigen
        
        // Function which chooses the smallest box size possible (number of levels), while
        // still satifying the necessary IFGF assumption
        void setup_IFGF_choose_levels(); // Delete

        void create_IFGF_object(); // Change
        
        void compute_new_order_points_RP(); // Change

        bool check_patch_in_neighbours(); // Change

        void initialize_indexes_HO(); // Change 

        void compute_intensities_patch(const long long npatch,
                                       const std::complex<double>* phi,
                                       std::complex<double>* intensities); // Change

        void compute_intensities(const VectorXcd& phi,
                                 std::vector<std::complex<double>>& intensities); // Eigen
       

        void compute_integral_acc(const VectorXcd& phi, const std::vector<std::complex<double>>& vec_coeffs,
                                  VectorXcd& integral); // Eigen
       
        
        
        inline void fct_4(const double x1, const double x2, const double x3,
                                    const double y1, const double y2, const double y3,
                                    const double normal1, const double normal2, const double normal3,
                                    const std::complex<double> coupling_parameter, const std::complex<double> wavenumber,
                                    const std::complex<double> density, 
                                    std::complex<double>& phi) 
        {

            std::complex<double> kernel;

            HH2(x1, x2, x3, 
                y1, y2, y3,
                normal1, normal2, normal3,
                coupling_parameter, wavenumber,
                kernel);

            const double kerreal = kernel.real();
            const double kerimag = kernel.imag();

            const double phireal = density.real() * kerreal - density.imag() * kerimag;
            const double phiimag = density.real() * kerimag + density.imag() * kerreal;

            phi = {phireal, phiimag};

        }

        inline void fac_1(const double distance, const std::complex<double> wavenumber, std::complex<double>& sol)
        {
            static std::complex<double> I(0.0,1.0);

            sol = std::exp(I * wavenumber * distance) / distance;

        }
        
        
        void iterator_function_acc(const VectorXcd& phi,
                                   VectorXcd& rhs); //Change
       
        void create_fftw_objects(); // Delete
       
        void compute_oversampling_M1(const std::vector<double>& psi, std::vector<double>& psi_overs); // Change
       

        void compute_oversampling_M2(const std::vector<double>& psi, std::vector<double>& psi_overs); // Delete
        
        void compute_singular_points(); // Change
        
        void int_near_overs(const double r_0, const double r_1, const double r_2,
                            const long long q, const int flag_u_loc, const int flag_v_loc, 
                            const double u_a_loc, const double u_b_loc, const double v_a_loc, const double v_b_loc,
                            const double ubar_loc, const double vbar_loc,
                            const double* psi_overs,
                            std::complex<double>& solution);  // Change
        

        void write_psi_M1(const VectorXcd& phi, std::vector<double>& psi); // Change
        

        void write_psi_M2(const long long patch_num, const std::complex<double>* phi, std::vector<double>& psi); // delete
       

        void compute_integral_overs_M1(const VectorXcd& phi,
                                       const std::vector<double>& psi_overs,
                                       VectorXcd& integral); // change
       
        void compute_integral_overs_M2(const VectorXcd& phi,
                                       VectorXcd& integral); // delete
        
        void compute_integral_overs_acc_M1(const VectorXcd& phi,
                                           const std::vector<double>& psi_overs,
                                           VectorXcd& integral); // change
        
        void compute_integral_overs_acc_M2(const VectorXcd& phi,
                                           VectorXcd& integral); // delete
       
        void iterator_function_overs(const VectorXcd& phi,
                                     VectorXcd& rhs); // change

        void iterator_function_overs_acc(const VectorXcd& phi,
                                         VectorXcd& rhs); // change

        void iterator_function(const VectorXcd& phi,
                               VectorXcd& rhs); // change
        
        void setup(bool timing); // merge
       
        VectorXcd solve(const VectorXcd& rhs); // change
        

        VectorXcd compute_incident_field(); // change

        VectorXcd solve_u_inc(bool timing); // change

        void int_far_field(const double xVers_0, const double xVers_1, const double xVers_2,
                           const long long npatch, 
                           const std::complex<double>* phi,
                           std::complex<double>& solution); // Change  

        std::complex<double> compute_far_field_approx(const VectorXcd& phi, 
                                                      const double xVers_0, const double xVers_1, const double xVers_2); // Change
       

        std::complex<double> compute_far_field_exact(const double xVers_0, const double xVers_1, const double xVers_2); 
        
        void compute_far_field_error(const bool timing, const VectorXcd& phi); // change
        
        //std::complex<double> compute_incident_field(double x, double y, double z) 
        

        void int_near_field(const double x_0, const double x_1, const double x_2,
                            const long long npatch, 
                            const std::complex<double>* phi,
                            std::complex<double>& solution); // change                     
       

        std::complex<double> compute_near_field_approx(const VectorXcd& phi, 
                                                       const double x_0, const double x_1, const double x_2); // change
        
        /*
        Custom function to evaluate the near field, based off of a small change from compute
        near_field_approx
        */
        std::complex<double> compute_near_field(const VectorXcd& phi, 
                                                       const double x_0, const double x_1, const double x_2); // merge
       
        // void compute_near_field_error(const bool timing, const VectorXcd& phi)
       
        std::complex<double> iterator_function_Manuel(const VectorXcd& u, const VectorXcd& v); // change / delete 
       

        ////////////////////// Getters/ Setters //////////////////////////////////////////
        MPI_Comm get_mpi_comm() const { return mpi_comm_;}

        // Get the number of patches origonally
        long long get_num_unknowns() const {return Q_ * Qx_ * Qy_ * Nv_int_ * Nu_int_;}

        int get_patch_split_x() const {return Qx_;}
        int get_patch_split_y() const {return Qy_;}
        int get_nlevels_IFGF()  const {return nlevels_;}
        std::complex<double> get_coup_param() const {return coupling_parameter_;}

        std::vector<double> get_disc_points_x() const {return  disc_points_x_all_;}
        std::vector<double> get_disc_points_y() const {return  disc_points_y_all_;}
        std::vector<double> get_disc_points_z() const {return  disc_points_z_all_;}

        void set_eq_form(int form) {EQUATION_FORMULATION = form;}
        void set_int_ext(int int_ext) {INT_EXT = int_ext;}
        void set_coup_param(std::complex<double> cp) {coupling_parameter_ = cp; init_compute_coupling_param_ = false;}
        void set_num_wl_per_patch(double num_wl) {num_wl_per_patch = num_wl;}
        //////////////////////////////////////////////////////////////////////////////////

        // basic code to initialize the solver, after calling the constructor, which
        // just constructs the geometry.
        // In between you can set problem specific paramters
        void init_solver(const bool timing, 
        const std::complex<double> k, 
        const int* n_pts_per_patch,
        const int* n_split_per_patch,
        const int* n_pts_sing_int,
        const double proximity_box_size,
        const int n_levels_IFGF, 
        const MPI_Comm& mpi_comm);



        // Constructor for the sphere
        Solver(double sphere_radius = 1.0,
               double sphere_centerX = 0.0, double sphere_centerY = 0.0, double sphere_centerZ = 0.0);

        // Constructor for other geometries. This assumes that all files for a certain
        // geometry live in one directory, and that they are the only files in that directory
        // Files should be named "file_prefix1,file_prefix2, etc";
        // Pass in n_split_per_patch = {-1,-1} for automatic patch splitting based on the 
        // wave number
        Solver( const std::string directory,
                const std::string file_prefix);
 

        ~Solver();

};


#endif