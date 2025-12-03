#ifndef SOLVER2_H
#define SOLVER2_H

#include <limits>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <fstream>
#include <algorithm>

#include "edge_cv.h"
#include "rp_cv.h"
#include "parametrization.h"
#include "GreenFunctions.hpp"

#include "BoxTree.h"

#include "mkl.h"
#include <omp.h>
#include "mpi.h"

#include <iomanip>
#include <cmath>

#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>





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

        // TODO: MOVE THESE to the initialization
        long long Nu_int_, Nv_int_, Nu_prec_, Nv_prec_;

        long long Q_, Qx_, Qy_;

        std::unordered_map<long long, InterpPatch> interp_surface_;

        std::vector<double> fejer_nodes_u_int_, fejer_weights_u_int_;
        std::vector<double> fejer_nodes_v_int_, fejer_weights_v_int_;
        std::vector<double> fejer_nodes_u_prec_, fejer_weights_u_prec_;
        std::vector<double> fejer_nodes_v_prec_, fejer_weights_v_prec_;
        std::vector<double> fejer_weights_u_v_int_;

        std::vector<std::complex<double>> Tn_;
        std::vector<std::complex<double>> Tm_;

        int world_rank_;
        int world_size_;

        MPI_Comm mpi_comm_;

        std::vector<long long> split_points_;
        std::vector<long long> split_points_2_;

        std::vector<MPI_Count> recv_counts_;
        std::vector<MPI_Aint> displs_;
        std::vector<MPI_Count> recv_counts_2_;
        std::vector<MPI_Aint> displs_2_;

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

        double coupling_parameter_;

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

        std::unordered_map<long long, int> flags_domain_u_not_in_rank_;
        std::unordered_map<long long, int> flags_domain_v_not_in_rank_;

        std::unordered_map<long long, std::vector<double>> disc_points_x_not_in_rank_;
        std::unordered_map<long long, std::vector<double>> disc_points_y_not_in_rank_;
        std::unordered_map<long long, std::vector<double>> disc_points_z_not_in_rank_;
        
        std::unordered_map<long long, std::vector<double>> dsdtjac_not_in_rank_;

    public:
        //////////////////////////////////// Parameters //////////////////////////////////
        int G_SEARCH_MAX_ITER = 50;
        double G_SEARCH_TOL = 1E-12;


        // Choose integral equation formulation:

        // Single Layer = 1
        // Double Layer = 2
        // Combined Layer = 3 
        // 4
        int EQUATION_FORMULATION = 3;

        // Discretization parameters:

        // Number of points used per dimension per patch in every patch
        // TODO: MOVED CONSTEXPR FROM HERE
        // long long N_PTS_PER_PATCH[2] = {6, 6};
        // // Number of patch partitions along first and second dimensions
        // long long N_SPLIT_PER_PATCH[2] = {16, 16};
        
        // // Singular integration parameters:

        // // Number of points used with the singular integration
        // int N_PTS_SING_INT[2] = {40, 40};
        // Proximity distance that determines near singular integrals
        int DELTA_METHOD = 1;
        double PROXIMITY_BOX_SIZE = 0.1/16;
        double PERCENT_BOX_SIZE = 0.15;

        // GMRES options:
        int MAX_ITER = 100;
        double TOL_GMRES = 1E-4;

        
        // Geometry parameters:

        // Sphere = 0
        // Others = 1

        int GEOMETRY = 0;

        // Sphere
        double SPHERE_RADIUS = 1.0;
        std::array<double,3> SPHERE_CENTER = {0.0, 0.0, 0.0};

        // Currently edge geometries are unsupported
        std::vector<bool> EDGE_FLAG_U_A;
        std::vector<bool> EDGE_FLAG_U_B;
        std::vector<bool> EDGE_FLAG_V_A;
        std::vector<bool> EDGE_FLAG_V_B;

        std::string DIRECTORY = "";
        std::string FILE_NAME = "";

        // Incident field parameters:
        // 0 = PLANE WAVE, 1 = POINT SOURCE
        int PLANE_OR_POINT = 0; 

        // Incident plane wave parameters (PLANE_OR_POINT = 0):
        // kx = k*cos(the)*sin(phi)
        // ky = k*sin(the)*sin(phi)
        // kz = k*cos(phi)

        // TODO: Move lambda out
        double LAMBDA = 2.0 * SPHERE_RADIUS / 8.0; // 2.0 * M_PI / WAVE_NUMBER
        double WAVE_NUMBER = 2.0 * M_PI / LAMBDA;
        double PLANE_WAVE_THE = 0.0; // in [0, 2pi)
        double PLANE_WAVE_PHI = M_PI; // in [0, pi]

        // Incident source points (PLANE_OR_POINT = 1):
        int NUM_POINT_SOURCES = 2;
        std::vector<std::vector<double>> POINT_SOURCE_CENTER = {{0.0, 0.0, 3.3},
                                                                    {0.0, 3.3, 0.0}};


        ///////////////////////////////// Methods ////////////////////////////////////////

        void compute_parallel_parameters(); // ConstructorsandSetup.cpp
      
        void load_interpolated_surface(); // ConstructorsandSetup.cpp

        void compute_fejer_nodes_and_weights(); // ConstructorsandSetup.cpp

        void compute_chebyshev_evaluations(); // ConstructorsandSetup.cpp

        void compute_flags_domain(); // ConstructorsandSetup.cpp
        
        void compute_discretization_domain(); // ConstructorsandSetup.cpp   

        void compute_coupling_parameter(); // ConstructorsandSetup.cpp

        void compute_near_singular_patches_estimate(); // ConstructorsandSetup.cpp

        void beta(const double r_0, const double r_1, const double r_2,
                  const long long q, const int flag_u_loc, const int flag_v_loc, 
                  const double u_a_loc, const double u_b_loc, const double v_a_loc, const double v_b_loc,
                  const double ubar_loc, const double vbar_loc,
                  std::vector<std::complex<double>>& prec); // RP.cpp
        
        void compute_precomputations(); // RP.cpp
        
        void get_coeffs(const long long q, const std::complex<double>* phi,
                        std::complex<double>* coeffs); // RP.cpp
       
        void compute_coeffs(std::complex<double>* phi,
                            std::unordered_map<long long, std::vector<std::complex<double>>>& phi_not_in_rank,
                            std::vector<std::complex<double>>& vec_coeffs); // RP.cpp
       

        void int_near(const std::complex<double>* coeffs,
                      const std::complex<double>* precomputations,
                      std::complex<double>& solution); //RP.cpp

        void int_far(const double r_0, const double r_1, const double r_2,
            const long long npatch, 
            const std::complex<double>* phi,
            std::complex<double>& solution);

        void create_IFGF_object(); // IFGF.cpp
           
        void set_precomputations_data_IFGF(); // IFGF.cpp
        
        void compute_new_order_points_RP(); // IFGF.cpp 

        void check_patch_in_neighbours(); // IFGF.cpp
        

        inline void static fct_4(const double x1, const double x2, const double x3,
                                 const double y1, const double y2, const double y3,
                                 const double normal1, const double normal2, const double normal3,
                                 const double coupling_parameter, const double wavenumber, const int equation_formulation,
                                 const std::complex<double> density, 
                                 std::complex<double>& phi) 
        {

            std::complex<double> kernel;

            HH2(x1, x2, x3, 
                y1, y2, y3,
                normal1, normal2, normal3,
                coupling_parameter, wavenumber, equation_formulation,
                kernel);

            const double kerreal = kernel.real();
            const double kerimag = kernel.imag();

            const double phireal = density.real() * kerreal - density.imag() * kerimag;
            const double phiimag = density.real() * kerimag + density.imag() * kerreal;

            phi = {phireal, phiimag};

        }

        // TODO: 
        inline void static fac_1(const double distance, double wavenumber, std::complex<double>& sol)
        {

            const double re = std::cos(wavenumber * distance) / distance;
            const double im = std::sin(wavenumber * distance) / distance;

            sol = {re, im};

        }

        void compute_intensities_patch(const long long npatch,
                                       const std::complex<double>* phi, 
                                       std::complex<double>* intensities); // RP.cpp
       

        void compute_intensities(const std::complex<double>* phi,
                                 std::vector<std::complex<double>>& intensities); // RP.cpp
        
        void compute_integral_acc(const std::complex<double>* phi,
                                  const std::vector<std::complex<double>>& vec_coeffs,
                                  std::complex<double>* integral); // SolveandEval.cpp
      

        void redistribute_data_RP(const std::vector<std::complex<double>>& intensities, std::complex<double>* rhs); // RP.cpp
        

        void iterator_function_acc(std::complex<double>* phi, std::complex<double>* rhs); // SolveandEval.cpp
        

        void compute_phi_not_in_rank(const std::complex<double>* phi, std::unordered_map<long long, 
            std::vector<std::complex<double>>>& phi_not_in_rank); // SolveandEval.cpp
       

        void iterator_function(std::complex<double>* phi, std::complex<double>* rhs); // SolveandEval.cpp


        void setup(bool timing); // ConstructorsandSetup.cpp
        
        static PetscErrorCode iterator_function_2(Mat A, Vec x, Vec y); // SolveandEval.cpp
       

        std::vector<std::complex<double>> solve(const std::vector<std::complex<double>>& rhs); // SolveandEval.cpp 
        

        std::vector<std::complex<double>> compute_incident_field(); // SolveandEval.cpp 
        

        std::complex<double> compute_incident_field(double x, double y, double z);  // SolveandEval.cpp 
               

        std::vector<std::complex<double>> solve_u_inc(bool timing); 
       
        void int_far_field(const double xVers_0, const double xVers_1, const double xVers_2,
                           const long long npatch, 
                           const std::complex<double>* phi,
                           std::complex<double>& solution); // SolveandEval.cpp 

        std::complex<double> compute_far_field_approx(const std::vector<std::complex<double>>& phi, 
                                                      const double xVers_0, const double xVers_1, const double xVers_2); // SolveandEval.cpp    
        

        void int_near_field(const double x_0, const double x_1, const double x_2,
                            const long long npatch, 
                            const std::complex<double>* phi,
                            std::complex<double>& solution); // SolveandEval.cpp                      


        std::complex<double> compute_near_field_approx(const std::vector<std::complex<double>>& phi, 
                                                       const double x_0, const double x_1, const double x_2); // SolveandEval.cpp  


        std::complex<double> compute_far_field_exact(const double xVers_0, const double xVers_1, const double xVers_2); // SolveandEval.cpp   
       

        void compute_far_field_error(const bool timing, const std::vector<std::complex<double>>& phi); // SolveandEval.cpp 
      

        void compute_near_field(const bool timing, const std::vector<std::complex<double>>& phi); // SolveandEval.cpp 
       

        void init_solver(const bool timing, 
        const std::complex<double> k, 
        const int* n_pts_per_patch,
        const int* n_split_per_patch,
        const int* n_pts_sing_int,
        const double proximity_box_size,
        const int n_levels_IFGF, 
        const MPI_Comm& mpi_comm);

        void init_solver(const bool timing, 
            const double k);



        // Constructor for the sphere
        Solver(double sphere_radius = 1.0,
               double sphere_centerX = 0.0, double sphere_centerY = 0.0, double sphere_centerZ = 0.0);

        // Constructor for other geometries. This assumes that all files for a certain
        // geometry live in one directory, and that they are the only files in that directory
        // Files should be named "file_prefix1,file_prefix2, etc";
        Solver( const std::string directory,
                const std::string file_prefix);


        // Solver(bool timing)
        // {

        //     setup(timing);

        // }

        ~Solver() {           

        }    

};

#endif