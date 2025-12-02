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

#include "BoxTree.h"

#include "mkl.h"
#include <omp.h>
#include "mpi.h"

#include <iomanip>
#include <cmath>

#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>



inline void HH(
    const double p1x, const double p1y, const double p1z,
    const double p2x, const double p2y, const double p2z,
    const double nx, const double ny, const double nz,
    const double dxdsx, const double dxdsy, const double dxdsz,
    const double dxdtx, const double dxdty, const double dxdtz,
    const double dxdsdsx, const double dxdsdsy, const double dxdsdsz,
    const double dxdsdtx, const double dxdsdty, const double dxdsdtz,
    const double dxdtdtx, const double dxdtdty, const double dxdtdtz,
    const double coupling_parameter,
    std::complex<double>& solution)
{
    
    const double INV_4_PI = 1.0 / (4.0 * M_PI);

    const double dx = p2x - p1x;
    const double dy = p2y - p1y;
    const double dz = p2z - p1z;
    
    const double norm_diff_sq = dx * dx + dy * dy + dz * dz;

    if (norm_diff_sq < 1e-20) { 

        solution = std::complex<double>(0.0, 0.0);
        return;

    } else {
        
        const double norm_diff = std::sqrt(norm_diff_sq);
        const double inv_norm_diff = 1.0 / norm_diff;
        const double inv_norm_diff_sq = inv_norm_diff * inv_norm_diff;

        const double val_cos = std::cos(WAVE_NUMBER * norm_diff);
        const double val_sin = std::sin(WAVE_NUMBER * norm_diff);

        double solution_real, solution_imag;
        
        switch (EQUATION_FORMULATION) {
            
            case 1: {

                solution_real = val_cos * INV_4_PI * inv_norm_diff;
                solution_imag = val_sin * INV_4_PI * inv_norm_diff;
                break;

            }

            case 2: {
                
                // Compute beta = < r, n > * |r|^-2

                const double rDotNorm = dx * nx + dy * ny + dz * nz;
                double beta = rDotNorm * inv_norm_diff_sq;

                if ((norm_diff_sq < 1.0E-10) && (std::abs(rDotNorm) < 1.0E-6)) {

                    beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, 
                                nx, ny, nz, dxdsx, dxdsy, dxdsz, 
                                dxdtx, dxdty, dxdtz, dxdsdsx, dxdsdsy, dxdsdsz, 
                                dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

                }
                
                const double term1 = val_cos * inv_norm_diff;
                const double term2 = WAVE_NUMBER * val_sin;
                
                solution_real = -INV_4_PI * beta * (term1 + term2);

                const double term3 = val_sin * inv_norm_diff;
                const double term4 = WAVE_NUMBER * val_cos;
                
                solution_imag = -INV_4_PI * beta * (term3 - term4);
                break;

            }
            
            case 3: {

                const double solution_real_S = val_cos * INV_4_PI * inv_norm_diff;
                const double solution_imag_S = val_sin * INV_4_PI * inv_norm_diff;

                const double rDotNorm = dx * nx + dy * ny + dz * nz;
                double beta = rDotNorm * inv_norm_diff_sq;

                if ((norm_diff_sq < 1.0E-10) && (std::abs(rDotNorm) < 1.0E-6)) {

                    beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, 
                                nx, ny, nz, dxdsx, dxdsy, dxdsz, 
                                dxdtx, dxdty, dxdtz, dxdsdsx, dxdsdsy, dxdsdsz, 
                                dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

                }

                const double term1 = val_cos * inv_norm_diff;
                const double term2 = WAVE_NUMBER * val_sin;
                const double solution_real_D = -INV_4_PI * beta * (term1 + term2);

                const double term3 = val_sin * inv_norm_diff;
                const double term4 = WAVE_NUMBER * val_cos;
                const double solution_imag_D = -INV_4_PI * beta * (term3 - term4);

                solution_real = solution_real_D + coupling_parameter * solution_imag_S;
                solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

                break;

            }

            default: {
                
                solution_real = 0.0;
                solution_imag = 0.0;
                break;

            }

        }

        solution = std::complex<double>{solution_real, solution_imag};

    }

}

inline void HH2(const double p1x, const double p1y, const double p1z,
         const double p2x, const double p2y, const double p2z,
         const double nx, const double ny, const double nz,
         double coupling_parameter,
         std::complex<double>& solution)
{

    double solution_real, solution_imag;

    const double norm_diff = std::sqrt((p2x - p1x)*(p2x - p1x) + (p2y - p1y)*(p2y - p1y) + (p2z - p1z)*(p2z - p1z));

    if (norm_diff < 1e-14) {

        solution_real = 0.0;
        solution_imag = 0.0;

    } else {   

        const double val_cos = std::cos(WAVE_NUMBER * norm_diff);
        const double val_sin = std::sin(WAVE_NUMBER * norm_diff);

        if (EQUATION_FORMULATION == 1) {

            solution_real = val_cos / (4.0 * M_PI * norm_diff);
            solution_imag = val_sin / (4.0 * M_PI * norm_diff);

        } else if (EQUATION_FORMULATION == 2) {

            // Compute beta = < r, n > / |r|^2
            
            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            solution_real = -(1.0 / (4.0 * M_PI)) * beta * (val_cos / norm_diff + WAVE_NUMBER * val_sin);
            solution_imag = -(1.0 / (4.0 * M_PI)) * beta * (val_sin / norm_diff - WAVE_NUMBER * val_cos);        

            if (USE_ACCELERATOR) {

                solution_real *= -1.0;
                solution_imag *= -1.0;

            }

        } else {

            // EQUATION_FORMULATION == 3 

            const double solution_real_S = val_cos / (4.0 * M_PI * norm_diff);
            const double solution_imag_S = val_sin / (4.0 * M_PI * norm_diff);

            // Compute beta = < r, n > / |r|^2

            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            double solution_real_D = -(1.0 / (4.0 * M_PI)) * beta * (val_cos / norm_diff + WAVE_NUMBER * val_sin);
            double solution_imag_D = -(1.0 / (4.0 * M_PI)) * beta * (val_sin / norm_diff - WAVE_NUMBER * val_cos);        

            if (USE_ACCELERATOR) {

                solution_real_D *= -1.0;
                solution_imag_D *= -1.0;

            }

            solution_real = solution_real_D + coupling_parameter * solution_imag_S;
            solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

        }

        solution = std::complex<double>{solution_real, solution_imag};

    }

}

inline void HH_far(const double xVers_0, const double xVers_1, const double xVers_2,
            const double y_0, const double y_1, const double y_2,
            const double n_0, const double n_1, const double n_2,
            const double coupling_parameter,
            std::complex<double>& solution)
{

    double solution_real, solution_imag;

    const double inner_prod = - WAVE_NUMBER * (xVers_0 * y_0 + xVers_1 * y_1 + xVers_2 * y_2);
    const double inner_prod_2 = - WAVE_NUMBER * (xVers_0 * n_0 + xVers_1 * n_1 + xVers_2 * n_2);

    const double S_real = std::cos(inner_prod);
    const double S_imag = std::sin(inner_prod);
    const double D_real = -inner_prod_2 * S_imag;
    const double D_imag = inner_prod_2 * S_real;

    if (EQUATION_FORMULATION == 1) {

        solution_real = S_real;
        solution_imag = S_imag;

    } else if (EQUATION_FORMULATION == 2) {

        solution_real = D_real;
        solution_imag = D_imag;

    } else {

        // EQUATION_FORMULATION == 3

        solution_real = D_real + coupling_parameter * S_imag;
        solution_imag = D_imag - coupling_parameter * S_real;

    } 

    solution = std::complex<double>{solution_real, solution_imag};

}

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

        static constexpr long long Nu_int_ = N_PTS_PER_PATCH[0], Nv_int_ = N_PTS_PER_PATCH[1];
        static constexpr long long Nu_prec_ = N_PTS_SING_INT[0], Nv_prec_ = N_PTS_SING_INT[1];

        static constexpr long long Q_ = N_PATCHES_ORIG, Qx_ = N_SPLIT_PER_PATCH[0], Qy_ = N_SPLIT_PER_PATCH[1];

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

        inline static void fct_4(const double x1, const double x2, const double x3,
                                 const double y1, const double y2, const double y3,
                                 const double normal1, const double normal2, const double normal3,
                                 const double coupling_parameter,
                                 const std::complex<double> density, 
                                 std::complex<double>& phi) 
        {

            std::complex<double> kernel;

            HH2(x1, x2, x3, 
                y1, y2, y3,
                normal1, normal2, normal3,
                coupling_parameter,
                kernel);

            const double kerreal = kernel.real();
            const double kerimag = kernel.imag();

            const double phireal = density.real() * kerreal - density.imag() * kerimag;
            const double phiimag = density.real() * kerimag + density.imag() * kerreal;

            phi = {phireal, phiimag};

        }

        inline static void fac_1(const double distance, std::complex<double>& sol)
        {

            const double re = std::cos(WAVE_NUMBER * distance) / distance;
            const double im = std::sin(WAVE_NUMBER * distance) / distance;

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
       

        Solver(bool timing)
        {

            setup(timing);

        }

        ~Solver() {           

        }    

};

#endif