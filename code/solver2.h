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
                      std::complex<double>& solution) //RP.cpp
        {         

            solution = std::complex<double>(0.0, 0.0);

            #pragma omp simd
            for (int i = 0; i < Nu_int_*Nv_int_; i++) {

                solution += precomputations[i] * coeffs[i];

            }

        }

        void int_far(const double r_0, const double r_1, const double r_2,
            const long long npatch, 
            const std::complex<double>* phi,
            std::complex<double>& solution)
            
        {

        solution = std::complex<double>(0.0, 0.0);

        const long long size = Nu_int_*Nv_int_;

        const double* px_ptr = &disc_points_x_all_[npatch * size];
        const double* py_ptr = &disc_points_y_all_[npatch * size];
        const double* pz_ptr = &disc_points_z_all_[npatch * size];

        const double* nx_ptr = &norm_points_x_all_[npatch * size];
        const double* ny_ptr = &norm_points_y_all_[npatch * size];
        const double* nz_ptr = &norm_points_z_all_[npatch * size];

        const double* dsdtjac_ptr = &dsdtjac_all_[npatch * size];

        #pragma omp simd
        for (int i = 0; i < size; i++) {

            const double dsdtjac_loc = dsdtjac_ptr[i];
            const double constant = dsdtjac_loc * fejer_weights_u_v_int_[i];

            const double px = px_ptr[i];
            const double py = py_ptr[i];
            const double pz = pz_ptr[i];

            const double nx = nx_ptr[i];
            const double ny = ny_ptr[i];
            const double nz = nz_ptr[i];

            std::complex<double> kernel;
            HH2(r_0, r_1, r_2, px, py, pz, nx, ny, nz, coupling_parameter_, kernel);
            
            solution += constant * kernel * phi[i];

        }

        }

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
        

        std::vector<std::complex<double>> compute_incident_field() 
        {

            std::vector<std::complex<double>> u_inc(point_up_ - point_low_, 0);

            if (PLANE_OR_POINT == 0) {

                const double k_hat_0 = WAVE_NUMBER * std::cos(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
                const double k_hat_1 = WAVE_NUMBER * std::sin(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
                const double k_hat_2 = WAVE_NUMBER * std::cos(PLANE_WAVE_PHI);
                
                #pragma omp parallel for
                for (long long i = 0; i < point_up_ - point_low_; i++) {

                    long long idx = USE_ACCELERATOR ? i : i + point_low_;

                    const double x_0 = disc_points_x_all_[idx];
                    const double x_1 = disc_points_y_all_[idx];
                    const double x_2 = disc_points_z_all_[idx];               

                    const double inner_product = x_0 * k_hat_0 + x_1 * k_hat_1 + x_2 * k_hat_2;

                    const std::complex<double> value((-1.0) * std::cos(inner_product), (-1.0) * std::sin(inner_product));

                    u_inc[i] = value;

                }

            } else {

                const double k = WAVE_NUMBER;
                
                #pragma omp parallel for
                for (long long j = 0; j < point_up_ - point_low_; j++) {

                    long long idx = USE_ACCELERATOR ? j : j + point_low_;

                    const double y_0 = disc_points_x_all_[idx];
                    const double y_1 = disc_points_y_all_[idx];
                    const double y_2 = disc_points_z_all_[idx];  

                    for (int i = 0; i < NUM_POINT_SOURCES; i++) {

                        const double x_0 = POINT_SOURCE_CENTER[i][0];
                        const double x_1 = POINT_SOURCE_CENTER[i][1];
                        const double x_2 = POINT_SOURCE_CENTER[i][2];                                 

                        const double diff = std::sqrt((x_0-y_0)*(x_0-y_0) + (x_1-y_1)*(x_1-y_1) + (x_2-y_2)*(x_2-y_2));

                        const std::complex<double> value((-1.0) * std::cos(k * diff) / diff, (-1.0) * std::sin(k * diff) / diff);

                        u_inc[j] += value;

                    }

                }

            }

            return u_inc;

        }          

        std::complex<double> compute_incident_field(double x, double y, double z) 
        {

            std::complex<double> u_inc = {0.0, 0.0};

            if (PLANE_OR_POINT == 0) {

                const double k_hat_0 = WAVE_NUMBER * std::cos(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
                const double k_hat_1 = WAVE_NUMBER * std::sin(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
                const double k_hat_2 = WAVE_NUMBER * std::cos(PLANE_WAVE_PHI);            

                const double inner_product = x * k_hat_0 + y * k_hat_1 + z * k_hat_2;

                const std::complex<double> value((-1.0) * std::cos(inner_product), (-1.0) * std::sin(inner_product));

                u_inc = value;

            } else {

                const double k = WAVE_NUMBER;

                for (int i = 0; i < NUM_POINT_SOURCES; i++) {

                    const double x_0 = POINT_SOURCE_CENTER[i][0];
                    const double x_1 = POINT_SOURCE_CENTER[i][1];
                    const double x_2 = POINT_SOURCE_CENTER[i][2];                                 

                    const double diff = std::sqrt((x_0-x)*(x_0-x) + (x_1-y)*(x_1-y) + (x_2-z)*(x_2-z));

                    const std::complex<double> value((-1.0) * std::cos(k * diff) / diff, (-1.0) * std::sin(k * diff) / diff);

                    u_inc += value;

                }

            }

            return u_inc;

        }              

        std::vector<std::complex<double>> solve_u_inc(bool timing) 
        {
            
            MPI_Barrier(MPI_COMM_WORLD);

            double start_1 = MPI_Wtime();

            std::vector<std::complex<double>> rhs = compute_incident_field();

            double end_1 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time compute incident field: " << end_1 - start_1 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            double start_2 = MPI_Wtime();

            std::vector<std::complex<double>> solution = solve(rhs);

            //std::complex<double>* rhs_data = rhs.data();
            //std::complex<double> solution2[point_up_-point_low_];
            //iterator_function(rhs_data, solution2);

            double end_2 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time solve: " << end_2 - start_2 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            return solution;

        }

        void int_far_field(const double xVers_0, const double xVers_1, const double xVers_2,
                           const long long npatch, 
                           const std::complex<double>* phi,
                           std::complex<double>& solution)                      
        {

            solution = std::complex<double>(0.0, 0.0);

            const long long size = Nu_int_*Nv_int_;

            const double* px_ptr = &disc_points_x_all_[npatch * size];
            const double* py_ptr = &disc_points_y_all_[npatch * size];
            const double* pz_ptr = &disc_points_z_all_[npatch * size];

            const double* nx_ptr = &norm_points_x_all_[npatch * size];
            const double* ny_ptr = &norm_points_y_all_[npatch * size];
            const double* nz_ptr = &norm_points_z_all_[npatch * size];

            const double* dsdtjac_ptr = &dsdtjac_all_[npatch * size];

            #pragma omp simd
            for (int i = 0; i < size; i++) {

                const double dsdtjac_loc = dsdtjac_ptr[i];
                const double constant = dsdtjac_loc * fejer_weights_u_v_int_[i];

                const double px = px_ptr[i];
                const double py = py_ptr[i];
                const double pz = pz_ptr[i];

                const double nx = nx_ptr[i];
                const double ny = ny_ptr[i];
                const double nz = nz_ptr[i];

                std::complex<double> kernel;
                HH_far(xVers_0, xVers_1, xVers_2, px, py, pz, nx, ny, nz, coupling_parameter_, kernel);
                
                solution += constant * kernel * phi[i];

            }

        }

        std::complex<double> compute_far_field_approx(const std::vector<std::complex<double>>& phi, 
                                                      const double xVers_0, const double xVers_1, const double xVers_2) 
        {
            
            std::complex<double> solution_rank(0.0, 0.0);

            for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

                std::complex<double> solution_loc;

                const long long patch = USE_ACCELERATOR ? (patch_num - patch_low_) : patch_num;

                int_far_field(xVers_0, xVers_1, xVers_2,
                              patch,
                              &phi[(patch_num - patch_low_) * Nu_int_*Nv_int_],
                              solution_loc);                        
                
                solution_rank += solution_loc;

            }

            solution_rank *= 1.0 / (4.0 * M_PI);

            return solution_rank;

        }

        void int_near_field(const double x_0, const double x_1, const double x_2,
                            const long long npatch, 
                            const std::complex<double>* phi,
                            std::complex<double>& solution)                      
        {

            int_far(x_0, x_1, x_2,
                    npatch,
                    phi,
                    solution);

        }

        std::complex<double> compute_near_field_approx(const std::vector<std::complex<double>>& phi, 
                                                       const double x_0, const double x_1, const double x_2) 
        {
            
            std::complex<double> solution_rank(0.0, 0.0);

            for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

                std::complex<double> solution_loc;

                const long long patch = USE_ACCELERATOR ? (patch_num - patch_low_) : patch_num;

                int_near_field(x_0, x_1, x_2,
                               patch,
                               &phi[patch_num * Nu_int_*Nv_int_],
                               solution_loc);                        
                
                solution_rank += solution_loc;

            }

            return solution_rank;

        }

        std::complex<double> compute_far_field_exact(const double xVers_0, const double xVers_1, const double xVers_2) 
        {

            std::complex<double> solution(0.0, 0.0);

            const std::complex<double> I(0.0, 1.0);

            long long nterms = 0;

            std::complex<double> yn = std::sph_neumann(nterms, WAVE_NUMBER * SPHERE_RADIUS);
            std::complex<double> jn = std::sph_bessel(nterms, WAVE_NUMBER * SPHERE_RADIUS);
            std::complex<double> hn = jn + I * yn;

            while (std::abs(hn) * std::abs(hn) < 1.0e14) {

                nterms++;

                yn = std::sph_neumann(nterms, WAVE_NUMBER * SPHERE_RADIUS);
                jn = std::sph_bessel(nterms, WAVE_NUMBER * SPHERE_RADIUS);
                hn = jn + I * yn;

            }

            const double kVers[3] = {std::cos(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI), std::sin(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI), std::cos(PLANE_WAVE_PHI)};
            const double x = xVers_0 * kVers[0] + xVers_1 * kVers[1] + xVers_2 * kVers[2];
            
            for (long long n = 0; n <= nterms; n++) {

                const std::complex<double> yn = std::sph_neumann(n, WAVE_NUMBER * SPHERE_RADIUS);
                const std::complex<double> jn = std::sph_bessel(n, WAVE_NUMBER * SPHERE_RADIUS);
                const std::complex<double> hn = jn + I * yn;

                const double P = std::legendre(n, x);

                solution += (2.0 * n + 1.0) * jn/hn * P;

            }

            solution *= I / (WAVE_NUMBER * SPHERE_RADIUS);

            return solution;

        }

        void compute_far_field_error(const bool timing, const std::vector<std::complex<double>>& phi)
        {           

            MPI_Barrier(MPI_COMM_WORLD);
            double start = MPI_Wtime();

            const double deltaPhi = (M_PI - 2.0 * 1.0e-5) / (N_FAR_PTS[0] - 1);
            const double deltaTheta = 2.0 * M_PI / (N_FAR_PTS[1] - 1);

            const long long total_pts = N_FAR_PTS[0] * N_FAR_PTS[1];

            std::vector<long long> split_points_far_field(world_size_ + 1);
            const long long points_per_rank = total_pts / world_size_;
            long long remaining_points = total_pts % world_size_;

            split_points_far_field[0] = 0;
            split_points_far_field[world_size_] = total_pts;

            for (int i = 1; i < world_size_; i++) {

                split_points_far_field[i] = split_points_far_field[i-1] + points_per_rank;

                if (remaining_points > 0) {

                    split_points_far_field[i]++;
                    remaining_points--;

                }

            }

            std::vector<int> recv_counts_far_field(world_size_);
            std::vector<int> displs_far_field(world_size_, 0);

            for (int i = 0; i < world_size_; i++) {

                recv_counts_far_field[i] = split_points_far_field[i+1] - split_points_far_field[i];

                if (i != 0) {

                    displs_far_field[i] = displs_far_field[i-1] + recv_counts_far_field[i-1];

                }

            }

            long long point_low_far_field = split_points_far_field[world_rank_];
            long long point_up_far_field = split_points_far_field[world_rank_+1];

            std::vector<double> xVers_0_loc(point_up_far_field-point_low_far_field);
            std::vector<double> xVers_1_loc(point_up_far_field-point_low_far_field);
            std::vector<double> xVers_2_loc(point_up_far_field-point_low_far_field);

            std::vector<std::complex<double>> far_field_exact_loc(point_up_far_field-point_low_far_field);

            #pragma omp parallel for
            for (long long point = 0; point < point_up_far_field-point_low_far_field; point++) {

                const long long m = (point + point_low_far_field) / N_FAR_PTS[1];
                const long long n = (point + point_low_far_field) % N_FAR_PTS[1];

                const double phi_m = 1.0e-5 + m * deltaPhi;
                const double theta_n = n * deltaTheta;

                xVers_0_loc[point] = std::sin(phi_m) * std::cos(theta_n);
                xVers_1_loc[point] = std::sin(phi_m) * std::sin(theta_n);
                xVers_2_loc[point] = std::cos(phi_m);

                far_field_exact_loc[point] = compute_far_field_exact(xVers_0_loc[point], xVers_1_loc[point], xVers_2_loc[point]);
                
            }

            std::vector<double> xVers_0_all(total_pts);
            std::vector<double> xVers_1_all(total_pts);
            std::vector<double> xVers_2_all(total_pts);
            
            MPI_Allgatherv(&xVers_0_loc[0], point_up_far_field-point_low_far_field, MPI_DOUBLE, &xVers_0_all[0], &recv_counts_far_field[0], &displs_far_field[0], MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgatherv(&xVers_1_loc[0], point_up_far_field-point_low_far_field, MPI_DOUBLE, &xVers_1_all[0], &recv_counts_far_field[0], &displs_far_field[0], MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgatherv(&xVers_2_loc[0], point_up_far_field-point_low_far_field, MPI_DOUBLE, &xVers_2_all[0], &recv_counts_far_field[0], &displs_far_field[0], MPI_DOUBLE, MPI_COMM_WORLD);

            std::vector<double>().swap(xVers_0_loc);
            std::vector<double>().swap(xVers_1_loc);
            std::vector<double>().swap(xVers_2_loc);

            std::vector<std::complex<double>> far_field_approx_loc(total_pts);

            #pragma omp parallel for
            for (long long point = 0; point < total_pts; point++) {
                
                far_field_approx_loc[point] = compute_far_field_approx(phi, xVers_0_all[point], xVers_1_all[point], xVers_2_all[point]);

            }

            std::vector<std::complex<double>> far_field_approx_all(total_pts);

            MPI_Allreduce(&far_field_approx_loc[0], &far_field_approx_all[0], total_pts, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

            std::vector<std::complex<double>>().swap(far_field_approx_loc);

            double error_1_loc = 0.0;
            double error_2_loc = 0.0;
            
            for (long long point = point_low_far_field; point < point_up_far_field; point++) {

                const std::complex<double> far_field_approx = far_field_approx_all[point];
                const std::complex<double> far_field_exact = far_field_exact_loc[point - point_low_far_field];
                
                const double far_field_approx_norm = std::abs(far_field_approx);
                const double far_field_exact_norm = std::abs(far_field_exact);
                    
                const double value1 = std::abs(far_field_exact_norm - far_field_approx_norm);
                const double value2 = std::abs(far_field_exact_norm);

                error_1_loc = std::max(error_1_loc, value1);
                error_2_loc = std::max(error_2_loc, value2);

            }
             
            double error_1;
            double error_2;

            MPI_Allreduce(&error_1_loc, &error_1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&error_2_loc, &error_2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            double end = MPI_Wtime();
            MPI_Barrier(MPI_COMM_WORLD);

            if (world_rank_ == 0) {

                if (timing) {

                    std::cout << "Time compute far field error: " << end - start << " seconds\n";
                    print_max_RSS();

                }

                std::cout << "Far field error: " << std::setprecision(15) << std::fixed << error_1 / error_2 << "\n";

            }

        }

        void compute_near_field(const bool timing, const std::vector<std::complex<double>>& phi)
        {

            const int nNearZones = N_NEAR_PTS.size();

            for (int i = 0; i < nNearZones; i++) {

                MPI_Barrier(MPI_COMM_WORLD);
                double start = MPI_Wtime();

                const double xmin = NEAR_FIELD_LIMITS[i][0];
                const double xmax = NEAR_FIELD_LIMITS[i][1];
                const double ymin = NEAR_FIELD_LIMITS[i][2];
                const double ymax = NEAR_FIELD_LIMITS[i][3];
                const double zmin = NEAR_FIELD_LIMITS[i][4];
                const double zmax = NEAR_FIELD_LIMITS[i][5];

                const int imax = N_NEAR_PTS[i][0];
                const int jmax = N_NEAR_PTS[i][1];
                const int kmax = N_NEAR_PTS[i][2];

                double deltaX = (xmax - xmin) / (imax - 1);
                double deltaY = (ymax - ymin) / (jmax - 1);
                double deltaZ = (zmax - zmin) / (kmax - 1);

                if (imax == 1) deltaX = 0.0;
                if (jmax == 1) deltaY = 0.0;
                if (kmax == 1) deltaZ = 0.0;

                const long long total_pts = imax * jmax * kmax;

                std::vector<long long> split_points_near_field(world_size_ + 1);
                const long long points_per_rank = total_pts / world_size_;
                long long remaining_points = total_pts % world_size_;

                split_points_near_field[0] = 0;
                split_points_near_field[world_size_] = total_pts;

                for (int l = 1; l < world_size_; l++) {

                    split_points_near_field[l] = split_points_near_field[l-1] + points_per_rank;

                    if (remaining_points > 0) {

                        split_points_near_field[l]++;
                        remaining_points--;

                    }

                }

                std::vector<int> recv_counts_near_field(world_size_);
                std::vector<int> displs_near_field(world_size_, 0);

                for (int l = 0; l < world_size_; l++) {

                    recv_counts_near_field[l] = split_points_near_field[l+1] - split_points_near_field[l];

                    if (l != 0) {

                        displs_near_field[l] = displs_near_field[l-1] + recv_counts_near_field[l-1];

                    }

                }

                long long point_low_near_field = split_points_near_field[world_rank_];
                long long point_up_near_field = split_points_near_field[world_rank_+1];

                std::vector<double> x_loc(point_up_near_field-point_low_near_field);
                std::vector<double> y_loc(point_up_near_field-point_low_near_field);
                std::vector<double> z_loc(point_up_near_field-point_low_near_field);

                #pragma omp parallel for
                for (long long point = 0; point < point_up_near_field-point_low_near_field; point++) {

                    const long long ii = (point + point_low_near_field) / (jmax * kmax);
                    const long long jj = ((point + point_low_near_field) % (jmax * kmax)) / kmax;
                    const long long kk = ((point + point_low_near_field) % (jmax * kmax)) % kmax;

                    x_loc[point] = xmin + ii * deltaX;
                    y_loc[point] = ymin + jj * deltaY;
                    z_loc[point] = zmin + kk * deltaZ;

                }

                std::vector<double> x_all(total_pts);
                std::vector<double> y_all(total_pts);
                std::vector<double> z_all(total_pts);
                
                MPI_Allgatherv(&x_loc[0], point_up_near_field-point_low_near_field, MPI_DOUBLE, &x_all[0], &recv_counts_near_field[0], &displs_near_field[0], MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Allgatherv(&y_loc[0], point_up_near_field-point_low_near_field, MPI_DOUBLE, &y_all[0], &recv_counts_near_field[0], &displs_near_field[0], MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Allgatherv(&z_loc[0], point_up_near_field-point_low_near_field, MPI_DOUBLE, &z_all[0], &recv_counts_near_field[0], &displs_near_field[0], MPI_DOUBLE, MPI_COMM_WORLD);

                std::vector<double>().swap(x_loc);
                std::vector<double>().swap(y_loc);
                std::vector<double>().swap(z_loc);

                std::vector<std::complex<double>> u_inc_loc(total_pts);
                std::vector<std::complex<double>> u_scat_loc(total_pts);
                std::vector<std::complex<double>> u_total_loc(total_pts);

                #pragma omp parallel for
                for (long long point = 0; point < total_pts; point++) {

                    u_inc_loc[point] = -compute_incident_field(x_all[point], y_all[point], z_all[point]);

                    u_scat_loc[point] = compute_near_field_approx(phi, x_all[point], y_all[point], z_all[point]);

                    u_total_loc[point] = u_inc_loc[point] + u_scat_loc[point];

                }  
                
                std::vector<std::complex<double>> u_inc_all(total_pts);
                std::vector<std::complex<double>> u_scat_all(total_pts);
                std::vector<std::complex<double>> u_total_all(total_pts);

                MPI_Allreduce(&u_inc_loc[0], &u_inc_all[0], total_pts, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&u_scat_loc[0], &u_scat_all[0], total_pts, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&u_total_loc[0], &u_total_all[0], total_pts, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

                std::vector<std::complex<double>>().swap(u_inc_loc);
                std::vector<std::complex<double>>().swap(u_scat_loc);
                std::vector<std::complex<double>>().swap(u_total_loc);

                double end = MPI_Wtime();
                MPI_Barrier(MPI_COMM_WORLD);

                if (world_rank_ == 0 && timing) {

                    std::cout << "Time compute near field: " << end - start << " seconds\n";
                    print_max_RSS();

                }

                std::ofstream file_i_u_inc, file_i_u_scat, file_i_u_total;

                if (world_rank_ == 0) {

                    file_i_u_inc.open("u_inc_" + std::to_string(i) + ".txt");
                    file_i_u_scat.open("u_scat_" + std::to_string(i) + ".txt");
                    file_i_u_total.open("u_total_" + std::to_string(i) + ".txt");

                    for (int i = 0; i < total_pts; i++) {

                        file_i_u_inc << std::fixed << std::setprecision(16) << x_all[i] << " " << y_all[i] << " " << z_all[i] << " " << u_inc_all[i] << "\n";
                        file_i_u_scat << std::fixed << std::setprecision(16) << x_all[i] << " " << y_all[i] << " " << z_all[i] << " " << u_scat_all[i] << "\n";
                        file_i_u_total << std::fixed << std::setprecision(16) << x_all[i] << " " << y_all[i] << " " << z_all[i] << " " << u_total_all[i] << "\n";

                    }

                    file_i_u_inc.close();
                    file_i_u_scat.close();
                    file_i_u_total.close();

                }

                std::vector<double>().swap(x_all);
                std::vector<double>().swap(y_all);
                std::vector<double>().swap(z_all);

                std::vector<std::complex<double>>().swap(u_inc_all);
                std::vector<std::complex<double>>().swap(u_scat_all);
                std::vector<std::complex<double>>().swap(u_total_all);

            }

        }

        Solver(bool timing)
        {

            setup(timing);

        }

        ~Solver() {           

        }    

};

#endif