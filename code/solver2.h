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
        
        void compute_precomputations()
        {
            
            std::vector<double> vec_mean_size;

            if (DELTA_METHOD == 2) {

                vec_mean_size.resize(patch_up_ - patch_low_);

                #pragma omp parallel for
                for (long long npatch = patch_low_; npatch < patch_up_; npatch++) {

                    std::vector<double> min_vec(3, std::numeric_limits<double>::max());
                    std::vector<double> max_vec(3, std::numeric_limits<double>::lowest());

                    for (int ii = 0; ii < Nu_int_; ii++) {
                        for (int jj = 0; jj < Nv_int_; jj++) {

                            const long long position = (npatch - patch_low_) * Nu_int_*Nv_int_ + ii * Nv_int_ + jj;

                            const double pointx = disc_points_x_all_[position];
                            const double pointy = disc_points_y_all_[position];
                            const double pointz = disc_points_z_all_[position];

                            if (min_vec[0] > pointx) min_vec[0] = pointx;
                            if (min_vec[1] > pointy) min_vec[1] = pointy;
                            if (min_vec[2] > pointz) min_vec[2] = pointz;
                            
                            if (max_vec[0] < pointx) max_vec[0] = pointx;
                            if (max_vec[1] < pointy) max_vec[1] = pointy;
                            if (max_vec[2] < pointz) max_vec[2] = pointz;  

                        }
                    }

                    vec_mean_size[npatch - patch_low_] = ((max_vec[0] - min_vec[0]) +
                                                          (max_vec[1] - min_vec[1]) + 
                                                          (max_vec[2] - min_vec[2])) / 3.0;

                }

            }

            #pragma omp parallel
            {

            std::vector<std::complex<double>> precomputations_thread;
            std::vector<long long> point_precomputations_thread;
            std::vector<long long> patch_num_precomputations_thread; 

            precomputations_thread.reserve(Nu_int_*Nv_int_ * (patch_up_-patch_low_) * 9 * Nu_int_*Nv_int_ / NTHREADS);
            point_precomputations_thread.reserve(Nu_int_*Nv_int_ * (patch_up_-patch_low_) * 9 / NTHREADS);
            patch_num_precomputations_thread.reserve(Nu_int_*Nv_int_ * (patch_up_-patch_low_) * 9 / NTHREADS);
            
            #pragma omp for
            for (long long npoint = point_low_; npoint < point_up_; npoint++) {

                const long long idx_point = USE_ACCELERATOR ? (npoint - point_low_) : npoint;

                const double r_0 = disc_points_x_all_[idx_point];
                const double r_1 = disc_points_y_all_[idx_point];
                const double r_2 = disc_points_z_all_[idx_point];

                const long long patch_num_own = npoint / (Nu_int_*Nv_int_);

                const long long start = start_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];
                const long long total = size_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];

                for (long long i = 0; i < total; i++) {

                    const long long patch_num = sing_and_near_sing_patches_estimate_[start + i];

                    double min_dist_loc;
                    long long indx_point_min_dist_loc;

                    bool sing_near_sing_flag = false;

                    if (patch_num_own == patch_num) {

                        min_dist_loc = 0.0;
                        indx_point_min_dist_loc = npoint % (Nu_int_*Nv_int_);  

                        sing_near_sing_flag = true; 

                    } else {  

                        min_dist_loc = std::numeric_limits<double>::max();

                        for (int ii = 0; ii < Nu_int_; ii++) {
                            for (int jj = 0; jj < Nv_int_; jj++) {

                                double x_loc;
                                double y_loc;
                                double z_loc;

                                if (USE_ACCELERATOR) {

                                    if ((patch_num >= patch_low_) && (patch_num < patch_up_)) {

                                        const long long position = (patch_num - patch_low_) * Nu_int_*Nv_int_ + ii * Nv_int_ + jj;

                                        x_loc = disc_points_x_all_[position];
                                        y_loc = disc_points_y_all_[position];
                                        z_loc = disc_points_z_all_[position];

                                    } else {

                                        x_loc = disc_points_x_not_in_rank_[patch_num][ii * Nv_int_ + jj];
                                        y_loc = disc_points_y_not_in_rank_[patch_num][ii * Nv_int_ + jj];
                                        z_loc = disc_points_z_not_in_rank_[patch_num][ii * Nv_int_ + jj];

                                    }

                                } else {

                                    const long long position = patch_num * Nu_int_*Nv_int_ + ii * Nv_int_ + jj;

                                    x_loc = disc_points_x_all_[position];
                                    y_loc = disc_points_y_all_[position];
                                    z_loc = disc_points_z_all_[position];

                                }

                                const double diff_0 = x_loc - r_0;
                                const double diff_1 = y_loc - r_1;
                                const double diff_2 = z_loc - r_2;

                                const double dist_aux = diff_0*diff_0 + diff_1*diff_1 + diff_2*diff_2;

                                if (dist_aux < min_dist_loc) {

                                    min_dist_loc = dist_aux;
                                    indx_point_min_dist_loc = ii * Nv_int_ + jj;

                                }

                                if (!sing_near_sing_flag && DELTA_METHOD == 1) {

                                    const bool proximity = (std::abs(std::floor(r_0/PROXIMITY_BOX_SIZE) - std::floor(x_loc/PROXIMITY_BOX_SIZE)) <= 1) &&
                                                           (std::abs(std::floor(r_1/PROXIMITY_BOX_SIZE) - std::floor(y_loc/PROXIMITY_BOX_SIZE)) <= 1) &&
                                                           (std::abs(std::floor(r_2/PROXIMITY_BOX_SIZE) - std::floor(z_loc/PROXIMITY_BOX_SIZE)) <= 1);

                                    if (proximity) sing_near_sing_flag = true;

                                }

                            }
                        }

                        min_dist_loc = std::sqrt(min_dist_loc);

                    }

                    if ((DELTA_METHOD == 2) && (min_dist_loc <= vec_mean_size[patch_num_own - patch_low_] * PERCENT_BOX_SIZE)) {

                        sing_near_sing_flag = true;

                    } 
               
                    if (!sing_near_sing_flag) continue;

                    const long long nu_loc = indx_point_min_dist_loc / Nv_int_;
                    const long long nv_loc = indx_point_min_dist_loc % Nv_int_;  

                    int flag_u_loc;
                    int flag_v_loc;

                    if (USE_ACCELERATOR) {

                        if ((patch_num >= patch_low_) && (patch_num < patch_up_)) {

                            flag_u_loc = flags_domain_u_all_[patch_num - patch_low_];
                            flag_v_loc = flags_domain_v_all_[patch_num - patch_low_];                               

                        } else {

                            flag_u_loc = flags_domain_u_not_in_rank_[patch_num];
                            flag_v_loc = flags_domain_v_not_in_rank_[patch_num];
                            
                        }

                    } else {

                        flag_u_loc = flags_domain_u_all_[patch_num];
                        flag_v_loc = flags_domain_v_all_[patch_num];

                    }                        

                    const long long q = patch_num / (Qx_*Qy_);
                    const long long q_x = (patch_num / Qy_) % Qx_;
                    const long long q_y = patch_num % Qy_;

                    const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
                    const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
                    const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
                    const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_; 

                    double ubar_loc, vbar_loc;
                                
                    if (min_dist_loc < 1E-12) {
                        
                        ubar_loc = eta(fejer_nodes_u_int_[nu_loc], flag_u_loc);
                        vbar_loc = eta(fejer_nodes_v_int_[nv_loc], flag_v_loc);

                    } else {

                        double u_a_min, u_b_min, v_a_min, v_b_min;

                        if (nu_loc == 0) {

                            u_a_min = eta(fejer_nodes_u_int_[nu_loc+1], flag_u_loc);
                            u_b_min = 1.0;

                        } else if (nu_loc == Nu_int_-1) {

                            u_a_min = -1.0;
                            u_b_min = eta(fejer_nodes_u_int_[nu_loc-1], flag_u_loc);

                        } else {

                            u_a_min = eta(fejer_nodes_u_int_[nu_loc+1], flag_u_loc);
                            u_b_min = eta(fejer_nodes_u_int_[nu_loc-1], flag_u_loc);

                        }

                        if (nv_loc == 0) {

                            v_a_min = eta(fejer_nodes_v_int_[nv_loc+1], flag_v_loc);
                            v_b_min = 1.0;

                        } else if (nv_loc == Nv_int_-1) {

                            v_a_min = -1.0;
                            v_b_min = eta(fejer_nodes_v_int_[nv_loc-1], flag_v_loc);

                        } else {

                            v_a_min = eta(fejer_nodes_v_int_[nv_loc+1], flag_v_loc);
                            v_b_min = eta(fejer_nodes_v_int_[nv_loc-1], flag_v_loc);

                        } 

                        if (GEOMETRY == 0) {
                        
                            auto dist_func = [r_0, r_1, r_2, q, u_a_loc, u_b_loc, v_a_loc, v_b_loc](double s, double t) -> double {const double ss = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, s);
                                                                                                                                    const double tt = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, t);
                                                                                                                                    double rp_0, rp_1, rp_2;
                                                                                                                                    parametrization_q(ss, tt, q, rp_0, rp_1, rp_2);
                                                                                                                                    return ((r_0 - rp_0)*(r_0 - rp_0) + (r_1 - rp_1)*(r_1 - rp_1) + (r_2 - rp_2)*(r_2 - rp_2));};
                    
                            golden_search(dist_func, G_SEARCH_MAX_ITER, G_SEARCH_TOL, u_a_min, u_b_min, v_a_min, v_b_min, ubar_loc, vbar_loc); 

                        } else {

                            InterpPatch patch = interp_surface_[q];

                            auto dist_func = [r_0, r_1, r_2, q, u_a_loc, u_b_loc, v_a_loc, v_b_loc, patch](double s, double t) -> double {const std::vector<double> ss = {ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, s)};
                                                                                                                                            const std::vector<double> tt = {ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, t)};
                                                                                                                                            double rp_0, rp_1, rp_2;
                                                                                                                                            lagrange_interpolation_2D(patch.uNodes, patch.vNodes, patch.uWeights, patch.vWeights,
                                                                                                                                                                    patch.x, ss, tt, &rp_0);
                                                                                                                                            lagrange_interpolation_2D(patch.uNodes, patch.vNodes, patch.uWeights, patch.vWeights,
                                                                                                                                                                    patch.y, ss, tt, &rp_1);
                                                                                                                                            lagrange_interpolation_2D(patch.uNodes, patch.vNodes, patch.uWeights, patch.vWeights,
                                                                                                                                                                    patch.z, ss, tt, &rp_2);
                                                                                                                                            return ((r_0 - rp_0)*(r_0 - rp_0) + (r_1 - rp_1)*(r_1 - rp_1) + (r_2 - rp_2)*(r_2 - rp_2));};
                    
                            golden_search(dist_func, G_SEARCH_MAX_ITER, G_SEARCH_TOL, u_a_min, u_b_min, v_a_min, v_b_min, ubar_loc, vbar_loc); 
                            
                        }

                    } 

                    std::vector<std::complex<double>> prec_loc(Nu_int_*Nv_int_);
                            
                    beta(r_0, r_1, r_2,
                        q, flag_u_loc, flag_v_loc,
                        u_a_loc, u_b_loc, v_a_loc, v_b_loc,
                        ubar_loc, vbar_loc, 
                        prec_loc);

                    precomputations_thread.insert(precomputations_thread.end(), std::make_move_iterator(prec_loc.begin()), std::make_move_iterator(prec_loc.end()));
                    point_precomputations_thread.push_back(npoint);
                    patch_num_precomputations_thread.push_back(patch_num);

                }

            }

            #pragma omp for ordered
            for (int i = 0; i < NTHREADS; i++) {

                #pragma omp ordered
                {

                    precomputations_.insert(precomputations_.end(), std::make_move_iterator(precomputations_thread.begin()), std::make_move_iterator(precomputations_thread.end()));
                    point_precomputations_.insert(point_precomputations_.end(), std::make_move_iterator(point_precomputations_thread.begin()), std::make_move_iterator(point_precomputations_thread.end()));
                    patch_num_precomputations_.insert(patch_num_precomputations_.end(), std::make_move_iterator(patch_num_precomputations_thread.begin()), std::make_move_iterator(patch_num_precomputations_thread.end()));
                    
                    std::vector<std::complex<double>>().swap(precomputations_thread);
                    std::vector<long long>().swap(point_precomputations_thread);
                    std::vector<long long>().swap(patch_num_precomputations_thread);

                }

            }

            }

            std::vector<long long> patch_num_coeffs_aux = patch_num_precomputations_;
            std::sort(patch_num_coeffs_aux.begin(), patch_num_coeffs_aux.end());
            auto last_unique = std::unique(patch_num_coeffs_aux.begin(), patch_num_coeffs_aux.end());
            patch_num_coeffs_aux.erase(last_unique, patch_num_coeffs_aux.end());
            patch_num_coeffs_ = std::move(patch_num_coeffs_aux);

            std::vector<long long>().swap(sing_and_near_sing_patches_estimate_);
            std::vector<long long>().swap(start_sing_and_near_sing_patches_estimate_);
            std::vector<long long>().swap(size_sing_and_near_sing_patches_estimate_);

        }

        void get_coeffs(const long long q, const std::complex<double>* phi,
                        std::complex<double>* coeffs)
        {

            const int size = Nu_int_ * Nv_int_;

            const double* jac_ptr = nullptr;

            if (USE_ACCELERATOR) {

                if ((q >= patch_low_) && (q < patch_up_)) {

                    jac_ptr = &dsdtjac_all_[(q - patch_low_) * size];

                } else {
                    
                    jac_ptr = dsdtjac_not_in_rank_[q].data(); 

                }

            } else {

                jac_ptr = &dsdtjac_all_[q * size];

            }

            std::vector<std::complex<double>> evals(size);

            #pragma omp simd
            for (int i = 0; i < size; i++) {

                evals[i] = jac_ptr[i] * phi[i];

            }

            std::vector<std::complex<double>> matprod1(Nu_int_*Nv_int_);

            std::complex<double> one(1.0, 0.0);
            std::complex<double> zero(0.0, 0.0);

            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        Nu_int_, Nv_int_, Nu_int_, 
                        &one, &Tn_[0], Nu_int_, &evals[0], Nv_int_, &zero, &matprod1[0], Nv_int_);
            
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        Nu_int_, Nv_int_, Nv_int_,
                        &one, &matprod1[0], Nv_int_, &Tm_[0], Nv_int_, &zero, &coeffs[0], Nv_int_);
            
            const double inv_val = 1.0 / (double)(Nu_int_ * Nv_int_);
            const double scale = 4.0 * inv_val;
        
            #pragma omp simd
            for (int i = 0; i < size; i++) {

                coeffs[i] *= scale;

            }
            
            #pragma omp simd
            for (int m = 0; m < Nv_int_; m++) {

                coeffs[m] *= 0.5;

            }
        
            for (int n = 0; n < Nu_int_; n++) {

                coeffs[n * Nv_int_] *= 0.5;

            }

        }

        void compute_coeffs(std::complex<double>* phi,
                            std::unordered_map<long long, std::vector<std::complex<double>>>& phi_not_in_rank,
                            std::vector<std::complex<double>>& vec_coeffs)
        {
            
            #pragma omp parallel for schedule(dynamic)
            for (long long i = 0; i < patch_num_coeffs_.size(); i++) {

                const long long patch_num = patch_num_coeffs_[i];

                std::complex<double>* phi_loc;

                if ((patch_num >= patch_low_) && (patch_num < patch_up_)) {
                    phi_loc = &phi[(patch_num - patch_low_) * Nu_int_*Nv_int_];
                } else {
                    phi_loc = phi_not_in_rank.at(patch_num).data();
                }

                get_coeffs(patch_num, phi_loc, &vec_coeffs[i * Nu_int_*Nv_int_]);

            }

        }

        void int_near(const std::complex<double>* coeffs,
                      const std::complex<double>* precomputations,
                      std::complex<double>& solution)
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

        void create_IFGF_object() 
        {

            boxes_.InitializeObject(disc_points_x_all_.begin(), disc_points_x_all_.end(),
                                    disc_points_y_all_.begin(), disc_points_y_all_.end(),
                                    disc_points_z_all_.begin(), disc_points_z_all_.end(),
                                    norm_points_x_all_.begin(), norm_points_x_all_.end(),
                                    norm_points_y_all_.begin(), norm_points_y_all_.end(),
                                    norm_points_z_all_.begin(), norm_points_z_all_.end(),
                                    split_points_2_, recv_counts_2_, displs_2_,
                                    coupling_parameter_,
                                    new_order_points_IFGF_);

        }
           
        void set_precomputations_data_IFGF() 
        {

            std::unordered_map<long long, std::unordered_set<long long>> patch_to_sing_point;
            std::vector<std::vector<long long>> points_not_in_rank(world_size_);

            for (long long i = 0; i < point_up_ - point_low_; i++) {

                const long long point = new_order_points_IFGF_[i];

                auto it = std::upper_bound(split_points_2_.begin(), split_points_2_.end(), point);

                int rank = static_cast<int>(std::distance(split_points_2_.begin(), it)) - 1;

                if (rank < 0) {
                
                    rank = 0;
                
                } else if (rank >= world_size_) {
                
                    rank = world_size_ - 1; 
                
                }

                if (world_rank_ == rank) {

                    const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point);
                    const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point);

                    const size_t start_idx = std::distance(point_precomputations_.begin(), point_begin);
                    const size_t end_idx = std::distance(point_precomputations_.begin(), point_end);

                    for (size_t j = start_idx; j < end_idx; ++j) {
                 
                        patch_to_sing_point[patch_num_precomputations_[j]].insert(point);
            
                    }

                } else {

                    points_not_in_rank[rank].push_back(point);

                }

            }

            std::vector<MPI_Count> send_counts_points(world_size_);
            std::vector<MPI_Aint> sdispls_points(world_size_);
            MPI_Count total_send_points = 0;    

            for (int i = 0; i < world_size_; ++i) {

                send_counts_points[i] = static_cast<MPI_Count>(points_not_in_rank[i].size());
                sdispls_points[i] = total_send_points;
                total_send_points += send_counts_points[i];

            }

            std::vector<long long> flat_send_points_buffer(static_cast<size_t>(total_send_points));
            
            MPI_Aint current_point_pos = 0;
    
            for (int i = 0; i < world_size_; ++i) {

                std::copy(points_not_in_rank[i].begin(), points_not_in_rank[i].end(),
                          flat_send_points_buffer.begin() + current_point_pos);
                current_point_pos += send_counts_points[i];

            }

            points_not_in_rank.clear();

            std::vector<MPI_Count> recv_counts_points(world_size_);
    
            MPI_Alltoall(send_counts_points.data(), 1, MPI_COUNT, recv_counts_points.data(), 1, MPI_COUNT, MPI_COMM_WORLD);

            std::vector<MPI_Aint> rdispls_points(world_size_);
            MPI_Count total_recv_points = 0;

            for (int i = 0; i < world_size_; ++i) {

                rdispls_points[i] = total_recv_points;
                total_recv_points += recv_counts_points[i];

            }

            std::vector<long long> flat_recv_points_buffer(static_cast<size_t>(total_recv_points));
            
            MPI_Alltoallv_c(flat_send_points_buffer.data(), send_counts_points.data(), sdispls_points.data(), MPI_LONG_LONG,
                            flat_recv_points_buffer.data(), recv_counts_points.data(), rdispls_points.data(), MPI_LONG_LONG,
                            MPI_COMM_WORLD);

            std::vector<std::vector<MPI_Count>> patch_counts_to_send_back(world_size_);
            std::vector<std::vector<long long>> patches_data_to_send_back(world_size_);

            MPI_Aint current_flat_recv_points_idx = 0;
    
            for (int sender_rank = 0; sender_rank < world_size_; ++sender_rank) {
        
                MPI_Count num_points_from_this_sender = recv_counts_points[sender_rank];

                for (MPI_Count i = 0; i < num_points_from_this_sender; ++i) {
            
                    const long long point = flat_recv_points_buffer[current_flat_recv_points_idx + i];

                    const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point);
                    const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point);

                    const size_t start_idx = std::distance(point_precomputations_.begin(), point_begin);
                    const size_t end_idx = std::distance(point_precomputations_.begin(), point_end);

                    MPI_Count num_patches_for_this_point = static_cast<MPI_Count>(end_idx - start_idx);

                    patch_counts_to_send_back[sender_rank].push_back(num_patches_for_this_point);
            
                    for (size_t j = start_idx; j < end_idx; ++j) {
                
                        patches_data_to_send_back[sender_rank].push_back(patch_num_precomputations_[j]);
                
                    }

                }

                current_flat_recv_points_idx += num_points_from_this_sender;

            }            

            flat_recv_points_buffer.clear();

            std::vector<MPI_Count> send_counts_patch_counts(world_size_);
            std::vector<MPI_Aint> sdispls_patch_counts(world_size_);
            MPI_Count total_send_patch_counts_flat = 0;

            for (int i = 0; i < world_size_; ++i) {

                send_counts_patch_counts[i] = static_cast<MPI_Count>(patch_counts_to_send_back[i].size());
                sdispls_patch_counts[i] = total_send_patch_counts_flat;
                total_send_patch_counts_flat += send_counts_patch_counts[i];

            }

            std::vector<MPI_Count> flat_send_patch_counts_buffer(static_cast<size_t>(total_send_patch_counts_flat));
    
            MPI_Count current_flat_send_pos = 0;

            for (int i = 0; i < world_size_; ++i) {

                std::copy(patch_counts_to_send_back[i].begin(), patch_counts_to_send_back[i].end(),
                        flat_send_patch_counts_buffer.begin() + current_flat_send_pos);
                current_flat_send_pos += send_counts_patch_counts[i];

            }

            patch_counts_to_send_back.clear();

            std::vector<MPI_Count> recv_counts_patch_counts(world_size_);
    
            MPI_Alltoall(send_counts_patch_counts.data(), 1, MPI_COUNT,
                         recv_counts_patch_counts.data(), 1, MPI_COUNT, MPI_COMM_WORLD);

            MPI_Count total_recv_patch_counts_flat = 0;
            std::vector<MPI_Aint> rdispls_patch_counts(world_size_);

            for (int i = 0; i < world_size_; ++i) {

                rdispls_patch_counts[i] = total_recv_patch_counts_flat;
                total_recv_patch_counts_flat += recv_counts_patch_counts[i];

            }

            std::vector<MPI_Count> flat_recv_patch_counts_buffer(static_cast<size_t>(total_recv_patch_counts_flat));

            MPI_Alltoallv_c(flat_send_patch_counts_buffer.data(), send_counts_patch_counts.data(), sdispls_patch_counts.data(), MPI_COUNT,
                            flat_recv_patch_counts_buffer.data(), recv_counts_patch_counts.data(), rdispls_patch_counts.data(), MPI_COUNT,
                            MPI_COMM_WORLD);
    
            flat_send_patch_counts_buffer.clear();

            std::vector<MPI_Count> send_counts_patches_data(world_size_);
            std::vector<MPI_Aint> sdispls_patches_data(world_size_);
            MPI_Count total_send_patches_data_flat = 0;

            for (int i = 0; i < world_size_; ++i) {

                send_counts_patches_data[i] = static_cast<MPI_Count>(patches_data_to_send_back[i].size());
                sdispls_patches_data[i] = total_send_patches_data_flat;
                total_send_patches_data_flat += send_counts_patches_data[i];

            }

            std::vector<long long> flat_send_patches_data_buffer(static_cast<size_t>(total_send_patches_data_flat));
            current_flat_send_pos = 0;

            for (int i = 0; i < world_size_; ++i) {

                std::copy(patches_data_to_send_back[i].begin(), patches_data_to_send_back[i].end(),
                        flat_send_patches_data_buffer.begin() + current_flat_send_pos);
                current_flat_send_pos += send_counts_patches_data[i];

            }

            patches_data_to_send_back.clear();

            std::vector<MPI_Count> recv_counts_patches_data(world_size_);

            MPI_Alltoall(send_counts_patches_data.data(), 1, MPI_COUNT,
                         recv_counts_patches_data.data(), 1, MPI_COUNT, MPI_COMM_WORLD);

            MPI_Count total_recv_patches_data_flat = 0;
            std::vector<MPI_Aint> rdispls_patches_data(world_size_);

            for (int i = 0; i < world_size_; ++i) {

                rdispls_patches_data[i] = total_recv_patches_data_flat;
                total_recv_patches_data_flat += recv_counts_patches_data[i];

            }

            std::vector<long long> flat_recv_patches_data_buffer(static_cast<size_t>(total_recv_patches_data_flat));

            MPI_Alltoallv_c(flat_send_patches_data_buffer.data(), send_counts_patches_data.data(), sdispls_patches_data.data(), MPI_LONG_LONG,
                            flat_recv_patches_data_buffer.data(), recv_counts_patches_data.data(), rdispls_patches_data.data(), MPI_LONG_LONG,
                            MPI_COMM_WORLD);

            flat_send_patches_data_buffer.clear();

            MPI_Count count = 0;

            for (MPI_Count i = 0; i < total_send_points; i++) {

                const long long point = flat_send_points_buffer[i];
                const long long size = flat_recv_patch_counts_buffer[i];

                for (long long j = 0; j < size; j++) {

                    const long long patch = flat_recv_patches_data_buffer[count];

                    patch_to_sing_point[patch].insert(point);

                    count++;

                }

            }

            flat_send_points_buffer.clear();
            flat_recv_patch_counts_buffer.clear();
            flat_recv_patches_data_buffer.clear();

            boxes_.set_precomputations_data(std::move(patch_to_sing_point));

        }

        void compute_new_order_points_RP() {

            std::vector<std::vector<long long>> points_not_in_rank(world_size_), order_points_not_in_rank(world_size_);

            new_order_points_RP_ = std::vector<long long>(point_up_-point_low_);

            for (long long i = point_low_; i < point_up_; i++) {

                const long long point = new_order_points_IFGF_[i - point_low_];

                auto it = std::upper_bound(split_points_2_.begin(), split_points_2_.end(), point);

                int rank = std::clamp<int>(static_cast<int>(std::distance(split_points_2_.begin(), it)) - 1, 0, world_size_-1);

                if (world_rank_ == rank) {

                    new_order_points_RP_[point - point_low_] = i;

                } else {

                    points_not_in_rank[rank].push_back(point);
                    order_points_not_in_rank[rank].push_back(i);

                }

            }

            std::vector<MPI_Count> send_counts_points(world_size_);
            std::vector<MPI_Aint> sdispls_points(world_size_);
            MPI_Count total_send_points = 0;    

            for (int i = 0; i < world_size_; ++i) {

                send_counts_points[i] = static_cast<MPI_Count>(points_not_in_rank[i].size());
                sdispls_points[i] = total_send_points;
                total_send_points += send_counts_points[i];

            }

            std::vector<long long> flat_send_points_buffer(static_cast<size_t>(total_send_points));
            std::vector<long long> flat_send_orders_buffer(static_cast<size_t>(total_send_points));
    
            MPI_Aint current_point_pos = 0;
    
            for (int i = 0; i < world_size_; ++i) {

                std::copy(points_not_in_rank[i].begin(), points_not_in_rank[i].end(),
                          flat_send_points_buffer.begin() + current_point_pos);
                std::copy(order_points_not_in_rank[i].begin(), order_points_not_in_rank[i].end(),
                          flat_send_orders_buffer.begin() + current_point_pos);
                current_point_pos += send_counts_points[i];

            }

            points_not_in_rank.clear();
            order_points_not_in_rank.clear();

            std::vector<MPI_Count> recv_counts_points(world_size_);
    
            MPI_Alltoall(send_counts_points.data(), 1, MPI_COUNT, recv_counts_points.data(), 1, MPI_COUNT, MPI_COMM_WORLD);

            std::vector<MPI_Aint> rdispls_points(world_size_);
            MPI_Count total_recv_points = 0;

            for (int i = 0; i < world_size_; ++i) {

                rdispls_points[i] = total_recv_points;
                total_recv_points += recv_counts_points[i];

            }

            std::vector<long long> flat_recv_points_buffer(static_cast<size_t>(total_recv_points));
            std::vector<long long> flat_recv_orders_buffer(static_cast<size_t>(total_recv_points));

            MPI_Alltoallv_c(flat_send_points_buffer.data(), send_counts_points.data(), sdispls_points.data(), MPI_LONG_LONG,
                            flat_recv_points_buffer.data(), recv_counts_points.data(), rdispls_points.data(), MPI_LONG_LONG,
                            MPI_COMM_WORLD);

            MPI_Alltoallv_c(flat_send_orders_buffer.data(), send_counts_points.data(), sdispls_points.data(), MPI_LONG_LONG, 
                            flat_recv_orders_buffer.data(), recv_counts_points.data(), rdispls_points.data(), MPI_LONG_LONG, 
                            MPI_COMM_WORLD);

            flat_send_points_buffer.clear();
            flat_send_orders_buffer.clear();

            for (long long i = 0; i < total_recv_points; ++i) {
        
                const long long point_global_id = flat_recv_points_buffer[i];
                const long long original_global_index = flat_recv_orders_buffer[i];
        
                new_order_points_RP_[point_global_id - point_low_] = original_global_index;
        
            }

            flat_recv_points_buffer.clear();
            flat_recv_orders_buffer.clear();

        }

        void check_patch_in_neighbours() {

            // Compute point to box

            std::vector<long long> point_to_box(point_up_-point_low_);

            #pragma omp parallel for
            for (long long i = 0; i < point_up_-point_low_; i++) {

                point_to_box[i] = boxes_.get_box(disc_points_x_all_[i], disc_points_y_all_[i], disc_points_z_all_[i]);

            }

            // Compute patch to box

            std::vector<std::vector<long long>> patch_to_box(patch_up_-patch_low_);           

            for (long long patch_num = 0; patch_num < patch_up_-patch_low_; patch_num++) {

                auto &vec = patch_to_box[patch_num];

                vec.assign(point_to_box.begin() + patch_num * Nu_int_*Nv_int_, point_to_box.begin() + (patch_num + 1) * Nu_int_*Nv_int_);

                std::sort(vec.begin(), vec.end());
                auto last_unique = std::unique(vec.begin(), vec.end());
                vec.erase(last_unique, vec.end());

            }

            // Send patch_to_box map to all ranks

            std::vector<long long> patches_loc;
            std::vector<long long> size_loc;
            std::vector<long long> boxes_loc;

            patches_loc.reserve(patch_up_-patch_low_);
            size_loc.reserve(patch_up_-patch_low_);
            
            for (long long patch_num = 0; patch_num < patch_up_-patch_low_; patch_num++) {

                patches_loc.push_back(patch_num + patch_low_);
                size_loc.push_back(patch_to_box[patch_num].size());
                boxes_loc.insert(boxes_loc.end(), patch_to_box[patch_num].begin(), patch_to_box[patch_num].end());

            }

            patch_to_box.clear();

            MPI_Count total_1 = patches_loc.size();
            MPI_Count total_2 = boxes_loc.size();

            std::vector<MPI_Count> recv_counts_1(world_size_);
            std::vector<MPI_Count> recv_counts_2(world_size_);

            MPI_Allgather(&total_1, 1, MPI_COUNT, &recv_counts_1[0], 1, MPI_COUNT, MPI_COMM_WORLD);
            MPI_Allgather(&total_2, 1, MPI_COUNT, &recv_counts_2[0], 1, MPI_COUNT, MPI_COMM_WORLD);

            std::vector<MPI_Aint> displs_1(world_size_, 0);
            std::vector<MPI_Aint> displs_2(world_size_, 0);

            MPI_Count total_1_all = 0;
            MPI_Count total_2_all = 0;

            for (int i = 0; i < world_size_; i++) {

                displs_1[i] = total_1_all;
                displs_2[i] = total_2_all;

                total_1_all += recv_counts_1[i];
                total_2_all += recv_counts_2[i];

            }

            std::vector<long long> patches_all(total_1_all);
            std::vector<long long> size_all(total_1_all);
            std::vector<long long> boxes_all(total_2_all);

            MPI_Allgatherv_c(&patches_loc[0], total_1, MPI_LONG_LONG, &patches_all[0], &recv_counts_1[0], &displs_1[0], MPI_LONG_LONG, MPI_COMM_WORLD);
            MPI_Allgatherv_c(&size_loc[0], total_1, MPI_LONG_LONG, &size_all[0], &recv_counts_1[0], &displs_1[0], MPI_LONG_LONG, MPI_COMM_WORLD);
            MPI_Allgatherv_c(&boxes_loc[0], total_2, MPI_LONG_LONG, &boxes_all[0], &recv_counts_2[0], &displs_2[0], MPI_LONG_LONG, MPI_COMM_WORLD);

            std::vector<long long>().swap(patches_loc);
            std::vector<long long>().swap(size_loc);
            std::vector<long long>().swap(boxes_loc);

            std::vector<std::vector<long long>> patch_to_box_all(Q_ * Qx_*Qy_);   

            long long counter = 0;

            for (long long i = 0; i < Q_ * Qx_*Qy_; i++) {

                long long patch = patches_all[i];
                long long size_boxes = size_all[i];

                auto &vec = patch_to_box_all[i];

                vec.assign(boxes_all.begin() + counter, boxes_all.begin() + counter + size_boxes);

                counter += size_boxes;

            }

            std::vector<long long>().swap(patches_all);
            std::vector<long long>().swap(size_all);
            std::vector<long long>().swap(boxes_all);

            // Compute relevant morton boxes

            std::unordered_set<long long> mortonidofrelboxes;

            for (long long i = 0; i < Q_ * Qx_*Qy_; i++) {

                mortonidofrelboxes.insert(patch_to_box_all[i].begin(), patch_to_box_all[i].end());

            }    

            // Get size boxes level D

            std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_all;

            if (USE_ADAPTIVITY) {

                mortonbox2discretizationpoints_all = boxes_.get_mortonbox2discretizationpoints_all();

            }

            // Check if patch is in neighbour union

            #pragma omp parallel
            {

            long long actual_npoint = -1;
            std::vector<long long> neighbours_box_npoint;

            #pragma omp for
            for (long long i = 0; i < point_precomputations_.size(); i++) {

                const long long npoint = point_precomputations_[i];

                if (actual_npoint != npoint) {
                
                    const long long box_npoint = point_to_box[npoint - point_low_];
                    neighbours_box_npoint = boxes_.get_neighbours_box(box_npoint);

                    std::vector<long long> relevant_neighbours;
                    for (const auto & elem : neighbours_box_npoint) {

                        if (mortonidofrelboxes.find(elem) != mortonidofrelboxes.end()) {

                            relevant_neighbours.push_back(elem);

                        }

                    }

                    neighbours_box_npoint = std::move(relevant_neighbours);

                    actual_npoint = npoint;

                }

                const long long patch_num_singular = patch_num_precomputations_[i];
                const std::vector<long long> boxes_patch_num_singular = patch_to_box_all[patch_num_singular];

                for (const auto & box : boxes_patch_num_singular) {

                    bool is_in_adaptive_boxes = true;
                    
                    if (USE_ADAPTIVITY) {

                        is_in_adaptive_boxes = (mortonbox2discretizationpoints_all[box][1] <= MAX_ELEMS_LEAF);

                    }                    
                    
                    bool is_in_neighbour_union = (std::find(neighbours_box_npoint.begin(), neighbours_box_npoint.end(), box) != neighbours_box_npoint.end());

                    if (!(is_in_neighbour_union || is_in_adaptive_boxes)) {

                        std::cout << "Point " << npoint << " and singular / near singular patches are not contained in neighbor union.\n";
                        std::exit(0);

                    }

                }

            }

            }

        }

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
                                       std::complex<double>* intensities) 
        {        

            const long long total_points = Nu_int_ * Nv_int_;
            const long long patch_offset = npatch * total_points;

            const double* dsdtjac_ptr = &dsdtjac_all_[patch_offset];
            const double* w_ptr   = fejer_weights_u_v_int_.data();

            #pragma omp simd
            for (long long k = 0; k < total_points; k++) {

                double scalar = dsdtjac_ptr[k] * w_ptr[k];
                intensities[k] = scalar * phi[k];

            }   

        }

        void compute_intensities(const std::complex<double>* phi,
                                 std::vector<std::complex<double>>& intensities)
        {

            const long long num_local_points = point_up_ - point_low_;
            const long long num_local_patches = patch_up_ - patch_low_;
            const int patch_size = Nu_int_ * Nv_int_;

            intensities.resize(num_local_points);

            const long long* points_ptr = new_order_points_RP_.data();

            int max_threads = omp_get_max_threads();

            std::vector<std::vector<MPI_Count>> thread_counts(max_threads, std::vector<MPI_Count>(world_size_, 0));

            #pragma omp parallel
            {

            int tid = omp_get_thread_num();
            auto& my_counts = thread_counts[tid];

            #pragma omp for
            for (long long patch_idx = 0; patch_idx < num_local_patches; ++patch_idx) {

                long long start_idx = patch_idx * patch_size;

                for (int k = 0; k < patch_size; ++k) {

                    long long point_idx = start_idx + k;

                    long long point = points_ptr[point_idx];

                    auto it = std::upper_bound(split_points_2_.begin(), split_points_2_.end(), point);
                    int rank = static_cast<int>(std::distance(split_points_2_.begin(), it)) - 1;
                
                    if (rank < 0) rank = 0;
                    if (rank >= world_size_) rank = world_size_ - 1;
                        
                    my_counts[rank]++;

                }

            }

            }

            std::vector<MPI_Count> global_send_counts(world_size_, 0);
            std::vector<MPI_Aint> sdispls(world_size_);
            std::vector<std::vector<MPI_Aint>> thread_offsets(max_threads, std::vector<MPI_Aint>(world_size_));

            MPI_Count total_send_points = 0;

            for (int r = 0; r < world_size_; r++) {
            
                sdispls[r] = total_send_points;
            
                for (int t = 0; t < max_threads; t++) {

                    thread_offsets[t][r] = total_send_points;
                    global_send_counts[r] += thread_counts[t][r];
                    total_send_points += thread_counts[t][r];

                }
            
            }

            MPI_Count local_count = global_send_counts[world_rank_];
            global_send_counts[world_rank_] = 0;

            std::vector<long long> flat_send_points(total_send_points);
            std::vector<std::complex<double>> flat_send_data(total_send_points);

            long long* send_pts_ptr = flat_send_points.data();
            std::complex<double>* send_data_ptr = flat_send_data.data();
            std::complex<double>* intensities_ptr = intensities.data();


            #pragma omp parallel
            {

            int tid = omp_get_thread_num();
            std::vector<std::complex<double>> patch_buffer(patch_size);
            std::complex<double>* patch_buf_ptr = patch_buffer.data();
            auto& my_offsets = thread_offsets[tid];

            #pragma omp for
            for (long long patch_idx = 0; patch_idx < num_local_patches; ++patch_idx) {

                compute_intensities_patch(patch_idx, 
                                        &phi[patch_idx * patch_size],
                                        patch_buf_ptr);

                long long start_idx = patch_idx * patch_size;

                for (int k = 0; k < patch_size; ++k) {
                
                    long long point_idx = start_idx + k;

                    long long point = points_ptr[point_idx];
                
                    std::complex<double> val = patch_buf_ptr[k];

                    auto it = std::upper_bound(split_points_2_.begin(), split_points_2_.end(), point);
                    int rank = static_cast<int>(std::distance(split_points_2_.begin(), it)) - 1;
                    if (rank < 0) rank = 0;
                    if (rank >= world_size_) rank = world_size_ - 1;
                        
                    if (rank == world_rank_) {

                        intensities_ptr[point - point_low_] = val;

                    } else {
                        
                        MPI_Aint pos = my_offsets[rank]++;
                        send_pts_ptr[pos] = point;
                        send_data_ptr[pos] = val;

                    }
                    
                }
                    
            }
            
            }

            std::vector<MPI_Count> recv_counts(world_size_);
            MPI_Alltoall(global_send_counts.data(), 1, MPI_COUNT, recv_counts.data(), 1, MPI_COUNT, MPI_COMM_WORLD);

            std::vector<MPI_Aint> rdispls(world_size_);
            MPI_Count total_recv_points = 0;

            for (int i = 0; i < world_size_; ++i) {

                rdispls[i] = total_recv_points;
                total_recv_points += recv_counts[i];

            }
            
            std::vector<long long> flat_recv_points(total_recv_points);
            std::vector<std::complex<double>> flat_recv_data(total_recv_points);

            MPI_Alltoallv_c(flat_send_points.data(), global_send_counts.data(), sdispls.data(), MPI_LONG_LONG,
                            flat_recv_points.data(), recv_counts.data(), rdispls.data(), MPI_LONG_LONG,
                            MPI_COMM_WORLD);
            std::vector<long long>().swap(flat_send_points);

            MPI_Alltoallv_c(flat_send_data.data(), global_send_counts.data(), sdispls.data(), MPI_DOUBLE_COMPLEX,
                            flat_recv_data.data(), recv_counts.data(), rdispls.data(), MPI_DOUBLE_COMPLEX,
                            MPI_COMM_WORLD);
            std::vector<std::complex<double>>().swap(flat_send_data);

            long long* recv_pts_ptr = flat_recv_points.data();
            std::complex<double>* recv_data_ptr = flat_recv_data.data();

            #pragma omp parallel for
            for (MPI_Count i = 0; i < total_recv_points; i++) {

                intensities_ptr[recv_pts_ptr[i] - point_low_] = recv_data_ptr[i];
            
            }

        }

        void compute_integral_acc(const std::complex<double>* phi,
                                  const std::vector<std::complex<double>>& vec_coeffs,
                                  std::complex<double>* integral)
        {

            #pragma omp parallel for schedule(dynamic)
            for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

                const double r_0 = disc_points_x_all_[npoint];
                const double r_1 = disc_points_y_all_[npoint];
                const double r_2 = disc_points_z_all_[npoint];

                integral[npoint] = std::complex<double>(0.0, 0.0);

                const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);
                const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);                      

                const long long start_idx_precomp = std::distance(point_precomputations_.begin(), point_begin);
                const long long end_idx_precomp = std::distance(point_precomputations_.begin(), point_end);

                for (long long i = start_idx_precomp; i < end_idx_precomp; ++i) {

                    const long long patch_num = patch_num_precomputations_[i];

                    std::complex<double> solution(0.0, 0.0);    

                    const long long start_precomputations = i * Nu_int_*Nv_int_;
                    
                    const auto pos_coeffs = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
                    const auto idx_coeffs = std::distance(patch_num_coeffs_.begin(), pos_coeffs);

                    const long long start_coeffs = idx_coeffs * Nu_int_*Nv_int_;

                    int_near(&vec_coeffs[start_coeffs], 
                             &precomputations_[start_precomputations], 
                             solution);                    

                    integral[npoint] += solution;

                }

                if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

                    integral[npoint] += 0.5 * phi[npoint];

                }

            }

        }

        void redistribute_data_RP(const std::vector<std::complex<double>>& intensities, std::complex<double>* rhs)
        {

            std::vector<std::vector<long long>> points_not_in_rank(world_size_);
            std::vector<std::vector<std::complex<double>>> intensities_not_in_rank(world_size_);

            for (long long i = 0; i < point_up_ - point_low_; i++) {

                const long long point = new_order_points_IFGF_[i];

                auto it = std::upper_bound(split_points_2_.begin(), split_points_2_.end(), point);

                int rank = std::clamp<int>(static_cast<int>(std::distance(split_points_2_.begin(), it)) - 1, 0, world_size_-1);

                if (world_rank_ == rank) {

                    rhs[point - point_low_] += intensities[i];

                } else {

                    points_not_in_rank[rank].push_back(point);
                    intensities_not_in_rank[rank].push_back(intensities[i]);

                }

            }

            std::vector<MPI_Count> send_counts_points(world_size_);
            std::vector<MPI_Aint> sdispls_points(world_size_);
            MPI_Count total_send_points = 0;    

            for (int i = 0; i < world_size_; ++i) {

                send_counts_points[i] = static_cast<MPI_Count>(points_not_in_rank[i].size());
                sdispls_points[i] = total_send_points;
                total_send_points += send_counts_points[i];

            }

            std::vector<long long> flat_send_points_buffer(static_cast<size_t>(total_send_points));
            std::vector<std::complex<double>> flat_send_intensities_buffer(static_cast<size_t>(total_send_points));
            
            MPI_Aint current_point_pos = 0;
    
            for (int i = 0; i < world_size_; ++i) {

                std::copy(points_not_in_rank[i].begin(), points_not_in_rank[i].end(),
                          flat_send_points_buffer.begin() + current_point_pos);
                std::copy(intensities_not_in_rank[i].begin(), intensities_not_in_rank[i].end(),
                          flat_send_intensities_buffer.begin() + current_point_pos);
                current_point_pos += send_counts_points[i];

            }

            points_not_in_rank.clear();
            intensities_not_in_rank.clear();

            std::vector<MPI_Count> recv_counts_points(world_size_);
    
            MPI_Alltoall(send_counts_points.data(), 1, MPI_COUNT, recv_counts_points.data(), 1, MPI_COUNT, MPI_COMM_WORLD);

            std::vector<MPI_Aint> rdispls_points(world_size_);
            MPI_Count total_recv_points = 0;

            for (int i = 0; i < world_size_; ++i) {

                rdispls_points[i] = total_recv_points;
                total_recv_points += recv_counts_points[i];

            }

            std::vector<long long> flat_recv_points_buffer(static_cast<size_t>(total_recv_points));
            std::vector<std::complex<double>> flat_recv_intensities_buffer(static_cast<size_t>(total_recv_points));
            
            MPI_Alltoallv_c(flat_send_points_buffer.data(), send_counts_points.data(), sdispls_points.data(), MPI_LONG_LONG,
                            flat_recv_points_buffer.data(), recv_counts_points.data(), rdispls_points.data(), MPI_LONG_LONG,
                            MPI_COMM_WORLD);

            MPI_Alltoallv_c(flat_send_intensities_buffer.data(), send_counts_points.data(), sdispls_points.data(), MPI_DOUBLE_COMPLEX,
                            flat_recv_intensities_buffer.data(), recv_counts_points.data(), rdispls_points.data(), MPI_DOUBLE_COMPLEX,
                            MPI_COMM_WORLD);

            std::vector<long long>().swap(flat_send_points_buffer);
            std::vector<std::complex<double>>().swap(flat_send_intensities_buffer);

            #pragma omp parallel for
            for (MPI_Count i = 0; i < total_recv_points; i++) {

                long long point = flat_recv_points_buffer[i];
                std::complex<double> intensity = flat_recv_intensities_buffer[i];

                rhs[point - point_low_] += intensity;

            }

            std::vector<long long>().swap(flat_recv_points_buffer);
            std::vector<std::complex<double>>().swap(flat_recv_intensities_buffer);

        }

        void iterator_function_acc(std::complex<double>* phi, std::complex<double>* rhs)
        {

            std::unordered_map<long long, std::vector<std::complex<double>>> phi_not_in_rank;

            compute_phi_not_in_rank(phi, phi_not_in_rank);

            std::vector<std::complex<double>> vec_coeffs(patch_num_coeffs_.size() * Nu_int_*Nv_int_);

            compute_coeffs(phi, phi_not_in_rank, vec_coeffs);

            phi_not_in_rank.clear();

            compute_integral_acc(phi, vec_coeffs, rhs);      

            std::vector<std::complex<double>>().swap(vec_coeffs);

            std::vector<std::complex<double>> intensities;

            compute_intensities(phi, intensities);

            boxes_.Solve<&fct_4, &fac_1>(intensities);

            redistribute_data_RP(intensities, rhs);

            std::vector<std::complex<double>>().swap(intensities);

        }

        void compute_phi_not_in_rank(const std::complex<double>* phi, std::unordered_map<long long, std::vector<std::complex<double>>>& phi_not_in_rank)
        {

            std::vector<std::vector<long long>> patches_not_in_rank(world_size_);

            for (long long patch_num : patch_num_coeffs_) {

                auto it = std::upper_bound(split_points_.begin(), split_points_.end(), patch_num);

                int rank = static_cast<int>(std::distance(split_points_.begin(), it)) - 1;

                if (rank < 0) rank = 0;                
                if (rank >= world_size_) rank = world_size_ - 1; 

                if (world_rank_ != rank) {

                    patches_not_in_rank[rank].push_back(patch_num);

                }

            }

            std::vector<MPI_Count> send_counts_patches(world_size_);
            std::vector<MPI_Aint> sdispls_patches(world_size_);
            MPI_Count total_send_patches = 0;    

            for (int i = 0; i < world_size_; ++i) {

                send_counts_patches[i] = static_cast<MPI_Count>(patches_not_in_rank[i].size());
                sdispls_patches[i] = total_send_patches;
                total_send_patches += send_counts_patches[i];

            }

            std::vector<long long> flat_send_patches_buffer(static_cast<size_t>(total_send_patches));
        
            MPI_Aint current_patch_pos = 0;
    
            for (int i = 0; i < world_size_; ++i) {

                std::copy(patches_not_in_rank[i].begin(), patches_not_in_rank[i].end(),
                          flat_send_patches_buffer.begin() + current_patch_pos);
                current_patch_pos += send_counts_patches[i];

            }
            
            patches_not_in_rank.clear();

            std::vector<MPI_Count> recv_counts_patches(world_size_);

            MPI_Alltoall(send_counts_patches.data(), 1, MPI_COUNT, recv_counts_patches.data(), 1, MPI_COUNT, MPI_COMM_WORLD);

            std::vector<MPI_Aint> rdispls_patches(world_size_);
            MPI_Count total_recv_patches = 0;

            for (int i = 0; i < world_size_; ++i) {

                rdispls_patches[i] = total_recv_patches;
                total_recv_patches += recv_counts_patches[i];

            }

            std::vector<long long> flat_recv_patches_buffer(static_cast<size_t>(total_recv_patches));
        
            MPI_Request request_patches;
            MPI_Ialltoallv_c(flat_send_patches_buffer.data(), send_counts_patches.data(), sdispls_patches.data(), MPI_LONG_LONG,
                             flat_recv_patches_buffer.data(), recv_counts_patches.data(), rdispls_patches.data(), MPI_LONG_LONG,
                             MPI_COMM_WORLD, &request_patches);

            const size_t patch_data_size = Nu_int_ * Nv_int_; 
            MPI_Count total_send_data_back = total_recv_patches * static_cast<MPI_Count>(patch_data_size);

            std::vector<std::complex<double>> flat_send_data_back_buffer(static_cast<size_t>(total_send_data_back));
    
            std::vector<MPI_Count> send_counts_data_back(world_size_);
            std::vector<MPI_Aint> sdispls_data_back(world_size_);
            MPI_Aint current_data_pos = 0;

            for (int sender_rank = 0; sender_rank < world_size_; ++sender_rank) {

                MPI_Count num_patches = recv_counts_patches[sender_rank];
                MPI_Count data_size = num_patches * static_cast<MPI_Count>(patch_data_size);

                send_counts_data_back[sender_rank] = data_size;
                sdispls_data_back[sender_rank] = current_data_pos;
                current_data_pos += data_size;

            }

            MPI_Wait(&request_patches, MPI_STATUS_IGNORE); 

            MPI_Aint current_patch_idx = 0;
    
            for (int sender_rank = 0; sender_rank < world_size_; ++sender_rank) {

                MPI_Count num_patches_from_this_sender = recv_counts_patches[sender_rank];

                size_t dest_start_offset = sdispls_data_back[sender_rank]; 

                for (MPI_Count i = 0; i < num_patches_from_this_sender; ++i) {

                    const long long patch = flat_recv_patches_buffer[current_patch_idx + i];
                    const long long patch_scaled = patch - patch_low_;
                    
                    size_t source_start_idx = static_cast<size_t>(patch_scaled) * patch_data_size;
                    size_t dest_current_idx = dest_start_offset + i * patch_data_size;
                    
                    std::copy(phi + source_start_idx, 
                              phi + source_start_idx + patch_data_size,
                              flat_send_data_back_buffer.begin() + dest_current_idx);

                }

                current_patch_idx += num_patches_from_this_sender;

            }

            flat_recv_patches_buffer.clear(); 

            std::vector<MPI_Count> recv_counts_data_back(world_size_);

            MPI_Alltoall(send_counts_data_back.data(), 1, MPI_COUNT,
                         recv_counts_data_back.data(), 1, MPI_COUNT, MPI_COMM_WORLD);

            MPI_Count total_recv_data_back = 0;
            std::vector<MPI_Aint> rdispls_data_back(world_size_);

            for (int i = 0; i < world_size_; ++i) {

                rdispls_data_back[i] = total_recv_data_back;
                total_recv_data_back += recv_counts_data_back[i];

            }

            std::vector<std::complex<double>> flat_recv_data_buffer_back(static_cast<size_t>(total_recv_data_back));

            MPI_Request request_data;
            MPI_Ialltoallv_c(flat_send_data_back_buffer.data(), send_counts_data_back.data(), sdispls_data_back.data(), MPI_DOUBLE_COMPLEX,
                             flat_recv_data_buffer_back.data(), recv_counts_data_back.data(), rdispls_data_back.data(), MPI_DOUBLE_COMPLEX,
                             MPI_COMM_WORLD, &request_data); 

            flat_send_data_back_buffer.clear();

            MPI_Wait(&request_data, MPI_STATUS_IGNORE);

            MPI_Count count = 0;
            
            for (MPI_Count i = 0; i < total_send_patches; i++) {

                const long long patch = flat_send_patches_buffer[i];

                size_t data_start_index = static_cast<size_t>(i) * patch_data_size;

                std::vector<std::complex<double>> phi_local(
                    flat_recv_data_buffer_back.begin() + data_start_index,
                    flat_recv_data_buffer_back.begin() + data_start_index + patch_data_size
                );

                phi_not_in_rank[patch] = std::move(phi_local); 

            }

            flat_send_patches_buffer.clear();
            flat_recv_data_buffer_back.clear();   

        }

        void iterator_function(std::complex<double>* phi, std::complex<double>* rhs)
        {

            iterator_function_acc(phi, rhs);

        }

        void setup(bool timing)
        {

            MPI_Barrier(MPI_COMM_WORLD);

            double start_1 = MPI_Wtime();

            compute_parallel_parameters();

            double end_1 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time compute parallel parameters: " << end_1 - start_1 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            double start_2 = MPI_Wtime();

            if (GEOMETRY != 0) {

                load_interpolated_surface();

            }

            double end_2 = MPI_Wtime();

            if (GEOMETRY != 0 && timing && world_rank_ == 0) {

                std::cout << "Time load interpolated surface: " << end_2 - start_2 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            double start_3 = MPI_Wtime();

            compute_fejer_nodes_and_weights();

            double end_3 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time compute Fejer nodes and weights: " << end_3 - start_3 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            double start_4 = MPI_Wtime();

            compute_chebyshev_evaluations();

            double end_4 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time compute Chebyshev evaluations: " << end_4 - start_4 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            double start_5 = MPI_Wtime();

            compute_flags_domain();

            double end_5 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time compute flags domain: " << end_5 - start_5 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            double start_6 = MPI_Wtime();

            compute_discretization_domain();

            double end_6 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time compute discretization domain: " << end_6 - start_6 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            double start_7 = MPI_Wtime();

            compute_coupling_parameter();  

            double end_7 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time compute coupling parameter: " << end_7 - start_7 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);  

            double start_8 = MPI_Wtime();   

            compute_near_singular_patches_estimate();   

            double end_8 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time compute near singular patches estimate: " << end_8 - start_8 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);  

            double start_9 = MPI_Wtime();

            compute_precomputations();

            double end_9 = MPI_Wtime();

            if (timing && world_rank_ == 0) {

                std::cout << "Time compute precomputations data: " << end_9 - start_9 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            double start_10 = MPI_Wtime();

            if (USE_ACCELERATOR) {

                create_IFGF_object();

            }

            double end_10 = MPI_Wtime();

            if (USE_ACCELERATOR && timing && world_rank_ == 0) {

                std::cout << "Time create IFGF object: " << end_10 - start_10 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

            double start_11 = MPI_Wtime();

            if (USE_ACCELERATOR) {

                set_precomputations_data_IFGF();

            }

            double end_11 = MPI_Wtime();

            if (USE_ACCELERATOR && timing && world_rank_ == 0) {

                std::cout << "Time set precomputation data in IFGF: " << end_11 - start_11 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD); 

            double start_12 = MPI_Wtime();

            if (USE_ACCELERATOR) {

                compute_new_order_points_RP();                

            }

            double end_12 = MPI_Wtime();

            if (USE_ACCELERATOR && timing && world_rank_ == 0) {

                std::cout << "Time compute new order points RP: " << end_12 - start_12 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD); 

            double start_13 = MPI_Wtime();

            if (USE_ACCELERATOR) {

                check_patch_in_neighbours();

            }

            double end_13 = MPI_Wtime();

            if (USE_ACCELERATOR && timing && world_rank_ == 0) {

                std::cout << "Time check patch in neighbours: " << end_13 - start_13 << " seconds\n";
                print_max_RSS();

            }

            MPI_Barrier(MPI_COMM_WORLD);

        }

        struct UserContext{
            std::function<void(std::complex<double>*, std::complex<double>*)> func;
            long long n;
        };

        static PetscErrorCode iterator_function_2(Mat A, Vec x, Vec y) 
        {

            PetscFunctionBeginUser;
            
            UserContext *ctx;
            PetscErrorCode ierr;          
    
            ierr = MatShellGetContext(A, &ctx);
            CHKERRQ(ierr);

            PetscScalar *x_array;
            PetscScalar *y_array;

            ierr = VecGetArrayRead(x, &x_array); 
            CHKERRQ(ierr);
            ierr = VecGetArray(y, &y_array); 
            CHKERRQ(ierr);

            std::complex<double> *complex_x_ptr = reinterpret_cast<std::complex<double>*>(x_array);
            std::complex<double> *complex_y_ptr = reinterpret_cast<std::complex<double>*>(y_array);
    
            ctx->func(complex_x_ptr, complex_y_ptr);

            ierr = VecRestoreArrayRead(x, &x_array); 
            CHKERRQ(ierr);
            ierr = VecRestoreArray(y, &y_array); 
            CHKERRQ(ierr);

            PetscFunctionReturn(0);

        }

        std::vector<std::complex<double>> solve(const std::vector<std::complex<double>>& rhs) 
        {

            KSP ksp;
            Mat A;
            Vec x, b;
            PetscErrorCode ierr; 

            long long N = static_cast<long long>(Q_) * Qx_*Qy_ * Nu_int_*Nv_int_;
            long long local_size = point_up_ - point_low_;

            UserContext user_ctx;

            user_ctx.func = [this](std::complex<double>* x, std::complex<double>* solution) -> void {iterator_function(x, solution);};
            user_ctx.n = N;

            ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, local_size, N, (const PetscScalar*)rhs.data(), &b);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);

            ierr = VecDuplicate(b, &x);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            
            ierr = MatCreateShell(PETSC_COMM_WORLD, local_size, local_size, N, N, &user_ctx, &A);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)(void))iterator_function_2);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);  
            
            ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = KSPSetOperators(ksp, A, A); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = KSPSetType(ksp, KSPGMRES); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = KSPSetTolerances(ksp, TOL_GMRES, PETSC_DEFAULT, PETSC_DEFAULT, MAX_ITER); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = KSPSetFromOptions(ksp); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            
            ierr = KSPSolve(ksp, b, x); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);

            KSPConvergedReason reason;
            ierr = KSPGetConvergedReason(ksp, &reason); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            if (reason < 0) {
                PetscPrintf(PETSC_COMM_WORLD, "GMRES diverged with reason %d\n", reason);
            } else {
                PetscInt its;
                PetscReal rnorm;
                ierr = KSPGetIterationNumber(ksp, &its); 
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
                ierr = KSPGetResidualNorm(ksp, &rnorm); 
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
                PetscPrintf(PETSC_COMM_WORLD, "GMRES converged in %d iterations with residual norm %g\n", its, (double)rnorm);
            }
            
            const PetscScalar *x_array;
            ierr = VecGetArrayRead(x, &x_array); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);

            const std::complex<double> *complex_ptr = reinterpret_cast<const std::complex<double>*>(x_array);
            std::vector<std::complex<double>> x_sol(complex_ptr, complex_ptr + local_size);

            ierr = VecRestoreArrayRead(x, &x_array); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);

            VecDestroy(&x);
            VecDestroy(&b);
            MatDestroy(&A);
            KSPDestroy(&ksp);

            return x_sol;

        }

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