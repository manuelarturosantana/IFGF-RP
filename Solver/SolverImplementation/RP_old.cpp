#include "../solver2.h"

using namespace Eigen;

void Solver::compute_precomputations()
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

                    const long long position = npatch * Nu_int_*Nv_int_ + ii * Nv_int_ + jj;

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

        long long idx_point;
        idx_point = npoint;

        const double r_0 = disc_points_x_all_[idx_point];
        const double r_1 = disc_points_y_all_[idx_point];
        const double r_2 = disc_points_z_all_[idx_point];

        const long long patch_num_own = npoint / (Nu_int_*Nv_int_);

        const long long start = start_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];
        const long long total = size_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];

        for (long long i = 0; i < total; i++) {

            const long long patch_num = sing_and_near_sing_patches_estimate_[start + i];

            double min_dist_loc;
            int indx_point_min_dist_loc;

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

                        const long long position = patch_num * Nu_int_*Nv_int_ + ii * Nv_int_ + jj;

                        x_loc = disc_points_x_all_[position];
                        y_loc = disc_points_y_all_[position];
                        z_loc = disc_points_z_all_[position];

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

                            if (proximity) {

                                sing_near_sing_flag = true;

                            }

                        }

                    }
                }

                min_dist_loc = std::sqrt(min_dist_loc);

            }

            if ((DELTA_METHOD == 2) && (min_dist_loc <= vec_mean_size[patch_num_own - patch_low_] * PERCENT_BOX_SIZE)) {

                sing_near_sing_flag = true;

            }     
        
            if (sing_near_sing_flag) {  

                const int nu_loc = indx_point_min_dist_loc / Nv_int_;
                const int nv_loc = indx_point_min_dist_loc % Nv_int_;  

                int flag_u_loc;
                int flag_v_loc;

                flag_u_loc = flags_domain_u_all_[patch_num];
                flag_v_loc = flags_domain_v_all_[patch_num];

                const long long q = patch_num / (Qx_*Qy_);
                const int q_x = (patch_num / Qy_) % Qx_;
                const int q_y = patch_num % Qy_;

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
                    
                        auto dist_func = [r_0, r_1, r_2, q, u_a_loc, u_b_loc, v_a_loc, v_b_loc, this](double s, double t) -> double {const double ss = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, s);
                                                                                                                                const double tt = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, t);
                                                                                                                                double rp_0, rp_1, rp_2;
                                                                                                                                parametrization_q(this->SPHERE_RADIUS, this->SPHERE_CENTER,ss, tt, q, rp_0, rp_1, rp_2);
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

                precomputations_thread.insert(precomputations_thread.end(), prec_loc.begin(), prec_loc.end());
                point_precomputations_thread.push_back(npoint);
                patch_num_precomputations_thread.push_back(patch_num);

            }

        }

    }

    #pragma omp for ordered
    for (int i = 0; i < NTHREADS; i++) {

        #pragma omp ordered
        {

            precomputations_.insert(precomputations_.end(), precomputations_thread.begin(), precomputations_thread.end());
            point_precomputations_.insert(point_precomputations_.end(), point_precomputations_thread.begin(), point_precomputations_thread.end());
            patch_num_precomputations_.insert(patch_num_precomputations_.end(), patch_num_precomputations_thread.begin(), patch_num_precomputations_thread.end());
            
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

void Solver::get_coeffs(const long long q, const std::complex<double>* phi,
                std::complex<double>* coeffs)
{

    std::vector<std::complex<double>> evals(Nu_int_*Nv_int_);

    for (int i = 0; i < Nu_int_; i++) {
        for (int j = 0; j < Nv_int_; j++) {

            evals[i*Nv_int_ + j] = dsdtjac_all_[q * Nu_int_*Nv_int_ + i*Nv_int_ + j] * phi[i*Nv_int_ + j];

        }
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
    
    for (int n = 0; n < Nu_int_; n++) {

        double alpha_n;

        if (n == 0) {
            alpha_n = 1;
        } else {
            alpha_n = 2;
        }

        for (int m = 0; m < Nv_int_; m++) {

            double alpha_m;

            if (m == 0) {
                alpha_m = 1;
            } else {
                alpha_m = 2;
            }

            coeffs[n*Nv_int_ + m] *= alpha_n * alpha_m / (Nu_int_*Nv_int_);

        }

    }

}

void Solver::compute_coeffs(const VectorXcd& phi,
                    std::vector<std::complex<double>>& vec_coeffs)
{
    
    #pragma omp parallel for
    for (long long i = 0; i < patch_num_coeffs_.size(); i++) {

        const long long patch_num = patch_num_coeffs_[i];

        get_coeffs(patch_num, &phi[patch_num * Nu_int_*Nv_int_], &vec_coeffs[i * Nu_int_*Nv_int_]);

    }

}

void Solver::int_near(const std::complex<double>* coeffs,
                const std::complex<double>* precomputations,
                std::complex<double>& solution)
{         

    solution = std::complex<double>(0.0, 0.0);

    for (int i = 0; i < Nu_int_; i++) {
        for (int j = 0; j < Nv_int_; j++) {

            solution += precomputations[i*Nv_int_ + j] * coeffs[i*Nv_int_ + j];

        }
    }

}

void Solver::int_far(const double r_0, const double r_1, const double r_2,
                const long long npatch, 
                const std::complex<double>* phi,
                std::complex<double>& solution) 
                
{

    solution = std::complex<double>(0.0, 0.0);

    for (int i = 0; i < Nu_int_; i++) {
        for (int j = 0; j < Nv_int_; j++) {

            const double dsdtjac_loc = dsdtjac_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double constant = dsdtjac_loc * fejer_weights_u_int_[i] * fejer_weights_v_int_[j];

            const double px = disc_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double py = disc_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double pz = disc_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

            const double nx = norm_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double ny = norm_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double nz = norm_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

            std::complex<double> kernel;
            HH2(r_0, r_1, r_2, px, py, pz, nx, ny, nz, coupling_parameter_, wavenumber_, kernel);
            
            solution += constant * kernel * phi[i*Nv_int_ + j];

        }
    }

}

void Solver::compute_integral(const VectorXcd& phi, const std::vector<std::complex<double>>& vec_coeffs,
                        VectorXcd& integral)
{

    #pragma omp parallel for
    for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

        const double r_0 = disc_points_x_all_[point_low_ + npoint];
        const double r_1 = disc_points_y_all_[point_low_ + npoint];
        const double r_2 = disc_points_z_all_[point_low_ + npoint];

        const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);
        const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);

        integral[npoint] = std::complex<double>(0.0, 0.0);

        for (long long q = 0; q < Q_; q++) {             

            for (int q_x = 0; q_x < Qx_; q_x++) {
                for (int q_y = 0; q_y < Qy_; q_y++) {

                    const long long patch_num = q * Qx_*Qy_ + q_x * Qy_ + q_y;

                    const auto patch_begin = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_begin);
                    const auto patch_end = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_end);                            
                    const auto pos = std::lower_bound(patch_begin, patch_end, patch_num);

                    const auto idx = std::distance(patch_num_precomputations_.begin(), pos);
                    
                    std::complex<double> solution(0.0,0.0);

                    if ((point_precomputations_[idx] == (point_low_ + npoint)) && (patch_num_precomputations_[idx] == patch_num)) {

                        const long long start_precomputations = idx * Nu_int_*Nv_int_;
                        
                        const auto pos_coeffs = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
                        const auto idx_coeffs = std::distance(patch_num_coeffs_.begin(), pos_coeffs);

                        const long long start_coeffs = idx_coeffs * Nu_int_*Nv_int_;

                        int_near(&vec_coeffs[start_coeffs], 
                                    &precomputations_[start_precomputations], 
                                    solution);
                        
                    } else {

                        int_far(r_0, r_1, r_2,
                                patch_num, 
                                &phi[patch_num * Nu_int_*Nv_int_],
                                solution);

                    }

                    integral[npoint] += solution;

                }
            }

        }

        if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

            integral[npoint] += 0.5 * phi[point_low_ + npoint];

        }                

    }

}

void Solver::compute_intensities_patch(const long long npatch,
                                const std::complex<double>* phi,
                                std::complex<double>* intensities)
{        

    for (int i = 0; i < Nu_int_; i++) {
        for (int j = 0; j < Nv_int_; j++) {

            intensities[i * Nv_int_ + j] = dsdtjac_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j] * fejer_weights_u_int_[i] * fejer_weights_v_int_[j] * phi[i * Nv_int_ + j];

        }
    }   

}

void Solver::compute_intensities(const VectorXcd& phi,
                            std::vector<std::complex<double>>& intensities)
{

    intensities = std::vector<std::complex<double>>(point_up_-point_low_);

    #pragma omp parallel for
    for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

        compute_intensities_patch(patch_num, 
                                    &phi[patch_num * Nu_int_*Nv_int_],
                                    &intensities[(patch_num - patch_low_) * Nu_int_*Nv_int_]);

    }                    

    std::vector<std::complex<double>> intensities_all(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
    
    MPI_Allgatherv(&intensities[0], point_up_-point_low_, MPI_DOUBLE_COMPLEX, &intensities_all[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE_COMPLEX, mpi_comm_);
    
    intensities = std::vector<std::complex<double>>(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
    
    #pragma omp parallel for
    for (long long i = 0; i < Q_ * Qx_*Qy_ * Nu_int_*Nv_int_; i++) {

        intensities[i] = intensities_all[new_order_points_IFGF_[i]];

    }

}

void Solver::compute_integral_acc(const VectorXcd& phi, const std::vector<std::complex<double>>& vec_coeffs,
                            VectorXcd& integral)
{

    integral.setZero();

    #pragma omp parallel for
    for (long long i = 0; i < point_precomputations_.size(); i++) {

        const long long npoint = point_precomputations_[i];
        const long long patch_num = patch_num_precomputations_[i];
        
        const long long start_precomputations = i * Nu_int_*Nv_int_;

        const auto pos_coeffs = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
        const auto idx_coeffs = std::distance(patch_num_coeffs_.begin(), pos_coeffs);

        const long long start_coeffs = idx_coeffs * Nu_int_*Nv_int_;

        std::complex<double> solution;

        int_near(&vec_coeffs[start_coeffs], 
                    &precomputations_[start_precomputations], 
                    solution);

        integral[npoint - point_low_] += solution;

    }

    if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

        #pragma omp parallel for
        for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

            integral[npoint] += 0.5 * phi[point_low_ + npoint];

        }

    }

} 

void Solver::create_fftw_objects()
{

    if (METHOD_FFTW == 1) {

        fftw_init_threads();

        fftw_plan_with_nthreads(NTHREADS);

        int rank = 2;

        int n_orig[2] = {Nu_int_, Nv_int_};
        int n_overs[2] = {Nu_int_*P_over_, Nv_int_*P_over_};

        original_patch_ = (double*) fftw_malloc(sizeof(double) * n_orig[0]*n_orig[1] * patch_num_coeffs_.size() * 2);
        oversampled_patch_ = (double*) fftw_malloc(sizeof(double) * n_overs[0]*n_overs[1] * patch_num_coeffs_.size() * 2);   
            
        int istride_orig = 1;
        int idist_orig = n_orig[0]*n_orig[1];
        
        int istride_overs = 1;
        int idist_overs = n_overs[0]*n_overs[1];

        fftw_r2r_kind kind_orig[2] = {FFTW_REDFT10, FFTW_REDFT10};
        fftw_r2r_kind kind_overs[2] = {FFTW_REDFT01, FFTW_REDFT01};

        if (TYPE_FFTW == 1) {

            plan_original_patch_ = fftw_plan_many_r2r(rank, n_orig, patch_num_coeffs_.size() * 2, original_patch_, n_orig, istride_orig, idist_orig, original_patch_, n_orig, istride_orig, idist_orig, kind_orig, FFTW_PATIENT);
            plan_oversampled_patch_ = fftw_plan_many_r2r(rank, n_overs, patch_num_coeffs_.size() * 2, oversampled_patch_, n_overs, istride_overs, idist_overs, oversampled_patch_, n_overs, istride_overs, idist_overs, kind_overs, FFTW_PATIENT);

        } else {

            fftw_iodim dims_orig[rank];
            dims_orig[1].n = n_orig[1];
            dims_orig[1].is = istride_orig;
            dims_orig[1].os = istride_orig;
            dims_orig[0].n = n_orig[0];
            dims_orig[0].is = n_orig[1] * dims_orig[1].is;
            dims_orig[0].os = n_orig[1] * dims_orig[1].os;

            fftw_iodim dims_overs[rank];
            dims_overs[1].n = n_overs[1];
            dims_overs[1].is = istride_overs;
            dims_overs[1].os = istride_overs;
            dims_overs[0].n = n_overs[0];
            dims_overs[0].is = n_overs[1] * dims_overs[1].is;
            dims_overs[0].os = n_overs[1] * dims_overs[1].os;

            int howmany_rank = 1;

            fftw_iodim howmany_dims_orig[howmany_rank]; 
            howmany_dims_orig[0].n = patch_num_coeffs_.size() * 2;
            howmany_dims_orig[0].is = idist_orig;
            howmany_dims_orig[0].os = idist_orig;

            fftw_iodim howmany_dims_overs[howmany_rank]; 
            howmany_dims_overs[0].n = patch_num_coeffs_.size() * 2;
            howmany_dims_overs[0].is = idist_overs;
            howmany_dims_overs[0].os = idist_overs;

            plan_original_patch_ = fftw_plan_guru_r2r(rank, dims_orig, howmany_rank, howmany_dims_orig, original_patch_, original_patch_, kind_orig, FFTW_PATIENT);
            plan_oversampled_patch_ = fftw_plan_guru_r2r(rank, dims_overs, howmany_rank, howmany_dims_overs, oversampled_patch_, oversampled_patch_, kind_overs, FFTW_PATIENT);

        }

    } else {

        int rank = 2;

        int n_orig[2] = {Nu_int_, Nv_int_};
        int n_overs[2] = {Nu_int_*P_over_, Nv_int_*P_over_};
        
        original_patch_ = (double*) fftw_malloc(sizeof(double) * n_orig[0]*n_orig[1] * 2);
        oversampled_patch_ = (double*) fftw_malloc(sizeof(double) * n_overs[0]*n_overs[1] * 2);   

        int istride_orig = 1;
        int idist_orig = n_orig[0]*n_orig[1];
        
        int istride_overs = 1;
        int idist_overs = n_overs[0]*n_overs[1];

        fftw_r2r_kind kind_orig[2] = {FFTW_REDFT10, FFTW_REDFT10};
        fftw_r2r_kind kind_overs[2] = {FFTW_REDFT01, FFTW_REDFT01};

        if (TYPE_FFTW == 1) {

            plan_original_patch_ = fftw_plan_many_r2r(rank, n_orig, 2, original_patch_, n_orig, istride_orig, idist_orig, original_patch_, n_orig, istride_orig, idist_orig, kind_orig, FFTW_PATIENT);
            plan_oversampled_patch_ = fftw_plan_many_r2r(rank, n_overs, 2, oversampled_patch_, n_overs, istride_overs, idist_overs, oversampled_patch_, n_overs, istride_overs, idist_overs, kind_overs, FFTW_PATIENT);

        } else {

            fftw_iodim dims_orig[rank];
            dims_orig[1].n = n_orig[1];
            dims_orig[1].is = istride_orig;
            dims_orig[1].os = istride_orig;
            dims_orig[0].n = n_orig[0];
            dims_orig[0].is = n_orig[1] * dims_orig[1].is;
            dims_orig[0].os = n_orig[1] * dims_orig[1].os;

            fftw_iodim dims_overs[rank];
            dims_overs[1].n = n_overs[1];
            dims_overs[1].is = istride_overs;
            dims_overs[1].os = istride_overs;
            dims_overs[0].n = n_overs[0];
            dims_overs[0].is = n_overs[1] * dims_overs[1].is;
            dims_overs[0].os = n_overs[1] * dims_overs[1].os;

            int howmany_rank = 1;

            fftw_iodim howmany_dims_orig[howmany_rank]; 
            howmany_dims_orig[0].n = 2;
            howmany_dims_orig[0].is = idist_orig;
            howmany_dims_orig[0].os = idist_orig;

            fftw_iodim howmany_dims_overs[howmany_rank]; 
            howmany_dims_overs[0].n = 2;
            howmany_dims_overs[0].is = idist_overs;
            howmany_dims_overs[0].os = idist_overs;

            plan_original_patch_ = fftw_plan_guru_r2r(rank, dims_orig, howmany_rank, howmany_dims_orig, original_patch_, original_patch_, kind_orig, FFTW_PATIENT);
            plan_oversampled_patch_ = fftw_plan_guru_r2r(rank, dims_overs, howmany_rank, howmany_dims_overs, oversampled_patch_, oversampled_patch_, kind_overs, FFTW_PATIENT);

        }

    }

}

void Solver::compute_oversampling_M1(const std::vector<double>& psi, std::vector<double>& psi_overs) 
{

    #pragma omp parallel for
    for (long long i = 0; i < psi.size(); i++) {

        original_patch_[i] = psi[i] / (Nu_int_ * Nv_int_);

    }

    fftw_execute(plan_original_patch_);

    std::memset(oversampled_patch_, 0.0, patch_num_coeffs_.size() * Nu_int_*Nv_int_ * P_over_*P_over_ * 2 * sizeof(double));

    #pragma omp parallel for
    for (long long k = 0; k < patch_num_coeffs_.size(); k++) {

        for (int i = 0; i < Nu_int_; i++) {
            for (int j = 0; j < Nv_int_; j++) {

                oversampled_patch_[k * 2 * Nu_int_*Nv_int_ * P_over_*P_over_ + i * Nv_int_ * P_over_ + j] = original_patch_[k * 2 * Nu_int_*Nv_int_ + i * Nv_int_ + j] * 0.5 * 0.5;
                oversampled_patch_[k * 2 * Nu_int_*Nv_int_ * P_over_*P_over_ + Nu_int_*Nv_int_ * P_over_*P_over_ + i * Nv_int_ * P_over_ + j] = original_patch_[k * 2 * Nu_int_*Nv_int_ + Nu_int_*Nv_int_ + i * Nv_int_ + j] * 0.5 * 0.5;

            }
        }

    }

    fftw_execute(plan_oversampled_patch_);

    psi_overs = std::vector<double>(oversampled_patch_, oversampled_patch_ + 2 * patch_num_coeffs_.size() * Nu_int_*Nv_int_ * P_over_*P_over_);

}

void Solver::compute_oversampling_M2(const std::vector<double>& psi, std::vector<double>& psi_overs) 
{

    for (long long i = 0; i < psi.size(); i++) {

        original_patch_[i] = psi[i] / (Nu_int_ * Nv_int_);

    }

    fftw_execute(plan_original_patch_);

    std::memset(oversampled_patch_, 0.0, Nu_int_*Nv_int_ * P_over_*P_over_ * 2 * sizeof(double));

    for (int i = 0; i < Nu_int_; i++) {
        for (int j = 0; j < Nv_int_; j++) {

            oversampled_patch_[i * Nv_int_ * P_over_ + j] = original_patch_[i * Nv_int_ + j] * 0.5 * 0.5;
            oversampled_patch_[Nu_int_*Nv_int_ * P_over_*P_over_ + i * Nv_int_ * P_over_ + j] = original_patch_[Nu_int_*Nv_int_ + i * Nv_int_ + j] * 0.5 * 0.5;

        }
    }

    fftw_execute(plan_oversampled_patch_);

    psi_overs = std::vector<double>(oversampled_patch_, oversampled_patch_ + 2 * Nu_int_*Nv_int_ * P_over_*P_over_);

}

void Solver::compute_singular_points()
{

    #pragma omp parallel
    {

    std::vector<long long> point_precomputations_thread;
    std::vector<long long> patch_num_precomputations_thread;
    std::vector<double> argmin_precomputations_u_thread;
    std::vector<double> argmin_precomputations_v_thread;

    #pragma omp for      
    for (long long npoint = point_low_; npoint < point_up_; npoint++) {

        const double r_0 = disc_points_x_all_[npoint];
        const double r_1 = disc_points_y_all_[npoint];
        const double r_2 = disc_points_z_all_[npoint];

        const long long patch_num_own = npoint / (Nu_int_*Nv_int_);

        const long long start = start_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];
        const long long total = size_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];

        for (long long i = 0; i < total; i++) {

            const long long patch_num = sing_and_near_sing_patches_estimate_[start + i];

            double min_dist_loc;
            int indx_point_min_dist_loc;

            bool sing_near_sing_flag = false;

            if (patch_num_own == patch_num) {

                min_dist_loc = 0.0;
                indx_point_min_dist_loc = npoint % (Nu_int_*Nv_int_);   

                sing_near_sing_flag = true; 

            } else {  

                min_dist_loc = std::numeric_limits<double>::max();

                for (int ii = 0; ii < Nu_int_; ii++) {
                    for (int jj = 0; jj < Nv_int_; jj++) {

                        const long long position = patch_num * Nu_int_*Nv_int_ + ii * Nv_int_ + jj;

                        const double diff_0 = disc_points_x_all_[position] - r_0;
                        const double diff_1 = disc_points_y_all_[position] - r_1;
                        const double diff_2 = disc_points_z_all_[position] - r_2;

                        const double dist_aux = diff_0*diff_0 + diff_1*diff_1 + diff_2*diff_2;

                        if (dist_aux < min_dist_loc) {

                            min_dist_loc = dist_aux;
                            indx_point_min_dist_loc = ii * Nv_int_ + jj;

                        }

                        const bool proximity = (std::abs(std::floor(r_0/proximity_) - std::floor(disc_points_x_all_[position]/proximity_)) <= 1) &&
                                                (std::abs(std::floor(r_1/proximity_) - std::floor(disc_points_y_all_[position]/proximity_)) <= 1) &&
                                                (std::abs(std::floor(r_2/proximity_) - std::floor(disc_points_z_all_[position]/proximity_)) <= 1);

                        if (!sing_near_sing_flag && proximity) {

                            sing_near_sing_flag = true;

                        }

                    }
                }

                min_dist_loc = std::sqrt(min_dist_loc);

            }
        
            if (sing_near_sing_flag) {   

                const int nu_loc = indx_point_min_dist_loc / Nv_int_;
                const int nv_loc = indx_point_min_dist_loc % Nv_int_;  

                const int flag_u_loc = flags_domain_u_all_[patch_num];
                const int flag_v_loc = flags_domain_v_all_[patch_num];

                const long long q = patch_num / (Qx_ * Qy_);
                const int q_x = (patch_num / Qy_) % Qx_;
                const int q_y = patch_num % Qy_;

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
                    
                        auto dist_func = [r_0, r_1, r_2, q, u_a_loc, u_b_loc, v_a_loc, v_b_loc, this](double s, double t) -> double {const double ss = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, s);
                                                                                                                                const double tt = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, t);
                                                                                                                                double rp_0, rp_1, rp_2;
                                                                                                                                parametrization_q(this->SPHERE_RADIUS,this->SPHERE_CENTER,ss, tt, q, rp_0, rp_1, rp_2);
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

                point_precomputations_thread.push_back(npoint);
                patch_num_precomputations_thread.push_back(patch_num);
                argmin_precomputations_u_thread.push_back(ubar_loc);
                argmin_precomputations_v_thread.push_back(vbar_loc);


            }

        }

    }

    #pragma omp for ordered
    for (int i = 0; i < NTHREADS; i++) {

        #pragma omp ordered
        {

            point_precomputations_.insert(point_precomputations_.end(), point_precomputations_thread.begin(), point_precomputations_thread.end());
            patch_num_precomputations_.insert(patch_num_precomputations_.end(), patch_num_precomputations_thread.begin(), patch_num_precomputations_thread.end());
            argmin_precomputations_u_.insert(argmin_precomputations_u_.end(), argmin_precomputations_u_thread.begin(), argmin_precomputations_u_thread.end());
            argmin_precomputations_v_.insert(argmin_precomputations_v_.end(), argmin_precomputations_v_thread.begin(), argmin_precomputations_v_thread.end());

            std::vector<long long>().swap(point_precomputations_thread);
            std::vector<long long>().swap(patch_num_precomputations_thread);
            std::vector<double>().swap(argmin_precomputations_u_thread);
            std::vector<double>().swap(argmin_precomputations_v_thread);

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

void Solver::int_near_overs(const double r_0, const double r_1, const double r_2,
                    const long long q, const int flag_u_loc, const int flag_v_loc, 
                    const double u_a_loc, const double u_b_loc, const double v_a_loc, const double v_b_loc,
                    const double ubar_loc, const double vbar_loc,
                    const double* psi_overs,
                    std::complex<double>& solution) 
{

    solution = std::complex<double>(0.0, 0.0);

    std::vector<double> s(Nu_int_*P_over_ + 2);
    std::vector<double> t(Nv_int_*P_over_ + 2);

    s[0] = u_a_loc;

    for (int i = 0; i < Nu_int_*P_over_; i++) {

        const double xi = std::cos(M_PI * (2.0 * (Nu_int_*P_over_ - 1 - i) + 1.0) / (2.0 * Nu_int_*P_over_));
        s[i+1] = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, eta(xi, flag_u_loc));

    }          

    s[Nu_int_*P_over_ + 1] = u_b_loc;

    t[0] = v_a_loc;
    
    for (int i = 0; i < Nv_int_*P_over_; i++) {

        const double xi = std::cos(M_PI * (2.0 * (Nv_int_*P_over_ - 1 - i) + 1.0) / (2.0 * Nv_int_*P_over_));
        t[i+1] = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, eta(xi, flag_v_loc));

    }

    t[Nv_int_*P_over_ + 1] = v_b_loc;                

    std::vector<double> ui(Nu_prec_), si(Nu_prec_), muwi(Nu_prec_);
    std::vector<double> vj(Nv_prec_), tj(Nv_prec_), muwj(Nv_prec_);

    for (int i = 0; i < Nu_prec_; i++) {

        ui[i] = xi(ubar_loc, fejer_nodes_u_prec_[i]);
        si[i] = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, eta(ui[i], flag_u_loc));
        muwi[i] = der_xi(ubar_loc, fejer_nodes_u_prec_[i]) * fejer_weights_u_prec_[i];
    
    }

    for (int j = 0; j < Nv_prec_; j++) {

        vj[j] = xi(vbar_loc, fejer_nodes_v_prec_[j]);
        tj[j] = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, eta(vj[j], flag_v_loc));
        muwj[j] = der_xi(vbar_loc, fejer_nodes_v_prec_[j]) * fejer_weights_v_prec_[j];
    
    }

    std::array<double, N_LOC_X> s_loc;
    std::array<double, N_LOC_Y> t_loc;
    std::array<double, N_LOC_X*N_LOC_Y> psi_overs_loc_real, psi_overs_loc_imag;

    std::vector<std::complex<double>> psi_overs_loc(Nu_prec_*Nv_prec_);
    
    for (int i = 0; i < Nu_prec_; i++) {

        const double ss = si[i];
        auto idx_s = std::lower_bound(s.begin(), s.end(), ss);
        int idx_ss = idx_s - s.begin();

        if (idx_ss == 0) {
            idx_ss = 1;
        } else if (idx_ss == Nu_int_*P_over_ + 1) {
            idx_ss = Nu_int_*P_over_;
        }

        if (idx_ss - (N_LOC_X - 1) / 2 < 1) {
            idx_ss = (N_LOC_X - 1) / 2 + 1;
        } else if (idx_ss + (N_LOC_X - 1) / 2 > Nu_int_*P_over_) {
            idx_ss = Nu_int_*P_over_ - (N_LOC_X - 1) / 2;
        }

        for (int k = - (N_LOC_X - 1) / 2; k <= (N_LOC_X - 1) / 2; k++) {
            s_loc[k + (N_LOC_X - 1) / 2] = s[idx_ss + k];
        }

        for (int j = 0; j < Nv_prec_; j++) {

            const double tt = tj[j];
            auto idx_t = std::lower_bound(t.begin(), t.end(), tt);
            int idx_tt = idx_t - t.begin();

            if (idx_tt == 0) {
                idx_tt = 1;
            } else if (idx_tt == Nv_int_*P_over_ + 1) {
                idx_tt = Nv_int_*P_over_;
            }

            if (idx_tt - (N_LOC_Y - 1) / 2 < 1) {
                idx_tt = (N_LOC_Y - 1) / 2 + 1;
            } else if (idx_tt + (N_LOC_Y - 1) / 2 > Nv_int_*P_over_) {
                idx_tt = Nv_int_*P_over_ - (N_LOC_Y - 1) / 2;
            }

            for (int l = - (N_LOC_Y - 1) / 2; l <= (N_LOC_Y - 1) / 2; l++) {
                t_loc[l + (N_LOC_Y - 1) / 2] = t[idx_tt + l];
            }

            for (int k = - (N_LOC_X - 1) / 2; k <= (N_LOC_X - 1) / 2; k++) {
                for (int l = - (N_LOC_Y - 1) / 2; l <= (N_LOC_Y - 1) / 2; l++) {

                    const double real_val = psi_overs[(Nu_int_*P_over_ - idx_ss - k) * Nv_int_*P_over_ + (Nv_int_*P_over_ - idx_tt - l)];
                    const double imag_val = psi_overs[Nu_int_*Nv_int_ * P_over_*P_over_ + (Nu_int_*P_over_ - idx_ss - k) * Nv_int_*P_over_ + (Nv_int_*P_over_ - idx_tt - l)];

                    psi_overs_loc_real[(k + (N_LOC_X - 1) / 2) * N_LOC_Y + (l + (N_LOC_Y - 1) / 2)] = real_val;
                    psi_overs_loc_imag[(k + (N_LOC_X - 1) / 2) * N_LOC_Y + (l + (N_LOC_Y - 1) / 2)] = imag_val;

                }
            } 

            const double interp_val_real = LocalInterpNewton2D<double, N_LOC_X, N_LOC_Y>(s_loc, t_loc, psi_overs_loc_real, ss, tt);
            const double interp_val_imag = LocalInterpNewton2D<double, N_LOC_X, N_LOC_Y>(s_loc, t_loc, psi_overs_loc_imag, ss, tt);
            const std::complex<double> interp_val(interp_val_real, interp_val_imag);

            psi_overs_loc[i * Nv_prec_ + j] = interp_val;

        }

    }

    if (GEOMETRY == 0) {

        for (int i = 0; i < Nu_prec_; i++) {
            for (int j = 0; j < Nv_prec_; j++) {       

                const double ss = si[i];
                const double tt = tj[j];                 

                double px, py, pz;
                double nx, ny, nz;
                double dxdsx, dxdsy, dxdsz;
                double dxdtx, dxdty, dxdtz;
                double dxdsdsx, dxdsdsy, dxdsdsz;
                double dxdsdtx, dxdsdty, dxdsdtz;
                double dxdtdtx, dxdtdty, dxdtdtz;

                parametrization_q(SPHERE_RADIUS, SPHERE_CENTER, ss, tt, q, px, py, pz);
                normal_q(SPHERE_RADIUS,ss, tt, q, nx, ny, nz);
                dxds_q(SPHERE_RADIUS,ss, tt, q, dxdsx, dxdsy, dxdsz);
                dxdt_q(SPHERE_RADIUS,ss, tt, q, dxdtx, dxdty, dxdtz);
                dxdsds_q(SPHERE_RADIUS,ss, tt, q, dxdsdsx, dxdsdsy, dxdsdsz);
                dxdsdt_q(SPHERE_RADIUS,ss, tt, q, dxdsdtx, dxdsdty, dxdsdtz);
                dxdtdt_q(SPHERE_RADIUS,ss, tt, q, dxdtdtx, dxdtdty, dxdtdtz);

                dxdsx *= 0.5 * (u_b_loc - u_a_loc);
                dxdsy *= 0.5 * (u_b_loc - u_a_loc);
                dxdsz *= 0.5 * (u_b_loc - u_a_loc);

                dxdtx *= 0.5 * (v_b_loc - v_a_loc);
                dxdty *= 0.5 * (v_b_loc - v_a_loc);
                dxdtz *= 0.5 * (v_b_loc - v_a_loc);

                dxdsdsx *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsdsy *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsdsz *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (u_b_loc - u_a_loc);

                dxdsdtx *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdsdty *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdsdtz *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (v_b_loc - v_a_loc);

                dxdtdtx *= 0.5 * (v_b_loc - v_a_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdtdty *= 0.5 * (v_b_loc - v_a_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdtdtz *= 0.5 * (v_b_loc - v_a_loc) * 0.5 * (v_b_loc - v_a_loc);

                std::complex<double> kernel;
                HH(r_0, r_1, r_2,
                px, py, pz,
                nx, ny, nz,
                dxdsx, dxdsy, dxdsz,
                dxdtx, dxdty, dxdtz,
                dxdsdsx, dxdsdsy, dxdsdsz,
                dxdsdtx, dxdsdty, dxdsdtz,
                dxdtdtx, dxdtdty, dxdtdtz,
                coupling_parameter_, wavenumber_,
                kernel);

                solution += kernel * muwi[i] * muwj[j] * psi_overs_loc[i * Nv_prec_ + j];

            }
        }

    } else {

        InterpPatch elem = interp_surface_[q];

        std::vector<double> px(Nu_prec_*Nv_prec_), py(Nu_prec_*Nv_prec_), pz(Nu_prec_*Nv_prec_);
        std::vector<double> nx(Nu_prec_*Nv_prec_), ny(Nu_prec_*Nv_prec_), nz(Nu_prec_*Nv_prec_);
        std::vector<double> dxdsx(Nu_prec_*Nv_prec_), dxdsy(Nu_prec_*Nv_prec_), dxdsz(Nu_prec_*Nv_prec_);
        std::vector<double> dxdtx(Nu_prec_*Nv_prec_), dxdty(Nu_prec_*Nv_prec_), dxdtz(Nu_prec_*Nv_prec_);

        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.x, si, tj, &px[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.y, si, tj, &py[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.z, si, tj, &pz[0]);

        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.nuX, si, tj, &nx[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.nuY, si, tj, &ny[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.nuZ, si, tj, &nz[0]);

        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dxdu, si, tj, &dxdsx[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dydu, si, tj, &dxdsy[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dzdu, si, tj, &dxdsz[0]);

        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dxdv, si, tj, &dxdtx[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dydv, si, tj, &dxdty[0]);
        lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
                                    elem.dzdv, si, tj, &dxdtz[0]);

        for (int i = 0; i < Nu_prec_; i++) {
            for (int j = 0; j < Nv_prec_; j++) {

                std::complex<double> kernel;

                HH(r_0, r_1, r_2,
                px[i*Nv_prec_+j], py[i*Nv_prec_+j], pz[i*Nv_prec_+j],
                nx[i*Nv_prec_+j], ny[i*Nv_prec_+j], nz[i*Nv_prec_+j],
                dxdsx[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc), dxdsy[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc), dxdsz[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc),
                dxdtx[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc), dxdty[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc), dxdtz[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc),
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                coupling_parameter_, wavenumber_,
                kernel);

                solution += kernel * muwi[i] * muwj[j] * psi_overs_loc[i * Nv_prec_ + j];

            }
        }

    }

}

void Solver::write_psi_M1(const VectorXcd& phi, std::vector<double>& psi)
{

    #pragma omp parallel for
    for (long long i = 0; i < patch_num_coeffs_.size(); i++) {

        const long long patch_num = patch_num_coeffs_[i];                      

        for (int nu = 0; nu < Nu_int_; nu++) {
            for (int nv = 0; nv < Nv_int_; nv++) {

                const long long position = patch_num * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;
                const long long position_real = i * 2 * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;
                const long long position_imag = i * 2 * Nu_int_*Nv_int_ + Nu_int_*Nv_int_ + nu * Nv_int_ + nv;

                psi[position_real] = phi[position].real() * dsdtjac_all_[position];
                psi[position_imag] = phi[position].imag() * dsdtjac_all_[position];

            }
        }

    }

}

void Solver::write_psi_M2(const long long patch_num, const std::complex<double>* phi, std::vector<double>& psi)
{

    for (int nu = 0; nu < Nu_int_; nu++) {
        for (int nv = 0; nv < Nv_int_; nv++) {

            const int position = nu * Nv_int_ + nv;
            const int position_real = nu * Nv_int_ + nv;
            const int position_imag = Nu_int_*Nv_int_ + nu * Nv_int_ + nv;

            psi[position_real] = phi[position].real() * dsdtjac_all_[patch_num * Nu_int_*Nv_int_ + position];
            psi[position_imag] = phi[position].imag() * dsdtjac_all_[patch_num * Nu_int_*Nv_int_ + position];

        }
    }

}

void Solver::compute_integral_overs_M1(const VectorXcd& phi,
                                const std::vector<double>& psi_overs,
                                VectorXcd& integral)
{
    
    #pragma omp parallel for
    for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

        const double r_0 = disc_points_x_all_[point_low_ + npoint];
        const double r_1 = disc_points_y_all_[point_low_ + npoint];
        const double r_2 = disc_points_z_all_[point_low_ + npoint];

        const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);
        const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);

        integral[npoint] = std::complex<double>(0.0, 0.0);

        for (long long q = 0; q < Q_; q++) {             

            for (int q_x = 0; q_x < Qx_; q_x++) {
                for (int q_y = 0; q_y < Qy_; q_y++) {

                    const long long patch_num = q * Qx_*Qy_ + q_x * Qy_ + q_y;

                    const auto patch_begin = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_begin);
                    const auto patch_end = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_end);                            
                    const auto pos = std::lower_bound(patch_begin, patch_end, patch_num);

                    const auto idx = std::distance(patch_num_precomputations_.begin(), pos);
                    
                    std::complex<double> solution;

                    if ((point_precomputations_[idx] == (point_low_ + npoint)) && (patch_num_precomputations_[idx] == patch_num)) {

                        const long long flag_u_loc = flags_domain_u_all_[patch_num];
                        const long long flag_v_loc = flags_domain_v_all_[patch_num];

                        const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
                        const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
                        const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
                        const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

                        const double ubar_loc = argmin_precomputations_u_[idx];
                        const double vbar_loc = argmin_precomputations_v_[idx];

                        const auto patch_begin_2 = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
                        const auto idx_2 = std::distance(patch_num_coeffs_.begin(), patch_begin_2);
                        
                        int_near_overs(r_0, r_1, r_2,
                                        q, flag_u_loc, flag_v_loc,
                                        u_a_loc, u_b_loc, v_a_loc, v_b_loc,
                                        ubar_loc, vbar_loc,
                                        &psi_overs[idx_2 * 2 * Nu_int_*Nv_int_ * P_over_*P_over_],
                                        solution);

                    } else {                            

                        int_far(r_0, r_1, r_2,
                                patch_num, 
                                &phi[patch_num * Nu_int_*Nv_int_],
                                solution);
                        
                    }

                    integral[npoint] += solution;

                }
            }

        }

        if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

            integral[npoint] += 0.5 * phi[point_low_ + npoint];

        }

    }

}

void Solver::compute_integral_overs_M2(const VectorXcd& phi,
                                VectorXcd& integral)
{

    integral.setZero();   

    #pragma omp parallel for
    for (long long q = 0; q < Q_; q++) {

        for (int q_x = 0; q_x < Qx_; q_x++) {
            for (int q_y = 0; q_y < Qy_; q_y++) {

                const long long patch_num = q * Qx_*Qy_ + q_x * Qy_ + q_y;

                const int flag_u_loc = flags_domain_u_all_[patch_num];
                const int flag_v_loc = flags_domain_v_all_[patch_num];

                const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
                const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
                const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
                const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

                const auto patch_begin_2 = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
                const auto idx_2 = std::distance(patch_num_coeffs_.begin(), patch_begin_2);
                
                std::vector<double> psi(2 * Nu_int_*Nv_int_);            
                std::vector<double> psi_overs;
                
                #pragma omp critical
                {

                if (patch_num_coeffs_[idx_2] == patch_num) {

                    write_psi_M2(patch_num, &phi[patch_num * Nu_int_*Nv_int_], psi);
                    
                    compute_oversampling_M2(psi, psi_overs);

                }  

                }      
                
                for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

                    const double r_0 = disc_points_x_all_[point_low_ + npoint];
                    const double r_1 = disc_points_y_all_[point_low_ + npoint];
                    const double r_2 = disc_points_z_all_[point_low_ + npoint];

                    const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);
                    const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);

                    const auto patch_begin = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_begin);
                    const auto patch_end = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_end);                            
                    const auto pos = std::lower_bound(patch_begin, patch_end, patch_num);

                    const auto idx = std::distance(patch_num_precomputations_.begin(), pos);
                    
                    std::complex<double> solution;

                    if ((point_precomputations_[idx] == (point_low_ + npoint)) && (patch_num_precomputations_[idx] == patch_num)) {

                        const double ubar_loc = argmin_precomputations_u_[idx];
                        const double vbar_loc = argmin_precomputations_v_[idx];
                        
                        int_near_overs(r_0, r_1, r_2,
                                        q, flag_u_loc, flag_v_loc,
                                        u_a_loc, u_b_loc, v_a_loc, v_b_loc,
                                        ubar_loc, vbar_loc,
                                        &psi_overs[0],
                                        solution);

                    } else {                            

                        int_far(r_0, r_1, r_2,
                                patch_num, 
                                &phi[patch_num * Nu_int_*Nv_int_],
                                solution);
                        
                    }

                    integral[npoint] += solution;

                }

            }
        }

    } 

    if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

        #pragma omp parallel for
        for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

            integral[npoint] += 0.5 * phi[point_low_ + npoint];

        }

    }  

}

void Solver::compute_integral_overs_acc_M1(const VectorXcd& phi,
                                    const std::vector<double>& psi_overs,
                                    VectorXcd& integral)
{          

    integral.setZero();

    #pragma omp parallel for
    for (long long i = 0; i < point_precomputations_.size(); i++) {

        const long long npoint = point_precomputations_[i];
        const long long patch_num = patch_num_precomputations_[i];

        const double r_0 = disc_points_x_all_[npoint];
        const double r_1 = disc_points_y_all_[npoint];
        const double r_2 = disc_points_z_all_[npoint];

        const int flag_u_loc = flags_domain_u_all_[patch_num];
        const int flag_v_loc = flags_domain_v_all_[patch_num];

        const long long q = patch_num / (Qx_*Qy_);
        const int q_x = (patch_num % (Qx_*Qy_)) / Qy_;
        const int q_y = (patch_num % (Qx_*Qy_)) % Qy_;

        const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
        const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
        const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
        const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

        const double ubar_loc = argmin_precomputations_u_[i];
        const double vbar_loc = argmin_precomputations_v_[i];

        const auto patch_begin = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
        const auto idx_2 = std::distance(patch_num_coeffs_.begin(), patch_begin);

        std::complex<double> solution;

        int_near_overs(r_0, r_1, r_2,
                        q, flag_u_loc, flag_v_loc,
                        u_a_loc, u_b_loc, v_a_loc, v_b_loc,
                        ubar_loc, vbar_loc,
                        &psi_overs[idx_2 * 2 * Nu_int_*Nv_int_ * P_over_*P_over_],
                        solution);

        integral[npoint - point_low_] += solution;

    }

    if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

        #pragma omp parallel for
        for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

            integral[npoint] += 0.5 * phi[point_low_ + npoint];

        }

    }        

}

void Solver::compute_integral_overs_acc_M2(const VectorXcd& phi,
                                    VectorXcd& integral)
{

    integral.setZero();           

    #pragma omp parallel for 
    for (long long patch_num : patch_num_coeffs_) {

        const long long q = patch_num / (Qx_*Qy_);
        const int q_x = (patch_num % (Qx_*Qy_)) / Qy_;
        const int q_y = (patch_num % (Qx_*Qy_)) % Qy_;

        const int flag_u_loc = flags_domain_u_all_[patch_num];
        const int flag_v_loc = flags_domain_v_all_[patch_num];

        const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
        const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
        const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
        const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

        std::vector<double> psi(2 * Nu_int_*Nv_int_);            
        std::vector<double> psi_overs;

        #pragma omp critical
        {

        write_psi_M2(patch_num, &phi[patch_num * Nu_int_*Nv_int_], psi);

        compute_oversampling_M2(psi, psi_overs);

        }

        std::vector<long long> points; // = patch_num_FFTW_[patch_num];

        for (long long npoint : points) {

            const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), npoint);
            const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), npoint);
            const auto patch_begin = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_begin);
            const auto patch_end = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_end);                            
            const auto pos = std::lower_bound(patch_begin, patch_end, patch_num);
            const auto idx = std::distance(patch_num_precomputations_.begin(), pos);
                    
            std::complex<double> solution;

            const double r_0 = disc_points_x_all_[npoint];
            const double r_1 = disc_points_y_all_[npoint];
            const double r_2 = disc_points_z_all_[npoint];

            const double ubar_loc = argmin_precomputations_u_[idx];
            const double vbar_loc = argmin_precomputations_v_[idx];

            int_near_overs(r_0, r_1, r_2,
                            q, flag_u_loc, flag_v_loc,
                            u_a_loc, u_b_loc, v_a_loc, v_b_loc,
                            ubar_loc, vbar_loc,
                            &psi_overs[0],
                            solution);

            integral[npoint - point_low_] += solution;

        }

    }

    if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

        #pragma omp parallel for
        for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

            integral[npoint] += 0.5 * phi[point_low_ + npoint];

        }

    }   

}