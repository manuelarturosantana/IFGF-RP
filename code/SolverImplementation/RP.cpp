
#include "../solver2.h"

void Solver::beta(const double r_0, const double r_1, const double r_2,
                  const long long q, const int flag_u_loc, const int flag_v_loc, 
                  const double u_a_loc, const double u_b_loc, const double v_a_loc, const double v_b_loc,
                  const double ubar_loc, const double vbar_loc,
                  std::vector<std::complex<double>>& prec)
{

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
    
    std::vector<std::complex<double>> H(Nu_prec_*Nv_prec_);
    
    if (GEOMETRY == 0) {

        for (int i = 0; i < Nu_prec_; i++) {
            for (int j = 0; j < Nv_prec_; j++) {                        
                
                const double s = si[i];
                const double t = tj[j];

                double px, py, pz;
                double nx, ny, nz;
                double dxdsx, dxdsy, dxdsz;
                double dxdtx, dxdty, dxdtz;
                double dxdsdsx, dxdsdsy, dxdsdsz;
                double dxdsdtx, dxdsdty, dxdsdtz;
                double dxdtdtx, dxdtdty, dxdtdtz;

                parametrization_q(s, t, q, px, py, pz);
                normal_q(s, t, q, nx, ny, nz);
                dxds_q(s, t, q, dxdsx, dxdsy, dxdsz);
                dxdt_q(s, t, q, dxdtx, dxdty, dxdtz);
                dxdsds_q(s, t, q, dxdsdsx, dxdsdsy, dxdsdsz);
                dxdsdt_q(s, t, q, dxdsdtx, dxdsdty, dxdsdtz);
                dxdtdt_q(s, t, q, dxdtdtx, dxdtdty, dxdtdtz);                        

                dxdsx *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsy *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsz *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);

                dxdtx *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdty *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdtz *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);

                dxdsdsx *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsdsy *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);
                dxdsdsz *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc);

                dxdsdtx *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdsdty *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdsdtz *= der_eta(s, flag_u_loc) * 0.5 * (u_b_loc - u_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);

                dxdtdtx *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdtdty *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);
                dxdtdtz *= der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc) * der_eta(t, flag_v_loc) * 0.5 * (v_b_loc - v_a_loc);

                HH(r_0, r_1, r_2,
                px, py, pz,
                nx, ny, nz,
                dxdsx, dxdsy, dxdsz,
                dxdtx, dxdty, dxdtz,
                dxdsdsx, dxdsdsy, dxdsdsz,
                dxdsdtx, dxdsdty, dxdsdtz,
                dxdtdtx, dxdtdty, dxdtdtz,
                coupling_parameter_,
                H[i*Nv_prec_+j]);

                
                H[i*Nv_prec_+j] *= muwi[i]*muwj[j];

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

                HH(r_0, r_1, r_2,
                px[i*Nv_prec_+j], py[i*Nv_prec_+j], pz[i*Nv_prec_+j],
                nx[i*Nv_prec_+j], ny[i*Nv_prec_+j], nz[i*Nv_prec_+j],
                dxdsx[i*Nv_prec_+j] * der_eta(si[i], flag_u_loc) * 0.5 * (u_b_loc - u_a_loc), dxdsy[i*Nv_prec_+j] * der_eta(si[i], flag_u_loc) * 0.5 * (u_b_loc - u_a_loc), dxdsz[i*Nv_prec_+j] * der_eta(si[i], flag_u_loc) * 0.5 * (u_b_loc - u_a_loc),
                dxdtx[i*Nv_prec_+j] * der_eta(tj[j], flag_v_loc) * 0.5 * (v_b_loc - v_a_loc), dxdty[i*Nv_prec_+j] * der_eta(tj[j], flag_v_loc) * 0.5 * (v_b_loc - v_a_loc), dxdtz[i*Nv_prec_+j] * der_eta(tj[j], flag_v_loc) * 0.5 * (v_b_loc - v_a_loc),
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                coupling_parameter_,
                H[i*Nv_prec_+j]);
                
                H[i*Nv_prec_+j] *= muwi[i]*muwj[j];

            }
        }

    }

    std::vector<std::complex<double>> Tn_mat(Nu_prec_*Nu_int_), Tm_mat(Nv_prec_*Nv_int_);
    
    cheb_evals<Nu_prec_, Nu_int_>(ui, Tn_mat);
    cheb_evals<Nv_prec_, Nv_int_>(vj, Tm_mat);

    std::vector<std::complex<double>> matprod1(Nv_prec_*Nu_int_, {0.0, 0.0});

    std::complex<double> one(1.0, 0.0);
    std::complex<double> zero(0.0, 0.0);

    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                Nu_int_, Nv_prec_, Nu_prec_, 
                &one, &Tn_mat[0], Nu_prec_, &H[0], Nv_prec_, &zero, &matprod1[0], Nv_prec_);

    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                Nu_int_, Nv_int_, Nv_prec_,
                &one, &matprod1[0], Nv_prec_, &Tm_mat[0], Nv_prec_, &zero, &prec[0], Nv_int_);
    
}


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


void Solver::get_coeffs(const long long q, const std::complex<double>* phi, 
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


void Solver::compute_coeffs(std::complex<double>* phi,
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

void Solver::compute_intensities_patch(const long long npatch,
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

void Solver::compute_intensities(const std::complex<double>* phi,
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

void Solver::redistribute_data_RP(const std::vector<std::complex<double>>& intensities, 
                                    std::complex<double>* rhs)
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



