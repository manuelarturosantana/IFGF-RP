#include "../solver2.h"

#include <filesystem> // Required for std::filesystem
#include <algorithm>  // Required for std::count_if

void Solver::compute_parallel_parameters()
{

    int flag;

    MPI_Initialized(&flag);

    if (!flag) {
        throw std::logic_error("Cannot use this code without MPI initialization\n");
    }

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size_);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank_);

    split_points_ = std::vector<long long>(comm_size_ + 1);
    split_points_2_ = std::vector<long long>(comm_size_ + 1);

    const long long patches_per_rank = (Q_ * Qx_*Qy_) / comm_size_;
    long long remaining_patches = (Q_ * Qx_*Qy_) % comm_size_;

    split_points_[0] = 0;
    split_points_[comm_size_] = Q_ * Qx_*Qy_;

    split_points_2_[0] = 0;
    split_points_2_[comm_size_] = Q_ * Qx_*Qy_ * Nu_int_*Nv_int_;

    for (int i = 1; i < comm_size_; i++) {

        split_points_[i] = split_points_[i-1] + patches_per_rank;

        if (remaining_patches > 0) {
            split_points_[i]++;
            remaining_patches--;
        }

        split_points_2_[i] = split_points_[i] * Nu_int_*Nv_int_;

    }

    recv_counts_ = std::vector<MPI_Count>(comm_size_);
    displs_ = std::vector<MPI_Aint>(comm_size_, 0);

    recv_counts_2_ = std::vector<MPI_Count>(comm_size_);
    displs_2_ = std::vector<MPI_Aint>(comm_size_, 0);

    for (int rank = 0; rank < comm_size_; rank++) {

        recv_counts_[rank] = split_points_[rank + 1] - split_points_[rank];
        recv_counts_2_[rank] = split_points_2_[rank + 1] - split_points_2_[rank];

        if (rank != 0) {

            displs_[rank] = displs_[rank - 1] + recv_counts_[rank - 1];
            displs_2_[rank] = displs_2_[rank - 1] + recv_counts_2_[rank - 1];

        }

    }

    patch_low_ = split_points_[comm_rank_];
    patch_up_ = split_points_[comm_rank_ + 1];

    orig_patch_low_ = patch_low_ / (Qx_*Qy_);
    orig_patch_up_ = (patch_up_ - 1) / (Qx_*Qy_);

    point_low_ = patch_low_ * Nu_int_*Nv_int_;
    point_up_ = patch_up_ * Nu_int_*Nv_int_;

    omp_set_num_threads(NTHREADS);   
    
}

void Solver::load_interpolated_surface() 
{
    
    #pragma omp parallel for
    for (long long i = orig_patch_low_; i <= orig_patch_up_; i++) {

        std::ifstream fin(DIRECTORY + FILE_NAME + std::to_string(i+1) + ".txt");

        InterpPatch elem;
        
        fin >> elem.imax >> elem.jmax >> elem.kmax >> elem.zoneID;

        elem.uNodes = std::vector<double>(elem.imax);

        for (int j = 0; j < elem.imax; j++) {

            fin >> elem.uNodes[j];

        }

        elem.vNodes = std::vector<double>(elem.jmax);

        for (int j = 0; j < elem.jmax; j++) {

            fin >> elem.vNodes[j];

        }

        elem.uWeights = barycentric_weights(elem.uNodes);

        elem.vWeights = barycentric_weights(elem.vNodes);

        elem.x = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (int j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.x[j];

        }

        elem.y = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (int j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.y[j];

        }

        elem.z = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.z[j];

        }

        elem.mask = std::vector<int>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.mask[j];

        }

        elem.dxdu = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.dxdu[j];

        }

        elem.dydu = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.dydu[j];

        }

        elem.dzdu = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.dzdu[j];

        }

        elem.dxdv = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.dxdv[j];

        }

        elem.dydv = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.dydv[j];

        }

        elem.dzdv = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.dzdv[j];

        }

        elem.dS = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.dS[j];

        }

        elem.nuX = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.nuX[j];

        }

        elem.nuY = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.nuY[j];

        }

        elem.nuZ = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

        for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

            fin >> elem.nuZ[j];

        }

        interp_surface_[i] = std::move(elem);

        fin.close();

    }

}

void Solver::compute_fejer_nodes_and_weights()
{

    fejer_nodes_u_int_ = std::vector<double>(Nu_int_);
    fejer_weights_u_int_ = std::vector<double>(Nu_int_);

    fejer_nodes_v_int_ = std::vector<double>(Nv_int_);
    fejer_weights_v_int_ = std::vector<double>(Nv_int_);

    fejer_nodes_u_prec_ = std::vector<double>(Nu_prec_);
    fejer_weights_u_prec_ = std::vector<double>(Nu_prec_);

    fejer_nodes_v_prec_ = std::vector<double>(Nv_prec_);
    fejer_weights_v_prec_ = std::vector<double>(Nv_prec_);

    fejerquadrature1(fejer_nodes_u_int_, fejer_weights_u_int_, Nu_int_);
    fejerquadrature1(fejer_nodes_v_int_, fejer_weights_v_int_, Nv_int_);
    fejerquadrature1(fejer_nodes_u_prec_, fejer_weights_u_prec_, Nu_prec_);
    fejerquadrature1(fejer_nodes_v_prec_, fejer_weights_v_prec_, Nv_prec_);

    fejer_weights_u_v_int_ = std::vector<double>(Nu_int_*Nv_int_);
    
    for (int i = 0; i < Nu_int_; i++) {
        for (int j = 0; j < Nv_int_; j++) {

            fejer_weights_u_v_int_[i * Nv_int_ + j] = fejer_weights_u_int_[i] * fejer_weights_v_int_[j];

        }
    }

}

void Solver::compute_chebyshev_evaluations()
{

    Tn_ = std::vector<std::complex<double>>(Nu_int_*Nu_int_);
    Tm_ = std::vector<std::complex<double>>(Nv_int_*Nv_int_);

    cheb_evals(fejer_nodes_u_int_, Tn_, Nu_int_, Nu_int_);
    cheb_evals(fejer_nodes_v_int_, Tm_, Nv_int_, Nv_int_); 

}   

void Solver::compute_flags_domain() 
{

    const long long num_patches = patch_up_ - patch_low_;

    std::vector<int> flags_domain_u(num_patches, 0);
    std::vector<int> flags_domain_v(num_patches, 0);

    #pragma omp parallel for
    for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

        const long long q = patch_num / (Qx_*Qy_);
        const long long q_x = (patch_num % (Qx_*Qy_)) / Qy_;
        const long long q_y = (patch_num % (Qx_*Qy_)) % Qy_;

        const long long pos = patch_num - patch_low_;

        const double u_a = -1.0 + q_x * 2.0 / Qx_;
        const double u_b = -1.0 + (q_x + 1) * 2.0 / Qx_;
        const double v_a = -1.0 + q_y * 2.0 / Qy_;
        const double v_b = -1.0 + (q_y + 1) * 2.0 / Qy_;

        if (u_a == -1.0 && EDGE_FLAG_U_A[q]) {
            flags_domain_u[pos] += 1; 
        }

        if (u_b == 1.0 && EDGE_FLAG_U_B[q]) {
            flags_domain_u[pos] += 2;
        }

        if (v_a == -1.0 && EDGE_FLAG_V_A[q]) {
            flags_domain_v[pos] += 1;
        }

        if (v_b == 1.0 && EDGE_FLAG_V_B[q]) {
            flags_domain_v[pos] += 2; 
        }

    }

    if (USE_ACCELERATOR) {

        flags_domain_u_all_ = std::move(flags_domain_u);
        flags_domain_v_all_ = std::move(flags_domain_v);

    } else {

        flags_domain_u_all_ = std::vector<int>(Q_ * Qx_*Qy_);
        flags_domain_v_all_ = std::vector<int>(Q_ * Qx_*Qy_);

        MPI_Allgatherv_c(flags_domain_u.data(), num_patches, MPI_INT, flags_domain_u_all_.data(), &recv_counts_[0], &displs_[0], MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv_c(flags_domain_v.data(), num_patches, MPI_INT, flags_domain_v_all_.data(), &recv_counts_[0], &displs_[0], MPI_INT, MPI_COMM_WORLD);

    } 

}

void Solver::compute_discretization_domain() 
{

    const long long num_points = point_up_ - point_low_;

    std::vector<double> disc_points_x(num_points);
    std::vector<double> disc_points_y(num_points);
    std::vector<double> disc_points_z(num_points);

    std::vector<double> norm_points_x(num_points);
    std::vector<double> norm_points_y(num_points);
    std::vector<double> norm_points_z(num_points);

    std::vector<double> dsdtjac(num_points);
    
    #pragma omp parallel for
    for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

        const long long q = patch_num / (Qx_*Qy_);
        const long long q_x = (patch_num % (Qx_*Qy_)) / Qy_;
        const long long q_y = (patch_num % (Qx_*Qy_)) % Qy_;

        const long long pos = patch_num - patch_low_;

        const long long pos_flags = USE_ACCELERATOR ? pos : patch_num;

        const int flag_u_loc = flags_domain_u_all_[pos_flags];
        const int flag_v_loc = flags_domain_v_all_[pos_flags];

        const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
        const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
        const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
        const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

        std::vector<double> vec_si(Nu_int_), dsi(Nu_int_);
        std::vector<double> vec_tj(Nv_int_), dtj(Nv_int_);

        for (int i = 0; i < Nu_int_; i++) {

            double si = eta(fejer_nodes_u_int_[i], flag_u_loc);
            si = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, si);
            vec_si[i] = si;

            dsi[i] = der_eta(fejer_nodes_u_int_[i], flag_u_loc) * (u_b_loc - u_a_loc) * 0.5;

        }

        for (int j = 0; j < Nv_int_; j++) {

            double tj = eta(fejer_nodes_v_int_[j], flag_v_loc);                            
            tj = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, tj);
            vec_tj[j] = tj;

            dtj[j] = der_eta(fejer_nodes_v_int_[j], flag_v_loc) * (v_b_loc - v_a_loc) * 0.5;

        }

        if (GEOMETRY == 0) {

            for (int nu = 0; nu < Nu_int_; nu++) {
                for (int nv = 0; nv < Nv_int_; nv++) {

                    const long long position = pos * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;

                    const double si = vec_si[nu];
                    const double tj = vec_tj[nv];

                    parametrization_q(SPHERE_RADIUS, SPHERE_CENTER, si, tj, q, disc_points_x[position], disc_points_y[position], disc_points_z[position]);
                    normal_q(SPHERE_RADIUS, si, tj, q, norm_points_x[position], norm_points_y[position], norm_points_z[position]);
                    dsdtjac[position] = dsi[nu] * dtj[nv] * jacobian_q(SPHERE_RADIUS, si, tj, q);

                }
            }

        } else {             

            InterpPatch interpolation_patch = interp_surface_[q];

            lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
                                    interpolation_patch.uWeights, interpolation_patch.vWeights,
                                    interpolation_patch.x, vec_si, vec_tj, &disc_points_x[pos * Nu_int_*Nv_int_]);
            lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
                                    interpolation_patch.uWeights, interpolation_patch.vWeights,
                                    interpolation_patch.y, vec_si, vec_tj, &disc_points_y[pos * Nu_int_*Nv_int_]);
            lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
                                    interpolation_patch.uWeights, interpolation_patch.vWeights,
                                    interpolation_patch.z, vec_si, vec_tj, &disc_points_z[pos * Nu_int_*Nv_int_]);

            lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
                                    interpolation_patch.uWeights, interpolation_patch.vWeights,
                                    interpolation_patch.nuX, vec_si, vec_tj, &norm_points_x[pos * Nu_int_*Nv_int_]);
            lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
                                    interpolation_patch.uWeights, interpolation_patch.vWeights,
                                    interpolation_patch.nuY, vec_si, vec_tj, &norm_points_y[pos * Nu_int_*Nv_int_]);
            lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
                                    interpolation_patch.uWeights, interpolation_patch.vWeights,
                                    interpolation_patch.nuZ, vec_si, vec_tj, &norm_points_z[pos * Nu_int_*Nv_int_]);

            lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
                                    interpolation_patch.uWeights, interpolation_patch.vWeights,
                                    interpolation_patch.dS, vec_si, vec_tj, &dsdtjac[pos * Nu_int_*Nv_int_]);

            for (int nu = 0; nu < Nu_int_; nu++) {
                for (int nv = 0; nv < Nv_int_; nv++) {

                    const long long position = pos * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;
                    
                    dsdtjac[position] *= dsi[nu] * dtj[nv];

                }
            }               

        }  

    }

    if (USE_ACCELERATOR) {

        disc_points_x_all_ = std::move(disc_points_x);
        disc_points_y_all_ = std::move(disc_points_y);
        disc_points_z_all_ = std::move(disc_points_z);

        norm_points_x_all_ = std::move(norm_points_x);
        norm_points_y_all_ = std::move(norm_points_y);
        norm_points_z_all_ = std::move(norm_points_z);

        dsdtjac_all_ = std::move(dsdtjac);

    } else {

        disc_points_x_all_ = std::vector<double>(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
        disc_points_y_all_ = std::vector<double>(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
        disc_points_z_all_ = std::vector<double>(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

        norm_points_x_all_ = std::vector<double>(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
        norm_points_y_all_ = std::vector<double>(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
        norm_points_z_all_ = std::vector<double>(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

        dsdtjac_all_ = std::vector<double>(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

        MPI_Allgatherv_c(&disc_points_x[0], point_up_-point_low_, MPI_DOUBLE, &disc_points_x_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv_c(&disc_points_y[0], point_up_-point_low_, MPI_DOUBLE, &disc_points_y_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv_c(&disc_points_z[0], point_up_-point_low_, MPI_DOUBLE, &disc_points_z_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, MPI_COMM_WORLD);

        MPI_Allgatherv_c(&norm_points_x[0], point_up_-point_low_, MPI_DOUBLE, &norm_points_x_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv_c(&norm_points_y[0], point_up_-point_low_, MPI_DOUBLE, &norm_points_y_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv_c(&norm_points_z[0], point_up_-point_low_, MPI_DOUBLE, &norm_points_z_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, MPI_COMM_WORLD);

        MPI_Allgatherv_c(&dsdtjac[0], point_up_-point_low_, MPI_DOUBLE, &dsdtjac_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, MPI_COMM_WORLD);

    }

}

void Solver::compute_coupling_parameter() 
{

    double loc_max_x = std::numeric_limits<double>::lowest();
    double loc_max_y = std::numeric_limits<double>::lowest();
    double loc_max_z = std::numeric_limits<double>::lowest();

    double loc_min_x = std::numeric_limits<double>::max();
    double loc_min_y = std::numeric_limits<double>::max();
    double loc_min_z = std::numeric_limits<double>::max();

    #pragma omp parallel for reduction(max:loc_max_x,loc_max_y,loc_max_z) \
                                reduction(min:loc_min_x,loc_min_y,loc_min_z)
    for (long long i = point_low_; i < point_up_; i++) {

        long long idx = USE_ACCELERATOR ? i - point_low_ : i;

        const double pointx = disc_points_x_all_[idx];
        const double pointy = disc_points_y_all_[idx];
        const double pointz = disc_points_z_all_[idx];

        if (loc_max_x < pointx) loc_max_x = pointx;
        if (loc_max_y < pointy) loc_max_y = pointy;
        if (loc_max_z < pointz) loc_max_z = pointz;
        
        if (loc_min_x > pointx) loc_min_x = pointx;
        if (loc_min_y > pointy) loc_min_y = pointy;
        if (loc_min_z > pointz) loc_min_z = pointz;

    }

    double global_max_x;
    double global_max_y;
    double global_max_z;

    double global_min_x;
    double global_min_y;
    double global_min_z;

    MPI_Allreduce(&loc_max_x, &global_max_x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_max_y, &global_max_y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_max_z, &global_max_z, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_min_x, &global_min_x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_min_y, &global_min_y, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_min_z, &global_min_z, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    const double scatDiam = std::max({global_max_x - global_min_x, global_max_y - global_min_y, global_max_z - global_min_z});
    
    coupling_parameter_ = std::max(3.0, scatDiam / LAMBDA);

}


void Solver::compute_near_singular_patches_estimate()
{

    std::vector<double> x_min(patch_up_ - patch_low_, std::numeric_limits<double>::max());
    std::vector<double> y_min(patch_up_ - patch_low_, std::numeric_limits<double>::max());
    std::vector<double> z_min(patch_up_ - patch_low_, std::numeric_limits<double>::max());

    std::vector<double> x_max(patch_up_ - patch_low_, std::numeric_limits<double>::lowest());
    std::vector<double> y_max(patch_up_ - patch_low_, std::numeric_limits<double>::lowest());
    std::vector<double> z_max(patch_up_ - patch_low_, std::numeric_limits<double>::lowest());

    double loc_min_x = std::numeric_limits<double>::max();
    double loc_min_y = std::numeric_limits<double>::max();
    double loc_min_z = std::numeric_limits<double>::max();
    
    double loc_max_x = std::numeric_limits<double>::lowest();
    double loc_max_y = std::numeric_limits<double>::lowest();
    double loc_max_z = std::numeric_limits<double>::lowest();   
    
    double cell_size_x_loc = 0.0;
    double cell_size_y_loc = 0.0;
    double cell_size_z_loc = 0.0;
    
    #pragma omp parallel for reduction(max:loc_max_x,loc_max_y,loc_max_z) \
                                reduction(min:loc_min_x,loc_min_y,loc_min_z) \
                                reduction(+:cell_size_x_loc,cell_size_y_loc,cell_size_z_loc)
    for (long long patch_num = 0; patch_num < patch_up_ - patch_low_; patch_num++) {

        const long long idx = USE_ACCELERATOR ? patch_num : (patch_num + patch_low_);

        for (int nu = 0; nu < Nu_int_; nu++) {
            for (int nv = 0; nv < Nv_int_; nv++) {                        

                const long long i = idx * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;               

                const double pointx = disc_points_x_all_[i];
                const double pointy = disc_points_y_all_[i];
                const double pointz = disc_points_z_all_[i];

                if (x_min[patch_num] > pointx) x_min[patch_num] = pointx;
                if (y_min[patch_num] > pointy) y_min[patch_num] = pointy;
                if (z_min[patch_num] > pointz) z_min[patch_num] = pointz;
                
                if (x_max[patch_num] < pointx) x_max[patch_num] = pointx;
                if (y_max[patch_num] < pointy) y_max[patch_num] = pointy;
                if (z_max[patch_num] < pointz) z_max[patch_num] = pointz;  

            }
        } 

        x_max[patch_num] += 1e-10;
        y_max[patch_num] += 1e-10;
        z_max[patch_num] += 1e-10;

        x_min[patch_num] -= 1e-10;
        y_min[patch_num] -= 1e-10;
        z_min[patch_num] -= 1e-10;
        
        if (loc_min_x > x_min[patch_num]) loc_min_x = x_min[patch_num];
        if (loc_min_y > y_min[patch_num]) loc_min_y = y_min[patch_num];
        if (loc_min_z > z_min[patch_num]) loc_min_z = z_min[patch_num];
        
        if (loc_max_x < x_max[patch_num]) loc_max_x = x_max[patch_num];
        if (loc_max_y < y_max[patch_num]) loc_max_y = y_max[patch_num];
        if (loc_max_z < z_max[patch_num]) loc_max_z = z_max[patch_num];   
        
        cell_size_x_loc += (x_max[patch_num] - x_min[patch_num]);
        cell_size_y_loc += (y_max[patch_num] - y_min[patch_num]);
        cell_size_z_loc += (z_max[patch_num] - z_min[patch_num]);

    }

    double global_max_x;
    double global_max_y;
    double global_max_z;

    double global_min_x;
    double global_min_y;
    double global_min_z;

    double cell_size_x;
    double cell_size_y;
    double cell_size_z;

    MPI_Allreduce(&loc_max_x, &global_max_x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_max_y, &global_max_y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_max_z, &global_max_z, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    MPI_Allreduce(&loc_min_x, &global_min_x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_min_y, &global_min_y, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_min_z, &global_min_z, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Allreduce(&cell_size_x_loc, &cell_size_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&cell_size_y_loc, &cell_size_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&cell_size_z_loc, &cell_size_z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    cell_size_x /= (Q_ * Qx_*Qy_);
    cell_size_y /= (Q_ * Qx_*Qy_);
    cell_size_z /= (Q_ * Qx_*Qy_);

    const long long Nx = std::max(1LL, static_cast<long long>(std::ceil((global_max_x - global_min_x) / cell_size_x)));
    const long long Ny = std::max(1LL, static_cast<long long>(std::ceil((global_max_y - global_min_y) / cell_size_y)));
    const long long Nz = std::max(1LL, static_cast<long long>(std::ceil((global_max_z - global_min_z) / cell_size_z)));

    std::vector<std::unordered_map<long long, std::vector<long long>>> grids_omp(omp_get_max_threads());

    #pragma omp parallel
    {

        #pragma omp for
        for (long long patch_num = 0; patch_num < patch_up_ - patch_low_; patch_num++) {

            long long ix_min = std::clamp(static_cast<long long>((x_min[patch_num] - global_min_x) / cell_size_x), 0LL, Nx - 1);
            long long ix_max = std::clamp(static_cast<long long>((x_max[patch_num] - global_min_x) / cell_size_x), 0LL, Nx - 1);
            long long iy_min = std::clamp(static_cast<long long>((y_min[patch_num] - global_min_y) / cell_size_y), 0LL, Ny - 1);
            long long iy_max = std::clamp(static_cast<long long>((y_max[patch_num] - global_min_y) / cell_size_y), 0LL, Ny - 1);
            long long iz_min = std::clamp(static_cast<long long>((z_min[patch_num] - global_min_z) / cell_size_z), 0LL, Nz - 1);
            long long iz_max = std::clamp(static_cast<long long>((z_max[patch_num] - global_min_z) / cell_size_z), 0LL, Nz - 1);

            for (long long ix = ix_min; ix <= ix_max; ++ix) {
                for (long long iy = iy_min; iy <= iy_max; ++iy) {
                    for (long long iz = iz_min; iz <= iz_max; ++iz) {

                        long long key = ((ix * Ny + iy) * Nz + iz);

                        grids_omp[omp_get_thread_num()][key].push_back(patch_num + patch_low_);

                    }
                }
            }


        }

    }

    std::map<long long, std::vector<long long>> grid;

    for (auto & local_grid : grids_omp) {

        for (auto & [key, vec] : local_grid) {

            auto & target = grid[key];
            target.insert(target.end(), vec.begin(), vec.end());

        }

        local_grid.clear();

    }

    grids_omp.clear();

    std::vector<long long> keys_loc, starts_loc, values_loc;

    keys_loc.reserve(grid.size());
    starts_loc.reserve(grid.size());

    long long total_values = 0;

    for (const auto & [key, vec] : grid) {

        keys_loc.push_back(key);
        starts_loc.push_back(total_values);
        values_loc.insert(values_loc.end(), vec.begin(), vec.end());

        total_values += vec.size();

    }

    std::vector<MPI_Count> recv_counts(comm_size_);
    std::vector<MPI_Aint> recv_displs(comm_size_);

    MPI_Count total_keys = grid.size();

    MPI_Allgather(&total_keys, 1, MPI_COUNT, recv_counts.data(), 1, MPI_COUNT, MPI_COMM_WORLD);
    
    MPI_Aint total_keys_all = 0;

    for (int i = 0; i < comm_size_; i++) {

        recv_displs[i] = total_keys_all;

        total_keys_all += recv_counts[i];

    }

    std::vector<long long> keys_all(total_keys_all), starts_all(total_keys_all);

    MPI_Allgatherv_c(keys_loc.data(), total_keys, MPI_LONG_LONG, keys_all.data(), recv_counts.data(), recv_displs.data(), MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgatherv_c(starts_loc.data(), total_keys, MPI_LONG_LONG, starts_all.data(), recv_counts.data(), recv_displs.data(), MPI_LONG_LONG, MPI_COMM_WORLD);
    
    std::vector<long long>().swap(keys_loc);
    std::vector<long long>().swap(starts_loc);

    std::vector<long long> total_values_all(comm_size_);

    MPI_Allgather(&total_values, 1, MPI_LONG_LONG, total_values_all.data(), 1, MPI_LONG_LONG, MPI_COMM_WORLD);

    MPI_Win win_values;

    MPI_Win_create(values_loc.data(), total_values * sizeof(long long), sizeof(long long), MPI_INFO_NULL, MPI_COMM_WORLD, &win_values);
    
    MPI_Win_lock_all(MPI_MODE_NOCHECK, win_values);

    struct RMARequest {
        int rank;
        long long key;
        long long start;
        long long size;
        long long local_offset;
    };

    std::vector<RMARequest> requests;
    requests.reserve(grid.size() * (comm_size_ - 1));

    for (int rank = 0; rank < comm_size_; rank++) {

        if (rank == comm_rank_) continue;

        const auto recv_start = recv_displs[rank];
        const auto recv_end = recv_displs[rank] + recv_counts[rank];

        const auto keys_begin = keys_all.begin() + recv_start;
        const auto keys_end = keys_all.begin() + recv_end;

        for (const auto & [key, vec] : grid) {

            auto it = std::lower_bound(keys_begin, keys_end, key);

            if (it == keys_end || *it != key) continue;

            const long long idx = std::distance(keys_begin, it);

            const long long start = starts_all[recv_displs[rank] + idx];
            const long long next_start = (idx == (recv_counts[rank] - 1)) ? total_values_all[rank] : starts_all[recv_start + idx + 1];
            const long long size = next_start - start;

            const long long local_offset = grid[key].size();
            grid[key].resize(local_offset + size);

            requests.push_back({rank, key, start, size, local_offset});

        }

    }

    std::vector<MPI_Request> mpi_reqs;
    mpi_reqs.reserve(requests.size());

    for (const auto & req : requests) {

        MPI_Request mpi_req;

        long long *target_ptr = grid[req.key].data() + req.local_offset;

        MPI_Rget(target_ptr, req.size, MPI_LONG_LONG, req.rank, req.start, req.size, MPI_LONG_LONG, win_values, &mpi_req);

        mpi_reqs.push_back(mpi_req);

    }

    MPI_Waitall(static_cast<int>(mpi_reqs.size()), mpi_reqs.data(), MPI_STATUSES_IGNORE);

    MPI_Win_unlock_all(win_values);

    MPI_Win_free(&win_values);    

    std::vector<RMARequest>().swap(requests);
    std::vector<MPI_Request>().swap(mpi_reqs);

    std::vector<long long>().swap(keys_all);
    std::vector<long long>().swap(starts_all);

    std::vector<long long>().swap(values_loc);    

    start_sing_and_near_sing_patches_estimate_ = std::vector<long long>(patch_up_ - patch_low_);
    size_sing_and_near_sing_patches_estimate_ = std::vector<long long>(patch_up_ - patch_low_);

    std::vector<std::unordered_set<long long>> patches_not_in_rank_acc(comm_size_);
    std::vector<std::unordered_set<long long>> patches_not_in_rank(comm_size_);

    #pragma omp parallel
    {

    std::vector<std::unordered_set<long long>> patches_not_in_rank_acc_thread(comm_size_);
    std::vector<std::unordered_set<long long>> patches_not_in_rank_thread(comm_size_);

    std::vector<long long> sing_and_near_sing_patches_estimate_thread;
    
    #pragma omp for
    for (long long patch_num = 0; patch_num < patch_up_ - patch_low_; patch_num++) {  

        long long total = 0;

        long long ix_min = std::clamp(static_cast<long long>((x_min[patch_num] - global_min_x) / cell_size_x), 1LL, Nx - 1);
        long long ix_max = std::clamp(static_cast<long long>((x_max[patch_num] - global_min_x) / cell_size_x), 0LL, Nx - 2);
        long long iy_min = std::clamp(static_cast<long long>((y_min[patch_num] - global_min_y) / cell_size_y), 1LL, Ny - 1);
        long long iy_max = std::clamp(static_cast<long long>((y_max[patch_num] - global_min_y) / cell_size_y), 0LL, Ny - 2);
        long long iz_min = std::clamp(static_cast<long long>((z_min[patch_num] - global_min_z) / cell_size_z), 1LL, Nz - 1);
        long long iz_max = std::clamp(static_cast<long long>((z_max[patch_num] - global_min_z) / cell_size_z), 0LL, Nz - 2);

        std::vector<long long> patches_unique;

        for (long long ix = ix_min - 1; ix <= ix_max + 1; ++ix) {
            for (long long iy = iy_min - 1; iy <= iy_max + 1; ++iy) {
                for (long long iz = iz_min - 1; iz <= iz_max + 1; ++iz) {

                    long long key = ((ix * Ny + iy) * Nz + iz);
                    
                    if (grid.count(key) == 0) continue;

                    patches_unique.insert(patches_unique.end(), grid[key].begin(), grid[key].end());            
                    
                }
            }
        }

        std::sort(patches_unique.begin(), patches_unique.end());
        auto last_unique = std::unique(patches_unique.begin(), patches_unique.end());
        patches_unique.erase(last_unique, patches_unique.end());

        for (const auto & patch_num_2 : patches_unique) {

            sing_and_near_sing_patches_estimate_thread.push_back(patch_num_2);
            total++;

            auto it = std::upper_bound(split_points_.begin(), split_points_.end(), patch_num_2);        
            int rank = static_cast<int>(std::distance(split_points_.begin(), it)) - 1;        
            rank = std::max(0, std::min(rank, comm_size_ - 1));

            const long long q = patch_num_2 / (Qx_ * Qy_);

            if (USE_ACCELERATOR && (rank != comm_rank_)) {

                patches_not_in_rank_acc_thread[rank].insert(patch_num_2);

            }

            if ((GEOMETRY != 0) && (interp_surface_.count(q) == 0))  {
                
                patches_not_in_rank_thread[rank].insert(q);           

            }

        }    

        size_sing_and_near_sing_patches_estimate_[patch_num] = total;

    }

    #pragma omp for ordered
    for (int i = 0; i < NTHREADS; i++) {

        #pragma omp ordered
        {

            sing_and_near_sing_patches_estimate_.insert(sing_and_near_sing_patches_estimate_.end(), sing_and_near_sing_patches_estimate_thread.begin(), sing_and_near_sing_patches_estimate_thread.end());

            for (int ii = 0; ii < comm_size_; ii++) {

                patches_not_in_rank_acc[ii].insert(patches_not_in_rank_acc_thread[ii].begin(), patches_not_in_rank_acc_thread[ii].end());
                patches_not_in_rank[ii].insert(patches_not_in_rank_thread[ii].begin(), patches_not_in_rank_thread[ii].end());

            }

            std::vector<std::unordered_set<long long>>().swap(patches_not_in_rank_acc_thread);
            std::vector<std::unordered_set<long long>>().swap(patches_not_in_rank_thread);
            std::vector<long long>().swap(sing_and_near_sing_patches_estimate_thread);

        }

    }

    }

    std::vector<double>().swap(x_min);
    std::vector<double>().swap(y_min);
    std::vector<double>().swap(z_min);

    std::vector<double>().swap(x_max);
    std::vector<double>().swap(y_max);
    std::vector<double>().swap(z_max);

    grid.clear();

    long long index = 0;

    for (long long patch_num = 0; patch_num < patch_up_ - patch_low_; patch_num++) { 

        start_sing_and_near_sing_patches_estimate_[patch_num] = index;
        index += size_sing_and_near_sing_patches_estimate_[patch_num];

    }

    if (USE_ACCELERATOR) {

        const long long num_points = Nu_int_ * Nv_int_;
        const long long num_patches = patch_up_ - patch_low_;

        const size_t elems_per_patch = 2 + 4 * num_points;
        const size_t bytes_per_patch = sizeof(int) * 2 + sizeof(double) * 4 * num_points;
    
        std::vector<char> patch_buffer_local(bytes_per_patch * num_patches);
        char* base_ptr = patch_buffer_local.data();

        for (long long i = 0; i < num_patches; ++i) {

            char* ptr = base_ptr + i * bytes_per_patch;

            int* flag_u_ptr = reinterpret_cast<int*>(ptr);
            int* flag_v_ptr = flag_u_ptr + 1;

            double* x_ptr = reinterpret_cast<double*>(flag_v_ptr + 1);
            double* y_ptr = x_ptr + num_points;
            double* z_ptr = y_ptr + num_points;

            double* dsdtjac_ptr = z_ptr + num_points;
    
            *flag_u_ptr = flags_domain_u_all_[i];
            *flag_v_ptr = flags_domain_v_all_[i];
    
            std::copy_n(&disc_points_x_all_[i * num_points], num_points, x_ptr);
            std::copy_n(&disc_points_y_all_[i * num_points], num_points, y_ptr);
            std::copy_n(&disc_points_z_all_[i * num_points], num_points, z_ptr);

            std::copy_n(&dsdtjac_all_[i * num_points], num_points, dsdtjac_ptr);

        }

        MPI_Win win;
        MPI_Info info;
        MPI_Info_create(&info);
        MPI_Info_set(info, "accumulate_ordering", "none");
        MPI_Info_set(info, "same_size", "true");
        MPI_Info_set(info, "same_disp_unit", "true");

        MPI_Win_create(patch_buffer_local.data(), patch_buffer_local.size(), sizeof(char), info, MPI_COMM_WORLD, &win);

        MPI_Info_free(&info);
        
        MPI_Win_lock_all(MPI_MODE_NOCHECK, win);

        std::vector<char> recv_patch(bytes_per_patch);

        for (int rank = 0; rank < comm_size_; ++rank) {

            if (rank == comm_rank_) continue;

            for (const auto& patch : patches_not_in_rank_acc[rank]) {

                MPI_Aint disp = static_cast<MPI_Aint>(patch - split_points_[rank]) * bytes_per_patch;

                MPI_Get(recv_patch.data(), bytes_per_patch, MPI_BYTE, rank, disp, bytes_per_patch, MPI_BYTE, win);

                MPI_Win_flush(rank, win);

                const int* flag_u_ptr = reinterpret_cast<int*>(recv_patch.data());
                const int* flag_v_ptr = flag_u_ptr + 1;

                const double* x_ptr = reinterpret_cast<const double*>(flag_v_ptr + 1);
                const double* y_ptr = x_ptr + num_points;
                const double* z_ptr = y_ptr + num_points;
                
                const double* dsdtjac_ptr = z_ptr + num_points;

                flags_domain_u_not_in_rank_[patch] = *flag_u_ptr;
                flags_domain_v_not_in_rank_[patch] = *flag_v_ptr;

                std::vector<double> x(x_ptr, x_ptr + num_points);
                std::vector<double> y(y_ptr, y_ptr + num_points);
                std::vector<double> z(z_ptr, z_ptr + num_points);

                std::vector<double> dsdtjac(dsdtjac_ptr, dsdtjac_ptr + num_points);

                disc_points_x_not_in_rank_[patch] = std::move(x);
                disc_points_y_not_in_rank_[patch] = std::move(y);
                disc_points_z_not_in_rank_[patch] = std::move(z);

                dsdtjac_not_in_rank_[patch] = std::move(dsdtjac);

            }

        }

        MPI_Win_unlock_all(win);
        MPI_Win_free(&win);

    }

    std::vector<std::unordered_set<long long>>().swap(patches_not_in_rank_acc);
    
    if (GEOMETRY != 0) {

        for (int rank = 0; rank < comm_size_; rank++) {

            const long long size = patches_not_in_rank[rank].size();
            const std::vector<long long> elements(patches_not_in_rank[rank].begin(), patches_not_in_rank[rank].end());

            std::vector<long long> size_all(comm_size_);

            MPI_Gather(&size, 1, MPI_LONG_LONG, &size_all[0], 1, MPI_LONG_LONG, rank, MPI_COMM_WORLD);

            std::vector<long long> patches;
            
            std::vector<int> recv_counts(comm_size_);
            std::vector<int> displs(comm_size_, 0);

            if (comm_rank_ == rank) {

                long long total_elements = 0;

                for (int i = 0; i < comm_size_; i++) {

                    recv_counts[i] = size_all[i];

                    if (i != 0) {

                        displs[i] = displs[i-1] + recv_counts[i-1];

                    }

                    total_elements += size_all[i];

                }

                patches.resize(total_elements);

            }

            MPI_Gatherv(&elements[0], size, MPI_LONG_LONG, &patches[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, rank, MPI_COMM_WORLD);

            std::vector<int> vector_imax;
            std::vector<int> vector_jmax;
            std::vector<int> vector_kmax;

            std::vector<long long> vector_zoneID;

            std::vector<double> vector_uNodes;
            std::vector<double> vector_vNodes;

            std::vector<double> vector_uWeights;
            std::vector<double> vector_vWeights;

            std::vector<double> vector_x, vector_y, vector_z;
            std::vector<int> vector_mask;
            std::vector<double> vector_dxdu, vector_dydu, vector_dzdu;
            std::vector<double> vector_dxdv, vector_dydv, vector_dzdv;
            std::vector<double> vector_dS;
            std::vector<double> vector_nuX, vector_nuY, vector_nuZ;

            if (comm_rank_ == rank) {
                
                for (long long i = 0; i < patches.size(); i++) {

                    InterpPatch elem = interp_surface_[patches[i]];

                    vector_imax.push_back(elem.imax);
                    vector_jmax.push_back(elem.jmax);
                    vector_kmax.push_back(elem.kmax);

                    vector_zoneID.push_back(elem.zoneID);

                    vector_uNodes.insert(vector_uNodes.end(), elem.uNodes.begin(), elem.uNodes.end());
                    vector_vNodes.insert(vector_vNodes.end(), elem.vNodes.begin(), elem.vNodes.end());

                    vector_uWeights.insert(vector_uWeights.end(), elem.uWeights.begin(), elem.uWeights.end());
                    vector_vWeights.insert(vector_vWeights.end(), elem.vWeights.begin(), elem.vWeights.end());

                    vector_x.insert(vector_x.end(), elem.x.begin(), elem.x.end());
                    vector_y.insert(vector_y.end(), elem.y.begin(), elem.y.end());
                    vector_z.insert(vector_z.end(), elem.z.begin(), elem.z.end());

                    vector_mask.insert(vector_mask.end(), elem.mask.begin(), elem.mask.end());

                    vector_dxdu.insert(vector_dxdu.end(), elem.dxdu.begin(), elem.dxdu.end());
                    vector_dydu.insert(vector_dydu.end(), elem.dydu.begin(), elem.dydu.end());
                    vector_dzdu.insert(vector_dzdu.end(), elem.dzdu.begin(), elem.dzdu.end());

                    vector_dxdv.insert(vector_dxdv.end(), elem.dxdv.begin(), elem.dxdv.end());
                    vector_dydv.insert(vector_dydv.end(), elem.dydv.begin(), elem.dydv.end());
                    vector_dzdv.insert(vector_dzdv.end(), elem.dzdv.begin(), elem.dzdv.end());

                    vector_dS.insert(vector_dS.end(), elem.dS.begin(), elem.dS.end());

                    vector_nuX.insert(vector_nuX.end(), elem.nuX.begin(), elem.nuX.end());
                    vector_nuY.insert(vector_nuY.end(), elem.nuY.begin(), elem.nuY.end());
                    vector_nuZ.insert(vector_nuZ.end(), elem.nuZ.begin(), elem.nuZ.end());

                }

            }
            
            std::vector<int> vector_imax_loc(size);
            std::vector<int> vector_jmax_loc(size);
            std::vector<int> vector_kmax_loc(size);

            std::vector<long long> vector_zoneID_loc(size);

            MPI_Scatterv(&vector_imax[0], &recv_counts[0], &displs[0], MPI_INT, &vector_imax_loc[0], size, MPI_INT, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_jmax[0], &recv_counts[0], &displs[0], MPI_INT, &vector_jmax_loc[0], size, MPI_INT, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_kmax[0], &recv_counts[0], &displs[0], MPI_INT, &vector_kmax_loc[0], size, MPI_INT, rank, MPI_COMM_WORLD);

            MPI_Scatterv(&vector_zoneID[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, &vector_zoneID_loc[0], size, MPI_LONG_LONG, rank, MPI_COMM_WORLD);

            int imax_total = 0;
            int jmax_total = 0;
            int imaxjmaxkmax_total = 0;

            for (long long i = 0; i < size; i++) {

                imax_total += vector_imax_loc[i];
                jmax_total += vector_jmax_loc[i];
                imaxjmaxkmax_total += vector_imax_loc[i] * vector_jmax_loc[i] * vector_kmax_loc[i];

            }

            std::vector<int> recv_counts_imax(comm_size_);
            std::vector<int> recv_counts_jmax(comm_size_);
            std::vector<int> recv_counts_imaxjmaxkmax(comm_size_);

            MPI_Gather(&imax_total, 1, MPI_INT, &recv_counts_imax[0], 1, MPI_INT, rank, MPI_COMM_WORLD);
            MPI_Gather(&jmax_total, 1, MPI_INT, &recv_counts_jmax[0], 1, MPI_INT, rank, MPI_COMM_WORLD);
            MPI_Gather(&imaxjmaxkmax_total, 1, MPI_INT, &recv_counts_imaxjmaxkmax[0], 1, MPI_INT, rank, MPI_COMM_WORLD);                

            std::vector<int> displs_imax(comm_size_, 0);
            std::vector<int> displs_jmax(comm_size_, 0);
            std::vector<int> displs_imaxjmaxkmax(comm_size_, 0);

            if (comm_rank_ == rank) {

                for (int i = 1; i < comm_size_; i++) {

                    displs_imax[i] = displs_imax[i-1] + recv_counts_imax[i-1];
                    displs_jmax[i] = displs_jmax[i-1] + recv_counts_jmax[i-1];
                    displs_imaxjmaxkmax[i] = displs_imaxjmaxkmax[i-1] + recv_counts_imaxjmaxkmax[i-1];

                }

            }

            std::vector<double> vector_uNodes_loc(imax_total);
            std::vector<double> vector_vNodes_loc(jmax_total);

            std::vector<double> vector_uWeights_loc(imax_total);
            std::vector<double> vector_vWeights_loc(jmax_total);

            std::vector<double> vector_x_loc(imaxjmaxkmax_total), vector_y_loc(imaxjmaxkmax_total), vector_z_loc(imaxjmaxkmax_total);
            std::vector<int> vector_mask_loc(imaxjmaxkmax_total);
            std::vector<double> vector_dxdu_loc(imaxjmaxkmax_total), vector_dydu_loc(imaxjmaxkmax_total), vector_dzdu_loc(imaxjmaxkmax_total);
            std::vector<double> vector_dxdv_loc(imaxjmaxkmax_total), vector_dydv_loc(imaxjmaxkmax_total), vector_dzdv_loc(imaxjmaxkmax_total);
            std::vector<double> vector_dS_loc(imaxjmaxkmax_total);
            std::vector<double> vector_nuX_loc(imaxjmaxkmax_total), vector_nuY_loc(imaxjmaxkmax_total), vector_nuZ_loc(imaxjmaxkmax_total);

            MPI_Scatterv(&vector_uNodes[0], &recv_counts_imax[0], &displs_imax[0], MPI_DOUBLE, &vector_uNodes_loc[0], imax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_vNodes[0], &recv_counts_jmax[0], &displs_jmax[0], MPI_DOUBLE, &vector_vNodes_loc[0], jmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            MPI_Scatterv(&vector_uWeights[0], &recv_counts_imax[0], &displs_imax[0], MPI_DOUBLE, &vector_uWeights_loc[0], imax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_vWeights[0], &recv_counts_jmax[0], &displs_jmax[0], MPI_DOUBLE, &vector_vWeights_loc[0], jmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            MPI_Scatterv(&vector_x[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_x_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_y[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_y_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_z[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_z_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            
            MPI_Scatterv(&vector_mask[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_INT, &vector_mask_loc[0], imaxjmaxkmax_total, MPI_INT, rank, MPI_COMM_WORLD);
            
            MPI_Scatterv(&vector_dxdu[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dxdu_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_dydu[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dydu_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_dzdu[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dzdu_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            
            MPI_Scatterv(&vector_dxdv[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dxdv_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_dydv[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dydv_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_dzdv[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dzdv_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            MPI_Scatterv(&vector_dS[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dS_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            MPI_Scatterv(&vector_nuX[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_nuX_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_nuY[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_nuY_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            MPI_Scatterv(&vector_nuZ[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_nuZ_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, MPI_COMM_WORLD);

            long long start_imax = 0;
            long long start_jmax = 0;
            long long start_imaxjmaxkmax = 0;

            for (long long i = 0; i < size; i++) {

                InterpPatch new_elem;

                new_elem.imax = vector_imax_loc[i];
                new_elem.jmax = vector_jmax_loc[i];
                new_elem.kmax = vector_kmax_loc[i];

                new_elem.zoneID = vector_zoneID_loc[i];

                new_elem.uNodes.insert(new_elem.uNodes.end(), vector_uNodes_loc.begin() + start_imax, vector_uNodes_loc.begin() + start_imax + new_elem.imax);
                new_elem.vNodes.insert(new_elem.vNodes.end(), vector_vNodes_loc.begin() + start_jmax, vector_vNodes_loc.begin() + start_jmax + new_elem.jmax);

                new_elem.uWeights.insert(new_elem.uWeights.end(), vector_uWeights_loc.begin() + start_imax, vector_uWeights_loc.begin() + start_imax + new_elem.imax);
                new_elem.vWeights.insert(new_elem.vWeights.end(), vector_vWeights_loc.begin() + start_jmax, vector_vWeights_loc.begin() + start_jmax + new_elem.jmax);

                new_elem.x.insert(new_elem.x.end(), vector_x_loc.begin() + start_imaxjmaxkmax, vector_x_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
                new_elem.y.insert(new_elem.y.end(), vector_y_loc.begin() + start_imaxjmaxkmax, vector_y_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
                new_elem.z.insert(new_elem.z.end(), vector_z_loc.begin() + start_imaxjmaxkmax, vector_z_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

                new_elem.mask.insert(new_elem.mask.end(), vector_mask_loc.begin() + start_imaxjmaxkmax, vector_mask_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

                new_elem.dxdu.insert(new_elem.dxdu.end(), vector_dxdu_loc.begin() + start_imaxjmaxkmax, vector_dxdu_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
                new_elem.dydu.insert(new_elem.dydu.end(), vector_dydu_loc.begin() + start_imaxjmaxkmax, vector_dydu_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
                new_elem.dzdu.insert(new_elem.dzdu.end(), vector_dzdu_loc.begin() + start_imaxjmaxkmax, vector_dzdu_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

                new_elem.dxdv.insert(new_elem.dxdv.end(), vector_dxdv_loc.begin() + start_imaxjmaxkmax, vector_dxdv_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
                new_elem.dydv.insert(new_elem.dydv.end(), vector_dydv_loc.begin() + start_imaxjmaxkmax, vector_dydv_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
                new_elem.dzdv.insert(new_elem.dzdv.end(), vector_dzdv_loc.begin() + start_imaxjmaxkmax, vector_dzdv_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

                new_elem.dS.insert(new_elem.dS.end(), vector_dS_loc.begin() + start_imaxjmaxkmax, vector_dS_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

                new_elem.nuX.insert(new_elem.nuX.end(), vector_nuX_loc.begin() + start_imaxjmaxkmax, vector_nuX_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
                new_elem.nuY.insert(new_elem.nuY.end(), vector_nuY_loc.begin() + start_imaxjmaxkmax, vector_nuY_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
                new_elem.nuZ.insert(new_elem.nuZ.end(), vector_nuZ_loc.begin() + start_imaxjmaxkmax, vector_nuZ_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

                interp_surface_[elements[i]] = new_elem;

                start_imax += new_elem.imax;
                start_jmax += new_elem.jmax;
                start_imaxjmaxkmax += new_elem.imax*new_elem.jmax*new_elem.kmax;

            }

        }

    }

    std::vector<std::unordered_set<long long>>().swap(patches_not_in_rank);
            
}

void Solver::setup(bool timing)
{

    MPI_Barrier(MPI_COMM_WORLD);

    double start_1 = MPI_Wtime();

    compute_parallel_parameters();

    double end_1 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute parallel parameters: " << end_1 - start_1 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    double start_2 = MPI_Wtime();

    if (GEOMETRY != 0) {

        load_interpolated_surface();

    }

    double end_2 = MPI_Wtime();

    if (GEOMETRY != 0 && timing && comm_rank_ == 0) {

        std::cout << "Time load interpolated surface: " << end_2 - start_2 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    double start_3 = MPI_Wtime();

    compute_fejer_nodes_and_weights();

    double end_3 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute Fejer nodes and weights: " << end_3 - start_3 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    double start_4 = MPI_Wtime();

    compute_chebyshev_evaluations();

    double end_4 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute Chebyshev evaluations: " << end_4 - start_4 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    double start_5 = MPI_Wtime();

    compute_flags_domain();

    double end_5 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute flags domain: " << end_5 - start_5 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    double start_6 = MPI_Wtime();

    compute_discretization_domain();

    double end_6 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute discretization domain: " << end_6 - start_6 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    double start_7 = MPI_Wtime();

    compute_coupling_parameter();  

    double end_7 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute coupling parameter: " << end_7 - start_7 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);  

    double start_8 = MPI_Wtime();   

    compute_near_singular_patches_estimate();   

    double end_8 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute near singular patches estimate: " << end_8 - start_8 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);  

    double start_9 = MPI_Wtime();

    compute_precomputations();

    double end_9 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute precomputations data: " << end_9 - start_9 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    double start_10 = MPI_Wtime();

    if (USE_ACCELERATOR) {

        create_IFGF_object();

    }

    double end_10 = MPI_Wtime();

    if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

        std::cout << "Time create IFGF object: " << end_10 - start_10 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);

    double start_11 = MPI_Wtime();

    if (USE_ACCELERATOR) {

        set_precomputations_data_IFGF();

    }

    double end_11 = MPI_Wtime();

    if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

        std::cout << "Time set precomputation data in IFGF: " << end_11 - start_11 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD); 

    double start_12 = MPI_Wtime();

    if (USE_ACCELERATOR) {

        compute_new_order_points_RP();                

    }

    double end_12 = MPI_Wtime();

    if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

        std::cout << "Time compute new order points RP: " << end_12 - start_12 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD); 

    double start_13 = MPI_Wtime();

    if (USE_ACCELERATOR) {

        check_patch_in_neighbours();

    }

    double end_13 = MPI_Wtime();

    if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

        std::cout << "Time check patch in neighbours: " << end_13 - start_13 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(MPI_COMM_WORLD);

}


Solver::Solver(
        double sphere_radius,
        double sphere_centerX, double sphere_centerY, double sphere_centerZ)
{

    SPHERE_RADIUS = sphere_radius;
    SPHERE_CENTER[0] = sphere_centerX; SPHERE_CENTER[1] = sphere_centerY; SPHERE_CENTER[2] = sphere_centerZ;
    // Q_ is the number of patches before splitting.
    Q_ = 6;
    EDGE_FLAG_U_A = std::vector<bool>(Q_, false);
    EDGE_FLAG_U_B = std::vector<bool>(Q_, false);
    EDGE_FLAG_V_A = std::vector<bool>(Q_, false);
    EDGE_FLAG_V_B = std::vector<bool>(Q_, false);

}



Solver::Solver( const std::string directory,
                const std::string file_prefix)
{

    DIRECTORY = directory;
    FILE_NAME = file_prefix;
    Q_ = 0;

    GEOMETRY = 1;

     try {
        // Create a directory_iterator for the specified path
        std::filesystem::directory_iterator dirIter(DIRECTORY);

        // Count regular files using std::count_if and a lambda expression
        Q_ = std::count_if(begin(dirIter), end(dirIter),
                                      [](const std::filesystem::directory_entry& entry) {
                                          return entry.is_regular_file();
                                      });

    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Error accessing directory: " << ex.what() << std::endl;
        
    }

    EDGE_FLAG_U_A = std::vector<bool>(Q_, false);
    EDGE_FLAG_U_B = std::vector<bool>(Q_, false);
    EDGE_FLAG_V_A = std::vector<bool>(Q_, false);
    EDGE_FLAG_V_B = std::vector<bool>(Q_, false);
               
}

void Solver::init_solver(const bool timing, const double k, MPI_Comm mpi_comm) {

    WAVE_NUMBER = k;
    mpi_comm_   = mpi_comm;
    
    Nu_int_ = N_PTS_PER_PATCH[0]; Nv_int_ = N_PTS_PER_PATCH[1];
    Nu_prec_ = N_PTS_SING_INT[0]; Nv_prec_ = N_PTS_SING_INT[1];
    
    Qx_ = N_SPLIT_PER_PATCH[0], Qy_ = N_SPLIT_PER_PATCH[1];


    setup(timing);

}


