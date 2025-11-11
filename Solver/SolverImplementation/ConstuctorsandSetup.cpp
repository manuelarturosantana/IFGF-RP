#include "../solver2.h"

using namespace Eigen;

void Solver::compute_parallel_parameters() 
{

    int flag;

    MPI_Initialized(&flag);

    if (!flag) {
        throw std::logic_error("Cannot use this code without MPI initialization\n");
    }

    MPI_Comm_size(mpi_comm_, &comm_size_);
    MPI_Comm_rank(mpi_comm_, &comm_rank_);

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

    recv_counts_ = std::vector<int>(comm_size_);
    displs_ = std::vector<int>(comm_size_, 0);

    recv_counts_2_ = std::vector<int>(comm_size_);
    displs_2_ = std::vector<int>(comm_size_, 0);

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

        if (!fin.is_open()) {
            std::cerr << "Error: could not open file "
                    << DIRECTORY + FILE_NAME + std::to_string(i+1) + ".txt" << "\n";
            // handle error
            std::exit(1);
        }

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

        interp_surface_[i] = elem;

        fin.close();

    }

}

// Helper function for compute patch splits
double inline norm3(double x1, double x2, double x3) {
    return std::sqrt(x1*x1 + x2*x2 + x3*x3);
}

void Solver::compute_number_patch_splits() {

    // std::cout << "Before (Qx, Qy) : (" << Qx_ << ", " << Qy_ << ")" << std::endl;

    // Track umaxlen, vmaxlen
    double umaxlen = 0.0; double vmaxlen = 0.0;

    // track number of fejer nodes/weights
    int num_u_nodes = 0; int num_v_nodes = 0;
    // init vector of fejer nodes/weights to a guess (say 200?)
    std::vector<double> fejer_weights_u; fejer_weights_u.reserve(200);
    std::vector<double>fejer_weights_v; fejer_weights_v.reserve(200);
    // Track the length of the four sides of the patch
    double lenum1, lenup1, lenvm1, lenvp1;

    // loop over patches for this rank assuming the setup for the parallel parameters has
    // already be called;
    for (long long i = orig_patch_low_; i <= orig_patch_up_; i++) {
         // set m1,u1,vm1,v1 to 0.0
        lenum1 = 0.0; lenup1 = 0.0; lenvm1 = 0.0; lenvp1 = 0.0;

        // Alias for the current patch
        InterpPatch& cp = interp_surface_[i];
            // for both u/v 
    //    if the number of nodes is different, recompute fejer weights
        if (cp.imax != fejer_weights_u.size()) {
            fejer_weights_u.resize(cp.imax);
            fejerquadrature1(cp.uNodes, fejer_weights_u, cp.imax);
        }

        if (cp.jmax != fejer_weights_v.size()) {
            fejer_weights_v.resize(cp.jmax);
            fejerquadrature1(cp.vNodes, fejer_weights_v, cp.jmax);
        }

        // loop and grab the gradient, the use it to compute the length of the four sides
        for (int uind = 0; uind < cp.imax; uind++) {
            int indm1 = uind;
            int indp1 = cp.imax * (cp.jmax - 1) + uind;
            
            lenum1 += fejer_weights_u[uind] * norm3(cp.dxdu[indm1], cp.dydu[indm1], cp.dzdu[indm1]);
            lenup1 += fejer_weights_u[uind] * norm3(cp.dxdu[indp1], cp.dydu[indp1], cp.dzdu[indp1]);

        }

        for (int vind = 0; vind < cp.jmax; vind++) {
            int indm1 = vind * cp.imax;
            int indp1 = vind * cp.imax + (cp.imax - 1);
            
            lenvm1 += fejer_weights_v[vind] * norm3(cp.dxdu[indm1], cp.dydu[indm1], cp.dzdu[indm1]);
            lenvp1 += fejer_weights_v[vind] * norm3(cp.dxdu[indp1], cp.dydu[indp1], cp.dzdu[indp1]);

        }

        umaxlen = std::max({umaxlen, lenum1, lenup1});
        vmaxlen = std::max({vmaxlen, lenvm1, lenvp1});
    }

    MPI_Allreduce(MPI_IN_PLACE, &umaxlen, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
    MPI_Allreduce(MPI_IN_PLACE, &vmaxlen, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);

    // reduction on the maxes 
    // Compute number of patch splits per numer of wave lengths

    double lambda = (2.0 * M_PI) / std::real(wavenumber_);
    Qx_ = std::ceil(umaxlen / (lambda * num_wl_per_patch));
    Qy_ = std::ceil(vmaxlen / (lambda * num_wl_per_patch));

    // std::cout << "After (Qx, Qy) : (" << Qx_ << ", " << Qy_ << ")" << std::endl;
    // std::exit(1);

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

    std::vector<int> flags_domain_u(patch_up_ - patch_low_);
    std::vector<int> flags_domain_v(patch_up_ - patch_low_);
    
    #pragma omp parallel for
    for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

        const long long q = patch_num / (Qx_*Qy_);
        const int q_x = (patch_num % (Qx_*Qy_)) / Qy_;
        const int q_y = (patch_num % (Qx_*Qy_)) % Qy_;

        const long long pos = patch_num - patch_low_;

        flags_domain_u[pos] = 0; 
        flags_domain_v[pos] = 0; 

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

    flags_domain_u_all_.resize(Q_ * Qx_*Qy_);
    flags_domain_v_all_.resize(Q_ * Qx_*Qy_);

    MPI_Allgatherv(&flags_domain_u[0], patch_up_-patch_low_, MPI_INT, &flags_domain_u_all_[0], &recv_counts_[0], &displs_[0], MPI_INT, mpi_comm_);
    MPI_Allgatherv(&flags_domain_v[0], patch_up_-patch_low_, MPI_INT, &flags_domain_v_all_[0], &recv_counts_[0], &displs_[0], MPI_INT, mpi_comm_);

}


void Solver::compute_discretization_domain() 
{

    std::vector<double> disc_points_x(point_up_ - point_low_);
    std::vector<double> disc_points_y(point_up_ - point_low_);
    std::vector<double> disc_points_z(point_up_ - point_low_);

    std::vector<double> norm_points_x(point_up_ - point_low_);
    std::vector<double> norm_points_y(point_up_ - point_low_);
    std::vector<double> norm_points_z(point_up_ - point_low_);

    std::vector<double> dsdtjac(point_up_ - point_low_);
    
    #pragma omp parallel for
    for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

        const long long q = patch_num / (Qx_*Qy_);
        const int q_x = (patch_num % (Qx_*Qy_)) / Qy_;
        const int q_y = (patch_num % (Qx_*Qy_)) % Qy_;

        const long long pos = patch_num - patch_low_;

        const int flag_u_loc = flags_domain_u_all_[patch_num];
        const int flag_v_loc = flags_domain_v_all_[patch_num];

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
                    dsdtjac[position] = dsi[nu] * dtj[nv] * jacobian_q(SPHERE_RADIUS,si, tj, q);

                }
            }

        } else {             

            InterpPatch interpolation_patch = interp_surface_[q];

            // std::cout << "made it to load interpolated surface " << std::endl;
            // std::exit(1);

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

    disc_points_x_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
    disc_points_y_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
    disc_points_z_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

    norm_points_x_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
    norm_points_y_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
    norm_points_z_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

    dsdtjac_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

    MPI_Allgatherv(&disc_points_x[0], point_up_-point_low_, MPI_DOUBLE, &disc_points_x_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&disc_points_y[0], point_up_-point_low_, MPI_DOUBLE, &disc_points_y_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&disc_points_z[0], point_up_-point_low_, MPI_DOUBLE, &disc_points_z_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);

    MPI_Allgatherv(&norm_points_x[0], point_up_-point_low_, MPI_DOUBLE, &norm_points_x_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&norm_points_y[0], point_up_-point_low_, MPI_DOUBLE, &norm_points_y_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&norm_points_z[0], point_up_-point_low_, MPI_DOUBLE, &norm_points_z_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);

    MPI_Allgatherv(&dsdtjac[0], point_up_-point_low_, MPI_DOUBLE, &dsdtjac_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);

}            

void Solver::compute_coupling_parameter() 
{

    double loc_max_x = std::numeric_limits<double>::lowest();
    double loc_max_y = std::numeric_limits<double>::lowest();
    double loc_max_z = std::numeric_limits<double>::lowest();

    double loc_min_x = std::numeric_limits<double>::max();
    double loc_min_y = std::numeric_limits<double>::max();
    double loc_min_z = std::numeric_limits<double>::max();

    for (long long i = point_low_; i < point_up_; i++) {

        const double pointx = disc_points_x_all_[i];
        const double pointy = disc_points_y_all_[i];
        const double pointz = disc_points_z_all_[i];

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

    MPI_Allreduce(&loc_max_x, &global_max_x, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
    MPI_Allreduce(&loc_max_y, &global_max_y, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
    MPI_Allreduce(&loc_max_z, &global_max_z, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
    MPI_Allreduce(&loc_min_x, &global_min_x, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_);
    MPI_Allreduce(&loc_min_y, &global_min_y, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_);
    MPI_Allreduce(&loc_min_z, &global_min_z, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_);

    const double scatDiam = std::max({global_max_x - global_min_x, global_max_y - global_min_y, global_max_z - global_min_z});
    
    coupling_parameter_ = std::max(3.0, scatDiam / std::sqrt(lambda_.real()*lambda_.real() + lambda_.imag()*lambda_.imag()));
    
    if (use_flipped_eta && wavenumber_.real() > 0) {
        coupling_parameter_ *= -1;
    }

    

}

void Solver::compute_near_singular_patches_estimate()
{

    std::vector<double> x_min(patch_up_ - patch_low_, std::numeric_limits<double>::max());
    std::vector<double> y_min(patch_up_ - patch_low_, std::numeric_limits<double>::max());
    std::vector<double> z_min(patch_up_ - patch_low_, std::numeric_limits<double>::max());

    std::vector<double> x_max(patch_up_ - patch_low_, std::numeric_limits<double>::lowest());
    std::vector<double> y_max(patch_up_ - patch_low_, std::numeric_limits<double>::lowest());
    std::vector<double> z_max(patch_up_ - patch_low_, std::numeric_limits<double>::lowest());
    
    #pragma omp parallel for
    for (long long patch_num = 0; patch_num < patch_up_ - patch_low_; patch_num++) {

        for (int nu = 0; nu < Nu_int_; nu++) {
            for (int nv = 0; nv < Nv_int_; nv++) {

                const long long i = (patch_num + patch_low_) * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;               

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

        const double Lx = x_max[patch_num] - x_min[patch_num];
        const double Ly = y_max[patch_num] - y_min[patch_num];
        const double Lz = z_max[patch_num] - z_min[patch_num];

        const double Lmax = std::max({Lx, Ly, Lz});
        const double factor = Lmax * 2.0E-2; 

        x_min[patch_num] -= factor;
        x_max[patch_num] += factor;
        y_min[patch_num] -= factor;
        y_max[patch_num] += factor;
        z_min[patch_num] -= factor;
        z_max[patch_num] += factor;

    }    

    MPI_Barrier(mpi_comm_);        

    std::vector<double> x_min_all(Q_ * Qx_*Qy_);
    std::vector<double> y_min_all(Q_ * Qx_*Qy_);
    std::vector<double> z_min_all(Q_ * Qx_*Qy_);

    std::vector<double> x_max_all(Q_ * Qx_*Qy_);
    std::vector<double> y_max_all(Q_ * Qx_*Qy_);
    std::vector<double> z_max_all(Q_ * Qx_*Qy_);

    MPI_Allgatherv(&x_min[0], patch_up_-patch_low_, MPI_DOUBLE, &x_min_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&y_min[0], patch_up_-patch_low_, MPI_DOUBLE, &y_min_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&z_min[0], patch_up_-patch_low_, MPI_DOUBLE, &z_min_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);

    MPI_Allgatherv(&x_max[0], patch_up_-patch_low_, MPI_DOUBLE, &x_max_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&y_max[0], patch_up_-patch_low_, MPI_DOUBLE, &y_max_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&z_max[0], patch_up_-patch_low_, MPI_DOUBLE, &z_max_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);

    start_sing_and_near_sing_patches_estimate_ = std::vector<long long>(patch_up_ - patch_low_);
    size_sing_and_near_sing_patches_estimate_ = std::vector<long long>(patch_up_ - patch_low_);

    std::vector<std::unordered_set<long long>> patches_not_in_rank(comm_size_);

    long long index = 0;
    
    for (long long patch_num = 0; patch_num < patch_up_ - patch_low_; patch_num++) {

        long long total = 0;
        
        for (long long patch_num_2 = 0; patch_num_2 < Q_ * Qx_*Qy_; patch_num_2++) {

            if ((x_max[patch_num] >= x_min_all[patch_num_2]) && (x_max_all[patch_num_2] >= x_min[patch_num]) &&
                (y_max[patch_num] >= y_min_all[patch_num_2]) && (y_max_all[patch_num_2] >= y_min[patch_num]) &&
                (z_max[patch_num] >= z_min_all[patch_num_2]) && (z_max_all[patch_num_2] >= z_min[patch_num])) {

                sing_and_near_sing_patches_estimate_.push_back(patch_num_2);
                total++;

                if (GEOMETRY != 0) {
                
                    const long long q = patch_num_2 / (Qx_ * Qy_);

                    if (interp_surface_.count(q) == 0) {
                    
                        int rank;

                        for (int i = 0; i < comm_size_+1; i++) {

                            if ((patch_num_2 >= split_points_[i]) && (patch_num_2 < split_points_[i+1])) {

                                rank = i;
                                break;

                            }

                        }
                    
                        patches_not_in_rank[rank].insert(q); 

                    }               

                }

            }

        }

        start_sing_and_near_sing_patches_estimate_[patch_num] = index;
        size_sing_and_near_sing_patches_estimate_[patch_num] = total;
        index += total;

    }

    // This trick frees the memory of the vectors
    std::vector<double>().swap(x_min);
    std::vector<double>().swap(y_min);
    std::vector<double>().swap(z_min);

    std::vector<double>().swap(x_min_all);
    std::vector<double>().swap(y_min_all);
    std::vector<double>().swap(z_min_all);

    std::vector<double>().swap(x_max);
    std::vector<double>().swap(y_max);
    std::vector<double>().swap(z_max);

    std::vector<double>().swap(x_max_all);
    std::vector<double>().swap(y_max_all);
    std::vector<double>().swap(z_max_all);

    MPI_Barrier(mpi_comm_);
    
    if (GEOMETRY != 0) {

        for (int rank = 0; rank < comm_size_; rank++) {

            const long long size = patches_not_in_rank[rank].size();
            const std::vector<long long> elements(patches_not_in_rank[rank].begin(), patches_not_in_rank[rank].end());

            std::vector<long long> size_all(comm_size_);

            MPI_Gather(&size, 1, MPI_LONG_LONG, &size_all[0], 1, MPI_LONG_LONG, rank, mpi_comm_);

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

            MPI_Gatherv(&elements[0], size, MPI_LONG_LONG, &patches[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, rank, mpi_comm_);

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

            MPI_Scatterv(&vector_imax[0], &recv_counts[0], &displs[0], MPI_INT, &vector_imax_loc[0], size, MPI_INT, rank, mpi_comm_);
            MPI_Scatterv(&vector_jmax[0], &recv_counts[0], &displs[0], MPI_INT, &vector_jmax_loc[0], size, MPI_INT, rank, mpi_comm_);
            MPI_Scatterv(&vector_kmax[0], &recv_counts[0], &displs[0], MPI_INT, &vector_kmax_loc[0], size, MPI_INT, rank, mpi_comm_);

            MPI_Scatterv(&vector_zoneID[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, &vector_zoneID_loc[0], size, MPI_LONG_LONG, rank, mpi_comm_);

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

            MPI_Gather(&imax_total, 1, MPI_INT, &recv_counts_imax[0], 1, MPI_INT, rank, mpi_comm_);
            MPI_Gather(&jmax_total, 1, MPI_INT, &recv_counts_jmax[0], 1, MPI_INT, rank, mpi_comm_);
            MPI_Gather(&imaxjmaxkmax_total, 1, MPI_INT, &recv_counts_imaxjmaxkmax[0], 1, MPI_INT, rank, mpi_comm_);                

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

            MPI_Scatterv(&vector_uNodes[0], &recv_counts_imax[0], &displs_imax[0], MPI_DOUBLE, &vector_uNodes_loc[0], imax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_vNodes[0], &recv_counts_jmax[0], &displs_jmax[0], MPI_DOUBLE, &vector_vNodes_loc[0], jmax_total, MPI_DOUBLE, rank, mpi_comm_);

            MPI_Scatterv(&vector_uWeights[0], &recv_counts_imax[0], &displs_imax[0], MPI_DOUBLE, &vector_uWeights_loc[0], imax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_vWeights[0], &recv_counts_jmax[0], &displs_jmax[0], MPI_DOUBLE, &vector_vWeights_loc[0], jmax_total, MPI_DOUBLE, rank, mpi_comm_);

            MPI_Scatterv(&vector_x[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_x_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_y[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_y_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_z[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_z_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            
            MPI_Scatterv(&vector_mask[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_INT, &vector_mask_loc[0], imaxjmaxkmax_total, MPI_INT, rank, mpi_comm_);
            
            MPI_Scatterv(&vector_dxdu[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dxdu_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_dydu[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dydu_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_dzdu[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dzdu_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            
            MPI_Scatterv(&vector_dxdv[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dxdv_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_dydv[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dydv_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_dzdv[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dzdv_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);

            MPI_Scatterv(&vector_dS[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dS_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);

            MPI_Scatterv(&vector_nuX[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_nuX_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_nuY[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_nuY_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            MPI_Scatterv(&vector_nuZ[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_nuZ_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);

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
            
}

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

                parametrization_q(SPHERE_RADIUS,SPHERE_CENTER,s, t, q, px, py, pz);
                normal_q(SPHERE_RADIUS,s, t, q, nx, ny, nz);
                dxds_q(SPHERE_RADIUS,s, t, q, dxdsx, dxdsy, dxdsz);
                dxdt_q(SPHERE_RADIUS,s, t, q, dxdtx, dxdty, dxdtz);
                dxdsds_q(SPHERE_RADIUS,s, t, q, dxdsdsx, dxdsdsy, dxdsdsz);
                dxdsdt_q(SPHERE_RADIUS,s, t, q, dxdsdtx, dxdsdty, dxdsdtz);
                dxdtdt_q(SPHERE_RADIUS,s, t, q, dxdtdtx, dxdtdty, dxdtdtz);

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

                HH(r_0, r_1, r_2,
                px, py, pz,
                nx, ny, nz,
                dxdsx, dxdsy, dxdsz,
                dxdtx, dxdty, dxdtz,
                dxdsdsx, dxdsdsy, dxdsdsz,
                dxdsdtx, dxdsdty, dxdsdtz,
                dxdtdtx, dxdtdty, dxdtdtz,
                coupling_parameter_,
                wavenumber_,
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
                dxdsx[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc), dxdsy[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc), dxdsz[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc),
                dxdtx[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc), dxdty[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc), dxdtz[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc),
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                coupling_parameter_,
                wavenumber_,
                H[i*Nv_prec_+j]);

                H[i*Nv_prec_+j] *= muwi[i]*muwj[j];

            }
        }

    }

    std::vector<std::complex<double>> Tn_mat(Nu_prec_*Nu_int_), Tm_mat(Nv_prec_*Nv_int_);
    
    cheb_evals(ui, Tn_mat, Nu_prec_, Nu_int_);
    cheb_evals(vj, Tm_mat, Nv_prec_, Nv_int_);

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

void Solver::setup(bool timing) 
{

    MPI_Barrier(mpi_comm_);

    double start_1 = MPI_Wtime();


    compute_parallel_parameters();

    double end_1 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute parallel parameters: "  << end_1 - start_1 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);

    double start_2 = MPI_Wtime();

    if (GEOMETRY != 0) {
        load_interpolated_surface();
    }

    double end_2 = MPI_Wtime();

    if (GEOMETRY != 0 && timing && comm_rank_ == 0) {
        std::cout << "Time load interpolated surface: " << end_2 - start_2 << " seconds\n";
       
        print_max_RSS();

    }

    double start_3 = MPI_Wtime();

    if (split_patch_by_wavenumber_ && GEOMETRY!=0) {
        compute_number_patch_splits();
        compute_parallel_parameters();
        // There is a gotcha! here. interp_surface is not broadcast across the rankes
        // orig patch low and orig patch high can change as Q_x, Q_y change,
        // due to extra patches. To fix this so we need to read in the surface again after
        // computing parallel paramter.s
        interp_surface_.clear();
        load_interpolated_surface();

    }

    double end_3 = MPI_Wtime();

    if (split_patch_by_wavenumber_ && timing && comm_rank_ == 0) {
        std::cout << "Time compute number patch splits: " << end_3 - start_3 << " seconds\n";
        std::cout << "(Qx, Qy) : (" << Qx_ << ", " << Qy_ << ")" << std::endl;
       
        print_max_RSS();
    }

    MPI_Barrier(mpi_comm_);

    double start_4 = MPI_Wtime();

    compute_fejer_nodes_and_weights();

    double end_4 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute Fejer nodes and weights: " << end_4 - start_4 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);

    double start_5 = MPI_Wtime();

    if (!USE_OVERSAMPLING) {

        compute_chebyshev_evaluations();

    }

    double end_5 = MPI_Wtime();

    if (!USE_OVERSAMPLING && timing && comm_rank_ == 0) {

        std::cout << "Time compute Chebyshev evaluations: " << end_5 - start_5 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);

    double start_6 = MPI_Wtime();

    compute_flags_domain();

    double end_6 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute flags domain: " << end_6 - start_6 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);

    double start_7 = MPI_Wtime();

    compute_discretization_domain();

    double end_7 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute discretization domain: " << end_7 - start_7 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);

    double start_8 = MPI_Wtime();

    if (init_compute_coupling_param_) {
        compute_coupling_parameter();  
    }

    double end_8 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute coupling parameter: " << end_8 - start_8 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);   

    double start_9 = MPI_Wtime();   

    compute_near_singular_patches_estimate();   

    double end_9 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute near singular patches estimate: " << end_9 - start_9 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);   

    double start_10 = MPI_Wtime();

    if (USE_OVERSAMPLING) {

        compute_singular_points();

    } else {

        compute_precomputations();

    }

    double end_10 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute precomputations data: " << end_10 - start_10 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);  

    double start_11 = MPI_Wtime();

    if (USE_ACCELERATOR) {

        setup_IFGF_choose_levels();

        // if (comm_rank_ == 0) {
        //     std::cout << "The number of levels chosen by the algorithm is " << nlevels_ << std::endl;
        // }

        // create_IFGF_object();

    }

    double end_11 = MPI_Wtime();

    if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

        std::cout << "Time create IFGF object: " << end_11 - start_11 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);  


    // double start_12 = MPI_Wtime();

    // if (USE_ACCELERATOR) {

    //     compute_new_order_points_RP();                

    // }

    // double end_12 = MPI_Wtime();

    // if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

    //     std::cout << "Time compute new order points RP: " << end_12 - start_12 << " seconds\n";
    //     print_max_RSS();

    // }

    // MPI_Barrier(mpi_comm_);  


    // double start_13 = MPI_Wtime();

    // if (USE_ACCELERATOR) {

    //     check_patch_in_neighbours();

    // }

    // double end_13 = MPI_Wtime();

    // if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

    //     std::cout << "Time check patch in neighbours: " << end_13 - start_13 << " seconds\n";
    //     print_max_RSS();

    // }

    // MPI_Barrier(mpi_comm_); 


    double start_14 = MPI_Wtime();

    if (USE_ACCELERATOR && USE_HIGH_ORDER) {

        initialize_indexes_HO();

    }

    double end_14 = MPI_Wtime();

    if (USE_ACCELERATOR && USE_HIGH_ORDER && timing && comm_rank_ == 0) {

        std::cout << "Time initialize indexes HO: " << end_14 - start_14 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_); 

    
    double start_15 = MPI_Wtime();

    if (USE_OVERSAMPLING) {

        create_fftw_objects();

    }

    double end_15 = MPI_Wtime();

    if (USE_OVERSAMPLING && timing && comm_rank_ == 0) {

        std::cout << "Time create FFTW objects: " << end_15 - start_15 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_); 

}

void Solver::init_solver(const bool timing, 
        const std::complex<double> k, 
        const int* n_pts_per_patch,
        const int* n_split_per_patch,
        const int* n_pts_sing_int,
        const double proximity_box_size,
        const int n_levels_IFGF, 
        const MPI_Comm& mpi_comm) {
    
    

    Nu_int_ = n_pts_per_patch[0];
    Nv_int_ = n_pts_per_patch[1];

    Nu_prec_ = n_pts_sing_int[0];
    Nv_prec_ = n_pts_sing_int[1];
    
    // Flag for splitting the patch by wavenumber
    if (n_split_per_patch[0] == -1) {
        split_patch_by_wavenumber_ = true;
        Qx_ = 1;
        Qy_ = 1;
    } else {
        Qx_ = n_split_per_patch[0];
        Qy_ = n_split_per_patch[1];
    }



    wavenumber_ = k;
    lambda_ = 2.0 * M_PI / wavenumber_;

    proximity_ = proximity_box_size;
    nlevels_ = n_levels_IFGF;

    mpi_comm_ = mpi_comm;



    setup(timing);

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



Solver::~Solver() {

    if (USE_OVERSAMPLING) {

        fftw_destroy_plan(plan_original_patch_);
        fftw_destroy_plan(plan_oversampled_patch_);

        fftw_free(original_patch_);
        fftw_free(oversampled_patch_);

    }

}    