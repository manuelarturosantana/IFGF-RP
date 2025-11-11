// #include "solver2.h"
// Things are split into files now so that compilation doesn't take forever while making changes


// using namespace Eigen;

// // int G_SEARCH_MAX_ITER = 50;
// // double G_SEARCH_TOL = 1E-12;

// void Solver::HH(const double p1x, const double p1y, const double p1z,
//         const double p2x, const double p2y, const double p2z,
//         const double nx, const double ny, const double nz,
//         const double dxdsx, const double dxdsy, const double dxdsz,
//         const double dxdtx, const double dxdty, const double dxdtz,
//         const double dxdsdsx, const double dxdsdsy, const double dxdsdsz,
//         const double dxdsdtx, const double dxdsdty, const double dxdsdtz,
//         const double dxdtdtx, const double dxdtdty, const double dxdtdtz,
//         const double coupling_parameter, std::complex<double> wavenumber,
//         std::complex<double>& solution)
// {

//     double solution_real, solution_imag;

//     const double norm_diff = std::sqrt((p2x - p1x)*(p2x - p1x) + (p2y - p1y)*(p2y - p1y) + (p2z - p1z)*(p2z - p1z));

//     if (norm_diff < 1e-14) {

//         solution_real = 0.0;
//         solution_imag = 0.0;

//     } else {   

//         const double val_cos = std::cos(wavenumber.real() * norm_diff);
//         const double val_sin = std::sin(wavenumber.real() * norm_diff);
//         const double val_const = std::exp(-wavenumber.imag() * norm_diff);

//         if (EQUATION_FORMULATION == 1) {

//             solution_real = val_const * val_cos / (4.0 * M_PI * norm_diff);
//             solution_imag = val_const * val_sin / (4.0 * M_PI * norm_diff);

//         } else if (EQUATION_FORMULATION == 2) {

//             // Compute beta = < r, n > / |r|^2
            
//             const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
//             double beta = rDotNorm / (norm_diff * norm_diff);

//             if ((norm_diff < 1.0E-5) && (std::abs(rDotNorm) < 1.0E-6)) {

//                 beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, nx, ny, nz, dxdsx, dxdsy, dxdsz, dxdtx, dxdty, dxdtz, 
//                                                 dxdsdsx, dxdsdsy, dxdsdsz, dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

//             }            

//             solution_real = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_cos / norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
//             solution_imag = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_sin / norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

//         } else {

//             // EQUATION_FORMULATION == 3 

//             const double solution_real_S = val_const * val_cos / (4.0 * M_PI * norm_diff);
//             const double solution_imag_S = val_const * val_sin / (4.0 * M_PI * norm_diff);

//             // Compute beta = < r, n > / |r|^2

//             const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
//             double beta = rDotNorm / (norm_diff * norm_diff);

//             if ((norm_diff < 1.0E-5) && (std::abs(rDotNorm) < 1.0E-6)) {

//                 beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, nx, ny, nz, dxdsx, dxdsy, dxdsz, dxdtx, dxdty, dxdtz, 
//                                                 dxdsdsx, dxdsdsy, dxdsdsz, dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

//             } 

//             const double solution_real_D = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_cos / norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
//             const double solution_imag_D = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_sin / norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

//             solution_real = solution_real_D + coupling_parameter * solution_imag_S;
//             solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

//         }

//         solution = std::complex<double>{solution_real, solution_imag};

//     }

// }

// void Solver::HH2(const double p1x, const double p1y, const double p1z,
//          const double p2x, const double p2y, const double p2z,
//          const double nx, const double ny, const double nz,
//          double coupling_parameter, std::complex<double> wavenumber,
//          std::complex<double>& solution)
// {

//     double solution_real, solution_imag;

//     const double norm_diff = std::sqrt((p2x - p1x)*(p2x - p1x) + (p2y - p1y)*(p2y - p1y) + (p2z - p1z)*(p2z - p1z));

//     if (norm_diff < 1e-14) {

//         solution_real = 0.0;
//         solution_imag = 0.0;

//     } else {   

//         const double val_cos = std::cos(wavenumber.real() * norm_diff);
//         const double val_sin = std::sin(wavenumber.real() * norm_diff);
//         const double val_const = std::exp(-wavenumber.imag() * norm_diff);

//         if (EQUATION_FORMULATION == 1) {

//             solution_real = val_const * val_cos / (4.0 * M_PI * norm_diff);
//             solution_imag = val_const * val_sin / (4.0 * M_PI * norm_diff);

//         } else if (EQUATION_FORMULATION == 2) {

//             // Compute beta = < r, n > / |r|^2
            
//             const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
//             double beta = rDotNorm / (norm_diff * norm_diff);

//             solution_real = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_cos / norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
//             solution_imag = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_sin / norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

//             if (USE_ACCELERATOR) {

//                 solution_real *= -1.0;
//                 solution_imag *= -1.0;

//             }

//         } else {

//             // EQUATION_FORMULATION == 3 

//             const double solution_real_S = val_const * val_cos / (4.0 * M_PI * norm_diff);
//             const double solution_imag_S = val_const * val_sin / (4.0 * M_PI * norm_diff);

//             // Compute beta = < r, n > / |r|^2

//             const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
//             double beta = rDotNorm / (norm_diff * norm_diff);

//             double solution_real_D = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_cos / norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
//             double solution_imag_D = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_sin / norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

//             if (USE_ACCELERATOR) {

//                 solution_real_D *= -1.0;
//                 solution_imag_D *= -1.0;

//             }

//             solution_real = solution_real_D + coupling_parameter * solution_imag_S;
//             solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

//         }

//         solution = std::complex<double>{solution_real, solution_imag};

//     }

// }

// void Solver::HH_far(const double xVers_0, const double xVers_1, const double xVers_2,
//             const double y_0, const double y_1, const double y_2,
//             const double n_0, const double n_1, const double n_2,
//             const double coupling_parameter, const double wavenumber,
//             std::complex<double>& solution)
// {

//     double solution_real, solution_imag;

//     const double inner_prod = - wavenumber * (xVers_0 * y_0 + xVers_1 * y_1 + xVers_2 * y_2);
//     const double inner_prod_2 = - wavenumber * (xVers_0 * n_0 + xVers_1 * n_1 + xVers_2 * n_2);

//     const double S_real = std::cos(inner_prod);
//     const double S_imag = std::sin(inner_prod);
//     const double D_real = -inner_prod_2 * S_imag;
//     const double D_imag = inner_prod_2 * S_real;

//     if (EQUATION_FORMULATION == 1) {

//         solution_real = S_real;
//         solution_imag = S_imag;

//     } else if (EQUATION_FORMULATION == 2) {

//         solution_real = D_real;
//         solution_imag = D_imag;

//     } else {

//         // EQUATION_FORMULATION == 3

//         solution_real = D_real + coupling_parameter * S_imag;
//         solution_imag = D_imag - coupling_parameter * S_real;

//     } 

//     solution = std::complex<double>{solution_real, solution_imag};

// }




// void Solver::compute_parallel_parameters() 
// {

//     int flag;

//     MPI_Initialized(&flag);

//     if (!flag) {
//         throw std::logic_error("Cannot use this code without MPI initialization\n");
//     }

//     MPI_Comm_size(mpi_comm_, &comm_size_);
//     MPI_Comm_rank(mpi_comm_, &comm_rank_);

//     split_points_ = std::vector<long long>(comm_size_ + 1);
//     split_points_2_ = std::vector<long long>(comm_size_ + 1);

//     const long long patches_per_rank = (Q_ * Qx_*Qy_) / comm_size_;
//     long long remaining_patches = (Q_ * Qx_*Qy_) % comm_size_;

//     split_points_[0] = 0;
//     split_points_[comm_size_] = Q_ * Qx_*Qy_;

//     split_points_2_[0] = 0;
//     split_points_2_[comm_size_] = Q_ * Qx_*Qy_ * Nu_int_*Nv_int_;

//     for (int i = 1; i < comm_size_; i++) {

//         split_points_[i] = split_points_[i-1] + patches_per_rank;

//         if (remaining_patches > 0) {
//             split_points_[i]++;
//             remaining_patches--;
//         }

//         split_points_2_[i] = split_points_[i] * Nu_int_*Nv_int_;

//     }

//     recv_counts_ = std::vector<int>(comm_size_);
//     displs_ = std::vector<int>(comm_size_, 0);

//     recv_counts_2_ = std::vector<int>(comm_size_);
//     displs_2_ = std::vector<int>(comm_size_, 0);

//     for (int rank = 0; rank < comm_size_; rank++) {

//         recv_counts_[rank] = split_points_[rank + 1] - split_points_[rank];
//         recv_counts_2_[rank] = split_points_2_[rank + 1] - split_points_2_[rank];

//         if (rank != 0) {

//             displs_[rank] = displs_[rank - 1] + recv_counts_[rank - 1];
//             displs_2_[rank] = displs_2_[rank - 1] + recv_counts_2_[rank - 1];

//         }

//     }

//     patch_low_ = split_points_[comm_rank_];
//     patch_up_ = split_points_[comm_rank_ + 1];

//     orig_patch_low_ = patch_low_ / (Qx_*Qy_);
//     orig_patch_up_ = (patch_up_ - 1) / (Qx_*Qy_);

//     point_low_ = patch_low_ * Nu_int_*Nv_int_;
//     point_up_ = patch_up_ * Nu_int_*Nv_int_;

//     omp_set_num_threads(NTHREADS);       

// }

// void Solver::load_interpolated_surface() 
// {
    
//     #pragma omp parallel for
//     for (long long i = orig_patch_low_; i <= orig_patch_up_; i++) {

      

//         std::ifstream fin(DIRECTORY + FILE_NAME + std::to_string(i+1) + ".txt");

//         InterpPatch elem;
        
//         fin >> elem.imax >> elem.jmax >> elem.kmax >> elem.zoneID;

//         elem.uNodes = std::vector<double>(elem.imax);

//         for (int j = 0; j < elem.imax; j++) {

//             fin >> elem.uNodes[j];

//         }

//         elem.vNodes = std::vector<double>(elem.jmax);

//         for (int j = 0; j < elem.jmax; j++) {

//             fin >> elem.vNodes[j];

//         }

//         elem.uWeights = barycentric_weights(elem.uNodes);

//         elem.vWeights = barycentric_weights(elem.vNodes);

//         elem.x = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (int j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.x[j];

//         }

//         elem.y = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (int j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.y[j];

//         }

//         elem.z = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.z[j];

//         }

//         elem.mask = std::vector<int>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.mask[j];

//         }

//         elem.dxdu = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.dxdu[j];

//         }

//         elem.dydu = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.dydu[j];

//         }

//         elem.dzdu = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.dzdu[j];

//         }

//         elem.dxdv = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.dxdv[j];

//         }

//         elem.dydv = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.dydv[j];

//         }

//         elem.dzdv = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.dzdv[j];

//         }

//         elem.dS = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.dS[j];

//         }

//         elem.nuX = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.nuX[j];

//         }

//         elem.nuY = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.nuY[j];

//         }

//         elem.nuZ = std::vector<double>(elem.imax * elem.jmax * elem.kmax);

//         for (long long j = 0; j < elem.imax * elem.jmax * elem.kmax; j++) {

//             fin >> elem.nuZ[j];

//         }

//         interp_surface_[i] = elem;

//                   // debug 
//         std::cout << "Rank " << comm_rank_ << " loading patch " << i << std::endl;

//         fin.close();

//     }

// }

// void Solver::compute_fejer_nodes_and_weights() 
// {

//     fejer_nodes_u_int_ = std::vector<double>(Nu_int_);
//     fejer_weights_u_int_ = std::vector<double>(Nu_int_);

//     fejer_nodes_v_int_ = std::vector<double>(Nv_int_);
//     fejer_weights_v_int_ = std::vector<double>(Nv_int_);

//     fejer_nodes_u_prec_ = std::vector<double>(Nu_prec_);
//     fejer_weights_u_prec_ = std::vector<double>(Nu_prec_);

//     fejer_nodes_v_prec_ = std::vector<double>(Nv_prec_);
//     fejer_weights_v_prec_ = std::vector<double>(Nv_prec_);

//     fejerquadrature1(fejer_nodes_u_int_, fejer_weights_u_int_, Nu_int_);
//     fejerquadrature1(fejer_nodes_v_int_, fejer_weights_v_int_, Nv_int_);
//     fejerquadrature1(fejer_nodes_u_prec_, fejer_weights_u_prec_, Nu_prec_);
//     fejerquadrature1(fejer_nodes_v_prec_, fejer_weights_v_prec_, Nv_prec_);

// }

// void Solver::compute_chebyshev_evaluations()
// {

//     Tn_ = std::vector<std::complex<double>>(Nu_int_*Nu_int_);
//     Tm_ = std::vector<std::complex<double>>(Nv_int_*Nv_int_);

//     cheb_evals(fejer_nodes_u_int_, Tn_, Nu_int_, Nu_int_);
//     cheb_evals(fejer_nodes_v_int_, Tm_, Nv_int_, Nv_int_); 

// }
    
// void Solver::compute_flags_domain()
// {

//     std::vector<int> flags_domain_u(patch_up_ - patch_low_);
//     std::vector<int> flags_domain_v(patch_up_ - patch_low_);
    
//     #pragma omp parallel for
//     for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

//         const long long q = patch_num / (Qx_*Qy_);
//         const int q_x = (patch_num % (Qx_*Qy_)) / Qy_;
//         const int q_y = (patch_num % (Qx_*Qy_)) % Qy_;

//         const long long pos = patch_num - patch_low_;

//         flags_domain_u[pos] = 0; 
//         flags_domain_v[pos] = 0; 

//         const double u_a = -1.0 + q_x * 2.0 / Qx_;
//         const double u_b = -1.0 + (q_x + 1) * 2.0 / Qx_;
//         const double v_a = -1.0 + q_y * 2.0 / Qy_;
//         const double v_b = -1.0 + (q_y + 1) * 2.0 / Qy_;

//         if (u_a == -1.0 && EDGE_FLAG_U_A[q]) {
//             flags_domain_u[pos] += 1; 
//         }

//         if (u_b == 1.0 && EDGE_FLAG_U_B[q]) {
//             flags_domain_u[pos] += 2;
//         }

//         if (v_a == -1.0 && EDGE_FLAG_V_A[q]) {
//             flags_domain_v[pos] += 1;
//         }

//         if (v_b == 1.0 && EDGE_FLAG_V_B[q]) {
//             flags_domain_v[pos] += 2; 
//         }

//     }

//     flags_domain_u_all_.resize(Q_ * Qx_*Qy_);
//     flags_domain_v_all_.resize(Q_ * Qx_*Qy_);

//     MPI_Allgatherv(&flags_domain_u[0], patch_up_-patch_low_, MPI_INT, &flags_domain_u_all_[0], &recv_counts_[0], &displs_[0], MPI_INT, mpi_comm_);
//     MPI_Allgatherv(&flags_domain_v[0], patch_up_-patch_low_, MPI_INT, &flags_domain_v_all_[0], &recv_counts_[0], &displs_[0], MPI_INT, mpi_comm_);

// }

// void Solver::compute_discretization_domain() 
// {

//     std::vector<double> disc_points_x(point_up_ - point_low_);
//     std::vector<double> disc_points_y(point_up_ - point_low_);
//     std::vector<double> disc_points_z(point_up_ - point_low_);

//     std::vector<double> norm_points_x(point_up_ - point_low_);
//     std::vector<double> norm_points_y(point_up_ - point_low_);
//     std::vector<double> norm_points_z(point_up_ - point_low_);

//     std::vector<double> dsdtjac(point_up_ - point_low_);
    
//     #pragma omp parallel for
//     for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

//         const long long q = patch_num / (Qx_*Qy_);
//         const int q_x = (patch_num % (Qx_*Qy_)) / Qy_;
//         const int q_y = (patch_num % (Qx_*Qy_)) % Qy_;

//         const long long pos = patch_num - patch_low_;

//         const int flag_u_loc = flags_domain_u_all_[patch_num];
//         const int flag_v_loc = flags_domain_v_all_[patch_num];

//         const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
//         const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
//         const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
//         const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

//         std::vector<double> vec_si(Nu_int_), dsi(Nu_int_);
//         std::vector<double> vec_tj(Nv_int_), dtj(Nv_int_);

//         for (int i = 0; i < Nu_int_; i++) {

//             double si = eta(fejer_nodes_u_int_[i], flag_u_loc);
//             si = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, si);
//             vec_si[i] = si;

//             dsi[i] = der_eta(fejer_nodes_u_int_[i], flag_u_loc) * (u_b_loc - u_a_loc) * 0.5;

//         }

//         for (int j = 0; j < Nv_int_; j++) {

//             double tj = eta(fejer_nodes_v_int_[j], flag_v_loc);                            
//             tj = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, tj);
//             vec_tj[j] = tj;

//             dtj[j] = der_eta(fejer_nodes_v_int_[j], flag_v_loc) * (v_b_loc - v_a_loc) * 0.5;

//         }

//         if (GEOMETRY == 0) {

//             for (int nu = 0; nu < Nu_int_; nu++) {
//                 for (int nv = 0; nv < Nv_int_; nv++) {

//                     const long long position = pos * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;

//                     const double si = vec_si[nu];
//                     const double tj = vec_tj[nv];

//                     parametrization_q(SPHERE_RADIUS, SPHERE_CENTER, si, tj, q, disc_points_x[position], disc_points_y[position], disc_points_z[position]);
//                     normal_q(SPHERE_RADIUS, si, tj, q, norm_points_x[position], norm_points_y[position], norm_points_z[position]);
//                     dsdtjac[position] = dsi[nu] * dtj[nv] * jacobian_q(SPHERE_RADIUS,si, tj, q);

//                 }
//             }

//         } else {             

//             InterpPatch interpolation_patch = interp_surface_[q];

//             lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
//                                         interpolation_patch.uWeights, interpolation_patch.vWeights,
//                                         interpolation_patch.x, vec_si, vec_tj, &disc_points_x[pos * Nu_int_*Nv_int_]);
//             lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
//                                         interpolation_patch.uWeights, interpolation_patch.vWeights,
//                                         interpolation_patch.y, vec_si, vec_tj, &disc_points_y[pos * Nu_int_*Nv_int_]);
//             lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
//                                         interpolation_patch.uWeights, interpolation_patch.vWeights,
//                                         interpolation_patch.z, vec_si, vec_tj, &disc_points_z[pos * Nu_int_*Nv_int_]);

//             lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
//                                         interpolation_patch.uWeights, interpolation_patch.vWeights,
//                                         interpolation_patch.nuX, vec_si, vec_tj, &norm_points_x[pos * Nu_int_*Nv_int_]);
//             lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
//                                         interpolation_patch.uWeights, interpolation_patch.vWeights,
//                                         interpolation_patch.nuY, vec_si, vec_tj, &norm_points_y[pos * Nu_int_*Nv_int_]);
//             lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
//                                         interpolation_patch.uWeights, interpolation_patch.vWeights,
//                                         interpolation_patch.nuZ, vec_si, vec_tj, &norm_points_z[pos * Nu_int_*Nv_int_]);

//             lagrange_interpolation_2D(interpolation_patch.uNodes, interpolation_patch.vNodes,
//                                         interpolation_patch.uWeights, interpolation_patch.vWeights,
//                                         interpolation_patch.dS, vec_si, vec_tj, &dsdtjac[pos * Nu_int_*Nv_int_]);

//             for (int nu = 0; nu < Nu_int_; nu++) {
//                 for (int nv = 0; nv < Nv_int_; nv++) {

//                     const long long position = pos * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;
                    
//                     dsdtjac[position] *= dsi[nu] * dtj[nv];

//                 }
//             }               

//         }  

//     }

//     disc_points_x_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
//     disc_points_y_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
//     disc_points_z_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

//     norm_points_x_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
//     norm_points_y_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
//     norm_points_z_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

//     dsdtjac_all_.resize(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

//     MPI_Allgatherv(&disc_points_x[0], point_up_-point_low_, MPI_DOUBLE, &disc_points_x_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&disc_points_y[0], point_up_-point_low_, MPI_DOUBLE, &disc_points_y_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&disc_points_z[0], point_up_-point_low_, MPI_DOUBLE, &disc_points_z_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);

//     MPI_Allgatherv(&norm_points_x[0], point_up_-point_low_, MPI_DOUBLE, &norm_points_x_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&norm_points_y[0], point_up_-point_low_, MPI_DOUBLE, &norm_points_y_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&norm_points_z[0], point_up_-point_low_, MPI_DOUBLE, &norm_points_z_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);

//     MPI_Allgatherv(&dsdtjac[0], point_up_-point_low_, MPI_DOUBLE, &dsdtjac_all_[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE, mpi_comm_);

// }            

// void Solver::compute_coupling_parameter() 
// {

//     double loc_max_x = std::numeric_limits<double>::lowest();
//     double loc_max_y = std::numeric_limits<double>::lowest();
//     double loc_max_z = std::numeric_limits<double>::lowest();

//     double loc_min_x = std::numeric_limits<double>::max();
//     double loc_min_y = std::numeric_limits<double>::max();
//     double loc_min_z = std::numeric_limits<double>::max();

//     for (long long i = point_low_; i < point_up_; i++) {

//         const double pointx = disc_points_x_all_[i];
//         const double pointy = disc_points_y_all_[i];
//         const double pointz = disc_points_z_all_[i];

//         if (loc_max_x < pointx) loc_max_x = pointx;
//         if (loc_max_y < pointy) loc_max_y = pointy;
//         if (loc_max_z < pointz) loc_max_z = pointz;
        
//         if (loc_min_x > pointx) loc_min_x = pointx;
//         if (loc_min_y > pointy) loc_min_y = pointy;
//         if (loc_min_z > pointz) loc_min_z = pointz;

//     }

//     double global_max_x;
//     double global_max_y;
//     double global_max_z;

//     double global_min_x;
//     double global_min_y;
//     double global_min_z;

//     MPI_Allreduce(&loc_max_x, &global_max_x, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
//     MPI_Allreduce(&loc_max_y, &global_max_y, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
//     MPI_Allreduce(&loc_max_z, &global_max_z, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
//     MPI_Allreduce(&loc_min_x, &global_min_x, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_);
//     MPI_Allreduce(&loc_min_y, &global_min_y, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_);
//     MPI_Allreduce(&loc_min_z, &global_min_z, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_);

//     const double scatDiam = std::max({global_max_x - global_min_x, global_max_y - global_min_y, global_max_z - global_min_z});
    
//     coupling_parameter_ = std::max(3.0, scatDiam / std::sqrt(lambda_.real()*lambda_.real() + lambda_.imag()*lambda_.imag()));

// }

// void Solver::compute_near_singular_patches_estimate()
// {

//     std::vector<double> x_min(patch_up_ - patch_low_, std::numeric_limits<double>::max());
//     std::vector<double> y_min(patch_up_ - patch_low_, std::numeric_limits<double>::max());
//     std::vector<double> z_min(patch_up_ - patch_low_, std::numeric_limits<double>::max());

//     std::vector<double> x_max(patch_up_ - patch_low_, std::numeric_limits<double>::lowest());
//     std::vector<double> y_max(patch_up_ - patch_low_, std::numeric_limits<double>::lowest());
//     std::vector<double> z_max(patch_up_ - patch_low_, std::numeric_limits<double>::lowest());
    
//     #pragma omp parallel for
//     for (long long patch_num = 0; patch_num < patch_up_ - patch_low_; patch_num++) {

//         for (int nu = 0; nu < Nu_int_; nu++) {
//             for (int nv = 0; nv < Nv_int_; nv++) {

//                 const long long i = (patch_num + patch_low_) * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;               

//                 const double pointx = disc_points_x_all_[i];
//                 const double pointy = disc_points_y_all_[i];
//                 const double pointz = disc_points_z_all_[i];

//                 if (x_min[patch_num] > pointx) x_min[patch_num] = pointx;
//                 if (y_min[patch_num] > pointy) y_min[patch_num] = pointy;
//                 if (z_min[patch_num] > pointz) z_min[patch_num] = pointz;
                
//                 if (x_max[patch_num] < pointx) x_max[patch_num] = pointx;
//                 if (y_max[patch_num] < pointy) y_max[patch_num] = pointy;
//                 if (z_max[patch_num] < pointz) z_max[patch_num] = pointz;  

//             }
//         }                

//         const double Lx = x_max[patch_num] - x_min[patch_num];
//         const double Ly = y_max[patch_num] - y_min[patch_num];
//         const double Lz = z_max[patch_num] - z_min[patch_num];

//         const double Lmax = std::max({Lx, Ly, Lz});
//         const double factor = Lmax * 2.0E-2; 

//         x_min[patch_num] -= factor;
//         x_max[patch_num] += factor;
//         y_min[patch_num] -= factor;
//         y_max[patch_num] += factor;
//         z_min[patch_num] -= factor;
//         z_max[patch_num] += factor;

//     }    

//     MPI_Barrier(mpi_comm_);        

//     std::vector<double> x_min_all(Q_ * Qx_*Qy_);
//     std::vector<double> y_min_all(Q_ * Qx_*Qy_);
//     std::vector<double> z_min_all(Q_ * Qx_*Qy_);

//     std::vector<double> x_max_all(Q_ * Qx_*Qy_);
//     std::vector<double> y_max_all(Q_ * Qx_*Qy_);
//     std::vector<double> z_max_all(Q_ * Qx_*Qy_);

//     MPI_Allgatherv(&x_min[0], patch_up_-patch_low_, MPI_DOUBLE, &x_min_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&y_min[0], patch_up_-patch_low_, MPI_DOUBLE, &y_min_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&z_min[0], patch_up_-patch_low_, MPI_DOUBLE, &z_min_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);

//     MPI_Allgatherv(&x_max[0], patch_up_-patch_low_, MPI_DOUBLE, &x_max_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&y_max[0], patch_up_-patch_low_, MPI_DOUBLE, &y_max_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&z_max[0], patch_up_-patch_low_, MPI_DOUBLE, &z_max_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);

//     start_sing_and_near_sing_patches_estimate_ = std::vector<long long>(patch_up_ - patch_low_);
//     size_sing_and_near_sing_patches_estimate_ = std::vector<long long>(patch_up_ - patch_low_);

//     std::vector<std::unordered_set<long long>> patches_not_in_rank(comm_size_);

//     long long index = 0;
    
//     for (long long patch_num = 0; patch_num < patch_up_ - patch_low_; patch_num++) {

//         long long total = 0;
        
//         for (long long patch_num_2 = 0; patch_num_2 < Q_ * Qx_*Qy_; patch_num_2++) {

//             if ((x_max[patch_num] >= x_min_all[patch_num_2]) && (x_max_all[patch_num_2] >= x_min[patch_num]) &&
//                 (y_max[patch_num] >= y_min_all[patch_num_2]) && (y_max_all[patch_num_2] >= y_min[patch_num]) &&
//                 (z_max[patch_num] >= z_min_all[patch_num_2]) && (z_max_all[patch_num_2] >= z_min[patch_num])) {

//                 sing_and_near_sing_patches_estimate_.push_back(patch_num_2);
//                 total++;

//                 if (GEOMETRY != 0) {
                
//                     const long long q = patch_num_2 / (Qx_ * Qy_);

//                     if (interp_surface_.count(q) == 0) {
                    
//                         int rank;

//                         for (int i = 0; i < comm_size_+1; i++) {

//                             if ((patch_num_2 >= split_points_[i]) && (patch_num_2 < split_points_[i+1])) {

//                                 rank = i;
//                                 break;

//                             }

//                         }
                    
//                         patches_not_in_rank[rank].insert(q); 

//                     }               

//                 }

//             }

//         }

//         start_sing_and_near_sing_patches_estimate_[patch_num] = index;
//         size_sing_and_near_sing_patches_estimate_[patch_num] = total;
//         index += total;

//     }

//     // This trick frees the memory of the vectors
//     std::vector<double>().swap(x_min);
//     std::vector<double>().swap(y_min);
//     std::vector<double>().swap(z_min);

//     std::vector<double>().swap(x_min_all);
//     std::vector<double>().swap(y_min_all);
//     std::vector<double>().swap(z_min_all);

//     std::vector<double>().swap(x_max);
//     std::vector<double>().swap(y_max);
//     std::vector<double>().swap(z_max);

//     std::vector<double>().swap(x_max_all);
//     std::vector<double>().swap(y_max_all);
//     std::vector<double>().swap(z_max_all);

//     MPI_Barrier(mpi_comm_);
    
//     if (GEOMETRY != 0) {

//         for (int rank = 0; rank < comm_size_; rank++) {

//             const long long size = patches_not_in_rank[rank].size();
//             const std::vector<long long> elements(patches_not_in_rank[rank].begin(), patches_not_in_rank[rank].end());

//             std::vector<long long> size_all(comm_size_);

//             MPI_Gather(&size, 1, MPI_LONG_LONG, &size_all[0], 1, MPI_LONG_LONG, rank, mpi_comm_);

//             std::vector<long long> patches;
            
//             std::vector<int> recv_counts(comm_size_);
//             std::vector<int> displs(comm_size_, 0);

//             if (comm_rank_ == rank) {

//                 long long total_elements = 0;

//                 for (int i = 0; i < comm_size_; i++) {

//                     recv_counts[i] = size_all[i];

//                     if (i != 0) {

//                         displs[i] = displs[i-1] + recv_counts[i-1];

//                     }

//                     total_elements += size_all[i];

//                 }

//                 patches.resize(total_elements);

//             }

//             MPI_Gatherv(&elements[0], size, MPI_LONG_LONG, &patches[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, rank, mpi_comm_);

//             std::vector<int> vector_imax;
//             std::vector<int> vector_jmax;
//             std::vector<int> vector_kmax;

//             std::vector<long long> vector_zoneID;

//             std::vector<double> vector_uNodes;
//             std::vector<double> vector_vNodes;

//             std::vector<double> vector_uWeights;
//             std::vector<double> vector_vWeights;

//             std::vector<double> vector_x, vector_y, vector_z;
//             std::vector<int> vector_mask;
//             std::vector<double> vector_dxdu, vector_dydu, vector_dzdu;
//             std::vector<double> vector_dxdv, vector_dydv, vector_dzdv;
//             std::vector<double> vector_dS;
//             std::vector<double> vector_nuX, vector_nuY, vector_nuZ;

//             if (comm_rank_ == rank) {
                
//                 for (long long i = 0; i < patches.size(); i++) {

//                     InterpPatch elem = interp_surface_[patches[i]];

//                     vector_imax.push_back(elem.imax);
//                     vector_jmax.push_back(elem.jmax);
//                     vector_kmax.push_back(elem.kmax);

//                     vector_zoneID.push_back(elem.zoneID);

//                     vector_uNodes.insert(vector_uNodes.end(), elem.uNodes.begin(), elem.uNodes.end());
//                     vector_vNodes.insert(vector_vNodes.end(), elem.vNodes.begin(), elem.vNodes.end());

//                     vector_uWeights.insert(vector_uWeights.end(), elem.uWeights.begin(), elem.uWeights.end());
//                     vector_vWeights.insert(vector_vWeights.end(), elem.vWeights.begin(), elem.vWeights.end());

//                     vector_x.insert(vector_x.end(), elem.x.begin(), elem.x.end());
//                     vector_y.insert(vector_y.end(), elem.y.begin(), elem.y.end());
//                     vector_z.insert(vector_z.end(), elem.z.begin(), elem.z.end());

//                     vector_mask.insert(vector_mask.end(), elem.mask.begin(), elem.mask.end());

//                     vector_dxdu.insert(vector_dxdu.end(), elem.dxdu.begin(), elem.dxdu.end());
//                     vector_dydu.insert(vector_dydu.end(), elem.dydu.begin(), elem.dydu.end());
//                     vector_dzdu.insert(vector_dzdu.end(), elem.dzdu.begin(), elem.dzdu.end());

//                     vector_dxdv.insert(vector_dxdv.end(), elem.dxdv.begin(), elem.dxdv.end());
//                     vector_dydv.insert(vector_dydv.end(), elem.dydv.begin(), elem.dydv.end());
//                     vector_dzdv.insert(vector_dzdv.end(), elem.dzdv.begin(), elem.dzdv.end());

//                     vector_dS.insert(vector_dS.end(), elem.dS.begin(), elem.dS.end());

//                     vector_nuX.insert(vector_nuX.end(), elem.nuX.begin(), elem.nuX.end());
//                     vector_nuY.insert(vector_nuY.end(), elem.nuY.begin(), elem.nuY.end());
//                     vector_nuZ.insert(vector_nuZ.end(), elem.nuZ.begin(), elem.nuZ.end());

//                 }

//             }
            
//             std::vector<int> vector_imax_loc(size);
//             std::vector<int> vector_jmax_loc(size);
//             std::vector<int> vector_kmax_loc(size);

//             std::vector<long long> vector_zoneID_loc(size);

//             MPI_Scatterv(&vector_imax[0], &recv_counts[0], &displs[0], MPI_INT, &vector_imax_loc[0], size, MPI_INT, rank, mpi_comm_);
//             MPI_Scatterv(&vector_jmax[0], &recv_counts[0], &displs[0], MPI_INT, &vector_jmax_loc[0], size, MPI_INT, rank, mpi_comm_);
//             MPI_Scatterv(&vector_kmax[0], &recv_counts[0], &displs[0], MPI_INT, &vector_kmax_loc[0], size, MPI_INT, rank, mpi_comm_);

//             MPI_Scatterv(&vector_zoneID[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, &vector_zoneID_loc[0], size, MPI_LONG_LONG, rank, mpi_comm_);

//             int imax_total = 0;
//             int jmax_total = 0;
//             int imaxjmaxkmax_total = 0;

//             for (long long i = 0; i < size; i++) {

//                 imax_total += vector_imax_loc[i];
//                 jmax_total += vector_jmax_loc[i];
//                 imaxjmaxkmax_total += vector_imax_loc[i] * vector_jmax_loc[i] * vector_kmax_loc[i];

//             }

//             std::vector<int> recv_counts_imax(comm_size_);
//             std::vector<int> recv_counts_jmax(comm_size_);
//             std::vector<int> recv_counts_imaxjmaxkmax(comm_size_);

//             MPI_Gather(&imax_total, 1, MPI_INT, &recv_counts_imax[0], 1, MPI_INT, rank, mpi_comm_);
//             MPI_Gather(&jmax_total, 1, MPI_INT, &recv_counts_jmax[0], 1, MPI_INT, rank, mpi_comm_);
//             MPI_Gather(&imaxjmaxkmax_total, 1, MPI_INT, &recv_counts_imaxjmaxkmax[0], 1, MPI_INT, rank, mpi_comm_);                

//             std::vector<int> displs_imax(comm_size_, 0);
//             std::vector<int> displs_jmax(comm_size_, 0);
//             std::vector<int> displs_imaxjmaxkmax(comm_size_, 0);

//             if (comm_rank_ == rank) {

//                 for (int i = 1; i < comm_size_; i++) {

//                     displs_imax[i] = displs_imax[i-1] + recv_counts_imax[i-1];
//                     displs_jmax[i] = displs_jmax[i-1] + recv_counts_jmax[i-1];
//                     displs_imaxjmaxkmax[i] = displs_imaxjmaxkmax[i-1] + recv_counts_imaxjmaxkmax[i-1];

//                 }

//             }

//             std::vector<double> vector_uNodes_loc(imax_total);
//             std::vector<double> vector_vNodes_loc(jmax_total);

//             std::vector<double> vector_uWeights_loc(imax_total);
//             std::vector<double> vector_vWeights_loc(jmax_total);

//             std::vector<double> vector_x_loc(imaxjmaxkmax_total), vector_y_loc(imaxjmaxkmax_total), vector_z_loc(imaxjmaxkmax_total);
//             std::vector<int> vector_mask_loc(imaxjmaxkmax_total);
//             std::vector<double> vector_dxdu_loc(imaxjmaxkmax_total), vector_dydu_loc(imaxjmaxkmax_total), vector_dzdu_loc(imaxjmaxkmax_total);
//             std::vector<double> vector_dxdv_loc(imaxjmaxkmax_total), vector_dydv_loc(imaxjmaxkmax_total), vector_dzdv_loc(imaxjmaxkmax_total);
//             std::vector<double> vector_dS_loc(imaxjmaxkmax_total);
//             std::vector<double> vector_nuX_loc(imaxjmaxkmax_total), vector_nuY_loc(imaxjmaxkmax_total), vector_nuZ_loc(imaxjmaxkmax_total);

//             MPI_Scatterv(&vector_uNodes[0], &recv_counts_imax[0], &displs_imax[0], MPI_DOUBLE, &vector_uNodes_loc[0], imax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_vNodes[0], &recv_counts_jmax[0], &displs_jmax[0], MPI_DOUBLE, &vector_vNodes_loc[0], jmax_total, MPI_DOUBLE, rank, mpi_comm_);

//             MPI_Scatterv(&vector_uWeights[0], &recv_counts_imax[0], &displs_imax[0], MPI_DOUBLE, &vector_uWeights_loc[0], imax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_vWeights[0], &recv_counts_jmax[0], &displs_jmax[0], MPI_DOUBLE, &vector_vWeights_loc[0], jmax_total, MPI_DOUBLE, rank, mpi_comm_);

//             MPI_Scatterv(&vector_x[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_x_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_y[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_y_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_z[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_z_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            
//             MPI_Scatterv(&vector_mask[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_INT, &vector_mask_loc[0], imaxjmaxkmax_total, MPI_INT, rank, mpi_comm_);
            
//             MPI_Scatterv(&vector_dxdu[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dxdu_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_dydu[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dydu_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_dzdu[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dzdu_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
            
//             MPI_Scatterv(&vector_dxdv[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dxdv_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_dydv[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dydv_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_dzdv[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dzdv_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);

//             MPI_Scatterv(&vector_dS[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_dS_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);

//             MPI_Scatterv(&vector_nuX[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_nuX_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_nuY[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_nuY_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);
//             MPI_Scatterv(&vector_nuZ[0], &recv_counts_imaxjmaxkmax[0], &displs_imaxjmaxkmax[0], MPI_DOUBLE, &vector_nuZ_loc[0], imaxjmaxkmax_total, MPI_DOUBLE, rank, mpi_comm_);

//             long long start_imax = 0;
//             long long start_jmax = 0;
//             long long start_imaxjmaxkmax = 0;

//             for (long long i = 0; i < size; i++) {

//                 InterpPatch new_elem;

//                 new_elem.imax = vector_imax_loc[i];
//                 new_elem.jmax = vector_jmax_loc[i];
//                 new_elem.kmax = vector_kmax_loc[i];

//                 new_elem.zoneID = vector_zoneID_loc[i];

//                 new_elem.uNodes.insert(new_elem.uNodes.end(), vector_uNodes_loc.begin() + start_imax, vector_uNodes_loc.begin() + start_imax + new_elem.imax);
//                 new_elem.vNodes.insert(new_elem.vNodes.end(), vector_vNodes_loc.begin() + start_jmax, vector_vNodes_loc.begin() + start_jmax + new_elem.jmax);

//                 new_elem.uWeights.insert(new_elem.uWeights.end(), vector_uWeights_loc.begin() + start_imax, vector_uWeights_loc.begin() + start_imax + new_elem.imax);
//                 new_elem.vWeights.insert(new_elem.vWeights.end(), vector_vWeights_loc.begin() + start_jmax, vector_vWeights_loc.begin() + start_jmax + new_elem.jmax);

//                 new_elem.x.insert(new_elem.x.end(), vector_x_loc.begin() + start_imaxjmaxkmax, vector_x_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
//                 new_elem.y.insert(new_elem.y.end(), vector_y_loc.begin() + start_imaxjmaxkmax, vector_y_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
//                 new_elem.z.insert(new_elem.z.end(), vector_z_loc.begin() + start_imaxjmaxkmax, vector_z_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

//                 new_elem.mask.insert(new_elem.mask.end(), vector_mask_loc.begin() + start_imaxjmaxkmax, vector_mask_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

//                 new_elem.dxdu.insert(new_elem.dxdu.end(), vector_dxdu_loc.begin() + start_imaxjmaxkmax, vector_dxdu_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
//                 new_elem.dydu.insert(new_elem.dydu.end(), vector_dydu_loc.begin() + start_imaxjmaxkmax, vector_dydu_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
//                 new_elem.dzdu.insert(new_elem.dzdu.end(), vector_dzdu_loc.begin() + start_imaxjmaxkmax, vector_dzdu_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

//                 new_elem.dxdv.insert(new_elem.dxdv.end(), vector_dxdv_loc.begin() + start_imaxjmaxkmax, vector_dxdv_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
//                 new_elem.dydv.insert(new_elem.dydv.end(), vector_dydv_loc.begin() + start_imaxjmaxkmax, vector_dydv_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
//                 new_elem.dzdv.insert(new_elem.dzdv.end(), vector_dzdv_loc.begin() + start_imaxjmaxkmax, vector_dzdv_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

//                 new_elem.dS.insert(new_elem.dS.end(), vector_dS_loc.begin() + start_imaxjmaxkmax, vector_dS_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

//                 new_elem.nuX.insert(new_elem.nuX.end(), vector_nuX_loc.begin() + start_imaxjmaxkmax, vector_nuX_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
//                 new_elem.nuY.insert(new_elem.nuY.end(), vector_nuY_loc.begin() + start_imaxjmaxkmax, vector_nuY_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);
//                 new_elem.nuZ.insert(new_elem.nuZ.end(), vector_nuZ_loc.begin() + start_imaxjmaxkmax, vector_nuZ_loc.begin() + start_imaxjmaxkmax + new_elem.imax*new_elem.jmax*new_elem.kmax);

//                 interp_surface_[elements[i]] = new_elem;

//                 start_imax += new_elem.imax;
//                 start_jmax += new_elem.jmax;
//                 start_imaxjmaxkmax += new_elem.imax*new_elem.jmax*new_elem.kmax;

//             }

//         }

//     }
            
// }

// void Solver::beta(const double r_0, const double r_1, const double r_2,
//             const long long q, const int flag_u_loc, const int flag_v_loc, 
//             const double u_a_loc, const double u_b_loc, const double v_a_loc, const double v_b_loc,
//             const double ubar_loc, const double vbar_loc,
//             std::vector<std::complex<double>>& prec)
// {

//     std::vector<double> ui(Nu_prec_), si(Nu_prec_), muwi(Nu_prec_);
//     std::vector<double> vj(Nv_prec_), tj(Nv_prec_), muwj(Nv_prec_);

//     for (int i = 0; i < Nu_prec_; i++) {

//         ui[i] = xi(ubar_loc, fejer_nodes_u_prec_[i]);
//         si[i] = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, eta(ui[i], flag_u_loc));
//         muwi[i] = der_xi(ubar_loc, fejer_nodes_u_prec_[i]) * fejer_weights_u_prec_[i];
        
//     }

//     for (int j = 0; j < Nv_prec_; j++) {

//         vj[j] = xi(vbar_loc, fejer_nodes_v_prec_[j]);
//         tj[j] = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, eta(vj[j], flag_v_loc));
//         muwj[j] = der_xi(vbar_loc, fejer_nodes_v_prec_[j]) * fejer_weights_v_prec_[j];
        
//     }

//     std::vector<std::complex<double>> H(Nu_prec_*Nv_prec_);
    
//     if (GEOMETRY == 0) {

//         for (int i = 0; i < Nu_prec_; i++) {
//             for (int j = 0; j < Nv_prec_; j++) {                        
                
//                 const double s = si[i];
//                 const double t = tj[j];

//                 double px, py, pz;
//                 double nx, ny, nz;
//                 double dxdsx, dxdsy, dxdsz;
//                 double dxdtx, dxdty, dxdtz;
//                 double dxdsdsx, dxdsdsy, dxdsdsz;
//                 double dxdsdtx, dxdsdty, dxdsdtz;
//                 double dxdtdtx, dxdtdty, dxdtdtz;

//                 parametrization_q(SPHERE_RADIUS,SPHERE_CENTER,s, t, q, px, py, pz);
//                 normal_q(SPHERE_RADIUS,s, t, q, nx, ny, nz);
//                 dxds_q(SPHERE_RADIUS,s, t, q, dxdsx, dxdsy, dxdsz);
//                 dxdt_q(SPHERE_RADIUS,s, t, q, dxdtx, dxdty, dxdtz);
//                 dxdsds_q(SPHERE_RADIUS,s, t, q, dxdsdsx, dxdsdsy, dxdsdsz);
//                 dxdsdt_q(SPHERE_RADIUS,s, t, q, dxdsdtx, dxdsdty, dxdsdtz);
//                 dxdtdt_q(SPHERE_RADIUS,s, t, q, dxdtdtx, dxdtdty, dxdtdtz);

//                 dxdsx *= 0.5 * (u_b_loc - u_a_loc);
//                 dxdsy *= 0.5 * (u_b_loc - u_a_loc);
//                 dxdsz *= 0.5 * (u_b_loc - u_a_loc);

//                 dxdtx *= 0.5 * (v_b_loc - v_a_loc);
//                 dxdty *= 0.5 * (v_b_loc - v_a_loc);
//                 dxdtz *= 0.5 * (v_b_loc - v_a_loc);

//                 dxdsdsx *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (u_b_loc - u_a_loc);
//                 dxdsdsy *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (u_b_loc - u_a_loc);
//                 dxdsdsz *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (u_b_loc - u_a_loc);

//                 dxdsdtx *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (v_b_loc - v_a_loc);
//                 dxdsdty *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (v_b_loc - v_a_loc);
//                 dxdsdtz *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (v_b_loc - v_a_loc);

//                 dxdtdtx *= 0.5 * (v_b_loc - v_a_loc) * 0.5 * (v_b_loc - v_a_loc);
//                 dxdtdty *= 0.5 * (v_b_loc - v_a_loc) * 0.5 * (v_b_loc - v_a_loc);
//                 dxdtdtz *= 0.5 * (v_b_loc - v_a_loc) * 0.5 * (v_b_loc - v_a_loc);

//                 HH(r_0, r_1, r_2,
//                 px, py, pz,
//                 nx, ny, nz,
//                 dxdsx, dxdsy, dxdsz,
//                 dxdtx, dxdty, dxdtz,
//                 dxdsdsx, dxdsdsy, dxdsdsz,
//                 dxdsdtx, dxdsdty, dxdsdtz,
//                 dxdtdtx, dxdtdty, dxdtdtz,
//                 coupling_parameter_,
//                 wavenumber_,
//                 H[i*Nv_prec_+j]);

//                 H[i*Nv_prec_+j] *= muwi[i]*muwj[j];

//             }
//         }

//     } else {

//         InterpPatch elem = interp_surface_[q];

//         std::vector<double> px(Nu_prec_*Nv_prec_), py(Nu_prec_*Nv_prec_), pz(Nu_prec_*Nv_prec_);
//         std::vector<double> nx(Nu_prec_*Nv_prec_), ny(Nu_prec_*Nv_prec_), nz(Nu_prec_*Nv_prec_);
//         std::vector<double> dxdsx(Nu_prec_*Nv_prec_), dxdsy(Nu_prec_*Nv_prec_), dxdsz(Nu_prec_*Nv_prec_);
//         std::vector<double> dxdtx(Nu_prec_*Nv_prec_), dxdty(Nu_prec_*Nv_prec_), dxdtz(Nu_prec_*Nv_prec_);

//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.x, si, tj, &px[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.y, si, tj, &py[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.z, si, tj, &pz[0]);

//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.nuX, si, tj, &nx[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.nuY, si, tj, &ny[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.nuZ, si, tj, &nz[0]);

//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dxdu, si, tj, &dxdsx[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dydu, si, tj, &dxdsy[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dzdu, si, tj, &dxdsz[0]);

//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dxdv, si, tj, &dxdtx[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dydv, si, tj, &dxdty[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dzdv, si, tj, &dxdtz[0]);

//         for (int i = 0; i < Nu_prec_; i++) {
//             for (int j = 0; j < Nv_prec_; j++) {

//                 HH(r_0, r_1, r_2,
//                 px[i*Nv_prec_+j], py[i*Nv_prec_+j], pz[i*Nv_prec_+j],
//                 nx[i*Nv_prec_+j], ny[i*Nv_prec_+j], nz[i*Nv_prec_+j],
//                 dxdsx[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc), dxdsy[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc), dxdsz[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc),
//                 dxdtx[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc), dxdty[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc), dxdtz[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc),
//                 0.0, 0.0, 0.0,
//                 0.0, 0.0, 0.0,
//                 0.0, 0.0, 0.0,
//                 coupling_parameter_,
//                 wavenumber_,
//                 H[i*Nv_prec_+j]);

//                 H[i*Nv_prec_+j] *= muwi[i]*muwj[j];

//             }
//         }

//     }

//     std::vector<std::complex<double>> Tn_mat(Nu_prec_*Nu_int_), Tm_mat(Nv_prec_*Nv_int_);
    
//     cheb_evals(ui, Tn_mat, Nu_prec_, Nu_int_);
//     cheb_evals(vj, Tm_mat, Nv_prec_, Nv_int_);

//     std::vector<std::complex<double>> matprod1(Nv_prec_*Nu_int_, {0.0, 0.0});

//     std::complex<double> one(1.0, 0.0);
//     std::complex<double> zero(0.0, 0.0);

//     cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//                 Nu_int_, Nv_prec_, Nu_prec_, 
//                 &one, &Tn_mat[0], Nu_prec_, &H[0], Nv_prec_, &zero, &matprod1[0], Nv_prec_);

//     cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
//                 Nu_int_, Nv_int_, Nv_prec_,
//                 &one, &matprod1[0], Nv_prec_, &Tm_mat[0], Nv_prec_, &zero, &prec[0], Nv_int_);
    
// }
    
// void Solver::compute_precomputations()
// {
    
//     std::vector<double> vec_mean_size;

//     if (DELTA_METHOD == 2) {

//         vec_mean_size.resize(patch_up_ - patch_low_);

//         #pragma omp parallel for
//         for (long long npatch = patch_low_; npatch < patch_up_; npatch++) {

//             std::vector<double> min_vec(3, std::numeric_limits<double>::max());
//             std::vector<double> max_vec(3, std::numeric_limits<double>::lowest());

//             for (int ii = 0; ii < Nu_int_; ii++) {
//                 for (int jj = 0; jj < Nv_int_; jj++) {

//                     const long long position = npatch * Nu_int_*Nv_int_ + ii * Nv_int_ + jj;

//                     const double pointx = disc_points_x_all_[position];
//                     const double pointy = disc_points_y_all_[position];
//                     const double pointz = disc_points_z_all_[position];

//                     if (min_vec[0] > pointx) min_vec[0] = pointx;
//                     if (min_vec[1] > pointy) min_vec[1] = pointy;
//                     if (min_vec[2] > pointz) min_vec[2] = pointz;
                    
//                     if (max_vec[0] < pointx) max_vec[0] = pointx;
//                     if (max_vec[1] < pointy) max_vec[1] = pointy;
//                     if (max_vec[2] < pointz) max_vec[2] = pointz;  

//                 }
//             }

//             vec_mean_size[npatch - patch_low_] = ((max_vec[0] - min_vec[0]) +
//                                                     (max_vec[1] - min_vec[1]) + 
//                                                     (max_vec[2] - min_vec[2])) / 3.0;

//         }

//     }

//     #pragma omp parallel
//     {

//     std::vector<std::complex<double>> precomputations_thread;
//     std::vector<long long> point_precomputations_thread;
//     std::vector<long long> patch_num_precomputations_thread; 

//     precomputations_thread.reserve(Nu_int_*Nv_int_ * (patch_up_-patch_low_) * 9 * Nu_int_*Nv_int_ / NTHREADS);
//     point_precomputations_thread.reserve(Nu_int_*Nv_int_ * (patch_up_-patch_low_) * 9 / NTHREADS);
//     patch_num_precomputations_thread.reserve(Nu_int_*Nv_int_ * (patch_up_-patch_low_) * 9 / NTHREADS);
    
//     #pragma omp for
//     for (long long npoint = point_low_; npoint < point_up_; npoint++) {

//         long long idx_point;
//         idx_point = npoint;

//         const double r_0 = disc_points_x_all_[idx_point];
//         const double r_1 = disc_points_y_all_[idx_point];
//         const double r_2 = disc_points_z_all_[idx_point];

//         const long long patch_num_own = npoint / (Nu_int_*Nv_int_);

//         const long long start = start_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];
//         const long long total = size_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];

//         for (long long i = 0; i < total; i++) {

//             const long long patch_num = sing_and_near_sing_patches_estimate_[start + i];

//             double min_dist_loc;
//             int indx_point_min_dist_loc;

//             bool sing_near_sing_flag = false;

//             if (patch_num_own == patch_num) {

//                 min_dist_loc = 0.0;
//                 indx_point_min_dist_loc = npoint % (Nu_int_*Nv_int_);  

//                 sing_near_sing_flag = true; 

//             } else {  

//                 min_dist_loc = std::numeric_limits<double>::max();

//                 for (int ii = 0; ii < Nu_int_; ii++) {
//                     for (int jj = 0; jj < Nv_int_; jj++) {

//                         double x_loc;
//                         double y_loc;
//                         double z_loc;

//                         const long long position = patch_num * Nu_int_*Nv_int_ + ii * Nv_int_ + jj;

//                         x_loc = disc_points_x_all_[position];
//                         y_loc = disc_points_y_all_[position];
//                         z_loc = disc_points_z_all_[position];

//                         const double diff_0 = x_loc - r_0;
//                         const double diff_1 = y_loc - r_1;
//                         const double diff_2 = z_loc - r_2;

//                         const double dist_aux = diff_0*diff_0 + diff_1*diff_1 + diff_2*diff_2;

//                         if (dist_aux < min_dist_loc) {

//                             min_dist_loc = dist_aux;
//                             indx_point_min_dist_loc = ii * Nv_int_ + jj;

//                         }

//                         if (!sing_near_sing_flag && DELTA_METHOD == 1) {

//                             const bool proximity = (std::abs(std::floor(r_0/PROXIMITY_BOX_SIZE) - std::floor(x_loc/PROXIMITY_BOX_SIZE)) <= 1) &&
//                                                     (std::abs(std::floor(r_1/PROXIMITY_BOX_SIZE) - std::floor(y_loc/PROXIMITY_BOX_SIZE)) <= 1) &&
//                                                     (std::abs(std::floor(r_2/PROXIMITY_BOX_SIZE) - std::floor(z_loc/PROXIMITY_BOX_SIZE)) <= 1);

//                             if (proximity) {

//                                 sing_near_sing_flag = true;

//                             }

//                         }

//                     }
//                 }

//                 min_dist_loc = std::sqrt(min_dist_loc);

//             }

//             if ((DELTA_METHOD == 2) && (min_dist_loc <= vec_mean_size[patch_num_own - patch_low_] * PERCENT_BOX_SIZE)) {

//                 sing_near_sing_flag = true;

//             }     
        
//             if (sing_near_sing_flag) {  

//                 const int nu_loc = indx_point_min_dist_loc / Nv_int_;
//                 const int nv_loc = indx_point_min_dist_loc % Nv_int_;  

//                 int flag_u_loc;
//                 int flag_v_loc;

//                 flag_u_loc = flags_domain_u_all_[patch_num];
//                 flag_v_loc = flags_domain_v_all_[patch_num];

//                 const long long q = patch_num / (Qx_*Qy_);
//                 const int q_x = (patch_num / Qy_) % Qx_;
//                 const int q_y = patch_num % Qy_;

//                 const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
//                 const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
//                 const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
//                 const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_; 

//                 double ubar_loc, vbar_loc;
                        
//                 if (min_dist_loc < 1E-12) {
                    
//                     ubar_loc = eta(fejer_nodes_u_int_[nu_loc], flag_u_loc);
//                     vbar_loc = eta(fejer_nodes_v_int_[nv_loc], flag_v_loc);

//                 } else {

//                     double u_a_min, u_b_min, v_a_min, v_b_min;

//                     if (nu_loc == 0) {

//                         u_a_min = eta(fejer_nodes_u_int_[nu_loc+1], flag_u_loc);
//                         u_b_min = 1.0;

//                     } else if (nu_loc == Nu_int_-1) {

//                         u_a_min = -1.0;
//                         u_b_min = eta(fejer_nodes_u_int_[nu_loc-1], flag_u_loc);

//                     } else {

//                         u_a_min = eta(fejer_nodes_u_int_[nu_loc+1], flag_u_loc);
//                         u_b_min = eta(fejer_nodes_u_int_[nu_loc-1], flag_u_loc);

//                     }

//                     if (nv_loc == 0) {

//                         v_a_min = eta(fejer_nodes_v_int_[nv_loc+1], flag_v_loc);
//                         v_b_min = 1.0;

//                     } else if (nv_loc == Nv_int_-1) {

//                         v_a_min = -1.0;
//                         v_b_min = eta(fejer_nodes_v_int_[nv_loc-1], flag_v_loc);

//                     } else {

//                         v_a_min = eta(fejer_nodes_v_int_[nv_loc+1], flag_v_loc);
//                         v_b_min = eta(fejer_nodes_v_int_[nv_loc-1], flag_v_loc);

//                     } 

//                     if (GEOMETRY == 0) {
                    
//                         auto dist_func = [r_0, r_1, r_2, q, u_a_loc, u_b_loc, v_a_loc, v_b_loc, this](double s, double t) -> double {const double ss = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, s);
//                                                                                                                                 const double tt = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, t);
//                                                                                                                                 double rp_0, rp_1, rp_2;
//                                                                                                                                 parametrization_q(this->SPHERE_RADIUS, this->SPHERE_CENTER,ss, tt, q, rp_0, rp_1, rp_2);
//                                                                                                                                 return ((r_0 - rp_0)*(r_0 - rp_0) + (r_1 - rp_1)*(r_1 - rp_1) + (r_2 - rp_2)*(r_2 - rp_2));};
                
//                         golden_search(dist_func, G_SEARCH_MAX_ITER, G_SEARCH_TOL, u_a_min, u_b_min, v_a_min, v_b_min, ubar_loc, vbar_loc); 

//                     } else {

//                         InterpPatch patch = interp_surface_[q];

//                         auto dist_func = [r_0, r_1, r_2, q, u_a_loc, u_b_loc, v_a_loc, v_b_loc, patch](double s, double t) -> double {const std::vector<double> ss = {ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, s)};
//                                                                                                                                         const std::vector<double> tt = {ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, t)};
//                                                                                                                                         double rp_0, rp_1, rp_2;
//                                                                                                                                         lagrange_interpolation_2D(patch.uNodes, patch.vNodes, patch.uWeights, patch.vWeights,
//                                                                                                                                                                 patch.x, ss, tt, &rp_0);
//                                                                                                                                         lagrange_interpolation_2D(patch.uNodes, patch.vNodes, patch.uWeights, patch.vWeights,
//                                                                                                                                                                 patch.y, ss, tt, &rp_1);
//                                                                                                                                         lagrange_interpolation_2D(patch.uNodes, patch.vNodes, patch.uWeights, patch.vWeights,
//                                                                                                                                                                 patch.z, ss, tt, &rp_2);
//                                                                                                                                         return ((r_0 - rp_0)*(r_0 - rp_0) + (r_1 - rp_1)*(r_1 - rp_1) + (r_2 - rp_2)*(r_2 - rp_2));};
                
//                         golden_search(dist_func, G_SEARCH_MAX_ITER, G_SEARCH_TOL, u_a_min, u_b_min, v_a_min, v_b_min, ubar_loc, vbar_loc); 
                        
//                     }

//                 }    

//                 std::vector<std::complex<double>> prec_loc(Nu_int_*Nv_int_);
                        
//                 beta(r_0, r_1, r_2,
//                     q, flag_u_loc, flag_v_loc,
//                     u_a_loc, u_b_loc, v_a_loc, v_b_loc,
//                     ubar_loc, vbar_loc, 
//                     prec_loc);

//                 precomputations_thread.insert(precomputations_thread.end(), prec_loc.begin(), prec_loc.end());
//                 point_precomputations_thread.push_back(npoint);
//                 patch_num_precomputations_thread.push_back(patch_num);

//             }

//         }

//     }

//     #pragma omp for ordered
//     for (int i = 0; i < NTHREADS; i++) {

//         #pragma omp ordered
//         {

//             precomputations_.insert(precomputations_.end(), precomputations_thread.begin(), precomputations_thread.end());
//             point_precomputations_.insert(point_precomputations_.end(), point_precomputations_thread.begin(), point_precomputations_thread.end());
//             patch_num_precomputations_.insert(patch_num_precomputations_.end(), patch_num_precomputations_thread.begin(), patch_num_precomputations_thread.end());
            
//             std::vector<std::complex<double>>().swap(precomputations_thread);
//             std::vector<long long>().swap(point_precomputations_thread);
//             std::vector<long long>().swap(patch_num_precomputations_thread);

//         }

//     }

//     }

//     std::vector<long long> patch_num_coeffs_aux = patch_num_precomputations_;
//     std::sort(patch_num_coeffs_aux.begin(), patch_num_coeffs_aux.end());
//     auto last_unique = std::unique(patch_num_coeffs_aux.begin(), patch_num_coeffs_aux.end());
//     patch_num_coeffs_aux.erase(last_unique, patch_num_coeffs_aux.end());
//     patch_num_coeffs_ = std::move(patch_num_coeffs_aux);

//     std::vector<long long>().swap(sing_and_near_sing_patches_estimate_);
//     std::vector<long long>().swap(start_sing_and_near_sing_patches_estimate_);
//     std::vector<long long>().swap(size_sing_and_near_sing_patches_estimate_);

// }

// void Solver::get_coeffs(const long long q, const std::complex<double>* phi,
//                 std::complex<double>* coeffs)
// {

//     std::vector<std::complex<double>> evals(Nu_int_*Nv_int_);

//     for (int i = 0; i < Nu_int_; i++) {
//         for (int j = 0; j < Nv_int_; j++) {

//             evals[i*Nv_int_ + j] = dsdtjac_all_[q * Nu_int_*Nv_int_ + i*Nv_int_ + j] * phi[i*Nv_int_ + j];

//         }
//     }

//     std::vector<std::complex<double>> matprod1(Nu_int_*Nv_int_);

//     std::complex<double> one(1.0, 0.0);
//     std::complex<double> zero(0.0, 0.0);

//     cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//                 Nu_int_, Nv_int_, Nu_int_, 
//                 &one, &Tn_[0], Nu_int_, &evals[0], Nv_int_, &zero, &matprod1[0], Nv_int_);
    
//     cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
//                 Nu_int_, Nv_int_, Nv_int_,
//                 &one, &matprod1[0], Nv_int_, &Tm_[0], Nv_int_, &zero, &coeffs[0], Nv_int_);
    
//     for (int n = 0; n < Nu_int_; n++) {

//         double alpha_n;

//         if (n == 0) {
//             alpha_n = 1;
//         } else {
//             alpha_n = 2;
//         }

//         for (int m = 0; m < Nv_int_; m++) {

//             double alpha_m;

//             if (m == 0) {
//                 alpha_m = 1;
//             } else {
//                 alpha_m = 2;
//             }

//             coeffs[n*Nv_int_ + m] *= alpha_n * alpha_m / (Nu_int_*Nv_int_);

//         }

//     }

// }

// void Solver::compute_coeffs(const VectorXcd& phi,
//                     std::vector<std::complex<double>>& vec_coeffs)
// {
    
//     #pragma omp parallel for
//     for (long long i = 0; i < patch_num_coeffs_.size(); i++) {

//         const long long patch_num = patch_num_coeffs_[i];

//         get_coeffs(patch_num, &phi[patch_num * Nu_int_*Nv_int_], &vec_coeffs[i * Nu_int_*Nv_int_]);

//     }

// }

// void Solver::int_near(const std::complex<double>* coeffs,
//                 const std::complex<double>* precomputations,
//                 std::complex<double>& solution)
// {         

//     solution = std::complex<double>(0.0, 0.0);

//     for (int i = 0; i < Nu_int_; i++) {
//         for (int j = 0; j < Nv_int_; j++) {

//             solution += precomputations[i*Nv_int_ + j] * coeffs[i*Nv_int_ + j];

//         }
//     }

// }

// void Solver::int_far(const double r_0, const double r_1, const double r_2,
//                 const long long npatch, 
//                 const std::complex<double>* phi,
//                 std::complex<double>& solution) 
                
// {

//     solution = std::complex<double>(0.0, 0.0);

//     for (int i = 0; i < Nu_int_; i++) {
//         for (int j = 0; j < Nv_int_; j++) {

//             const double dsdtjac_loc = dsdtjac_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double constant = dsdtjac_loc * fejer_weights_u_int_[i] * fejer_weights_v_int_[j];

//             const double px = disc_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double py = disc_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double pz = disc_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

//             const double nx = norm_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double ny = norm_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double nz = norm_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

//             std::complex<double> kernel;
//             HH2(r_0, r_1, r_2, px, py, pz, nx, ny, nz, coupling_parameter_, wavenumber_, kernel);
            
//             solution += constant * kernel * phi[i*Nv_int_ + j];

//         }
//     }

// }

// void Solver::compute_integral(const VectorXcd& phi, const std::vector<std::complex<double>>& vec_coeffs,
//                         VectorXcd& integral)
// {

//     #pragma omp parallel for
//     for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

//         const double r_0 = disc_points_x_all_[point_low_ + npoint];
//         const double r_1 = disc_points_y_all_[point_low_ + npoint];
//         const double r_2 = disc_points_z_all_[point_low_ + npoint];

//         const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);
//         const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);

//         integral[npoint] = std::complex<double>(0.0, 0.0);

//         for (long long q = 0; q < Q_; q++) {             

//             for (int q_x = 0; q_x < Qx_; q_x++) {
//                 for (int q_y = 0; q_y < Qy_; q_y++) {

//                     const long long patch_num = q * Qx_*Qy_ + q_x * Qy_ + q_y;

//                     const auto patch_begin = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_begin);
//                     const auto patch_end = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_end);                            
//                     const auto pos = std::lower_bound(patch_begin, patch_end, patch_num);

//                     const auto idx = std::distance(patch_num_precomputations_.begin(), pos);
                    
//                     std::complex<double> solution(0.0,0.0);

//                     if ((point_precomputations_[idx] == (point_low_ + npoint)) && (patch_num_precomputations_[idx] == patch_num)) {

//                         const long long start_precomputations = idx * Nu_int_*Nv_int_;
                        
//                         const auto pos_coeffs = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
//                         const auto idx_coeffs = std::distance(patch_num_coeffs_.begin(), pos_coeffs);

//                         const long long start_coeffs = idx_coeffs * Nu_int_*Nv_int_;

//                         int_near(&vec_coeffs[start_coeffs], 
//                                     &precomputations_[start_precomputations], 
//                                     solution);
                        
//                     } else {

//                         int_far(r_0, r_1, r_2,
//                                 patch_num, 
//                                 &phi[patch_num * Nu_int_*Nv_int_],
//                                 solution);

//                     }

//                     integral[npoint] += solution;

//                 }
//             }

//         }

//         if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

//             integral[npoint] += 0.5 * phi[point_low_ + npoint];

//         }                

//     }

// }

// void Solver::iterator_function_unacc(const VectorXcd& phi,
//                                 VectorXcd& rhs)
// {

//     std::vector<std::complex<double>> vec_coeffs(patch_num_coeffs_.size() * Nu_int_*Nv_int_);

//     compute_coeffs(phi, vec_coeffs);
                
//     compute_integral(phi, vec_coeffs, rhs);         

// }

// void Solver::setup_IFGF_choose_levels() {

//     int  nlevels_old = nlevels_;
//     bool loc_cond_satisfied;
//     bool glob_cond_satisfied;

//     while (true) {
//         create_IFGF_object();
//         compute_new_order_points_RP();
//         loc_cond_satisfied = check_patch_in_neighbours();

//         MPI_Allreduce(&loc_cond_satisfied, &glob_cond_satisfied, 1, MPI_C_BOOL, MPI_LAND, mpi_comm_);

//         // Case where we go from to many levels to just the right amount
//         if (nlevels_ <= 3) {
//             nlevels_ = 3;
//             break;
//         } else if(glob_cond_satisfied && nlevels_old > nlevels_) {

//              break;

//         // Case where on the first try it worked, so we try to go deeper
//         } else if (glob_cond_satisfied && nlevels_old <= nlevels_) {
//             nlevels_old = nlevels_;
//             nlevels_ += 1;

//         // Case where the first try it didn't work, so we move down a level
//         } else if (!glob_cond_satisfied) {
//             nlevels_old = nlevels_;
//             nlevels_ -= 1;
//         } 

//     }


// }

// void Solver::create_IFGF_object() 
// {

//     std::vector<double> x, y, z;
//     std::vector<double> normal_x, normal_y, normal_z;

//     x.insert(x.end(), disc_points_x_all_.begin() + point_low_, disc_points_x_all_.begin() + point_up_);
//     y.insert(y.end(), disc_points_y_all_.begin() + point_low_, disc_points_y_all_.begin() + point_up_);
//     z.insert(z.end(), disc_points_z_all_.begin() + point_low_, disc_points_z_all_.begin() + point_up_);

//     normal_x.insert(normal_x.end(), norm_points_x_all_.begin() + point_low_, norm_points_x_all_.begin() + point_up_);
//     normal_y.insert(normal_y.end(), norm_points_y_all_.begin() + point_low_, norm_points_y_all_.begin() + point_up_);
//     normal_z.insert(normal_z.end(), norm_points_z_all_.begin() + point_low_, norm_points_z_all_.begin() + point_up_);

//     boxes_.InitializeObject(x, y, z, normal_x, normal_y, normal_z, 
//                             split_points_2_, recv_counts_2_, displs_2_,
//                             coupling_parameter_,
//                             new_order_points_IFGF_,
//                             nlevels_, wavenumber_,
//                             mpi_comm_);

//     std::vector<long long> point_precomputations_all;
//     std::vector<long long> patch_num_precomputations_all;

//     int size_loc = point_precomputations_.size();

//     std::vector<int> recv_counts(comm_size_);
//     std::vector<int> displs(comm_size_, 0);

//     MPI_Allgather(&size_loc, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, mpi_comm_);

//     int size_all = recv_counts[0];
//     for (int i = 1; i < comm_size_; i++) {
//         size_all += recv_counts[i];
//         displs[i] = displs[i-1] + recv_counts[i-1];
//     }

//     point_precomputations_all.resize(size_all);
//     patch_num_precomputations_all.resize(size_all);

//     MPI_Allgatherv(&point_precomputations_[0], size_loc, MPI_LONG_LONG, &point_precomputations_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
//     MPI_Allgatherv(&patch_num_precomputations_[0], size_loc, MPI_LONG_LONG, &patch_num_precomputations_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

//     std::unordered_map<long long, std::unordered_set<long long>> precomputations_data_all;

//     for (long long i = 0; i < size_all; i++) {

//         const long long point = point_precomputations_all[i];
//         const long long patch = patch_num_precomputations_all[i];

//         precomputations_data_all[point].insert(patch);

//     }

//     std::vector<long long>().swap(point_precomputations_all);
//     std::vector<long long>().swap(patch_num_precomputations_all);

//     std::unordered_map<long long, std::unordered_set<long long>> precomputations_data_loc;

//     for (long long i = point_low_; i < point_up_; i++) {

//         const long long point = new_order_points_IFGF_[i];

//         precomputations_data_loc[point] = precomputations_data_all[point];

//     }

//     precomputations_data_all.clear();

//     boxes_.set_precomputations_data(precomputations_data_loc);

//     precomputations_data_loc.clear();

//     boxes_.set_n_pts_per_patch(Nu_int_, Nv_int_);

// }

// void Solver::compute_new_order_points_RP() {

//     std::vector<std::vector<long long>> points_not_in_rank(comm_size_), order_points_not_in_rank(comm_size_);

//     new_order_points_RP_ = std::vector<long long>(point_up_-point_low_);

//     for (long long i = point_low_; i < point_up_; i++) {

//         const long long point = new_order_points_IFGF_[i];

//         int rank;

//         for (int j = 0; j < comm_size_+1; j++) {

//             if ((point >= split_points_2_[j]) && (point < split_points_2_[j+1])) {

//                 rank = j;
//                 break;

//             }

//         }

//         if (comm_rank_ == rank) {

//             new_order_points_RP_[point-point_low_] = i;

//         } else {

//             points_not_in_rank[rank].push_back(point);
//             order_points_not_in_rank[rank].push_back(i);

//         }

//     }

//     for (int rank = 0; rank < comm_size_; rank++) {

//         std::vector<long long> points = points_not_in_rank[rank];
//         std::vector<long long> orders = order_points_not_in_rank[rank];

//         std::vector<long long> points_all;
//         std::vector<long long> orders_all;

//         int size_points = points.size();
//         int size_total;
        
//         std::vector<int> recv_counts(comm_size_);
//         std::vector<int> displs(comm_size_, 0);

//         MPI_Gather(&size_points, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, rank, mpi_comm_);

//         if (rank == comm_rank_) {

//             size_total = recv_counts[0];

//             for (int i = 1; i < comm_size_; i++) {

//                 size_total += recv_counts[i];

//                 displs[i] = displs[i-1] + recv_counts[i-1];

//             }

//             points_all.resize(size_total);
//             orders_all.resize(size_total);

//         }

//         MPI_Gatherv(&points[0], size_points, MPI_LONG_LONG, &points_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, rank, mpi_comm_);
//         MPI_Gatherv(&orders[0], size_points, MPI_LONG_LONG, &orders_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, rank, mpi_comm_);

//         std::vector<long long>().swap(points);
//         std::vector<long long>().swap(orders);

//         if (rank == comm_rank_) {

//             for (long long i = 0; i < size_total; i++) {

//                 new_order_points_RP_[points_all[i]-point_low_] = orders_all[i];

//             }

//         }

//         std::vector<long long>().swap(orders_all);
//         std::vector<long long>().swap(points_all);

//     }

// }

// bool Solver::check_patch_in_neighbours() {

//     // Compute point to box

//     std::vector<long long> point_to_box(point_up_-point_low_);

//     for (long long i = 0; i < point_up_-point_low_; i++) {

//         point_to_box[i] = boxes_.get_box(new_order_points_RP_[i]);

//     }

//     // Compute patch to box

//     std::unordered_map<long long, std::set<long long>> patch_to_box;

//     std::vector<long long> patches_loc;
//     std::vector<long long> boxes_loc;
//     std::vector<long long> total_boxes_loc;

//     for (long long patch_num = 0; patch_num < patch_up_-patch_low_; patch_num++) {

//         std::unordered_set<long long> boxes(point_to_box.begin() + patch_num * Nu_int_*Nv_int_, point_to_box.begin() + (patch_num+1) * Nu_int_*Nv_int_);

//         patches_loc.push_back(patch_num + patch_low_);
//         total_boxes_loc.push_back(boxes.size());
//         boxes_loc.insert(boxes_loc.end(), boxes.begin(), boxes.end());

//     }

//     int size_1_loc = patches_loc.size();
//     int size_2_loc = boxes_loc.size();

//     std::vector<int> recv_counts_1(comm_size_);
//     std::vector<int> recv_counts_2(comm_size_);
//     std::vector<int> displs_1(comm_size_, 0);
//     std::vector<int> displs_2(comm_size_, 0);

//     MPI_Allgather(&size_1_loc, 1, MPI_INT, &recv_counts_1[0], 1, MPI_INT, mpi_comm_);
//     MPI_Allgather(&size_2_loc, 1, MPI_INT, &recv_counts_2[0], 1, MPI_INT, mpi_comm_);

//     int size_1_all = recv_counts_1[0];
//     int size_2_all = recv_counts_2[0];

//     for (int i = 1; i < comm_size_; i++) {

//         size_1_all += recv_counts_1[i];
//         size_2_all += recv_counts_2[i];

//         displs_1[i] = displs_1[i-1] + recv_counts_1[i-1];
//         displs_2[i] = displs_2[i-1] + recv_counts_2[i-1];

//     }

//     std::vector<long long> patches_all(size_1_all);
//     std::vector<long long> total_boxes_all(size_1_all);
//     std::vector<long long> boxes_all(size_2_all);

//     MPI_Allgatherv(&patches_loc[0], size_1_loc, MPI_LONG_LONG, &patches_all[0], &recv_counts_1[0], &displs_1[0], MPI_LONG_LONG, mpi_comm_);
//     MPI_Allgatherv(&total_boxes_loc[0], size_1_loc, MPI_LONG_LONG, &total_boxes_all[0], &recv_counts_1[0], &displs_1[0], MPI_LONG_LONG, mpi_comm_);
//     MPI_Allgatherv(&boxes_loc[0], size_2_loc, MPI_LONG_LONG, &boxes_all[0], &recv_counts_2[0], &displs_2[0], MPI_LONG_LONG, mpi_comm_);

//     std::vector<long long>().swap(patches_loc);
//     std::vector<long long>().swap(total_boxes_loc);
//     std::vector<long long>().swap(boxes_loc);

//     long long counter = 0;

//     for (long long i = 0; i < size_1_all; i++) {

//         long long patch = patches_all[i];
//         long long size_boxes = total_boxes_all[i];

//         patch_to_box[patch].insert(boxes_all.begin() + counter, boxes_all.begin() + counter + size_boxes);

//         counter += size_boxes;

//     }

//     std::vector<long long>().swap(patches_all);
//     std::vector<long long>().swap(total_boxes_all);
//     std::vector<long long>().swap(boxes_all);
    
//     // Compute relevant morton boxes

//     std::unordered_set<long long> mortonidofrelboxes;
//     boxes_.get_mortonidofrelboxes(mortonidofrelboxes);

//     std::vector<long long> mortonidofrelboxes_loc(mortonidofrelboxes.begin(), mortonidofrelboxes.end());
//     mortonidofrelboxes.clear();

//     int size_loc = mortonidofrelboxes_loc.size();

//     std::vector<int> recv_counts(comm_size_);
//     std::vector<int> displs(comm_size_, 0);

//     MPI_Allgather(&size_loc, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, mpi_comm_);

//     int size_all = recv_counts[0];

//     for (int i = 1; i < comm_size_; i++) {

//         displs[i] = displs[i-1] + recv_counts[i-1];

//         size_all += recv_counts[i];

//     }

//     std::vector<long long> mortonidofrelboxes_all(size_all);

//     MPI_Allgatherv(&mortonidofrelboxes_loc[0], size_loc, MPI_LONG_LONG, &mortonidofrelboxes_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

//     std::vector<long long>().swap(mortonidofrelboxes_loc);

//     for (long long i = 0; i < size_all; i++) {

//         mortonidofrelboxes.insert(mortonidofrelboxes_all[i]);

//     }

//     std::vector<long long>().swap(mortonidofrelboxes_all);

//     // Check if patch is in neighbour union

//     long long actual_npoint = -1;
//     std::vector<long long> neighbours_box_npoint;

//     for (long long i = 0; i < point_precomputations_.size(); i++) {

//         const long long npoint = point_precomputations_[i];

//         if (actual_npoint != npoint) {
        
//             const long long box_npoint = point_to_box[npoint - point_low_];
//             neighbours_box_npoint = boxes_.get_neighbours_box(box_npoint);

//             std::vector<long long> relevant_neighbours;
//             for (long long j = 0; j < neighbours_box_npoint.size(); j++) {

//                 if (mortonidofrelboxes.count(neighbours_box_npoint[j]) != 0) {
//                     relevant_neighbours.push_back(neighbours_box_npoint[j]);
//                 }

//             }

//             neighbours_box_npoint = relevant_neighbours;

//             actual_npoint = npoint;

//         }

//         const long long patch_num_singular = patch_num_precomputations_[i];
//         const std::set<long long> boxes_patch_num_singular = patch_to_box[patch_num_singular];

//         if (!(std::includes(neighbours_box_npoint.begin(), neighbours_box_npoint.end(), boxes_patch_num_singular.begin(), boxes_patch_num_singular.end()))) {

//             //std::cout << "Point " << npoint << " and singular / near singular patches are not contained in neighbor union.\n";
//             // std::exit(0);
//             return false;

//         }

//     }

//     return true;

// }



// void Solver::initialize_indexes_HO() {

//     boxes_.InitializeIndexes(
//     [&](double a1,double a2,double a3,double a4,double a5,double a6,
//         double a7,double a8,double a9,double a10,
//         std::complex<double> c1,std::complex<double> c2,std::complex<double>& out) {
//         Solver::fct_4(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,c1,c2,out);
//     }
//     );

//     // boxes_.InitializeIndexes<fct_4>();

// }

// void Solver::compute_intensities_patch(const long long npatch,
//                                 const std::complex<double>* phi,
//                                 std::complex<double>* intensities)
// {        

//     for (int i = 0; i < Nu_int_; i++) {
//         for (int j = 0; j < Nv_int_; j++) {

//             intensities[i * Nv_int_ + j] = dsdtjac_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j] * fejer_weights_u_int_[i] * fejer_weights_v_int_[j] * phi[i * Nv_int_ + j];

//         }
//     }   

// }

// void Solver::compute_intensities(const VectorXcd& phi,
//                             std::vector<std::complex<double>>& intensities)
// {

//     intensities = std::vector<std::complex<double>>(point_up_-point_low_);

//     #pragma omp parallel for
//     for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

//         compute_intensities_patch(patch_num, 
//                                     &phi[patch_num * Nu_int_*Nv_int_],
//                                     &intensities[(patch_num - patch_low_) * Nu_int_*Nv_int_]);

//     }                    

//     std::vector<std::complex<double>> intensities_all(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
    
//     MPI_Allgatherv(&intensities[0], point_up_-point_low_, MPI_DOUBLE_COMPLEX, &intensities_all[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE_COMPLEX, mpi_comm_);
    
//     intensities = std::vector<std::complex<double>>(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);
    
//     #pragma omp parallel for
//     for (long long i = 0; i < Q_ * Qx_*Qy_ * Nu_int_*Nv_int_; i++) {

//         intensities[i] = intensities_all[new_order_points_IFGF_[i]];

//     }

// }

// void Solver::compute_integral_acc(const VectorXcd& phi, const std::vector<std::complex<double>>& vec_coeffs,
//                             VectorXcd& integral)
// {

//     integral.setZero();

//     #pragma omp parallel for
//     for (long long i = 0; i < point_precomputations_.size(); i++) {

//         const long long npoint = point_precomputations_[i];
//         const long long patch_num = patch_num_precomputations_[i];
        
//         const long long start_precomputations = i * Nu_int_*Nv_int_;

//         const auto pos_coeffs = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
//         const auto idx_coeffs = std::distance(patch_num_coeffs_.begin(), pos_coeffs);

//         const long long start_coeffs = idx_coeffs * Nu_int_*Nv_int_;

//         std::complex<double> solution;

//         int_near(&vec_coeffs[start_coeffs], 
//                     &precomputations_[start_precomputations], 
//                     solution);

//         integral[npoint - point_low_] += solution;

//     }

//     if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

//         #pragma omp parallel for
//         for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

//             integral[npoint] += 0.5 * phi[point_low_ + npoint];

//         }

//     }

// } 

// inline void Solver::fct_4(const double x1, const double x2, const double x3,
//                             const double y1, const double y2, const double y3,
//                             const double normal1, const double normal2, const double normal3,
//                             const double coupling_parameter, const std::complex<double> wavenumber,
//                             const std::complex<double> density, 
//                             std::complex<double>& phi) 
// {

//     std::complex<double> kernel;

//     HH2(x1, x2, x3, 
//         y1, y2, y3,
//         normal1, normal2, normal3,
//         coupling_parameter, wavenumber,
//         kernel);

//     const double kerreal = kernel.real();
//     const double kerimag = kernel.imag();

//     const double phireal = density.real() * kerreal - density.imag() * kerimag;
//     const double phiimag = density.real() * kerimag + density.imag() * kerreal;

//     phi = {phireal, phiimag};

// }

// inline void Solver::fac_1(const double distance, const std::complex<double> wavenumber, std::complex<double>& sol)
// {

//     const double re = std::cos(wavenumber.real() * distance) / distance;
//     const double im = std::sin(wavenumber.real() * distance) / distance;

//     sol = {re, im};

// }

// void Solver::iterator_function_acc(const VectorXcd& phi,
//                             VectorXcd& rhs)
// {

//     std::vector<std::complex<double>> intensities; 

//     compute_intensities(phi, intensities);

//     std::vector<std::complex<double>> vec_coeffs(patch_num_coeffs_.size() * Nu_int_*Nv_int_);

//     compute_coeffs(phi, vec_coeffs);

//     compute_integral_acc(phi, vec_coeffs, rhs);  

//     std::vector<std::complex<double>> solution_1;
    
//     solution_1 = intensities;

//     // boxes_.Solve<&fct_4, &fac_1>(solution_1); 
//     boxes_.Solve(
//     solution_1,
//     [&](double a1,double a2,double a3,double a4,double a5,double a6,
//         double a7,double a8,double a9,double a10,
//         std::complex<double> c1,std::complex<double> c2,std::complex<double>& out) {
//         Solver::fct_4(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,c1,c2,out);
//     },
//     [&](double a1,std::complex<double> c1,std::complex<double>& out) {
//         Solver::fac_1(a1,c1,out);
//     }
// );

//     #pragma omp parallel for
//     for (long long i = 0; i < point_up_-point_low_; i++) {

//         rhs[i] += solution_1[new_order_points_RP_[i]];
        
//     }   

// }

// void Solver::create_fftw_objects()
// {

//     if (METHOD_FFTW == 1) {

//         fftw_init_threads();

//         fftw_plan_with_nthreads(NTHREADS);

//         int rank = 2;

//         int n_orig[2] = {Nu_int_, Nv_int_};
//         int n_overs[2] = {Nu_int_*P_over_, Nv_int_*P_over_};

//         original_patch_ = (double*) fftw_malloc(sizeof(double) * n_orig[0]*n_orig[1] * patch_num_coeffs_.size() * 2);
//         oversampled_patch_ = (double*) fftw_malloc(sizeof(double) * n_overs[0]*n_overs[1] * patch_num_coeffs_.size() * 2);   
            
//         int istride_orig = 1;
//         int idist_orig = n_orig[0]*n_orig[1];
        
//         int istride_overs = 1;
//         int idist_overs = n_overs[0]*n_overs[1];

//         fftw_r2r_kind kind_orig[2] = {FFTW_REDFT10, FFTW_REDFT10};
//         fftw_r2r_kind kind_overs[2] = {FFTW_REDFT01, FFTW_REDFT01};

//         if (TYPE_FFTW == 1) {

//             plan_original_patch_ = fftw_plan_many_r2r(rank, n_orig, patch_num_coeffs_.size() * 2, original_patch_, n_orig, istride_orig, idist_orig, original_patch_, n_orig, istride_orig, idist_orig, kind_orig, FFTW_PATIENT);
//             plan_oversampled_patch_ = fftw_plan_many_r2r(rank, n_overs, patch_num_coeffs_.size() * 2, oversampled_patch_, n_overs, istride_overs, idist_overs, oversampled_patch_, n_overs, istride_overs, idist_overs, kind_overs, FFTW_PATIENT);

//         } else {

//             fftw_iodim dims_orig[rank];
//             dims_orig[1].n = n_orig[1];
//             dims_orig[1].is = istride_orig;
//             dims_orig[1].os = istride_orig;
//             dims_orig[0].n = n_orig[0];
//             dims_orig[0].is = n_orig[1] * dims_orig[1].is;
//             dims_orig[0].os = n_orig[1] * dims_orig[1].os;

//             fftw_iodim dims_overs[rank];
//             dims_overs[1].n = n_overs[1];
//             dims_overs[1].is = istride_overs;
//             dims_overs[1].os = istride_overs;
//             dims_overs[0].n = n_overs[0];
//             dims_overs[0].is = n_overs[1] * dims_overs[1].is;
//             dims_overs[0].os = n_overs[1] * dims_overs[1].os;

//             int howmany_rank = 1;

//             fftw_iodim howmany_dims_orig[howmany_rank]; 
//             howmany_dims_orig[0].n = patch_num_coeffs_.size() * 2;
//             howmany_dims_orig[0].is = idist_orig;
//             howmany_dims_orig[0].os = idist_orig;

//             fftw_iodim howmany_dims_overs[howmany_rank]; 
//             howmany_dims_overs[0].n = patch_num_coeffs_.size() * 2;
//             howmany_dims_overs[0].is = idist_overs;
//             howmany_dims_overs[0].os = idist_overs;

//             plan_original_patch_ = fftw_plan_guru_r2r(rank, dims_orig, howmany_rank, howmany_dims_orig, original_patch_, original_patch_, kind_orig, FFTW_PATIENT);
//             plan_oversampled_patch_ = fftw_plan_guru_r2r(rank, dims_overs, howmany_rank, howmany_dims_overs, oversampled_patch_, oversampled_patch_, kind_overs, FFTW_PATIENT);

//         }

//     } else {

//         int rank = 2;

//         int n_orig[2] = {Nu_int_, Nv_int_};
//         int n_overs[2] = {Nu_int_*P_over_, Nv_int_*P_over_};
        
//         original_patch_ = (double*) fftw_malloc(sizeof(double) * n_orig[0]*n_orig[1] * 2);
//         oversampled_patch_ = (double*) fftw_malloc(sizeof(double) * n_overs[0]*n_overs[1] * 2);   

//         int istride_orig = 1;
//         int idist_orig = n_orig[0]*n_orig[1];
        
//         int istride_overs = 1;
//         int idist_overs = n_overs[0]*n_overs[1];

//         fftw_r2r_kind kind_orig[2] = {FFTW_REDFT10, FFTW_REDFT10};
//         fftw_r2r_kind kind_overs[2] = {FFTW_REDFT01, FFTW_REDFT01};

//         if (TYPE_FFTW == 1) {

//             plan_original_patch_ = fftw_plan_many_r2r(rank, n_orig, 2, original_patch_, n_orig, istride_orig, idist_orig, original_patch_, n_orig, istride_orig, idist_orig, kind_orig, FFTW_PATIENT);
//             plan_oversampled_patch_ = fftw_plan_many_r2r(rank, n_overs, 2, oversampled_patch_, n_overs, istride_overs, idist_overs, oversampled_patch_, n_overs, istride_overs, idist_overs, kind_overs, FFTW_PATIENT);

//         } else {

//             fftw_iodim dims_orig[rank];
//             dims_orig[1].n = n_orig[1];
//             dims_orig[1].is = istride_orig;
//             dims_orig[1].os = istride_orig;
//             dims_orig[0].n = n_orig[0];
//             dims_orig[0].is = n_orig[1] * dims_orig[1].is;
//             dims_orig[0].os = n_orig[1] * dims_orig[1].os;

//             fftw_iodim dims_overs[rank];
//             dims_overs[1].n = n_overs[1];
//             dims_overs[1].is = istride_overs;
//             dims_overs[1].os = istride_overs;
//             dims_overs[0].n = n_overs[0];
//             dims_overs[0].is = n_overs[1] * dims_overs[1].is;
//             dims_overs[0].os = n_overs[1] * dims_overs[1].os;

//             int howmany_rank = 1;

//             fftw_iodim howmany_dims_orig[howmany_rank]; 
//             howmany_dims_orig[0].n = 2;
//             howmany_dims_orig[0].is = idist_orig;
//             howmany_dims_orig[0].os = idist_orig;

//             fftw_iodim howmany_dims_overs[howmany_rank]; 
//             howmany_dims_overs[0].n = 2;
//             howmany_dims_overs[0].is = idist_overs;
//             howmany_dims_overs[0].os = idist_overs;

//             plan_original_patch_ = fftw_plan_guru_r2r(rank, dims_orig, howmany_rank, howmany_dims_orig, original_patch_, original_patch_, kind_orig, FFTW_PATIENT);
//             plan_oversampled_patch_ = fftw_plan_guru_r2r(rank, dims_overs, howmany_rank, howmany_dims_overs, oversampled_patch_, oversampled_patch_, kind_overs, FFTW_PATIENT);

//         }

//     }

// }

// void Solver::compute_oversampling_M1(const std::vector<double>& psi, std::vector<double>& psi_overs) 
// {

//     #pragma omp parallel for
//     for (long long i = 0; i < psi.size(); i++) {

//         original_patch_[i] = psi[i] / (Nu_int_ * Nv_int_);

//     }

//     fftw_execute(plan_original_patch_);

//     std::memset(oversampled_patch_, 0.0, patch_num_coeffs_.size() * Nu_int_*Nv_int_ * P_over_*P_over_ * 2 * sizeof(double));

//     #pragma omp parallel for
//     for (long long k = 0; k < patch_num_coeffs_.size(); k++) {

//         for (int i = 0; i < Nu_int_; i++) {
//             for (int j = 0; j < Nv_int_; j++) {

//                 oversampled_patch_[k * 2 * Nu_int_*Nv_int_ * P_over_*P_over_ + i * Nv_int_ * P_over_ + j] = original_patch_[k * 2 * Nu_int_*Nv_int_ + i * Nv_int_ + j] * 0.5 * 0.5;
//                 oversampled_patch_[k * 2 * Nu_int_*Nv_int_ * P_over_*P_over_ + Nu_int_*Nv_int_ * P_over_*P_over_ + i * Nv_int_ * P_over_ + j] = original_patch_[k * 2 * Nu_int_*Nv_int_ + Nu_int_*Nv_int_ + i * Nv_int_ + j] * 0.5 * 0.5;

//             }
//         }

//     }

//     fftw_execute(plan_oversampled_patch_);

//     psi_overs = std::vector<double>(oversampled_patch_, oversampled_patch_ + 2 * patch_num_coeffs_.size() * Nu_int_*Nv_int_ * P_over_*P_over_);

// }

// void Solver::compute_oversampling_M2(const std::vector<double>& psi, std::vector<double>& psi_overs) 
// {

//     for (long long i = 0; i < psi.size(); i++) {

//         original_patch_[i] = psi[i] / (Nu_int_ * Nv_int_);

//     }

//     fftw_execute(plan_original_patch_);

//     std::memset(oversampled_patch_, 0.0, Nu_int_*Nv_int_ * P_over_*P_over_ * 2 * sizeof(double));

//     for (int i = 0; i < Nu_int_; i++) {
//         for (int j = 0; j < Nv_int_; j++) {

//             oversampled_patch_[i * Nv_int_ * P_over_ + j] = original_patch_[i * Nv_int_ + j] * 0.5 * 0.5;
//             oversampled_patch_[Nu_int_*Nv_int_ * P_over_*P_over_ + i * Nv_int_ * P_over_ + j] = original_patch_[Nu_int_*Nv_int_ + i * Nv_int_ + j] * 0.5 * 0.5;

//         }
//     }

//     fftw_execute(plan_oversampled_patch_);

//     psi_overs = std::vector<double>(oversampled_patch_, oversampled_patch_ + 2 * Nu_int_*Nv_int_ * P_over_*P_over_);

// }

// void Solver::compute_singular_points()
// {

//     #pragma omp parallel
//     {

//     std::vector<long long> point_precomputations_thread;
//     std::vector<long long> patch_num_precomputations_thread;
//     std::vector<double> argmin_precomputations_u_thread;
//     std::vector<double> argmin_precomputations_v_thread;

//     #pragma omp for      
//     for (long long npoint = point_low_; npoint < point_up_; npoint++) {

//         const double r_0 = disc_points_x_all_[npoint];
//         const double r_1 = disc_points_y_all_[npoint];
//         const double r_2 = disc_points_z_all_[npoint];

//         const long long patch_num_own = npoint / (Nu_int_*Nv_int_);

//         const long long start = start_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];
//         const long long total = size_sing_and_near_sing_patches_estimate_[patch_num_own - patch_low_];

//         for (long long i = 0; i < total; i++) {

//             const long long patch_num = sing_and_near_sing_patches_estimate_[start + i];

//             double min_dist_loc;
//             int indx_point_min_dist_loc;

//             bool sing_near_sing_flag = false;

//             if (patch_num_own == patch_num) {

//                 min_dist_loc = 0.0;
//                 indx_point_min_dist_loc = npoint % (Nu_int_*Nv_int_);   

//                 sing_near_sing_flag = true; 

//             } else {  

//                 min_dist_loc = std::numeric_limits<double>::max();

//                 for (int ii = 0; ii < Nu_int_; ii++) {
//                     for (int jj = 0; jj < Nv_int_; jj++) {

//                         const long long position = patch_num * Nu_int_*Nv_int_ + ii * Nv_int_ + jj;

//                         const double diff_0 = disc_points_x_all_[position] - r_0;
//                         const double diff_1 = disc_points_y_all_[position] - r_1;
//                         const double diff_2 = disc_points_z_all_[position] - r_2;

//                         const double dist_aux = diff_0*diff_0 + diff_1*diff_1 + diff_2*diff_2;

//                         if (dist_aux < min_dist_loc) {

//                             min_dist_loc = dist_aux;
//                             indx_point_min_dist_loc = ii * Nv_int_ + jj;

//                         }

//                         const bool proximity = (std::abs(std::floor(r_0/proximity_) - std::floor(disc_points_x_all_[position]/proximity_)) <= 1) &&
//                                                 (std::abs(std::floor(r_1/proximity_) - std::floor(disc_points_y_all_[position]/proximity_)) <= 1) &&
//                                                 (std::abs(std::floor(r_2/proximity_) - std::floor(disc_points_z_all_[position]/proximity_)) <= 1);

//                         if (!sing_near_sing_flag && proximity) {

//                             sing_near_sing_flag = true;

//                         }

//                     }
//                 }

//                 min_dist_loc = std::sqrt(min_dist_loc);

//             }
        
//             if (sing_near_sing_flag) {   

//                 const int nu_loc = indx_point_min_dist_loc / Nv_int_;
//                 const int nv_loc = indx_point_min_dist_loc % Nv_int_;  

//                 const int flag_u_loc = flags_domain_u_all_[patch_num];
//                 const int flag_v_loc = flags_domain_v_all_[patch_num];

//                 const long long q = patch_num / (Qx_ * Qy_);
//                 const int q_x = (patch_num / Qy_) % Qx_;
//                 const int q_y = patch_num % Qy_;

//                 const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
//                 const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
//                 const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
//                 const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_; 

//                 double ubar_loc, vbar_loc;
                        
//                 if (min_dist_loc < 1E-12) {
                    
//                     ubar_loc = eta(fejer_nodes_u_int_[nu_loc], flag_u_loc);
//                     vbar_loc = eta(fejer_nodes_v_int_[nv_loc], flag_v_loc);

//                 } else {

//                     double u_a_min, u_b_min, v_a_min, v_b_min;

//                     if (nu_loc == 0) {

//                         u_a_min = eta(fejer_nodes_u_int_[nu_loc+1], flag_u_loc);
//                         u_b_min = 1.0;

//                     } else if (nu_loc == Nu_int_-1) {

//                         u_a_min = -1.0;
//                         u_b_min = eta(fejer_nodes_u_int_[nu_loc-1], flag_u_loc);

//                     } else {

//                         u_a_min = eta(fejer_nodes_u_int_[nu_loc+1], flag_u_loc);
//                         u_b_min = eta(fejer_nodes_u_int_[nu_loc-1], flag_u_loc);

//                     }

//                     if (nv_loc == 0) {

//                         v_a_min = eta(fejer_nodes_v_int_[nv_loc+1], flag_v_loc);
//                         v_b_min = 1.0;

//                     } else if (nv_loc == Nv_int_-1) {

//                         v_a_min = -1.0;
//                         v_b_min = eta(fejer_nodes_v_int_[nv_loc-1], flag_v_loc);

//                     } else {

//                         v_a_min = eta(fejer_nodes_v_int_[nv_loc+1], flag_v_loc);
//                         v_b_min = eta(fejer_nodes_v_int_[nv_loc-1], flag_v_loc);

//                     } 

//                     if (GEOMETRY == 0) {
                    
//                         auto dist_func = [r_0, r_1, r_2, q, u_a_loc, u_b_loc, v_a_loc, v_b_loc, this](double s, double t) -> double {const double ss = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, s);
//                                                                                                                                 const double tt = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, t);
//                                                                                                                                 double rp_0, rp_1, rp_2;
//                                                                                                                                 parametrization_q(this->SPHERE_RADIUS,this->SPHERE_CENTER,ss, tt, q, rp_0, rp_1, rp_2);
//                                                                                                                                 return ((r_0 - rp_0)*(r_0 - rp_0) + (r_1 - rp_1)*(r_1 - rp_1) + (r_2 - rp_2)*(r_2 - rp_2));};
                
//                         golden_search(dist_func, G_SEARCH_MAX_ITER, G_SEARCH_TOL, u_a_min, u_b_min, v_a_min, v_b_min, ubar_loc, vbar_loc); 

//                     } else {

//                         InterpPatch patch = interp_surface_[q];

//                         auto dist_func = [r_0, r_1, r_2, q, u_a_loc, u_b_loc, v_a_loc, v_b_loc, patch](double s, double t) -> double {const std::vector<double> ss = {ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, s)};
//                                                                                                                                         const std::vector<double> tt = {ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, t)};
//                                                                                                                                         double rp_0, rp_1, rp_2;
//                                                                                                                                         lagrange_interpolation_2D(patch.uNodes, patch.vNodes, patch.uWeights, patch.vWeights,
//                                                                                                                                                                 patch.x, ss, tt, &rp_0);
//                                                                                                                                         lagrange_interpolation_2D(patch.uNodes, patch.vNodes, patch.uWeights, patch.vWeights,
//                                                                                                                                                                 patch.y, ss, tt, &rp_1);
//                                                                                                                                         lagrange_interpolation_2D(patch.uNodes, patch.vNodes, patch.uWeights, patch.vWeights,
//                                                                                                                                                                 patch.z, ss, tt, &rp_2);
//                                                                                                                                         return ((r_0 - rp_0)*(r_0 - rp_0) + (r_1 - rp_1)*(r_1 - rp_1) + (r_2 - rp_2)*(r_2 - rp_2));};
                
//                         golden_search(dist_func, G_SEARCH_MAX_ITER, G_SEARCH_TOL, u_a_min, u_b_min, v_a_min, v_b_min, ubar_loc, vbar_loc); 

//                     }

//                 }    

//                 point_precomputations_thread.push_back(npoint);
//                 patch_num_precomputations_thread.push_back(patch_num);
//                 argmin_precomputations_u_thread.push_back(ubar_loc);
//                 argmin_precomputations_v_thread.push_back(vbar_loc);


//             }

//         }

//     }

//     #pragma omp for ordered
//     for (int i = 0; i < NTHREADS; i++) {

//         #pragma omp ordered
//         {

//             point_precomputations_.insert(point_precomputations_.end(), point_precomputations_thread.begin(), point_precomputations_thread.end());
//             patch_num_precomputations_.insert(patch_num_precomputations_.end(), patch_num_precomputations_thread.begin(), patch_num_precomputations_thread.end());
//             argmin_precomputations_u_.insert(argmin_precomputations_u_.end(), argmin_precomputations_u_thread.begin(), argmin_precomputations_u_thread.end());
//             argmin_precomputations_v_.insert(argmin_precomputations_v_.end(), argmin_precomputations_v_thread.begin(), argmin_precomputations_v_thread.end());

//             std::vector<long long>().swap(point_precomputations_thread);
//             std::vector<long long>().swap(patch_num_precomputations_thread);
//             std::vector<double>().swap(argmin_precomputations_u_thread);
//             std::vector<double>().swap(argmin_precomputations_v_thread);

//         }

//     }

//     }

//     std::vector<long long> patch_num_coeffs_aux = patch_num_precomputations_;
//     std::sort(patch_num_coeffs_aux.begin(), patch_num_coeffs_aux.end());
//     auto last_unique = std::unique(patch_num_coeffs_aux.begin(), patch_num_coeffs_aux.end());
//     patch_num_coeffs_aux.erase(last_unique, patch_num_coeffs_aux.end());
//     patch_num_coeffs_ = std::move(patch_num_coeffs_aux);

//     std::vector<long long>().swap(sing_and_near_sing_patches_estimate_);
//     std::vector<long long>().swap(start_sing_and_near_sing_patches_estimate_);
//     std::vector<long long>().swap(size_sing_and_near_sing_patches_estimate_);

// }

// void Solver::int_near_overs(const double r_0, const double r_1, const double r_2,
//                     const long long q, const int flag_u_loc, const int flag_v_loc, 
//                     const double u_a_loc, const double u_b_loc, const double v_a_loc, const double v_b_loc,
//                     const double ubar_loc, const double vbar_loc,
//                     const double* psi_overs,
//                     std::complex<double>& solution) 
// {

//     solution = std::complex<double>(0.0, 0.0);

//     std::vector<double> s(Nu_int_*P_over_ + 2);
//     std::vector<double> t(Nv_int_*P_over_ + 2);

//     s[0] = u_a_loc;

//     for (int i = 0; i < Nu_int_*P_over_; i++) {

//         const double xi = std::cos(M_PI * (2.0 * (Nu_int_*P_over_ - 1 - i) + 1.0) / (2.0 * Nu_int_*P_over_));
//         s[i+1] = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, eta(xi, flag_u_loc));

//     }          

//     s[Nu_int_*P_over_ + 1] = u_b_loc;

//     t[0] = v_a_loc;
    
//     for (int i = 0; i < Nv_int_*P_over_; i++) {

//         const double xi = std::cos(M_PI * (2.0 * (Nv_int_*P_over_ - 1 - i) + 1.0) / (2.0 * Nv_int_*P_over_));
//         t[i+1] = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, eta(xi, flag_v_loc));

//     }

//     t[Nv_int_*P_over_ + 1] = v_b_loc;                

//     std::vector<double> ui(Nu_prec_), si(Nu_prec_), muwi(Nu_prec_);
//     std::vector<double> vj(Nv_prec_), tj(Nv_prec_), muwj(Nv_prec_);

//     for (int i = 0; i < Nu_prec_; i++) {

//         ui[i] = xi(ubar_loc, fejer_nodes_u_prec_[i]);
//         si[i] = ab2cd(-1.0, 1.0, u_a_loc, u_b_loc, eta(ui[i], flag_u_loc));
//         muwi[i] = der_xi(ubar_loc, fejer_nodes_u_prec_[i]) * fejer_weights_u_prec_[i];
    
//     }

//     for (int j = 0; j < Nv_prec_; j++) {

//         vj[j] = xi(vbar_loc, fejer_nodes_v_prec_[j]);
//         tj[j] = ab2cd(-1.0, 1.0, v_a_loc, v_b_loc, eta(vj[j], flag_v_loc));
//         muwj[j] = der_xi(vbar_loc, fejer_nodes_v_prec_[j]) * fejer_weights_v_prec_[j];
    
//     }

//     std::array<double, N_LOC_X> s_loc;
//     std::array<double, N_LOC_Y> t_loc;
//     std::array<double, N_LOC_X*N_LOC_Y> psi_overs_loc_real, psi_overs_loc_imag;

//     std::vector<std::complex<double>> psi_overs_loc(Nu_prec_*Nv_prec_);
    
//     for (int i = 0; i < Nu_prec_; i++) {

//         const double ss = si[i];
//         auto idx_s = std::lower_bound(s.begin(), s.end(), ss);
//         int idx_ss = idx_s - s.begin();

//         if (idx_ss == 0) {
//             idx_ss = 1;
//         } else if (idx_ss == Nu_int_*P_over_ + 1) {
//             idx_ss = Nu_int_*P_over_;
//         }

//         if (idx_ss - (N_LOC_X - 1) / 2 < 1) {
//             idx_ss = (N_LOC_X - 1) / 2 + 1;
//         } else if (idx_ss + (N_LOC_X - 1) / 2 > Nu_int_*P_over_) {
//             idx_ss = Nu_int_*P_over_ - (N_LOC_X - 1) / 2;
//         }

//         for (int k = - (N_LOC_X - 1) / 2; k <= (N_LOC_X - 1) / 2; k++) {
//             s_loc[k + (N_LOC_X - 1) / 2] = s[idx_ss + k];
//         }

//         for (int j = 0; j < Nv_prec_; j++) {

//             const double tt = tj[j];
//             auto idx_t = std::lower_bound(t.begin(), t.end(), tt);
//             int idx_tt = idx_t - t.begin();

//             if (idx_tt == 0) {
//                 idx_tt = 1;
//             } else if (idx_tt == Nv_int_*P_over_ + 1) {
//                 idx_tt = Nv_int_*P_over_;
//             }

//             if (idx_tt - (N_LOC_Y - 1) / 2 < 1) {
//                 idx_tt = (N_LOC_Y - 1) / 2 + 1;
//             } else if (idx_tt + (N_LOC_Y - 1) / 2 > Nv_int_*P_over_) {
//                 idx_tt = Nv_int_*P_over_ - (N_LOC_Y - 1) / 2;
//             }

//             for (int l = - (N_LOC_Y - 1) / 2; l <= (N_LOC_Y - 1) / 2; l++) {
//                 t_loc[l + (N_LOC_Y - 1) / 2] = t[idx_tt + l];
//             }

//             for (int k = - (N_LOC_X - 1) / 2; k <= (N_LOC_X - 1) / 2; k++) {
//                 for (int l = - (N_LOC_Y - 1) / 2; l <= (N_LOC_Y - 1) / 2; l++) {

//                     const double real_val = psi_overs[(Nu_int_*P_over_ - idx_ss - k) * Nv_int_*P_over_ + (Nv_int_*P_over_ - idx_tt - l)];
//                     const double imag_val = psi_overs[Nu_int_*Nv_int_ * P_over_*P_over_ + (Nu_int_*P_over_ - idx_ss - k) * Nv_int_*P_over_ + (Nv_int_*P_over_ - idx_tt - l)];

//                     psi_overs_loc_real[(k + (N_LOC_X - 1) / 2) * N_LOC_Y + (l + (N_LOC_Y - 1) / 2)] = real_val;
//                     psi_overs_loc_imag[(k + (N_LOC_X - 1) / 2) * N_LOC_Y + (l + (N_LOC_Y - 1) / 2)] = imag_val;

//                 }
//             } 

//             const double interp_val_real = LocalInterpNewton2D<double, N_LOC_X, N_LOC_Y>(s_loc, t_loc, psi_overs_loc_real, ss, tt);
//             const double interp_val_imag = LocalInterpNewton2D<double, N_LOC_X, N_LOC_Y>(s_loc, t_loc, psi_overs_loc_imag, ss, tt);
//             const std::complex<double> interp_val(interp_val_real, interp_val_imag);

//             psi_overs_loc[i * Nv_prec_ + j] = interp_val;

//         }

//     }

//     if (GEOMETRY == 0) {

//         for (int i = 0; i < Nu_prec_; i++) {
//             for (int j = 0; j < Nv_prec_; j++) {       

//                 const double ss = si[i];
//                 const double tt = tj[j];                 

//                 double px, py, pz;
//                 double nx, ny, nz;
//                 double dxdsx, dxdsy, dxdsz;
//                 double dxdtx, dxdty, dxdtz;
//                 double dxdsdsx, dxdsdsy, dxdsdsz;
//                 double dxdsdtx, dxdsdty, dxdsdtz;
//                 double dxdtdtx, dxdtdty, dxdtdtz;

//                 parametrization_q(SPHERE_RADIUS, SPHERE_CENTER, ss, tt, q, px, py, pz);
//                 normal_q(SPHERE_RADIUS,ss, tt, q, nx, ny, nz);
//                 dxds_q(SPHERE_RADIUS,ss, tt, q, dxdsx, dxdsy, dxdsz);
//                 dxdt_q(SPHERE_RADIUS,ss, tt, q, dxdtx, dxdty, dxdtz);
//                 dxdsds_q(SPHERE_RADIUS,ss, tt, q, dxdsdsx, dxdsdsy, dxdsdsz);
//                 dxdsdt_q(SPHERE_RADIUS,ss, tt, q, dxdsdtx, dxdsdty, dxdsdtz);
//                 dxdtdt_q(SPHERE_RADIUS,ss, tt, q, dxdtdtx, dxdtdty, dxdtdtz);

//                 dxdsx *= 0.5 * (u_b_loc - u_a_loc);
//                 dxdsy *= 0.5 * (u_b_loc - u_a_loc);
//                 dxdsz *= 0.5 * (u_b_loc - u_a_loc);

//                 dxdtx *= 0.5 * (v_b_loc - v_a_loc);
//                 dxdty *= 0.5 * (v_b_loc - v_a_loc);
//                 dxdtz *= 0.5 * (v_b_loc - v_a_loc);

//                 dxdsdsx *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (u_b_loc - u_a_loc);
//                 dxdsdsy *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (u_b_loc - u_a_loc);
//                 dxdsdsz *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (u_b_loc - u_a_loc);

//                 dxdsdtx *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (v_b_loc - v_a_loc);
//                 dxdsdty *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (v_b_loc - v_a_loc);
//                 dxdsdtz *= 0.5 * (u_b_loc - u_a_loc) * 0.5 * (v_b_loc - v_a_loc);

//                 dxdtdtx *= 0.5 * (v_b_loc - v_a_loc) * 0.5 * (v_b_loc - v_a_loc);
//                 dxdtdty *= 0.5 * (v_b_loc - v_a_loc) * 0.5 * (v_b_loc - v_a_loc);
//                 dxdtdtz *= 0.5 * (v_b_loc - v_a_loc) * 0.5 * (v_b_loc - v_a_loc);

//                 std::complex<double> kernel;
//                 HH(r_0, r_1, r_2,
//                 px, py, pz,
//                 nx, ny, nz,
//                 dxdsx, dxdsy, dxdsz,
//                 dxdtx, dxdty, dxdtz,
//                 dxdsdsx, dxdsdsy, dxdsdsz,
//                 dxdsdtx, dxdsdty, dxdsdtz,
//                 dxdtdtx, dxdtdty, dxdtdtz,
//                 coupling_parameter_, wavenumber_,
//                 kernel);

//                 solution += kernel * muwi[i] * muwj[j] * psi_overs_loc[i * Nv_prec_ + j];

//             }
//         }

//     } else {

//         InterpPatch elem = interp_surface_[q];

//         std::vector<double> px(Nu_prec_*Nv_prec_), py(Nu_prec_*Nv_prec_), pz(Nu_prec_*Nv_prec_);
//         std::vector<double> nx(Nu_prec_*Nv_prec_), ny(Nu_prec_*Nv_prec_), nz(Nu_prec_*Nv_prec_);
//         std::vector<double> dxdsx(Nu_prec_*Nv_prec_), dxdsy(Nu_prec_*Nv_prec_), dxdsz(Nu_prec_*Nv_prec_);
//         std::vector<double> dxdtx(Nu_prec_*Nv_prec_), dxdty(Nu_prec_*Nv_prec_), dxdtz(Nu_prec_*Nv_prec_);

//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.x, si, tj, &px[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.y, si, tj, &py[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.z, si, tj, &pz[0]);

//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.nuX, si, tj, &nx[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.nuY, si, tj, &ny[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.nuZ, si, tj, &nz[0]);

//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dxdu, si, tj, &dxdsx[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dydu, si, tj, &dxdsy[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dzdu, si, tj, &dxdsz[0]);

//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dxdv, si, tj, &dxdtx[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dydv, si, tj, &dxdty[0]);
//         lagrange_interpolation_2D(elem.uNodes, elem.vNodes, elem.uWeights, elem.vWeights,
//                                     elem.dzdv, si, tj, &dxdtz[0]);

//         for (int i = 0; i < Nu_prec_; i++) {
//             for (int j = 0; j < Nv_prec_; j++) {

//                 std::complex<double> kernel;

//                 HH(r_0, r_1, r_2,
//                 px[i*Nv_prec_+j], py[i*Nv_prec_+j], pz[i*Nv_prec_+j],
//                 nx[i*Nv_prec_+j], ny[i*Nv_prec_+j], nz[i*Nv_prec_+j],
//                 dxdsx[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc), dxdsy[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc), dxdsz[i*Nv_prec_+j] * 0.5 * (u_b_loc - u_a_loc),
//                 dxdtx[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc), dxdty[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc), dxdtz[i*Nv_prec_+j] * 0.5 * (v_b_loc - v_a_loc),
//                 0.0, 0.0, 0.0,
//                 0.0, 0.0, 0.0,
//                 0.0, 0.0, 0.0,
//                 coupling_parameter_, wavenumber_,
//                 kernel);

//                 solution += kernel * muwi[i] * muwj[j] * psi_overs_loc[i * Nv_prec_ + j];

//             }
//         }

//     }

// }

// void Solver::write_psi_M1(const VectorXcd& phi, std::vector<double>& psi)
// {

//     #pragma omp parallel for
//     for (long long i = 0; i < patch_num_coeffs_.size(); i++) {

//         const long long patch_num = patch_num_coeffs_[i];                      

//         for (int nu = 0; nu < Nu_int_; nu++) {
//             for (int nv = 0; nv < Nv_int_; nv++) {

//                 const long long position = patch_num * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;
//                 const long long position_real = i * 2 * Nu_int_*Nv_int_ + nu * Nv_int_ + nv;
//                 const long long position_imag = i * 2 * Nu_int_*Nv_int_ + Nu_int_*Nv_int_ + nu * Nv_int_ + nv;

//                 psi[position_real] = phi[position].real() * dsdtjac_all_[position];
//                 psi[position_imag] = phi[position].imag() * dsdtjac_all_[position];

//             }
//         }

//     }

// }

// void Solver::write_psi_M2(const long long patch_num, const std::complex<double>* phi, std::vector<double>& psi)
// {

//     for (int nu = 0; nu < Nu_int_; nu++) {
//         for (int nv = 0; nv < Nv_int_; nv++) {

//             const int position = nu * Nv_int_ + nv;
//             const int position_real = nu * Nv_int_ + nv;
//             const int position_imag = Nu_int_*Nv_int_ + nu * Nv_int_ + nv;

//             psi[position_real] = phi[position].real() * dsdtjac_all_[patch_num * Nu_int_*Nv_int_ + position];
//             psi[position_imag] = phi[position].imag() * dsdtjac_all_[patch_num * Nu_int_*Nv_int_ + position];

//         }
//     }

// }

// void Solver::compute_integral_overs_M1(const VectorXcd& phi,
//                                 const std::vector<double>& psi_overs,
//                                 VectorXcd& integral)
// {
    
//     #pragma omp parallel for
//     for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

//         const double r_0 = disc_points_x_all_[point_low_ + npoint];
//         const double r_1 = disc_points_y_all_[point_low_ + npoint];
//         const double r_2 = disc_points_z_all_[point_low_ + npoint];

//         const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);
//         const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);

//         integral[npoint] = std::complex<double>(0.0, 0.0);

//         for (long long q = 0; q < Q_; q++) {             

//             for (int q_x = 0; q_x < Qx_; q_x++) {
//                 for (int q_y = 0; q_y < Qy_; q_y++) {

//                     const long long patch_num = q * Qx_*Qy_ + q_x * Qy_ + q_y;

//                     const auto patch_begin = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_begin);
//                     const auto patch_end = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_end);                            
//                     const auto pos = std::lower_bound(patch_begin, patch_end, patch_num);

//                     const auto idx = std::distance(patch_num_precomputations_.begin(), pos);
                    
//                     std::complex<double> solution;

//                     if ((point_precomputations_[idx] == (point_low_ + npoint)) && (patch_num_precomputations_[idx] == patch_num)) {

//                         const long long flag_u_loc = flags_domain_u_all_[patch_num];
//                         const long long flag_v_loc = flags_domain_v_all_[patch_num];

//                         const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
//                         const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
//                         const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
//                         const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

//                         const double ubar_loc = argmin_precomputations_u_[idx];
//                         const double vbar_loc = argmin_precomputations_v_[idx];

//                         const auto patch_begin_2 = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
//                         const auto idx_2 = std::distance(patch_num_coeffs_.begin(), patch_begin_2);
                        
//                         int_near_overs(r_0, r_1, r_2,
//                                         q, flag_u_loc, flag_v_loc,
//                                         u_a_loc, u_b_loc, v_a_loc, v_b_loc,
//                                         ubar_loc, vbar_loc,
//                                         &psi_overs[idx_2 * 2 * Nu_int_*Nv_int_ * P_over_*P_over_],
//                                         solution);

//                     } else {                            

//                         int_far(r_0, r_1, r_2,
//                                 patch_num, 
//                                 &phi[patch_num * Nu_int_*Nv_int_],
//                                 solution);
                        
//                     }

//                     integral[npoint] += solution;

//                 }
//             }

//         }

//         if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

//             integral[npoint] += 0.5 * phi[point_low_ + npoint];

//         }

//     }

// }

// void Solver::compute_integral_overs_M2(const VectorXcd& phi,
//                                 VectorXcd& integral)
// {

//     integral.setZero();   

//     #pragma omp parallel for
//     for (long long q = 0; q < Q_; q++) {

//         for (int q_x = 0; q_x < Qx_; q_x++) {
//             for (int q_y = 0; q_y < Qy_; q_y++) {

//                 const long long patch_num = q * Qx_*Qy_ + q_x * Qy_ + q_y;

//                 const int flag_u_loc = flags_domain_u_all_[patch_num];
//                 const int flag_v_loc = flags_domain_v_all_[patch_num];

//                 const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
//                 const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
//                 const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
//                 const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

//                 const auto patch_begin_2 = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
//                 const auto idx_2 = std::distance(patch_num_coeffs_.begin(), patch_begin_2);
                
//                 std::vector<double> psi(2 * Nu_int_*Nv_int_);            
//                 std::vector<double> psi_overs;
                
//                 #pragma omp critical
//                 {

//                 if (patch_num_coeffs_[idx_2] == patch_num) {

//                     write_psi_M2(patch_num, &phi[patch_num * Nu_int_*Nv_int_], psi);
                    
//                     compute_oversampling_M2(psi, psi_overs);

//                 }  

//                 }      
                
//                 for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

//                     const double r_0 = disc_points_x_all_[point_low_ + npoint];
//                     const double r_1 = disc_points_y_all_[point_low_ + npoint];
//                     const double r_2 = disc_points_z_all_[point_low_ + npoint];

//                     const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);
//                     const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point_low_ + npoint);

//                     const auto patch_begin = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_begin);
//                     const auto patch_end = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_end);                            
//                     const auto pos = std::lower_bound(patch_begin, patch_end, patch_num);

//                     const auto idx = std::distance(patch_num_precomputations_.begin(), pos);
                    
//                     std::complex<double> solution;

//                     if ((point_precomputations_[idx] == (point_low_ + npoint)) && (patch_num_precomputations_[idx] == patch_num)) {

//                         const double ubar_loc = argmin_precomputations_u_[idx];
//                         const double vbar_loc = argmin_precomputations_v_[idx];
                        
//                         int_near_overs(r_0, r_1, r_2,
//                                         q, flag_u_loc, flag_v_loc,
//                                         u_a_loc, u_b_loc, v_a_loc, v_b_loc,
//                                         ubar_loc, vbar_loc,
//                                         &psi_overs[0],
//                                         solution);

//                     } else {                            

//                         int_far(r_0, r_1, r_2,
//                                 patch_num, 
//                                 &phi[patch_num * Nu_int_*Nv_int_],
//                                 solution);
                        
//                     }

//                     integral[npoint] += solution;

//                 }

//             }
//         }

//     } 

//     if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

//         #pragma omp parallel for
//         for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

//             integral[npoint] += 0.5 * phi[point_low_ + npoint];

//         }

//     }  

// }

// void Solver::compute_integral_overs_acc_M1(const VectorXcd& phi,
//                                     const std::vector<double>& psi_overs,
//                                     VectorXcd& integral)
// {          

//     integral.setZero();

//     #pragma omp parallel for
//     for (long long i = 0; i < point_precomputations_.size(); i++) {

//         const long long npoint = point_precomputations_[i];
//         const long long patch_num = patch_num_precomputations_[i];

//         const double r_0 = disc_points_x_all_[npoint];
//         const double r_1 = disc_points_y_all_[npoint];
//         const double r_2 = disc_points_z_all_[npoint];

//         const int flag_u_loc = flags_domain_u_all_[patch_num];
//         const int flag_v_loc = flags_domain_v_all_[patch_num];

//         const long long q = patch_num / (Qx_*Qy_);
//         const int q_x = (patch_num % (Qx_*Qy_)) / Qy_;
//         const int q_y = (patch_num % (Qx_*Qy_)) % Qy_;

//         const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
//         const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
//         const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
//         const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

//         const double ubar_loc = argmin_precomputations_u_[i];
//         const double vbar_loc = argmin_precomputations_v_[i];

//         const auto patch_begin = std::lower_bound(patch_num_coeffs_.begin(), patch_num_coeffs_.end(), patch_num);
//         const auto idx_2 = std::distance(patch_num_coeffs_.begin(), patch_begin);

//         std::complex<double> solution;

//         int_near_overs(r_0, r_1, r_2,
//                         q, flag_u_loc, flag_v_loc,
//                         u_a_loc, u_b_loc, v_a_loc, v_b_loc,
//                         ubar_loc, vbar_loc,
//                         &psi_overs[idx_2 * 2 * Nu_int_*Nv_int_ * P_over_*P_over_],
//                         solution);

//         integral[npoint - point_low_] += solution;

//     }

//     if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

//         #pragma omp parallel for
//         for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

//             integral[npoint] += 0.5 * phi[point_low_ + npoint];

//         }

//     }        

// }

// void Solver::compute_integral_overs_acc_M2(const VectorXcd& phi,
//                                     VectorXcd& integral)
// {

//     integral.setZero();           

//     #pragma omp parallel for 
//     for (long long patch_num : patch_num_coeffs_) {

//         const long long q = patch_num / (Qx_*Qy_);
//         const int q_x = (patch_num % (Qx_*Qy_)) / Qy_;
//         const int q_y = (patch_num % (Qx_*Qy_)) % Qy_;

//         const int flag_u_loc = flags_domain_u_all_[patch_num];
//         const int flag_v_loc = flags_domain_v_all_[patch_num];

//         const double u_a_loc = -1.0 + q_x * 2.0 / Qx_;
//         const double u_b_loc = -1.0 + (q_x + 1) * 2.0 / Qx_;
//         const double v_a_loc = -1.0 + q_y * 2.0 / Qy_;
//         const double v_b_loc = -1.0 + (q_y + 1) * 2.0 / Qy_;

//         std::vector<double> psi(2 * Nu_int_*Nv_int_);            
//         std::vector<double> psi_overs;

//         #pragma omp critical
//         {

//         write_psi_M2(patch_num, &phi[patch_num * Nu_int_*Nv_int_], psi);

//         compute_oversampling_M2(psi, psi_overs);

//         }

//         std::vector<long long> points; // = patch_num_FFTW_[patch_num];

//         for (long long npoint : points) {

//             const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), npoint);
//             const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), npoint);
//             const auto patch_begin = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_begin);
//             const auto patch_end = patch_num_precomputations_.begin() + std::distance(point_precomputations_.begin(), point_end);                            
//             const auto pos = std::lower_bound(patch_begin, patch_end, patch_num);
//             const auto idx = std::distance(patch_num_precomputations_.begin(), pos);
                    
//             std::complex<double> solution;

//             const double r_0 = disc_points_x_all_[npoint];
//             const double r_1 = disc_points_y_all_[npoint];
//             const double r_2 = disc_points_z_all_[npoint];

//             const double ubar_loc = argmin_precomputations_u_[idx];
//             const double vbar_loc = argmin_precomputations_v_[idx];

//             int_near_overs(r_0, r_1, r_2,
//                             q, flag_u_loc, flag_v_loc,
//                             u_a_loc, u_b_loc, v_a_loc, v_b_loc,
//                             ubar_loc, vbar_loc,
//                             &psi_overs[0],
//                             solution);

//             integral[npoint - point_low_] += solution;

//         }

//     }

//     if (EQUATION_FORMULATION == 2 || EQUATION_FORMULATION == 3) {

//         #pragma omp parallel for
//         for (long long npoint = 0; npoint < point_up_-point_low_; npoint++) {

//             integral[npoint] += 0.5 * phi[point_low_ + npoint];

//         }

//     }   

// }

// void Solver::iterator_function_overs(const VectorXcd& phi,
//                                 VectorXcd& rhs)
// {
                
//     if (METHOD_FFTW == 1) {

//         std::vector<double> psi(2 * patch_num_coeffs_.size() * Nu_int_*Nv_int_);            
//         std::vector<double> psi_overs;
    
//         write_psi_M1(phi, psi);
        
//         compute_oversampling_M1(psi, psi_overs);

//         compute_integral_overs_M1(phi, psi_overs, rhs);   

//     } else {

//         compute_integral_overs_M2(phi, rhs);  

//     }   

// }

// void Solver::iterator_function_overs_acc(const VectorXcd& phi,
//                                     VectorXcd& rhs)
// {
    
//     if (METHOD_FFTW == 1) {

//         std::vector<double> psi(2 * patch_num_coeffs_.size() * Nu_int_*Nv_int_);            
//         std::vector<double> psi_overs;
    
//         write_psi_M1(phi, psi);

//         compute_oversampling_M1(psi, psi_overs);

//         compute_integral_overs_acc_M1(phi, psi_overs, rhs);        

//     } else {

//         compute_integral_overs_acc_M2(phi, rhs);  

//     }               

//     std::vector<std::complex<double>> intensities; 

//     compute_intensities(phi, intensities);

//     std::vector<std::complex<double>> solution_1;

//     solution_1 = intensities;

//     // boxes_.Solve<&fct_4, &fac_1>(solution_1); 
//         boxes_.Solve(
//     solution_1,
//     [&](double a1,double a2,double a3,double a4,double a5,double a6,
//         double a7,double a8,double a9,double a10,
//         std::complex<double> c1,std::complex<double> c2,std::complex<double>& out) {
//         Solver::fct_4(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,c1,c2,out);
//     },
//     [&](double a1,std::complex<double> c1,std::complex<double>& out) {
//         Solver::fac_1(a1,c1,out);
//     }
// ); 

//     #pragma omp parallel for
//     for (long long i = 0; i < point_up_-point_low_; i++) {

//         rhs[i] += solution_1[new_order_points_RP_[i]];

//     }  
    
// }

// void Solver::iterator_function(const VectorXcd& phi,
//                         VectorXcd& rhs)
// {

//     VectorXcd rhs_loc(point_up_ - point_low_);

//     if (USE_ACCELERATOR) {

//         if (USE_OVERSAMPLING) {

//             iterator_function_overs_acc(phi, rhs_loc);

//         } else {

//             iterator_function_acc(phi, rhs_loc);

//         }

//     } else {            
    
//         if (USE_OVERSAMPLING) {

//             iterator_function_overs(phi, rhs_loc);

//         } else {

//             iterator_function_unacc(phi, rhs_loc);

//         }

//     }

//     rhs = VectorXcd(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

//     MPI_Allgatherv(&rhs_loc[0], point_up_-point_low_, MPI_DOUBLE_COMPLEX, &rhs[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE_COMPLEX, mpi_comm_);

// }

// void Solver::setup(bool timing) 
// {

//     MPI_Barrier(mpi_comm_);

//     double start_1 = MPI_Wtime();


//     compute_parallel_parameters();

//     double end_1 = MPI_Wtime();

//     if (timing && comm_rank_ == 0) {

//         std::cout << "Time compute parallel parameters: "  << end_1 - start_1 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);

//     double start_2 = MPI_Wtime();

//     if (GEOMETRY != 0) {
//         load_interpolated_surface();
        

//     }

//     double end_2 = MPI_Wtime();

//     if (GEOMETRY != 0 && timing && comm_rank_ == 0) {
//         std::cout << "Time load interpolated surface: " << end_2 - start_2 << " seconds\n";
       
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);

//     double start_4 = MPI_Wtime();

//     compute_fejer_nodes_and_weights();

//     double end_4 = MPI_Wtime();

//     if (timing && comm_rank_ == 0) {

//         std::cout << "Time compute Fejer nodes and weights: " << end_4 - start_4 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);

//     double start_5 = MPI_Wtime();

//     if (!USE_OVERSAMPLING) {

//         compute_chebyshev_evaluations();

//     }

//     double end_5 = MPI_Wtime();

//     if (!USE_OVERSAMPLING && timing && comm_rank_ == 0) {

//         std::cout << "Time compute Chebyshev evaluations: " << end_5 - start_5 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);

//     double start_6 = MPI_Wtime();

//     compute_flags_domain();

//     double end_6 = MPI_Wtime();

//     if (timing && comm_rank_ == 0) {

//         std::cout << "Time compute flags domain: " << end_6 - start_6 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);

//     double start_7 = MPI_Wtime();

//     compute_discretization_domain();

//     double end_7 = MPI_Wtime();

//     if (timing && comm_rank_ == 0) {

//         std::cout << "Time compute discretization domain: " << end_7 - start_7 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);

//     double start_8 = MPI_Wtime();

//     compute_coupling_parameter();  

//     double end_8 = MPI_Wtime();

//     if (timing && comm_rank_ == 0) {

//         std::cout << "Time compute coupling parameter: " << end_8 - start_8 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);   

//     double start_9 = MPI_Wtime();   

//     compute_near_singular_patches_estimate();   

//     double end_9 = MPI_Wtime();

//     if (timing && comm_rank_ == 0) {

//         std::cout << "Time compute near singular patches estimate: " << end_9 - start_9 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);   

//     double start_10 = MPI_Wtime();

//     if (USE_OVERSAMPLING) {

//         compute_singular_points();

//     } else {

//         compute_precomputations();

//     }

//     double end_10 = MPI_Wtime();

//     if (timing && comm_rank_ == 0) {

//         std::cout << "Time compute precomputations data: " << end_10 - start_10 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);  

//     double start_11 = MPI_Wtime();

//     if (USE_ACCELERATOR) {

//         setup_IFGF_choose_levels();

//         if (comm_rank_ == 0) {
//             std::cout << "The number of levels chosen by the algorithm is " << nlevels_ << std::endl;
//         }

//         // create_IFGF_object();

//     }

//     double end_11 = MPI_Wtime();

//     if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

//         std::cout << "Time create IFGF object: " << end_11 - start_11 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);  


//     // double start_12 = MPI_Wtime();

//     // if (USE_ACCELERATOR) {

//     //     compute_new_order_points_RP();                

//     // }

//     // double end_12 = MPI_Wtime();

//     // if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

//     //     std::cout << "Time compute new order points RP: " << end_12 - start_12 << " seconds\n";
//     //     print_max_RSS();

//     // }

//     // MPI_Barrier(mpi_comm_);  


//     // double start_13 = MPI_Wtime();

//     // if (USE_ACCELERATOR) {

//     //     check_patch_in_neighbours();

//     // }

//     // double end_13 = MPI_Wtime();

//     // if (USE_ACCELERATOR && timing && comm_rank_ == 0) {

//     //     std::cout << "Time check patch in neighbours: " << end_13 - start_13 << " seconds\n";
//     //     print_max_RSS();

//     // }

//     // MPI_Barrier(mpi_comm_); 


//     double start_14 = MPI_Wtime();

//     if (USE_ACCELERATOR && USE_HIGH_ORDER) {

//         initialize_indexes_HO();

//     }

//     double end_14 = MPI_Wtime();

//     if (USE_ACCELERATOR && USE_HIGH_ORDER && timing && comm_rank_ == 0) {

//         std::cout << "Time initialize indexes HO: " << end_14 - start_14 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_); 

    
//     double start_15 = MPI_Wtime();

//     if (USE_OVERSAMPLING) {

//         create_fftw_objects();

//     }

//     double end_15 = MPI_Wtime();

//     if (USE_OVERSAMPLING && timing && comm_rank_ == 0) {

//         std::cout << "Time create FFTW objects: " << end_15 - start_15 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_); 

// }

// VectorXcd Solver::solve(const VectorXcd& rhs) 
// {

//     long long N = rhs.size();            
    
//     auto A = [this](const VectorXcd& x, VectorXcd& solution) -> void {iterator_function(x, solution);};

//     LinearOperator B(N, N, A);

//     opGMRES<LinearOperator> func(B);
//     func.setMaxIterations(MAX_ITER);
//     func.set_restart(MAX_ITER);
//     func.setTolerance(TOL_GMRES);
//     func.set_mpi_comm(mpi_comm_);

//     VectorXcd x = func.solve(rhs);

//     int rank;
//     MPI_Comm_rank(mpi_comm_, &rank);
//     if (rank == 0) {
//         std::cout << "GMRES converged in " << func.iterations() << " iterations.\n";
//     }

//     return x;

// }

// VectorXcd Solver::compute_incident_field() 
// {

//     std::vector<std::complex<double>> u_inc(point_up_ - point_low_, 0);

//     const double k_hat_0 = wavenumber_.real() * std::cos(0.0) * std::sin(M_PI);
//     const double k_hat_1 = wavenumber_.real() * std::sin(0.0) * std::sin(M_PI);
//     const double k_hat_2 = wavenumber_.real() * std::cos(M_PI);
    
//     #pragma omp parallel for
//     for (long long i = 0; i < point_up_ - point_low_; i++) {

//         const double x_0 = disc_points_x_all_[i + point_low_];
//         const double x_1 = disc_points_y_all_[i + point_low_];
//         const double x_2 = disc_points_z_all_[i + point_low_];               

//         const double inner_product = x_0 * k_hat_0 + x_1 * k_hat_1 + x_2 * k_hat_2;

//         const std::complex<double> value((-1.0) * std::cos(inner_product), (-1.0) * std::sin(inner_product));

//         u_inc[i] = value;

//     }

//     VectorXcd u_inc_all(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

//     MPI_Allgatherv(&u_inc[0], point_up_-point_low_, MPI_DOUBLE_COMPLEX, &u_inc_all[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE_COMPLEX, mpi_comm_);

//     return u_inc_all;

// }     
    
// //     VectorXcd Solver::compute_incident_field() 
// // {

// //     std::vector<std::complex<double>> u_inc(point_up_ - point_low_, 0);

// //     if (PLANE_OR_POINT == 0) {

// //         const double k_hat_0 = WAVE_NUMBER * std::cos(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
// //         const double k_hat_1 = WAVE_NUMBER * std::sin(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
// //         const double k_hat_2 = WAVE_NUMBER * std::cos(PLANE_WAVE_PHI);
        
// //         #pragma omp parallel for
// //         for (long long i = 0; i < point_up_ - point_low_; i++) {

// //             const double x_0 = disc_points_x_all_[i + point_low_];
// //             const double x_1 = disc_points_y_all_[i + point_low_];
// //             const double x_2 = disc_points_z_all_[i + point_low_];               

// //             const double inner_product = x_0 * k_hat_0 + x_1 * k_hat_1 + x_2 * k_hat_2;

// //             const std::complex<double> value((-1.0) * std::cos(inner_product), (-1.0) * std::sin(inner_product));

// //             u_inc[i] = value;

// //         }

// //     } else {

// //         const double k = WAVE_NUMBER;
        
// //         #pragma omp parallel for
// //         for (long long j = 0; j < point_up_ - point_low_; j++) {

// //             const double y_0 = disc_points_x_all_[j + point_low_];
// //             const double y_1 = disc_points_y_all_[j + point_low_];
// //             const double y_2 = disc_points_z_all_[j + point_low_];  

// //             for (int i = 0; i < NUM_POINT_SOURCES; i++) {

// //                 const double x_0 = POINT_SOURCE_CENTER[i][0];
// //                 const double x_1 = POINT_SOURCE_CENTER[i][1];
// //                 const double x_2 = POINT_SOURCE_CENTER[i][2];                                 

// //                 const double diff = std::sqrt((x_0-y_0)*(x_0-y_0) + (x_1-y_1)*(x_1-y_1) + (x_2-y_2)*(x_2-y_2));

// //                 const std::complex<double> value((-1.0) * std::cos(k *diff) / diff, (-1.0) * std::sin(k * diff) / diff);

// //                 u_inc[j] += value;

// //             }

// //         }

// //     }

// //     VectorXcd u_inc_all(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

// //     MPI_Allgatherv(&u_inc[0], point_up_-point_low_, MPI_DOUBLE_COMPLEX, &u_inc_all[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE_COMPLEX, mpi_comm_);

// //     return u_inc_all;

// // }     
        
    

// VectorXcd Solver::solve_u_inc(bool timing) 
// {
    
//     MPI_Barrier(mpi_comm_);


//     double start_1 = MPI_Wtime();

//     VectorXcd rhs = compute_incident_field();

//     double end_1 = MPI_Wtime();

//     if (timing && comm_rank_ == 0) {

//         std::cout << "Time compute incident field: " << end_1 - start_1 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);

    
//     double start_2 = MPI_Wtime();

//     VectorXcd solution = solve(rhs);

//     double end_2 = MPI_Wtime();

//     if (timing && comm_rank_ == 0) {

//         std::cout << "Time solve: " << end_2 - start_2 << " seconds\n";
//         print_max_RSS();

//     }

//     MPI_Barrier(mpi_comm_);


//     return solution;

// }

// void Solver::int_far_field(const double xVers_0, const double xVers_1, const double xVers_2,
//                     const long long npatch, 
//                     const std::complex<double>* phi,
//                     std::complex<double>& solution)                      
// {

//     solution = std::complex<double>(0.0, 0.0);

//     for (int i = 0; i < Nu_int_; i++) {
//         for (int j = 0; j < Nv_int_; j++) {

//             const double dsdtjac_loc = dsdtjac_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double constant = dsdtjac_loc * fejer_weights_u_int_[i] * fejer_weights_v_int_[j];

//             const double px = disc_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double py = disc_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double pz = disc_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

//             const double nx = norm_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double ny = norm_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double nz = norm_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

//             std::complex<double> kernel;
//             HH_far(xVers_0, xVers_1, xVers_2, px, py, pz, nx, ny, nz, coupling_parameter_, wavenumber_.real(),kernel);

//             solution += constant * kernel * phi[i*Nv_int_ + j];

//         }
//     }

// }

// std::complex<double> Solver::compute_far_field_approx(const VectorXcd& phi, 
//                                                 const double xVers_0, const double xVers_1, const double xVers_2) 
// {
    
//     std::complex<double> solution_rank(0.0, 0.0);

//     for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

//         std::complex<double> solution_loc;

//         int_far_field(xVers_0, xVers_1, xVers_2,
//                         patch_num,
//                         &phi[patch_num * Nu_int_*Nv_int_],
//                         solution_loc);                        
        
//         solution_rank += solution_loc;

//     }

//     solution_rank *= 1.0 / (4.0 * M_PI);

//     return solution_rank;

// }

// std::complex<double> Solver::compute_far_field_exact(const double xVers_0, const double xVers_1, const double xVers_2) 
// {

//     std::complex<double> solution(0.0, 0.0);

//     const std::complex<double> I(0.0, 1.0);

//     int nterms = 0;

//     std::complex<double> yn = std::sph_neumann(nterms, wavenumber_.real() * SPHERE_RADIUS);
//     std::complex<double> jn = std::sph_bessel(nterms, wavenumber_.real() * SPHERE_RADIUS);
//     std::complex<double> hn = jn + I * yn;

//     while (std::abs(hn) * std::abs(hn) < 1.0e14) {

//         nterms++;

//         yn = std::sph_neumann(nterms, wavenumber_.real() * SPHERE_RADIUS);
//         jn = std::sph_bessel(nterms, wavenumber_.real() * SPHERE_RADIUS);
//         hn = jn + I * yn;

//     }

//     const double kVers[3] = {std::cos(0.0) * std::sin(M_PI), std::sin(0.0) * std::sin(M_PI), std::cos(M_PI)};
//     const double x = xVers_0 * kVers[0] + xVers_1 * kVers[1] + xVers_2 * kVers[2];

//     for (int n = 0; n <= nterms; n++) {

//         const std::complex<double> yn = std::sph_neumann(n, wavenumber_.real() * SPHERE_RADIUS);
//         const std::complex<double> jn = std::sph_bessel(n, wavenumber_.real() * SPHERE_RADIUS);
//         const std::complex<double> hn = jn + I * yn;

//         const double P = std::legendre(n, x);

//         solution += (2.0 * n + 1.0) * jn/hn * P;

//     }

//     solution *= I / (wavenumber_.real() * SPHERE_RADIUS);

//     return solution;

// }

// void Solver::compute_far_field_error(const bool timing, const VectorXcd& phi)
// {

//     MPI_Barrier(mpi_comm_);
//     double start = MPI_Wtime();

//     const double deltaPhi = M_PI / (200 - 1);
//     const double deltaTheta = 2.0 * M_PI / (200 - 1);

//     const long long total_pts = 200 * 200;

//     std::vector<long long> split_points_far_field(comm_size_ + 1);
//     const long long points_per_rank = total_pts / comm_size_;
//     long long remaining_points = total_pts % comm_size_;

//     split_points_far_field[0] = 0;
//     split_points_far_field[comm_size_] = total_pts;

//     for (int i = 1; i < comm_size_; i++) {

//         split_points_far_field[i] = split_points_far_field[i-1] + points_per_rank;

//         if (remaining_points > 0) {

//             split_points_far_field[i]++;
//             remaining_points--;

//         }

//     }

//     std::vector<int> recv_counts_far_field(comm_size_);
//     std::vector<int> displs_far_field(comm_size_, 0);

//     for (int i = 0; i < comm_size_; i++) {

//         recv_counts_far_field[i] = split_points_far_field[i+1] - split_points_far_field[i];

//         if (i != 0) {

//             displs_far_field[i] = displs_far_field[i-1] + recv_counts_far_field[i-1];

//         }

//     }

//     long long point_low_far_field = split_points_far_field[comm_rank_];
//     long long point_up_far_field = split_points_far_field[comm_rank_+1];

//     std::vector<double> xVers_0_loc(point_up_far_field-point_low_far_field);
//     std::vector<double> xVers_1_loc(point_up_far_field-point_low_far_field);
//     std::vector<double> xVers_2_loc(point_up_far_field-point_low_far_field);

//     std::vector<std::complex<double>> far_field_exact_loc(point_up_far_field-point_low_far_field);

//     #pragma omp parallel for
//     for (long long point = 0; point < point_up_far_field-point_low_far_field; point++) {

//         const long long m = (point + point_low_far_field) / 200;
//         const long long n = (point + point_low_far_field) % 200;

//         const double phi_m = (m - 1) * deltaPhi;
//         const double theta_n = (n - 1) * deltaTheta;

//         xVers_0_loc[point] = std::sin(phi_m) * std::cos(theta_n);
//         xVers_1_loc[point] = std::sin(phi_m) * std::sin(theta_n);
//         xVers_2_loc[point] = std::cos(phi_m);

//         far_field_exact_loc[point] = compute_far_field_exact(xVers_0_loc[point], xVers_1_loc[point], xVers_2_loc[point]);

//     }

//     std::vector<double> xVers_0_all(total_pts);
//     std::vector<double> xVers_1_all(total_pts);
//     std::vector<double> xVers_2_all(total_pts);
    
//     MPI_Allgatherv(&xVers_0_loc[0], point_up_far_field-point_low_far_field, MPI_DOUBLE, &xVers_0_all[0], &recv_counts_far_field[0], &displs_far_field[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&xVers_1_loc[0], point_up_far_field-point_low_far_field, MPI_DOUBLE, &xVers_1_all[0], &recv_counts_far_field[0], &displs_far_field[0], MPI_DOUBLE, mpi_comm_);
//     MPI_Allgatherv(&xVers_2_loc[0], point_up_far_field-point_low_far_field, MPI_DOUBLE, &xVers_2_all[0], &recv_counts_far_field[0], &displs_far_field[0], MPI_DOUBLE, mpi_comm_);

//     std::vector<double>().swap(xVers_0_loc);
//     std::vector<double>().swap(xVers_1_loc);
//     std::vector<double>().swap(xVers_2_loc);

//     std::vector<std::complex<double>> far_field_approx_loc(total_pts);

//     #pragma omp parallel for
//     for (long long point = 0; point < total_pts; point++) {
        
//         far_field_approx_loc[point] = compute_far_field_approx(phi, xVers_0_all[point], xVers_1_all[point], xVers_2_all[point]);

//     }

//     std::vector<std::complex<double>> far_field_approx_all(total_pts);

//     MPI_Allreduce(&far_field_approx_loc[0], &far_field_approx_all[0], total_pts, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_);

//     std::vector<std::complex<double>>().swap(far_field_approx_loc);

//     double error_1_loc = 0.0;
//     double error_2_loc = 0.0;

//     for (long long point = point_low_far_field; point < point_up_far_field; point++) {

//         const std::complex<double> far_field_approx = far_field_approx_all[point];
//         const std::complex<double> far_field_exact = far_field_exact_loc[point - point_low_far_field];
        
//         const double far_field_approx_norm = std::sqrt(far_field_approx.real()*far_field_approx.real() + far_field_approx.imag()*far_field_approx.imag());
//         const double far_field_exact_norm = std::sqrt(far_field_exact.real()*far_field_exact.real() + far_field_exact.imag()*far_field_exact.imag());
            
//         const double value1 = std::abs(far_field_exact_norm - far_field_approx_norm);
//         const double value2 = std::abs(far_field_exact_norm);

//         if (error_1_loc < value1) error_1_loc = value1;
//         if (error_2_loc < value2) error_2_loc = value2;

//     }
        
//     double error_1;
//     double error_2;

//     MPI_Allreduce(&error_1_loc, &error_1, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
//     MPI_Allreduce(&error_2_loc, &error_2, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);

//     double end = MPI_Wtime();
//     MPI_Barrier(mpi_comm_);

//     if (comm_rank_ == 0) {

//         if (timing) {

//             std::cout << "Time compute far field error: " << end - start << " seconds\n";
//             print_max_RSS();

//         }

//         std::cout << "Far field error: " << error_1 / error_2 << "\n";

//     }

// }

//     // /*std::complex<double> compute_incident_field(double x, double y, double z) 
//     // {

//     //     std::complex<double> u_inc = {0.0, 0.0};

//     //     if (PLANE_OR_POINT == 0) {

//     //         const double k_hat_0 = WAVE_NUMBER * std::cos(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
//     //         const double k_hat_1 = WAVE_NUMBER * std::sin(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
//     //         const double k_hat_2 = WAVE_NUMBER * std::cos(PLANE_WAVE_PHI);            

//     //         const double inner_product = x * k_hat_0 + y * k_hat_1 + z * k_hat_2;

//     //         const std::complex<double> value((-1.0) * std::cos(inner_product), (-1.0) * std::sin(inner_product));

//     //         u_inc = value;

//     //     } else {

//     //         const double k = WAVE_NUMBER;

//     //         for (int i = 0; i < NUM_POINT_SOURCES; i++) {

//     //             const double x_0 = POINT_SOURCE_CENTER[i][0];
//     //             const double x_1 = POINT_SOURCE_CENTER[i][1];
//     //             const double x_2 = POINT_SOURCE_CENTER[i][2];                                 

//     //             const double diff = std::sqrt((x_0-x)*(x_0-x) + (x_1-y)*(x_1-y) + (x_2-z)*(x_2-z));

//     //             const std::complex<double> value((-1.0) * std::cos(k * diff) / diff, (-1.0) * std::sin(k * diff) / diff);

//     //             u_inc += value;

//     //         }

//     //     }

//     //     return u_inc;

//     // }  */

// void Solver::int_near_field(const double x_0, const double x_1, const double x_2,
//                     const long long npatch, 
//                     const std::complex<double>* phi,
//                     std::complex<double>& solution)                      
// {

//     solution = std::complex<double>(0.0, 0.0);

//     for (int i = 0; i < Nu_int_; i++) {
//         for (int j = 0; j < Nv_int_; j++) {

//             const double dsdtjac_loc = dsdtjac_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double constant = dsdtjac_loc * fejer_weights_u_int_[i] * fejer_weights_v_int_[j];

//             const double px = disc_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double py = disc_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double pz = disc_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

//             const double nx = norm_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double ny = norm_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
//             const double nz = norm_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

//             std::complex<double> kernel;
//             HH2(x_0, x_1, x_2, px, py, pz, nx, ny, nz, coupling_parameter_, wavenumber_, kernel);

//             solution += constant * kernel * phi[i*Nv_int_ + j];

//         }
//     }

// }

// std::complex<double> Solver::compute_near_field_approx(const VectorXcd& phi, 
//                                                 const double x_0, const double x_1, const double x_2) 
// {
    
//     std::complex<double> solution_rank(0.0, 0.0);

//     for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

//         std::complex<double> solution_loc;

//         int_near_field(x_0, x_1, x_2,
//                         patch_num,
//                         &phi[patch_num * Nu_int_*Nv_int_],
//                         solution_loc);                        
        
//         solution_rank += solution_loc;

//     }

//     return solution_rank;
// }
// /*
// Custom function to evaluate the near field, based off of a small change from compute
// near_field_approx
// */
// std::complex<double> Solver::compute_near_field(const VectorXcd& phi, 
//                                                 const double x_0, const double x_1, const double x_2) 
// {
    
//     std::complex<double> solution_rank(0.0, 0.0);
//     if (comm_rank_ == 0) {
//         #pragma omp parallel for
//         for (long long patch_num = 0; patch_num < Q_*Qx_*Qy_; patch_num++) {

//             std::complex<double> solution_loc;

//             int_near_field(x_0, x_1, x_2,
//                         patch_num,
//                         &phi[patch_num * Nu_int_*Nv_int_],
//                         solution_loc);                        
            
//             solution_rank += solution_loc;
//         }
//     }
//     // Make sure every process has the solution :)
//     MPI_Bcast(&solution_rank, 1, MPI_DOUBLE_COMPLEX, 0, mpi_comm_);

//     return solution_rank;

// }

//     // /*
//     // void Solver::compute_near_field_error(const bool timing, const VectorXcd& phi)
//     // {

//     //     double loc_max_x = std::numeric_limits<double>::lowest();
//     //     double loc_max_y = std::numeric_limits<double>::lowest();
//     //     double loc_max_z = std::numeric_limits<double>::lowest();

//     //     double loc_min_x = std::numeric_limits<double>::max();
//     //     double loc_min_y = std::numeric_limits<double>::max();
//     //     double loc_min_z = std::numeric_limits<double>::max();

//     //     std::vector<double> bb_x_min(patch_up_-patch_low_, std::numeric_limits<double>::max()), 
//     //                         bb_x_max(patch_up_-patch_low_, std::numeric_limits<double>::lowest()), 
//     //                         bb_y_min(patch_up_-patch_low_, std::numeric_limits<double>::max()), 
//     //                         bb_y_max(patch_up_-patch_low_, std::numeric_limits<double>::lowest()), 
//     //                         bb_z_min(patch_up_-patch_low_, std::numeric_limits<double>::max()), 
//     //                         bb_z_max(patch_up_-patch_low_, std::numeric_limits<double>::lowest());

//     //     for (long long i = patch_low_; i < patch_up_; i++) {

//     //         for (long long j = 0; j < Nu_int_; j++) {
//     //             for (long long k = 0; k < Nv_int_; k++) {

//     //                 const long long idx = i * Nu_int_*Nv_int_ + j * Nv_int_ + k;

//     //                 const double pointx = disc_points_x_all_[idx];
//     //                 const double pointy = disc_points_y_all_[idx];
//     //                 const double pointz = disc_points_z_all_[idx];

//     //                 if (loc_max_x < pointx) loc_max_x = pointx;
//     //                 if (loc_max_y < pointy) loc_max_y = pointy;
//     //                 if (loc_max_z < pointz) loc_max_z = pointz;
                    
//     //                 if (loc_min_x > pointx) loc_min_x = pointx;
//     //                 if (loc_min_y > pointy) loc_min_y = pointy;
//     //                 if (loc_min_z > pointz) loc_min_z = pointz;

//     //                 if (bb_x_max[i - patch_low_] < pointx) bb_x_max[i - patch_low_] = pointx;
//     //                 if (bb_y_max[i - patch_low_] < pointy) bb_y_max[i - patch_low_] = pointy;
//     //                 if (bb_z_max[i - patch_low_] < pointz) bb_z_max[i - patch_low_] = pointz;
                    
//     //                 if (bb_x_min[i - patch_low_] > pointx) bb_x_min[i - patch_low_] = pointx;
//     //                 if (bb_y_min[i - patch_low_] > pointy) bb_y_min[i - patch_low_] = pointy;
//     //                 if (bb_z_min[i - patch_low_] > pointz) bb_z_min[i - patch_low_] = pointz;

//     //             }
//     //         }  

//     //         const double Lx = bb_x_max[i - patch_low_] - bb_x_min[i - patch_low_];
//     //         const double Ly = bb_y_max[i - patch_low_] - bb_y_min[i - patch_low_];
//     //         const double Lz = bb_z_max[i - patch_low_] - bb_z_min[i - patch_low_];

//     //         const double Lmax = std::max({Lx, Ly, Lz});
//     //         const double factor = Lmax * 2.0E-2; 

//     //         bb_x_max[i - patch_low_] += factor;
//     //         bb_x_min[i - patch_low_] -= factor;
//     //         bb_y_max[i - patch_low_] += factor;
//     //         bb_y_min[i - patch_low_] -= factor;
//     //         bb_z_max[i - patch_low_] += factor;
//     //         bb_z_min[i - patch_low_] -= factor;

//     //     }

//     //     double global_max_x;
//     //     double global_max_y;
//     //     double global_max_z;

//     //     double global_min_x;
//     //     double global_min_y;
//     //     double global_min_z;

//     //     MPI_Allreduce(&loc_max_x, &global_max_x, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
//     //     MPI_Allreduce(&loc_max_y, &global_max_y, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
//     //     MPI_Allreduce(&loc_max_z, &global_max_z, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
//     //     MPI_Allreduce(&loc_min_x, &global_min_x, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_);
//     //     MPI_Allreduce(&loc_min_y, &global_min_y, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_);
//     //     MPI_Allreduce(&loc_min_z, &global_min_z, 1, MPI_DOUBLE, MPI_MIN, mpi_comm_);

//     //     std::vector<double> bb_x_min_all(Q_ * Qx_*Qy_), 
//     //                         bb_x_max_all(Q_ * Qx_*Qy_), 
//     //                         bb_y_min_all(Q_ * Qx_*Qy_), 
//     //                         bb_y_max_all(Q_ * Qx_*Qy_), 
//     //                         bb_z_min_all(Q_ * Qx_*Qy_), 
//     //                         bb_z_max_all(Q_ * Qx_*Qy_);

//     //     MPI_Allgatherv(&bb_x_min[0], patch_up_-patch_low_, MPI_DOUBLE, &bb_x_min_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
//     //     MPI_Allgatherv(&bb_y_min[0], patch_up_-patch_low_, MPI_DOUBLE, &bb_y_min_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
//     //     MPI_Allgatherv(&bb_z_min[0], patch_up_-patch_low_, MPI_DOUBLE, &bb_z_min_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
//     //     MPI_Allgatherv(&bb_x_max[0], patch_up_-patch_low_, MPI_DOUBLE, &bb_x_max_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
//     //     MPI_Allgatherv(&bb_y_max[0], patch_up_-patch_low_, MPI_DOUBLE, &bb_y_max_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);
//     //     MPI_Allgatherv(&bb_z_max[0], patch_up_-patch_low_, MPI_DOUBLE, &bb_z_max_all[0], &recv_counts_[0], &displs_[0], MPI_DOUBLE, mpi_comm_);

//     //     std::vector<double>().swap(bb_x_min);
//     //     std::vector<double>().swap(bb_y_min);
//     //     std::vector<double>().swap(bb_z_min);
//     //     std::vector<double>().swap(bb_x_max);
//     //     std::vector<double>().swap(bb_y_max);
//     //     std::vector<double>().swap(bb_z_max);

//     //     const int nNearZones = N_NEAR_PTS.size();

//     //     for (int i = 0; i < nNearZones; i++) {

//     //         MPI_Barrier(mpi_comm_);
//     //         double start = MPI_Wtime();

//     //         const double xmin = NEAR_FIELD_LIMITS[i][0];
//     //         const double xmax = NEAR_FIELD_LIMITS[i][1];
//     //         const double ymin = NEAR_FIELD_LIMITS[i][2];
//     //         const double ymax = NEAR_FIELD_LIMITS[i][3];
//     //         const double zmin = NEAR_FIELD_LIMITS[i][4];
//     //         const double zmax = NEAR_FIELD_LIMITS[i][5];

//     //         const int imax = N_NEAR_PTS[i][0];
//     //         const int jmax = N_NEAR_PTS[i][1];
//     //         const int kmax = N_NEAR_PTS[i][2];

//     //         double deltaX = (xmax - xmin) / (imax - 1);
//     //         double deltaY = (ymax - ymin) / (jmax - 1);
//     //         double deltaZ = (zmax - zmin) / (kmax - 1);

//     //         if (imax == 1) deltaX = 0.0;
//     //         if (jmax == 1) deltaY = 0.0;
//     //         if (kmax == 1) deltaZ = 0.0;

//     //         const long long total_pts = imax * jmax * kmax;

//     //         std::vector<long long> split_points_near_field(comm_size_ + 1);
//     //         const long long points_per_rank = total_pts / comm_size_;
//     //         long long remaining_points = total_pts % comm_size_;

//     //         split_points_near_field[0] = 0;
//     //         split_points_near_field[comm_size_] = total_pts;

//     //         for (int l = 1; l < comm_size_; l++) {

//     //             split_points_near_field[l] = split_points_near_field[l-1] + points_per_rank;

//     //             if (remaining_points > 0) {

//     //                 split_points_near_field[l]++;
//     //                 remaining_points--;

//     //             }

//     //         }

//     //         std::vector<int> recv_counts_near_field(comm_size_);
//     //         std::vector<int> displs_near_field(comm_size_, 0);

//     //         for (int l = 0; l < comm_size_; l++) {

//     //             recv_counts_near_field[l] = split_points_near_field[l+1] - split_points_near_field[l];

//     //             if (l != 0) {

//     //                 displs_near_field[l] = displs_near_field[l-1] + recv_counts_near_field[l-1];

//     //             }

//     //         }

//     //         long long point_low_near_field = split_points_near_field[comm_rank_];
//     //         long long point_up_near_field = split_points_near_field[comm_rank_+1];

//     //         std::vector<double> x_loc(point_up_near_field-point_low_near_field);
//     //         std::vector<double> y_loc(point_up_near_field-point_low_near_field);
//     //         std::vector<double> z_loc(point_up_near_field-point_low_near_field);

//     //         std::vector<int> mask_loc(point_up_near_field-point_low_near_field, 2);

//     //         #pragma omp parallel for
//     //         for (long long point = 0; point < point_up_near_field-point_low_near_field; point++) {

//     //             const long long ii = (point + point_low_near_field) / (jmax * kmax);
//     //             const long long jj = ((point + point_low_near_field) % (jmax * kmax)) / kmax;
//     //             const long long kk = ((point + point_low_near_field) % (jmax * kmax)) % kmax;

//     //             x_loc[point] = xmin + ii * deltaX;
//     //             y_loc[point] = ymin + jj * deltaY;
//     //             z_loc[point] = zmin + kk * deltaZ;

//     //             if (x_loc[point] < global_min_x ||
//     //                 x_loc[point] > global_max_x ||
//     //                 y_loc[point] < global_min_y ||
//     //                 y_loc[point] > global_max_y ||
//     //                 z_loc[point] < global_min_z ||
//     //                 z_loc[point] > global_max_z) {

//     //                     mask_loc[point] = 1;

//     //             }

//     //         }

//     //         if (imax == 1) {



//     //         }

//     //         if (jmax == 1) {



//     //         }

//     //         if (kmax == 1) {



//     //         }

//     //         std::vector<std::complex<double>> u_inc_loc(point_up_near_field-point_low_near_field);
//     //         std::vector<std::complex<double>> u_scat_loc(point_up_near_field-point_low_near_field);
//     //         std::vector<std::complex<double>> u_total_loc(point_up_near_field-point_low_near_field);

//     //         #pragma omp parallel for
//     //         for (long long point = 0; point < point_up_near_field-point_low_near_field; point++) {

//     //             if (mask_loc[point] == 2) {

//     //                 u_inc_loc[point] = {0.0, 0.0};
//     //                 u_scat_loc[point] = {0.0, 0.0};
//     //                 u_total_loc[point] = {0.0, 0.0};

//     //                 continue;

//     //             }

//     //             if (mask_loc[point] <= 0) {

//     //                 u_inc_loc[point] = {0.0, 0.0};
//     //                 u_scat_loc[point] = {0.0, 0.0};
//     //                 u_total_loc[point] = {0.0, 0.0};

//     //             } else {

//     //                 u_inc_loc[point] = -compute_incident_field(x_loc[point], y_loc[point], z_loc[point]);

//     //             }

//     //             u_scat_loc[point] = compute_near_field_approx(phi, x_loc[point], y_loc[point], z_loc[point]);

//     //             u_total_loc[point] = u_inc_loc[point] + u_scat_loc[point];

//     //         }  
            
//     //     }

//     // }*/

//     std::complex<double> Solver::iterator_function_Manuel(const VectorXcd& u, const VectorXcd& v) 
//     {

//         VectorXcd y = solve(v);

//         return (u.conjugate().transpose() * y);

//     }

// void Solver::init_solver(const bool timing, 
//         const std::complex<double> k, 
//         const int* n_pts_per_patch,
//         const int* n_split_per_patch,
//         const int* n_pts_sing_int,
//         const double proximity_box_size,
//         const int n_levels_IFGF, 
//         const MPI_Comm& mpi_comm) {
    
//     Nu_int_ = n_pts_per_patch[0];
//     Nv_int_ = n_pts_per_patch[1];

//     Nu_prec_ = n_pts_sing_int[0];
//     Nv_prec_ = n_pts_sing_int[1];

//     Qx_ = n_split_per_patch[0];
//     Qy_ = n_split_per_patch[1];

//     wavenumber_ = k;
//     lambda_ = 2.0 * M_PI / wavenumber_;

//     proximity_ = proximity_box_size;
//     nlevels_ = n_levels_IFGF;

//     mpi_comm_ = mpi_comm;



//     setup(timing);

//         }

// Solver::Solver(const bool timing, 
//         const std::complex<double> k, 
//         const int* n_pts_per_patch,
//         const int* n_split_per_patch,
//         const int* n_pts_sing_int,
//         const double proximity_box_size,
//         const int n_levels_IFGF, 
//         const MPI_Comm& mpi_comm, double sphere_radius,
//         double sphere_centerX, double sphere_centerY, double sphere_centerZ)
// {

//     SPHERE_RADIUS = sphere_radius;
//     SPHERE_CENTER[0] = sphere_centerX; SPHERE_CENTER[1] = sphere_centerY; SPHERE_CENTER[2] = sphere_centerZ;
//     // Q_ is the number of patches before splitting.
//     Q_ = 6;
//     EDGE_FLAG_U_A = std::vector<bool>(Q_, false);
//     EDGE_FLAG_U_B = std::vector<bool>(Q_, false);
//     EDGE_FLAG_V_A = std::vector<bool>(Q_, false);
//     EDGE_FLAG_V_B = std::vector<bool>(Q_, false);

//     init_solver(timing, k, n_pts_per_patch, n_split_per_patch, n_pts_sing_int, proximity_box_size,
//     n_levels_IFGF, mpi_comm);
// }

// Solver::Solver(const bool timing, 
//                const std::complex<double> k, 
//                const int* n_pts_per_patch,
//                const int* n_split_per_patch,
//                const int* n_pts_sing_int,
//                const double proximity_box_size,
//                const int n_levels_IFGF, 
//                const MPI_Comm& mpi_comm, 
//                const std::string directory,
//                const std::string file_prefix) 
// {

//     DIRECTORY = directory;
//     FILE_NAME = file_prefix;
//     Q_ = 0;

//     GEOMETRY = 1;

//      try {
//         // Create a directory_iterator for the specified path
//         std::filesystem::directory_iterator dirIter(DIRECTORY);

//         // Count regular files using std::count_if and a lambda expression
//         Q_ = std::count_if(begin(dirIter), end(dirIter),
//                                       [](const std::filesystem::directory_entry& entry) {
//                                           return entry.is_regular_file();
//                                       });

//     } catch (const std::filesystem::filesystem_error& ex) {
//         std::cerr << "Error accessing directory: " << ex.what() << std::endl;
        
//     }

//     EDGE_FLAG_U_A = std::vector<bool>(Q_, false);
//     EDGE_FLAG_U_B = std::vector<bool>(Q_, false);
//     EDGE_FLAG_V_A = std::vector<bool>(Q_, false);
//     EDGE_FLAG_V_B = std::vector<bool>(Q_, false);

//     init_solver(timing, k, n_pts_per_patch, n_split_per_patch, n_pts_sing_int, proximity_box_size,
//     n_levels_IFGF, mpi_comm);
               

// }

// Solver::~Solver() {

//     if (USE_OVERSAMPLING) {

//         fftw_destroy_plan(plan_original_patch_);
//         fftw_destroy_plan(plan_oversampled_patch_);

//         fftw_free(original_patch_);
//         fftw_free(oversampled_patch_);

//     }

// }    

