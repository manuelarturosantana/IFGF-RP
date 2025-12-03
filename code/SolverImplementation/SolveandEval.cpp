#include "../solver2.h"

 void Solver::compute_integral_acc(const std::complex<double>* phi,
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

void Solver::iterator_function_acc(std::complex<double>* phi, std::complex<double>* rhs)
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

void Solver::compute_phi_not_in_rank(const std::complex<double>* phi, std::unordered_map<long long, 
    std::vector<std::complex<double>>>& phi_not_in_rank)
 {

    std::vector<std::vector<long long>> patches_not_in_rank(comm_size_);

    for (long long patch_num : patch_num_coeffs_) {

        auto it = std::upper_bound(split_points_.begin(), split_points_.end(), patch_num);

        int rank = static_cast<int>(std::distance(split_points_.begin(), it)) - 1;

        if (rank < 0) rank = 0;                
        if (rank >= comm_size_) rank = comm_size_ - 1; 

        if (comm_rank_ != rank) {

            patches_not_in_rank[rank].push_back(patch_num);

        }

    }

    std::vector<MPI_Count> send_counts_patches(comm_size_);
    std::vector<MPI_Aint> sdispls_patches(comm_size_);
    MPI_Count total_send_patches = 0;    

    for (int i = 0; i < comm_size_; ++i) {

        send_counts_patches[i] = static_cast<MPI_Count>(patches_not_in_rank[i].size());
        sdispls_patches[i] = total_send_patches;
        total_send_patches += send_counts_patches[i];

    }

    std::vector<long long> flat_send_patches_buffer(static_cast<size_t>(total_send_patches));

    MPI_Aint current_patch_pos = 0;

    for (int i = 0; i < comm_size_; ++i) {

        std::copy(patches_not_in_rank[i].begin(), patches_not_in_rank[i].end(),
                    flat_send_patches_buffer.begin() + current_patch_pos);
        current_patch_pos += send_counts_patches[i];

    }
    
    patches_not_in_rank.clear();

    std::vector<MPI_Count> recv_counts_patches(comm_size_);

    MPI_Alltoall(send_counts_patches.data(), 1, MPI_COUNT, recv_counts_patches.data(), 1, MPI_COUNT, mpi_comm_);

    std::vector<MPI_Aint> rdispls_patches(comm_size_);
    MPI_Count total_recv_patches = 0;

    for (int i = 0; i < comm_size_; ++i) {

        rdispls_patches[i] = total_recv_patches;
        total_recv_patches += recv_counts_patches[i];

    }

    std::vector<long long> flat_recv_patches_buffer(static_cast<size_t>(total_recv_patches));

    MPI_Request request_patches;
    MPI_Ialltoallv_c(flat_send_patches_buffer.data(), send_counts_patches.data(), sdispls_patches.data(), MPI_LONG_LONG,
                        flat_recv_patches_buffer.data(), recv_counts_patches.data(), rdispls_patches.data(), MPI_LONG_LONG,
                        mpi_comm_, &request_patches);

    const size_t patch_data_size = Nu_int_ * Nv_int_; 
    MPI_Count total_send_data_back = total_recv_patches * static_cast<MPI_Count>(patch_data_size);

    std::vector<std::complex<double>> flat_send_data_back_buffer(static_cast<size_t>(total_send_data_back));

    std::vector<MPI_Count> send_counts_data_back(comm_size_);
    std::vector<MPI_Aint> sdispls_data_back(comm_size_);
    MPI_Aint current_data_pos = 0;

    for (int sender_rank = 0; sender_rank < comm_size_; ++sender_rank) {

        MPI_Count num_patches = recv_counts_patches[sender_rank];
        MPI_Count data_size = num_patches * static_cast<MPI_Count>(patch_data_size);

        send_counts_data_back[sender_rank] = data_size;
        sdispls_data_back[sender_rank] = current_data_pos;
        current_data_pos += data_size;

    }

    MPI_Wait(&request_patches, MPI_STATUS_IGNORE); 

    MPI_Aint current_patch_idx = 0;

    for (int sender_rank = 0; sender_rank < comm_size_; ++sender_rank) {

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

    std::vector<MPI_Count> recv_counts_data_back(comm_size_);

    MPI_Alltoall(send_counts_data_back.data(), 1, MPI_COUNT,
                    recv_counts_data_back.data(), 1, MPI_COUNT, mpi_comm_);

    MPI_Count total_recv_data_back = 0;
    std::vector<MPI_Aint> rdispls_data_back(comm_size_);

    for (int i = 0; i < comm_size_; ++i) {

        rdispls_data_back[i] = total_recv_data_back;
        total_recv_data_back += recv_counts_data_back[i];

    }

    std::vector<std::complex<double>> flat_recv_data_buffer_back(static_cast<size_t>(total_recv_data_back));

    MPI_Request request_data;
    MPI_Ialltoallv_c(flat_send_data_back_buffer.data(), send_counts_data_back.data(), sdispls_data_back.data(), MPI_DOUBLE_COMPLEX,
                        flat_recv_data_buffer_back.data(), recv_counts_data_back.data(), rdispls_data_back.data(), MPI_DOUBLE_COMPLEX,
                        mpi_comm_, &request_data); 

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

void Solver::iterator_function(std::complex<double>* phi, std::complex<double>* rhs)
{

    iterator_function_acc(phi, rhs);

}

struct UserContext{
    std::function<void(std::complex<double>*, std::complex<double>*)> func;
    long long n;
};

PetscErrorCode Solver::iterator_function_2(Mat A, Vec x, Vec y)
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

std::vector<std::complex<double>> Solver::solve(const std::vector<std::complex<double>>& rhs)
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




std::vector<std::complex<double>> Solver::compute_incident_field() 
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


std::complex<double> Solver::compute_incident_field(double x, double y, double z) 
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


std::vector<std::complex<double>> Solver::solve_u_inc(bool timing) 
{
    
    MPI_Barrier(mpi_comm_);

    double start_1 = MPI_Wtime();

    std::vector<std::complex<double>> rhs = compute_incident_field();

    double end_1 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute incident field: " << end_1 - start_1 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);

    double start_2 = MPI_Wtime();

    std::vector<std::complex<double>> solution = solve(rhs);

    //std::complex<double>* rhs_data = rhs.data();
    //std::complex<double> solution2[point_up_-point_low_];
    //iterator_function(rhs_data, solution2);

    double end_2 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time solve: " << end_2 - start_2 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);

    return solution;

}

void Solver::int_far_field(const double xVers_0, const double xVers_1, const double xVers_2,
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

    double sol_re = 0.0, sol_im = 0.0;

    #pragma omp simd reduction(+:sol_re, sol_im)
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
        HH_far(xVers_0, xVers_1, xVers_2, px, py, pz, nx, ny, nz, coupling_parameter_, 
              WAVE_NUMBER, EQUATION_FORMULATION, kernel);

        sol_re += constant * (kernel.real() * phi[i].real() - kernel.imag() * phi[i].imag());
        sol_im += constant * (kernel.real() * phi[i].imag() + kernel.imag() * phi[i].real());
        
        // solution += constant * kernel * phi[i];

    }

    solution = std::complex<double>(sol_re, sol_im);

}

 std::complex<double> Solver::compute_far_field_approx(const std::vector<std::complex<double>>& phi, 
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

void Solver::int_near_field(const double x_0, const double x_1, const double x_2,
                            const long long npatch, 
                            const std::complex<double>* phi,
                            std::complex<double>& solution)
{

    int_far(x_0, x_1, x_2,
            npatch,
            phi,
            solution);

}

std::complex<double> Solver::compute_near_field_approx(const std::vector<std::complex<double>>& phi, 
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

std::complex<double> Solver::compute_far_field_exact(const double xVers_0, const double xVers_1, const double xVers_2)
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

void Solver::compute_far_field_error(const bool timing, const std::vector<std::complex<double>>& phi)
{           

    MPI_Barrier(mpi_comm_);
    double start = MPI_Wtime();

    const double deltaPhi = (M_PI - 2.0 * 1.0e-5) / (N_FAR_PTS[0] - 1);
    const double deltaTheta = 2.0 * M_PI / (N_FAR_PTS[1] - 1);

    const long long total_pts = N_FAR_PTS[0] * N_FAR_PTS[1];

    std::vector<long long> split_points_far_field(comm_size_ + 1);
    const long long points_per_rank = total_pts / comm_size_;
    long long remaining_points = total_pts % comm_size_;

    split_points_far_field[0] = 0;
    split_points_far_field[comm_size_] = total_pts;

    for (int i = 1; i < comm_size_; i++) {

        split_points_far_field[i] = split_points_far_field[i-1] + points_per_rank;

        if (remaining_points > 0) {

            split_points_far_field[i]++;
            remaining_points--;

        }

    }

    std::vector<int> recv_counts_far_field(comm_size_);
    std::vector<int> displs_far_field(comm_size_, 0);

    for (int i = 0; i < comm_size_; i++) {

        recv_counts_far_field[i] = split_points_far_field[i+1] - split_points_far_field[i];

        if (i != 0) {

            displs_far_field[i] = displs_far_field[i-1] + recv_counts_far_field[i-1];

        }

    }

    long long point_low_far_field = split_points_far_field[comm_rank_];
    long long point_up_far_field = split_points_far_field[comm_rank_+1];

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

    MPI_Allgatherv(&xVers_0_loc[0], point_up_far_field-point_low_far_field, MPI_DOUBLE, &xVers_0_all[0], &recv_counts_far_field[0], &displs_far_field[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&xVers_1_loc[0], point_up_far_field-point_low_far_field, MPI_DOUBLE, &xVers_1_all[0], &recv_counts_far_field[0], &displs_far_field[0], MPI_DOUBLE, mpi_comm_);
    MPI_Allgatherv(&xVers_2_loc[0], point_up_far_field-point_low_far_field, MPI_DOUBLE, &xVers_2_all[0], &recv_counts_far_field[0], &displs_far_field[0], MPI_DOUBLE, mpi_comm_);

    std::vector<double>().swap(xVers_0_loc);
    std::vector<double>().swap(xVers_1_loc);
    std::vector<double>().swap(xVers_2_loc);

    std::vector<std::complex<double>> far_field_approx_loc(total_pts);

    #pragma omp parallel for
    for (long long point = 0; point < total_pts; point++) {
        
        far_field_approx_loc[point] = compute_far_field_approx(phi, xVers_0_all[point], xVers_1_all[point], xVers_2_all[point]);

    }

    std::vector<std::complex<double>> far_field_approx_all(total_pts);

    MPI_Allreduce(&far_field_approx_loc[0], &far_field_approx_all[0], total_pts, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_);

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

    MPI_Allreduce(&error_1_loc, &error_1, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);
    MPI_Allreduce(&error_2_loc, &error_2, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_);

    double end = MPI_Wtime();
    MPI_Barrier(mpi_comm_);

    if (comm_rank_ == 0) {

        if (timing) {

            std::cout << "Time compute far field error: " << end - start << " seconds\n";
            print_max_RSS();

        }

        std::cout << "Far field error: " << std::setprecision(15) << std::fixed << error_1 / error_2 << "\n";

    }

}

void Solver::compute_near_field(const bool timing, const std::vector<std::complex<double>>& phi)
{

    const int nNearZones = N_NEAR_PTS.size();

    for (int i = 0; i < nNearZones; i++) {

        MPI_Barrier(mpi_comm_);
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

        std::vector<long long> split_points_near_field(comm_size_ + 1);
        const long long points_per_rank = total_pts / comm_size_;
        long long remaining_points = total_pts % comm_size_;

        split_points_near_field[0] = 0;
        split_points_near_field[comm_size_] = total_pts;

        for (int l = 1; l < comm_size_; l++) {

            split_points_near_field[l] = split_points_near_field[l-1] + points_per_rank;

            if (remaining_points > 0) {

                split_points_near_field[l]++;
                remaining_points--;

            }

        }

        std::vector<int> recv_counts_near_field(comm_size_);
        std::vector<int> displs_near_field(comm_size_, 0);

        for (int l = 0; l < comm_size_; l++) {

            recv_counts_near_field[l] = split_points_near_field[l+1] - split_points_near_field[l];

            if (l != 0) {

                displs_near_field[l] = displs_near_field[l-1] + recv_counts_near_field[l-1];

            }

        }

        long long point_low_near_field = split_points_near_field[comm_rank_];
        long long point_up_near_field = split_points_near_field[comm_rank_+1];

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
        
        MPI_Allgatherv(&x_loc[0], point_up_near_field-point_low_near_field, MPI_DOUBLE, &x_all[0], &recv_counts_near_field[0], &displs_near_field[0], MPI_DOUBLE, mpi_comm_);
        MPI_Allgatherv(&y_loc[0], point_up_near_field-point_low_near_field, MPI_DOUBLE, &y_all[0], &recv_counts_near_field[0], &displs_near_field[0], MPI_DOUBLE, mpi_comm_);
        MPI_Allgatherv(&z_loc[0], point_up_near_field-point_low_near_field, MPI_DOUBLE, &z_all[0], &recv_counts_near_field[0], &displs_near_field[0], MPI_DOUBLE, mpi_comm_);

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

        MPI_Allreduce(&u_inc_loc[0], &u_inc_all[0], total_pts, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_);
        MPI_Allreduce(&u_scat_loc[0], &u_scat_all[0], total_pts, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_);
        MPI_Allreduce(&u_total_loc[0], &u_total_all[0], total_pts, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_);

        std::vector<std::complex<double>>().swap(u_inc_loc);
        std::vector<std::complex<double>>().swap(u_scat_loc);
        std::vector<std::complex<double>>().swap(u_total_loc);

        double end = MPI_Wtime();
        MPI_Barrier(mpi_comm_);

        if (comm_rank_ == 0 && timing) {

            std::cout << "Time compute near field: " << end - start << " seconds\n";
            print_max_RSS();

        }

        std::ofstream file_i_u_inc, file_i_u_scat, file_i_u_total;

        if (comm_rank_ == 0) {

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


                                                    



