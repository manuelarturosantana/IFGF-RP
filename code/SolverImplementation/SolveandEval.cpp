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