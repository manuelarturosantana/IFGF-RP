#include <iostream>
#include <cstdlib>
#include <random>
#include <iomanip>
#include "solver2.h"
#include <chrono>

int main(int argc, char* argv[]) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE) {
        throw std::logic_error("Cannot provide threaded MPI\n");
    }

    try {

        int N_PATCHES_ORIG = 6;
        int world_rank;

        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        std::default_random_engine generator(1);
        // std::normal_distribution<double> distribution(0.0, 1.0);
         std::uniform_real_distribution<double> distribution(-1,1);

        int n_iters = 1;

        std::vector<std::complex<double>> vec_k(n_iters);
        std::vector<std::complex<double>> vec_val(n_iters);

        for (int i = 0; i < n_iters; i++) {

            // Wavenumber which approximates an eigenvalue
            std::complex<double> k (2, 0);

            // Number of points used per dimension per patch in every patch
            int n_pts_per_patch[2] = {10, 10};
            // Number of patch partitions along first and second dimensions
            int n_split_per_patch[2] = {(1 << (i + 3)), (1 << (i + 3))};
            
            // Number of points used with the singular integration
            int n_pts_sing_int[2] = {40, 40};
            // Proximity distance that determines near singular integrals
            double proximity_box_size = 0.0125 / (1 << i);

            // IFGF number of levels
            // int n_levels_IFGF = 4 + i;
            int n_levels_IFGF = 3;

            // Total points
            long long N = N_PATCHES_ORIG * n_split_per_patch[0]*n_split_per_patch[1] * n_pts_per_patch[0]*n_pts_per_patch[1];

            // Random vectors
            Eigen::VectorXcd u(N), v(N);

            // Set unique vectors for all ranks
            if (world_rank == 0) {

                // auto start = std::chrono::high_resolution_clock::now();
                #pragma omp parallel for
                for (int i = 0; i < N; i++) {

                    u[i] = std::complex<double>(distribution(generator), distribution(generator));
                    v[i] = std::complex<double>(distribution(generator), distribution(generator));

                }

                

                // auto end = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> elapsed = end - start;
                // std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
            }

            MPI_Bcast(&u[0], N, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&v[0], N, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

            // Create object
            auto precomp_start = std::chrono::high_resolution_clock::now();
            Solver new_object;
            new_object.set_eq_form(3);
            new_object.init_solver(true, k, n_pts_per_patch, n_split_per_patch, n_pts_sing_int, 
                proximity_box_size, n_levels_IFGF, MPI_COMM_WORLD);

            auto precomp_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> precomp_elapsed = precomp_end - precomp_start;

            PetscErrorCode ierr;  

            ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
        
            if (world_rank == 0) {
                std::cout << "Coupling param" << " " << new_object.get_coup_param() << std::endl;
            }

            // Full solve test
            VectorXcd vv = new_object.compute_incident_field();
            VectorXcd ww = new_object.solve(vv);
            // VectorXcd ww = new_object.solve(u);
            new_object.compute_far_field_error(false, ww);
            
            // Eigenfunction test
            // new_object.solve_spherical_harmonics();

            ierr = PetscFinalize(); 
            CHKERRABORT(PETSC_COMM_WORLD, ierr);

        }

    }

    catch (std::exception& e) {

        std::cout << e.what() << std::endl;
        return 1;

    }

    MPI_Finalize();

    return 0;
    
}