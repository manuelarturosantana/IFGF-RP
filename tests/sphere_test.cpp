#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include "../src/solver2.h"

int main(int argc, char* argv[]) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE) {
        throw std::logic_error("Cannot provide threaded MPI\n");
    }

    try {

        // Change the interpolation order
        Solver<3,5> new_object(1,0,0,0);

        new_object.set_coup_param(std::complex<double>(1.0, 0.0));

        // Solver new_object(1,0,0,0);
        double LAMBDA = 2.0 * 1 / 8.0; // 2.0 * M_PI / WAVE_NUMBER
        std::complex<double> WAVE_NUMBER(2.0 * M_PI / LAMBDA,0);

        new_object.init_solver(true, WAVE_NUMBER, MPI_COMM_WORLD);

        PetscErrorCode ierr;  

        ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);

        std::vector<std::complex<double>> phi = new_object.solve_u_inc(true);

        ierr = PetscFinalize(); 
        CHKERRABORT(PETSC_COMM_WORLD, ierr);

        new_object.compute_far_field_error(true, phi);
        
    }

    catch (std::exception& e) {

        std::cout << e.what() << std::endl;
        return 1;

    }

    MPI_Finalize();

    return 0;
    
}