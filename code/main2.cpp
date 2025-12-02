#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include "solver2.h"

int main(int argc, char* argv[]) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE) {
        throw std::logic_error("Cannot provide threaded MPI\n");
    }

    try {

        Solver new_object(true);

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