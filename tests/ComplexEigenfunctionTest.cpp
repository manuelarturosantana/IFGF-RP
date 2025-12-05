#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include "../src/solver2.h"
#include "spherical_harmonics.hpp"

/*
Note because I am lazy for this to compile one must go make everything privatve in solver2.h public
*/

int main(int argc, char* argv[]) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE) {
        throw std::logic_error("Cannot provide threaded MPI\n");
    }

    try {

        // Change the interpolation order
        Solver<3,5> new_object(1,0,0,0);
      

        // Solver new_object(1,0,0,0);
        double LAMBDA = 4.0 * 1 / 8.0; // 2.0 * M_PI / WAVE_NUMBER
        std::complex<double> WAVE_NUMBER(2.0 * M_PI / LAMBDA, -1.5);

        new_object.init_solver(true, WAVE_NUMBER, MPI_COMM_WORLD);

        solve_spherical_harmonics(new_object);
        
    }

    catch (std::exception& e) {

        std::cout << e.what() << std::endl;
        return 1;

    }

    MPI_Finalize();

    return 0;
    
}