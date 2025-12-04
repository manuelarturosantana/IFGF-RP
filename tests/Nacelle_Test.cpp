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

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

        // Change the interpolation order
        Solver<3,5> new_object("../Nacelle/", "Nacelle-");

        std::complex<double> WAVE_NUMBER(16, 0.0);
        new_object.set_n_pts_per_patch(12,12);

        new_object.init_solver(true, WAVE_NUMBER, MPI_COMM_WORLD);
        
        if (rank == 0) {
            std::cout << "The number of unknowns is " << new_object.get_num_unknowns() << std::endl;
            std::cout << "The number of patch splits is (" << new_object.get_patch_split_x() << ", " << new_object.get_patch_split_y() << ") " << std::endl;
        }
        

    }

    catch (std::exception& e) {

        std::cout << e.what() << std::endl;
        return 1;

    }

    MPI_Finalize();

    return 0;
    
}