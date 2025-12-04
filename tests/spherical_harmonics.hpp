#include "../src/solver2.h"
#include "complex_bessel-master/include/complex_bessel.h"
/*
Note because I am lazy and don't want to write setters and getters for everything,
to make this work one would have to go make all the private members public in solver.
*/
template<int PS, int PT>
std::vector<std::complex<double>> compute_spherical_harmonics(Solver<PS, PT>& S)
{

    std::vector<std::complex<double>> u_inc(S.point_up_ - S.point_low_, 0);

    #pragma omp parallel for
    for (long long i = 0; i < S.point_up_ - S.point_low_; i++) {

        long long idx;
        if (S.USE_ACCELERATOR) {
            idx = i;
        } else {
            idx = i + S.point_low_;
        }

        const double x = S.disc_points_x_all_[idx];
        const double y = S.disc_points_y_all_[idx];
        const double z = S.disc_points_z_all_[idx];      

        std::complex<double> I(0.0, 1.0);

        // Spherical harmonic (5, 2)
        u_inc[i] = (1.0/8.0) * std::sqrt(1155.0/(2.0*M_PI)) * (x + I * y) * (x + I * y) * (3.0 * z*z*z - z);

    }

    return u_inc;

}

template<int PS, int PT>
std::vector<std::complex<double>> compute_spherical_harmonics_solution_exact(Solver<PS,PT>& S)
{

    std::vector<std::complex<double>> phi(S.point_up_ - S.point_low_);
    std::vector<std::complex<double>> real_rhs(S.point_up_ - S.point_low_);

    std::complex<double> k = S.WAVE_NUMBER;
    // double k = wavenumber_.real();
    std::complex<double> I(0.0, 1.0);
    //  std::complex<double> k = wavenumber_;

    double nu = 5.5;
    auto factor = std::sqrt(M_PI / (2.0 * k));

    std::complex<double> j_5  = factor * sp_bessel::besselJ(nu, k);
    std::complex<double> dj_5 = factor * sp_bessel::besselJp(nu, k, 1)
                            - 0.5 * sp_bessel::besselJ(nu, k) * factor / k;

    // std::complex<double> h_5  = factor * sp_bessel::hankelH1(nu, k);
    std::complex<double> h_5 = sp_bessel::sph_hankelH1(5,k);
    std::complex<double> dh_5 = factor * sp_bessel::hankelH1p(nu, k, 1)
                            - 0.5 * sp_bessel::hankelH1(nu, k) * factor / k;
    
      // FOR SOME REASON THE HANKEL FUNCTION IS I TIMES THIS. FOLLOWING THE CONVENTION OF THE RP PAPER
    h_5 *= I;
    dh_5 *= I;
                           
 
    #pragma omp parallel for
    for (long long i = 0; i < S.point_up_ - S.point_low_; i++) {

        long long idx;
        if (S.USE_ACCELERATOR) {
            idx = i;
        } else {
            idx = i + S.point_low_;
        }

        const double x = S.disc_points_x_all_[idx];
        const double y = S.disc_points_y_all_[idx];
        const double z = S.disc_points_z_all_[idx];      

        phi[i] = (1.0/8.0) * std::sqrt(1155.0/(2.0*M_PI)) * (x + I * y) * (x + I * y) * (3.0 * z*z*z - z);
        
    }

    if (S.EQUATION_FORMULATION == 1) {

        #pragma omp parallel for
        for (long long i = 0; i < S.point_up_ - S.point_low_; i++) {
        
            real_rhs[i] = k * j_5 * h_5 * phi[i];

        }

    } else if (S.EQUATION_FORMULATION == 2) {

        #pragma omp parallel for
        for (long long i = 0; i < S.point_up_ - S.point_low_; i++) {
        
            real_rhs[i] = 0.5 * phi[i] + (k*k) * 0.5 * (j_5 * dh_5 + h_5 * dj_5) * phi[i];

        }

    } else {

        #pragma omp parallel for
        for (long long i = 0; i < S.point_up_ - S.point_low_; i++) {
        
            real_rhs[i] = 0.5 * phi[i] + ((k*k) * 0.5 * (j_5 * dh_5 + h_5 * dj_5) * phi[i]) - I * S.get_coup_param() * (k * j_5 * h_5 * phi[i]);

        }

    }  

    return real_rhs;

}

template<int PS, int PT>
void  solve_spherical_harmonics(Solver<PS,PT>& S)
{

    std::vector<std::complex<double>> rhs = compute_spherical_harmonics(S);

    
    std::vector<std::complex<double>> approx_sol(S.get_num_unknowns());

    S.iterator_function(rhs.data(), approx_sol.data());

    std::vector<std::complex<double>> exact_sol = compute_spherical_harmonics_solution_exact(S);
    // Compute MSE

    double error = 0.0;

    for (long long i = 0; i < S.point_up_ - S.point_low_; i++) {

        error += std::abs(exact_sol[i] - approx_sol[i]) * std::abs(exact_sol[i] - approx_sol[i]);

    }

    double error_total = 0.0;

    MPI_Allreduce(&error, &error_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (S.comm_rank_ == 0) {
    
        std::cout << "MSE: " << error_total / (S.get_num_unknowns()) << "\n";

    }

}

