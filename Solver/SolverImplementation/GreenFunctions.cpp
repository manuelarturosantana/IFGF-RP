// #include "../../complex_bessel-master/include/complex_bessel.h"
#include "../solver2.h"


using namespace Eigen;

// int G_SEARCH_MAX_ITER = 50;
// double G_SEARCH_TOL = 1E-12;

void Solver::HH(const double p1x, const double p1y, const double p1z,
        const double p2x, const double p2y, const double p2z,
        const double nx, const double ny, const double nz,
        const double dxdsx, const double dxdsy, const double dxdsz,
        const double dxdtx, const double dxdty, const double dxdtz,
        const double dxdsdsx, const double dxdsdsy, const double dxdsdsz,
        const double dxdsdtx, const double dxdsdty, const double dxdsdtz,
        const double dxdtdtx, const double dxdtdty, const double dxdtdtz,
        const std::complex<double> coupling_parameter, std::complex<double> wavenumber,
        std::complex<double>& solution)
{

    double solution_real, solution_imag;

    const double norm_diff = std::sqrt((p2x - p1x)*(p2x - p1x) + (p2y - p1y)*(p2y - p1y) + (p2z - p1z)*(p2z - p1z));

    if (norm_diff < 1e-14) {

        solution_real = 0.0;
        solution_imag = 0.0;

    } else {   

        const double val_cos = std::cos(wavenumber.real() * norm_diff);
        const double val_sin = std::sin(wavenumber.real() * norm_diff);
        const double val_const = std::exp(-wavenumber.imag() * norm_diff);

        if (EQUATION_FORMULATION == 1) {

            solution_real = val_const * val_cos / (4.0 * M_PI * norm_diff);
            solution_imag = val_const * val_sin / (4.0 * M_PI * norm_diff);

        } else if (EQUATION_FORMULATION == 2) {

            // Compute beta = < r, n > / |r|^2
            
            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            if ((norm_diff < 1.0E-5) && (std::abs(rDotNorm) < 1.0E-6)) {

                beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, nx, ny, nz, dxdsx, dxdsy, dxdsz, dxdtx, dxdty, dxdtz, 
                                                dxdsdsx, dxdsdsy, dxdsdsz, dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

            }            

            solution_real = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_cos / norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
            solution_imag = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_sin / norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

        } else {

            // EQUATION_FORMULATION == 3 

            const double solution_real_S = val_const * val_cos / (4.0 * M_PI * norm_diff);
            const double solution_imag_S = val_const * val_sin / (4.0 * M_PI * norm_diff);

            // Compute beta = < r, n > / |r|^2

            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            if ((norm_diff < 1.0E-5) && (std::abs(rDotNorm) < 1.0E-6)) {

                beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, nx, ny, nz, dxdsx, dxdsy, dxdsz, dxdtx, dxdty, dxdtz, 
                                                dxdsdsx, dxdsdsy, dxdsdsz, dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

            } 

            const double solution_real_D = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_cos / norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
            const double solution_imag_D = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_sin / norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

            std::complex<double> solution_D = std::complex<double>(solution_real_D, solution_imag_D);
            std::complex<double> solution_S = std::complex<double>(solution_real_S, solution_imag_S);

            std::complex<double> I(0.0,1.0);
            std::complex<double> sol = solution_D - I * coupling_parameter * solution_S;
            solution_real = sol.real();
            solution_imag = sol.imag();

            // solution_real = solution_real_D + coupling_parameter * solution_imag_S;
            // solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

        }

        solution = std::complex<double>{solution_real, solution_imag};

    }

}

void Solver::HH2(const double p1x, const double p1y, const double p1z,
         const double p2x, const double p2y, const double p2z,
         const double nx, const double ny, const double nz,
         std::complex<double> coupling_parameter, std::complex<double> wavenumber,
         std::complex<double>& solution)
{

    double solution_real, solution_imag;

    const double norm_diff = std::sqrt((p2x - p1x)*(p2x - p1x) + (p2y - p1y)*(p2y - p1y) + (p2z - p1z)*(p2z - p1z));

    if (norm_diff < 1e-14) {

        solution_real = 0.0;
        solution_imag = 0.0;

    } else {   

        const double val_cos = std::cos(wavenumber.real() * norm_diff);
        const double val_sin = std::sin(wavenumber.real() * norm_diff);
        const double val_const = std::exp(-wavenumber.imag() * norm_diff);

        if (EQUATION_FORMULATION == 1) {

            solution_real = val_const * val_cos / (4.0 * M_PI * norm_diff);
            solution_imag = val_const * val_sin / (4.0 * M_PI * norm_diff);

        } else if (EQUATION_FORMULATION == 2) {

            // Compute beta = < r, n > / |r|^2
            
            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            solution_real = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_cos / norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
            solution_imag = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_sin / norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

            if (USE_ACCELERATOR) {

                solution_real *= -1.0;
                solution_imag *= -1.0;

            }

        } else {

            // EQUATION_FORMULATION == 3 

            const double solution_real_S = val_const * val_cos / (4.0 * M_PI * norm_diff);
            const double solution_imag_S = val_const * val_sin / (4.0 * M_PI * norm_diff);

            // Compute beta = < r, n > / |r|^2

            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            double solution_real_D = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_cos / norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
            double solution_imag_D = -(1.0 / (4.0 * M_PI)) * INT_EXT * beta * val_const * (val_sin / norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

            if (USE_ACCELERATOR) {

                solution_real_D *= -1.0;
                solution_imag_D *= -1.0;

            }
            std::complex<double> solution_D = std::complex<double>(solution_real_D, solution_imag_D);
            std::complex<double> solution_S = std::complex<double>(solution_real_S, solution_imag_S);

            std::complex<double> I(0.0,1.0);
            std::complex<double> sol = solution_D - I * coupling_parameter * solution_S;
            solution_real = sol.real();
            solution_imag = sol.imag();

            // solution_real = solution_real_D + coupling_parameter * solution_imag_S;
            // solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

        }

        solution = std::complex<double>{solution_real, solution_imag};

    }

}

void Solver::HH_far(const double xVers_0, const double xVers_1, const double xVers_2,
            const double y_0, const double y_1, const double y_2,
            const double n_0, const double n_1, const double n_2,
            const std::complex<double> coupling_parameter, const double wavenumber,
            std::complex<double>& solution)
{

    double solution_real, solution_imag;

    const double inner_prod = - wavenumber * (xVers_0 * y_0 + xVers_1 * y_1 + xVers_2 * y_2);
    const double inner_prod_2 = - wavenumber * (xVers_0 * n_0 + xVers_1 * n_1 + xVers_2 * n_2);

    const double S_real = std::cos(inner_prod);
    const double S_imag = std::sin(inner_prod);
    const double D_real = -inner_prod_2 * S_imag;
    const double D_imag = inner_prod_2 * S_real;

    if (EQUATION_FORMULATION == 1) {

        solution_real = S_real;
        solution_imag = S_imag;

    } else if (EQUATION_FORMULATION == 2) {

        solution_real = D_real;
        solution_imag = D_imag;

    } else {

        // EQUATION_FORMULATION == 3

        std::complex<double> solution_D = std::complex<double>(D_real, D_imag);
        std::complex<double> solution_S = std::complex<double>(S_real, S_imag);

        std::complex<double> I(0.0,1.0);
        std::complex<double> sol = solution_D - I * coupling_parameter * solution_S;
        solution_real = sol.real();
        solution_imag = sol.imag();

        // solution_real = D_real + coupling_parameter * S_imag;
        // solution_imag = D_imag - coupling_parameter * S_real;

    } 

    solution = std::complex<double>{solution_real, solution_imag};

}

// Spherical Harmoics test. Doesn't really belong in Greens function, but it was the best I 
// could come up with without creating a new file
// WARNING: Due to some compatability issues this doesn't work
// for more than one processor on MPI
// std::vector<std::complex<double>> Solver::compute_spherical_harmonics()
// {

//     std::vector<std::complex<double>> u_inc(point_up_ - point_low_, 0);

//     #pragma omp parallel for
//     for (long long i = 0; i < point_up_ - point_low_; i++) {

//         long long idx;
//         if (USE_ACCELERATOR) {
//             idx = i;
//         } else {
//             idx = i + point_low_;
//         }

//         const double x = disc_points_x_all_[idx];
//         const double y = disc_points_y_all_[idx];
//         const double z = disc_points_z_all_[idx];      

//         std::complex<double> I(0.0, 1.0);

//         // Spherical harmonic (5, 2)
//         u_inc[i] = (1.0/8.0) * std::sqrt(1155.0/(2.0*M_PI)) * (x + I * y) * (x + I * y) * (3.0 * z*z*z - z);

//     }

//     return u_inc;

// }

// std::vector<std::complex<double>> Solver::compute_spherical_harmonics_solution_exact()
// {

//     std::vector<std::complex<double>> phi(point_up_ - point_low_);
//     std::vector<std::complex<double>> real_rhs(point_up_ - point_low_);

//     std::complex<double> k = wavenumber_;
//     // double k = wavenumber_.real();
//     std::complex<double> I(0.0, 1.0);
//     //  std::complex<double> k = wavenumber_;

//     double nu = 5.5;
//     auto factor = std::sqrt(M_PI / (2.0 * k));

//     std::complex<double> j_5  = factor * sp_bessel::besselJ(nu, k);
//     std::complex<double> dj_5 = factor * sp_bessel::besselJp(nu, k, 1)
//                             - 0.5 * sp_bessel::besselJ(nu, k) * factor / k;

//     // std::complex<double> h_5  = factor * sp_bessel::hankelH1(nu, k);
//     std::complex<double> h_5 = sp_bessel::sph_hankelH1(5,k);
//     std::complex<double> dh_5 = factor * sp_bessel::hankelH1p(nu, k, 1)
//                             - 0.5 * sp_bessel::hankelH1(nu, k) * factor / k;
    
//       // FOR SOME REASON THE HANKEL FUNCTION IS I TIMES THIS. FOLLOWING THE CONVENTION OF THE RP PAPER
//     h_5 *= I;
//     dh_5 *= I;
                           
//     // std::complex<double> j_5  = std::sqrt((M_PI / (2.0 * k)))  * sp_bessel::besselJ(5.5, k);
//     // std::complex<double> dj_5 = std::sqrt((M_PI / (2.0 * k))) * sp_bessel::besselJp(5.5, k,1) + sp_bessel::besselJ(5.5, k) * std::sqrt(M_PI / 2.0) * (-0.5) * std::pow(k, -3.0/2.0);
    
//     // std::complex<double> h_5  = std::sqrt((M_PI / (2.0 * k))) * sp_bessel::hankelH1(5.5, k);
//     // std::complex<double> dh_5 = std::sqrt((M_PI / (2.0 * k))) * sp_bessel::hankelH1p(5.5, k,1) + sp_bessel::hankelH1(5.5, k) * std::sqrt(M_PI / 2.0) * (-0.5) * std::pow(k, -3.0/2.0);
    
//     // FOR SOME REASON THE HANKEL FUNCTION IS I TIMES THIS. FOLLOWING THE CONVENTION OF THE RP PAPER
//     // double j_52 = (15 * (std::pow(k, 4) - 28 * k*k + 63) * sin(k)) / std::pow(k, 6) + ((-std::pow(k, 4) + 105 * k*k - 945) * cos(k)) / std::pow(k, 5);
//     // double dj_52 = ((16 * std::pow(k, 4) - 735 * k*k + 5670) * cos(k)) / std::pow(k, 6) + ((std::pow(k, 6) - 135 * std::pow(k, 4) + 2625 * k*k - 5670) * sin(k)) / std::pow(k, 7);
//     // std::complex<double> h_52((k*(std::pow(k,4)-105*std::pow(k,2)+945)*sin(k)+15*(std::pow(k,4)-28*k*k+63)*cos(k))/std::pow(k,6), (15*(std::pow(k,4)-28*k*k+63)*sin(k)-k*(std::pow(k,4)-105*k*k+945)*cos(k))/std::pow(k,6));
//     // std::complex<double> dh_52((k*(-16.0*std::pow(k,4.0)+735*k*k-5670)*sin(k)+(std::pow(k,6)-135.0*std::pow(k,4)+2625*k*k-5670)*cos(k))/std::pow(k,7), (k*(16.0*std::pow(k,4)-735*k*k+5670)*cos(k)+(std::pow(k,6)-135.0*std::pow(k,4)+22625.0*k*k-5670.0)*sin(k))/std::pow(k,7));
//     // // check the difference and print
//     // std::cout << "Difference in j_5: " << std::abs(j_5 - j_52) << std::endl;
//     // std::cout << "Difference in dj_5: " << std::abs(dj_5 - dj_52) << std::endl;
//     // std::cout << "Difference in h_5: " << std::abs(h_5 - h_52) << std:: endl;
//     // std::cout << "Difference in dh_5: " << std::abs(dh_5 - dh_52) << std:: endl;
//     // std::cout << "h_5: " << h_5 << std:: endl;
//     // std::cout << "dh_5: " << dh_5 << std:: endl;
//     // std::cout << "h_52: " << h_52 << std:: endl;
//     // std::cout << "dh_52: " << dh_52 << std:: endl;


//     #pragma omp parallel for
//     for (long long i = 0; i < point_up_ - point_low_; i++) {

//         long long idx;
//         if (USE_ACCELERATOR) {
//             idx = i;
//         } else {
//             idx = i + point_low_;
//         }

//         const double x = disc_points_x_all_[idx];
//         const double y = disc_points_y_all_[idx];
//         const double z = disc_points_z_all_[idx];      

//         phi[i] = (1.0/8.0) * std::sqrt(1155.0/(2.0*M_PI)) * (x + I * y) * (x + I * y) * (3.0 * z*z*z - z);
        
//     }

//     if (EQUATION_FORMULATION == 1) {

//         #pragma omp parallel for
//         for (long long i = 0; i < point_up_ - point_low_; i++) {
        
//             real_rhs[i] = k * j_5 * h_5 * phi[i];

//         }

//     } else if (EQUATION_FORMULATION == 2) {

//         #pragma omp parallel for
//         for (long long i = 0; i < point_up_ - point_low_; i++) {
        
//             real_rhs[i] = 0.5 * phi[i] + (k*k) * 0.5 * (j_5 * dh_5 + h_5 * dj_5) * phi[i];

//         }

//     } else {

//         #pragma omp parallel for
//         for (long long i = 0; i < point_up_ - point_low_; i++) {
        
//             real_rhs[i] = 0.5 * phi[i] + ((k*k) * 0.5 * (j_5 * dh_5 + h_5 * dj_5) * phi[i]) - I * coupling_parameter_ * (k * j_5 * h_5 * phi[i]);

//         }

//     }  

//     return real_rhs;

// }

// void  Solver::solve_spherical_harmonics()
// {

//     std::vector<std::complex<double>> rhs = compute_spherical_harmonics();

//     // Eigen::VectorXcd approx_sol(point_up_ - point_low_);
//     Eigen::VectorXcd approx_sol(this->get_num_unknowns());

//     Eigen::Map<Eigen::VectorXcd> eigen_rhs(rhs.data(), rhs.size());

//     std::cout << "Made it here " << std::endl;
//     iterator_function(eigen_rhs, approx_sol);

//     std::vector<std::complex<double>> exact_sol = compute_spherical_harmonics_solution_exact();
//     std::cout << "Made it here2 " << std::endl;
//     // Compute MSE

//     double error = 0.0;

//     for (long long i = 0; i < point_up_ - point_low_; i++) {

//         error += std::abs(exact_sol[i] - approx_sol(i)) * std::abs(exact_sol[i] - approx_sol(i));

//     }

//     double error_total = 0.0;

//     MPI_Allreduce(&error, &error_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//     if (comm_rank_ == 0) {
    
//         std::cout << "MSE: " << error_total / (Q_ * Qx_ * Qy_ * Nu_int_ * Nv_int_) << "\n";

//     }

// }
