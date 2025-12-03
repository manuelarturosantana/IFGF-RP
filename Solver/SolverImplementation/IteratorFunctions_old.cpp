#include "../solver2.h"


using namespace Eigen;

// int G_SEARCH_MAX_ITER = 50;
// double G_SEARCH_TOL = 1E-12;





void Solver::iterator_function_unacc(const VectorXcd& phi,
                                VectorXcd& rhs)
{

    std::vector<std::complex<double>> vec_coeffs(patch_num_coeffs_.size() * Nu_int_*Nv_int_);

    compute_coeffs(phi, vec_coeffs);
                
    compute_integral(phi, vec_coeffs, rhs);         

}


void Solver::iterator_function_acc(const VectorXcd& phi,
                            VectorXcd& rhs)
{

    std::vector<std::complex<double>> intensities; 

    compute_intensities(phi, intensities);

    std::vector<std::complex<double>> vec_coeffs(patch_num_coeffs_.size() * Nu_int_*Nv_int_);

    compute_coeffs(phi, vec_coeffs);

    compute_integral_acc(phi, vec_coeffs, rhs);  

    std::vector<std::complex<double>> solution_1;
    
    solution_1 = intensities;

    // boxes_.Solve<&fct_4, &fac_1>(solution_1); 
    boxes_.Solve(
    solution_1,
    [&](double a1,double a2,double a3,double a4,double a5,double a6,
        double a7,double a8,double a9,std::complex<double> a10,
        std::complex<double> c1,std::complex<double> c2,std::complex<double>& out) {
        Solver::fct_4(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,c1,c2,out);
    },
    [&](double a1,std::complex<double> c1,std::complex<double>& out) {
        Solver::fac_1(a1,c1,out);
    }
);

    #pragma omp parallel for
    for (long long i = 0; i < point_up_-point_low_; i++) {

        rhs[i] += solution_1[new_order_points_RP_[i]];
        
    }   

}



void Solver::iterator_function_overs(const VectorXcd& phi,
                                VectorXcd& rhs)
{
                
    if (METHOD_FFTW == 1) {

        std::vector<double> psi(2 * patch_num_coeffs_.size() * Nu_int_*Nv_int_);            
        std::vector<double> psi_overs;
    
        write_psi_M1(phi, psi);
        
        compute_oversampling_M1(psi, psi_overs);

        compute_integral_overs_M1(phi, psi_overs, rhs);   

    } else {

        compute_integral_overs_M2(phi, rhs);  

    }   

}

void Solver::iterator_function_overs_acc(const VectorXcd& phi,
                                    VectorXcd& rhs)
{
    
    if (METHOD_FFTW == 1) {

        std::vector<double> psi(2 * patch_num_coeffs_.size() * Nu_int_*Nv_int_);            
        std::vector<double> psi_overs;
    
        write_psi_M1(phi, psi);

        compute_oversampling_M1(psi, psi_overs);

        compute_integral_overs_acc_M1(phi, psi_overs, rhs);        

    } else {

        compute_integral_overs_acc_M2(phi, rhs);  

    }               

    std::vector<std::complex<double>> intensities; 

    compute_intensities(phi, intensities);

    std::vector<std::complex<double>> solution_1;

    solution_1 = intensities;

    // boxes_.Solve<&fct_4, &fac_1>(solution_1); 
        boxes_.Solve(
    solution_1,
    [&](double a1,double a2,double a3,double a4,double a5,double a6,
        double a7,double a8,double a9,std::complex<double> a10,
        std::complex<double> c1,std::complex<double> c2,std::complex<double>& out) {
        Solver::fct_4(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,c1,c2,out);
    },
    [&](double a1,std::complex<double> c1,std::complex<double>& out) {
        Solver::fac_1(a1,c1,out);
    }
); 

    #pragma omp parallel for
    for (long long i = 0; i < point_up_-point_low_; i++) {

        rhs[i] += solution_1[new_order_points_RP_[i]];

    }  
    
}

void Solver::iterator_function(const VectorXcd& phi,
                        VectorXcd& rhs)
{
    double start = MPI_Wtime();
    VectorXcd rhs_loc(point_up_ - point_low_);

    if (USE_ACCELERATOR) {

        if (USE_OVERSAMPLING) {

            iterator_function_overs_acc(phi, rhs_loc);

        } else {

            iterator_function_acc(phi, rhs_loc);

        }

    } else {            
    
        if (USE_OVERSAMPLING) {

            iterator_function_overs(phi, rhs_loc);

        } else {

            iterator_function_unacc(phi, rhs_loc);

        }

    }

    rhs = VectorXcd(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

    MPI_Allgatherv(&rhs_loc[0], point_up_-point_low_, MPI_DOUBLE_COMPLEX, &rhs[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE_COMPLEX, mpi_comm_);
    
    double end = MPI_Wtime();
    if (comm_rank_ == 0) {
        std::cout << "Iterator function time: " << end - start << " seconds." << std::endl;
    }
}



VectorXcd Solver::solve(const VectorXcd& rhs) 
{

    long long N = rhs.size();            
    
    auto A = [this](const VectorXcd& x, VectorXcd& solution) -> void {iterator_function(x, solution);};

    LinearOperator B(N, N, A);

    opGMRES<LinearOperator> func(B);
    func.setMaxIterations(MAX_GMRES_ITER);
    func.set_restart(MAX_GMRES_ITER);
    func.setTolerance(TOL_GMRES);
    func.set_mpi_comm(mpi_comm_);

    VectorXcd x = func.solve(rhs);

    int rank;
    MPI_Comm_rank(mpi_comm_, &rank);
    if (rank == 0) {
        std::cout << "GMRES converged in " << func.iterations() << " iterations.\n";
        std::cout << "The error was " << func.error() << std::endl;
    }

    return x;

}

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
    
VectorXcd Solver::compute_incident_field() 
{

    std::vector<std::complex<double>> u_inc(point_up_ - point_low_, 0);

    if (PLANE_OR_POINT == 0) {

        const double k_hat_0 = wavenumber_.real() * std::cos(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
        const double k_hat_1 = wavenumber_.real() * std::sin(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
        const double k_hat_2 = wavenumber_.real() * std::cos(PLANE_WAVE_PHI);
        
        #pragma omp parallel for
        for (long long i = 0; i < point_up_ - point_low_; i++) {

            const double x_0 = disc_points_x_all_[i + point_low_];
            const double x_1 = disc_points_y_all_[i + point_low_];
            const double x_2 = disc_points_z_all_[i + point_low_];               

            const double inner_product = x_0 * k_hat_0 + x_1 * k_hat_1 + x_2 * k_hat_2;

            const std::complex<double> value((-1.0) * std::cos(inner_product), (-1.0) * std::sin(inner_product));

            u_inc[i] = value;

        }

    } else {

        const double k = wavenumber_.real();
        
        #pragma omp parallel for
        for (long long j = 0; j < point_up_ - point_low_; j++) {

            const double y_0 = disc_points_x_all_[j + point_low_];
            const double y_1 = disc_points_y_all_[j + point_low_];
            const double y_2 = disc_points_z_all_[j + point_low_];  

            for (int i = 0; i < NUM_POINT_SOURCES; i++) {

                const double x_0 = POINT_SOURCE_CENTER[i][0];
                const double x_1 = POINT_SOURCE_CENTER[i][1];
                const double x_2 = POINT_SOURCE_CENTER[i][2];                                 

                const double diff = std::sqrt((x_0-y_0)*(x_0-y_0) + (x_1-y_1)*(x_1-y_1) + (x_2-y_2)*(x_2-y_2));

                const std::complex<double> value((-1.0) * std::cos(k *diff) / diff, (-1.0) * std::sin(k * diff) / diff);

                u_inc[j] += value;

            }

        }

    }

    VectorXcd u_inc_all(Q_ * Qx_*Qy_ * Nu_int_*Nv_int_);

    MPI_Allgatherv(&u_inc[0], point_up_-point_low_, MPI_DOUBLE_COMPLEX, &u_inc_all[0], &recv_counts_2_[0], &displs_2_[0], MPI_DOUBLE_COMPLEX, mpi_comm_);

    return u_inc_all;

}     
        
    

VectorXcd Solver::solve_u_inc(bool timing) 
{
    
    MPI_Barrier(mpi_comm_);


    double start_1 = MPI_Wtime();

    VectorXcd rhs = compute_incident_field();

    double end_1 = MPI_Wtime();

    if (timing && comm_rank_ == 0) {

        std::cout << "Time compute incident field: " << end_1 - start_1 << " seconds\n";
        print_max_RSS();

    }

    MPI_Barrier(mpi_comm_);

    
    double start_2 = MPI_Wtime();

    VectorXcd solution = solve(rhs);

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

    for (int i = 0; i < Nu_int_; i++) {
        for (int j = 0; j < Nv_int_; j++) {

            const double dsdtjac_loc = dsdtjac_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double constant = dsdtjac_loc * fejer_weights_u_int_[i] * fejer_weights_v_int_[j];

            const double px = disc_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double py = disc_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double pz = disc_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

            const double nx = norm_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double ny = norm_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double nz = norm_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

            std::complex<double> kernel;
            HH_far(xVers_0, xVers_1, xVers_2, px, py, pz, nx, ny, nz, coupling_parameter_, wavenumber_.real(),kernel);

            solution += constant * kernel * phi[i*Nv_int_ + j];

        }
    }

}

std::complex<double> Solver::compute_far_field_approx(const VectorXcd& phi, 
                                                const double xVers_0, const double xVers_1, const double xVers_2) 
{
    
    std::complex<double> solution_rank(0.0, 0.0);

    for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

        std::complex<double> solution_loc;

        int_far_field(xVers_0, xVers_1, xVers_2,
                        patch_num,
                        &phi[patch_num * Nu_int_*Nv_int_],
                        solution_loc);                        
        
        solution_rank += solution_loc;

    }

    solution_rank *= 1.0 / (4.0 * M_PI);

    return solution_rank;

}

std::complex<double> Solver::compute_far_field_exact(const double xVers_0, const double xVers_1, const double xVers_2) 
{

    std::complex<double> solution(0.0, 0.0);

    const std::complex<double> I(0.0, 1.0);

    int nterms = 0;

    std::complex<double> yn = std::sph_neumann(nterms, wavenumber_.real() * SPHERE_RADIUS);
    std::complex<double> jn = std::sph_bessel(nterms, wavenumber_.real() * SPHERE_RADIUS);
    std::complex<double> hn = jn + I * yn;

    while (std::abs(hn) * std::abs(hn) < 1.0e14) {

        nterms++;

        yn = std::sph_neumann(nterms, wavenumber_.real() * SPHERE_RADIUS);
        jn = std::sph_bessel(nterms, wavenumber_.real() * SPHERE_RADIUS);
        hn = jn + I * yn;

    }

    const double kVers[3] = {std::cos(0.0) * std::sin(M_PI), std::sin(0.0) * std::sin(M_PI), std::cos(M_PI)};
    const double x = xVers_0 * kVers[0] + xVers_1 * kVers[1] + xVers_2 * kVers[2];

    for (int n = 0; n <= nterms; n++) {

        const std::complex<double> yn = std::sph_neumann(n, wavenumber_.real() * SPHERE_RADIUS);
        const std::complex<double> jn = std::sph_bessel(n, wavenumber_.real() * SPHERE_RADIUS);
        const std::complex<double> hn = jn + I * yn;

        const double P = std::legendre(n, x);

        solution += (2.0 * n + 1.0) * jn/hn * P;

    }

    solution *= I / (wavenumber_.real() * SPHERE_RADIUS);

    return solution;

}

void Solver::compute_far_field_error(const bool timing, const VectorXcd& phi)
{

    MPI_Barrier(mpi_comm_);
    double start = MPI_Wtime();

    const double deltaPhi = M_PI / (200 - 1);
    const double deltaTheta = 2.0 * M_PI / (200 - 1);

    const long long total_pts = 200 * 200;

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

        const long long m = (point + point_low_far_field) / 200;
        const long long n = (point + point_low_far_field) % 200;

        const double phi_m = (m - 1) * deltaPhi;
        const double theta_n = (n - 1) * deltaTheta;

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
        
        const double far_field_approx_norm = std::sqrt(far_field_approx.real()*far_field_approx.real() + far_field_approx.imag()*far_field_approx.imag());
        const double far_field_exact_norm = std::sqrt(far_field_exact.real()*far_field_exact.real() + far_field_exact.imag()*far_field_exact.imag());
            
        const double value1 = std::abs(far_field_exact_norm - far_field_approx_norm);
        const double value2 = std::abs(far_field_exact_norm);

        if (error_1_loc < value1) error_1_loc = value1;
        if (error_2_loc < value2) error_2_loc = value2;

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

        std::cout << "Far field error: " << error_1 / error_2 << "\n";

    }

}

    // /*std::complex<double> compute_incident_field(double x, double y, double z) 
    // {

    //     std::complex<double> u_inc = {0.0, 0.0};

    //     if (PLANE_OR_POINT == 0) {

    //         const double k_hat_0 = WAVE_NUMBER * std::cos(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
    //         const double k_hat_1 = WAVE_NUMBER * std::sin(PLANE_WAVE_THE) * std::sin(PLANE_WAVE_PHI);
    //         const double k_hat_2 = WAVE_NUMBER * std::cos(PLANE_WAVE_PHI);            

    //         const double inner_product = x * k_hat_0 + y * k_hat_1 + z * k_hat_2;

    //         const std::complex<double> value((-1.0) * std::cos(inner_product), (-1.0) * std::sin(inner_product));

    //         u_inc = value;

    //     } else {

    //         const double k = WAVE_NUMBER;

    //         for (int i = 0; i < NUM_POINT_SOURCES; i++) {

    //             const double x_0 = POINT_SOURCE_CENTER[i][0];
    //             const double x_1 = POINT_SOURCE_CENTER[i][1];
    //             const double x_2 = POINT_SOURCE_CENTER[i][2];                                 

    //             const double diff = std::sqrt((x_0-x)*(x_0-x) + (x_1-y)*(x_1-y) + (x_2-z)*(x_2-z));

    //             const std::complex<double> value((-1.0) * std::cos(k * diff) / diff, (-1.0) * std::sin(k * diff) / diff);

    //             u_inc += value;

    //         }

    //     }

    //     return u_inc;

    // }  */

void Solver::int_near_field(const double x_0, const double x_1, const double x_2,
                    const long long npatch, 
                    const std::complex<double>* phi,
                    std::complex<double>& solution)                      
{

    solution = std::complex<double>(0.0, 0.0);

    for (int i = 0; i < Nu_int_; i++) {
        for (int j = 0; j < Nv_int_; j++) {

            const double dsdtjac_loc = dsdtjac_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double constant = dsdtjac_loc * fejer_weights_u_int_[i] * fejer_weights_v_int_[j];

            const double px = disc_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double py = disc_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double pz = disc_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

            const double nx = norm_points_x_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double ny = norm_points_y_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];
            const double nz = norm_points_z_all_[npatch * Nu_int_*Nv_int_ + i * Nv_int_ + j];

            std::complex<double> kernel;
            HH2(x_0, x_1, x_2, px, py, pz, nx, ny, nz, coupling_parameter_, wavenumber_, kernel);

            solution += constant * kernel * phi[i*Nv_int_ + j];

        }
    }

}

std::complex<double> Solver::compute_near_field_approx(const VectorXcd& phi, 
                                                const double x_0, const double x_1, const double x_2) 
{
    
    std::complex<double> solution_rank(0.0, 0.0);

    for (long long patch_num = patch_low_; patch_num < patch_up_; patch_num++) {

        std::complex<double> solution_loc;

        int_near_field(x_0, x_1, x_2,
                        patch_num,
                        &phi[patch_num * Nu_int_*Nv_int_],
                        solution_loc);                        
        
        solution_rank += solution_loc;

    }

    return solution_rank;
}
/*
Custom function to evaluate the near field, based off of a small change from compute
near_field_approx
*/
std::complex<double> Solver::compute_near_field(const VectorXcd& phi, 
                                                const double x_0, const double x_1, const double x_2) 
{
    
    std::complex<double> solution_rank(0.0, 0.0);
    if (comm_rank_ == 0) {
        #pragma omp parallel for
        for (long long patch_num = 0; patch_num < Q_*Qx_*Qy_; patch_num++) {

            std::complex<double> solution_loc;

            int_near_field(x_0, x_1, x_2,
                        patch_num,
                        &phi[patch_num * Nu_int_*Nv_int_],
                        solution_loc);                        
            
            solution_rank += solution_loc;
        }
    }
    // Make sure every process has the solution :)
    MPI_Bcast(&solution_rank, 1, MPI_DOUBLE_COMPLEX, 0, mpi_comm_);

    return solution_rank;

}


std::complex<double> Solver::iterator_function_Manuel(const VectorXcd& u, const VectorXcd& v) 
{

    VectorXcd y = solve(v);

    return (u.conjugate().transpose() * y);

}



