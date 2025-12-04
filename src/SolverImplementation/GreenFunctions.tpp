#pragma once

#include "../solver2.h"

// template<int PS, int PT>
// void Solver<PS,PT>::HH(
//     const double p1x, const double p1y, const double p1z,
//     const double p2x, const double p2y, const double p2z,
//     const double nx, const double ny, const double nz,
//     const double dxdsx, const double dxdsy, const double dxdsz,
//     const double dxdtx, const double dxdty, const double dxdtz,
//     const double dxdsdsx, const double dxdsdsy, const double dxdsdsz,
//     const double dxdsdtx, const double dxdsdty, const double dxdsdtz,
//     const double dxdtdtx, const double dxdtdty, const double dxdtdtz,
//     const double coupling_parameter, double wavenumber,
//     std::complex<double>& solution)
// {
    
//     const double INV_4_PI = 1.0 / (4.0 * M_PI);

//     const double dx = p2x - p1x;
//     const double dy = p2y - p1y;
//     const double dz = p2z - p1z;
    
//     const double norm_diff_sq = dx * dx + dy * dy + dz * dz;

//     if (norm_diff_sq < 1e-20) { 

//         solution = std::complex<double>(0.0, 0.0);
//         return;

//     } else {
        
//         const double norm_diff = std::sqrt(norm_diff_sq);
//         const double inv_norm_diff = 1.0 / norm_diff;
//         const double inv_norm_diff_sq = inv_norm_diff * inv_norm_diff;

//         const double val_cos = std::cos(wavenumber * norm_diff);
//         const double val_sin = std::sin(wavenumber * norm_diff);

//         double solution_real, solution_imag;
        
//         switch (EQUATION_FORMULATION) {
            
//             case 1: {

//                 solution_real = val_cos * INV_4_PI * inv_norm_diff;
//                 solution_imag = val_sin * INV_4_PI * inv_norm_diff;
//                 break;

//             }

//             case 2: {
                
//                 // Compute beta = < r, n > * |r|^-2

//                 const double rDotNorm = dx * nx + dy * ny + dz * nz;
//                 double beta = rDotNorm * inv_norm_diff_sq;

//                 if ((norm_diff_sq < 1.0E-10) && (std::abs(rDotNorm) < 1.0E-6)) {

//                     beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, 
//                                 nx, ny, nz, dxdsx, dxdsy, dxdsz, 
//                                 dxdtx, dxdty, dxdtz, dxdsdsx, dxdsdsy, dxdsdsz, 
//                                 dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

//                 }
                
//                 const double term1 = val_cos * inv_norm_diff;
//                 const double term2 = wavenumber * val_sin;
                
//                 solution_real = -INV_4_PI * INT_EXT *beta * (term1 + term2);

//                 const double term3 = val_sin * inv_norm_diff;
//                 const double term4 = wavenumber * val_cos;
                
//                 solution_imag = -INV_4_PI * INT_EXT * beta * (term3 - term4);
//                 break;

//             }
            
//             case 3: {

//                 const double solution_real_S = val_cos * INV_4_PI * inv_norm_diff;
//                 const double solution_imag_S = val_sin * INV_4_PI * inv_norm_diff;

//                 const double rDotNorm = dx * nx + dy * ny + dz * nz;
//                 double beta = rDotNorm * inv_norm_diff_sq;

//                 if ((norm_diff_sq < 1.0E-10) && (std::abs(rDotNorm) < 1.0E-6)) {

//                     beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, 
//                                 nx, ny, nz, dxdsx, dxdsy, dxdsz, 
//                                 dxdtx, dxdty, dxdtz, dxdsdsx, dxdsdsy, dxdsdsz, 
//                                 dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

//                 }

//                 const double term1 = val_cos * inv_norm_diff;
//                 const double term2 = wavenumber * val_sin;
//                 const double solution_real_D = -INV_4_PI * INT_EXT * beta * (term1 + term2);

//                 const double term3 = val_sin * inv_norm_diff;
//                 const double term4 = wavenumber * val_cos;
//                 const double solution_imag_D = -INV_4_PI * INT_EXT * beta * (term3 - term4);

//                 solution_real = solution_real_D + coupling_parameter * solution_imag_S;
//                 solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

//                 break;

//             }

//             default: {
                
//                 solution_real = 0.0;
//                 solution_imag = 0.0;
//                 break;

//             }

//         }

//         solution = std::complex<double>{solution_real, solution_imag};

//     }

// }



template<int PS, int PT>
void Solver<PS,PT>::HH(
    const double p1x, const double p1y, const double p1z,
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
    
    static const double INV_4_PI = 1.0 / (4.0 * M_PI);

    const double dx = p2x - p1x;
    const double dy = p2y - p1y;
    const double dz = p2z - p1z;
    
    const double norm_diff_sq = dx * dx + dy * dy + dz * dz;

    if (norm_diff_sq < 1e-20) { 

        solution = std::complex<double>(0.0, 0.0);
        return;

    } else {
        
        const double norm_diff = std::sqrt(norm_diff_sq);
        const double inv_norm_diff = 1.0 / norm_diff;
        const double inv_norm_diff_sq = inv_norm_diff * inv_norm_diff;

        const double val_cos = std::cos(wavenumber.real() * norm_diff);
        const double val_sin = std::sin(wavenumber.real() * norm_diff);
        const double val_eimagk = std::exp(-wavenumber.imag() * norm_diff);

        double solution_real, solution_imag;
        
        switch (EQUATION_FORMULATION) {
            
            case 1: {

                solution_real = val_eimagk * val_cos * INV_4_PI * inv_norm_diff;
                solution_imag = val_eimagk * val_sin * INV_4_PI * inv_norm_diff;
                break;

            }

            case 2: {
                
                // Compute beta = < r, n > * |r|^-2

                const double rDotNorm = dx * nx + dy * ny + dz * nz;
                double beta = rDotNorm * inv_norm_diff_sq;

                if ((norm_diff_sq < 1.0E-10) && (std::abs(rDotNorm) < 1.0E-6)) {

                    beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, 
                                nx, ny, nz, dxdsx, dxdsy, dxdsz, 
                                dxdtx, dxdty, dxdtz, dxdsdsx, dxdsdsy, dxdsdsz, 
                                dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

                }
                
                const double term1 = val_cos * inv_norm_diff;
                const double term2 =  wavenumber.real() * val_sin + wavenumber.imag() * val_cos;
                
                solution_real = -INV_4_PI * INT_EXT * beta * val_eimagk * (term1 + term2);

                const double term3 = val_sin * inv_norm_diff;
                const double term4 = wavenumber.real() * val_cos + wavenumber.imag();
                
                solution_imag = -INV_4_PI * INT_EXT * beta * val_eimagk * (term3 - term4);
                break;

            }
            
            case 3: {

                const double solution_real_S = val_eimagk * val_cos * INV_4_PI * inv_norm_diff;
                const double solution_imag_S = val_eimagk * val_sin * INV_4_PI * inv_norm_diff;

                const double rDotNorm = dx * nx + dy * ny + dz * nz;
                double beta = rDotNorm * inv_norm_diff_sq;

                if ((norm_diff_sq < 1.0E-10) && (std::abs(rDotNorm) < 1.0E-6)) {

                    beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, 
                                nx, ny, nz, dxdsx, dxdsy, dxdsz, 
                                dxdtx, dxdty, dxdtz, dxdsdsx, dxdsdsy, dxdsdsz, 
                                dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

                }

                const double term1 = val_cos * inv_norm_diff;
                const double term2 = wavenumber.real() * val_sin + wavenumber.imag() * val_cos;;
                const double solution_real_D = -INV_4_PI * INT_EXT * beta * (term1 + term2);

                const double term3 = val_sin * inv_norm_diff;
                const double term4 = wavenumber.real() * val_cos + wavenumber.imag();
                const double solution_imag_D = -INV_4_PI * INT_EXT * beta * (term3 - term4);

                // For clarity, use shorter names within this block:
                const double S_R = solution_real_S;
                const double S_I = solution_imag_S;
                const double D_R = solution_real_D;
                const double D_I = solution_imag_D;
                
                
                const double coupling_R = coupling_parameter.real();
                const double coupling_I = coupling_parameter.imag();


                // 1. Calculate the components of the product P = lambda * S
                // P_R = coupling_R * S_R - coupling_I * S_I
                // P_I = coupling_R * S_I + coupling_I * S_R
                const double product_R = coupling_R * S_R - coupling_I * S_I;
                const double product_I = coupling_R * S_I + coupling_I * S_R;

                // 2. Calculate final solution sol = D - i * P 
                // Since i * P = (-P_I + i * P_R), then D - i * P = (D_R - (-P_I)) + i * (D_I - P_R)
                // sol_R = D_R + P_I
                // sol_I = D_I - P_R
                solution_real = D_R + product_I;
                solution_imag = D_I - product_R;
                break;

            }

            default: {
                
                solution_real = 0.0;
                solution_imag = 0.0;
                break;

            }

        }

        solution = std::complex<double>{solution_real, solution_imag};

    }

}



template<int PS, int PT>
void Solver<PS,PT>::HH2(const double p1x, const double p1y, const double p1z,
         const double p2x, const double p2y, const double p2z,
         const double nx, const double ny, const double nz,
         std::complex<double> coupling_parameter, std::complex<double> wavenumber,
         std::complex<double>& solution)
{

    static const double INV_4_PI = 1.0 / (4.0 * M_PI);

    double solution_real, solution_imag;

    const double norm_diff = std::sqrt((p2x - p1x)*(p2x - p1x) + (p2y - p1y)*(p2y - p1y) + (p2z - p1z)*(p2z - p1z));

    const double inv_norm_diff = 1.0 / norm_diff;

    if (norm_diff < 1e-14) {

        solution_real = 0.0;
        solution_imag = 0.0;

    } else {   

        const double val_cos   = std::cos(wavenumber.real() * norm_diff);
        const double val_sin   = std::sin(wavenumber.real() * norm_diff);
        const double val_eimagk = std::exp(-wavenumber.imag() * norm_diff);

        if (EQUATION_FORMULATION == 1) {

            solution_real = val_eimagk * val_cos * INV_4_PI * inv_norm_diff;
            solution_imag = val_eimagk * val_sin * INV_4_PI * inv_norm_diff;

        } else if (EQUATION_FORMULATION == 2) {

            // Compute beta = < r, n > / |r|^2
            
            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            solution_real = -INV_4_PI * beta * INT_EXT  * val_eimagk * (val_cos * inv_norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
            solution_imag = -INV_4_PI * beta * INT_EXT  * val_eimagk * (val_sin * inv_norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

            if (USE_ACCELERATOR) {

                solution_real *= -1.0;
                solution_imag *= -1.0;

            }

        } else {

            // EQUATION_FORMULATION == 3 

            const double solution_real_S = val_eimagk * val_cos * INV_4_PI * inv_norm_diff;
            const double solution_imag_S = val_eimagk * val_sin * INV_4_PI * inv_norm_diff;

            // Compute beta = < r, n > / |r|^2

            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            double solution_real_D = -INV_4_PI * beta * INT_EXT * val_eimagk * (val_cos * inv_norm_diff + wavenumber.real() * val_sin + wavenumber.imag() * val_cos);
            double solution_imag_D = -INV_4_PI * beta * INT_EXT * val_eimagk * (val_sin * inv_norm_diff - wavenumber.real() * val_cos + wavenumber.imag() * val_sin);        

            if (USE_ACCELERATOR) {

                solution_real_D *= -1.0;
                solution_imag_D *= -1.0;

            }

            // For clarity, use shorter names within this block:
            const double S_R = solution_real_S;
            const double S_I = solution_imag_S;
            const double D_R = solution_real_D;
            const double D_I = solution_imag_D;
            
            
            const double coupling_R = coupling_parameter.real();
            const double coupling_I = coupling_parameter.imag();


            // 1. Calculate the components of the product P = lambda * S
            // P_R = coupling_R * S_R - coupling_I * S_I
            // P_I = coupling_R * S_I + coupling_I * S_R
            const double product_R = coupling_R * S_R - coupling_I * S_I;
            const double product_I = coupling_R * S_I + coupling_I * S_R;

            // 2. Calculate final solution sol = D - i * P 
            // Since i * P = (-P_I + i * P_R), then D - i * P = (D_R - (-P_I)) + i * (D_I - P_R)
            // sol_R = D_R + P_I
            // sol_I = D_I - P_R
            solution_real = D_R + product_I;
            solution_imag = D_I - product_R;

            // solution_real = solution_real_D + coupling_parameter * solution_imag_S;
            // solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

        }

        solution = std::complex<double>{solution_real, solution_imag};

    }

}

template<int PS, int PT>
void Solver<PS,PT>::HH_far(const double xVers_0, const double xVers_1, const double xVers_2,
            const double y_0, const double y_1, const double y_2,
            const double n_0, const double n_1, const double n_2,
            const double coupling_parameter, double wavenumber,
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

        solution_real = D_real + coupling_parameter * S_imag;
        solution_imag = D_imag - coupling_parameter * S_real;

    } 

    solution = std::complex<double>{solution_real, solution_imag};

}

template<int PS, int PT>
void Solver<PS,PT>::fct_4(const double x1, const double x2, const double x3,
                                 const double y1, const double y2, const double y3,
                                 const double normal1, const double normal2, const double normal3,
                                 const std::complex<double> coupling_parameter, const std::complex<double> wavenumber,
                                 const std::complex<double> density, 
                                 std::complex<double>& phi)
{

    std::complex<double> kernel;

    HH2(x1, x2, x3, 
        y1, y2, y3,
        normal1, normal2, normal3,
        coupling_parameter, wavenumber,
        kernel);

    const double kerreal = kernel.real();
    const double kerimag = kernel.imag();

    const double phireal = density.real() * kerreal - density.imag() * kerimag;
    const double phiimag = density.real() * kerimag + density.imag() * kerreal;

    phi = {phireal, phiimag};

}

template<int PS, int PT>
void Solver<PS,PT>::fac_1(const double distance, std::complex<double> wavenumber, std::complex<double>& sol) 
{

    const double dist_inv = 1.0 / distance;
    const double val_eimak = std::exp(-wavenumber.imag() * distance);
    const double re = val_eimak * std::cos(wavenumber.real() * distance) * dist_inv;
    const double im = val_eimak * std::sin(wavenumber.real() * distance) * dist_inv;

    sol = {re, im};

}

