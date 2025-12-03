#pragma once

#include "../solver2.h"


void Solver::HH(
    const double p1x, const double p1y, const double p1z,
    const double p2x, const double p2y, const double p2z,
    const double nx, const double ny, const double nz,
    const double dxdsx, const double dxdsy, const double dxdsz,
    const double dxdtx, const double dxdty, const double dxdtz,
    const double dxdsdsx, const double dxdsdsy, const double dxdsdsz,
    const double dxdsdtx, const double dxdsdty, const double dxdsdtz,
    const double dxdtdtx, const double dxdtdty, const double dxdtdtz,
    const double coupling_parameter, double wavenumber,
    std::complex<double>& solution)
{
    
    const double INV_4_PI = 1.0 / (4.0 * M_PI);

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

        const double val_cos = std::cos(wavenumber * norm_diff);
        const double val_sin = std::sin(wavenumber * norm_diff);

        double solution_real, solution_imag;
        
        switch (EQUATION_FORMULATION) {
            
            case 1: {

                solution_real = val_cos * INV_4_PI * inv_norm_diff;
                solution_imag = val_sin * INV_4_PI * inv_norm_diff;
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
                const double term2 = wavenumber * val_sin;
                
                solution_real = -INV_4_PI * beta * (term1 + term2);

                const double term3 = val_sin * inv_norm_diff;
                const double term4 = wavenumber * val_cos;
                
                solution_imag = -INV_4_PI * beta * (term3 - term4);
                break;

            }
            
            case 3: {

                const double solution_real_S = val_cos * INV_4_PI * inv_norm_diff;
                const double solution_imag_S = val_sin * INV_4_PI * inv_norm_diff;

                const double rDotNorm = dx * nx + dy * ny + dz * nz;
                double beta = rDotNorm * inv_norm_diff_sq;

                if ((norm_diff_sq < 1.0E-10) && (std::abs(rDotNorm) < 1.0E-6)) {

                    beta = -beta_double_layer_close(p1x, p1y, p1z, p2x, p2y, p2z, 
                                nx, ny, nz, dxdsx, dxdsy, dxdsz, 
                                dxdtx, dxdty, dxdtz, dxdsdsx, dxdsdsy, dxdsdsz, 
                                dxdsdtx, dxdsdty, dxdsdtz, dxdtdtx, dxdtdty, dxdtdtz);

                }

                const double term1 = val_cos * inv_norm_diff;
                const double term2 = wavenumber * val_sin;
                const double solution_real_D = -INV_4_PI * beta * (term1 + term2);

                const double term3 = val_sin * inv_norm_diff;
                const double term4 = wavenumber * val_cos;
                const double solution_imag_D = -INV_4_PI * beta * (term3 - term4);

                solution_real = solution_real_D + coupling_parameter * solution_imag_S;
                solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

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

void Solver::HH2(const double p1x, const double p1y, const double p1z,
         const double p2x, const double p2y, const double p2z,
         const double nx, const double ny, const double nz,
         double coupling_parameter, double wavenumber,
         std::complex<double>& solution)
{

    double solution_real, solution_imag;

    const double norm_diff = std::sqrt((p2x - p1x)*(p2x - p1x) + (p2y - p1y)*(p2y - p1y) + (p2z - p1z)*(p2z - p1z));

    if (norm_diff < 1e-14) {

        solution_real = 0.0;
        solution_imag = 0.0;

    } else {   

        const double val_cos = std::cos(wavenumber * norm_diff);
        const double val_sin = std::sin(wavenumber * norm_diff);

        if (EQUATION_FORMULATION == 1) {

            solution_real = val_cos / (4.0 * M_PI * norm_diff);
            solution_imag = val_sin / (4.0 * M_PI * norm_diff);

        } else if (EQUATION_FORMULATION == 2) {

            // Compute beta = < r, n > / |r|^2
            
            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            solution_real = -(1.0 / (4.0 * M_PI)) * beta * (val_cos / norm_diff + wavenumber * val_sin);
            solution_imag = -(1.0 / (4.0 * M_PI)) * beta * (val_sin / norm_diff - wavenumber * val_cos);        

            if (USE_ACCELERATOR) {

                solution_real *= -1.0;
                solution_imag *= -1.0;

            }

        } else {

            // EQUATION_FORMULATION == 3 

            const double solution_real_S = val_cos / (4.0 * M_PI * norm_diff);
            const double solution_imag_S = val_sin / (4.0 * M_PI * norm_diff);

            // Compute beta = < r, n > / |r|^2

            const double rDotNorm = (p2x - p1x) * nx + (p2y - p1y) * ny + (p2z - p1z) * nz;
            double beta = rDotNorm / (norm_diff * norm_diff);

            double solution_real_D = -(1.0 / (4.0 * M_PI)) * beta * (val_cos / norm_diff + wavenumber * val_sin);
            double solution_imag_D = -(1.0 / (4.0 * M_PI)) * beta * (val_sin / norm_diff - wavenumber * val_cos);        

            if (USE_ACCELERATOR) {

                solution_real_D *= -1.0;
                solution_imag_D *= -1.0;

            }

            solution_real = solution_real_D + coupling_parameter * solution_imag_S;
            solution_imag = solution_imag_D - coupling_parameter * solution_real_S;

        }

        solution = std::complex<double>{solution_real, solution_imag};

    }

}

void Solver::HH_far(const double xVers_0, const double xVers_1, const double xVers_2,
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

void Solver::fct_4(const double x1, const double x2, const double x3,
                                 const double y1, const double y2, const double y3,
                                 const double normal1, const double normal2, const double normal3,
                                 const double coupling_parameter, const double wavenumber,
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

void Solver::fac_1(const double distance, double wavenumber, std::complex<double>& sol) 
{

    const double re = std::cos(wavenumber * distance) / distance;
    const double im = std::sin(wavenumber * distance) / distance;

    sol = {re, im};

}

