#ifndef GLOBAL_H
#define GLOBAL_H

#include <functional>
#include <cmath>
#include <array>
#include <iostream>
#include <cstdlib>
#include <complex>
#include <sys/time.h>
#include <sys/resource.h>
#include <fstream>
#include <unistd.h>

#include <mkl.h>
#include <omp.h>
#include <mpi.h>

namespace {
    double EQUAL_TOL = 1.0e-12;
}


// Beta double layer close

inline double beta_double_layer_close(const double p1x, const double p1y, const double p1z,
                                const double p2x, const double p2y, const double p2z,
                                const double nx, const double ny, const double nz,
                                const double dxdsx, const double dxdsy, const double dxdsz,
                                const double dxdtx, const double dxdty, const double dxdtz,
                                const double dxdsdsx, const double dxdsdsy, const double dxdsdsz,
                                const double dxdsdtx, const double dxdsdty, const double dxdsdtz,
                                const double dxdtdtx, const double dxdtdty, const double dxdtdtz)
{

    // Compute cos(theta) / |x-x'|

    const double gamma_u = dxdsx * (p1x - p2x) + dxdsy * (p1y - p2y) + dxdsz * (p1z - p2z);
    const double gamma_v = dxdtx * (p1x - p2x) + dxdty * (p1y - p2y) + dxdtz * (p1z - p2z);

    const double numerator = gamma_v * (dxdsx * dxdsx + dxdsy * dxdsy + dxdsz * dxdsz) - gamma_u * (dxdsx * dxdtx + dxdsy * dxdty + dxdsz * dxdtz);
    const double denominator = gamma_u * (dxdtx * dxdtx + dxdty * dxdty + dxdtz * dxdtz) - gamma_v * (dxdsx * dxdtx + dxdsy * dxdty + dxdsz * dxdtz);

    const double phi = std::atan2(numerator, denominator);

    double cos_phi, sin_phi;
    sin_phi = std::sin(phi);
    cos_phi = std::cos(phi);

    const double numerator2 = nx * (cos_phi * cos_phi * dxdsdsx + 2.0 * cos_phi * sin_phi * dxdsdtx + sin_phi * sin_phi * dxdtdtx) + 
                              ny * (cos_phi * cos_phi * dxdsdsy + 2.0 * cos_phi * sin_phi * dxdsdty + sin_phi * sin_phi * dxdtdty) + 
                              nz * (cos_phi * cos_phi * dxdsdsz + 2.0 * cos_phi * sin_phi * dxdsdtz + sin_phi * sin_phi * dxdtdtz);
    const double denominator2 = (cos_phi * dxdsx + sin_phi * dxdtx) * (cos_phi * dxdsx + sin_phi * dxdtx) + 
                                (cos_phi * dxdsy + sin_phi * dxdty) * (cos_phi * dxdsy + sin_phi * dxdty) + 
                                (cos_phi * dxdsz + sin_phi * dxdtz) * (cos_phi * dxdsz + sin_phi * dxdtz);

    return 0.5 * numerator2 / denominator2;

}

// Sign function

inline double sign(double x) 
{

    return (x > 0.0) - (x < 0.0);

}

// Send x in [a, b] to [c, d]

inline double ab2cd(const double a, const double b, const double c, const double d, const double x)
{

    return ((d - c) / (b - a)) * (x - a) + c;

}

// Fejer 1st quadrature rule

template <int N>
inline void fejerquadrature1(std::vector<double>& nodes, std::vector<double>& weights)
{

    for (int i = 0; i < N; i++) {

        nodes[i] = std::cos(M_PI * (2.0 * i + 1.0) / (2.0 * N));
        weights[i] = 0.0;

        for (int j = 1; j <= std::floor(N * 0.5); j++) {

            weights[i] += std::cos(j * M_PI * (2.0 * i + 1.0) / N) / (4.0 * j*j - 1.0);

        }

        weights[i] = (2.0 / N) * (1.0 - 2.0 * weights[i]);

    }

}

// Colton-Kress change of variables
// 0 <= tau <= 2pi

inline double fun_v(double p, double tau)
{

    const double quotient = (M_PI - tau) / M_PI;

    return ((1.0 / p) - 0.5) * quotient*quotient*quotient - (1.0 / p) * quotient + 0.5;

}

inline double fun_dv(double p, double tau)
{

    const double quotient = (M_PI - tau) / M_PI;

    return ((1.0 / p) - 0.5) * (-3.0 / M_PI) * quotient*quotient + (1.0 / p) * (1.0 / M_PI);

}

inline double fun_w(double p, double tau) 
{        

    const double fun_v_p = std::pow(fun_v(p, tau), p);
    const double fun_v_p_2pi = std::pow(fun_v(p, 2.0 * M_PI - tau), p);

    return 2.0 * M_PI * fun_v_p / (fun_v_p + fun_v_p_2pi);

}

inline double fun_dw(double p, double tau)
{

    const double funv = fun_v(p, tau);
    const double funv2pi = fun_v(p, 2.0 * M_PI - tau);

    const double fun_v_p = std::pow(funv, p-1);
    const double fun_v_p_2pi = std::pow(funv2pi, p-1);

    return 2.0 * M_PI * (p * fun_v_p * fun_v_p_2pi * fun_dv(p, tau) * (funv + funv2pi)) / ((fun_v_p*funv + fun_v_p_2pi*funv2pi) * (fun_v_p*funv + fun_v_p_2pi*funv2pi));

}

// Chebyshev evaluations

template<int N, int M>
inline void cheb_evals(const std::vector<double>& x, std::vector<std::complex<double>>& evals)
{

    if (M == 1) {

        for (int i = 0; i < N; i++) {

            evals[i] = std::complex<double>(1.0, 0.0);

        }

    } else if (M == 2) {

        for (int i = 0; i < N; i++) {

            evals[i] = std::complex<double>(1.0, 0.0);
            evals[N + i] = std::complex<double>(x[i], 0.0);

        }

    } else {

        for (int i = 0; i < N; i++) {

            evals[i] = std::complex<double>(1.0, 0.0);
            evals[N + i] = std::complex<double>(x[i], 0.0);

            for (int j = 2; j < M; j++) {

                evals[j*N + i] = 2.0 * x[i] * evals[(j-1)*N + i] - evals[(j-2)*N + i];

            }

        }    

    }

}

// Golden search algorithm

inline void golden_search(const std::function<double (double, double)>& f, int max_iter, double tol, double u_a, double u_b, double v_a, double v_b, double& sol_x, double& sol_y)
{

    double g_ratio = (1.0 + sqrt(5.0)) / 2.0;

    double u_c, u_d, v_c, v_d;
    double minimum, values[4];
    int index, count = 0;
    int act_u, act_v;

    bool condition = (count < max_iter) && (std::abs(u_b - u_a) > tol) && (std::abs(v_b - v_a) > tol);

    while (condition) {
        
        u_c = u_b - (u_b - u_a) / g_ratio;
        u_d = u_a + (u_b - u_a) / g_ratio;

        v_c = v_b - (v_b - v_a) / g_ratio;
        v_d = v_a + (v_b - v_a) / g_ratio;

        values[0] = f(u_c, v_c);
        values[1] = f(u_c, v_d);
        values[2] = f(u_d, v_c);
        values[3] = f(u_d, v_d);

        minimum = values[0];
        index = 0;

        for (int i = 1; i < 4; i++) {

            if (values[i] < minimum) {

                minimum = values[i];
                index = i;

            }

        }

        act_u = index / 2;
        act_v = index % 2;

        u_a = act_u * u_c + (1 - act_u) * u_a;
        u_b = (1 - act_u) * u_d + act_u * u_b;
        v_a = act_v * v_c + (1 - act_v) * v_a;
        v_b = (1 - act_v) * v_d + act_v * v_b;

        count++;

        condition = (count < max_iter) && (std::abs(u_b - u_a) > tol) && (std::abs(v_b - v_a) > tol);
    
    }

    sol_x = 0.5 * (u_b + u_a); 
    sol_y = 0.5 * (v_b + v_a);

}

inline std::vector<double> barycentric_weights(const std::vector<double>& nodes)
{

    int N = nodes.size();

    std::vector<double> weights(N, 1.0);
    
    for (int j = 0; j < N; j++) {

        double w = 1.0;

        for (int i = 0; i < N; i++) {

            if (i != j) {

                w *= (nodes[j] - nodes[i]);

            }

        }

        weights[j] = 1.0 / w;

    }

    return weights;

}

inline void lagrange_interpolation_2D(const std::vector<double>& xNodes, const std::vector<double>& yNodes,
                               const std::vector<double>& xWeights, const std::vector<double>& yWeights,
                               const std::vector<double>& fNodes, 
                               const std::vector<double>& x, const std::vector<double>& y,
                               double * solution)
{

    const int N = xNodes.size();
    const int M = yNodes.size();

    const int P = x.size();
    const int Q = y.size();

    std::vector<double> lx(P, 1.0);
    std::vector<double> ly(Q, 1.0);

    std::vector<int> idx_equal_x(P, -1);
    std::vector<int> idx_equal_y(Q, -1);

    std::vector<double> vec_x(P*N);
    std::vector<double> vec_y(M*Q);

    for (int i = 0; i < P; ++i) {

        double prod = 1.0;
        int eq_idx = -1;

        for (int j = 0; j < N; ++j) {

            double diff = x[i] - xNodes[j];

            if (std::abs(diff) < EQUAL_TOL) {

                eq_idx = j;
                break;

            }

            vec_x[i * N + j] = xWeights[j] / diff;

            prod *= diff;

        }

        if (eq_idx != -1) {

            std::fill(vec_x.begin() + i*N, vec_x.begin() + (i+1)*N, 0.0);
            vec_x[i*N + eq_idx] = 1.0;
            lx[i] = 1.0;

        } else {

            lx[i] = prod;

        }

    }

    for (int i = 0; i < Q; ++i) {

        double prod = 1.0;
        int eq_idx = -1;

        for (int j = 0; j < M; ++j) {

            double diff = y[i] - yNodes[j];

            if (std::abs(diff) < EQUAL_TOL) {

                eq_idx = j;
                break;

            }

            vec_y[i*M + j] = yWeights[j] / diff;

            prod *= diff;
        }

        if (eq_idx != -1) {

            std::fill(vec_y.begin() + i*M, vec_y.begin() + (i+1)*M, 0.0);
            vec_y[i*M + eq_idx] = 1.0;
            ly[i] = 1.0;

        } else {

            ly[i] = prod;

        }

    }

    std::vector<double> mat_prod(P*M);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, P, M, N, 1.0, &vec_x[0], N, &fNodes[0], M, 0.0, &mat_prod[0], M);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, P, Q, M, 1.0, &mat_prod[0], M, &vec_y[0], M, 0.0, solution, Q);

    for (int i = 0; i < P; i++) {
        for (int j = 0; j < Q; j++) {

            solution[i * Q + j] *= lx[i]*ly[j];

        }
    }

}

inline void print_max_RSS() 
{
    
    rusage ru {};
    getrusage(RUSAGE_SELF, &ru);
    double peak_gb = (double)ru.ru_maxrss / (1024.0 * 1024.0);

    double current_gb;
    std::ifstream f("/proc/self/status");
    std::string key;
    long value_kb = 0;
    while (f >> key) {
        if (key == "VmRSS:") {
            f >> value_kb;
            current_gb = value_kb / (1024.0 * 1024.0);
        }
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    std::cout << "Peak RSS    = " << peak_gb    << " GB\n";
    std::cout << "Current RSS = " << current_gb << " GB\n";

}

inline double minusLog10JnInfty(int n, double x) {

    double ediv2 = 0.5 * std::exp(1.0);

    return 0.5 * std::log10(M_2_PI * x) - n * std::log10(ediv2 * x / n);

}

inline int bwdRecurrenceStartJn1(double x, int p) {

    double xabs = std::abs(x);

    int nkm1 = static_cast<int>(1.1 * xabs) + 1;
    
    double fkm1 = minusLog10JnInfty(nkm1, xabs) - p;

    int nk = nkm1 + 5;
    double fk = minusLog10JnInfty(nk, xabs) - p;

    int nkp1;

    for (int m = 1; m <= 20; m++) {

        nkp1 = nk - (nk - nkm1) / (1.0 - fkm1 / fk);

        double f = minusLog10JnInfty(nkp1, xabs) - p;

        if (std::abs(nkp1 - nk) < 1) {

            break;

        }

        nkm1 = nk;
        nk = nkp1;

        fkm1 = fk;
        fk = f;

    }

    return nkp1;

}

inline int bwdRecurrenceStartJn2(double x, int n, int p) {

    double xabs = std::abs(x);
    double halfP = 0.5 * p;

    double fn = minusLog10JnInfty(n, xabs);
    double q;

    if (fn <= halfP) {

        return bwdRecurrenceStartJn1(x, p) + 10;

    } else {

        q = halfP + fn;

    }

    int nkm1 = static_cast<int>(1.1 * xabs) + 1;    
    double fkm1 = minusLog10JnInfty(nkm1, xabs) - q;

    int nk = nkm1 + 5;
    double fk = minusLog10JnInfty(nk, xabs) - q;

    int nkp1;

    for (int m = 1; m <= 20; m++) {

        nkp1 = nk - (nk - nkm1) / (1.0 - fkm1 / fk);

        double f = minusLog10JnInfty(nkp1, xabs) - q;

        if (std::abs(nkp1 - nk) < 1) {

            break;

        }

        nkm1 = nk;
        nk = nkp1;

        fkm1 = fk;
        fk = f;

    }

    return nkp1 + 10;

}

inline void SphericalBesselYAndDerivative_Real(int n, double x, std::vector<double>& sY, std::vector<double>& dY) {

    int nm = n;

    if (x < 1.0e-60) {

        for (int k = 0; k <= n; k++) {

            sY[k] = -1.0e300;
            dY[k] = 1.0e300;

        }

        return;

    }

    sY[0] = -cos(x) / x;
    sY[1] = (sY[0] - sin(x)) / x;

    double f0 = sY[0];
    double f1 = sY[1];
    double f;

    for (int k = 2; k <= n; k++) {

        f = (2.0 * k - 1.0) * f1 / x - f0;
        sY[k] = f;
        nm = k-1;

        if (1.0e300 <= std::abs(f)) {

            break;

        }

        f0 = f1;
        f1 = f;

    }
    
    dY[0] = (sin(x) + cos(x) / x) / x;

    for (int k = 1; k <= nm; k++) {

        dY[k] = sY[k-1] - (k + 1.0) * sY[k] / x;

    }

}

inline void SphericalBesselJAndDerivative_Real(int n, double x, std::vector<double>& sj, std::vector<double>& dj) {

    int nm = n;

    if (std::abs(x) <= 1e-100) {

        for (int k = 0; k <= n; k++) {

            sj[k] = 0.0;
            dj[k] = 0.0;

        }

        sj[0] = 1.0;
        dj[1] = 1.0 / 3.0;

        return;

    }

    sj[0] = sin(x) / x;
    sj[1] = (sj[0] - cos(x)) / x;

    if (n >= 2) {

        double sa = sj[0];
        double sb = sj[1];
        int m = bwdRecurrenceStartJn1(x, 200);

        if (m < n) {

            nm = m;

        } else {

            m = bwdRecurrenceStartJn2(x, n, 15);

        }

        double f0 = 0.0;
        double f1 = 1.0e-100;
        double f;

        for (int k = m; k >= 0; k--) {

            f = (2.0 * k + 3.0) * f1 / x - f0;

            if (k <= nm) {

                sj[k] = f;

            }

            f0 = f1;
            f1 = f;

        }

        double cs;

        if (std::abs(sa) <= std::abs(sb)) {

            cs = sb / f0;

        } else {

            cs = sa / f;

        }

        for (int k = 0; k <= nm; k++) {

            sj[k] = cs * sj[k];

        }

    }

    dj[0] = (cos(x) - sin(x) / x) / x;

    for (int k = 1; k <= nm; k++) {

        dj[k] = sj[k-1] - (k + 1.0) * sj[k] / x;

    }

}

inline double sphericalBesselJn_Real(int n, double x) {

    double solution;

    if (std::abs(x) < 2.0 * std::numeric_limits<double>::epsilon()) {

        if (n == 0) {

            return 1.0;

        }

        double den = 1.0;

        for (int i = 1; i <= n; i++) {

            den *= (2.0 * i + 1.0);

        }

        solution = std::pow(x, n) / den;

    } else {

        if (n == 0) {

            solution = sin(x) / x;

        } else if (n == 1) {

            solution = sin(x) / (x*x) - cos(x) / x;

        } else {

            std::vector<double> sj(n+1), dj(n+1);

            SphericalBesselJAndDerivative_Real(n, x, sj, dj);

            solution = sj[n];

        }

    }

    return solution;

}

inline double sphericalBesselYn_Real(int n, double x) {

    double solution;

    if (std::abs(x) < 2.0 * std::numeric_limits<double>::epsilon()) {

        if (n == 0) {

            return -1.0 / x;

        }

        double num = 1.0;

        for (int i = 1; i <= n; i++) {

            num *= (2.0 * i - 1.0);

        }

        solution = -num / std::pow(x, n+1);

    } else {

        if (n == 0) {

            solution = -cos(x) / x;

        } else if (n == 1) {

            solution = -cos(x) / (x*x) - sin(x) / x;

        } else {

            std::vector<double> sj(n+1), dj(n+1);

            SphericalBesselYAndDerivative_Real(n, x, sj, dj);

            solution = sj[n];

        }

    }

    return solution;

}

inline std::complex<double> sphericalBesselH(int n, int k, double x) {

    std::complex<double> solution;

    if (k == 1) {

        if (n == 0) {

            solution = {sin(x) / x, -cos(x) / x};

        } else if (n == 1) {

            solution = {(sin(x) - x * cos(x)) / (x*x), - (x * sin(x) + cos(x)) / (x*x)};

        } else {

            solution = sphericalBesselJn_Real(n, x) + std::complex<double>(0.0, 1.0) * sphericalBesselYn_Real(n, x);

        }

    } else {

        if (n == 0) {

            solution = {sin(x) / x, cos(x) / x};

        } else if (n == 1) {

            solution = {(sin(x) - x * cos(x)) / (x*x), (x * sin(x) + cos(x)) / (x*x)};

        } else {

            solution = sphericalBesselJn_Real(n, x) - std::complex<double>(0.0, 1.0) * sphericalBesselYn_Real(n, x);

        }
        
    }   

    return solution;

}

#endif