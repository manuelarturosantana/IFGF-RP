#ifndef GLOBAL_INTERPOLATOR_H
#define GLOBAL_INTERPOLATOR_H

#include <array>
#include <cmath>

// Equivalent to FFTW_REDFT10 / N
template <int N>
std::array<double, N> DCT_II_1D(const std::array<double, N>& x)
{

    std::array<double, N> y;

    for (int i = 0; i < N; i++) {

        y[i] = 0.0;

        for (int j = 0; j < N; j++) {

            y[i] += x[j] * cos(i * (j + 0.5) * M_PI / N);

        }

        y[i] *= 2.0 / N;

    }

    return y;

}

// Equivalent to FFTW_REDFT10 / (N * M)
template <int N, int M>
std::array<double, N*M> DCT_II_2D(const std::array<double, N*M>& x)
{
    
    std::array<double, N*M> y;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {

            y[i*M + j] = 0.0;

            for (int n = 0; n < N; n++) {
                for (int m = 0; m < M; m++) {

                    y[i*M + j] += x[n*M + m] * cos(i * (n + 0.5) * M_PI / N) * cos(j * (m + 0.5) * M_PI / M);

                }
            }

            y[i*M + j] *= 4.0 / (N*M);

        }
    }

    return y;

}

// Equivalent to FFTW_REDFT10 / (N * M * O)
template <int N, int M, int O>
std::array<double, N*M*O> DCT_II_3D(const std::array<double, N*M*O>& x)
{
    
    std::array<double, N*M*O> y;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < O; k++) {

                y[i*M*O + j*O + k] = 0.0;

                for (int n = 0; n < N; n++) {
                    for (int m = 0; m < M; m++) {
                        for (int o = 0; o < O; o++) {

                            y[i*M*O + j*O + k] += x[n*M*O + m*O + o] * cos(i * (n + 0.5) * M_PI / N) * cos(j * (m + 0.5) * M_PI / M) * cos(k * (o + 0.5) * M_PI / O);

                        }
                    }
                }

                y[i*M*O + j*O + k] *= 8.0 / (N*M*O);

            }
        }
    }

    return y;

}

// Equivalent to FFTW_REDFT01 * 0.5
template <int N>
std::array<double, N> DCT_III_1D(const std::array<double, N>& x)
{

    std::array<double, N> y;

    double factor;

    for (int i = 0; i < N; i++) {

        y[i] = 0.0;

        for (int j = 0; j < N; j++) {

            if (j == 0) {
                factor = 0.5;
            } else {
                factor = 1.0;
            }

            y[i] += factor * x[j] * cos(j * (i + 0.5) * M_PI / N);

        }

    }

    return y;

}

// Equivalent to FFTW_REDFT01 * 0.5 * 0.5
template <int N, int M>
std::array<double, N*M> DCT_III_2D(const std::array<double, N*M>& x)
{

    std::array<double, N*M> y;

    double factor_n, factor_m;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {

            y[i*M + j] = 0.0;

            for (int n = 0; n < N; n++) {
                for (int m = 0; m < M; m++) {

                    if (n == 0) {
                        factor_n = 0.5;
                    } else {
                        factor_n = 1.0;
                    }

                    if (m == 0) {
                        factor_m = 0.5;
                    } else {
                        factor_m = 1.0;
                    }

                    y[i*M + j] += factor_n * factor_m * x[n*M + m] * cos(n * (i + 0.5) * M_PI / N) * cos(m * (j + 0.5) * M_PI / M);

                }
            }

        }
    }

    return y;

}

// Equivalent to FFTW_REDFT01 * 0.5 * 0.5 * 0.5
template <int N, int M, int O>
std::array<double, N*M*O> DCT_III_3D(const std::array<double, N*M*O>& x)
{

    std::array<double, N*M*O> y;

    double factor_n, factor_m, factor_o;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < O; k++) {

                y[i*M*O + j*O + k] = 0.0;

                for (int n = 0; n < N; n++) {
                    for (int m = 0; m < M; m++) {
                        for (int o = 0; o < O; o++) {

                            if (n == 0) {
                                factor_n = 0.5;
                            } else {
                                factor_n = 1.0;
                            }

                            if (m == 0) {
                                factor_m = 0.5;
                            } else {
                                factor_m = 1.0;
                            }

                            if (o == 0) {
                                factor_o = 0.5;
                            } else {
                                factor_o = 1.0;
                            }

                            y[i*M*O + j*O + k] += factor_n * factor_m * factor_o * x[n*M*O + m*O + o] * cos(n * (i + 0.5) * M_PI / N) * cos(m * (j + 0.5) * M_PI / M) * cos(o * (k + 0.5) * M_PI / O);

                        }
                    }
                }

            }
        }
    }

    return y;

}

template <typename T, int N>
void diffDiv(const std::array<double, N>& x, std::array<T, N>& f)
{

    for (int i = 1; i < N; i++) {
        for (int j = N-1; j >= i; j--) {
            f[j] = (f[j] - f[j-1]) / (x[j] - x[j-i]);
        }
    }

}

template <typename T, int N>
T NewtonInterp(const std::array<double, N>& x, const std::array<T, N>& coeffs, double X)
{

    double products = 1.0;
    T result = T(0);

    for (int i = 0; i < N; i++) {

        result += coeffs[i] * products;
        products *= (X - x[i]);

    }

    return result;

}

template <typename T, int N, int M>
T LocalInterpNewton2D(const std::array<double, N>& x, const std::array<double, M>& y, const std::array<T, N*M>& f, double x0, double y0) 
{

    std::array<T, M> coeffs_1;
    std::array<T, N> coeffs_2;

    for (int i = 0; i < N; i++) {

        for (int j = 0; j < M; j++) {

            coeffs_1[j] = f[i*M + j];

        }

        diffDiv<T, M>(y, coeffs_1);

        coeffs_2[i] = NewtonInterp<T, M>(y, coeffs_1, y0);

    }

    diffDiv<T, N>(x, coeffs_2);

    return NewtonInterp<T, N>(x, coeffs_2, x0);

}

#endif