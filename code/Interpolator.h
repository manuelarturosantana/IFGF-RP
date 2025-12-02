#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Chebyshev.h"
#include "parameters.h"
#include <complex>

enum InterpolatorScheme 
{
    Chebyshev = 0
};

template <InterpolatorScheme scheme>
class Interpolator
{

    private:

        static constexpr int Ps = PS;
        static constexpr int Pang = PT;

        const std::array<double, Ps> xs_;
        const std::array<double, Pang> xang_;

        const std::array<std::array<double, Ps>, Ps> Ts_; 
        const std::array<std::array<double, Pang>, Pang> Tang_;

        const int P_;

    public:

        constexpr Interpolator() : 
            xs_(Functions::SetupChebyshevPoints<Ps>()),
            xang_(Functions::SetupChebyshevPoints<Pang>()),
            Ts_(Functions::SetupMultipleTn<Ps>()),
            Tang_(Functions::SetupMultipleTn<Pang>()),
            P_(Ps*Pang*Pang)
        {
        }

        void GetInterpolationPoint(const int iter, double& x, double& y, double& z) const
        {
            x = xs_[iter % Ps];
            y = xang_[(iter / Ps) % Pang];
            z = xang_[iter / (Ps*Pang)];
        }

        constexpr int GetNInterpPoints() const {return P_;}

        constexpr int GetInterpId(const int is, const int itheta, const int iphi) const
        {
            return ( iphi*Pang + itheta )*Ps + is;
        }

        inline std::complex<double> Interpolate(const double x, const double y, const double z, const std::complex<double>* vals_begin) const
        {
            static thread_local std::array<double, Ps> TargetTs; 
            static thread_local std::array<double, Pang> TargetTphi;
            static thread_local std::array<double, Pang> TargetTtheta;

            Functions::SetupMultipleTnArr<Ps>(x, TargetTs);
            Functions::SetupMultipleTnArr<Pang>(y, TargetTtheta);
            Functions::SetupMultipleTnArr<Pang>(z, TargetTphi);

            std::complex<double> result = {0.0, 0.0};
            
            for (int sangiter2 = 0; sangiter2 < Pang; sangiter2++) {
                for (int sangiter1 = 0; sangiter1 < Pang; sangiter1++) {
                    for (int sraditer = 0; sraditer < Ps; sraditer++) {
                
                        const long long liniter = (static_cast<long long>(sangiter2)*Pang + sangiter1)*Ps + sraditer;
                        const double cheb = TargetTphi[sangiter2]*TargetTtheta[sangiter1]*TargetTs[sraditer];                
                
                        result += std::complex<double>(cheb * vals_begin[liniter].real(), cheb * vals_begin[liniter].imag());
                    }
                }
            }

            return result;
        }

        void GenerateInterpolant(std::complex<double>* const arr) const
        {
            static thread_local std::array<std::array<std::array<std::complex<double>, Pang>, Pang>, Ps> Fijo;
            static thread_local std::array<std::array<std::array<std::complex<double>, Pang>, Pang>, Ps> Gino;

            #pragma unenroll
            for (int i = 0; i < Ps; i++) {
                for (int j = 0; j < Pang; j++) {
                    for (int o = 0; o < Pang; o++) {
                        Fijo[i][j][o] = {0.0, 0.0};
                        Gino[i][j][o] = {0.0, 0.0};
                    }
                }
            }

            #pragma unenroll
            for (int i = 0; i < Ps; i++) {
                for (int j = 0; j < Pang; j++) {
                    for (int o = 0; o < Pang; o++) {
                        for (int k = 0; k < Pang; k++) {
                            Fijo[i][j][o] += std::complex<double>(arr[GetInterpId(i, j, k)].real()*Tang_[o][k], arr[GetInterpId(i, j, k)].imag()*Tang_[o][k]);
                        }
                    }
                }
            }

            #pragma unenroll
            for (int i = 0; i < Ps; i++) {
                for (int n = 0; n < Pang; n++) {
                    for (int o = 0; o < Pang; o++) {
                        for (int j = 0; j < Pang; j++) {
                            Gino[i][n][o] += std::complex<double>(Fijo[i][j][o].real()*Tang_[n][j], Fijo[i][j][o].imag()*Tang_[n][j]);
                        }
                    }
                }
            }

            #pragma unenroll
            for (int i = 0; i < GetNInterpPoints(); i++) {
                arr[i] = {0.0, 0.0};
            }

            #pragma unenroll
            for (int m = 0; m < Ps; m++) {
                for (int n = 0; n < Pang; n++) {
                    for (int o = 0; o < Pang; o++) {
                        const int id = GetInterpId(m, n, o);
                        for (int i = 0; i < Ps; i++) {
                            arr[id] += std::complex<double>(Gino[i][n][o].real()*Ts_[m][i], Gino[i][n][o].imag()*Ts_[m][i]);
                        }
                        double alphan = 2.0;
                        double alpham = 2.0;
                        double alphao = 2.0;
                        if (n == 0) alphan = 1.0;
                        if (m == 0) alpham = 1.0;
                        if (o == 0) alphao = 1.0;
                        arr[id] = {arr[id].real()*alpham*alphan*alphao/(Ps*Pang*Pang), arr[id].imag()*alpham*alphan*alphao/(Ps*Pang*Pang)};
                    }
                }
            }
        }

};

#endif