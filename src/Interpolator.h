#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Chebyshev.h"
#include <complex>


template <int PS, int PT>
class Interpolator
{

    private:

    //    int PS = PS;
    //    int PT = PT;

        const std::array<double, PS> xs_;
        const std::array<double, PT> xang_;

        const std::array<std::array<double, PS>, PS> Ts_; 
        const std::array<std::array<double, PT>, PT> Tang_;

        const int P_;

    public:

        constexpr Interpolator() : 
            xs_(Functions::SetupChebyshevPoints<PS>()),
            xang_(Functions::SetupChebyshevPoints<PT>()),
            Ts_(Functions::SetupMultipleTn<PS>()),
            Tang_(Functions::SetupMultipleTn<PT>()),
            // PS(PS),
            // PT(PT),
            P_(PS*PT*PT)
        {
        }

        void GetInterpolationPoint(const int iter, double& x, double& y, double& z) const
        {
            x = xs_[iter % PS];
            y = xang_[(iter / PS) % PT];
            z = xang_[iter / (PS*PT)];
        }

        constexpr int GetNInterpPoints() const {return P_;}

        constexpr int GetInterpId(const int is, const int itheta, const int iphi) const
        {
            return ( iphi*PT + itheta )*PS + is;
        }

        inline std::complex<double> Interpolate(const double x, const double y, const double z, const std::complex<double>* vals_begin) const
        {
            static thread_local std::array<double, PS> TargetTs; 
            static thread_local std::array<double, PT> TargetTphi;
            static thread_local std::array<double, PT> TargetTtheta;

            Functions::SetupMultipleTnArr<PS>(x, TargetTs);
            Functions::SetupMultipleTnArr<PT>(y, TargetTtheta);
            Functions::SetupMultipleTnArr<PT>(z, TargetTphi);

            std::complex<double> result = {0.0, 0.0};
            
            for (int sangiter2 = 0; sangiter2 < PT; sangiter2++) {
                for (int sangiter1 = 0; sangiter1 < PT; sangiter1++) {
                    for (int sraditer = 0; sraditer < PS; sraditer++) {
                
                        const long long liniter = (static_cast<long long>(sangiter2)*PT + sangiter1)*PS + sraditer;
                        const double cheb = TargetTphi[sangiter2]*TargetTtheta[sangiter1]*TargetTs[sraditer];                
                
                        result += std::complex<double>(cheb * vals_begin[liniter].real(), cheb * vals_begin[liniter].imag());
                    }
                }
            }

            return result;
        }

        void GenerateInterpolant(std::complex<double>* const arr) const
        {
            static thread_local std::array<std::array<std::array<std::complex<double>, PT>, PT>, PS> Fijo;
            static thread_local std::array<std::array<std::array<std::complex<double>, PT>, PT>, PS> Gino;

            #pragma unenroll
            for (int i = 0; i < PS; i++) {
                for (int j = 0; j < PT; j++) {
                    for (int o = 0; o < PT; o++) {
                        Fijo[i][j][o] = {0.0, 0.0};
                        Gino[i][j][o] = {0.0, 0.0};
                    }
                }
            }

            #pragma unenroll
            for (int i = 0; i < PS; i++) {
                for (int j = 0; j < PT; j++) {
                    for (int o = 0; o < PT; o++) {
                        for (int k = 0; k < PT; k++) {
                            Fijo[i][j][o] += std::complex<double>(arr[GetInterpId(i, j, k)].real()*Tang_[o][k], arr[GetInterpId(i, j, k)].imag()*Tang_[o][k]);
                        }
                    }
                }
            }

            #pragma unenroll
            for (int i = 0; i < PS; i++) {
                for (int n = 0; n < PT; n++) {
                    for (int o = 0; o < PT; o++) {
                        for (int j = 0; j < PT; j++) {
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
            for (int m = 0; m < PS; m++) {
                for (int n = 0; n < PT; n++) {
                    for (int o = 0; o < PT; o++) {
                        const int id = GetInterpId(m, n, o);
                        for (int i = 0; i < PS; i++) {
                            arr[id] += std::complex<double>(Gino[i][n][o].real()*Ts_[m][i], Gino[i][n][o].imag()*Ts_[m][i]);
                        }
                        double alphan = 2.0;
                        double alpham = 2.0;
                        double alphao = 2.0;
                        if (n == 0) alphan = 1.0;
                        if (m == 0) alpham = 1.0;
                        if (o == 0) alphao = 1.0;
                        arr[id] = {arr[id].real()*alpham*alphan*alphao/(PS*PT*PT), arr[id].imag()*alpham*alphan*alphao/(PS*PT*PT)};
                    }
                }
            }
        }

};

#endif