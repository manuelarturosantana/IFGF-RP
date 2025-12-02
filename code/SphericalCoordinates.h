#ifndef SPHERICALCOORDINATES_H
#define SPHERICALCOORDINATES_H

#include <cmath>

namespace Functions {

    inline double Atan2(const double y, const double x) noexcept {

        if (x > 0 && y >= 0)
            return atan(y/x);
        else if (x > 0 && y < 0)
            return atan(y/x) + 2.0*M_PI;
        else if (x < 0 && y > 0)
            return atan(y/x) + M_PI;
        else if (x < 0 && y == 0)
            return +M_PI;
        else if (x < 0 && y < 0)
            return atan(y/x) + M_PI;
        else if (x == 0 && y > 0)
            return M_PI_2;
        else if (x == 0 && y < 0)
            return 3.0*M_PI_2;
        else
            return 0;
            
    }

    inline void CartToSph(double& x, double& y, double& z) noexcept {

        const double rho = std::sqrt(x*x + y*y + z*z);

        if (rho == 0.0) {
            x = 0.0; y = 0.0; z = 0.0;
            return;
        }

        const double phi = Atan2(y, x);

        x = rho; // r
        y = acos(z/rho); // theta - 0 <= . <= pi
        z = phi; // phi - 0 <= . < 2pi

    }

    inline void SphToCart(double & r, double & theta, double & phi) noexcept {

        double sintheta, costheta;
        double sinphi, cosphi;
        sincos(theta, &sintheta, &costheta);
        sincos(phi, &sinphi, &cosphi);

        const double x = r * sintheta * cosphi;
        const double y = r * sintheta * sinphi;
        const double z = r * costheta;

        r = x;
        theta = y;
        phi = z;

    }

};

#endif