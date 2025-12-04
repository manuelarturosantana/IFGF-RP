#include <array>
#include <cmath>


void  parametrization_q(const double SPHERE_RADIUS, const std::array<double, 3>& SPHERE_CENTER,
    const double s, const double t, const int patch, double& x, double& y, double& z)
{

    const double denominator = 1.0 / std::sqrt(1.0 + s*s + t*t);

    if (patch == 0) {

        x = s * denominator;
        y = t * denominator;
        z = 1.0 * denominator;
        
    } else if (patch == 1) {

        x = s * denominator;
        y = 1.0 * denominator;
        z = t * denominator;

    } else if (patch == 2) {

        x = 1.0 * denominator;
        y = t * denominator;
        z = s * denominator;

    } else if (patch == 3) {

        x = s * denominator;
        y = t * denominator;
        z = -1.0 * denominator;

    } else if (patch == 4) {

        x = s * denominator;
        y = -1.0 * denominator;
        z = t * denominator;

    } else {

        x = -1.0 * denominator;
        y = t * denominator;
        z = s * denominator;

    }

    x *= SPHERE_RADIUS;
    y *= SPHERE_RADIUS;
    z *= SPHERE_RADIUS;

    x += SPHERE_CENTER[0];
    y += SPHERE_CENTER[1];
    z += SPHERE_CENTER[2];

}

// Cartesian product d/ds x d/dt for sphere

void  cartesian_product_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& nx, double& ny, double& nz) 
{

    const double denominator = 1.0 / std::sqrt(1.0 + s*s + t*t);
    const double denominator_3 = denominator * denominator * denominator;

    double dxds0, dxds1, dxds2, dxdt0, dxdt1, dxdt2;
    int normal_direction;

    if (patch == 0) {

        dxds0 = (1.0 + t*t) * denominator_3; 
        dxds1 = (- t*s) * denominator_3; 
        dxds2 = (- s) * denominator_3; 

        dxdt0 = (- s*t) * denominator_3; 
        dxdt1 = (1.0 + s*s) * denominator_3; 
        dxdt2 = (- t) * denominator_3; 

        normal_direction = 1;

    } else if (patch == 1) {

        dxds0 = (1.0 + t*t) * denominator_3; 
        dxds1 = (- s) * denominator_3; 
        dxds2 = (- s*t) * denominator_3; 

        dxdt0 = (- s*t) * denominator_3; 
        dxdt1 = (- t) * denominator_3; 
        dxdt2 = (1.0 + s*s) * denominator_3; 

        normal_direction = -1;

    } else if (patch == 2) {

        dxds0 = (- s) * denominator_3; 
        dxds1 = (- s*t) * denominator_3; 
        dxds2 = (1.0 + t*t) * denominator_3; 

        dxdt0 = (- t) * denominator_3; 
        dxdt1 = (1.0 + s*s) * denominator_3; 
        dxdt2 = (- s*t) * denominator_3; 

        normal_direction = -1;

    } else if (patch == 3) {

        dxds0 = (1.0 + t*t) * denominator_3; 
        dxds1 = (- s*t) * denominator_3; 
        dxds2 = s * denominator_3; 

        dxdt0 = (- s*t) * denominator_3; 
        dxdt1 = (1.0 + s*s) * denominator_3; 
        dxdt2 = t * denominator_3; 

        normal_direction = -1;

    } else if (patch == 4) {

        dxds0 = (1.0 + t*t) * denominator_3; 
        dxds1 = s * denominator_3; 
        dxds2 = (- s*t) * denominator_3; 

        dxdt0 = (- s*t) * denominator_3; 
        dxdt1 = t * denominator_3; 
        dxdt2 = (1.0 + s*s) * denominator_3; 

        normal_direction = 1;

    } else {

        dxds0 = s * denominator_3; 
        dxds1 = (- s*t) * denominator_3; 
        dxds2 = (1.0 + t*t) * denominator_3; 

        dxdt0 = t * denominator_3; 
        dxdt1 = (1.0 + s*s) * denominator_3; 
        dxdt2 = (- s*t) * denominator_3; 

        normal_direction = 1;

    }

    dxds0 *= SPHERE_RADIUS;
    dxds1 *= SPHERE_RADIUS;
    dxds2 *= SPHERE_RADIUS;

    dxdt0 *= SPHERE_RADIUS;
    dxdt1 *= SPHERE_RADIUS;
    dxdt2 *= SPHERE_RADIUS;

    nx = (dxdt2*dxds1 - dxdt1*dxds2) * normal_direction;
    ny = (dxdt0*dxds2 - dxdt2*dxds0) * normal_direction;
    nz = (dxdt1*dxds0 - dxdt0*dxds1) * normal_direction;

}

// Jacobian

double  jacobian_q(const double SPHERE_RADIUS, const double s, const double t, const int patch)
{

    double nx, ny, nz;

    cartesian_product_q(SPHERE_RADIUS,s, t, patch, nx, ny, nz); 

    return std::sqrt(nx*nx + ny*ny + nz*nz);

}

// Normal vector

void  normal_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& nx, double& ny, double& nz) 
{

    cartesian_product_q(SPHERE_RADIUS, s, t, patch, nx, ny, nz);

    const double norm = std::sqrt(nx*nx + ny*ny + nz*nz);

    nx /= norm;
    ny /= norm;
    nz /= norm;

}

// d/ds

void  dxds_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdsx, double& dxdsy, double& dxdsz) 
{

    const double denominator = 1.0 / std::sqrt(1.0 + s*s + t*t);
    const double denominator_3 = denominator * denominator * denominator;

    if (patch == 0) {

        dxdsx = (1.0 + t*t) * denominator_3;
        dxdsy = (-s*t) * denominator_3;
        dxdsz = (-s) * denominator_3;

    } else if (patch == 1) {

        dxdsx = (1.0 + t*t) * denominator_3;
        dxdsy = (-s) * denominator_3;
        dxdsz = (-s*t) * denominator_3;

    } else if (patch == 2) {

        dxdsx = (-s) * denominator_3;
        dxdsy = (-s*t) * denominator_3;
        dxdsz = (1.0 + t*t) * denominator_3;

    } else if (patch == 3) {

        dxdsx = (1.0 + t*t) * denominator_3;
        dxdsy = (-s*t) * denominator_3;
        dxdsz = s * denominator_3;

    } else if (patch == 4) {

        dxdsx = (1.0 + t*t) * denominator_3;
        dxdsy = s * denominator_3;
        dxdsz = (-s*t) * denominator_3;

    } else {

        dxdsx = s * denominator_3;
        dxdsy = (-s*t) * denominator_3;
        dxdsz = (1.0 + t*t) * denominator_3;

    }

    dxdsx *= SPHERE_RADIUS;
    dxdsy *= SPHERE_RADIUS;
    dxdsz *= SPHERE_RADIUS;

}

// d/dt

void  dxdt_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdtx, double& dxdty, double& dxdtz) 
{

    const double denominator = 1.0 / std::sqrt(1.0 + s*s + t*t);
    const double denominator_3 = denominator * denominator * denominator;

    if (patch == 0) {

        dxdtx = (-s*t) * denominator_3;
        dxdty = (1.0 + s*s) * denominator_3;
        dxdtz = (-t) * denominator_3;

    } else if (patch == 1) {

        dxdtx = (-s*t) * denominator_3;
        dxdty = (-t) * denominator_3;
        dxdtz = (1.0 + s*s) * denominator_3;

    } else if (patch == 2) {

        dxdtx = (-t) * denominator_3;
        dxdty = (1.0 + s*s) * denominator_3;
        dxdtz = (-s*t) * denominator_3;

    } else if (patch == 3) {

        dxdtx = (-s*t) * denominator_3;
        dxdty = (1.0 + s*s) * denominator_3;
        dxdtz = t * denominator_3;

    } else if (patch == 4) {

        dxdtx = (-s*t) * denominator_3;
        dxdty = t * denominator_3;
        dxdtz = (1.0 + s*s) * denominator_3;

    } else {

        dxdtx = t * denominator_3;
        dxdty = (1.0 + s*s) * denominator_3;
        dxdtz = (-s*t) * denominator_3;

    }

    dxdtx *= SPHERE_RADIUS;
    dxdty *= SPHERE_RADIUS;
    dxdtz *= SPHERE_RADIUS;

}

// d^2/ds^2

void  dxdsds_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdsdsx, double& dxdsdsy, double& dxdsdsz) 
{

    const double denominator = 1.0 / std::sqrt(1.0 + s*s + t*t);
    const double denominator_5 = denominator * denominator * denominator * denominator * denominator;

    if (patch == 0) {

        dxdsdsx = (-3.0 * s * (1.0 + t*t)) * denominator_5;
        dxdsdsy = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;
        dxdsdsz = (-1.0 + 2.0 * s*s - t*t) * denominator_5;

    } else if (patch == 1) {

        dxdsdsx = (-3.0 * s * (1.0 + t*t)) * denominator_5;
        dxdsdsy = (-1.0 + 2.0 * s*s - t*t) * denominator_5;
        dxdsdsz = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;

    } else if (patch == 2) {

        dxdsdsx = (-1.0 + 2.0 * s*s - t*t) * denominator_5;
        dxdsdsy = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;
        dxdsdsz = (-3.0 * s * (1.0 + t*t)) * denominator_5;

    } else if (patch == 3) {

        dxdsdsx = (-3.0 * s * (1.0 + t*t)) * denominator_5;
        dxdsdsy = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;
        dxdsdsz = -(-1.0 + 2.0 * s*s - t*t) * denominator_5;

    } else if (patch == 4) {

        dxdsdsx = (-3.0 * s * (1.0 + t*t)) * denominator_5;
        dxdsdsy = (-(-1.0 + 2.0 * s*s - t*t)) * denominator_5;
        dxdsdsz = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;

    } else {

        dxdsdsx = (-(-1.0 + 2.0 * s*s - t*t)) * denominator_5;
        dxdsdsy = (-(t * (1.0 - 2.0 * s*s + t*t))) * denominator_5;
        dxdsdsz = (-3.0 * s * (1.0 + t*t)) * denominator_5;

    }

    dxdsdsx *= SPHERE_RADIUS;
    dxdsdsy *= SPHERE_RADIUS;
    dxdsdsz *= SPHERE_RADIUS;

}

// d^2/dt^2

void  dxdtdt_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdtdtx, double& dxdtdty, double& dxdtdtz) 
{

    const double denominator = 1.0 / std::sqrt(1.0 + s*s + t*t);
    const double denominator_5 = denominator * denominator * denominator * denominator * denominator;

    if (patch == 0) {

        dxdtdtx = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;
        dxdtdty = (-3.0 * (1.0 + s*s) * t) * denominator_5;
        dxdtdtz = (-1.0 - s*s + 2.0 * t*t) * denominator_5;

    } else if (patch == 1) {

        dxdtdtx = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;
        dxdtdty = (-1.0 - s*s + 2.0 * t*t) * denominator_5;
        dxdtdtz = (-3.0 * (1.0 + s*s) * t) * denominator_5;

    } else if (patch == 2) {

        dxdtdtx = (-1.0 - s*s + 2.0 * t*t) * denominator_5;
        dxdtdty = (-3.0 * (1.0 + s*s) * t) * denominator_5;
        dxdtdtz = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;

    } else if (patch == 3) {

        dxdtdtx = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;
        dxdtdty = (-3.0 * (1.0 + s*s) * t) * denominator_5;
        dxdtdtz = -(-1.0 - s*s + 2.0 * t*t) * denominator_5;

    } else if (patch == 4) {

        dxdtdtx = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;
        dxdtdty = (-(-1.0 - s*s + 2.0 * t*t)) * denominator_5;
        dxdtdtz = (-3.0 * (1.0 + s*s) * t) * denominator_5;

    } else {

        dxdtdtx = (-(-1.0 - s*s + 2.0 * t*t)) * denominator_5;
        dxdtdty = (-3.0 * (1.0 + s*s) * t) * denominator_5;
        dxdtdtz = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;

    }

    dxdtdtx *= SPHERE_RADIUS;
    dxdtdty *= SPHERE_RADIUS;
    dxdtdtz *= SPHERE_RADIUS;

}

// d^2/dsdt

void  dxdsdt_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdsdtx, double& dxdsdty, double& dxdsdtz) 
{

    const double denominator = 1.0 / std::sqrt(1.0 + s*s + t*t);
    const double denominator_5 = denominator * denominator * denominator * denominator * denominator;

    if (patch == 0) {

        dxdsdtx = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;
        dxdsdty = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;
        dxdsdtz = (3.0 * s * t) * denominator_5;

    } else if (patch == 1) {

        dxdsdtx = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;
        dxdsdty = (3.0 * s * t) * denominator_5;
        dxdsdtz = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;

    } else if (patch == 2) {

        dxdsdtx = (3.0 * s * t) * denominator_5;
        dxdsdty = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;
        dxdsdtz = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;

    } else if (patch == 3) {

        dxdsdtx = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;
        dxdsdty = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;
        dxdsdtz = (-3.0 * s * t) * denominator_5;

    } else if (patch == 4) {

        dxdsdtx = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;
        dxdsdty = (-3.0 * s * t) * denominator_5;
        dxdsdtz = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;

    } else {

        dxdsdtx = (-3.0 * s * t) * denominator_5;
        dxdsdty = (-s * (1.0 + s*s - 2.0 * t*t)) * denominator_5;
        dxdsdtz = (-t * (1.0 - 2.0 * s*s + t*t)) * denominator_5;

    }

    dxdsdtx *= SPHERE_RADIUS;
    dxdsdty *= SPHERE_RADIUS;
    dxdsdtz *= SPHERE_RADIUS;

}