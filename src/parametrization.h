#ifndef PARAMETRIZATION_H
#define PARAMETRIZATION_H

#include <array>

// Parametrization for sphere

void parametrization_q(const double SPHERE_RADIUS, const std::array<double,3>& SPHERE_CENTER, const double s, const double t, const int patch, double& x, double& y, double& z);
// Cartesian product d/ds x d/dt for sphere

void cartesian_product_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& nx, double& ny, double& nz);

// Jacobian

double jacobian_q(const double SPHERE_RADIUS, const double s, const double t, const int patch);

// Normal vector

void normal_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& nx, double& ny, double& nz);

// d/ds

void dxds_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdsx, double& dxdsy, double& dxdsz); 

// d/dt

void dxdt_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdtx, double& dxdty, double& dxdtdz); 

// d^2/ds^2

void dxdsds_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdsdsx, double& dxdsdsy, double& dxdsdsz); 

// d^2/dt^2

void dxdtdt_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdtdtx, double& dxdtdty, double& dxdtdtz); 

// d^2/dsdt

void dxdsdt_q(const double SPHERE_RADIUS, const double s, const double t, const int patch, double& dxdsdtx, double& dxdsdty, double& dxdsdtz); 

#endif