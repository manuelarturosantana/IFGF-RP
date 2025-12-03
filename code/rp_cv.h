#ifndef RP_CV_H
#define RP_CV_H

#include "global.h"
#include "parameters.h"

namespace {
    double TOL_PROJ_EDGE = 1.0E-10;
    double P_CK_EDGE = 8.0;
}


// Rectangular polar change of variables
// tau and alpha are in [-1, 1]

double inline xi(double alpha, double tau) 
{

    double solution;

    if (std::abs(alpha - 1.0) < TOL_PROJ_EDGE) {

        solution = alpha - ((1.0 + alpha) / M_PI) * fun_w(P_CK_EDGE, M_PI * std::abs((tau - 1.0) / 2.0));

    } else if (std::abs(alpha + 1.0) < TOL_PROJ_EDGE) {

        solution = alpha + ((1.0 - alpha) / M_PI) * fun_w(P_CK_EDGE, M_PI * std::abs((tau + 1.0) / 2.0));

    } else {

        solution = alpha + ((sign(tau) - alpha) / M_PI) * fun_w(P_CK_EDGE, M_PI * std::abs(tau));

    }

    return solution;

}

double inline der_xi(double alpha, double tau) 
{

    double solution;

    if (std::abs(alpha - 1.0) < TOL_PROJ_EDGE) {

        solution = - ((1.0 + alpha) / 2.0) * fun_dw(P_CK_EDGE, M_PI * std::abs((tau - 1.0) / 2.0)) * sign(tau - 1.0);

    } else if (std::abs(alpha + 1.0) < TOL_PROJ_EDGE) {

        solution = ((1.0 - alpha) / 2.0) * fun_dw(P_CK_EDGE, M_PI * std::abs((tau + 1.0) / 2.0)) * sign(tau + 1.0);

    } else {

        solution = (sign(tau) - alpha) * fun_dw(P_CK_EDGE, M_PI * std::abs(tau)) * sign(tau);

    }

    return solution;

}

// Edge change of variables
// Takes u in [-1, 1] and maps it to [-1, 1] clustered
// 0 = no edge, 1 = edge at -1, 2 = edge at 1, 3 = edge at -/+ 1

inline double eta(double u, int flag) 
{

    double solution;

    if (flag == 0) {

        solution = u;

    } else if (flag == 1) {

        solution = -1.0 + (2.0 / M_PI) * fun_w(P_CK_EDGE, (M_PI / 2.0) * (u + 1.0));

    } else if (flag == 2) {

        solution = -3.0 + (2.0 / M_PI) * fun_w(P_CK_EDGE, M_PI + (M_PI / 2.0) * (u + 1.0));

    } else {

        solution = -1.0 + (1.0 / M_PI) * fun_w(P_CK_EDGE, M_PI * (u + 1.0));

    }

    return solution;

}

inline double der_eta(double u, int flag) 
{

    double solution;

    if (flag == 0) {

        solution = 1.0;

    } else if (flag == 1) {

        solution = fun_dw(P_CK_EDGE, (M_PI / 2.0) * (u + 1.0));

    } else if (flag == 2) {

        solution = fun_dw(P_CK_EDGE, M_PI + (M_PI / 2.0) * (u + 1.0));

    } else {

        solution = fun_dw(P_CK_EDGE, M_PI * (u + 1.0));

    }

    return solution;

} 


#endif