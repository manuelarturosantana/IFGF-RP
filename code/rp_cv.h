#ifndef RP_CV_H
#define RP_CV_H

#include "global.h"
#include "parameters.h"

namespace {
    double TOL_PROJ_EDGE = 1.0E-10;
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

#endif