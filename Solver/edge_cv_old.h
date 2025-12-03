#ifndef EDGE_CV_H
#define EDGE_CV_H

#include "global.h"
#include "parameters.h"

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