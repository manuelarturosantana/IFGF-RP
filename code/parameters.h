#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <string>
#include <vector>



// Set problem parameters


        constexpr long long N_PTS_PER_PATCH[2] = {6, 6};
        // Number of patch partitions along first and second dimensions
        constexpr long long N_SPLIT_PER_PATCH[2] = {16, 16};
        
        // Singular integration parameters:

        // Number of points used with the singular integration
        constexpr int N_PTS_SING_INT[2] = {40, 40};
        constexpr long long N_PATCHES_ORIG = 6;

// Nacelle

/*
constexpr int N_PATCHES_ORIG = 134;
const bool EDGE_FLAG_U_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_U_B[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_B[N_PATCHES_ORIG] = {false};

const std::string DIRECTORY = "./Grids/Nacelle/";
const std::string FILE_NAME = "Nacelle-";

const double SPHERE_RADIUS = 1.0;
const double SPHERE_CENTER[3] = {0.0};
*/

// Submarine

/*
constexpr int N_PATCHES_ORIG = 285;
const bool EDGE_FLAG_U_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_U_B[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_B[N_PATCHES_ORIG] = {false};

const std::string DIRECTORY = "./Grids/BigSubmarine/";
const std::string FILE_NAME = "BigSubmarine-";

const double SPHERE_RADIUS = 1.0;
const double SPHERE_CENTER[3] = {0.0};
*/

// Turbofan

/*
constexpr int N_PATCHES_ORIG = 198;
const bool EDGE_FLAG_U_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_U_B[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_B[N_PATCHES_ORIG] = {false};

const std::string DIRECTORY = "./Grids/Turbofan/";
const std::string FILE_NAME = "Turbofan-";
*/


// Far and near field points:
// TODO: Can we move these out as well
const int N_FAR_PTS[2] = {200, 200};
const std::vector<std::vector<int>> N_NEAR_PTS = {{10, 10, 1},
                                                  {10, 1, 10}};
const std::vector<std::vector<double>> NEAR_FIELD_LIMITS = {{-12.0, 12.0, -12.0, 12.0, -25.0, -25.0},
                                                            {-12.0, 12.0, 0.0, 0.0, -25.0, 15.0}};

// Auxiliary parameters:

const double P_CK_EDGE = 8.0;

// IFGF parameters:
// TODO: Does N LEVELS IFGF still need to be chosen with the adaptive method?
const bool USE_ACCELERATOR = true;
const int N_LEVELS_IFGF = 5;

const int PS = 3;
const int PT = 5;

const double BBSIZEOFFSET = 0.0;

// Adaptivity parameters:

const bool USE_ADAPTIVITY = false;
const long long MAX_ELEMS_LEAF = 100;

// OpenMP parameters:

const int NTHREADS = 16;

#endif