#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <string>
#include <vector>

// Set problem parameters

// Geometry parameters:

// Sphere = 0
// Others = 1

const int GEOMETRY = 0;

// Sphere


const double SPHERE_RADIUS = 1.0;
const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};
constexpr long long N_PATCHES_ORIG = 6;
const bool EDGE_FLAG_U_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_U_B[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_B[N_PATCHES_ORIG] = {false};

const std::string DIRECTORY = "";
const std::string FILE_NAME = "";


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

// Choose integral equation formulation:

// Single Layer = 1
// Double Layer = 2
// Combined Layer = 3 
// 4

const int EQUATION_FORMULATION = 3;

// Discretization parameters:

// Number of points used per dimension per patch in every patch
constexpr long long N_PTS_PER_PATCH[2] = {6, 6};
// Number of patch partitions along first and second dimensions
constexpr long long N_SPLIT_PER_PATCH[2] = {16, 16};
 
// Singular integration parameters:

// Number of points used with the singular integration
constexpr int N_PTS_SING_INT[2] = {40, 40};
// Proximity distance that determines near singular integrals
const int DELTA_METHOD = 1;
const double PROXIMITY_BOX_SIZE = 0.1/16;
const double PERCENT_BOX_SIZE = 0.15;

// GMRES options:

const int MAX_ITER = 100;
const double TOL_GMRES = 1E-4;

// Incident field parameters:
// 0 = PLANE WAVE, 1 = POINT SOURCE

const int PLANE_OR_POINT = 0; 

// Incident plane wave parameters (PLANE_OR_POINT = 0):
// kx = k*cos(the)*sin(phi)
// ky = k*sin(the)*sin(phi)
// kz = k*cos(phi)

const double LAMBDA = 2.0 * SPHERE_RADIUS / 8.0; // 2.0 * M_PI / WAVE_NUMBER
const double WAVE_NUMBER = 2.0 * M_PI / LAMBDA;
const double PLANE_WAVE_THE = 0.0; // in [0, 2pi)
const double PLANE_WAVE_PHI = M_PI; // in [0, pi]

// Incident source points (PLANE_OR_POINT = 1):

const int NUM_POINT_SOURCES = 2;
const std::vector<std::vector<double>> POINT_SOURCE_CENTER = {{0.0, 0.0, 3.3},
                                                              {0.0, 3.3, 0.0}};

// Far and near field points:

const int N_FAR_PTS[2] = {200, 200};
const std::vector<std::vector<int>> N_NEAR_PTS = {{10, 10, 1},
                                                  {10, 1, 10}};
const std::vector<std::vector<double>> NEAR_FIELD_LIMITS = {{-12.0, 12.0, -12.0, 12.0, -25.0, -25.0},
                                                            {-12.0, 12.0, 0.0, 0.0, -25.0, 15.0}};

// Auxiliary parameters:

const double P_CK_EDGE = 8.0;

// IFGF parameters:

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