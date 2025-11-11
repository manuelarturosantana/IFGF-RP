#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <string>
#include <vector>

// Set problem parameters

// Geometry parameters:

// Sphere = 0
// Others = 1

// const int GEOMETRY = 0;

// Sphere


// const double SPHERE_RADIUS = 1.0;
// const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};
// constexpr int N_PATCHES_ORIG = 6;
// const bool EDGE_FLAG_U_A[N_PATCHES_ORIG] = {false};
// const bool EDGE_FLAG_U_B[N_PATCHES_ORIG] = {false};
// const bool EDGE_FLAG_V_A[N_PATCHES_ORIG] = {false};
// const bool EDGE_FLAG_V_B[N_PATCHES_ORIG] = {false};

// const std::string DIRECTORY = "";
// const std::string FILE_NAME = "";


// Nacelle


// constexpr int N_PATCHES_ORIG = 134;
// const bool EDGE_FLAG_U_A[N_PATCHES_ORIG] = {false};
// const bool EDGE_FLAG_U_B[N_PATCHES_ORIG] = {false};
// const bool EDGE_FLAG_V_A[N_PATCHES_ORIG] = {false};
// const bool EDGE_FLAG_V_B[N_PATCHES_ORIG] = {false};

// const std::string DIRECTORY = "./Grids/Nacelle/";
// const std::string FILE_NAME = "Nacelle-";

// const double SPHERE_RADIUS = 0.0;
// const double SPHERE_CENTER[3] = {0.0};


// Submarine

/*
constexpr int N_PATCHES_ORIG = 285;
const bool EDGE_FLAG_U_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_U_B[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_A[N_PATCHES_ORIG] = {false};
const bool EDGE_FLAG_V_B[N_PATCHES_ORIG] = {false};

const std::string DIRECTORY = "./Grids/BigSubmarine/";
const std::string FILE_NAME = "BigSubmarine-";

const double SPHERE_RADIUS = 0.0;
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

const double SPHERE_RADIUS = 0.0;
const double SPHERE_CENTER[3] = {0.0};
*/

// Choose integral equation formulation:

// Single Layer = 1
// Double Layer = 2
// Combined Layer = 3
// INT_EXT = 1.0 or -1.0 for the exterior problem

// const int EQUATION_FORMULATION = 2;
// For the interior problem INT_EXT = -1.0
// For the exterior problem INT_EXT = 1.0
// const double INT_EXT = 1.0;


const int DELTA_METHOD = 2; // 1 uses PROXIMITY_BOX_SIZE and 2 uses PERCENT_BOX_SIZE
const double PROXIMITY_BOX_SIZE = 0.1/8;
const double PERCENT_BOX_SIZE = 0.15;

// GMRES options:

// const int MAX_ITER = 50;
// const double TOL_GMRES = 1E-6;

// Parameter related to the Colton-Kress change of variables
// only used in the change of variables file
const double P_CK_EDGE = 8.0;

// IFGF parameters:

// const bool USE_ACCELERATOR = true;

const int PS = 4; // Accuracy in the radial variable
const int PT = 6; // Accuracy in the theta variable
// const int PS = 6; // Accuracy in the radial variable
// const int PT = 8; // Accuracy in the theta variable

// OpenMP parameters:
const int NTHREADS = 1;

// The following parameters are the "experimental" part
// They are still being improved, you can use them if you want, but with the above things you are fine

// High order parameters:

const bool USE_HIGH_ORDER = false;
const int PS_II = 4; // Accuracy in theta
const int PT_II = 4; // Accuracy in phi
const double ACCURACY = 1E-17; // QR parameter (no need to change it)
const int TOTAL_COLS = -1; // Related with how many columns you want to take in QR (no need to change it)
const int TYPE_POINTS = 1; // 0 = Chebyshev, 1 = Uniform (no need to change it)

// Adaptivity parameters:

const bool USE_ADAPTIVITY = false;
const long long MAX_ELEMS_LEAF = 50;

// Oversampling parameters:
// (If you want to delete precomputations, but increases the time)
// Ideally good for very large problems (N >> 0)

const bool USE_OVERSAMPLING = false;
const int P_OVERSAMPLING = 2; // Size oversampling
const int METHOD_FFTW = 1; // No need to change this
const int TYPE_FFTW = 1; // No need to change this
constexpr int N_LOC_X = 5; // Size local interpolation in x
constexpr int N_LOC_Y = 5; // Size local interpolation in y

// Some parameters that are commented out in the code
// If you need to use them, uncomment the functions necessary for them


// Incident field parameters:
// 0 = PLANE WAVE, 1 = POINT SOURCE

// const int PLANE_OR_POINT = 0; 

//Incident plane wave parameters (PLANE_OR_POINT = 0):
// kx = k*cos(the)*sin(phi)
// ky = k*sin(the)*sin(phi)
// kz = k*cos(phi)


// const double PLANE_WAVE_THE = 0.0; // in [0, 2pi)
// const double PLANE_WAVE_PHI = M_PI; // in [0, pi]

/*
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
*/

#endif