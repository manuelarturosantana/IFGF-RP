#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <string>
#include <vector>





//const std::string DIRECTORY = "./Grids/Nacelle/";
//const std::string FILE_NAME = "Nacelle-";


//const std::string DIRECTORY = "./Grids/BigSubmarine/";
//const std::string FILE_NAME = "BigSubmarine-";

//const std::string DIRECTORY = "./Grids/Turbofan/";
//const std::string FILE_NAME = "Turbofan-";



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


const int PS = 3;
const int PT = 5;

const double BBSIZEOFFSET = 0.0;

// OpenMP parameters:

const int NTHREADS = 16;

#endif