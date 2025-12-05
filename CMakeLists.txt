set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
cmake_minimum_required(VERSION 3.20)
project(IFGF LANGUAGES CXX)

if(NOT DEFINED ENV{MKLROOT})
    message(FATAL_ERROR "Environment variable MKLROOT must be set.")
endif()

if(NOT DEFINED ENV{I_MPI_ROOT})
    message(FATAL_ERROR "Environment variable I_MPI_ROOT must be set.")
endif()

if (NOT DEFINED ENV{PETSC_DIR})
    message(FATAL_ERROR "Environment variable PETSC_DIR must be set.")
endif()

if (NOT DEFINED ENV{PETSC_ARCH})
    message(FATAL_ERROR "Environment variable PETSC_ARCH must be set.")
endif()

#-----------------------------------------
# User options
#-----------------------------------------
set(NP 64 CACHE STRING "Number of MPI ranks for the run target")
set(OMP_THREADS 1 CACHE STRING "OMP_NUM_THREADS for the run target")

#-----------------------------------------
# Variables (include dirs, link dirs, libs)
#-----------------------------------------
set(PETSC_DIR $ENV{PETSC_DIR})
set(PETSC_ARCH $ENV{PETSC_ARCH})

set(IFGF_INCLUDE
    ${PETSC_DIR}/include
    ${PETSC_DIR}/${PETSC_ARCH}/include
    /opt/intel/oneapi/mpi/latest/include
    /opt/intel/oneapi/mkl/latest/include
)


set(IFGF_LINK_DIRS
    ${PETSC_DIR}/${PETSC_ARCH}/lib
    /opt/intel/oneapi/mkl/latest/lib/intel64
)

set(IFGF_LIBS
    petsc
    mkl_intel_lp64
    mkl_sequential
    mkl_core
    iomp5
)
set(IFGF_SOURCE_DIR src)

set(IFGF_SOURCE_FILES
     ${IFGF_SOURCE_DIR}/parametrization.cpp
   )

#-----------------------------------------
# Far Field Sphere Test
#-----------------------------------------

add_executable(sphereff ${IFGF_SOURCE_FILES} tests/sphere_test.cpp)

target_compile_options(sphereff PRIVATE -O3 -qopenmp -w)

target_include_directories(sphereff PRIVATE ${IFGF_INCLUDE})

target_link_directories(sphereff PRIVATE ${IFGF_LINK_DIRS})

target_link_libraries(sphereff PRIVATE ${IFGF_LIBS})

# ------------------------------------------
# Eigen function test for complex wave numbers
# ------------------------------------------
# set(COMP_BESSEL_PATH "${CMAKE_SOURCE_DIR}/tests/complex_bessel-master/build/lib64/libcomplex_bessel.so")
# add_executable(complex_efunc ${IFGF_SOURCE_FILES} tests/ComplexEigenfunctionTest.cpp)
# target_compile_options(complex_efunc PRIVATE -O3 -qopenmp -w)
# target_include_directories(complex_efunc PRIVATE ${IFGF_INCLUDE})
# target_link_directories(complex_efunc PRIVATE ${IFGF_LINK_DIRS})
# target_link_libraries(complex_efunc PRIVATE ${IFGF_LIBS})
# target_link_libraries(complex_efunc PRIVATE ${COMP_BESSEL_PATH})

# ----------------------------------------------
# Nacelle Test
# -----------------------------------------------
add_executable(Nacelle ${IFGF_SOURCE_FILES} tests/Nacelle_Test.cpp)
target_compile_options(Nacelle PRIVATE -O3 -qopenmp -w)
target_include_directories(Nacelle PRIVATE ${IFGF_INCLUDE})
target_link_directories(Nacelle PRIVATE ${IFGF_LINK_DIRS})
target_link_libraries(Nacelle PRIVATE ${IFGF_LIBS})

#-----------------------------------------
# MPI (use same wrapper compiler as Makefile)
#-----------------------------------------
set(CMAKE_CXX_COMPILER mpiicpc)

#-----------------------------------------
# Run target (equivalent to Makefile MPI_RUN)
#-----------------------------------------
add_custom_target(run
    COMMAND ${CMAKE_COMMAND} -E env
        OMP_NUM_THREADS=${OMP_THREADS}
        mpirun -np ${NP}
        ./sphereff -ksp_monitor -ksp_type gmres
    DEPENDS main2
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running main2 with mpirun"
)
