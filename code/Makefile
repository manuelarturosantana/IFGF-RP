# Compiler
CXX = mpiicpc

# Compiler Flags
CXXFLAGS = -O3 -qopenmp -w

# Include Paths
INCLUDES = \
    -I${PETSC_DIR}/include \
    -I${PETSC_DIR}/${PETSC_ARCH}/include \
	-I/opt/intel/oneapi/mpi/latest/include \
	-I/opt/intel/oneapi/mkl/latest/include

# Library Paths and Linking
LDFLAGS = \
    -L${PETSC_DIR}/${PETSC_ARCH}/lib \
	-L/opt/intel/oneapi/mkl/latest/lib/intel64 \
    -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
    -lpetsc

# Source Files
SOURCES = main2.cpp

# Executable
TARGET = main2

# MPI Run Settings
NP = 8
MPI_RUN = OMP_NUM_THREADS=16 mpirun -np $(NP) ./$(TARGET) -ksp_monitor -ksp_type gmres

# Default target
all: $(TARGET)

# Build rule
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCES) -o $(TARGET) $(LDFLAGS)

# Run the executable with mpirun
run: $(TARGET)
	$(MPI_RUN)

# Clean build artifacts
clean:
	rm -f $(TARGET)

.PHONY: all run clean