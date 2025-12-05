# Define the C++ compiler
CXX = mpiicpx

# Define compiler flags
CXXFLAGS = -O3 -qopenmp -I$(MKLROOT)/include -I$(I_MPI_ROOT)/include -I/opt/Eigen/ 

# Define libraries
LDFLAGS = -L$(MKLROOT)/lib -L$(I_MPI_ROOT)/lib -lfftw3 -lfftw3_threads -lm -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# Define source files
SRCS = LinearOperator.cpp main.cpp

# Define the output executable name
TARGET = main

# Default number of MPI processes for running
NP ?= 3

.PHONY: all clean run

all: run

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

run: $(TARGET)
	mpirun -np $(NP) ./$(TARGET)

clean:
	rm -f $(TARGET) *.o