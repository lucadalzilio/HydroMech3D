# Compiler and flags
# PATH for required libraries
ARMADILLO_DIR = /opt/homebrew/opt/armadillo
ifneq ($(MKLROOT),)
    MKLROOT :=/opt/intel/oneapi/mkl/latest
endif
OPENBLAS_DIR=/opt/homebrew/opt/openblas
OPENMPI_DIR=/opt/homebrew/opt/openmpi
LIBOMP_DIR=/opt/homebrew/opt/libomp
#TBB_DIR=/home/zhenhuan001/Apps/oneTBB
CXX = mpicxx
USE_MKL=0
ifeq ($(USE_MKL), 1)
CXXFLAGS = -g -std=c++14 -Wall -Wextra -pedantic -O0 -m64 -DIL_MKL -DIL_BLAS

LDFLAGS  = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -L$(LIBOMP_DIR)/lib     -lomp -ldl -lpthread -lm -L$(ARMADILLO_DIR)/lib -larmadillo

INCLUDES := -I./JSON -I$(MKLROOT)/include -I./ -I$(ARMADILLO_DIR)/include -I    $(LIBOMP_DIR)/include
else 
CXXFLAGS = -g -std=c++14 -Wall -Wextra -pedantic -O0  

LDFLAGS  = -L$(OPENMPI_DIR)/lib -fopenmp -L$(OPENBLAS_DIR)/lib -lopenblas -L$(ARMADILLO_DIR)/lib -larmadillo

INCLUDES := -I./JSON -I$(OPENBLAS_DIR)/include -I./ -I$(ARMADILLO_DIR)/include -I$(OPENMPI_DIR)/include
endif
#CXXFLAGS = -g -std=c++17 -Wall -Wextra -pedantic -O0 -m64 -DIL_MKL -DIL_BLAS 
#LDFLAGS := -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -L$(LIBOMP_DIR)/lib -lomp -ldl -lpthread -lm -L$(ARMADILLO_DIR)/lib -larmadillo 
#LDFLAGS += $(shell hlibpro-3.1/bin/hpro-config --lflags) 
# Include directories
#INCLUDES += $(shell hlibpro-3.1/bin/hpro-config --cflags)
# Source files
SRC_DIR = src
SRC_FILES = $(SRC_DIR)/loadProgramArguments.cpp \
            $(SRC_DIR)/EQsolver.cpp \
            $(SRC_DIR)/ImportJsonInputData.cpp \
            $(SRC_DIR)/Mesh.cpp \
            $(SRC_DIR)/FullSpaceElasticity.cpp \
            $(SRC_DIR)/StressKernelsP0/StressKernelsDxP0.cpp \
            $(SRC_DIR)/StressKernelsP0/StressKernelsDyP0.cpp \
            $(SRC_DIR)/StressKernelsP0/StressKernelsDzP0.cpp \
            $(SRC_DIR)/Solution.cpp \
            $(SRC_DIR)/RK45.cpp \
            $(SRC_DIR)/Utils.cpp \
            $(SRC_DIR)/AssembleElasticityMatrix.cpp \
            main.cpp

# Object files
OBJ_FILES = $(SRC_FILES:.cpp=.o)

# Target executable
TARGET = 3DEQSim

# Rules
all: $(TARGET)

$(TARGET): $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) $(LDFLAGS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJ_FILES) $(TARGET)

.PHONY: all clean
