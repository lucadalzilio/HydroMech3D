cmake_minimum_required(VERSION 3.19)
project(3DEQSim)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++14")

################################################################################
# choose compiler: clang, gcc, ???
################################################################################

#set(CMAKE_C_COMPILER $ENV{HOME}/mpi/bin/mpicc)
#set(CMAKE_CXX_COMPILER $ENV{HOME}/mpi/bin/mpicxx)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)
MESSAGE("ID C compiler = ${CMAKE_C_COMPILER_ID}")
MESSAGE("ID C++ compiler = ${CMAKE_CXX_COMPILER_ID}")
MESSAGE("Build type = ${CMAKE_BUILD_TYPE}")

#####
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIR})
####

################################################################################
# set flags for compiler in Debug mode
################################################################################
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -O0")
        MESSAGE("Flags for building in Debug mode = ${CMAKE_CXX_FLAGS}")
    endif()
ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")


################################################################################
# math framework, i.e. Intel MKL
################################################################################

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
	    link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../../compiler/latest/mac/compiler/lib)
endif()
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
           set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
           set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
           link_directories($ENV{MKLROOT}/lib/intel64_lin)	   
endif()
################################################################################
#  FOR 3DEQSim
################################################################################
# Include InsideLoop library
include_directories(il)
#include_directories(hlibpro-3.1/include)

# Include JSON library
include_directories(JSON)
# Include Openblas library
include_directories(${OpenBLAS_INCLUDE_DIR})
#add_library(hpro STATIC IMPORTED)
#set_target_properties(hpro PROPERTIES IMPORTED_LOCATION /home/zhenhuan001/HydroMech3D/hlibpro-3.1/lib/libhpro.so)
# Set directory containing all the source files
set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

# Set the source files inside the directory previously set that will be compiled
set(SOURCE_FILES
        src/loadProgramArguments.cpp
        src/loadProgramArguments.h
        src/EQsolver.cpp
        src/EQsolver.h
        src/ImportJsonInputData.cpp
        src/ImportJsonInputData.h
        src/FluidProperties.h
        src/SolidMatrixProperties.h
        src/Mesh.cpp
        src/Mesh.h
        src/FaultInSituStress.h
        src/ElementData.h
        src/FullSpaceElasticity.cpp
        src/FullSpaceElasticity.h
        src/StressKernelsP0/StressKernelsDxP0.cpp
        src/StressKernelsP0/StressKernelsDxP0.h
        src/StressKernelsP0/StressKernelsDyP0.cpp
        src/StressKernelsP0/StressKernelsDyP0.h
        src/StressKernelsP0/StressKernelsDzP0.cpp
        src/StressKernelsP0/StressKernelsDzP0.h
        src/Solution.h
        src/Solution.cpp
        src/RK45.cpp
        src/RK45.h
        src/FrictionProperties.h
        src/PermeabilityProperties.h
        src/RightHandSidesODEs.h
        src/FaultProperties.h
        src/Injection.h
        src/SolverParameters.h
        src/Utils.cpp
        src/Utils.h
        src/AssembleElasticityMatrix.h
        src/AssembleElasticityMatrix.cpp)

# Create an executable called for 3DEQSim from the main.cpp file
add_executable(3DEQSim ${SOURCE_FILES} main.cpp)

# Link the MKL library to the executable 3DEQSim
# ${CMAKE_MKL_LIBRARIES} should specify the path in which the linker should
# search for libraries when linking 3DEQSim
#link_directories(hlibpro-3.1/lib)
#target_link_libraries(3DEQSim hpro ${Boost_LIBRARIES} ${TBB_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(3DEQSim ${CMAKE_MKL_LIBRARIES} ${ARMADILLO_LIBRARIES})
# Include directory to be compiled
target_include_directories(3DEQSim PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
