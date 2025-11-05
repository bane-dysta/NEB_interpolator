cmake_minimum_required(VERSION 3.10)
project(MolecularGeometryTools LANGUAGES CXX Fortran)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler-specific flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
endif()

# Option to force using system BLAS instead of MKL
option(USE_SYSTEM_BLAS "Force using system BLAS/LAPACK instead of MKL" OFF)

# Find required packages
# If USE_SYSTEM_BLAS is ON, limit search to avoid MKL
if(USE_SYSTEM_BLAS)
    # Save MKLROOT temporarily
    set(MKLROOT_SAVED $ENV{MKLROOT})
    # Set to empty string to prevent find_package from using MKL
    set(ENV{MKLROOT} "")
    # Set BLA_VENDOR to Generic to prefer standard BLAS and avoid MKL
    set(BLA_VENDOR "Generic" CACHE STRING "BLAS vendor to search for" FORCE)
    # Disable MKL search by setting MKL_ROOT to empty
    set(MKL_ROOT "" CACHE PATH "Intel MKL root directory" FORCE)
    message(STATUS "USE_SYSTEM_BLAS=ON: Limiting BLAS/LAPACK search to system libraries (skipping MKL)")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Restore MKLROOT if it was set
if(USE_SYSTEM_BLAS AND MKLROOT_SAVED)
    set(ENV{MKLROOT} ${MKLROOT_SAVED})
endif()

# Check for Intel MKL (preferred)
set(USE_MKL FALSE)
if(NOT USE_SYSTEM_BLAS AND DEFINED ENV{MKLROOT})
    message(STATUS "Intel MKL found at: $ENV{MKLROOT}")
    set(MKL_ROOT $ENV{MKLROOT})
    set(USE_MKL TRUE)
    
    # Intel MKL include directory
    set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include)
    
    # Intel MKL libraries
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(MKL_LIBRARIES
            ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
            ${MKL_ROOT}/lib/intel64/libmkl_intel_thread.so
            ${MKL_ROOT}/lib/intel64/libmkl_core.so
            -liomp5 -lpthread -lm -ldl
        )
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        set(MKL_LIBRARIES
            ${MKL_ROOT}/lib/libmkl_intel_lp64.dylib
            ${MKL_ROOT}/lib/libmkl_intel_thread.dylib
            ${MKL_ROOT}/lib/libmkl_core.dylib
            -liomp5 -lpthread -lm -ldl
        )
    endif()
    
    # Use MKL compiler flags if using Intel compilers
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl")
        set(FORTRAN_LIBS "")  # MKL linked via compiler flag
    else()
        set(FORTRAN_LIBS ${MKL_LIBRARIES})
    endif()
    
    if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl")
    endif()
    
else()
    if(USE_SYSTEM_BLAS)
        message(STATUS "USE_SYSTEM_BLAS=ON: Using system BLAS/LAPACK")
    else()
        message(STATUS "Intel MKL not found, using system BLAS/LAPACK")
    endif()
    set(MKL_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    set(MKL_INCLUDE_DIRS "")
    set(FORTRAN_LIBS ${MKL_LIBRARIES})
endif()

# Create directories
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

# Set source directory
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# Check if source files exist
if(NOT EXISTS ${SRC_DIR}/neb_interpolator.h)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/neb_interpolator.h")
endif()
if(NOT EXISTS ${SRC_DIR}/calc_rmsd_xyz.f90)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/calc_rmsd_xyz.f90")
endif()
if(NOT EXISTS ${SRC_DIR}/xyzgeom.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/xyzgeom.cpp")
endif()
if(NOT EXISTS ${SRC_DIR}/neb_interpolator.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/neb_interpolator.cpp")
endif()

# Copy header file to build directory
configure_file(${SRC_DIR}/neb_interpolator.h 
               ${CMAKE_BINARY_DIR}/include/neb_interpolator.h 
               COPYONLY)

# Include directories
include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${SRC_DIR})  # Also include src directory for headers
if(MKL_INCLUDE_DIRS)
    include_directories(${MKL_INCLUDE_DIRS})
endif()

# ========== Executables ==========

# 1. Fortran RMSD calculator executable
add_executable(calc_rmsd_xyz ${SRC_DIR}/calc_rmsd_xyz.f90)
target_link_libraries(calc_rmsd_xyz ${FORTRAN_LIBS})

# 2. C++ XYZ Geometry Analyzer (Enhanced version)
add_executable(xyzgeom ${SRC_DIR}/xyzgeom.cpp)
set_target_properties(xyzgeom PROPERTIES 
    COMPILE_FLAGS "-I${CMAKE_BINARY_DIR}/include -I${SRC_DIR}"
)

# 3. C++ NEB Interpolator standalone executable
add_executable(neb_interpolator ${SRC_DIR}/neb_interpolator.cpp)
set_target_properties(neb_interpolator PROPERTIES 
    COMPILE_FLAGS "-I${CMAKE_BINARY_DIR}/include -I${SRC_DIR}"
)

# Set output directory for all executables
set_target_properties(calc_rmsd_xyz xyzgeom neb_interpolator
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# ========== Installation ==========

install(TARGETS calc_rmsd_xyz xyzgeom neb_interpolator
    RUNTIME DESTINATION bin
)

install(FILES ${SRC_DIR}/neb_interpolator.h
    DESTINATION include
)

# ========== Build Information ==========

message(STATUS "")
message(STATUS "=== ${PROJECT_NAME} Build Configuration ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Source directory: ${SRC_DIR}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran flags: ${CMAKE_Fortran_FLAGS}")
if(USE_MKL)
    message(STATUS "BLAS/LAPACK: Intel MKL at ${MKL_ROOT}")
else()
    message(STATUS "BLAS/LAPACK: System libraries")
endif()
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
message(STATUS "Source files:")
message(STATUS "  ${SRC_DIR}/xyzgeom.cpp")
message(STATUS "  ${SRC_DIR}/calc_rmsd_xyz.f90")
message(STATUS "  ${SRC_DIR}/neb_interpolator.cpp")
message(STATUS "  ${SRC_DIR}/neb_interpolator.h")
message(STATUS "")
message(STATUS "Executables will be built in: ${CMAKE_BINARY_DIR}/bin")
message(STATUS "===========================")