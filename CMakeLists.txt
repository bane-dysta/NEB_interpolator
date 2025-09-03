cmake_minimum_required(VERSION 3.10)
project(NEB_RMSD LANGUAGES CXX Fortran)

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
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
endif()

# Find required packages
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Check for Intel MKL (preferred)
set(USE_MKL FALSE)
if(DEFINED ENV{MKLROOT})
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
    message(STATUS "Intel MKL not found, using system BLAS/LAPACK")
    set(MKL_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    set(MKL_INCLUDE_DIRS "")
    set(FORTRAN_LIBS ${MKL_LIBRARIES})
endif()

# Create bin directory
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Fortran RMSD calculator executable
add_executable(calc_rmsd_xyz src/calc_rmsd_xyz.f90)
target_link_libraries(calc_rmsd_xyz ${FORTRAN_LIBS})
if(MKL_INCLUDE_DIRS)
    target_include_directories(calc_rmsd_xyz PRIVATE ${MKL_INCLUDE_DIRS})
endif()

# C++ NEB interpolator executable
add_executable(neb_interpolator src/neb_interpolator.cpp)

# Set output directory for both executables
set_target_properties(calc_rmsd_xyz neb_interpolator
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Create a custom target to copy examples
add_custom_target(copy_examples ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/examples
        ${CMAKE_BINARY_DIR}/examples
    COMMENT "Copying example files to build directory"
)

# Installation
install(TARGETS calc_rmsd_xyz neb_interpolator
    RUNTIME DESTINATION bin
)

install(DIRECTORY examples/
    DESTINATION share/neb_rmsd/examples
    FILES_MATCHING PATTERN "*.xyz"
)

# Add a custom target for running tests
add_custom_target(test_run
    COMMAND ${CMAKE_BINARY_DIR}/bin/neb_interpolator 
        ${CMAKE_BINARY_DIR}/examples/molecule1.xyz 
        ${CMAKE_BINARY_DIR}/examples/molecule2.xyz
    DEPENDS neb_interpolator calc_rmsd_xyz copy_examples
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Running test with example molecules"
)

# Print build information
message(STATUS "=== Build Configuration ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
if(USE_MKL)
    message(STATUS "BLAS/LAPACK: Intel MKL at ${MKL_ROOT}")
else()
    message(STATUS "BLAS/LAPACK: System libraries")
endif()
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "===========================")

# Add help target
add_custom_target(help_build
    COMMAND echo "Available targets:"
    COMMAND echo "  all           - Build all executables"
    COMMAND echo "  calc_rmsd_xyz - Build Fortran RMSD calculator"
    COMMAND echo "  neb_interpolator - Build C++ NEB interpolator"
    COMMAND echo "  copy_examples - Copy example files"
    COMMAND echo "  test_run      - Run test with example molecules"
    COMMAND echo "  install       - Install executables"
    COMMAND echo ""
    COMMAND echo "Usage after build:"
    COMMAND echo "  cd bin && ./neb_interpolator ../examples/molecule1.xyz ../examples/molecule2.xyz"
)