cmake_minimum_required(VERSION 3.10)
project(MolecularGeometryTools LANGUAGES CXX Fortran)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
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

# Find required packages
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Python is only used for lightweight regression test scripts
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
else()
    set(PYTHON_EXECUTABLE python3)
endif()

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
add_executable(neb_interpolator 
    ${SRC_DIR}/neb_interpolator.cpp
    ${SRC_DIR}/internal_ic.cpp
    ${SRC_DIR}/zmat/covalent_radii.cpp
)
set_target_properties(neb_interpolator PROPERTIES 
    COMPILE_FLAGS "-I${CMAKE_BINARY_DIR}/include -I${SRC_DIR}"
)

# GCC < 9 may require linking stdc++fs explicitly for std::filesystem
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
        target_link_libraries(xyzgeom stdc++fs)
        target_link_libraries(neb_interpolator stdc++fs)
    endif()
endif()

# Set output directory for all executables
set_target_properties(calc_rmsd_xyz xyzgeom neb_interpolator
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# ========== Example Files ==========

# Create example XYZ files for testing
file(WRITE ${CMAKE_BINARY_DIR}/examples/water.xyz 
"3
Water molecule
O    0.000000    0.000000    0.000000
H    0.758602    0.000000    0.504284
H    0.758602    0.000000   -0.504284
")

file(WRITE ${CMAKE_BINARY_DIR}/examples/water_rotated.xyz 
"3
Water molecule (rotated)
O    0.000000    0.000000    0.000000
H    0.504284    0.758602    0.000000
H   -0.504284    0.758602    0.000000
")

file(WRITE ${CMAKE_BINARY_DIR}/examples/ammonia.xyz 
"4
Ammonia molecule
N    0.000000    0.000000    0.000000
H    0.000000    0.000000    1.012316
H    0.969771    0.000000   -0.290275
H   -0.390071   -0.887477   -0.290275
")

# Extra small molecules for regression tests
file(WRITE ${CMAKE_BINARY_DIR}/examples/range_ref.xyz
"3
Range test reference
H    0.000000    0.000000    0.000000
H    1.000000    0.000000    0.000000
H    0.000000    1.000000    0.000000
")

file(WRITE ${CMAKE_BINARY_DIR}/examples/range_mobile.xyz
"3
Range test mobile (atom 3 shifted)
H    0.000000    0.000000    0.000000
H    1.000000    0.000000    0.000000
H    0.000000    2.000000    0.000000
")

# Simple 3-atom bond-shift style example (H3 rearrangement-like)
file(WRITE ${CMAKE_BINARY_DIR}/examples/bondshift_init.xyz
"3
Bond shift init
H    0.000000    0.000000    0.000000
H    0.740000    0.000000    0.000000
H    3.000000    0.000000    0.000000
")

file(WRITE ${CMAKE_BINARY_DIR}/examples/bondshift_final.xyz
"3
Bond shift final
H    0.000000    0.000000    0.000000
H    2.260000    0.000000    0.000000
H    3.000000    0.000000    0.000000
")

# ========== Testing Targets ==========

# Test RMSD alignment
add_custom_target(test_rmsd
    COMMAND ${CMAKE_BINARY_DIR}/bin/calc_rmsd_xyz 
        ${CMAKE_BINARY_DIR}/examples/water.xyz 
        ${CMAKE_BINARY_DIR}/examples/water_rotated.xyz
    DEPENDS calc_rmsd_xyz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing RMSD alignment with water molecules"
)

# Test XYZ geometry analyzer
add_custom_target(test_xyzgeom
    COMMAND echo "11" | ${CMAKE_BINARY_DIR}/bin/xyzgeom 
        ${CMAKE_BINARY_DIR}/examples/water.xyz --print
    DEPENDS xyzgeom
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing XYZ geometry analyzer"
)

# Test NEB interpolation
add_custom_target(test_neb
    COMMAND ${CMAKE_BINARY_DIR}/bin/neb_interpolator 
        -n 3 -m liic
        ${CMAKE_BINARY_DIR}/examples/water.xyz 
        ${CMAKE_BINARY_DIR}/examples/water_rotated.xyz
    DEPENDS neb_interpolator calc_rmsd_xyz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing NEB interpolation"
)

# Test IIC interpolation (water)
add_custom_target(test_iic
    COMMAND ${CMAKE_BINARY_DIR}/bin/neb_interpolator
        --no-align -n 3 -m iic
        ${CMAKE_BINARY_DIR}/examples/water.xyz
        ${CMAKE_BINARY_DIR}/examples/water_rotated.xyz
    DEPENDS neb_interpolator
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing IIC interpolation with water molecules"
)

# Test DM interpolation (water)
add_custom_target(test_dm
    COMMAND ${CMAKE_BINARY_DIR}/bin/neb_interpolator
        --no-align -n 3 -m dm
        ${CMAKE_BINARY_DIR}/examples/water.xyz
        ${CMAKE_BINARY_DIR}/examples/water_rotated.xyz
    DEPENDS neb_interpolator
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing DM interpolation with water molecules"
)

# Test calc_rmsd_xyz range parsing (regression for a variable-reuse bug)
add_custom_target(test_rmsd_range
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_rmsd_range.py
        ${CMAKE_BINARY_DIR}/bin/calc_rmsd_xyz
        ${CMAKE_BINARY_DIR}/examples/range_ref.xyz
        ${CMAKE_BINARY_DIR}/examples/range_mobile.xyz
    DEPENDS calc_rmsd_xyz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing calc_rmsd_xyz with explicit atom ranges"
)

# Test IIC->DM fallback path (forced IIC failure via --ev-thresh)
add_custom_target(test_dm_fallback
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_dm_fallback.py
        ${CMAKE_BINARY_DIR}/bin/neb_interpolator
        ${CMAKE_BINARY_DIR}/examples/bondshift_init.xyz
        ${CMAKE_BINARY_DIR}/examples/bondshift_final.xyz
    DEPENDS neb_interpolator
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing DM fallback when IIC fails"
)

# Test external engine interface with a mock engine that outputs zero gradients
add_custom_target(test_external_engine
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_external_engine.py
        ${CMAKE_BINARY_DIR}/bin/neb_interpolator
        ${CMAKE_SOURCE_DIR}/tests/mock_engine_constant.py
        ${CMAKE_BINARY_DIR}/examples/water.xyz
        ${CMAKE_BINARY_DIR}/examples/water_rotated.xyz
    DEPENDS neb_interpolator
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing external engine interface"
)


# Run all tests
add_custom_target(test_all
    DEPENDS test_rmsd test_rmsd_range test_xyzgeom test_neb test_iic test_dm test_dm_fallback test_external_engine
    COMMENT "Running all tests"
)

# ========== Installation ==========

install(TARGETS calc_rmsd_xyz xyzgeom neb_interpolator
    RUNTIME DESTINATION bin
)

install(FILES ${SRC_DIR}/neb_interpolator.h
    DESTINATION include
)

install(DIRECTORY ${CMAKE_BINARY_DIR}/examples/
    DESTINATION share/${PROJECT_NAME}/examples
    FILES_MATCHING PATTERN "*.xyz"
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

# ========== Usage Target (替换原来的help) ==========

add_custom_target(usage
    COMMAND echo "Available targets:"
    COMMAND echo "  all            - Build all executables"
    COMMAND echo "  xyzgeom        - Build enhanced XYZ geometry analyzer"
    COMMAND echo "  calc_rmsd_xyz  - Build Fortran RMSD calculator"
    COMMAND echo "  neb_interpolator - Build standalone NEB interpolator"
    COMMAND echo "  test_all       - Run all tests"
    COMMAND echo "  test_rmsd      - Test RMSD alignment"
    COMMAND echo "  test_xyzgeom   - Test geometry analyzer"
    COMMAND echo "  test_neb       - Test NEB interpolation"
    COMMAND echo "  install        - Install executables and files"
    COMMAND echo "  clean          - Clean build files"
    COMMAND echo ""
    COMMAND echo "Usage examples after build:"
    COMMAND echo "  cd bin"
    COMMAND echo "  ./xyzgeom molecule.xyz              # Interactive mode"
    COMMAND echo "  ./xyzgeom molecule.xyz --print      # Print structure"
    COMMAND echo "  ./calc_rmsd_xyz ref.xyz mobile.xyz  # Align structures"
    COMMAND echo "  ./neb_interpolator initial.xyz final.xyz  # NEB interpolation"
    VERBATIM
)

# ========== Clean Target ==========

add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/include
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/examples
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/*.xyz 2>/dev/null || true
    COMMENT "Cleaning all generated files"
    VERBATIM
)