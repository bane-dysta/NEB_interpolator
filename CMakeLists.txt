cmake_minimum_required(VERSION 3.10)
project(NEB_interpo LANGUAGES CXX)

include(CheckLanguage)

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

# Optional toolchain/backends
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
        set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
        set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
    endif()
endif()

find_package(BLAS QUIET)
find_package(LAPACK QUIET)

# Python is only used for lightweight regression test scripts
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
else()
    set(PYTHON_EXECUTABLE python3)
endif()

# Intel MKL (optional)
option(USE_MKL "Use Intel MKL for BLAS/LAPACK" ON)

set(MKL_ROOT "")
set(MKL_INCLUDE_DIRS "")
set(MKL_LIBRARIES "")
set(FORTRAN_LIBS "")

set(MKL_FOUND FALSE)
set(HAVE_BLAS_LAPACK FALSE)
if(BLAS_FOUND AND LAPACK_FOUND)
    set(HAVE_BLAS_LAPACK TRUE)
endif()

if(USE_MKL AND DEFINED ENV{MKLROOT})
    set(MKL_ROOT $ENV{MKLROOT})
    if(EXISTS "${MKL_ROOT}/include/mkl.h")
        if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
            set(_MKL_LIB_DIR "${MKL_ROOT}/lib/intel64")
            if(EXISTS "${_MKL_LIB_DIR}/libmkl_core.so" AND
               EXISTS "${_MKL_LIB_DIR}/libmkl_intel_lp64.so" AND
               EXISTS "${_MKL_LIB_DIR}/libmkl_intel_thread.so")
                set(MKL_LIBRARIES
                    ${_MKL_LIB_DIR}/libmkl_intel_lp64.so
                    ${_MKL_LIB_DIR}/libmkl_intel_thread.so
                    ${_MKL_LIB_DIR}/libmkl_core.so
                    -liomp5 -lpthread -lm -ldl
                )
                set(MKL_FOUND TRUE)
            endif()
        elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
            set(_MKL_LIB_DIR "${MKL_ROOT}/lib")
            if(EXISTS "${_MKL_LIB_DIR}/libmkl_core.dylib" AND
               EXISTS "${_MKL_LIB_DIR}/libmkl_intel_lp64.dylib" AND
               EXISTS "${_MKL_LIB_DIR}/libmkl_intel_thread.dylib")
                set(MKL_LIBRARIES
                    ${_MKL_LIB_DIR}/libmkl_intel_lp64.dylib
                    ${_MKL_LIB_DIR}/libmkl_intel_thread.dylib
                    ${_MKL_LIB_DIR}/libmkl_core.dylib
                    -liomp5 -lpthread -lm -ldl
                )
                set(MKL_FOUND TRUE)
            endif()
        endif()
    endif()
endif()

if(MKL_FOUND)
    message(STATUS "Intel MKL enabled at: ${MKL_ROOT}")
    set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include)

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
    if(HAVE_BLAS_LAPACK)
        if(USE_MKL)
            message(STATUS "Intel MKL requested (USE_MKL=ON) but not found; using system BLAS/LAPACK")
        else()
            message(STATUS "Intel MKL disabled (USE_MKL=OFF); using system BLAS/LAPACK")
        endif()
        set(MKL_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        set(MKL_INCLUDE_DIRS "")
        set(FORTRAN_LIBS ${MKL_LIBRARIES})
    else()
        message(STATUS "BLAS/LAPACK not found; calc_rmsd_xyz and allign_lapack will not be built")
        set(MKL_LIBRARIES "")
        set(MKL_INCLUDE_DIRS "")
        set(FORTRAN_LIBS "")
    endif()
endif()

# Create directories
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/examples)

# Set source directory
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# Check if source files exist
if(NOT EXISTS ${SRC_DIR}/neb_driver.h)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/neb_driver.h")
endif()
if(NOT EXISTS ${SRC_DIR}/neb_driver.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/neb_driver.cpp")
endif()
if(NOT EXISTS ${SRC_DIR}/rmsd_align.h)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/rmsd_align.h")
endif()
if(NOT EXISTS ${SRC_DIR}/rmsd_align.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/rmsd_align.cpp")
endif()
if(NOT EXISTS ${SRC_DIR}/calc_rmsd_xyz.f90)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/calc_rmsd_xyz.f90")
endif()
if(NOT EXISTS ${SRC_DIR}/allign_lapack.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/allign_lapack.cpp")
endif()
if(NOT EXISTS ${SRC_DIR}/allign_eigen.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/allign_eigen.cpp")
endif()
if(NOT EXISTS ${SRC_DIR}/xyzgeom.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/xyzgeom.cpp")
endif()
if(NOT EXISTS ${SRC_DIR}/neb_interpolator.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/neb_interpolator.cpp")
endif()

if(NOT EXISTS ${SRC_DIR}/neb_engine.h)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/neb_engine.h")
endif()
if(NOT EXISTS ${SRC_DIR}/neb_engine.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/neb_engine.cpp")
endif()
if(NOT EXISTS ${SRC_DIR}/neb_core.h)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/neb_core.h")
endif()
if(NOT EXISTS ${SRC_DIR}/neb_core.cpp)
    message(FATAL_ERROR "Cannot find ${SRC_DIR}/neb_core.cpp")
endif()

# Copy header file to build directory
configure_file(${SRC_DIR}/neb_driver.h
               ${CMAKE_BINARY_DIR}/include/neb_driver.h
               COPYONLY)
configure_file(${SRC_DIR}/rmsd_align.h
               ${CMAKE_BINARY_DIR}/include/rmsd_align.h
               COPYONLY)

# Copy shared headers to build directory
configure_file(${SRC_DIR}/geometry.h
               ${CMAKE_BINARY_DIR}/include/geometry.h
               COPYONLY)
configure_file(${SRC_DIR}/util.h
               ${CMAKE_BINARY_DIR}/include/util.h
               COPYONLY)
configure_file(${SRC_DIR}/neb_engine.h
               ${CMAKE_BINARY_DIR}/include/neb_engine.h
               COPYONLY)
configure_file(${SRC_DIR}/neb_core.h
               ${CMAKE_BINARY_DIR}/include/neb_core.h
               COPYONLY)

# Include directories
include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${SRC_DIR})  # Also include src directory for headers
if(MKL_INCLUDE_DIRS)
    include_directories(${MKL_INCLUDE_DIRS})
endif()

# ========== Executables ==========

# Shared C++ library (driver + NEB core/engine + interpolation backends)
add_library(neb_lib
    ${SRC_DIR}/neb_driver.cpp
    ${SRC_DIR}/rmsd_align.cpp
    ${SRC_DIR}/neb_engine.cpp
    ${SRC_DIR}/neb_core.cpp
    ${SRC_DIR}/internal_ic.cpp
    ${SRC_DIR}/zmat/covalent_radii.cpp
)

set_target_properties(neb_lib PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

# 1. RMSD alignment backends
set(RMSD_BACKEND_TARGETS allign_eigen)

add_executable(allign_eigen ${SRC_DIR}/allign_eigen.cpp)
target_include_directories(allign_eigen PRIVATE ${CMAKE_SOURCE_DIR}/third_party/eigen-5.0.0)

if(HAVE_BLAS_LAPACK)
    add_executable(allign_lapack ${SRC_DIR}/allign_lapack.cpp)
    target_link_libraries(allign_lapack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    list(APPEND RMSD_BACKEND_TARGETS allign_lapack)
endif()

if(CMAKE_Fortran_COMPILER AND HAVE_BLAS_LAPACK)
    add_executable(calc_rmsd_xyz ${SRC_DIR}/calc_rmsd_xyz.f90)
    target_link_libraries(calc_rmsd_xyz ${FORTRAN_LIBS})
    list(APPEND RMSD_BACKEND_TARGETS calc_rmsd_xyz)
endif()

# 2. C++ XYZ Geometry Analyzer (Enhanced version)
add_executable(xyzgeom ${SRC_DIR}/xyzgeom.cpp)
target_link_libraries(xyzgeom neb_lib)
set_target_properties(xyzgeom PROPERTIES 
    COMPILE_FLAGS "-I${CMAKE_BINARY_DIR}/include -I${SRC_DIR}"
)

# 3. C++ NEB Interpolator standalone executable
add_executable(neb_interpolator ${SRC_DIR}/neb_interpolator.cpp)
target_link_libraries(neb_interpolator neb_lib)
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

# Set output directory for main executables
set(ALL_RUNTIME_TARGETS xyzgeom neb_interpolator ${RMSD_BACKEND_TARGETS})
set_target_properties(${ALL_RUNTIME_TARGETS}
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
if(TARGET calc_rmsd_xyz)
add_custom_target(test_rmsd
    COMMAND ${CMAKE_BINARY_DIR}/bin/calc_rmsd_xyz 
        ${CMAKE_BINARY_DIR}/examples/water.xyz 
        ${CMAKE_BINARY_DIR}/examples/water_rotated.xyz
    DEPENDS calc_rmsd_xyz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing RMSD alignment with water molecules"
)
endif()

# Test XYZ geometry analyzer
add_custom_target(test_xyzgeom
    COMMAND echo "11" | ${CMAKE_BINARY_DIR}/bin/xyzgeom 
        ${CMAKE_BINARY_DIR}/examples/water.xyz --print
    DEPENDS xyzgeom
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing XYZ geometry analyzer"
)

# Test LIC interpolation (Cartesian linear)
add_custom_target(test_lic
    COMMAND ${CMAKE_BINARY_DIR}/bin/neb_interpolator 
        -n 3 -m lic
        ${CMAKE_BINARY_DIR}/examples/water.xyz 
        ${CMAKE_BINARY_DIR}/examples/water_rotated.xyz
    DEPENDS neb_interpolator
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing LIC interpolation"
)

# Test LIIC interpolation (water)
add_custom_target(test_liic
    COMMAND ${CMAKE_BINARY_DIR}/bin/neb_interpolator
        --no-align -n 3 -m liic
        ${CMAKE_BINARY_DIR}/examples/water.xyz
        ${CMAKE_BINARY_DIR}/examples/water_rotated.xyz
    DEPENDS neb_interpolator
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing LIIC interpolation with water molecules"
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
if(TARGET calc_rmsd_xyz)
add_custom_target(test_rmsd_range
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_rmsd_range.py
        ${CMAKE_BINARY_DIR}/bin/calc_rmsd_xyz
        ${CMAKE_BINARY_DIR}/examples/range_ref.xyz
        ${CMAKE_BINARY_DIR}/examples/range_mobile.xyz
    DEPENDS calc_rmsd_xyz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing calc_rmsd_xyz with explicit atom ranges"
)
endif()

# Test LIIC->DM fallback path (forced LIIC failure via --ev-thresh)
add_custom_target(test_dm_fallback
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_dm_fallback.py
        ${CMAKE_BINARY_DIR}/bin/neb_interpolator
        ${CMAKE_BINARY_DIR}/examples/bondshift_init.xyz
        ${CMAKE_BINARY_DIR}/examples/bondshift_final.xyz
    DEPENDS neb_interpolator
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing DM fallback when LIIC fails"
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


# Minimal C++ unit test for NEB core (assembleNEBForces)
add_executable(neb_core_test
    ${CMAKE_SOURCE_DIR}/tests/test_neb_core.cpp
    ${SRC_DIR}/neb_core.cpp
)
set_target_properties(neb_core_test PROPERTIES
    COMPILE_FLAGS "-I${CMAKE_BINARY_DIR}/include -I${SRC_DIR}"
)
set_target_properties(neb_core_test
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

add_custom_target(test_neb_core
    COMMAND ${CMAKE_BINARY_DIR}/bin/neb_core_test
    DEPENDS neb_core_test
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Testing NEB core (assembleNEBForces)"
)

# Run all tests
set(TEST_ALL_DEPENDS test_xyzgeom test_lic test_liic test_dm test_dm_fallback test_external_engine test_neb_core)
if(TARGET calc_rmsd_xyz)
    list(APPEND TEST_ALL_DEPENDS test_rmsd test_rmsd_range)
endif()

add_custom_target(test_all
    DEPENDS ${TEST_ALL_DEPENDS}
    COMMENT "Running all tests"
)

# ========== Installation ==========

install(TARGETS xyzgeom neb_interpolator ${RMSD_BACKEND_TARGETS}
    RUNTIME DESTINATION bin
)

install(FILES ${SRC_DIR}/neb_driver.h
    DESTINATION include
)

install(FILES ${SRC_DIR}/rmsd_align.h
    DESTINATION include
)

install(FILES ${SRC_DIR}/neb_engine.h
    DESTINATION include
)

install(FILES ${SRC_DIR}/neb_core.h
    DESTINATION include
)

install(FILES ${SRC_DIR}/internal_ic.h
    DESTINATION include
)

install(FILES ${SRC_DIR}/geometry.h
    DESTINATION include
)

install(FILES ${SRC_DIR}/util.h
    DESTINATION include
)

install(FILES ${SRC_DIR}/zmat/covalent_radii.h
    DESTINATION include/zmat
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
if(CMAKE_Fortran_COMPILER)
    message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
    message(STATUS "Fortran flags: ${CMAKE_Fortran_FLAGS}")
else()
    message(STATUS "Fortran compiler: not found")
endif()
if(MKL_FOUND)
    message(STATUS "BLAS/LAPACK: Intel MKL at ${MKL_ROOT}")
elseif(HAVE_BLAS_LAPACK)
    message(STATUS "BLAS/LAPACK: System libraries")
else()
    message(STATUS "BLAS/LAPACK: not found")
endif()
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Executables will be built in: ${CMAKE_BINARY_DIR}/bin")
message(STATUS "Enabled RMSD backends: ${RMSD_BACKEND_TARGETS}")
message(STATUS "===========================")

# ========== Usage Target (替换原来的help) ==========

add_custom_target(usage
    COMMAND echo "Available targets:"
    COMMAND echo "  all            - Build all executables"
    COMMAND echo "  xyzgeom        - Build enhanced XYZ geometry analyzer"
    COMMAND echo "  calc_rmsd_xyz  - Build Fortran RMSD calculator (only when Fortran+BLAS/LAPACK are available)"
    COMMAND echo "  neb_interpolator - Build standalone NEB interpolator"
    COMMAND echo "  test_all       - Run all tests"
    COMMAND echo "  test_neb_core  - Run NEB core unit tests"
    COMMAND echo "  test_rmsd      - Test RMSD alignment"
    COMMAND echo "  test_xyzgeom   - Test geometry analyzer"
    COMMAND echo "  test_lic       - Test LIC interpolation"
    COMMAND echo "  test_liic      - Test LIIC interpolation"
    COMMAND echo "  test_dm        - Test DM interpolation"
    COMMAND echo "  test_dm_fallback - Test LIIC -> DM fallback"
    COMMAND echo "  test_external_engine - Test external engine interface"
    COMMAND echo "  install        - Install executables and files"
    COMMAND echo "  clean          - Clean build files"
    COMMAND echo ""
    COMMAND echo "Usage examples after build:"
    COMMAND echo "  cd bin"
    COMMAND echo "  ./xyzgeom molecule.xyz              # Interactive mode"
    COMMAND echo "  ./xyzgeom molecule.xyz --print      # Print structure"
    COMMAND echo "  ./calc_rmsd_xyz ref.xyz mobile.xyz  # Align structures (if built)"
    COMMAND echo "  ./neb_interpolator initial.xyz final.xyz  # NEB interpolation"
    VERBATIM
)

# ========== Clean Target ==========

add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/include
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/examples
    COMMAND ${CMAKE_COMMAND} -DBINARY_DIR=${CMAKE_BINARY_DIR} -P ${CMAKE_SOURCE_DIR}/cmake/clean_xyz.cmake
    COMMENT "Cleaning all generated files"
    VERBATIM
)