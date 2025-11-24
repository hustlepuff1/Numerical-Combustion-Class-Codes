cmake_minimum_required(VERSION 3.15)
project(CombustionHW6 LANGUAGES Fortran)

# --- Compiler Configuration ---
# Check if using Intel LLVM (ifx)
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    message(STATUS "Configuring for Intel ifx...")
    # -O3: Aggressive optimization
    # -g: Include debug info (useful even in release)
    # -traceback: precise line numbers if code crashes (CRITICAL for debugging)
    # -qmkl: Optional, if you ever decide to use MKL later for linear algebra
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -traceback")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -check all")
else()
    message(STATUS "Using generic Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# --- Source Files ---
# Order technically doesn't matter for CMake (it scans dependencies),
# but listing logical order is good practice.
set(SOURCES
    mod_global.f90
    mod_linalg.f90
    mod_chem.f90
    main.f90
)

# --- Targets ---
add_executable(solver ${SOURCES})

# Clean up module files (.mod) nicely
set_target_properties(solver PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)