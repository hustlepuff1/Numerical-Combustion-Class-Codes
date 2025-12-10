cmake_minimum_required(VERSION 3.10)
project(nasa_solver LANGUAGES Fortran)

# Optimization Flags + OpenMP
# Check if using GFortran or Intel (ifx)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -g -fopenmp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -g -qopenmp")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(nasa_solver 
    constants_mod.f90
    thermo_mod.f90
    reaction_mod.f90
    flux_mod.f90
    solver_mod.f90
    main.f90
)