cmake_minimum_required(VERSION 3.16)

# 1. Set compiler BEFORE project() if you want to force it, 
#    but better to let CMake detect what is available.
#    (We comment this out to allow auto-detection on ws104)
# set(CMAKE_Fortran_COMPILER ifx)

project(odwe_fortran LANGUAGES Fortran)

# 2. Find OpenMP (works for both Intel and GNU)
find_package(OpenMP REQUIRED)

set(SRC
  main.f90
  grid.f90
  state.f90
  flux.f90
  solver.f90
  bc.f90
  chemistry.f90
  post.f90
)

add_executable(odwe ${SRC})

# 3. Link OpenMP
target_link_libraries(odwe PRIVATE OpenMP::OpenMP_Fortran)

# 4. Apply flags based on the Detected Compiler
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    # Intel Flags (ifx / ifort)
    message(STATUS "Configuring for Intel Compiler...")
    target_compile_options(odwe PRIVATE -O3 -warn -traceback)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # GNU Flags (gfortran)
    message(STATUS "Configuring for GNU gfortran...")
    target_compile_options(odwe PRIVATE -O3 -Wall -fbacktrace -ffree-line-length-none)
else()
    # Fallback
    target_compile_options(odwe PRIVATE -O3)
endif()