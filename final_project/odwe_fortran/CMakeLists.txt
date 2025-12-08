cmake_minimum_required(VERSION 3.16)

project(odwe_fortran LANGUAGES Fortran)

# Use Intel oneAPI ifx
set(CMAKE_Fortran_COMPILER ifx)

# --- NEW: Enable OpenMP ---
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

# Link OpenMP to the target
target_link_libraries(odwe PRIVATE OpenMP::OpenMP_Fortran)

# Flags: Add -O3 (Optimization) and Traceback
# OpenMP flags are handled automatically by target_link_libraries usually,
# but we add them explicitly to compile options just in case.
target_compile_options(odwe PRIVATE ${OpenMP_Fortran_FLAGS} -O3 -warn -traceback)