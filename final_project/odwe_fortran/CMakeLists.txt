cmake_minimum_required(VERSION 3.16)

project(odwe_fortran LANGUAGES Fortran)

# Use Intel oneAPI ifx
set(CMAKE_Fortran_COMPILER ifx)

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

# Some reasonable flags
target_compile_options(odwe PRIVATE -O3 -warn -traceback)
