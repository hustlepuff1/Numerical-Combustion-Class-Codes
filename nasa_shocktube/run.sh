#!/bin/bash

# Stop on errors
set -e

# --- 1. Setup Directories ---
echo "-----------------------------------"
echo "   Setting up Project Directories  "
echo "-----------------------------------"

mkdir -p results_fortran
mkdir -p results_python
mkdir -p build

# --- 2. Build ---
echo ">>> Compiling..."
cd build
FC=ifx cmake ..
make
cd ..

# --- 3. Run Solver ---
echo " "
echo ">>> Running Solver (OpenMP)..."
export OMP_NUM_THREADS=8
./nasa_solver

# --- 4. Organize Output ---
echo ">>> Organizing Data..."
# Move all frame files and history
# 2>/dev/null suppresses error if no files found (e.g. if solver crashed early)
mv flow_*.dat results_fortran/ 2>/dev/null || true
mv history.dat results_fortran/ 2>/dev/null || true

# --- 5. Run Visualization ---
echo " "
echo ">>> Plotting Results..."

# Plot Fortran results(We can select which frame to plot by changing the number)
python3 plot_results.py 9

echo "DONE."