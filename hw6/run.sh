#!/bin/bash

# Stop on error
set -e

echo "========================================"
echo "   Combustion Homework Pipeline"
echo "========================================"

# 1. Setup Build Directory
if [ ! -d "build" ]; then
    echo "[Setup] Creating build directory..."
    mkdir build
fi

cd build

# 2. Configure CMake (Force ifx compiler)
# Note: Ensure you have sourced intel environment vars (source setvars.sh) before running this script
echo "[CMake] Configuring..."
FC=ifx cmake .. -DCMAKE_BUILD_TYPE=Release

# 3. Compile
echo "[Make] Compiling..."
make -j4

# 4. Run Solver
# We pipe "2" into the solver to auto-select 'Implicit Scheme'
# Change "2" to "1" if you want to run Explicit.
#echo "[Run] Running Explicit Solver..."
echo "[Run] Running Implicit Solver..."
echo "2" | ./solver

# Move result to main directory for the python script
mv results.csv ../results.csv
cd ..

# 5. Plotting
echo "[Plot] Generating Graphs..."
python3 plot_results.py

echo "========================================"
echo "   Done! Check the pop-up window."
echo "========================================"