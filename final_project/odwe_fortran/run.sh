#!/bin/bash
set -e

# --- NEW: Set Number of Threads (Change 12 to your CPU core count) ---
export OMP_NUM_THREADS=9

rm -rf build
mkdir -p build
cd build
cmake ..
cmake --build . -j
./odwe
cd ..

python3 plot_results.py