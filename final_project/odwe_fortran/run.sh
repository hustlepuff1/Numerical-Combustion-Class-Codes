#!/bin/bash
set -e

mkdir -p build
cd build
cmake ..
cmake --build . -j
./odwe
cd ..

python3 plot_results.py