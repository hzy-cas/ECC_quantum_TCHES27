## Overview
This C++ project uses the Shor algorithm based on recursive karatsuba  algorithm for circuit resource estimation.

## Build and Run
g++ main.cpp GateManager.cpp MatrixLoader.cpp QuantumBasic.cpp QuantumField.cpp QuantumInversion.cpp QuantumShor.cpp -o main -fopenmp -std=c++17 -O3
./main

## Notes
The parameters included are n_target and max_x_limit, the upper bound of the parallel width traversal, both of which are modified in the main function.