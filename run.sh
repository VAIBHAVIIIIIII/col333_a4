#!/bin/bash

# This script executes the compiled C++ program 'learner'.
# $1 is expected to be the .bif file (e.g., hailfinder.bif)
# $2 is expected to be the data file (e.g., records.dat)

# 2. Execute the compiled C++ program, passing the required two input files.
# The C++ program (main.cpp) is responsible for reading these and creating 
# the required output file: solved_hailfinder.bif
./starter $1 $2

echo "Execution complete. Check solved_hailfinder.bif for the learned network."
