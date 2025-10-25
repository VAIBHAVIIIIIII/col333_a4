#!/bin/bash

# $1 is expected to be the .bif file (e.g., hailfinder.bif)
# $2 is expected to be the data file (e.g., records.dat)
# Execute the compiled C++ program, passing the required two input files.
echo "Execution starting..."
./starter $1 $2
echo "Execution complete. Check solved_hailfinder.bif for the learned network."
