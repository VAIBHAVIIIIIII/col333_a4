#!/bin/bash
# compile.sh
g++ -O2 -std=c++17 -o learn_bn main.cpp
if [ $? -eq 0 ]; then
  echo "Compiled successfully: ./learn_bn"
else
  echo "Compilation failed"
fi
