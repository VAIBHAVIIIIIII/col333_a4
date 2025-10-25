#!/bin/bash
g++ -std=c++11 -O2 starter.cpp -o starter
if [ $? -eq 0 ]; then
  echo "Compiled successfully: ./starter"
else
  echo "Compilation failed"
fi
