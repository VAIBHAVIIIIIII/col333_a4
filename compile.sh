#!/bin/bash
#!/bin/bash

# Define the name of the main source file
SOURCE_FILE="starter.cpp"

# Define the output executable name
EXECUTABLE_NAME="starter"

# Compile the code using g++ with C++17 standard, optimization (-O3), 
# and all warnings enabled (-Wall).
g++ -std=c++17 -O3 -Wall $SOURCE_FILE -o $EXECUTABLE_NAME

if [ $? -eq 0 ]; then
  echo "Compiled successfully: ./starter"
else
  echo "Compilation failed"
fi
