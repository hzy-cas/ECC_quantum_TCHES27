#!/bin/bash

# --- Configuration Area ---

# Project root directory (i.e., the directory where this script is located)
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT_DIR="$PROJECT_DIR/input"
RESULTS_DIR="$PROJECT_DIR/results"
# The name of the executable file generated after compiling C++.
EXECUTABLE_NAME="untitled_parallel"

# --- Script begins ---

set -e

echo "# --- Script begins ----"

cd "$PROJECT_DIR"
echo "--- Compiling project... ---"
rm -rf build
mkdir build
cd build
cmake ..
make
echo "✅ Compilation successful."

cd "$PROJECT_DIR"
mkdir -p "$RESULTS_DIR"
echo "✅All results will be saved in: $RESULTS_DIR"


echo "--- Begin parallel processing of all files in the '$INPUT_DIR' directory ---"

run_single_file() {
    input_file="$1"
    base_name=$(basename "$input_file")
    
    echo "▶️ Processing: $base_name"
    
    "$PROJECT_DIR/build/$EXECUTABLE_NAME" "$input_file"
    
    mv "$PROJECT_DIR/result_$base_name" "$RESULTS_DIR/result_$base_name"
    
    echo "✅ Processing complete! Results saved to: $RESULTS_DIR/result_$base_name"
}

export -f run_single_file
export PROJECT_DIR
export EXECUTABLE_NAME
export RESULTS_DIR

find "$INPUT_DIR" -type f | xargs -I {} -P 1 bash -c 'run_single_file "{}"'
echo "-----------------------------------------------------"
echo "All files processed!" 
echo "--- Script execution finished ---"