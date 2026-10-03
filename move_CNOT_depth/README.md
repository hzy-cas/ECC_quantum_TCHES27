## Overview

This project aims to optimize the CNOT circuit depth of sparse linear layer matrices in batches.It automatically processes input matrix files and generates optimization results.

## File Structure

* `run_all.sh`: Script to batch compile and process all input files.
* `input/`: Directory for storing matrix files to be processed.
* `results/`: Directory where all processing results are saved.

## Dependencies

* CMake
* GCC or Clang (with OpenMP support)

## Build and Run

1. Place all input matrix files into the `input/` directory.
2. Run the following command in the project root directory:
```bash
./run_all.sh

```



The script will automatically compile the code and process all input files. The results will be saved in the `results/` directory.

## Notes

The entry point for the C++ main program is `main_openmp.cpp`, which supports multi-threaded optimization. You can modify `TOTAL_LOOPS` in the source code to adjust the number of optimization iterations.

## Configuration and Matrix Sizing

If the number of columns in your input matrix is , you must adjust the definitions in the source code.

Change the default values:

```cpp
#define SIZE 64*11
#define WORDS 11

```

To fit your matrix size:

```cpp
#define SIZE (( (n + 63) / 64 ) * 64 )
#define WORDS (( (n + 63) / 64 ))

```

**Explanation:**

*  is the actual width (number of columns) of your matrix.
* `SIZE` must be set to the smallest multiple of 64 that is greater than or equal to .
* `WORDS` corresponds to `ceil(n/64)`.
