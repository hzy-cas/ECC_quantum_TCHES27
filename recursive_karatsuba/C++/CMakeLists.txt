cmake_minimum_required(VERSION 3.16)
project(karatsuba_balanced LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_executable(karatsuba_balanced_estimator main.cpp GateManager.cpp MatrixLoader.cpp
    QuantumBasic.cpp QuantumField.cpp QuantumInversion.cpp QuantumShor.cpp)
target_compile_options(karatsuba_balanced_estimator PRIVATE -O3)
target_compile_definitions(karatsuba_balanced_estimator PRIVATE
    KARATSUBA_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
