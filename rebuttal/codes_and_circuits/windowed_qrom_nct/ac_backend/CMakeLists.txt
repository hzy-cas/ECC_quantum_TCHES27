cmake_minimum_required(VERSION 3.16)
project(ac_balanced_dw_search LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-O3 -Wall -Wextra)

set(CORE_SOURCES
    src/Infrastructure/GateManager.cpp
    src/Infrastructure/MatrixLoader.cpp
    src/BasicLogic/QuantumBasic.cpp
    src/FieldArithmetic/QuantumField.cpp
    src/FieldArithmetic/QuantumInversion.cpp
    src/AlgorithmFlow/BalancedQuantumShor.cpp
    src/AlgorithmFlow/WindowedQuantumShor.cpp
    src/AlgorithmFlow/BalancedResource.cpp)

add_executable(ac_balanced_estimator src/main.cpp ${CORE_SOURCES})
add_executable(test_schedule tests/test_schedule.cpp src/AlgorithmFlow/BalancedResource.cpp)

enable_testing()
add_test(NAME balanced_schedule COMMAND test_schedule)
