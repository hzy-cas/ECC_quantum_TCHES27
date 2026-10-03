cmake_minimum_required(VERSION 3.16)
project(windowed_qrom_toptimal LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_executable(retained_estimator
  backend/retained_main.cpp
  backend/src/Infrastructure/GateManager.cpp
  backend/src/Infrastructure/MatrixLoader.cpp
  backend/src/BasicLogic/QuantumBasic.cpp
  backend/src/FieldArithmetic/QuantumField.cpp
  backend/src/FieldArithmetic/QuantumInversion.cpp
  backend/src/AlgorithmFlow/BalancedQuantumShor.cpp
  backend/src/AlgorithmFlow/WindowedQuantumShor.cpp
  backend/src/AlgorithmFlow/BalancedWindowLayer.cpp
  backend/src/AlgorithmFlow/RetainedQuantumShor.cpp)
target_compile_options(retained_estimator PRIVATE -O3 -Wall -Wextra)
enable_testing()
add_executable(test_gate_manager tests/test_gate_manager.cpp backend/src/Infrastructure/GateManager.cpp)
add_test(NAME gate_manager COMMAND test_gate_manager)
