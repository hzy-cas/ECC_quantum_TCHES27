cmake_minimum_required(VERSION 3.10)
project(ShorEstimator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(ShorEstimator ${SOURCES})
target_compile_definitions(ShorEstimator PRIVATE AC_SOURCE_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
