cmake_minimum_required(VERSION 3.16)
project(untitled)

set(CMAKE_CXX_STANDARD 14)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

add_executable(untitled_parallel main_openmp.cpp strategy.cpp matrix.cpp reduce.cpp)

target_link_libraries(untitled_parallel Threads::Threads OpenMP::OpenMP_CXX)