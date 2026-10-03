cmake_minimum_required(VERSION 3.10)
project(DaCSynth)

set(CMAKE_CXX_STANDARD 17)
add_executable(dac_synth
    main.cpp
    utils.cpp
    hopcroft_karp.cpp
    dac_solver.cpp
    io_handler.cpp
)