// src/Infrastructure/GateManager.h
#pragma once
#include <vector>
#include <cstdint>
#include <tuple>
#include <array>
#include <algorithm>

class GateManager {
public:
    GateManager();
    
    // uint32_t inputs suffice because circuit width is normally below four billion.
    void add_X(uint32_t t);
    void add_CNOT(uint32_t c, uint32_t t);
    void add_Toffoli(uint32_t c1, uint32_t c2, uint32_t t);
    
    // Return-type upgrade:
    // (full_depth, current_depth, toffoli_depth, t_count, c_count)
    // use uint64_t throughout to support depths above 4.2 billion.
    std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> optimize_and_count(uint32_t width);
    
    void replay_reverse(size_t start_ptr, size_t end_ptr);
    size_t current_pointer() const;
    void clear();

private:
    // Ten million entries per chunk balance fragmentation and allocation count.
    static const int CHUNK_SIZE = 10000000; 
    std::vector<std::vector<uint64_t>> chunks;
    std::vector<uint64_t> current_chunk;

    // Bit packing constants
    static const int SHIFT_OP = 60;
    static const int SHIFT_ARG = 30;
    // Use 1ULL to ensure a 64-bit unsigned mask.
    static const uint64_t MASK_ARG = (1ULL << SHIFT_ARG) - 1;
    
    // Use uint32_t to avoid signed comparisons.
    static const uint32_t OP_DATA = 0;
    static const uint32_t OP_X = 1;
    static const uint32_t OP_CNOT = 2;
    static const uint32_t OP_TOFFOLI = 3;
    static const uint32_t OP_NOP = 15;
    
    void _append(uint64_t val);
};
