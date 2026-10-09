// src/Infrastructure/GateManager.h
#pragma once
#include <vector>
#include <cstdint>
#include <tuple>
#include <array>
#include <algorithm>

class GateManager {
public:
    explicit GateManager(size_t chunk_size = 10000000);
    

    void add_X(uint32_t t);
    void add_CNOT(uint32_t c, uint32_t t);
    void add_Toffoli(uint32_t c1, uint32_t c2, uint32_t t);
    

    // (full_depth, current_depth, toffoli_depth, t_count, c_count)

    std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> optimize_and_count(uint32_t width);
    
    void replay_reverse(size_t start_ptr, size_t end_ptr);
    size_t current_pointer() const;
    void clear();
    // Basis-state audit, including a literal inverse of the optimized stream.
    void apply_to_basis(std::vector<uint8_t>& state, bool reverse = false) const;

private:

    size_t chunk_size_; 
    std::vector<std::vector<uint64_t>> chunks;
    std::vector<uint64_t> current_chunk;

    // Bit packing constants
    static const int SHIFT_OP = 60;
    static const int SHIFT_ARG = 30;

    static const uint64_t MASK_ARG = (1ULL << SHIFT_ARG) - 1;
    

    static const uint32_t OP_DATA = 0;
    static const uint32_t OP_X = 1;
    static const uint32_t OP_CNOT = 2;
    static const uint32_t OP_TOFFOLI = 3;
    static const uint32_t OP_NOP = 15;
    
    void _ensure_room(size_t words);
    void _append(uint64_t val);
};
