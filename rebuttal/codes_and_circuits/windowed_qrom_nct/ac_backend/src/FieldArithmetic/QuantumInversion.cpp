// src/FieldArithmetic/QuantumInversion.cpp
#include "../BasicLogic/QuantumLib.h"
#include <iostream>
#include <tuple>

// ==========================================

// ==========================================
namespace {
    std::vector<int> concat(const std::vector<int>& a, const std::vector<int>& b) {
        std::vector<int> res;
        res.reserve(a.size() + b.size());
        res.insert(res.end(), a.begin(), a.end());
        res.insert(res.end(), b.begin(), b.end());
        return res;
    }
}

// ==========================================

// ==========================================
// S(x, p) -> Square(x, 2*n, p)
#define S_MACRO(x, p) Square(x, 2 * config.n, p)

// M(x, y, c_off) -> mul(..., Toffoli_qubits[c_off:], 0)

#define M_MACRO(x, y, c_off) mul(x, y, std::vector<int>(Toffoli_qubits.begin() + (c_off), Toffoli_qubits.end()), 0)

// C(x, y) -> CNOT_N(x, y)
#define C_MACRO(x, y) CNOT_N(x, y)


// ==========================================

// ==========================================
std::vector<int> QuantumContext::inversion_logic(const std::vector<int>& a, 
                                                 const std::vector<int>& sqr, 
                                                 const std::vector<int>& ancilla, 
                                                 const std::vector<int>& ancilla2, 
                                                 int& count) {
    if (config.n == 163) return inv_163(a, sqr, ancilla, ancilla2, count);
    if (config.n == 233) return inv_233(a, sqr, ancilla, ancilla2, count);
    if (config.n == 283) return inv_283(a, sqr, ancilla, ancilla2, count);
    if (config.n == 571) return inv_571(a, sqr, ancilla, ancilla2, count);
    
    std::cerr << "Error: Unsupported N for inversion: " << config.n << std::endl;
    return {};
}

// ==========================================
// GF(163) Inversion Implementation
// Ref: algorithms/implementations.py -> Shor163
// ==========================================
std::vector<int> QuantumContext::inv_163(std::vector<int> a, std::vector<int> sqr, 
                                         std::vector<int> ancilla, std::vector<int> ancilla2, int& count) {
    int n = 163;
    int mul_size = config.mul_cost_size; // 906
    int anc_limit = config.block_size;   // 743

    // Slices for sqr (simulating sqr[:163], sqr[163:326], etc.)
    std::vector<int> sqr1(sqr.begin(), sqr.begin() + n);
    std::vector<int> sqr2(sqr.begin() + n, sqr.begin() + 2 * n);
    std::vector<int> sqr3(sqr.begin() + 2 * n, sqr.begin() + 3 * n);
    std::vector<int> sqr4(sqr.begin() + 3 * n, sqr.begin() + 4 * n);

    // Ancilla split
    std::vector<int> anc_main_1(ancilla.begin(), ancilla.begin() + anc_limit);
    std::vector<int> anc_main_2(ancilla.begin() + anc_limit, ancilla.end());
    
    std::vector<int> anc_side_1(ancilla2.begin(), ancilla2.begin() + anc_limit);
    std::vector<int> anc_side_2(ancilla2.begin() + anc_limit, ancilla2.end());

    std::vector<int> A_1_1, A_2_0, A_2_2, A_4_0, A_4_4, A_8_0, A_8_8, A_16_0, A_16_16, A_32_0, A_32_32;
    std::vector<int> A_2_32, A_64_0, A_64_64, A_34_0, A_128_0, A_128_1, A_34_129, A_162_1;

    // 1 -> 2
    std::tie(a, A_1_1) = S_MACRO(concat(a, sqr1), 1);
    A_2_0 = M_MACRO(concat(a, anc_main_1), concat(A_1_1, anc_main_2), count);
    count += mul_size;

    // Chain 1 (2->4)
    std::tie(A_2_0, A_2_2) = S_MACRO(concat(A_2_0, sqr2), 2);
    std::tie(a, A_1_1) = S_MACRO(concat(a, A_1_1), 1); // clear sqr1
    sqr1 = A_1_1;
    A_4_0 = M_MACRO(concat(A_2_2, anc_main_1), concat(A_2_0, anc_main_2), count);
    count += mul_size;

    // 4 -> 8
    std::tie(A_2_0, A_2_2) = S_MACRO(concat(A_2_0, A_2_2), 2);
    sqr2 = A_2_2;
    std::tie(A_4_0, A_4_4) = S_MACRO(concat(A_4_0, sqr1), 4);
    A_8_0 = M_MACRO(concat(A_4_0, anc_main_1), concat(A_4_4, anc_main_2), count);
    count += mul_size;

    // 8 -> 16
    std::tie(A_4_0, A_4_4) = S_MACRO(concat(A_4_0, A_4_4), 4);
    sqr1 = A_4_4;
    std::tie(A_8_0, A_8_8) = S_MACRO(concat(A_8_0, sqr2), 8);
    A_16_0 = M_MACRO(concat(A_8_0, anc_main_1), concat(A_8_8, anc_main_2), count);
    count += mul_size;

    // 16 -> 32
    std::tie(A_8_0, A_8_8) = S_MACRO(concat(A_8_0, A_8_8), 8);
    sqr2 = A_8_8;
    std::tie(A_16_0, A_16_16) = S_MACRO(concat(A_16_0, sqr1), 16);
    A_32_0 = M_MACRO(concat(A_16_0, anc_main_1), concat(A_16_16, anc_main_2), count);
    count += mul_size;

    // S - Chain 2 (Prep 2->34 side branch)
    C_MACRO(A_2_0, sqr3);
    std::tie(sqr3, A_2_32) = S_MACRO(concat(sqr3, sqr4), 32);
    C_MACRO(A_2_0, sqr3);

    // Chain 1 (32 -> 64)
    std::tie(A_16_0, A_16_16) = S_MACRO(concat(A_16_0, A_16_16), 16);
    sqr1 = A_16_16;
    std::tie(A_32_0, A_32_32) = S_MACRO(concat(A_32_0, sqr2), 32);

    // Merge (32 + 2 = 34)
    C_MACRO(A_32_0, sqr3);
    // Main
    A_64_0 = M_MACRO(concat(A_32_0, anc_main_1), concat(A_32_32, anc_main_2), count);
    count += mul_size;
    // Side (A_2_32 * sqr3) -> A_34_0
    A_34_0 = M_MACRO(concat(A_2_32, anc_side_1), concat(sqr3, anc_side_2), count);
    count += mul_size;
    C_MACRO(A_32_0, sqr3);

    // 64 -> 128
    std::tie(A_32_0, A_32_32) = S_MACRO(concat(A_32_0, A_32_32), 32);
    sqr2 = A_32_32;
    std::tie(A_64_0, A_64_64) = S_MACRO(concat(A_64_0, sqr1), 64);
    A_128_0 = M_MACRO(concat(A_64_0, anc_main_1), concat(A_64_64, anc_main_2), count);
    count += mul_size;

    // Chain 1 Cleanup
    std::tie(A_64_0, A_64_64) = S_MACRO(concat(A_64_0, A_64_64), 64);
    sqr1 = A_64_64; // clear sqr1

    // Final Merge Prep (128->1, 34->129)
    std::tie(A_128_0, A_128_1) = S_MACRO(concat(A_128_0, sqr2), 1);
    std::tie(A_34_0, A_34_129) = S_MACRO(concat(A_34_0, sqr3), 129);

    // Final Mul
    A_162_1 = M_MACRO(concat(A_34_129, anc_side_1), concat(A_128_1, anc_side_2), count);
    count += mul_size;

    // Cleanup A_2_0
    std::tie(A_2_0, sqr4) = S_MACRO(concat(A_2_0, A_2_32), 32);

    return A_162_1;
}

// ==========================================
// GF(233) Inversion Implementation
// Ref: algorithms/implementations.py -> Shor233
// ==========================================
std::vector<int> QuantumContext::inv_233(std::vector<int> a, std::vector<int> sqr, 
                                         std::vector<int> ancilla, std::vector<int> ancilla2, int& count) {
    int n = 233;
    int mul_size = config.mul_cost_size;
    int anc_limit = config.block_size;

    std::vector<int> sqr1(sqr.begin(), sqr.begin() + n);
    std::vector<int> sqr2(sqr.begin() + n, sqr.begin() + 2 * n);
    std::vector<int> sqr3(sqr.begin() + 2 * n, sqr.begin() + 3 * n);
    std::vector<int> sqr4(sqr.begin() + 3 * n, sqr.begin() + 4 * n);

    std::vector<int> anc_main_1(ancilla.begin(), ancilla.begin() + anc_limit);
    std::vector<int> anc_main_2(ancilla.begin() + anc_limit, ancilla.end());
    std::vector<int> anc_side_1(ancilla2.begin(), ancilla2.begin() + anc_limit);
    std::vector<int> anc_side_2(ancilla2.begin() + anc_limit, ancilla2.end());

    std::vector<int> A_1_1, A_2_0, A_2_2, A_4_0, A_4_4, A_8_0, A_8_8, A_16_0, A_16_16;
    std::vector<int> A_32_0, A_32_32, A_8_32, A_64_0, A_64_64, A_40_0, A_40_64;
    std::vector<int> A_128_0, A_104_0, A_128_1, A_104_129, A_232_1;

    // 1 -> 2
    std::tie(a, A_1_1) = S_MACRO(concat(a, sqr1), 1);
    A_2_0 = M_MACRO(concat(a, anc_main_1), concat(A_1_1, anc_main_2), count);
    count += mul_size;

    // 2 -> 4
    std::tie(A_2_0, A_2_2) = S_MACRO(concat(A_2_0, sqr2), 2);
    std::tie(a, A_1_1) = S_MACRO(concat(a, A_1_1), 1);
    sqr1 = A_1_1;
    A_4_0 = M_MACRO(concat(A_2_2, anc_main_1), concat(A_2_0, anc_main_2), count);
    count += mul_size;

    // 4 -> 8
    std::tie(A_2_0, A_2_2) = S_MACRO(concat(A_2_0, A_2_2), 2);
    sqr2 = A_2_2;
    std::tie(A_4_0, A_4_4) = S_MACRO(concat(A_4_0, sqr1), 4);
    A_8_0 = M_MACRO(concat(A_4_0, anc_main_1), concat(A_4_4, anc_main_2), count);
    count += mul_size;

    // Chain 1 (8 -> 16)
    std::tie(A_4_0, A_4_4) = S_MACRO(concat(A_4_0, A_4_4), 4);
    sqr1 = A_4_4;
    std::tie(A_8_0, A_8_8) = S_MACRO(concat(A_8_0, sqr2), 8);
    A_16_0 = M_MACRO(concat(A_8_0, anc_main_1), concat(A_8_8, anc_main_2), count);
    count += mul_size;

    // 16 -> 32
    std::tie(A_8_0, A_8_8) = S_MACRO(concat(A_8_0, A_8_8), 8);
    sqr2 = A_8_8;
    std::tie(A_16_0, A_16_16) = S_MACRO(concat(A_16_0, sqr1), 16);
    A_32_0 = M_MACRO(concat(A_16_0, anc_main_1), concat(A_16_16, anc_main_2), count);
    count += mul_size;

    // S - Chain 2 (Prep 8 -> 40 side branch)
    C_MACRO(A_8_0, sqr3);
    std::tie(sqr3, A_8_32) = S_MACRO(concat(sqr3, sqr4), 32);
    C_MACRO(A_8_0, sqr3);

    // Chain 1 (32 -> 64)
    std::tie(A_16_0, A_16_16) = S_MACRO(concat(A_16_0, A_16_16), 16);
    sqr1 = A_16_16;
    std::tie(A_32_0, A_32_32) = S_MACRO(concat(A_32_0, sqr2), 32);

    // Merge (32 + 8 = 40)
    C_MACRO(A_32_0, sqr3);
    // Main
    A_64_0 = M_MACRO(concat(A_32_0, anc_main_1), concat(A_32_32, anc_main_2), count);
    count += mul_size;
    // Side (A_8_32 * sqr3) -> A_40_0
    A_40_0 = M_MACRO(concat(A_8_32, anc_side_1), concat(sqr3, anc_side_2), count);
    count += mul_size;
    C_MACRO(A_32_0, sqr3);

    // S - Chain 2 (Prep 40 -> 104 side branch)
    std::tie(A_8_0, sqr4) = S_MACRO(concat(A_8_0, A_8_32), 32); // clear sqr4
    std::tie(A_40_0, A_40_64) = S_MACRO(concat(A_40_0, sqr3), 64);

    // Chain 1 (64 -> 128)
    std::tie(A_32_0, A_32_32) = S_MACRO(concat(A_32_0, A_32_32), 32);
    sqr2 = A_32_32;
    std::tie(A_64_0, A_64_64) = S_MACRO(concat(A_64_0, sqr1), 64);

    // Merge (64 + 40 = 104)
    C_MACRO(A_64_0, sqr4);
    // Main
    A_128_0 = M_MACRO(concat(A_64_0, anc_main_1), concat(A_64_64, anc_main_2), count);
    count += mul_size;
    // Side
    A_104_0 = M_MACRO(concat(A_40_64, anc_side_1), concat(sqr4, anc_side_2), count);
    count += mul_size;
    C_MACRO(A_64_0, sqr4);

    // Chain 1 cleanup & Final Merge Prep
    std::tie(A_64_0, sqr1) = S_MACRO(concat(A_64_0, A_64_64), 64); // clear sqr1
    std::tie(A_128_0, A_128_1) = S_MACRO(concat(A_128_0, sqr2), 1);
    std::tie(A_104_0, A_104_129) = S_MACRO(concat(A_104_0, sqr4), 129);

    // Final Mul
    A_232_1 = M_MACRO(concat(A_104_129, anc_side_1), concat(A_128_1, anc_side_2), count);
    count += mul_size;

    std::tie(A_40_0, sqr3) = S_MACRO(concat(A_40_0, A_40_64), 64); // clear sqr3

    return A_232_1;
}

// ==========================================
// GF(283) Inversion Implementation
// Ref: algorithms/implementations.py -> Shor283
// ==========================================
std::vector<int> QuantumContext::inv_283(std::vector<int> a, std::vector<int> sqr, 
                                         std::vector<int> ancilla, std::vector<int> ancilla2, int& count) {
    int n = 283;
    int mul_size = config.mul_cost_size;
    int anc_limit = config.block_size;

    std::vector<int> sqr1(sqr.begin(), sqr.begin() + n);
    std::vector<int> sqr2(sqr.begin() + n, sqr.begin() + 2 * n);
    std::vector<int> sqr3(sqr.begin() + 2 * n, sqr.begin() + 3 * n);
    std::vector<int> sqr4(sqr.begin() + 3 * n, sqr.begin() + 4 * n);

    std::vector<int> anc_main_1(ancilla.begin(), ancilla.begin() + anc_limit);
    std::vector<int> anc_main_2(ancilla.begin() + anc_limit, ancilla.end());
    std::vector<int> anc_side_1(ancilla2.begin(), ancilla2.begin() + anc_limit);
    std::vector<int> anc_side_2(ancilla2.begin() + anc_limit, ancilla2.end());

    std::vector<int> A_1_1, A_2_0, A_2_2, A_4_0, A_4_4, A_8_0, A_8_8;
    std::vector<int> A_2_8, A_16_0, A_16_16, A_10_0, A_10_16, A_32_0, A_32_32, A_26_0;
    std::vector<int> A_64_0, A_64_64, A_128_0, A_128_128, A_256_0, A_256_1, A_26_257, A_282_1;

    // 1 -> 2
    std::tie(a, A_1_1) = S_MACRO(concat(a, sqr1), 1);
    A_2_0 = M_MACRO(concat(a, anc_main_1), concat(A_1_1, anc_main_2), count);
    count += mul_size;

    // Chain 1 (2 -> 4)
    std::tie(A_2_0, A_2_2) = S_MACRO(concat(A_2_0, sqr2), 2);
    std::tie(a, A_1_1) = S_MACRO(concat(a, A_1_1), 1);
    sqr1 = A_1_1;
    A_4_0 = M_MACRO(concat(A_2_2, anc_main_1), concat(A_2_0, anc_main_2), count);
    count += mul_size;

    // 4 -> 8
    std::tie(A_2_0, A_2_2) = S_MACRO(concat(A_2_0, A_2_2), 2);
    sqr2 = A_2_2;
    std::tie(A_4_0, A_4_4) = S_MACRO(concat(A_4_0, sqr1), 4);
    A_8_0 = M_MACRO(concat(A_4_0, anc_main_1), concat(A_4_4, anc_main_2), count);
    count += mul_size;

    // S - Chain 2 (Prep Side Branch)
    C_MACRO(A_2_0, sqr3);
    std::tie(sqr3, A_2_8) = S_MACRO(concat(sqr3, sqr4), 8);
    C_MACRO(A_2_0, sqr3);

    // Chain 1 (8 -> 16)
    std::tie(A_4_0, A_4_4) = S_MACRO(concat(A_4_0, A_4_4), 4);
    std::tie(A_8_0, A_8_8) = S_MACRO(concat(A_8_0, sqr2), 8);

    // Merge
    C_MACRO(A_8_0, sqr3);
    // Main Branch
    A_16_0 = M_MACRO(concat(A_8_0, anc_main_1), concat(A_8_8, anc_main_2), count);
    count += mul_size;
    // Side Branch
    A_10_0 = M_MACRO(concat(A_2_8, anc_side_1), concat(sqr3, anc_side_2), count);
    count += mul_size;
    C_MACRO(A_8_0, sqr3);

    // S - Chain 2
    std::tie(A_2_0, sqr4) = S_MACRO(concat(A_2_0, A_2_8), 8); // clear sqr4
    std::tie(A_10_0, A_10_16) = S_MACRO(concat(A_10_0, sqr3), 16);

    // Chain 1
    std::tie(A_8_0, sqr2) = S_MACRO(concat(A_8_0, A_8_8), 8);
    std::tie(A_16_0, A_16_16) = S_MACRO(concat(A_16_0, sqr1), 16);

    // Merge
    C_MACRO(A_16_0, sqr4);
    // Main Branch (16->32)
    A_32_0 = M_MACRO(concat(A_16_0, anc_main_1), concat(A_16_16, anc_main_2), count);
    count += mul_size;
    // Side Branch (10->26)
    A_26_0 = M_MACRO(concat(A_10_16, anc_side_1), concat(sqr4, anc_side_2), count);
    count += mul_size;
    C_MACRO(A_16_0, sqr4);

    // Chain 1 (32 -> 64)
    std::tie(A_16_0, A_16_16) = S_MACRO(concat(A_16_0, A_16_16), 16);
    sqr1 = A_16_16;
    std::tie(A_32_0, A_32_32) = S_MACRO(concat(A_32_0, sqr2), 32);
    A_64_0 = M_MACRO(concat(A_32_0, anc_main_1), concat(A_32_32, anc_main_2), count);
    count += mul_size;

    // Chain 1 (64 -> 128)
    std::tie(A_32_0, A_32_32) = S_MACRO(concat(A_32_0, A_32_32), 32);
    sqr2 = A_32_32;
    std::tie(A_64_0, A_64_64) = S_MACRO(concat(A_64_0, sqr1), 64);
    A_128_0 = M_MACRO(concat(A_64_0, anc_main_1), concat(A_64_64, anc_main_2), count);
    count += mul_size;

    // Chain 1 (128 -> 256)
    std::tie(A_64_0, A_64_64) = S_MACRO(concat(A_64_0, A_64_64), 64);
    sqr1 = A_64_64;
    std::tie(A_128_0, A_128_128) = S_MACRO(concat(A_128_0, sqr2), 128);
    A_256_0 = M_MACRO(concat(A_128_0, anc_main_1), concat(A_128_128, anc_main_2), count);
    count += mul_size;

    // Final Prep
    std::tie(A_128_0, sqr2) = S_MACRO(concat(A_128_0, A_128_128), 128); // clear sqr2
    std::tie(A_256_0, A_256_1) = S_MACRO(concat(A_256_0, sqr1), 1);
    std::tie(A_26_0, A_26_257) = S_MACRO(concat(A_26_0, sqr4), 257);

    // Final Mul (Combine 256 + 26 = 282)
    A_282_1 = M_MACRO(concat(A_26_257, anc_side_1), concat(A_256_1, anc_side_2), count);
    count += mul_size;

    // Cleanup sqr3
    std::tie(A_10_0, sqr3) = S_MACRO(concat(A_10_0, A_10_16), 16);

    return A_282_1;
}

// ==========================================
// GF(571) Inversion Implementation
// Ref: algorithms/implementations.py -> Shor571
// ==========================================
// src/FieldArithmetic/QuantumInversion.cpp



std::vector<int> QuantumContext::inv_571(std::vector<int> a, std::vector<int> sqr, std::vector<int> anc, std::vector<int> anc2, int& count) {

    int n = 571;
    int sz = 2998; // ANC_LIMIT
    int mul_size = 3569;

    auto sqr1 = std::vector<int>(sqr.begin(), sqr.begin() + n);
    auto sqr2 = std::vector<int>(sqr.begin() + n, sqr.begin() + 2 * n);
    auto sqr3 = std::vector<int>(sqr.begin() + 2 * n, sqr.begin() + 3 * n);
    auto sqr4 = std::vector<int>(sqr.begin() + 3 * n, sqr.begin() + 4 * n);

    auto anc_main_1 = std::vector<int>(anc.begin(), anc.begin() + sz);
    auto anc_main_2 = std::vector<int>(anc.begin() + sz, anc.end());
    

    auto anc_side_1 = std::vector<int>(anc2.begin(), anc2.begin() + sz);
    auto anc_side_2 = std::vector<int>(anc2.begin() + sz, anc2.end());

    // 1 -> 2
    auto [res_a, A_1_1] = S_MACRO(concat(a, sqr1), 1);
    a = res_a; // Update a
    auto A_2_0 = M_MACRO(concat(a, anc_main_1), concat(A_1_1, anc_main_2), count);
    count += mul_size;

    // Chain 1 (2 -> 4)
    std::vector<int> A_2_2;
    std::tie(A_2_0, A_2_2) = S_MACRO(concat(A_2_0, sqr2), 2);
    
    std::tie(a, sqr1) = S_MACRO(concat(a, A_1_1), 1); // Python: a, sqr1 = S(a+A_1_1, 1) (Updating sqr1!)
    
    auto A_4_0 = M_MACRO(concat(A_2_2, anc_main_1), concat(A_2_0, anc_main_2), count);
    count += mul_size;

    // 4 -> 8
    std::tie(A_2_0, sqr2) = S_MACRO(concat(A_2_0, A_2_2), 2); // Update sqr2
    
    std::vector<int> A_4_4;
    std::tie(A_4_0, A_4_4) = S_MACRO(concat(A_4_0, sqr1), 4);
    
    auto A_8_0 = M_MACRO(concat(A_4_0, anc_main_1), concat(A_4_4, anc_main_2), count);
    count += mul_size;

    // S - Chain 2 (Prep Side)
    C_MACRO(A_2_0, sqr3);
    std::vector<int> A_2_8;
    std::tie(sqr3, A_2_8) = S_MACRO(concat(sqr3, sqr4), 8); // Update sqr3
    C_MACRO(A_2_0, sqr3); // CNOT restore

    // Chain 1 (8 -> 16)
    std::tie(A_4_0, A_4_4) = S_MACRO(concat(A_4_0, A_4_4), 4);
    std::vector<int> A_8_8;
    std::tie(A_8_0, A_8_8) = S_MACRO(concat(A_8_0, sqr2), 8);

    // Merge
    C_MACRO(A_8_0, sqr3);
    // Main
    auto A_16_0 = M_MACRO(concat(A_8_0, anc_main_1), concat(A_8_8, anc_main_2), count);
    count += mul_size;
    // Side
    auto A_10_0 = M_MACRO(concat(A_2_8, anc_side_1), concat(sqr3, anc_side_2), count);
    count += mul_size;
    C_MACRO(A_8_0, sqr3);

    // S - Chain 2
    std::tie(A_2_0, sqr4) = S_MACRO(concat(A_2_0, A_2_8), 8); // clear sqr4
    std::vector<int> A_10_16;
    std::tie(A_10_0, A_10_16) = S_MACRO(concat(A_10_0, sqr3), 16);

    // Chain 1
    std::tie(A_8_0, A_8_8) = S_MACRO(concat(A_8_0, A_8_8), 8);
    sqr2 = A_8_8; // Update sqr2
    std::vector<int> A_16_16;
    std::tie(A_16_0, A_16_16) = S_MACRO(concat(A_16_0, sqr1), 16);

    // Merge
    C_MACRO(A_16_0, sqr4);
    // Main
    auto A_32_0 = M_MACRO(concat(A_16_0, anc_main_1), concat(A_16_16, anc_main_2), count);
    count += mul_size;
    // Side
    auto A_26_0 = M_MACRO(concat(A_10_16, anc_side_1), concat(sqr4, anc_side_2), count);
    count += mul_size;
    C_MACRO(A_16_0, sqr4);

    // S - Chain 2
    std::vector<int> A_26_32;
    std::tie(A_26_0, A_26_32) = S_MACRO(concat(A_26_0, sqr4), 32);
    std::tie(A_10_0, sqr3) = S_MACRO(concat(A_10_0, A_10_16), 16); // clear sqr3

    // Chain 1 (32 -> 64)
    std::tie(A_16_0, sqr1) = S_MACRO(concat(A_16_0, A_16_16), 16); // Update sqr1
    std::vector<int> A_32_32;
    std::tie(A_32_0, A_32_32) = S_MACRO(concat(A_32_0, sqr2), 32);

    // Merge
    C_MACRO(A_32_0, sqr3);
    // Main
    auto A_64_0 = M_MACRO(concat(A_32_0, anc_main_1), concat(A_32_32, anc_main_2), count);
    count += mul_size;
    // Side (Note: A_26_32 + sqr3)
    auto A_58_0 = M_MACRO(concat(A_26_32, anc_side_1), concat(sqr3, anc_side_2), count);
    count += mul_size;
    C_MACRO(A_32_0, sqr3);

    // Chain 1 (64 -> 128)
    std::tie(A_32_0, sqr2) = S_MACRO(concat(A_32_0, A_32_32), 32); // Update sqr2
    std::vector<int> A_64_64;
    std::tie(A_64_0, A_64_64) = S_MACRO(concat(A_64_0, sqr1), 64);
    auto A_128_0 = M_MACRO(concat(A_64_0, anc_main_1), concat(A_64_64, anc_main_2), count);
    count += mul_size;

    // Chain 1 (128 -> 256)
    std::tie(A_64_0, sqr1) = S_MACRO(concat(A_64_0, A_64_64), 64); // Update sqr1
    std::vector<int> A_128_128;
    std::tie(A_128_0, A_128_128) = S_MACRO(concat(A_128_0, sqr2), 128);
    auto A_256_0 = M_MACRO(concat(A_128_0, anc_main_1), concat(A_128_128, anc_main_2), count);
    count += mul_size;

    // Chain 1 (256 -> 512)
    std::tie(A_128_0, sqr2) = S_MACRO(concat(A_128_0, A_128_128), 128); // Update sqr2
    std::vector<int> A_256_256;
    std::tie(A_256_0, A_256_256) = S_MACRO(concat(A_256_0, sqr1), 256);
    auto A_512_0 = M_MACRO(concat(A_256_0, anc_main_1), concat(A_256_256, anc_main_2), count);
    count += mul_size;

    // Final Prep
    std::tie(A_256_0, sqr1) = S_MACRO(concat(A_256_0, A_256_256), 256); // clear sqr1
    
    std::vector<int> A_512_1;
    std::tie(A_512_0, A_512_1) = S_MACRO(concat(A_512_0, sqr2), 1);
    
    std::vector<int> A_58_513;

    std::tie(A_58_0, A_58_513) = S_MACRO(concat(A_58_0, sqr3), 513);

    // Final Mul (Combine 58 + 512 = 570)
    auto A_570_1 = M_MACRO(concat(A_58_513, anc_side_1), concat(A_512_1, anc_side_2), count);
    count += mul_size;

    // Cleanup sqr4
    std::tie(A_26_0, sqr4) = S_MACRO(concat(A_26_0, A_26_32), 32);

    return A_570_1;
}