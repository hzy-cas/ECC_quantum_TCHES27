#include "QuantumLib.h"
#include <iostream>
#include <vector>
#include "config.h" 
extern GFConfig current_config;

// ==========================================
// Inversion Dispatcher 
// ==========================================
QuantumContext::InversionResult QuantumContext::Inverison_Itoh_Tsujii_based(
    std::vector<int> a, int n, 
    std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4,
    int count, std::vector<int> ancilla1, std::vector<int> ancilla2) 
{
    if (current_config.N == 163) return Inversion_163(a, sqr1, sqr2, sqr3, sqr4, count, ancilla1, ancilla2);
    if (current_config.N == 233) return Inversion_233(a, sqr1, sqr2, sqr3, sqr4, count, ancilla1, ancilla2);
    if (current_config.N == 283) return Inversion_283(a, sqr1, sqr2, sqr3, sqr4, count, ancilla1, ancilla2);
    if (current_config.N == 571) return Inversion_571(a, sqr1, sqr2, sqr3, sqr4, count, ancilla1, ancilla2);
    
    std::cerr << "Inversion implementation for N=" << n << " missing." << std::endl;
    exit(1);
}

// --- Inversion 163 ---
QuantumContext::InversionResult QuantumContext::Inversion_163(std::vector<int> a, std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4, int count, std::vector<int> ancilla1, std::vector<int> ancilla2) 
{
    int n = current_config.N; 

    // a, A_1_1 = Square(a+sqr1, n, 1) 
    auto res = Square(concat(a, sqr1), n, 1);
    a = res.r1; std::vector<int> A_1_1 = res.r2;

    count = 0;
    // A_2_0, count, ancilla1 = recursive_karatsuba(a, A_1_1 , n, count, ancilla1) 
    auto k_res = recursive_karatsuba(a, A_1_1, n, count, ancilla1);
    std::vector<int> A_2_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # chain1
    // A_2_0, A_2_2 = Square(A_2_0+sqr2, n, 2) 
    res = Square(concat(A_2_0, sqr2), n, 2);
    A_2_0 = res.r1; std::vector<int> A_2_2 = res.r2;

    // a, A_1_1 = Square(a+A_1_1, n, 1)  #clear sqr1
    res = Square(concat(a, A_1_1), n, 1);
    a = res.r1; A_1_1 = res.r2;
    sqr1 = A_1_1;

    count = 0;
    // A_4_0, count, ancilla1 = recursive_karatsuba(A_2_2, A_2_0, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_2_2, A_2_0, n, count, ancilla1);
    std::vector<int> A_4_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_2_0, A_2_2 = Square(A_2_0+A_2_2, n, 2)
    res = Square(concat(A_2_0, A_2_2), n, 2);
    A_2_0 = res.r1; A_2_2 = res.r2;
    sqr2 = A_2_2;

    // A_4_0, A_4_4 = Square(A_4_0+sqr1, n, 4)
    res = Square(concat(A_4_0, sqr1), n, 4);
    A_4_0 = res.r1; std::vector<int> A_4_4 = res.r2;

    count = 0;
    // A_8_0, count, ancilla1 = recursive_karatsuba(A_4_0, A_4_4, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_4_0, A_4_4, n, count, ancilla1);
    std::vector<int> A_8_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_4_0, A_4_4 = Square(A_4_0+A_4_4, n, 4)
    res = Square(concat(A_4_0, A_4_4), n, 4);
    A_4_0 = res.r1; A_4_4 = res.r2;
    sqr1 = A_4_4;

    // A_8_0, A_8_8 = Square(A_8_0+sqr2, n, 8)
    res = Square(concat(A_8_0, sqr2), n, 8);
    A_8_0 = res.r1; std::vector<int> A_8_8 = res.r2;

    count = 0;
    // A_16_0, count, ancilla1 = recursive_karatsuba(A_8_0, A_8_8, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_8_0, A_8_8, n, count, ancilla1);
    std::vector<int> A_16_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_8_0, A_8_8 = Square(A_8_0+A_8_8, n, 8)
    res = Square(concat(A_8_0, A_8_8), n, 8);
    A_8_0 = res.r1; A_8_8 = res.r2;
    sqr2 = A_8_8;

    // A_16_0, A_16_16 = Square(A_16_0+sqr1, n, 16)
    res = Square(concat(A_16_0, sqr1), n, 16);
    A_16_0 = res.r1; std::vector<int> A_16_16 = res.r2;

    count = 0;
    // A_32_0, count, ancilla1 = recursive_karatsuba(A_16_0, A_16_16, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_16_0, A_16_16, n, count, ancilla1);
    std::vector<int> A_32_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;
    
    // # S
    // # chain2
    CNOT_n(A_2_0, sqr3);
    // sqr3, A_2_32 = Square(sqr3+sqr4, n, 32)
    res = Square(concat(sqr3, sqr4), n, 32);
    sqr3 = res.r1; std::vector<int> A_2_32 = res.r2;
    CNOT_n(A_2_0, sqr3);

    // # chain1
    // A_16_0, A_16_16 = Square(A_16_0+A_16_16, n, 16)
    res = Square(concat(A_16_0, A_16_16), n, 16);
    A_16_0 = res.r1; A_16_16 = res.r2;
    sqr1 = A_16_16;

    // A_32_0, A_32_32 = Square(A_32_0+sqr2, n, 32)
    res = Square(concat(A_32_0, sqr2), n, 32);
    A_32_0 = res.r1; std::vector<int> A_32_32 = res.r2;

    // # #M
    CNOT_n(A_32_0, sqr3);
    count = 0;
    // A_64_0, count, ancilla1 = recursive_karatsuba(A_32_0, A_32_32, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_32_0, A_32_32, n, count, ancilla1);
    std::vector<int> A_64_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    count = 0;
    // A_34_0, count, ancilla2 = recursive_karatsuba(A_2_32, sqr3, n, count, ancilla2) 
    auto k_res2 = recursive_karatsuba(A_2_32, sqr3, n, count, ancilla2);
    std::vector<int> A_34_0 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    CNOT_n(A_32_0, sqr3);

    // A_32_0, A_32_32 = Square(A_32_0+A_32_32, n, 32)
    res = Square(concat(A_32_0, A_32_32), n, 32);
    A_32_0 = res.r1; A_32_32 = res.r2;
    sqr2 = A_32_32;

    // A_64_0, A_64_64 = Square(A_64_0+sqr1, n, 64)
    res = Square(concat(A_64_0, sqr1), n, 64);
    A_64_0 = res.r1; std::vector<int> A_64_64 = res.r2;

    count = 0;
    // A_128_0, count, ancilla1 = recursive_karatsuba(A_64_0, A_64_64, n, count, ancilla1)  
    k_res = recursive_karatsuba(A_64_0, A_64_64, n, count, ancilla1);
    std::vector<int> A_128_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # chain1
    // A_64_0, A_64_64 = Square(A_64_0+A_64_64, n, 64)
    res = Square(concat(A_64_0, A_64_64), n, 64);
    A_64_0 = res.r1; A_64_64 = res.r2;
    sqr1 = A_64_64; // # clear sqr1

    // A_128_0, A_128_1 = Square(A_128_0+sqr2, n, 1)
    res = Square(concat(A_128_0, sqr2), n, 1);
    A_128_0 = res.r1; std::vector<int> A_128_1 = res.r2;

    // A_34_0, A_34_129 = Square(A_34_0+sqr3, n, 129)
    res = Square(concat(A_34_0, sqr3), n, 129);
    A_34_0 = res.r1; std::vector<int> A_34_129 = res.r2;

    count = 0;
    // A_162_1, count, ancilla2 = recursive_karatsuba(A_34_129, A_128_1, n, count, ancilla2) 
    k_res2 = recursive_karatsuba(A_34_129, A_128_1, n, count, ancilla2);
    std::vector<int> A_162_1 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    // A_2_0, sqr4 = Square(A_2_0+A_2_32, n, 32)
    res = Square(concat(A_2_0, A_2_32), n, 32);
    A_2_0 = res.r1; sqr4 = res.r2;
    return {A_162_1, ancilla1, ancilla2};
}

// --- Inversion 233 ---
QuantumContext::InversionResult QuantumContext::Inversion_233(std::vector<int> a, std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4, int count, std::vector<int> ancilla1, std::vector<int> ancilla2)
{
    int n = current_config.N; 

    // a, A_1_1 = Square(a+sqr1, n, 1) # a^2
    auto res = Square(concat(a, sqr1), n, 1);
    a = res.r1; std::vector<int> A_1_1 = res.r2;

    count = 0;
    // A_2_0, count, ancilla1 = recursive_karatsuba(a, A_1_1 , n, count, ancilla1) 
    auto k_res = recursive_karatsuba(a, A_1_1, n, count, ancilla1);
    std::vector<int> A_2_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # 2^2
    // A_2_0, A_2_2 = Square(A_2_0+sqr2, n, 2)
    res = Square(concat(A_2_0, sqr2), n, 2);
    A_2_0 = res.r1; std::vector<int> A_2_2 = res.r2;

    // a, A_1_1 = Square(a+A_1_1, n, 1) 
    res = Square(concat(a, A_1_1), n, 1);
    a = res.r1; A_1_1 = res.r2;
    sqr1 = A_1_1;

    count = 0;
    // A_4_0, count, ancilla1 = recursive_karatsuba(A_2_2, A_2_0, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_2_2, A_2_0, n, count, ancilla1);
    std::vector<int> A_4_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # 2^4
    // A_2_0, A_2_2 = Square(A_2_0+A_2_2, n, 2)
    res = Square(concat(A_2_0, A_2_2), n, 2);
    A_2_0 = res.r1; A_2_2 = res.r2;
    sqr2 = A_2_2;

    // A_4_0, A_4_4 = Square(A_4_0+sqr1, n, 4)
    res = Square(concat(A_4_0, sqr1), n, 4);
    A_4_0 = res.r1; std::vector<int> A_4_4 = res.r2;

    count = 0;
    // A_8_0, count, ancilla1 = recursive_karatsuba(A_4_0, A_4_4, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_4_0, A_4_4, n, count, ancilla1);
    std::vector<int> A_8_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # chain1
    // A_4_0, A_4_4 = Square(A_4_0+A_4_4, n, 4)
    res = Square(concat(A_4_0, A_4_4), n, 4);
    A_4_0 = res.r1; A_4_4 = res.r2;
    sqr1 = A_4_4;

    // A_8_0, A_8_8 = Square(A_8_0+sqr2, n, 8)
    res = Square(concat(A_8_0, sqr2), n, 8);
    A_8_0 = res.r1; std::vector<int> A_8_8 = res.r2;

    // #M
    count = 0;
    // A_16_0, count, ancilla1 = recursive_karatsuba(A_8_0, A_8_8, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_8_0, A_8_8, n, count, ancilla1);
    std::vector<int> A_16_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # 2^16
    // A_8_0, A_8_8 = Square(A_8_0+A_8_8, n, 8)
    res = Square(concat(A_8_0, A_8_8), n, 8);
    A_8_0 = res.r1; A_8_8 = res.r2;
    sqr2 = A_8_8;

    // A_16_0, A_16_16 = Square(A_16_0+sqr1, n, 16)
    res = Square(concat(A_16_0, sqr1), n, 16);
    A_16_0 = res.r1; std::vector<int> A_16_16 = res.r2;

    count = 0;
    // A_32_0, count, ancilla1 = recursive_karatsuba(A_16_0, A_16_16, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_16_0, A_16_16, n, count, ancilla1);
    std::vector<int> A_32_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # S
    // # chain2
    CNOT_n(A_8_0, sqr3);
    // sqr3, A_8_32 = Square(sqr3+sqr4, n, 32)
    res = Square(concat(sqr3, sqr4), n, 32);
    sqr3 = res.r1; std::vector<int> A_8_32 = res.r2;
    CNOT_n(A_8_0, sqr3);

    // # chain1
    // A_16_0, A_16_16 = Square(A_16_0+A_16_16, n, 16)
    res = Square(concat(A_16_0, A_16_16), n, 16);
    A_16_0 = res.r1; A_16_16 = res.r2;
    sqr1 = A_16_16;

    // A_32_0, A_32_32 = Square(A_32_0+sqr2, n, 32)
    res = Square(concat(A_32_0, sqr2), n, 32);
    A_32_0 = res.r1; std::vector<int> A_32_32 = res.r2;

    // #M
    CNOT_n(A_32_0, sqr3);
    count = 0;
    // A_64_0, count, ancilla1 = recursive_karatsuba(A_32_0, A_32_32, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_32_0, A_32_32, n, count, ancilla1);
    std::vector<int> A_64_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    count = 0;
    // A_40_0, count, ancilla2 = recursive_karatsuba(A_8_32, sqr3, n, count, ancilla2) 
    auto k_res2 = recursive_karatsuba(A_8_32, sqr3, n, count, ancilla2);
    std::vector<int> A_40_0 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    CNOT_n(A_32_0, sqr3);

    // # S
    // # chain2
    // A_8_0,sqr4 = Square(A_8_0+A_8_32, n, 32) # clear sqr4
    res = Square(concat(A_8_0, A_8_32), n, 32);
    A_8_0 = res.r1; sqr4 = res.r2;

    // A_40_0, A_40_64 = Square(A_40_0+sqr3, n, 64)
    res = Square(concat(A_40_0, sqr3), n, 64);
    A_40_0 = res.r1; std::vector<int> A_40_64 = res.r2;

    // # chain1
    // A_32_0, A_32_32 = Square(A_32_0+A_32_32, n, 32)
    res = Square(concat(A_32_0, A_32_32), n, 32);
    A_32_0 = res.r1; A_32_32 = res.r2;
    sqr2 = A_32_32;

    // A_64_0, A_64_64 = Square(A_64_0+sqr1, n, 64)
    res = Square(concat(A_64_0, sqr1), n, 64);
    A_64_0 = res.r1; std::vector<int> A_64_64 = res.r2;

    // #M
    CNOT_n(A_64_0, sqr4);
    count = 0;
    // A_128_0, count, ancilla1 = recursive_karatsuba(A_64_0, A_64_64, n, count, ancilla1)  
    k_res = recursive_karatsuba(A_64_0, A_64_64, n, count, ancilla1);
    std::vector<int> A_128_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    count = 0;
    // A_104_0, count, ancilla2 = recursive_karatsuba(A_40_64, sqr4, n, count, ancilla2) 
    k_res2 = recursive_karatsuba(A_40_64, sqr4, n, count, ancilla2);
    std::vector<int> A_104_0 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;
    
    CNOT_n(A_64_0, sqr4);
    
    // # chain1
    // A_64_0, sqr1 = Square(A_64_0+A_64_64, n, 64) # clear sqr1 
    res = Square(concat(A_64_0, A_64_64), n, 64);
    A_64_0 = res.r1; sqr1 = res.r2;

    // A_128_0, A_128_1 = Square(A_128_0+sqr2, n, 1)
    res = Square(concat(A_128_0, sqr2), n, 1);
    A_128_0 = res.r1; std::vector<int> A_128_1 = res.r2;

    // A_104_0, A_104_129 = Square(A_104_0+sqr4, n, 129)
    res = Square(concat(A_104_0, sqr4), n, 129);
    A_104_0 = res.r1; std::vector<int> A_104_129 = res.r2;

    count = 0;
    // A_232_1, count, ancilla2 = recursive_karatsuba(A_104_129, A_128_1, n, count, ancilla2) 
    k_res2 = recursive_karatsuba(A_104_129, A_128_1, n, count, ancilla2);
    std::vector<int> A_232_1 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    // A_40_0, sqr3 = Square(A_40_0+A_40_64, n, 64)# clear sqr3 
    res = Square(concat(A_40_0, A_40_64), n, 64);
    A_40_0 = res.r1; sqr3 = res.r2;

    return {A_232_1, ancilla1, ancilla2};
}

// --- Inversion 283 ---
QuantumContext::InversionResult QuantumContext::Inversion_283(std::vector<int> a, std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4, int count, std::vector<int> ancilla1, std::vector<int> ancilla2)
{
    int n = current_config.N;

    // a, A_1_1 = Square(a+sqr1, n, 1) # a^2
    auto res = Square(concat(a, sqr1), n, 1);
    a = res.r1; std::vector<int> A_1_1 = res.r2;

    count = 0;
    // A_2_0, count, ancilla1 = recursive_karatsuba(a, A_1_1 , n, count, ancilla1) 
    auto k_res = recursive_karatsuba(a, A_1_1, n, count, ancilla1);
    std::vector<int> A_2_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # chain1
    // A_2_0, A_2_2 = Square(A_2_0+sqr2, n, 2)
    res = Square(concat(A_2_0, sqr2), n, 2);
    A_2_0 = res.r1; std::vector<int> A_2_2 = res.r2;

    // a, A_1_1 = Square(a+A_1_1, n, 1) 
    res = Square(concat(a, A_1_1), n, 1);
    a = res.r1; A_1_1 = res.r2;
    sqr1 = A_1_1;

    count = 0;
    // A_4_0, count, ancilla1 = recursive_karatsuba(A_2_2, A_2_0, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_2_2, A_2_0, n, count, ancilla1);
    std::vector<int> A_4_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_2_0, A_2_2 = Square(A_2_0+A_2_2, n, 2)
    res = Square(concat(A_2_0, A_2_2), n, 2);
    A_2_0 = res.r1; A_2_2 = res.r2;
    sqr2 = A_2_2;

    // A_4_0, A_4_4 = Square(A_4_0+sqr1, n, 4)
    res = Square(concat(A_4_0, sqr1), n, 4);
    A_4_0 = res.r1; std::vector<int> A_4_4 = res.r2;

    count = 0;
    // A_8_0, count, ancilla1 = recursive_karatsuba(A_4_0, A_4_4, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_4_0, A_4_4, n, count, ancilla1);
    std::vector<int> A_8_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # S
    // # chain2
    CNOT_n(A_2_0, sqr3);
    // sqr3, A_2_8 = Square(sqr3+sqr4, n, 8)
    res = Square(concat(sqr3, sqr4), n, 8);
    sqr3 = res.r1; std::vector<int> A_2_8 = res.r2;
    CNOT_n(A_2_0, sqr3);

    // # chain1
    // A_4_0, A_4_4 = Square(A_4_0+A_4_4, n, 4)
    res = Square(concat(A_4_0, A_4_4), n, 4);
    A_4_0 = res.r1; A_4_4 = res.r2;

    // A_8_0, A_8_8 = Square(A_8_0+sqr2, n, 8)
    res = Square(concat(A_8_0, sqr2), n, 8);
    A_8_0 = res.r1; std::vector<int> A_8_8 = res.r2;

    // # #M
    CNOT_n(A_8_0, sqr3);
    count = 0;
    // A_16_0, count, ancilla1 = recursive_karatsuba(A_8_0, A_8_8, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_8_0, A_8_8, n, count, ancilla1);
    std::vector<int> A_16_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    count = 0;
    // A_10_0, count, ancilla2 = recursive_karatsuba(A_2_8, sqr3, n, count, ancilla2) 
    auto k_res2 = recursive_karatsuba(A_2_8, sqr3, n, count, ancilla2);
    std::vector<int> A_10_0 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    CNOT_n(A_8_0, sqr3);

    // # S
    // # chain2
    // A_2_0,sqr4 = Square(A_2_0+A_2_8, n, 8) # clear sqr4
    res = Square(concat(A_2_0, A_2_8), n, 8);
    A_2_0 = res.r1; sqr4 = res.r2;

    // A_10_0, A_10_16 = Square(A_10_0+sqr3, n, 16)
    res = Square(concat(A_10_0, sqr3), n, 16);
    A_10_0 = res.r1; std::vector<int> A_10_16 = res.r2;

    // # chain1    
    // A_8_0, sqr2 = Square(A_8_0+A_8_8, n, 8)
    res = Square(concat(A_8_0, A_8_8), n, 8);
    A_8_0 = res.r1; sqr2 = res.r2;

    // A_16_0, A_16_16 = Square(A_16_0+sqr1, n, 16)
    res = Square(concat(A_16_0, sqr1), n, 16);
    A_16_0 = res.r1; std::vector<int> A_16_16 = res.r2;

    // #M
    CNOT_n(A_16_0, sqr4);
    count = 0;
    // A_32_0, count, ancilla1 = recursive_karatsuba(A_16_0, A_16_16, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_16_0, A_16_16, n, count, ancilla1);
    std::vector<int> A_32_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    count = 0;
    // A_26_0, count, ancilla2 = recursive_karatsuba(A_10_16, sqr4, n, count, ancilla2) 
    k_res2 = recursive_karatsuba(A_10_16, sqr4, n, count, ancilla2);
    std::vector<int> A_26_0 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    CNOT_n(A_16_0, sqr4);

    // A_16_0, A_16_16 = Square(A_16_0+A_16_16, n, 16)
    res = Square(concat(A_16_0, A_16_16), n, 16);
    A_16_0 = res.r1; A_16_16 = res.r2;
    sqr1 = A_16_16;

    // A_32_0, A_32_32 = Square(A_32_0+sqr2, n, 32)
    res = Square(concat(A_32_0, sqr2), n, 32);
    A_32_0 = res.r1; std::vector<int> A_32_32 = res.r2;

    count = 0;
    // A_64_0, count, ancilla1 = recursive_karatsuba(A_32_0, A_32_32, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_32_0, A_32_32, n, count, ancilla1);
    std::vector<int> A_64_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_32_0, A_32_32 = Square(A_32_0+A_32_32, n, 32)
    res = Square(concat(A_32_0, A_32_32), n, 32);
    A_32_0 = res.r1; A_32_32 = res.r2;
    sqr2 = A_32_32;

    // A_64_0, A_64_64 = Square(A_64_0+sqr1, n, 64)
    res = Square(concat(A_64_0, sqr1), n, 64);
    A_64_0 = res.r1; std::vector<int> A_64_64 = res.r2;

    count = 0;
    // A_128_0, count, ancilla1 = recursive_karatsuba(A_64_0, A_64_64, n, count, ancilla1)  
    k_res = recursive_karatsuba(A_64_0, A_64_64, n, count, ancilla1);
    std::vector<int> A_128_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_64_0, A_64_64 = Square(A_64_0+A_64_64, n, 64)
    res = Square(concat(A_64_0, A_64_64), n, 64);
    A_64_0 = res.r1; A_64_64 = res.r2;
    sqr1 = A_64_64;

    // A_128_0, A_128_128 = Square(A_128_0+sqr2, n, 128)
    res = Square(concat(A_128_0, sqr2), n, 128);
    A_128_0 = res.r1; std::vector<int> A_128_128 = res.r2;

    count = 0;
    // A_256_0, count, ancilla1 = recursive_karatsuba(A_128_0, A_128_128, n, count, ancilla1)  
    k_res = recursive_karatsuba(A_128_0, A_128_128, n, count, ancilla1);
    std::vector<int> A_256_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_128_0, sqr2 = Square(A_128_0+A_128_128, n, 128)
    res = Square(concat(A_128_0, A_128_128), n, 128);
    A_128_0 = res.r1; sqr2 = res.r2;

    // A_256_0, A_256_1 = Square(A_256_0+sqr1, n, 1)
    res = Square(concat(A_256_0, sqr1), n, 1);
    A_256_0 = res.r1; std::vector<int> A_256_1 = res.r2;

    // A_26_0, A_26_257 = Square(A_26_0+sqr4, n, 257)
    res = Square(concat(A_26_0, sqr4), n, 257);
    A_26_0 = res.r1; std::vector<int> A_26_257 = res.r2;

    count = 0;
    // A_282_1, count, ancilla2 = recursive_karatsuba(A_26_257, A_256_1, n, count, ancilla2) 
    k_res2 = recursive_karatsuba(A_26_257, A_256_1, n, count, ancilla2);
    std::vector<int> A_282_1 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    // A_10_0, sqr3 = Square(A_10_0+A_10_16, n, 16)
    res = Square(concat(A_10_0, A_10_16), n, 16);
    A_10_0 = res.r1; sqr3 = res.r2;
    // # clear sqr2
    
    return {A_282_1, ancilla1, ancilla2};
}

// --- Inversion 571 ---
QuantumContext::InversionResult QuantumContext::Inversion_571(std::vector<int> a, std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4, int count, std::vector<int> ancilla1, std::vector<int> ancilla2)
{
    int n = current_config.N;

    // a, A_1_1 = Square(a+sqr1, n, 1) # a^2
    auto res = Square(concat(a, sqr1), n, 1);
    a = res.r1; std::vector<int> A_1_1 = res.r2;

    count = 0;
    // A_2_0, count, ancilla1 = recursive_karatsuba(a, A_1_1 , n, count, ancilla1) 
    auto k_res = recursive_karatsuba(a, A_1_1, n, count, ancilla1);
    std::vector<int> A_2_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # chain1
    // A_2_0, A_2_2 = Square(A_2_0+sqr2, n, 2)
    res = Square(concat(A_2_0, sqr2), n, 2);
    A_2_0 = res.r1; std::vector<int> A_2_2 = res.r2;

    // a, sqr1 = Square(a+A_1_1, n, 1) 
    res = Square(concat(a, A_1_1), n, 1);
    a = res.r1; sqr1 = res.r2;

    count = 0;
    // A_4_0, count, ancilla1 = recursive_karatsuba(A_2_2, A_2_0, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_2_2, A_2_0, n, count, ancilla1);
    std::vector<int> A_4_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_2_0, sqr2 = Square(A_2_0+A_2_2, n, 2)
    res = Square(concat(A_2_0, A_2_2), n, 2);
    A_2_0 = res.r1; sqr2 = res.r2;

    // A_4_0, A_4_4 = Square(A_4_0+sqr1, n, 4)
    res = Square(concat(A_4_0, sqr1), n, 4);
    A_4_0 = res.r1; std::vector<int> A_4_4 = res.r2;

    count = 0;
    // A_8_0, count, ancilla1 = recursive_karatsuba(A_4_0, A_4_4, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_4_0, A_4_4, n, count, ancilla1);
    std::vector<int> A_8_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // # S
    // # chain2
    CNOT_n(A_2_0, sqr3);
    // sqr3, A_2_8 = Square(sqr3+sqr4, n, 8)
    res = Square(concat(sqr3, sqr4), n, 8);
    sqr3 = res.r1; std::vector<int> A_2_8 = res.r2;
    CNOT_n(A_2_0, sqr3);

    // # chain1
    // A_4_0, A_4_4 = Square(A_4_0+A_4_4, n, 4)
    res = Square(concat(A_4_0, A_4_4), n, 4);
    A_4_0 = res.r1; A_4_4 = res.r2;

    // A_8_0, A_8_8 = Square(A_8_0+sqr2, n, 8)
    res = Square(concat(A_8_0, sqr2), n, 8);
    A_8_0 = res.r1; std::vector<int> A_8_8 = res.r2;

    // #M
    CNOT_n(A_8_0, sqr3);
    count = 0;
    // A_16_0, count, ancilla1 = recursive_karatsuba(A_8_0, A_8_8, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_8_0, A_8_8, n, count, ancilla1);
    std::vector<int> A_16_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    count = 0;
    // A_10_0, count, ancilla2 = recursive_karatsuba(A_2_8, sqr3, n, count, ancilla2) 
    auto k_res2 = recursive_karatsuba(A_2_8, sqr3, n, count, ancilla2);
    std::vector<int> A_10_0 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    CNOT_n(A_8_0, sqr3);

    // # S
    // # chain2
    // A_2_0,sqr4 = Square(A_2_0+A_2_8, n, 8) # clear sqr4
    res = Square(concat(A_2_0, A_2_8), n, 8);
    A_2_0 = res.r1; sqr4 = res.r2;

    // A_10_0, A_10_16 = Square(A_10_0+sqr3, n, 16)
    res = Square(concat(A_10_0, sqr3), n, 16);
    A_10_0 = res.r1; std::vector<int> A_10_16 = res.r2;

    // # chain1    
    // A_8_0, A_8_8 = Square(A_8_0+A_8_8, n, 8)
    res = Square(concat(A_8_0, A_8_8), n, 8);
    A_8_0 = res.r1; A_8_8 = res.r2;
    sqr2 = A_8_8;

    // A_16_0, A_16_16 = Square(A_16_0+sqr1, n, 16)
    res = Square(concat(A_16_0, sqr1), n, 16);
    A_16_0 = res.r1; std::vector<int> A_16_16 = res.r2;

    // #M
    CNOT_n(A_16_0, sqr4);
    count = 0;
    // A_32_0, count, ancilla1 = recursive_karatsuba(A_16_0, A_16_16, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_16_0, A_16_16, n, count, ancilla1);
    std::vector<int> A_32_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    count = 0;
    // A_26_0, count, ancilla2 = recursive_karatsuba(A_10_16, sqr4, n, count, ancilla2) 
    k_res2 = recursive_karatsuba(A_10_16, sqr4, n, count, ancilla2);
    std::vector<int> A_26_0 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    CNOT_n(A_16_0, sqr4);

    // # S
    // # chain2
    // A_26_0,A_26_32= Square(A_26_0+sqr4, n, 32) 
    res = Square(concat(A_26_0, sqr4), n, 32);
    A_26_0 = res.r1; std::vector<int> A_26_32 = res.r2;

    // A_10_0, sqr3 = Square(A_10_0+A_10_16, n, 16)# clear sqr3
    res = Square(concat(A_10_0, A_10_16), n, 16);
    A_10_0 = res.r1; sqr3 = res.r2;

    // # chain1
    // A_16_0, sqr1 = Square(A_16_0+A_16_16, n, 16)
    res = Square(concat(A_16_0, A_16_16), n, 16);
    A_16_0 = res.r1; sqr1 = res.r2;

    // A_32_0, A_32_32 = Square(A_32_0+sqr2, n, 32)
    res = Square(concat(A_32_0, sqr2), n, 32);
    A_32_0 = res.r1; std::vector<int> A_32_32 = res.r2;

    // #M
    CNOT_n(A_32_0, sqr3);
    count = 0;
    // A_64_0, count, ancilla1 = recursive_karatsuba(A_32_0, A_32_32, n, count, ancilla1) 
    k_res = recursive_karatsuba(A_32_0, A_32_32, n, count, ancilla1);
    std::vector<int> A_64_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    count = 0;
    // A_58_0, count, ancilla2 = recursive_karatsuba(A_26_32, sqr3, n, count, ancilla2) 
    k_res2 = recursive_karatsuba(A_26_32, sqr3, n, count, ancilla2);
    std::vector<int> A_58_0 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    CNOT_n(A_32_0, sqr3);

    // A_32_0, sqr2 = Square(A_32_0+A_32_32, n, 32)
    res = Square(concat(A_32_0, A_32_32), n, 32);
    A_32_0 = res.r1; sqr2 = res.r2;

    // A_64_0, A_64_64 = Square(A_64_0+sqr1, n, 64)
    res = Square(concat(A_64_0, sqr1), n, 64);
    A_64_0 = res.r1; std::vector<int> A_64_64 = res.r2;

    count = 0;
    // A_128_0, count, ancilla1 = recursive_karatsuba(A_64_0, A_64_64, n, count, ancilla1)  
    k_res = recursive_karatsuba(A_64_0, A_64_64, n, count, ancilla1);
    std::vector<int> A_128_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_64_0, sqr1 = Square(A_64_0+A_64_64, n, 64)
    res = Square(concat(A_64_0, A_64_64), n, 64);
    A_64_0 = res.r1; sqr1 = res.r2;

    // A_128_0, A_128_128 = Square(A_128_0+sqr2, n, 128)
    res = Square(concat(A_128_0, sqr2), n, 128);
    A_128_0 = res.r1; std::vector<int> A_128_128 = res.r2;

    count = 0;
    // A_256_0, count, ancilla1 = recursive_karatsuba(A_128_0, A_128_128, n, count, ancilla1)  
    k_res = recursive_karatsuba(A_128_0, A_128_128, n, count, ancilla1);
    std::vector<int> A_256_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_128_0, sqr2 = Square(A_128_0+A_128_128, n, 128)
    res = Square(concat(A_128_0, A_128_128), n, 128);
    A_128_0 = res.r1; sqr2 = res.r2;

    // A_256_0, A_256_256 = Square(A_256_0+sqr1, n, 256)
    res = Square(concat(A_256_0, sqr1), n, 256);
    A_256_0 = res.r1; std::vector<int> A_256_256 = res.r2;

    count = 0;
    // A_512_0, count, ancilla1 = recursive_karatsuba(A_256_0, A_256_256, n, count, ancilla1)  
    k_res = recursive_karatsuba(A_256_0, A_256_256, n, count, ancilla1);
    std::vector<int> A_512_0 = Reduction(k_res.res);
    ancilla1 = k_res.ancilla;

    // A_256_0, sqr1 = Square(A_256_0+A_256_256, n, 256)# clear sqr1
    res = Square(concat(A_256_0, A_256_256), n, 256);
    A_256_0 = res.r1; sqr1 = res.r2;

    // A_512_0, A_512_1 = Square(A_512_0+sqr2, n, 1)
    res = Square(concat(A_512_0, sqr2), n, 1);
    A_512_0 = res.r1; std::vector<int> A_512_1 = res.r2;

    // A_58_0, A_58_513 = Square(A_58_0+sqr3, n, 513)
    res = Square(concat(A_58_0, sqr3), n, 513);
    A_58_0 = res.r1; std::vector<int> A_58_513 = res.r2;

    count = 0;
    // A_570_1, count, ancilla2 = recursive_karatsuba(A_58_513, A_512_1, n, count, ancilla2) 
    k_res2 = recursive_karatsuba(A_58_513, A_512_1, n, count, ancilla2);
    std::vector<int> A_570_1 = Reduction(k_res2.res);
    ancilla2 = k_res2.ancilla;

    // A_26_0,sqr4= Square(A_26_0+A_26_32, n, 32) 
    res = Square(concat(A_26_0, A_26_32), n, 32);
    A_26_0 = res.r1; sqr4 = res.r2;

    return {A_570_1, ancilla1, ancilla2};
}