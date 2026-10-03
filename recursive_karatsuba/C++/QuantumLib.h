// QuantumLib.h
#ifndef QUANTUMLIB_H
#define QUANTUMLIB_H

#include <vector>
#include <string>
#include <cstdint> 
#include "GateManager.h"
#include "MatrixLoader.h"

struct ResourceStats {
    uint64_t full_depth;
    uint64_t current_depth;
    uint64_t toffoli_depth;
    uint64_t toffoli;
    uint64_t cnot;
};

class QuantumContext {
public:
    GateManager gm;
    std::vector<int> Toffoli_qubits;
    int cnt;

    QuantumContext() : cnt(0) {}

    // ==========================================
    // 1. Basic Gates & Utils (Implemented in QuantumBasic.cpp)
    // ==========================================
    std::vector<int> X(int a);
    std::vector<int> Toffoli_gate(int a, int b, int c);
    std::vector<int> CNOT(int a, int b);
    void CNOT_n(const std::vector<int>& a, const std::vector<int>& b);
    
    void CONST_ADD_n(const std::vector<int>& target, int val);
    void CONST_ADD_n(const std::vector<int>& target, const std::vector<int>& val_bits);
    
    void Round_constant_XOR(const std::vector<int>& k, int rc, int bit);
    void CSWAP(int c, int t1, int t2);
    std::vector<int> copy_parallel(int value, std::vector<int> ancillas, int n);
    std::vector<int> combine(const std::vector<int>& a, const std::vector<int>& b, const std::vector<int>& r, int n);

    // ==========================================
    // 2. Field Arithmetic (Implemented in QuantumField.cpp)
    // ==========================================
    std::vector<int> Reduction(const std::vector<int>& result);
    void Modular_small(const std::vector<int>& input, const std::vector<int>& result, int size);

    struct KaratsubaResult {
        std::vector<int> res;
        int count;
        std::vector<int> ancilla;
    };
    KaratsubaResult recursive_karatsuba(std::vector<int> a, std::vector<int> b, int n, int count, std::vector<int> ancilla);

    struct SquareResult {
        std::vector<int> r1;
        std::vector<int> r2;
    };
    SquareResult Squaring_plus_input(std::vector<int> x, int n);
    SquareResult Square(std::vector<int> x, int n, int power);

    // ==========================================
    // 3. Inversion (Implemented in QuantumInversion.cpp)
    // ==========================================
    struct InversionResult {
        std::vector<int> res;
        std::vector<int> anc1;
        std::vector<int> anc2;
    };

    InversionResult Inverison_Itoh_Tsujii_based(std::vector<int> a, int n, 
        std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4,
        int count, std::vector<int> ancilla1, std::vector<int> ancilla2);

private:
    InversionResult Inversion_163(std::vector<int> a, std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4, int count, std::vector<int> ancilla1, std::vector<int> ancilla2);
    InversionResult Inversion_233(std::vector<int> a, std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4, int count, std::vector<int> ancilla1, std::vector<int> ancilla2);
    InversionResult Inversion_283(std::vector<int> a, std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4, int count, std::vector<int> ancilla1, std::vector<int> ancilla2);
    InversionResult Inversion_571(std::vector<int> a, std::vector<int> sqr1, std::vector<int> sqr2, std::vector<int> sqr3, std::vector<int> sqr4, int count, std::vector<int> ancilla1, std::vector<int> ancilla2);

public:
    // ==========================================
    // 4. Shor High-Level Logic (Implemented in QuantumShor.cpp)
    // ==========================================
    struct ConsMulResult {
        std::vector<int> product;
        int count;
        std::vector<int> ancilla;
        std::vector<std::vector<std::vector<int>>> tree;
    };
    ConsMulResult cons_mul(std::vector<std::vector<int>> xlist, int count, std::vector<int> ancilla);

    struct MontgomeryResult {
        std::vector<std::vector<int>> lambda_list;
        std::vector<int> ancilla;
    };
    MontgomeryResult Montgomerytrick(std::vector<std::vector<int>> x1x2_list, std::vector<std::vector<int>> y1_list,
        int count, std::vector<int> ancilla, std::vector<int> ancilla2, std::vector<int> sqr, 
        std::vector<int> target_copy, std::vector<int> inv_res);

    void shor_level(int parallel_num, 
                    const std::vector<int>& q_list,
                    std::vector<std::vector<int>> x1_list, std::vector<std::vector<int>> y1_list,
                    const std::vector<int>& x2_vals, const std::vector<int>& y2_vals,
                    std::vector<int> ancilla, std::vector<std::vector<int>> x1x2_list,
                    std::vector<int> ancilla2, std::vector<int> sqr, int level,
                    std::vector<int> lambda_res, std::vector<int> mult_res, std::vector<int> inv_res,
                    std::vector<int> cswap_ancilla);

    ResourceStats shor(int parallel_num, int level, int& length_out);

    void shor_reduction_level(int parallel_num, 
                              std::vector<std::vector<int>> x1_list, std::vector<std::vector<int>> y1_list,
                              std::vector<std::vector<int>> x2_list, std::vector<std::vector<int>> y2_list,
                              std::vector<int> ancilla, std::vector<std::vector<int>> x1x2_list,
                              std::vector<int> ancilla2, std::vector<int> sqr,
                              std::vector<int> lambda_res, std::vector<int> mult_res, std::vector<int> inv_res);

    ResourceStats shor_reduction(int parallel_num, int level, int& length_out);
};

std::vector<int> range_vec(int start, int n);
std::vector<int> concat(const std::vector<int>& a, const std::vector<int>& b);
std::vector<int> slice(const std::vector<int>& v, int start, int len);

#endif