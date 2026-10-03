#pragma once
#include "../Infrastructure/GateManager.h"
#include "../Infrastructure/MatrixLoader.h"
#include "../Infrastructure/config.h"
#include <vector>
#include <string>

class QuantumContext {
public:
    QuantumContext(GFConfig cfg);

    // Basic Ops
    void X(int a);
    void CNOT(int c, int t);
    void Toffoli(int c1, int c2, int t);
    void CNOT_N(const std::vector<int>& a, const std::vector<int>& b);
    void CONST_ADD(const std::vector<int>& reg, const std::vector<bool>& val_bits);
    void CSWAP(int c, int a, int b);
    std::vector<int> copy_parallel(int ctrl, const std::vector<int>& ancilla, int count);
    
    // Matrix Logic
    std::vector<int> mul_matrix(const std::vector<int>& x, int n_rows, std::string filename);
    std::pair<std::vector<int>, std::vector<int>> Square(const std::vector<int>& x, int n_rows, int power);
    std::pair<std::vector<int>, std::vector<int>> Squaring_plus_input(const std::vector<int>& x, int n_rows);
    
    // Field Arithmetic
    std::vector<int> mul(const std::vector<int>& a0, const std::vector<int>& b0, const std::vector<int>& c0, int offset);
    
    // Inversion
    std::vector<int> inversion_logic(const std::vector<int>& a, 
                                     const std::vector<int>& sqr, 
                                     const std::vector<int>& ancilla, 
                                     const std::vector<int>& ancilla2, 
                                     int& count);

    // Shor Logic
    struct TreeResult {
        std::vector<int> product;
        int count;
        std::vector<int> ancilla;
        std::vector<std::vector<std::vector<int>>> tree;
    };
    
    TreeResult cons_mul(std::vector<std::vector<int>> xlist, int count, std::vector<int> ancilla);
    
    std::tuple<std::vector<std::vector<int>>, std::vector<int>, int> 
    Montgomerytrick(std::vector<std::vector<int>> x1x2_list, std::vector<std::vector<int>> y1_list, 
                    int count, std::vector<int> ancilla, std::vector<int> ancilla2, std::vector<int> sqr);

    void shor_accumulation_step(int parallel_num, 
                                const std::vector<int>& q_list,
                                std::vector<std::vector<int>>& x1_list,
                                std::vector<std::vector<int>>& y1_list,
                                const std::vector<std::vector<bool>>& x2_val_list,
                                const std::vector<std::vector<bool>>& y2_val_list,
                                std::vector<int>& ancilla,
                                std::vector<std::vector<int>>& x1x2_list,
                                std::vector<int>& ancilla2,
                                std::vector<int>& sqr,
                                std::vector<int>& cswap_ancilla);

    void shor_reduction_step(int parallel_num,
                             std::vector<std::vector<int>>& x1_list,
                             std::vector<std::vector<int>>& y1_list,
                             std::vector<std::vector<int>>& x2_list,
                             std::vector<std::vector<int>>& y2_list,
                             std::vector<int>& ancilla,
                             std::vector<std::vector<int>>& x1x2_list,
                             std::vector<int>& ancilla2,
                             std::vector<int>& sqr);

    std::tuple<int, uint64_t, uint64_t, uint32_t, uint32_t, uint32_t> 
    run_shor_logic(int parallel_num, std::string mode);
    
    void estimate_total_resources(int total_inputs);

private:
    GFConfig config;
    GateManager gm;
    MatrixCache matrix_loader;
    std::vector<int> Toffoli_qubits;
    
    std::string get_matrix_path(std::string folder, std::string filename);
    int get_mul_total_size();
    
    // Inversion sub-routines
    std::vector<int> inv_163(std::vector<int> a, std::vector<int> sqr, std::vector<int> anc, std::vector<int> anc2, int& cnt);
    std::vector<int> inv_233(std::vector<int> a, std::vector<int> sqr, std::vector<int> anc, std::vector<int> anc2, int& cnt);
    std::vector<int> inv_283(std::vector<int> a, std::vector<int> sqr, std::vector<int> anc, std::vector<int> anc2, int& cnt);
    std::vector<int> inv_571(std::vector<int> a, std::vector<int> sqr, std::vector<int> anc, std::vector<int> anc2, int& cnt);
};