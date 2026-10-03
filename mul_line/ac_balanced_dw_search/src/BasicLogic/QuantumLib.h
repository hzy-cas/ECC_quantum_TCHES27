#pragma once

#include "../Infrastructure/GateManager.h"
#include "../Infrastructure/MatrixLoader.h"
#include "../Infrastructure/config.h"

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

enum class PointAddMode { Accumulation, Reduction };

struct LayerStats {
    uint64_t qubits{};
    uint64_t toffoli{};
    uint64_t cnot{};
    uint64_t full_depth{};
    uint64_t current_depth{};
    uint64_t toffoli_depth{};
};

using PrimitiveStats = LayerStats;

class QuantumContext {
public:
    explicit QuantumContext(GFConfig cfg);

    void X(int a);
    void CNOT(int c, int t);
    void Toffoli(int c1, int c2, int t);
    void CNOT_N(const std::vector<int>& a, const std::vector<int>& b);
    void CONST_ADD(const std::vector<int>& reg, const std::vector<bool>& val_bits);
    void CSWAP(int c, int a, int b);
    std::vector<int> copy_parallel(int ctrl, const std::vector<int>& ancilla, int count);

    std::vector<int> mul_matrix(const std::vector<int>& x, int n_rows, std::string filename);
    std::pair<std::vector<int>, std::vector<int>> Square(const std::vector<int>& x, int n_rows, int power);
    std::pair<std::vector<int>, std::vector<int>> Squaring_plus_input(const std::vector<int>& x, int n_rows);
    std::vector<int> mul(const std::vector<int>& a0, const std::vector<int>& b0,
                         const std::vector<int>& c0, int offset);
    std::vector<int> inversion_logic(const std::vector<int>& a,
                                     const std::vector<int>& sqr,
                                     const std::vector<int>& ancilla,
                                     const std::vector<int>& ancilla2,
                                     int& count);

    struct TreeResult {
        std::vector<int> product;
        int count{};
        std::vector<std::vector<std::vector<int>>> tree;
    };

    struct MontgomeryResult {
        std::vector<std::vector<int>> lambdas;
    };

    TreeResult cons_mul(const std::vector<std::vector<int>>& inputs, int count,
                        const std::vector<int>& ancilla);
    MontgomeryResult montgomery_clean(
        const std::vector<std::vector<int>>& denominators,
        const std::vector<std::vector<int>>& numerators,
        const std::vector<int>& ancilla,
        const std::vector<int>& recovery_copy,
        const std::vector<int>& square_workspace,
        const std::vector<int>& lambda_outputs,
        const std::vector<int>& inverse_outputs);

    void accumulation_point_add(
        int lanes, const std::vector<int>& controls,
        std::vector<std::vector<int>> x1, std::vector<std::vector<int>> y1,
        const std::vector<std::vector<bool>>& x2_constants,
        const std::vector<std::vector<bool>>& y2_constants,
        const std::vector<int>& ancilla,
        std::vector<std::vector<int>> x1_plus_x2,
        const std::vector<int>& recovery_copy,
        const std::vector<int>& square_workspace,
        const std::vector<int>& lambda_outputs,
        const std::vector<int>& multiply_outputs,
        const std::vector<int>& inverse_outputs,
        const std::vector<int>& cswap_controls);

    void reduction_point_add(
        int lanes,
        std::vector<std::vector<int>> x1, std::vector<std::vector<int>> y1,
        std::vector<std::vector<int>> x2, std::vector<std::vector<int>> y2,
        const std::vector<int>& ancilla,
        std::vector<std::vector<int>> x1_plus_x2,
        const std::vector<int>& recovery_copy,
        const std::vector<int>& square_workspace,
        const std::vector<int>& lambda_outputs,
        const std::vector<int>& multiply_outputs,
        const std::vector<int>& inverse_outputs);

    LayerStats run_balanced_layer(int lanes, PointAddMode mode);
    PrimitiveStats estimate_multiplication();
    PrimitiveStats estimate_inversion();

private:
    GFConfig config;
    GateManager gm;
    MatrixCache matrix_loader;
    std::vector<int> Toffoli_qubits;

    std::string get_matrix_path(std::string folder, std::string filename);
    int get_mul_total_size() const { return config.mul_cost_size; }

    std::vector<int> inv_163(std::vector<int> a, std::vector<int> sqr,
                             std::vector<int> anc, std::vector<int> anc2, int& cnt);
    std::vector<int> inv_233(std::vector<int> a, std::vector<int> sqr,
                             std::vector<int> anc, std::vector<int> anc2, int& cnt);
    std::vector<int> inv_283(std::vector<int> a, std::vector<int> sqr,
                             std::vector<int> anc, std::vector<int> anc2, int& cnt);
    std::vector<int> inv_571(std::vector<int> a, std::vector<int> sqr,
                             std::vector<int> anc, std::vector<int> anc2, int& cnt);
};
