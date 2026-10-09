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

// Packed, little-endian QROM words.  Entry i occupies words_per_entry()
// consecutive uint64_t values and stores x | (y << n).
struct PackedQromTable {
    int address_bits{};
    int word_bits{};
    std::vector<uint64_t> words;

    std::size_t entries() const;
    std::size_t words_per_entry() const;
    bool bit(std::size_t entry, int bit_index) const;
};

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

    // A literal low-width unary-iteration QROM and its gate-by-gate inverse.
    // inverse=true emits the exact reverse word on the same physical lines.
    void coherent_unary_qrom(
        const std::vector<int>& address,
        const std::vector<int>& selector_ancillas,
        const std::vector<int>& data_fanout,
        const std::vector<int>& bus,
        const PackedQromTable& table,
        bool inverse);

    // Windowed QQ point addition with clean z=x_R+x_T and v=y_R+y_T.
    // The QROM cache is uncomputed before Montgomery inversion and its lines
    // are reused by the arithmetic workspace in run_window_layer().
    void early_unlookup_point_add(
        int lanes,
        std::vector<std::vector<int>> x1,
        std::vector<std::vector<int>> y1,
        const std::vector<std::vector<int>>& x_table,
        const std::vector<std::vector<int>>& y_table,
        const std::vector<std::vector<int>>& addresses,
        const std::vector<std::vector<int>>& selectors,
        const std::vector<std::vector<int>>& data_fanout,
        const std::vector<PackedQromTable>& tables,
        const std::vector<int>& ancilla,
        std::vector<std::vector<int>> z,
        std::vector<std::vector<int>> v,
        const std::vector<int>& recovery_copy,
        const std::vector<int>& square_workspace,
        const std::vector<int>& lambda_outputs,
        const std::vector<int>& multiply_outputs,
        const std::vector<int>& inverse_outputs,
        const std::vector<bool>& curve_a);

    // Arithmetic suffix after a QROM pair has produced clean z and v and
    // released all lookup-only lines.  Keeping this as its own entry point
    // lets the modular estimator synthesize the table-independent suffix once
    // for every (n, lanes).
    void early_point_add_tail(
        int lanes,
        std::vector<std::vector<int>> x1,
        std::vector<std::vector<int>> y1,
        std::vector<std::vector<int>> z,
        std::vector<std::vector<int>> v,
        const std::vector<int>& ancilla,
        const std::vector<int>& recovery_copy,
        const std::vector<int>& square_workspace,
        const std::vector<int>& lambda_outputs,
        const std::vector<int>& multiply_outputs,
        const std::vector<int>& inverse_outputs,
        const std::vector<bool>& curve_a);

    LayerStats run_balanced_layer(int lanes, PointAddMode mode);
    LayerStats run_window_layer(const std::vector<PackedQromTable>& tables,
                                const std::vector<bool>& curve_a);
    LayerStats run_window_initialization(const std::vector<PackedQromTable>& tables);
    LayerStats run_qrom_pair(const PackedQromTable& table);
    LayerStats run_early_tail(int lanes, const std::vector<bool>& curve_a);
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
