#include "../BasicLogic/QuantumLib.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace {

template <typename T>
std::vector<T> concat(const std::vector<T>& left, const std::vector<T>& right) {
    std::vector<T> out;
    out.reserve(left.size() + right.size());
    out.insert(out.end(), left.begin(), left.end());
    out.insert(out.end(), right.begin(), right.end());
    return out;
}

template <typename T>
std::vector<T> slice_exact(const std::vector<T>& values, std::size_t start,
                           std::size_t length, const char* label) {
    if (start > values.size() || length > values.size() - start) {
        throw std::out_of_range(std::string("slice outside ") + label);
    }
    return std::vector<T>(values.begin() + static_cast<std::ptrdiff_t>(start),
                          values.begin() + static_cast<std::ptrdiff_t>(start + length));
}

std::vector<int> allocate_range(uint64_t& cursor, uint64_t length) {
    if (cursor + length > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("qubit index exceeds signed-int range");
    }
    std::vector<int> out;
    out.reserve(static_cast<std::size_t>(length));
    for (uint64_t i = 0; i < length; ++i) {
        out.push_back(static_cast<int>(cursor + i));
    }
    cursor += length;
    return out;
}

std::vector<std::vector<int>> split_register(const std::vector<int>& flat, int width,
                                              int count, const char* label) {
    if (width < 0 || count < 0 ||
        flat.size() != static_cast<std::size_t>(width) * static_cast<std::size_t>(count)) {
        throw std::invalid_argument(std::string("invalid register layout: ") + label);
    }
    std::vector<std::vector<int>> out;
    out.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        out.push_back(slice_exact(flat, static_cast<std::size_t>(i) * width,
                                  static_cast<std::size_t>(width), label));
    }
    return out;
}

std::vector<bool> constant_bits(int n, uint64_t value) {
    std::vector<bool> bits(static_cast<std::size_t>(n), false);
    for (int i = 0; i < n && i < 64; ++i) {
        bits[static_cast<std::size_t>(i)] = ((value >> i) & 1U) != 0;
    }
    return bits;
}

LayerStats unpack_stats(uint64_t qubits,
                        const std::tuple<uint64_t, uint64_t, uint64_t,
                                         uint64_t, uint64_t>& stats) {
    return LayerStats{qubits,
                      std::get<3>(stats),
                      std::get<4>(stats),
                      std::get<0>(stats),
                      std::get<1>(stats),
                      std::get<2>(stats)};
}

} // namespace

QuantumContext::TreeResult QuantumContext::cons_mul(
    const std::vector<std::vector<int>>& inputs, int count,
    const std::vector<int>& ancilla) {
    if (inputs.empty()) {
        throw std::invalid_argument("Montgomery product tree needs at least one input");
    }

    const int operand_workspace = config.block_size;
    const int pair_workspace = 2 * operand_workspace;
    const int mul_size = config.mul_cost_size;

    std::vector<std::vector<std::vector<int>>> tree;
    tree.push_back(inputs);
    auto current = inputs;

    while (current.size() > 1) {
        std::vector<std::vector<int>> next;
        next.reserve((current.size() + 1) / 2);
        std::size_t pair = 0;

        for (std::size_t i = 0; i < current.size(); i += 2) {
            if (i + 1 == current.size()) {
                next.push_back(current[i]);
                continue;
            }

            const std::size_t base = pair * static_cast<std::size_t>(pair_workspace);
            auto left = concat(current[i],
                               slice_exact(ancilla, base, operand_workspace, "product-tree ancilla"));
            auto right = concat(current[i + 1],
                                slice_exact(ancilla, base + operand_workspace,
                                            operand_workspace, "product-tree ancilla"));
            auto targets = slice_exact(Toffoli_qubits, static_cast<std::size_t>(count),
                                       mul_size, "product-tree Toffoli targets");
            next.push_back(mul(left, right, targets, 0));
            count += mul_size;
            ++pair;
        }

        tree.push_back(next);
        current = std::move(next);
    }

    return TreeResult{tree.back().front(), count, std::move(tree)};
}

QuantumContext::MontgomeryResult QuantumContext::montgomery_clean(
    const std::vector<std::vector<int>>& denominators,
    const std::vector<std::vector<int>>& numerators,
    const std::vector<int>& ancilla,
    const std::vector<int>& recovery_copy,
    const std::vector<int>& square_workspace,
    const std::vector<int>& lambda_outputs,
    const std::vector<int>& inverse_outputs) {
    const int lanes = static_cast<int>(denominators.size());
    if (lanes <= 0 || numerators.size() != denominators.size()) {
        throw std::invalid_argument("invalid Montgomery lane lists");
    }

    const int n = config.n;
    const int block = config.block_size;
    const int pair_workspace = 2 * block;
    const int mul_size = config.mul_cost_size;
    int count = 0;




    const std::size_t inverse_begin = gm.current_pointer();
    std::vector<std::vector<int>> dirty_inverses;

    if (lanes == 1) {
        auto inv_main = slice_exact(ancilla, 0, pair_workspace, "inversion main workspace");
        auto inv_side = slice_exact(ancilla, ancilla.size() - pair_workspace,
                                    pair_workspace, "inversion side workspace");
        dirty_inverses.push_back(
            inversion_logic(denominators.front(), square_workspace,
                            inv_main, inv_side, count));
    } else {
        auto product = cons_mul(denominators, count, ancilla);
        count = product.count;

        auto inv_main = slice_exact(ancilla, 0, pair_workspace, "inversion main workspace");
        auto inv_side = slice_exact(ancilla, ancilla.size() - pair_workspace,
                                    pair_workspace, "inversion side workspace");
        auto root_inverse = inversion_logic(product.product, square_workspace,
                                            inv_main, inv_side, count);

        std::vector<std::vector<std::vector<int>>> inverse_tree;
        inverse_tree.push_back({root_inverse});

        for (int level = static_cast<int>(product.tree.size()) - 2; level >= 0; --level) {
            const auto& children = product.tree[static_cast<std::size_t>(level)];
            const auto& parent_inverses = inverse_tree.back();
            std::vector<std::vector<int>> next;
            next.reserve(children.size());

            std::size_t parent = 0;
            std::size_t copy_offset = 0;
            for (std::size_t i = 0; i < children.size(); i += 2, ++parent) {
                if (i + 1 == children.size()) {
                    next.push_back(parent_inverses[parent]);
                    continue;
                }

                const std::size_t copy_begin = gm.current_pointer();
                auto parent_copy = slice_exact(recovery_copy, copy_offset, n,
                                               "recovery-tree copy workspace");
                CNOT_N(parent_inverses[parent], parent_copy);
                const std::size_t copy_end = gm.current_pointer();
                copy_offset += static_cast<std::size_t>(n);

                const std::size_t left_base = i * static_cast<std::size_t>(pair_workspace);
                auto left_a = concat(parent_inverses[parent],
                                     slice_exact(ancilla, left_base, block,
                                                 "recovery-tree ancilla"));
                auto left_b = concat(children[i + 1],
                                     slice_exact(ancilla, left_base + block, block,
                                                 "recovery-tree ancilla"));
                auto left_target = slice_exact(Toffoli_qubits,
                                               static_cast<std::size_t>(count), mul_size,
                                               "recovery-tree Toffoli targets");
                auto left_inverse = mul(left_a, left_b, left_target, 0);
                count += mul_size;

                const std::size_t right_base = (i + 1) * static_cast<std::size_t>(pair_workspace);
                auto right_a = concat(parent_copy,
                                      slice_exact(ancilla, right_base, block,
                                                  "recovery-tree ancilla"));
                auto right_b = concat(children[i],
                                      slice_exact(ancilla, right_base + block, block,
                                                  "recovery-tree ancilla"));
                auto right_target = slice_exact(Toffoli_qubits,
                                                static_cast<std::size_t>(count), mul_size,
                                                "recovery-tree Toffoli targets");
                auto right_inverse = mul(right_a, right_b, right_target, 0);
                count += mul_size;

                gm.replay_reverse(copy_begin, copy_end);
                next.push_back(std::move(left_inverse));
                next.push_back(std::move(right_inverse));
            }
            inverse_tree.push_back(std::move(next));
        }
        dirty_inverses = inverse_tree.back();
    }

    const std::size_t inverse_end = gm.current_pointer();
    std::vector<std::vector<int>> clean_inverses;
    clean_inverses.reserve(static_cast<std::size_t>(lanes));
    for (int i = 0; i < lanes; ++i) {
        auto target = slice_exact(inverse_outputs, static_cast<std::size_t>(i) * n,
                                  n, "clean inverse outputs");
        CNOT_N(dirty_inverses[static_cast<std::size_t>(i)], target);
        clean_inverses.push_back(std::move(target));
    }
    gm.replay_reverse(inverse_begin, inverse_end);


    const std::size_t lambda_begin = gm.current_pointer();
    std::vector<std::vector<int>> dirty_lambdas;
    dirty_lambdas.reserve(static_cast<std::size_t>(lanes));
    for (int i = 0; i < lanes; ++i) {
        const std::size_t base = static_cast<std::size_t>(i) * pair_workspace;
        auto left = concat(numerators[static_cast<std::size_t>(i)],
                           slice_exact(ancilla, base, block, "lambda ancilla"));
        auto right = concat(clean_inverses[static_cast<std::size_t>(i)],
                            slice_exact(ancilla, base + block, block, "lambda ancilla"));
        auto targets = slice_exact(Toffoli_qubits, static_cast<std::size_t>(i) * mul_size,
                                   mul_size, "lambda Toffoli targets");
        dirty_lambdas.push_back(mul(left, right, targets, 0));
    }
    const std::size_t lambda_end = gm.current_pointer();

    std::vector<std::vector<int>> clean_lambdas;
    clean_lambdas.reserve(static_cast<std::size_t>(lanes));
    for (int i = 0; i < lanes; ++i) {
        auto target = slice_exact(lambda_outputs, static_cast<std::size_t>(i) * n,
                                  n, "clean lambda outputs");
        CNOT_N(dirty_lambdas[static_cast<std::size_t>(i)], target);
        clean_lambdas.push_back(std::move(target));
    }
    gm.replay_reverse(lambda_begin, lambda_end);

    return MontgomeryResult{std::move(clean_lambdas)};
}

void QuantumContext::accumulation_point_add(
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
    const std::vector<int>& cswap_controls) {
    const int n = config.n;
    const int block = config.block_size;
    const int pair_workspace = 2 * block;
    const int mul_size = config.mul_cost_size;
    const auto curve_a = constant_bits(n, 3);

    for (int i = 0; i < lanes; ++i) {
        CNOT_N(x1[static_cast<std::size_t>(i)], x1_plus_x2[static_cast<std::size_t>(i)]);
        CONST_ADD(y1[static_cast<std::size_t>(i)], y2_constants[static_cast<std::size_t>(i)]);
        CONST_ADD(x1_plus_x2[static_cast<std::size_t>(i)], x2_constants[static_cast<std::size_t>(i)]);
    }

    auto montgomery = montgomery_clean(x1_plus_x2, y1, ancilla, recovery_copy,
                                       square_workspace, lambda_outputs, inverse_outputs);

    for (int i = 0; i < lanes; ++i) {
        const std::size_t lane = static_cast<std::size_t>(i);
        CONST_ADD(y1[lane], y2_constants[lane]);
        CONST_ADD(x1_plus_x2[lane], x2_constants[lane]);
        CONST_ADD(x1_plus_x2[lane], curve_a);

        auto transformed = Squaring_plus_input(
            concat(montgomery.lambdas[lane], x1_plus_x2[lane]), 2 * n);
        montgomery.lambdas[lane] = std::move(transformed.first);
        x1_plus_x2[lane] = std::move(transformed.second);

        const std::size_t mul_begin = gm.current_pointer();
        const std::size_t base = lane * static_cast<std::size_t>(pair_workspace);
        auto left = concat(x1_plus_x2[lane],
                           slice_exact(ancilla, base, block, "coordinate-multiply ancilla"));
        auto right = concat(montgomery.lambdas[lane],
                            slice_exact(ancilla, base + block, block,
                                        "coordinate-multiply ancilla"));
        auto targets = slice_exact(Toffoli_qubits, lane * static_cast<std::size_t>(mul_size),
                                   mul_size, "coordinate-multiply Toffoli targets");
        auto dirty_y = mul(left, right, targets, 0);
        const std::size_t mul_end = gm.current_pointer();

        auto clean_y = slice_exact(multiply_outputs, lane * static_cast<std::size_t>(n),
                                   n, "coordinate-multiply outputs");
        CNOT_N(dirty_y, clean_y);
        gm.replay_reverse(mul_begin, mul_end);

        CONST_ADD(x1_plus_x2[lane], x2_constants[lane]);
        CONST_ADD(clean_y, y2_constants[lane]);
        CNOT_N(x1_plus_x2[lane], clean_y);

        auto fanout = slice_exact(cswap_controls,
                                  lane * static_cast<std::size_t>(2 * n - 1),
                                  2 * n - 1, "CSWAP fanout workspace");
        const std::size_t fanout_begin = gm.current_pointer();
        auto control_copies = copy_parallel(controls[lane], fanout, 2 * n);
        const std::size_t fanout_end = gm.current_pointer();
        for (int bit = 0; bit < n; ++bit) {
            CSWAP(control_copies[static_cast<std::size_t>(bit)],
                  x1[lane][static_cast<std::size_t>(bit)],
                  x1_plus_x2[lane][static_cast<std::size_t>(bit)]);
            CSWAP(control_copies[static_cast<std::size_t>(n + bit)],
                  y1[lane][static_cast<std::size_t>(bit)],
                  clean_y[static_cast<std::size_t>(bit)]);
        }
        gm.replay_reverse(fanout_begin, fanout_end);
    }
}

void QuantumContext::reduction_point_add(
    int lanes,
    std::vector<std::vector<int>> x1, std::vector<std::vector<int>> y1,
    std::vector<std::vector<int>> x2, std::vector<std::vector<int>> y2,
    const std::vector<int>& ancilla,
    std::vector<std::vector<int>> x1_plus_x2,
    const std::vector<int>& recovery_copy,
    const std::vector<int>& square_workspace,
    const std::vector<int>& lambda_outputs,
    const std::vector<int>& multiply_outputs,
    const std::vector<int>& inverse_outputs) {
    const int n = config.n;
    const int block = config.block_size;
    const int pair_workspace = 2 * block;
    const int mul_size = config.mul_cost_size;
    const auto curve_a = constant_bits(n, 3);

    for (int i = 0; i < lanes; ++i) {
        const std::size_t lane = static_cast<std::size_t>(i);
        CNOT_N(x1[lane], x1_plus_x2[lane]);
        CNOT_N(x2[lane], x1_plus_x2[lane]);
        CNOT_N(y2[lane], y1[lane]);
    }

    auto montgomery = montgomery_clean(x1_plus_x2, y1, ancilla, recovery_copy,
                                       square_workspace, lambda_outputs, inverse_outputs);

    for (int i = 0; i < lanes; ++i) {
        const std::size_t lane = static_cast<std::size_t>(i);
        CONST_ADD(x1_plus_x2[lane], curve_a);
        auto transformed = Squaring_plus_input(
            concat(montgomery.lambdas[lane], x1_plus_x2[lane]), 2 * n);
        montgomery.lambdas[lane] = std::move(transformed.first);
        x1_plus_x2[lane] = std::move(transformed.second);
        CNOT_N(x2[lane], x1_plus_x2[lane]);

        const std::size_t mul_begin = gm.current_pointer();
        const std::size_t base = lane * static_cast<std::size_t>(pair_workspace);
        auto left = concat(x1_plus_x2[lane],
                           slice_exact(ancilla, base, block, "coordinate-multiply ancilla"));
        auto right = concat(montgomery.lambdas[lane],
                            slice_exact(ancilla, base + block, block,
                                        "coordinate-multiply ancilla"));
        auto targets = slice_exact(Toffoli_qubits, lane * static_cast<std::size_t>(mul_size),
                                   mul_size, "coordinate-multiply Toffoli targets");
        auto dirty_y = mul(left, right, targets, 0);
        const std::size_t mul_end = gm.current_pointer();

        auto clean_y = slice_exact(multiply_outputs, lane * static_cast<std::size_t>(n),
                                   n, "coordinate-multiply outputs");
        CNOT_N(dirty_y, clean_y);
        gm.replay_reverse(mul_begin, mul_end);

        CNOT_N(y2[lane], clean_y);
        CNOT_N(x2[lane], x1_plus_x2[lane]);
        CNOT_N(x1_plus_x2[lane], clean_y);
    }
}

LayerStats QuantumContext::run_balanced_layer(int lanes, PointAddMode mode) {
    if (lanes <= 0) {
        throw std::invalid_argument("point-add layer must contain at least one lane");
    }
    gm.clear();

    const uint64_t k = static_cast<uint64_t>(lanes);
    const uint64_t n = static_cast<uint64_t>(config.n);
    const uint64_t block = static_cast<uint64_t>(config.block_size);
    const uint64_t mul_size = static_cast<uint64_t>(config.mul_cost_size);
    const uint64_t tree_and_inverse_blocks =
        3 * (k - 1) + static_cast<uint64_t>(config.inversion_mul_blocks);

    uint64_t cursor = 0;
    std::vector<int> controls;
    if (mode == PointAddMode::Accumulation) controls = allocate_range(cursor, k);
    auto x1_flat = allocate_range(cursor, n * k);
    std::vector<int> x2_flat;
    if (mode == PointAddMode::Reduction) x2_flat = allocate_range(cursor, n * k);
    auto x_sum_flat = allocate_range(cursor, n * k);
    auto y1_flat = allocate_range(cursor, n * k);
    std::vector<int> y2_flat;
    if (mode == PointAddMode::Reduction) y2_flat = allocate_range(cursor, n * k);

    auto ancilla = allocate_range(cursor, 4 * block * k);
    Toffoli_qubits = allocate_range(cursor, mul_size * tree_and_inverse_blocks);
    auto square_workspace = allocate_range(cursor, 4 * n);
    auto recovery_copy = allocate_range(cursor, n * (k / 2));
    auto lambda_outputs = allocate_range(cursor, n * k);
    auto multiply_outputs = allocate_range(cursor, n * k);
    auto inverse_outputs = allocate_range(cursor, n * k);
    std::vector<int> cswap_controls;
    if (mode == PointAddMode::Accumulation) {
        cswap_controls = allocate_range(cursor, (2 * n - 1) * k);
    }

    auto x1 = split_register(x1_flat, config.n, lanes, "x1");
    auto y1 = split_register(y1_flat, config.n, lanes, "y1");
    auto x_sum = split_register(x_sum_flat, config.n, lanes, "x1+x2");

    if (mode == PointAddMode::Accumulation) {
        const std::vector<bool> all_ones(static_cast<std::size_t>(config.n), true);
        std::vector<std::vector<bool>> x2_constants(static_cast<std::size_t>(lanes), all_ones);
        std::vector<std::vector<bool>> y2_constants(static_cast<std::size_t>(lanes), all_ones);
        accumulation_point_add(lanes, controls, std::move(x1), std::move(y1),
                               x2_constants, y2_constants, ancilla, std::move(x_sum),
                               recovery_copy, square_workspace, lambda_outputs,
                               multiply_outputs, inverse_outputs, cswap_controls);
    } else {
        auto x2 = split_register(x2_flat, config.n, lanes, "x2");
        auto y2 = split_register(y2_flat, config.n, lanes, "y2");
        reduction_point_add(lanes, std::move(x1), std::move(y1), std::move(x2),
                            std::move(y2), ancilla, std::move(x_sum), recovery_copy,
                            square_workspace, lambda_outputs, multiply_outputs,
                            inverse_outputs);
    }

    auto stats = unpack_stats(cursor, gm.optimize_and_count(static_cast<uint32_t>(cursor)));
    gm.clear();
    return stats;
}

PrimitiveStats QuantumContext::estimate_multiplication() {
    gm.clear();
    uint64_t cursor = 0;
    const uint64_t l = static_cast<uint64_t>(config.mul_cost_size);
    auto left = allocate_range(cursor, l);
    auto right = allocate_range(cursor, l);
    auto target = allocate_range(cursor, l);
    Toffoli_qubits = target;
    (void)mul(left, right, target, 0);
    auto stats = unpack_stats(cursor, gm.optimize_and_count(static_cast<uint32_t>(cursor)));
    gm.clear();
    return stats;
}

PrimitiveStats QuantumContext::estimate_inversion() {
    gm.clear();
    uint64_t cursor = 0;
    const uint64_t n = static_cast<uint64_t>(config.n);
    const uint64_t block = static_cast<uint64_t>(config.block_size);
    const uint64_t l = static_cast<uint64_t>(config.mul_cost_size);
    auto input = allocate_range(cursor, n);
    auto ancilla = allocate_range(cursor, 4 * block);
    Toffoli_qubits = allocate_range(
        cursor, l * static_cast<uint64_t>(config.inversion_mul_blocks));
    auto square_workspace = allocate_range(cursor, 4 * n);
    int count = 0;
    auto main = slice_exact(ancilla, 0, 2 * block, "inversion main workspace");
    auto side = slice_exact(ancilla, 2 * block, 2 * block, "inversion side workspace");
    (void)inversion_logic(input, square_workspace, main, side, count);
    auto stats = unpack_stats(cursor, gm.optimize_and_count(static_cast<uint32_t>(cursor)));
    gm.clear();
    return stats;
}
