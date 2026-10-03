#include "../BasicLogic/QuantumLib.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>

namespace {

template <typename T>
std::vector<T> concat_window(const std::vector<T>& left, const std::vector<T>& right) {
    std::vector<T> out;
    out.reserve(left.size() + right.size());
    out.insert(out.end(), left.begin(), left.end());
    out.insert(out.end(), right.begin(), right.end());
    return out;
}

template <typename T>
std::vector<T> slice_window(const std::vector<T>& values, std::size_t start,
                            std::size_t length, const char* label) {
    if (start > values.size() || length > values.size() - start) {
        throw std::out_of_range(std::string("slice outside ") + label);
    }
    return std::vector<T>(values.begin() + static_cast<std::ptrdiff_t>(start),
                          values.begin() + static_cast<std::ptrdiff_t>(start + length));
}

std::vector<int> allocate_window(uint64_t& cursor, uint64_t length) {
    if (cursor + length > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("qubit index exceeds signed-int range");
    }
    std::vector<int> out;
    out.reserve(static_cast<std::size_t>(length));
    for (uint64_t i = 0; i < length; ++i) out.push_back(static_cast<int>(cursor + i));
    cursor += length;
    return out;
}

std::vector<std::vector<int>> split_window(const std::vector<int>& flat, int width,
                                            int count, const char* label) {
    if (width < 0 || count < 0 ||
        flat.size() != static_cast<std::size_t>(width) * static_cast<std::size_t>(count)) {
        throw std::invalid_argument(std::string("invalid register layout: ") + label);
    }
    std::vector<std::vector<int>> out;
    out.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        out.push_back(slice_window(flat, static_cast<std::size_t>(i) * width,
                                   static_cast<std::size_t>(width), label));
    }
    return out;
}

} // namespace

std::size_t PackedQromTable::words_per_entry() const {
    if (word_bits <= 0) throw std::invalid_argument("QROM word_bits must be positive");
    return (static_cast<std::size_t>(word_bits) + 63U) / 64U;
}

std::size_t PackedQromTable::entries() const {
    if (address_bits <= 0 || address_bits >= 63) {
        throw std::invalid_argument("QROM address_bits outside supported range");
    }
    const std::size_t expected = std::size_t{1} << address_bits;
    if (words.size() != expected * words_per_entry()) {
        throw std::invalid_argument("packed QROM length does not match address_bits/word_bits");
    }
    return expected;
}

bool PackedQromTable::bit(std::size_t entry, int bit_index) const {
    if (bit_index < 0 || bit_index >= word_bits || entry >= entries()) {
        throw std::out_of_range("packed QROM bit outside table");
    }
    const std::size_t offset = entry * words_per_entry() +
                               static_cast<std::size_t>(bit_index) / 64U;
    return ((words[offset] >> (bit_index % 64)) & 1U) != 0;
}

void QuantumContext::coherent_unary_qrom(
    const std::vector<int>& address,
    const std::vector<int>& selector_ancillas,
    const std::vector<int>& data_fanout,
    const std::vector<int>& bus,
    const PackedQromTable& table,
    bool inverse) {
    const std::size_t entry_count = table.entries();
    if (address.size() != static_cast<std::size_t>(table.address_bits) ||
        selector_ancillas.size() != static_cast<std::size_t>(table.address_bits - 1) ||
        bus.size() != static_cast<std::size_t>(table.word_bits) ||
        data_fanout.size() + 1 < bus.size()) {
        throw std::invalid_argument("QROM register dimensions do not match table");
    }
    std::vector<int> selection(address.rbegin(), address.rend());

    auto emit_word = [&](int control, std::size_t entry, bool reverse) {
        std::vector<int> set_bits;
        set_bits.reserve(static_cast<std::size_t>(table.word_bits));
        for (int bit_index = 0; bit_index < table.word_bits; ++bit_index) {
            if (table.bit(entry, bit_index)) set_bits.push_back(bit_index);
        }
        if (set_bits.empty()) return;
        const std::size_t fanout_begin = gm.current_pointer();
        const auto controls = copy_parallel(
            control, slice_window(data_fanout, 0, set_bits.size() - 1, "QROM data fanout"),
            static_cast<int>(set_bits.size()));
        const std::size_t fanout_end = gm.current_pointer();
        if (!reverse) {
            for (std::size_t i = 0; i < set_bits.size(); ++i) {
                CNOT(controls[i], bus[static_cast<std::size_t>(set_bits[i])]);
            }
        } else {
            for (std::size_t i = set_bits.size(); i-- > 0;) {
                CNOT(controls[i], bus[static_cast<std::size_t>(set_bits[i])]);
            }
        }
        gm.replay_reverse(fanout_begin, fanout_end);
    };

    auto negative_and = [&](int control, int select, int target) {
        X(select);
        Toffoli(control, select, target);
        X(select);
    };

    std::function<void(int, int, std::size_t, std::size_t, bool)> subtree;
    subtree = [&](int control, int level, std::size_t left, std::size_t right, bool reverse) {
        if (right - left == 1) {
            emit_word(control, left, reverse);
            return;
        }
        const std::size_t middle = (left + right) / 2;
        const int select = selection[static_cast<std::size_t>(level + 1)];
        const int ancilla = selector_ancillas[static_cast<std::size_t>(level)];
        if (!reverse) {
            negative_and(control, select, ancilla);
            subtree(ancilla, level + 1, left, middle, false);
            CNOT(control, ancilla);
            subtree(ancilla, level + 1, middle, right, false);
            Toffoli(control, select, ancilla);
        } else {
            Toffoli(control, select, ancilla);
            subtree(ancilla, level + 1, middle, right, true);
            CNOT(control, ancilla);
            subtree(ancilla, level + 1, left, middle, true);
            negative_and(control, select, ancilla);
        }
    };

    const int root = selection.front();
    const std::size_t middle = entry_count / 2;
    if (!inverse) {
        X(root);
        subtree(root, 0, 0, middle, false);
        X(root);
        subtree(root, 0, middle, entry_count, false);
    } else {
        subtree(root, 0, middle, entry_count, true);
        X(root);
        subtree(root, 0, 0, middle, true);
        X(root);
    }
}

void QuantumContext::early_unlookup_point_add(
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
    const std::vector<bool>& curve_a) {
    if (lanes <= 0 || tables.size() != static_cast<std::size_t>(lanes) ||
        curve_a.size() != static_cast<std::size_t>(config.n)) {
        throw std::invalid_argument("invalid early-unlookup point-add arguments");
    }
    // All lane QROMs are emitted into one gate stream.  Since their registers
    // are disjoint, the per-line ASAP pass schedules them in parallel.
    for (int i = 0; i < lanes; ++i) {
        const std::size_t lane = static_cast<std::size_t>(i);
        coherent_unary_qrom(addresses[lane], selectors[lane], data_fanout[lane],
                            concat_window(x_table[lane], y_table[lane]), tables[lane], false);
    }
    for (int i = 0; i < lanes; ++i) {
        const std::size_t lane = static_cast<std::size_t>(i);
        CNOT_N(x1[lane], z[lane]);
        CNOT_N(x_table[lane], z[lane]);
        CNOT_N(y1[lane], v[lane]);
        CNOT_N(y_table[lane], v[lane]);
    }
    for (int i = lanes - 1; i >= 0; --i) {
        const std::size_t lane = static_cast<std::size_t>(i);
        coherent_unary_qrom(addresses[lane], selectors[lane], data_fanout[lane],
                            concat_window(x_table[lane], y_table[lane]), tables[lane], true);
    }

    early_point_add_tail(lanes, std::move(x1), std::move(y1), std::move(z), std::move(v),
                         ancilla, recovery_copy, square_workspace, lambda_outputs,
                         multiply_outputs, inverse_outputs, curve_a);
}

void QuantumContext::early_point_add_tail(
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
    const std::vector<bool>& curve_a) {
    if (lanes <= 0 || curve_a.size() != static_cast<std::size_t>(config.n) ||
        x1.size() != static_cast<std::size_t>(lanes) ||
        y1.size() != static_cast<std::size_t>(lanes) ||
        z.size() != static_cast<std::size_t>(lanes) ||
        v.size() != static_cast<std::size_t>(lanes)) {
        throw std::invalid_argument("invalid early point-add tail arguments");
    }
    const int n = config.n;
    const int block = config.block_size;
    const int pair_workspace = 2 * block;
    const int mul_size = config.mul_cost_size;

    auto montgomery = montgomery_clean(z, v, ancilla, recovery_copy,
                                       square_workspace, lambda_outputs, inverse_outputs);

    for (int i = 0; i < lanes; ++i) {
        const std::size_t lane = static_cast<std::size_t>(i);
        CONST_ADD(z[lane], curve_a);
        auto transformed = Squaring_plus_input(
            concat_window(montgomery.lambdas[lane], z[lane]), 2 * n);
        montgomery.lambdas[lane] = std::move(transformed.first);
        z[lane] = std::move(transformed.second); // x3

        // y3 = lambda*(x1+x3) + x3 + y1.  This form deliberately
        // avoids every post-unlookup reference to x_table/y_table.
        CNOT_N(x1[lane], z[lane]);
        const std::size_t multiply_begin = gm.current_pointer();
        const std::size_t base = lane * static_cast<std::size_t>(pair_workspace);
        auto left = concat_window(z[lane],
                                  slice_window(ancilla, base, block, "window coordinate ancilla"));
        auto right = concat_window(
            montgomery.lambdas[lane],
            slice_window(ancilla, base + block, block, "window coordinate ancilla"));
        auto targets = slice_window(Toffoli_qubits,
                                    lane * static_cast<std::size_t>(mul_size), mul_size,
                                    "window coordinate Toffoli targets");
        auto dirty_y = mul(left, right, targets, 0);
        const std::size_t multiply_end = gm.current_pointer();

        auto clean_y = slice_window(multiply_outputs,
                                    lane * static_cast<std::size_t>(n), n,
                                    "window coordinate outputs");
        CNOT_N(dirty_y, clean_y);
        gm.replay_reverse(multiply_begin, multiply_end);
        CNOT_N(x1[lane], z[lane]); // restore x3
        CNOT_N(z[lane], clean_y);
        CNOT_N(y1[lane], clean_y);
    }
}

LayerStats QuantumContext::run_qrom_pair(const PackedQromTable& table) {
    if (table.word_bits != 2 * config.n) {
        throw std::invalid_argument("QROM pair word width must be 2n");
    }
    gm.clear();
    const uint64_t n = static_cast<uint64_t>(config.n);
    uint64_t cursor = 0;
    auto address = allocate_window(cursor, static_cast<uint64_t>(table.address_bits));
    auto x1 = allocate_window(cursor, n);
    auto y1 = allocate_window(cursor, n);
    auto z = allocate_window(cursor, n);
    auto v = allocate_window(cursor, n);
    auto x_table = allocate_window(cursor, n);
    auto y_table = allocate_window(cursor, n);
    auto selectors = allocate_window(cursor, static_cast<uint64_t>(table.address_bits - 1));
    auto fanout = allocate_window(cursor, 2 * n - 1);
    auto bus = concat_window(x_table, y_table);

    coherent_unary_qrom(address, selectors, fanout, bus, table, false);
    CNOT_N(x1, z);
    CNOT_N(x_table, z);
    CNOT_N(y1, v);
    CNOT_N(y_table, v);
    coherent_unary_qrom(address, selectors, fanout, bus, table, true);

    const auto stats = gm.optimize_and_count(static_cast<uint32_t>(cursor));
    LayerStats result{cursor, std::get<3>(stats), std::get<4>(stats), std::get<0>(stats),
                      std::get<1>(stats), std::get<2>(stats)};
    gm.clear();
    return result;
}

LayerStats QuantumContext::run_early_tail(int lanes, const std::vector<bool>& curve_a) {
    if (lanes <= 0 || curve_a.size() != static_cast<std::size_t>(config.n)) {
        throw std::invalid_argument("invalid early-tail arguments");
    }
    gm.clear();
    const uint64_t k = static_cast<uint64_t>(lanes);
    const uint64_t n = static_cast<uint64_t>(config.n);
    const uint64_t block = static_cast<uint64_t>(config.block_size);
    const uint64_t mul_size = static_cast<uint64_t>(config.mul_cost_size);
    const uint64_t arithmetic_blocks =
        3 * (k - 1) + static_cast<uint64_t>(config.inversion_mul_blocks);
    uint64_t cursor = 0;

    auto x1_flat = allocate_window(cursor, n * k);
    auto y1_flat = allocate_window(cursor, n * k);
    auto z_flat = allocate_window(cursor, n * k);
    auto v_flat = allocate_window(cursor, n * k);
    auto ancilla = allocate_window(cursor, 4 * block * k);
    Toffoli_qubits = allocate_window(cursor, mul_size * arithmetic_blocks);
    auto square_workspace = allocate_window(cursor, 4 * n);
    auto recovery_copy = allocate_window(cursor, n * (k / 2));
    auto lambda_outputs = allocate_window(cursor, n * k);
    auto multiply_outputs = allocate_window(cursor, n * k);
    auto inverse_outputs = allocate_window(cursor, n * k);

    auto x1 = split_window(x1_flat, config.n, lanes, "early-tail x1");
    auto y1 = split_window(y1_flat, config.n, lanes, "early-tail y1");
    auto z = split_window(z_flat, config.n, lanes, "early-tail z");
    auto v = split_window(v_flat, config.n, lanes, "early-tail v");
    early_point_add_tail(lanes, std::move(x1), std::move(y1), std::move(z), std::move(v),
                         ancilla, recovery_copy, square_workspace, lambda_outputs,
                         multiply_outputs, inverse_outputs, curve_a);

    const auto stats = gm.optimize_and_count(static_cast<uint32_t>(cursor));
    LayerStats result{cursor, std::get<3>(stats), std::get<4>(stats), std::get<0>(stats),
                      std::get<1>(stats), std::get<2>(stats)};
    gm.clear();
    return result;
}

LayerStats QuantumContext::run_window_layer(const std::vector<PackedQromTable>& tables,
                                             const std::vector<bool>& curve_a) {
    const int lanes = static_cast<int>(tables.size());
    if (lanes <= 0) throw std::invalid_argument("window layer needs at least one table");
    gm.clear();

    const uint64_t k = static_cast<uint64_t>(lanes);
    const uint64_t n = static_cast<uint64_t>(config.n);
    const uint64_t block = static_cast<uint64_t>(config.block_size);
    const uint64_t mul_size = static_cast<uint64_t>(config.mul_cost_size);
    const uint64_t arithmetic_blocks =
        3 * (k - 1) + static_cast<uint64_t>(config.inversion_mul_blocks);

    uint64_t cursor = 0;
    std::vector<std::vector<int>> addresses;
    addresses.reserve(tables.size());
    for (const auto& table : tables) {
        if (table.word_bits != 2 * config.n) {
            throw std::invalid_argument("window table word width must be 2n");
        }
        addresses.push_back(allocate_window(cursor, static_cast<uint64_t>(table.address_bits)));
    }
    auto x1_flat = allocate_window(cursor, n * k);
    auto y1_flat = allocate_window(cursor, n * k);
    auto z_flat = allocate_window(cursor, n * k);
    auto v_flat = allocate_window(cursor, n * k);

    uint64_t selector_bits = 0;
    for (const auto& table : tables) selector_bits += static_cast<uint64_t>(table.address_bits - 1);
    const uint64_t data_fanout_bits = (2 * n - 1) * k;
    const uint64_t qrom_scratch = 2 * n * k + selector_bits + data_fanout_bits;
    const uint64_t arithmetic_scratch =
        4 * block * k + mul_size * arithmetic_blocks + 4 * n + n * (k / 2) + 3 * n * k;
    auto reusable = allocate_window(cursor, std::max(qrom_scratch, arithmetic_scratch));

    auto x_table_flat = slice_window(reusable, 0, n * k, "QROM x bus");
    auto y_table_flat = slice_window(reusable, n * k, n * k, "QROM y bus");
    std::vector<std::vector<int>> selectors;
    selectors.reserve(tables.size());
    std::size_t selector_cursor = static_cast<std::size_t>(2 * n * k);
    for (const auto& table : tables) {
        const std::size_t count = static_cast<std::size_t>(table.address_bits - 1);
        selectors.push_back(slice_window(reusable, selector_cursor, count, "QROM selectors"));
        selector_cursor += count;
    }
    std::vector<std::vector<int>> data_fanout;
    data_fanout.reserve(tables.size());
    std::size_t fanout_cursor = static_cast<std::size_t>(2 * n * k + selector_bits);
    for (int i = 0; i < lanes; ++i) {
        data_fanout.push_back(slice_window(reusable, fanout_cursor,
                                           static_cast<std::size_t>(2 * n - 1),
                                           "QROM data fanout"));
        fanout_cursor += static_cast<std::size_t>(2 * n - 1);
    }

    std::size_t work = 0;
    auto ancilla = slice_window(reusable, work, 4 * block * k, "Montgomery ancilla");
    work += static_cast<std::size_t>(4 * block * k);
    Toffoli_qubits = slice_window(reusable, work, mul_size * arithmetic_blocks,
                                  "Montgomery Toffoli targets");
    work += static_cast<std::size_t>(mul_size * arithmetic_blocks);
    auto square_workspace = slice_window(reusable, work, 4 * n, "square workspace");
    work += static_cast<std::size_t>(4 * n);
    auto recovery_copy = slice_window(reusable, work, n * (k / 2), "recovery copy");
    work += static_cast<std::size_t>(n * (k / 2));
    auto lambda_outputs = slice_window(reusable, work, n * k, "lambda outputs");
    work += static_cast<std::size_t>(n * k);
    auto multiply_outputs = slice_window(reusable, work, n * k, "multiply outputs");
    work += static_cast<std::size_t>(n * k);
    auto inverse_outputs = slice_window(reusable, work, n * k, "inverse outputs");

    auto x1 = split_window(x1_flat, config.n, lanes, "window x1");
    auto y1 = split_window(y1_flat, config.n, lanes, "window y1");
    auto z = split_window(z_flat, config.n, lanes, "window z");
    auto v = split_window(v_flat, config.n, lanes, "window v");
    auto x_table = split_window(x_table_flat, config.n, lanes, "window x table");
    auto y_table = split_window(y_table_flat, config.n, lanes, "window y table");

    early_unlookup_point_add(lanes, std::move(x1), std::move(y1), x_table, y_table,
                             addresses, selectors, data_fanout, tables, ancilla,
                             std::move(z), std::move(v),
                             recovery_copy, square_workspace, lambda_outputs, multiply_outputs,
                             inverse_outputs, curve_a);

    const auto stats = gm.optimize_and_count(static_cast<uint32_t>(cursor));
    LayerStats result{cursor, std::get<3>(stats), std::get<4>(stats), std::get<0>(stats),
                      std::get<1>(stats), std::get<2>(stats)};
    gm.clear();
    return result;
}

LayerStats QuantumContext::run_window_initialization(
    const std::vector<PackedQromTable>& tables) {
    const int lanes = static_cast<int>(tables.size());
    if (lanes <= 0) throw std::invalid_argument("window initialization needs at least one table");
    gm.clear();
    const uint64_t n = static_cast<uint64_t>(config.n);
    uint64_t cursor = 0;

    std::vector<std::vector<int>> addresses;
    std::vector<std::vector<int>> accumulators;
    std::vector<std::vector<int>> selectors;
    std::vector<std::vector<int>> fanouts;
    addresses.reserve(tables.size());
    accumulators.reserve(tables.size());
    selectors.reserve(tables.size());
    fanouts.reserve(tables.size());
    for (const auto& table : tables) {
        if (table.word_bits != 2 * config.n) {
            throw std::invalid_argument("initialization table word width must be 2n");
        }
        addresses.push_back(allocate_window(cursor, static_cast<uint64_t>(table.address_bits)));
        accumulators.push_back(allocate_window(cursor, 2 * n));
        selectors.push_back(
            allocate_window(cursor, static_cast<uint64_t>(table.address_bits - 1)));
        fanouts.push_back(allocate_window(cursor, 2 * n - 1));
    }
    for (int i = 0; i < lanes; ++i) {
        const auto lane = static_cast<std::size_t>(i);
        coherent_unary_qrom(addresses[lane], selectors[lane], fanouts[lane],
                            accumulators[lane], tables[lane], false);
    }
    const auto stats = gm.optimize_and_count(static_cast<uint32_t>(cursor));
    LayerStats result{cursor, std::get<3>(stats), std::get<4>(stats), std::get<0>(stats),
                      std::get<1>(stats), std::get<2>(stats)};
    gm.clear();
    return result;
}
