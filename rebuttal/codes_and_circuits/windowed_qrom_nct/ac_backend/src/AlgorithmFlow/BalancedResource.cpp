#include "BalancedResource.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

uint64_t ceil_log2(uint64_t value) {
    if (value <= 1) return 0;
    uint64_t result = 0;
    uint64_t power = 1;
    while (power < value) {
        if (power > std::numeric_limits<uint64_t>::max() / 2) {
            throw std::overflow_error("ceil_log2 overflow");
        }
        power <<= 1;
        ++result;
    }
    return result;
}

uint64_t checked_add(uint64_t a, uint64_t b, const char* label) {
    if (b > std::numeric_limits<uint64_t>::max() - a) {
        throw std::overflow_error(label);
    }
    return a + b;
}

uint64_t checked_twice_plus(uint64_t value, uint64_t extra, const char* label) {
    if (value > (std::numeric_limits<uint64_t>::max() - extra) / 2) {
        throw std::overflow_error(label);
    }
    return 2 * value + extra;
}

int64_t as_signed(uint64_t value, const char* label) {
    if (value > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
        throw std::overflow_error(label);
    }
    return static_cast<int64_t>(value);
}

} // namespace

std::vector<ScheduledLayer> build_balanced_schedule(int total_inputs, int w) {
    if (total_inputs <= 0 || w <= 0 || 2 * w > total_inputs) {
        throw std::invalid_argument("Balanced schedule requires 1 <= w <= total_inputs/2");
    }

    std::vector<ScheduledLayer> schedule;


    schedule.push_back({PointAddMode::Accumulation, w});
    int remaining = total_inputs - 2 * w;
    while (remaining > 0) {
        const int lanes = std::min(w, remaining);
        schedule.push_back({PointAddMode::Accumulation, lanes});
        remaining -= lanes;
    }


    int items = w;
    while (items > 1) {
        const int additions = items / 2;
        schedule.push_back({PointAddMode::Reduction, additions});
        items = additions + (items % 2);
    }
    return schedule;
}

BalancedResult estimate_balanced(const GFConfig& config, int w,
                                 const LayerProvider& layer_provider) {
    const int total_inputs = config.total_shor_inputs();
    const auto schedule = build_balanced_schedule(total_inputs, w);
    const int64_t n = config.n;

    uint64_t total_toffoli = 0;
    uint64_t total_cnot = 0;
    uint64_t total_full_depth = 0;
    uint64_t total_current_depth = 0;
    uint64_t total_toffoli_depth = 0;
    int64_t peak_qubits = 0;
    int64_t clear_qubits = 0;
    int accumulation_layers = 0;
    int reduction_layers = 0;

    for (std::size_t layer_index = 0; layer_index < schedule.size(); ++layer_index) {
        const auto layer = schedule[layer_index];
        const auto stats = layer_provider(layer.lanes, layer.mode);
        const int64_t length = as_signed(stats.qubits, "layer width exceeds int64");
        const int64_t k = layer.lanes;

        total_toffoli = checked_add(total_toffoli, stats.toffoli, "Toffoli count overflow");
        total_toffoli_depth = checked_add(total_toffoli_depth, stats.toffoli_depth,
                                          "Toffoli depth overflow");

        if (layer_index == 0) {
            if (layer.mode != PointAddMode::Accumulation || layer.lanes != w) {
                throw std::logic_error("invalid first Balanced layer");
            }
            const uint64_t fanout_depth = 2 * ceil_log2(static_cast<uint64_t>(2 * n));
            total_full_depth = checked_add(total_full_depth,
                                           fanout_depth + 1 + stats.full_depth,
                                           "full depth overflow");
            total_current_depth = checked_add(total_current_depth,
                                              fanout_depth + 1 + stats.current_depth,
                                              "current depth overflow");
            const uint64_t initialization_cnot =
                static_cast<uint64_t>(2 * k * (2 * n - 1) + 2 * n * k);
            total_cnot = checked_add(total_cnot, initialization_cnot + stats.cnot,
                                     "CNOT count overflow");

            peak_qubits = length - k;
            clear_qubits = length - 6 * n * k - k;
            ++accumulation_layers;
            continue;
        }

        total_full_depth = checked_add(total_full_depth, stats.full_depth,
                                       "full depth overflow");
        total_current_depth = checked_add(total_current_depth, stats.current_depth,
                                          "current depth overflow");
        total_cnot = checked_add(total_cnot, stats.cnot, "CNOT count overflow");

        if (layer.mode == PointAddMode::Accumulation) {
            const int64_t demand = length - 2 * n * k - k;
            if (demand > clear_qubits) {
                peak_qubits += demand - clear_qubits;
                clear_qubits = length - 6 * n * k - k;
            } else {
                const int64_t over = clear_qubits - demand;
                clear_qubits = length - 6 * n * k - k + over;
            }
            ++accumulation_layers;
        } else {
            const int64_t demand = length - 4 * n * k;
            if (demand > clear_qubits) {
                peak_qubits += demand - clear_qubits;
                clear_qubits = length - 8 * n * k;
            } else {
                const int64_t over = clear_qubits - demand;
                clear_qubits = length - 8 * n * k + over;
            }
            ++reduction_layers;
        }
    }

    if (peak_qubits < 0) throw std::logic_error("negative peak width");
    const uint64_t final_qubits = static_cast<uint64_t>(peak_qubits + 4 * n + 2);
    const uint64_t final_toffoli = checked_twice_plus(total_toffoli, 0, "Toffoli overflow");
    const uint64_t final_cnot = checked_twice_plus(total_cnot, 2 * n, "CNOT overflow");
    const uint64_t final_full_depth =
        checked_twice_plus(total_full_depth, 2, "full depth overflow");
    const uint64_t final_current_depth =
        checked_twice_plus(total_current_depth, 2, "current depth overflow");
    const uint64_t final_toffoli_depth =
        checked_twice_plus(total_toffoli_depth, 0, "Toffoli depth overflow");

    return BalancedResult{
        config.n,
        total_inputs,
        w,
        accumulation_layers,
        reduction_layers,
        final_toffoli,
        final_cnot,
        final_full_depth,
        final_current_depth,
        final_toffoli_depth,
        final_qubits,
        static_cast<uint128_t>(final_full_depth) * final_qubits,
        static_cast<uint128_t>(final_current_depth) * final_qubits,
        static_cast<uint128_t>(final_toffoli_depth) * final_qubits};
}

std::string uint128_to_string(uint128_t value) {
    if (value == 0) return "0";
    std::string out;
    while (value != 0) {
        const unsigned digit = static_cast<unsigned>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

long double uint128_log2(uint128_t value) {
    if (value == 0) return -std::numeric_limits<long double>::infinity();
    int exponent = 0;
    uint128_t shifted = value;
    while (shifted > static_cast<uint128_t>(std::numeric_limits<uint64_t>::max())) {
        shifted >>= 1;
        ++exponent;
    }
    return std::log2(static_cast<long double>(static_cast<uint64_t>(shifted))) + exponent;
}
