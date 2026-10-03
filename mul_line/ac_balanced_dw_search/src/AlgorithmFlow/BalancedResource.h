#pragma once

#include "../BasicLogic/QuantumLib.h"
#include "../Infrastructure/config.h"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

using uint128_t = unsigned __int128;

struct ScheduledLayer {
    PointAddMode mode{};
    int lanes{};
};

struct BalancedResult {
    int n{};
    int total_inputs{};
    int w{};
    int accumulation_layers{};
    int reduction_layers{};
    uint64_t toffoli{};
    uint64_t cnot{};
    uint64_t full_depth{};
    uint64_t current_depth{};
    uint64_t toffoli_depth{};
    uint64_t qubits{};
    uint128_t dw_full{};
    uint128_t dw_current{};
    uint128_t tdw{};
};

using LayerProvider = std::function<LayerStats(int, PointAddMode)>;

std::vector<ScheduledLayer> build_balanced_schedule(int total_inputs, int w);
BalancedResult estimate_balanced(const GFConfig& config, int w,
                                 const LayerProvider& layer_provider);
std::string uint128_to_string(uint128_t value);
long double uint128_log2(uint128_t value);
