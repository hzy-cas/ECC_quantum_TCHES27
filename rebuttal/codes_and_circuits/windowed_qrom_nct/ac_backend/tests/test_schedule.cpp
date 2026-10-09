#include "../src/AlgorithmFlow/BalancedResource.h"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    const auto config = GFConfig::setup(163, "unused-in-schedule-test");
    const auto config571 = GFConfig::setup(571, "unused-in-schedule-test");
    require(config.scan_w_max() == 128, "n=163 must remain capped at w=128");
    require(config571.scan_w_max() == 256, "n=571 must support w up to 256");
    const auto schedule571 = build_balanced_schedule(config571.total_shor_inputs(), 256);
    require(!schedule571.empty() && schedule571.front().lanes == 256,
            "n=571, w=256 schedule was not constructed");
    require(schedule571[4].mode == PointAddMode::Reduction && schedule571[4].lanes == 128,
            "n=571, w=256 must start reduction with 128 parallel additions");

    const auto schedule = build_balanced_schedule(config.total_shor_inputs(), 30);

    std::vector<int> accumulation;
    std::vector<int> reduction;
    for (const auto& layer : schedule) {
        if (layer.mode == PointAddMode::Accumulation) accumulation.push_back(layer.lanes);
        else reduction.push_back(layer.lanes);
    }

    require((accumulation == std::vector<int>{30, 30, 30, 30, 30, 30, 30, 30, 30, 28}),
            "unexpected accumulation schedule");
    require((reduction == std::vector<int>{15, 7, 4, 2, 1}),
            "unexpected reduction schedule");

    int loaded_after_initialization = 30;
    for (int lanes : accumulation) loaded_after_initialization += lanes;
    require(loaded_after_initialization == 328, "schedule does not consume all inputs");


    const auto serial = build_balanced_schedule(328, 1);
    require(serial.size() == 327, "w=1 should have N-1 layers");
    for (const auto& layer : serial) {
        require(layer.mode == PointAddMode::Accumulation, "w=1 must not reduce");
        require(layer.lanes == 1, "w=1 layer has wrong lane count");
    }


    const LayerStats fake{10000, 2, 3, 5, 7, 11};
    const auto result = estimate_balanced(config, 1,
        [&](int lanes, PointAddMode) {
            require(lanes == 1, "provider received wrong lane count");
            return fake;
        });
    require(result.toffoli == 1308, "Toffoli aggregation mismatch");
    require(result.cnot == 4240, "CNOT aggregation mismatch");
    require(result.full_depth == 3310, "full-depth aggregation mismatch");
    require(result.current_depth == 4618, "current-depth aggregation mismatch");
    require(result.toffoli_depth == 7194, "Toffoli-depth aggregation mismatch");
    require(result.qubits == 223205, "peak-width aggregation mismatch");
    require(result.dw_full == static_cast<uint128_t>(3310) * 223205, "DW mismatch");
    require(result.tdw == static_cast<uint128_t>(7194) * 223205, "TDW mismatch");

    std::cout << "schedule/resource aggregation tests passed\n";
    return 0;
}
