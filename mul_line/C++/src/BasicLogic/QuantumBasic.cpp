#include "QuantumLib.h"
#include <iostream>

QuantumContext::QuantumContext(GFConfig cfg) : config(cfg) {}

void QuantumContext::X(int a) { gm.add_X(a); }
void QuantumContext::CNOT(int c, int t) { gm.add_CNOT(c, t); }
void QuantumContext::Toffoli(int c1, int c2, int t) { gm.add_Toffoli(c1, c2, t); }

void QuantumContext::CNOT_N(const std::vector<int>& a, const std::vector<int>& b) {
    for(size_t i=0; i<config.n; ++i) gm.add_CNOT(a[i], b[i]);
}

void QuantumContext::CONST_ADD(const std::vector<int>& reg, const std::vector<bool>& val_bits) {
    size_t limit = std::min(reg.size(), val_bits.size());
    for(size_t i=0; i<limit; ++i) {
        if (val_bits[i]) {
            gm.add_X(reg[i]);
        }
    }
}

void QuantumContext::CSWAP(int c, int a, int b) {
    gm.add_CNOT(b, a);
    gm.add_Toffoli(c, a, b);
    gm.add_CNOT(b, a);
}

std::vector<int> QuantumContext::copy_parallel(int ctrl, const std::vector<int>& ancilla, int total_count) {
    std::vector<int> controls = {ctrl};
    size_t anc_idx = 0;
    while(controls.size() < total_count) {
        std::vector<int> current = controls;
        for(int src : current) {
            if(controls.size() >= total_count) break;
            if(anc_idx >= ancilla.size()) throw std::runtime_error("Not enough ancilla");
            int target = ancilla[anc_idx++];
            gm.add_CNOT(src, target);
            controls.push_back(target);
        }
    }
    return controls;
}

std::string QuantumContext::get_matrix_path(std::string folder, std::string filename) {
    return config.data_root + "/quantum_" + std::to_string(config.n) + "/" + folder + "/" + filename;
}