// src/Infrastructure/GateManager.cpp
#include "GateManager.h"
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include <vector>

GateManager::GateManager() {
    current_chunk.reserve(CHUNK_SIZE);
}

void GateManager::_ensure_room(size_t words) {
    if (current_chunk.size() + words > static_cast<size_t>(CHUNK_SIZE)) {
        chunks.push_back(std::move(current_chunk));
        current_chunk = std::vector<uint64_t>();
        current_chunk.reserve(CHUNK_SIZE);
    }
}

void GateManager::_append(uint64_t val) {
    _ensure_room(1);
    current_chunk.push_back(val);
}

void GateManager::add_X(uint32_t t) {
    uint64_t val = ((uint64_t)OP_X << SHIFT_OP) | t;
    _append(val);
}

void GateManager::add_CNOT(uint32_t c, uint32_t t) {
    uint64_t val = ((uint64_t)OP_CNOT << SHIFT_OP) | ((uint64_t)c << SHIFT_ARG) | t;
    _append(val);
}

void GateManager::add_Toffoli(uint32_t c1, uint32_t c2, uint32_t t) {
    uint64_t val1 = ((uint64_t)OP_TOFFOLI << SHIFT_OP) | ((uint64_t)c1 << SHIFT_ARG) | c2;
    // A Toffoli occupies two words. Keep them in the same chunk so its target
    // can never appear at offset zero without the corresponding header.
    _ensure_room(2);
    current_chunk.push_back(val1);
    current_chunk.push_back((uint64_t)t);
}

size_t GateManager::current_pointer() const {
    size_t total = 0;
    for (const auto& c : chunks) total += c.size();
    total += current_chunk.size();
    return total;
}

void GateManager::clear() {
    chunks.clear();
    current_chunk.clear();
    current_chunk.reserve(CHUNK_SIZE);
}

void GateManager::replay_reverse(size_t start_idx, size_t end_idx) {
    if (start_idx >= end_idx) return;

    std::vector<uint64_t> replay_buffer;
    

    std::vector<const std::vector<uint64_t>*> all_chunks_ptrs;
    all_chunks_ptrs.reserve(chunks.size() + 1);
    for (const auto& c : chunks) all_chunks_ptrs.push_back(&c);
    all_chunks_ptrs.push_back(&current_chunk);

    std::vector<size_t> chunk_starts;
    size_t pos = 0;
    for (const auto* chunk_ptr : all_chunks_ptrs) {
        chunk_starts.push_back(pos);
        pos += chunk_ptr->size();
    }
    size_t total_size = pos;

    if (end_idx > total_size) return;

    auto safe_get_value = [&](size_t abs_idx) -> std::pair<bool, uint64_t> {
        if (abs_idx >= total_size) return {false, 0};

        for (size_t c_i = 0; c_i < all_chunks_ptrs.size(); ++c_i) {
            size_t chunk_start = chunk_starts[c_i];
            size_t chunk_end = chunk_start + all_chunks_ptrs[c_i]->size();
            if (abs_idx >= chunk_start && abs_idx < chunk_end) {
                return {true, (*all_chunks_ptrs[c_i])[abs_idx - chunk_start]};
            }
        }
        return {false, 0};
    };

    struct ReplayOp {
        uint32_t op_type;
        uint32_t c1, c2, t;
    };
    std::vector<ReplayOp> replay_ops;

    for (size_t abs_idx = end_idx; abs_idx > start_idx; --abs_idx) {
        size_t idx = abs_idx - 1;
        auto [found, val] = safe_get_value(idx);
        if (!found) continue;
        

        uint32_t op = val >> SHIFT_OP;
        
        if (op == OP_NOP) continue;

        if (op == OP_DATA) {
            uint32_t t = (uint32_t)val;
            if (idx > 0) {
                auto [found_prev, val_prev] = safe_get_value(idx - 1);

                if (found_prev && ((val_prev >> SHIFT_OP) == OP_TOFFOLI)) {
                    uint32_t c2 = val_prev & MASK_ARG;
                    uint32_t c1 = (val_prev >> SHIFT_ARG) & MASK_ARG;
                    replay_ops.push_back({OP_TOFFOLI, c1, c2, t});
                }
            }
        } else if (op == OP_TOFFOLI) {
            continue;
        } else if (op == OP_CNOT) {
            uint32_t t = val & MASK_ARG;
            uint32_t c = (val >> SHIFT_ARG) & MASK_ARG;
            replay_ops.push_back({OP_CNOT, c, 0, t});
        } else if (op == OP_X) {
            uint32_t t = val & MASK_ARG;
            replay_ops.push_back({OP_X, 0, 0, t});
        }
    }

    for (const auto& op : replay_ops) {
        if (op.op_type == OP_X) add_X(op.t);
        else if (op.op_type == OP_CNOT) add_CNOT(op.c1, op.t);
        else if (op.op_type == OP_TOFFOLI) add_Toffoli(op.c1, op.c2, op.t);
    }
}

// =========================================================

// =========================================================
std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> GateManager::optimize_and_count(uint32_t width) {
    
    std::vector<std::vector<uint64_t>*> all_chunks_ptrs;
    all_chunks_ptrs.reserve(chunks.size() + 1);
    for(auto& c : chunks) all_chunks_ptrs.push_back(&c);
    all_chunks_ptrs.push_back(&current_chunk);

    uint64_t VAL_NOP = ((uint64_t)OP_NOP << SHIFT_OP) | MASK_ARG;


    std::vector<std::vector<uint64_t>> wire_stacks(width);
    
    struct PendingState {
        uint32_t c1, c2;
        uint64_t ptr;
        bool active = false;
    } pending;


    for (size_t c_idx = 0; c_idx < all_chunks_ptrs.size(); ++c_idx) {
        std::vector<uint64_t>& chunk = *all_chunks_ptrs[c_idx];
        

        uint64_t c_idx_shifted = (uint64_t)c_idx << 32;

        for (size_t i = 0; i < chunk.size(); ++i) {
            uint64_t val = chunk[i];
            

            uint32_t op = val >> SHIFT_OP;

            if (pending.active) {
                uint32_t t = val & MASK_ARG;
                uint64_t gate_ptr = pending.ptr;

                auto& s_c1 = wire_stacks[pending.c1];
                auto& s_c2 = wire_stacks[pending.c2];
                auto& s_t  = wire_stacks[t];

                if (!s_c1.empty() && !s_c2.empty() && !s_t.empty()) {
                    if (s_c1.back() == s_c2.back() && s_c2.back() == s_t.back()) {
                        uint64_t old_ptr = s_c1.back();
                        uint32_t oc = old_ptr >> 32;
                        uint32_t oo = old_ptr & 0xFFFFFFFF;

                        auto& old_chunk = *all_chunks_ptrs[oc];
                        uint64_t old_header = old_chunk[oo];
                        if ((old_header >> SHIFT_OP) == OP_TOFFOLI &&
                            oo + 1 < old_chunk.size() &&
                            (old_chunk[oo + 1] & MASK_ARG) == t) {
                            old_chunk[oo] = VAL_NOP;
                            old_chunk[oo + 1] = VAL_NOP;
                            chunk[i] = VAL_NOP;
                            
                            uint32_t pc = pending.ptr >> 32;
                            uint32_t po = pending.ptr & 0xFFFFFFFF;
                            (*all_chunks_ptrs[pc])[po] = VAL_NOP;

                            s_c1.pop_back(); s_c2.pop_back(); s_t.pop_back();
                            pending.active = false;
                            continue; 
                        }
                    }
                }
                s_c1.push_back(gate_ptr); s_c2.push_back(gate_ptr); s_t.push_back(gate_ptr);
                pending.active = false;
            }
            else if (op == OP_CNOT) {
                uint32_t t = val & MASK_ARG;
                uint32_t c = (val >> SHIFT_ARG) & MASK_ARG;
                auto& s_c = wire_stacks[c];
                auto& s_t = wire_stacks[t];

                if (!s_c.empty() && !s_t.empty() && s_c.back() == s_t.back()) {
                    uint64_t old_ptr = s_c.back();
                    uint32_t oc = old_ptr >> 32;
                    uint32_t oo = old_ptr & 0xFFFFFFFF;
                    uint64_t old_val = (*all_chunks_ptrs[oc])[oo];

                    if ((old_val >> SHIFT_OP) == OP_CNOT && (old_val & MASK_ARG) == t) {
                        (*all_chunks_ptrs[oc])[oo] = VAL_NOP;
                        chunk[i] = VAL_NOP;
                        s_c.pop_back(); s_t.pop_back();
                        continue;
                    }
                }
                uint64_t ptr = c_idx_shifted | i;
                s_c.push_back(ptr); s_t.push_back(ptr);
            }
            else if (op == OP_TOFFOLI) {
                pending.c1 = (val >> SHIFT_ARG) & MASK_ARG;
                pending.c2 = val & MASK_ARG;
                pending.ptr = c_idx_shifted | i;
                pending.active = true;
            }
            else if (op == OP_X) {
                uint32_t t = val & MASK_ARG;
                auto& s_t = wire_stacks[t];
                
                if (!s_t.empty()) {
                    uint64_t old_ptr = s_t.back();
                    uint32_t oc = old_ptr >> 32;
                    uint32_t oo = old_ptr & 0xFFFFFFFF;
                    uint64_t old_val = (*all_chunks_ptrs[oc])[oo];
                    
                    if ((old_val >> SHIFT_OP) == OP_X) {
                        (*all_chunks_ptrs[oc])[oo] = VAL_NOP;
                        chunk[i] = VAL_NOP;
                        s_t.pop_back();
                        continue;
                    }
                }
                uint64_t ptr = c_idx_shifted | i;
                s_t.push_back(ptr);
            }
        }
    }

    // Phase 2: Depth Calculation

    std::vector<uint64_t> wire_tf_depths(width, 0);
    std::vector<uint64_t> wires_f_depth(width, 0);
    
    uint64_t t_count = 0;
    uint64_t c_count = 0;
    
    std::vector<std::vector<uint64_t>> clifford_buckets;
    std::vector<std::set<uint32_t>> mt_qubit_buckets;
    
    auto ensure_buckets = [&](size_t d) {
        while (clifford_buckets.size() <= d) clifford_buckets.emplace_back();
        while (mt_qubit_buckets.size() <= d + 1) mt_qubit_buckets.emplace_back();
    };
    
    struct ProcessingState { uint32_t c1, c2; bool active = false; } p2;

    for (const auto* chunk_ptr : all_chunks_ptrs) {
        for (uint64_t val : *chunk_ptr) {
            uint32_t op = val >> SHIFT_OP; // use uint32_t
            if (op == OP_NOP) continue;
            
            if (p2.active) {
                uint32_t t = val & MASK_ARG;
                t_count++;
                

                uint64_t tf_d = std::max({wire_tf_depths[p2.c1], wire_tf_depths[p2.c2], wire_tf_depths[t]}) + 1;
                wire_tf_depths[p2.c1] = wire_tf_depths[p2.c2] = wire_tf_depths[t] = tf_d;
                
                ensure_buckets(tf_d);
                mt_qubit_buckets[tf_d].insert(p2.c1);
                mt_qubit_buckets[tf_d].insert(p2.c2);
                mt_qubit_buckets[tf_d].insert(t);
                
                uint64_t f_mx = std::max({wires_f_depth[p2.c1], wires_f_depth[p2.c2], wires_f_depth[t]}) + 1;
                wires_f_depth[p2.c1] = wires_f_depth[p2.c2] = wires_f_depth[t] = f_mx;
                
                p2.active = false;
            }
            else if (op == OP_CNOT) {
                c_count++;
                uint32_t t = val & MASK_ARG;
                uint32_t c = (val >> SHIFT_ARG) & MASK_ARG;
                
                uint64_t tf_d = std::max(wire_tf_depths[c], wire_tf_depths[t]);
                wire_tf_depths[c] = wire_tf_depths[t] = tf_d;
                
                ensure_buckets(tf_d);
                clifford_buckets[tf_d].push_back(val);
                
                uint64_t f_mx = std::max(wires_f_depth[c], wires_f_depth[t]) + 1;
                wires_f_depth[c] = wires_f_depth[t] = f_mx;
            }
            else if (op == OP_TOFFOLI) {
                p2.c1 = (val >> SHIFT_ARG) & MASK_ARG;
                p2.c2 = val & MASK_ARG;
                p2.active = true;
            }
            else if (op == OP_X) {
                uint32_t t = val & MASK_ARG;
                uint64_t tf_d = wire_tf_depths[t]; // use uint64_t
                
                ensure_buckets(tf_d);
                clifford_buckets[tf_d].push_back(val);
                
                wires_f_depth[t]++;
            }
        }
    }

    // Phase 3: Current Depth Reconstruction

    std::vector<uint64_t> wires_c_depth(width, 0);
    size_t limit = std::max(clifford_buckets.size(), mt_qubit_buckets.size());
    
    for (size_t d = 0; d < limit; ++d) {
        if (d < clifford_buckets.size()) {
            for (uint64_t val : clifford_buckets[d]) {
                uint32_t op = val >> SHIFT_OP;
                if (op == OP_CNOT) {
                    uint32_t t = val & MASK_ARG;
                    uint32_t c = (val >> SHIFT_ARG) & MASK_ARG;
                    uint64_t new_v = std::max(wires_c_depth[c], wires_c_depth[t]) + 1;
                    wires_c_depth[c] = wires_c_depth[t] = new_v;
                } else if (op == OP_X) {
                    wires_c_depth[val & MASK_ARG]++;
                }
            }
        }
        
        size_t target_mt = d + 1;
        if (target_mt < mt_qubit_buckets.size() && !mt_qubit_buckets[target_mt].empty()) {
            const auto& qubits = mt_qubit_buckets[target_mt];
            uint64_t d_max = 0; // use uint64_t
            for (auto q : qubits) d_max = std::max(d_max, wires_c_depth[q]);
            
            d_max += 1;
            for (auto q : qubits) wires_c_depth[q] = d_max;
        }
    }

    // Final Results
    uint64_t final_toffoli_depth = 0;
    uint64_t final_full_depth = 0;
    uint64_t final_curr_depth = 0;
    
    if (width > 0) {
        for(auto v : wire_tf_depths) final_toffoli_depth = std::max(final_toffoli_depth, v);
        for(auto v : wires_f_depth) final_full_depth = std::max(final_full_depth, v);
        for(auto v : wires_c_depth) final_curr_depth = std::max(final_curr_depth, v);
    }
    
    return {final_full_depth, final_curr_depth, final_toffoli_depth, t_count, c_count};
}
