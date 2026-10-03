#include "../BasicLogic/QuantumLib.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include <iomanip>

namespace {

    // Binary curve coefficient used by the point-addition implementation.
    constexpr uint64_t CURVE_A = 3;

    template <typename T>
    std::vector<T> concat(const std::vector<T>& a, const std::vector<T>& b) {
        std::vector<T> res;
        res.reserve(a.size() + b.size());
        res.insert(res.end(), a.begin(), a.end());
        res.insert(res.end(), b.begin(), b.end());
        return res;
    }

    template <typename T>
    std::vector<T> sub_vec(const std::vector<T>& v, int start, int end) {
        if (start >= v.size()) return {};
        if (end > v.size()) end = v.size();
        return std::vector<T>(v.begin() + start, v.begin() + end);
    }

    std::vector<std::vector<int>> slice_vector(const std::vector<int>& flat, int chunk_size, int count) {
        std::vector<std::vector<int>> result;
        result.reserve(count);
        for (int i = 0; i < count; ++i) {
            if ((i + 1) * chunk_size > flat.size()) break;
            result.emplace_back(flat.begin() + i * chunk_size, flat.begin() + (i + 1) * chunk_size);
        }
        return result;
    }

    std::vector<bool> int_to_bits(int n, uint64_t val) {
        std::vector<bool> bits(n, false);
        for(int i=0; i<n && i<64; ++i) {
            if((val >> i) & 1) bits[i] = true;
        }
        return bits;
    }

} 


// =========================================================
// Shor's algorithm logic implementation
// =========================================================

QuantumContext::TreeResult QuantumContext::cons_mul(std::vector<std::vector<int>> xlist, int count, std::vector<int> ancilla) {
    std::vector<std::vector<std::vector<int>>> tree;
    tree.push_back(xlist);
    int num = config.block_size * 2;
    int total_mul_size = config.mul_cost_size;
    
    auto current_level = xlist;
    while (tree.back().size() > 1) {
        std::vector<std::vector<int>> next_level;
        int i = 0, j = 0;
        while(i < current_level.size()) {
            if (i + 1 < current_level.size()) {
             
                auto op1 = concat(current_level[i], sub_vec(ancilla, num*j, num*j + config.block_size));
                auto op2 = concat(current_level[i+1], sub_vec(ancilla, num*j + config.block_size, num*j + config.block_size*2));
                
                auto mul_anc = sub_vec(Toffoli_qubits, count, count + total_mul_size);
                
                
                auto product = mul(op1, op2, std::vector<int>(Toffoli_qubits.begin() + count, Toffoli_qubits.end()), 0);
                
                count += total_mul_size;
                next_level.push_back(product);
                i += 2; j += 1;
            } else {
                next_level.push_back(current_level[i]);
                i += 1;
            }
        }
        tree.push_back(next_level);
        current_level = next_level;
    }
    return {tree.back()[0], count, ancilla, tree};
}

std::tuple<std::vector<std::vector<int>>, std::vector<int>, int> 
QuantumContext::Montgomerytrick(
    std::vector<std::vector<int>> x1x2_list, 
    std::vector<std::vector<int>> y1_list, 
    int count, 
    std::vector<int> ancilla, 
    std::vector<int> ancilla2, 
    std::vector<int> sqr) 
{
    int num = config.block_size * 2;
    int parallel_num = x1x2_list.size();
    int total_mul_size = config.mul_cost_size;

    // Case 0: Empty
    if (parallel_num == 0) return {{}, ancilla, count};

    // Case 1: Single element (Direct Inversion)
    if (parallel_num == 1) {
        auto anc_first = sub_vec(ancilla, 0, num);
        auto anc_last = sub_vec(ancilla, ancilla.size() - num, ancilla.size());
        
        auto x1x2_inv = inversion_logic(x1x2_list[0], sqr, anc_first, anc_last, count);
        
        auto mul_anc = sub_vec(Toffoli_qubits, count, count + total_mul_size);
        auto lambda_res = mul(
            concat(y1_list[0], sub_vec(ancilla, 0, config.block_size)),
            concat(x1x2_inv, sub_vec(ancilla, config.block_size, config.block_size * 2)),
            mul_anc, 0
        );
        count += total_mul_size;
        return {{lambda_res}, ancilla, count};
    }

    // Case 2: Parallel > 1 (Montgomery Trick)
    
    auto tree_res = cons_mul(x1x2_list, count, ancilla);
    std::vector<int> product_final = tree_res.product;
    count = tree_res.count;
    auto tree = tree_res.tree;

    auto anc_first = sub_vec(ancilla, 0, num);
    auto anc_last = sub_vec(ancilla, ancilla.size() - num, ancilla.size());
    
    auto product_inv = inversion_logic(product_final, sqr, anc_first, anc_last, count);

    std::vector<std::vector<std::vector<int>>> inv_tree;
    inv_tree.push_back({product_inv});

    for (int level = (int)tree.size() - 2; level >= 0; --level) {
        auto& current_level_nodes = tree[level];
        auto& parent_inv_nodes = inv_tree.back();
        std::vector<std::vector<int>> next_inv_level;
        
        int inv_idx = 0;
        int idx = 0;
        int i = 0;
        
        while (i < current_level_nodes.size()) {
            if (i + 1 < current_level_nodes.size()) {
                // Pair Logic
                size_t begin = gm.current_pointer();
                
                auto anc2_slice = sub_vec(ancilla2, idx, idx + config.n);
                CNOT_N(parent_inv_nodes[inv_idx], anc2_slice);
                auto copy_inv = anc2_slice;
                idx += config.n;
                
                size_t end = gm.current_pointer();

                // Calc Left Inverse
                auto op1_a = concat(parent_inv_nodes[inv_idx], sub_vec(ancilla, num*i, num*i + config.block_size));
                auto op2_b = concat(current_level_nodes[i+1], sub_vec(ancilla, num*i + config.block_size, num*i + config.block_size*2));
                auto mul_anc_a = sub_vec(Toffoli_qubits, count, count + total_mul_size);
                auto a_inv = mul(op1_a, op2_b, mul_anc_a, 0);
                count += total_mul_size;

                // Calc Right Inverse
                auto op1_b = concat(copy_inv, sub_vec(ancilla, num*(i+1), num*(i+1) + config.block_size));
                auto op2_a = concat(current_level_nodes[i], sub_vec(ancilla, num*(i+1) + config.block_size, num*(i+1) + config.block_size*2));
                auto mul_anc_b = sub_vec(Toffoli_qubits, count, count + total_mul_size);
                auto b_inv = mul(op1_b, op2_a, mul_anc_b, 0);
                count += total_mul_size;

                gm.replay_reverse(begin, end);
                
                next_inv_level.push_back(a_inv);
                next_inv_level.push_back(b_inv);
                
                i += 2;
                inv_idx++;
            } else {
                // Odd Element
                next_inv_level.push_back(parent_inv_nodes[inv_idx]);
                i += 1;
                inv_idx++;
            }
        }
        inv_tree.push_back(next_inv_level);
    }

    auto& inverses = inv_tree.back();
    std::vector<std::vector<int>> lambda_list;
    
    for (int i = 0; i < parallel_num; ++i) {
        auto op_y = concat(y1_list[i], sub_vec(ancilla, num*i, num*i + config.block_size));
        auto op_inv = concat(inverses[i], sub_vec(ancilla, num*i + config.block_size, num*i + config.block_size*2));
        auto mul_anc = sub_vec(Toffoli_qubits, count, count + total_mul_size);
        
        auto lambda_i = mul(op_y, op_inv, mul_anc, 0);
        count += total_mul_size;
        lambda_list.push_back(lambda_i);
    }

    return {lambda_list, ancilla, count};
}

void QuantumContext::shor_accumulation_step(
    int parallel_num, 
    const std::vector<int>& q_list,
    std::vector<std::vector<int>>& x1_list,
    std::vector<std::vector<int>>& y1_list,
    const std::vector<std::vector<bool>>& x2_val_list,
    const std::vector<std::vector<bool>>& y2_val_list,
    std::vector<int>& ancilla,
    std::vector<std::vector<int>>& x1x2_list,
    std::vector<int>& ancilla2,
    std::vector<int>& sqr,
    std::vector<int>& cswap_ancilla) 
{
    std::vector<bool> curve_a_bits = int_to_bits(config.n, CURVE_A);

    int num = config.block_size * 2;
    int total_mul_size = config.mul_cost_size;
    int n = config.n;

    for(int i=0; i<parallel_num; ++i) {
        CNOT_N(x1_list[i], x1x2_list[i]);
        CONST_ADD(y1_list[i], y2_val_list[i]);
        CONST_ADD(x1x2_list[i], x2_val_list[i]);
    }

    int count = 0;
    auto [lambda_list, ancilla_ret, count_ret] = Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr);
    count = count_ret;

    for(int i=0; i<parallel_num; ++i) {
        CONST_ADD(x1x2_list[i], x2_val_list[i]); 
        CONST_ADD(x1x2_list[i], curve_a_bits); // XOR curve_a
        CONST_ADD(y1_list[i], y2_val_list[i]); 

        auto sqr_in = concat(lambda_list[i], x1x2_list[i]);
        auto [sqr_res1, sqr_res2] = Squaring_plus_input(sqr_in, 2*n);
        lambda_list[i] = sqr_res1;
        x1x2_list[i] = sqr_res2;

        auto op_x = concat(x1x2_list[i], sub_vec(ancilla, num*i, num*i + config.block_size));
        auto op_l = concat(lambda_list[i], sub_vec(ancilla, num*i + config.block_size, num*i + config.block_size*2));
        auto mul_anc = sub_vec(Toffoli_qubits, count, count + total_mul_size);
        
        auto y_target = mul(op_x, op_l, mul_anc, 0);
        count += total_mul_size;

        CONST_ADD(x1x2_list[i], x2_val_list[i]);
        CONST_ADD(y_target, y2_val_list[i]);
        CNOT_N(x1x2_list[i], y_target);

        auto current_cswap_anc = sub_vec(cswap_ancilla, i * (2*n - 1), (i+1) * (2*n - 1));
        
        size_t copy_start = gm.current_pointer();
        auto control_bits = copy_parallel(q_list[i], current_cswap_anc, 2*n);
        size_t copy_end = gm.current_pointer();

        for(int b=0; b<n; ++b) {
            CSWAP(control_bits[b], x1_list[i][b], x1x2_list[i][b]);
        }
        for(int b=0; b<n; ++b) {
            CSWAP(control_bits[n+b], y1_list[i][b], y_target[b]);
        }

        gm.replay_reverse(copy_start, copy_end);
    }
}

void QuantumContext::shor_reduction_step(
    int parallel_num,
    std::vector<std::vector<int>>& x1_list,
    std::vector<std::vector<int>>& y1_list,
    std::vector<std::vector<int>>& x2_list,
    std::vector<std::vector<int>>& y2_list,
    std::vector<int>& ancilla,
    std::vector<std::vector<int>>& x1x2_list,
    std::vector<int>& ancilla2,
    std::vector<int>& sqr) 
{
    std::vector<bool> curve_a_bits = int_to_bits(config.n, CURVE_A);

    int num = config.block_size * 2;
    int total_mul_size = config.mul_cost_size;
    int n = config.n;

    for(int i=0; i<parallel_num; ++i) {
        CNOT_N(x1_list[i], x1x2_list[i]);
        CNOT_N(x2_list[i], x1x2_list[i]);
        CNOT_N(y2_list[i], y1_list[i]);
    }

    int count = 0;
    auto [lambda_list, ancilla_ret, count_ret] = Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr);
    count = count_ret;

    for(int i=0; i<parallel_num; ++i) {
        CONST_ADD(x1x2_list[i], curve_a_bits);
        
        auto sqr_in = concat(lambda_list[i], x1x2_list[i]);
        auto [sqr_res1, sqr_res2] = Squaring_plus_input(sqr_in, 2*n);
        lambda_list[i] = sqr_res1;
        x1x2_list[i] = sqr_res2;

        CNOT_N(x2_list[i], x1x2_list[i]);

        auto op_x = concat(x1x2_list[i], sub_vec(ancilla, num*i, num*i + config.block_size));
        auto op_l = concat(lambda_list[i], sub_vec(ancilla, num*i + config.block_size, num*i + config.block_size*2));
        auto mul_anc = sub_vec(Toffoli_qubits, count, count + total_mul_size);
        
        auto y_target = mul(op_x, op_l, mul_anc, 0);
        count += total_mul_size;

        CNOT_N(y2_list[i], y_target);
        
        CNOT_N(x2_list[i], x1x2_list[i]); 
        CNOT_N(x1x2_list[i], y_target);   
    }
}

std::tuple<int, uint64_t, uint64_t, uint32_t, uint32_t, uint32_t> 
QuantumContext::run_shor_logic(int parallel_num, std::string mode) {
    gm.clear();
    
    int n = config.n;
    int block_size = config.block_size;
    int total_mul_size = config.mul_cost_size;
    
    long long nums1 = (long long)block_size * 4 * parallel_num;
    long long nums2 = (long long)total_mul_size * (3 * (parallel_num - 1) + config.nums2_const + 2 * parallel_num);
    
    int current_idx = 0;
    
    std::vector<int> q_list;
    if (mode == "accumulation") {
        for (int i = 0; i < parallel_num; ++i) q_list.push_back(current_idx++);
    }
    
    std::vector<int> x1; 
    for(int i=0; i < n * parallel_num; ++i) x1.push_back(current_idx++);
    
    std::vector<int> x2;
    if (mode == "reduction") {
        for(int i=0; i < n * parallel_num; ++i) x2.push_back(current_idx++);
    }
    
    std::vector<int> x1x2;
    for(int i=0; i < n * parallel_num; ++i) x1x2.push_back(current_idx++);
    
    std::vector<int> y1;
    for(int i=0; i < n * parallel_num; ++i) y1.push_back(current_idx++);
    
    std::vector<int> y2;
    if (mode == "reduction") {
        for(int i=0; i < n * parallel_num; ++i) y2.push_back(current_idx++);
    }
    
    auto x1_list = slice_vector(x1, n, parallel_num);
    auto y1_list = slice_vector(y1, n, parallel_num);
    auto x1x2_list = slice_vector(x1x2, n, parallel_num);
    
    std::vector<std::vector<int>> x2_list_q;
    std::vector<std::vector<int>> y2_list_q;
    
    std::vector<std::vector<bool>> x2_val_list;
    std::vector<std::vector<bool>> y2_val_list;
    
    if (mode == "reduction") {
        x2_list_q = slice_vector(x2, n, parallel_num);
        y2_list_q = slice_vector(y2, n, parallel_num);
    } else {
        std::vector<bool> all_ones(n, true);
        for(int i=0; i<parallel_num; ++i) {
            x2_val_list.push_back(all_ones);
            y2_val_list.push_back(all_ones);
        }
    }
    
    std::vector<int> ancilla;
    for(int i=0; i < nums1; ++i) ancilla.push_back(current_idx++);
    
    this->Toffoli_qubits.clear();
    for(int i=0; i < nums2; ++i) this->Toffoli_qubits.push_back(current_idx++);
    
    std::vector<int> sqr;
    for(int i=0; i < 4 * n; ++i) sqr.push_back(current_idx++);
    
    std::vector<int> ancilla2;
    int anc2_size = n * (parallel_num / 2);
    for(int i=0; i < anc2_size; ++i) ancilla2.push_back(current_idx++);
    
    std::vector<int> cswap_ancilla;
    if (mode == "accumulation") {
        int size = (2 * n - 1) * parallel_num;
        for(int i=0; i < size; ++i) cswap_ancilla.push_back(current_idx++);
    }
    
    int total_qubits_used = current_idx;

    if (mode == "accumulation") {
        shor_accumulation_step(
            parallel_num, q_list, x1_list, y1_list,
            x2_val_list, y2_val_list,
            ancilla, x1x2_list, ancilla2, sqr, cswap_ancilla
        );
    } else {
        shor_reduction_step(
            parallel_num, x1_list, y1_list,
            x2_list_q, y2_list_q,
            ancilla, x1x2_list, ancilla2, sqr
        );
    }
    
    auto stats = gm.optimize_and_count(total_qubits_used);
    return {
        total_qubits_used, 
        std::get<3>(stats), 
        std::get<4>(stats), 
        std::get<0>(stats), 
        std::get<1>(stats), 
        std::get<2>(stats)
    };
}

void QuantumContext::estimate_total_resources(int total_inputs) {
    uint64_t total_Toffoli = 0;
    uint64_t total_CNOT = 0;
    uint64_t total_fDepth = 0;
    uint64_t total_cDepth = 0;
    uint64_t total_tDepth = 0;
    
    long long peak_qubits = 0;
    long long clear_qubits = 0;
    
    int n = config.n;
    int current = total_inputs;
    int level = 0;
    
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "         Shor Algorithm Resource Estimation Strategy" << std::endl;
    std::cout << "                 (Strict Rules Applied)" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    while (current > 1) {
        int ops = current / 2;
        
        if (level == 0) {
            auto [length, t, c, fd, cd, td] = run_shor_logic(ops, "accumulation");
            
            double log_val = std::ceil(std::log2(2.0 * n));
            int depth_overhead = (int)(2 * log_val + 1);
            
            total_fDepth += (depth_overhead + fd);
            total_cDepth += (depth_overhead + cd);
            total_tDepth += td;
            
            long long cnot_overhead = (long long)2 * ops * (2 * n - 1) + (long long)2 * n * ops;
            total_CNOT += (cnot_overhead + c);
            total_Toffoli += t;
            
            peak_qubits = length - ops;
            clear_qubits = (long long)2 * n + n * (ops / 2) + (long long)(2 * n - 1) * ops + (long long)config.block_size * 4 * ops;
            
            std::cout << "Level " << level << " (ops=" << ops << "): "
                      << "Depth=" << (depth_overhead + fd) 
                      << ", Qubits=" << peak_qubits 
                      << ", Clear=" << clear_qubits << std::endl;
        }
        else {
            auto [length, t, c, fd, cd, td] = run_shor_logic(ops, "reduction");
            
            total_fDepth += fd;
            total_cDepth += cd;
            total_tDepth += td;
            total_CNOT += c;
            total_Toffoli += t;
            
            long long demand = length - (long long)4 * n * ops;
            
            if (demand > clear_qubits) {
                peak_qubits += (demand - clear_qubits);
                clear_qubits = (long long)2 * n + n * (ops / 2) + (long long)config.block_size * 4 * ops;
            } else {
                long long over = clear_qubits - demand;
                clear_qubits = (long long)2 * n + n * (ops / 2) + over + (long long)config.block_size * 4 * ops;
            }
            
            std::cout << "Level " << level << " (ops=" << ops << "): "
                      << "cDepth=" << cd 
                      << ", Demand=" << demand 
                      << ", Qubits=" << peak_qubits 
                      << ", Clear=" << clear_qubits << std::endl;
        }
        
        current = ops + (current % 2);
        level++;
    }
    
    uint64_t final_Toffoli = 2 * total_Toffoli;
    uint64_t final_CNOT = 2 * total_CNOT + 2 * n;
    
    uint64_t final_cDepth = 2 * total_cDepth + 2;
    uint64_t final_fDepth = 2 * total_fDepth + 2;
    uint64_t final_tDepth = 2 * total_tDepth;
    
    long long final_Qubits = peak_qubits + 4 * n + 2;

    std::cout << "\n" << std::string(40, '=') << std::endl;
    std::cout << "       Final Results (N=" << config.n << ")" << std::endl;
    std::cout << std::string(40, '=') << std::endl;
    std::cout << std::left << std::setw(20) << "Metric" << " | " << std::right << std::setw(15) << "Value" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    std::cout << std::left << std::setw(20) << "Toffoli Gates" << " | " << std::right << std::setw(15) << final_Toffoli << std::endl;
    std::cout << std::left << std::setw(20) << "CNOT Gates" << " | " << std::right << std::setw(15) << final_CNOT << std::endl;
    std::cout << std::left << std::setw(20) << "Full Depth" << " | " << std::right << std::setw(15) << final_fDepth << std::endl;
    std::cout << std::left << std::setw(20) << "Current Depth" << " | " << std::right << std::setw(15) << final_cDepth << std::endl;
    std::cout << std::left << std::setw(20) << "Toffoli Depth" << " | " << std::right << std::setw(15) << final_tDepth << std::endl;
    std::cout << std::left << std::setw(20) << "Qubits" << " | " << std::right << std::setw(15) << final_Qubits << std::endl;
    std::cout << std::string(40, '-') << std::endl;
}