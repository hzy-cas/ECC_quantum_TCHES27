#include "QuantumLib.h"
#include <iostream>
#include <random>
#include "config.h"
extern GFConfig current_config;

// ==========================================
// ConsMul & Montgomery
// ==========================================
QuantumContext::ConsMulResult QuantumContext::cons_mul(std::vector<std::vector<int>> xlist, int count, std::vector<int> ancilla)
{
    std::vector<std::vector<std::vector<int>>> tree;
    tree.push_back(xlist);

    int num = current_config.KARATSUBA_ANC_BASE;
    std::vector<std::vector<int>> current_level = xlist;
    int n = current_config.N;

    while (tree.back().size() > 1)
    {
        std::vector<std::vector<int>> next_level;
        int i = 0;
        int j = 0;
        while (i < current_level.size())
        {
            if (i + 1 < current_level.size())
            {
                count = 0;
                std::vector<int> sub_ancilla = slice(ancilla, num * j, num);

                auto res = recursive_karatsuba(current_level[i], current_level[i + 1], n, count, sub_ancilla);

                std::vector<int> product = Reduction(res.res);
                next_level.push_back(product);
                i += 2;
                j += 1;
            }
            else
            {
                next_level.push_back(current_level[i]);
                i += 1;
            }
        }
        tree.push_back(next_level);
        current_level = next_level;
    }

    return {tree.back()[0], count, ancilla, tree};
}

QuantumContext::MontgomeryResult QuantumContext::Montgomerytrick(
    std::vector<std::vector<int>> x1x2_list, std::vector<std::vector<int>> y1_list,
    int count, std::vector<int> ancilla, std::vector<int> ancilla2, std::vector<int> sqr,
    std::vector<int> target_copy, std::vector<int> inv_res)
{
    int n = current_config.N;
    size_t start_ptr_inv = gm.current_pointer();
    int num = current_config.KARATSUBA_ANC_BASE;
    int parallel_num = x1x2_list.size();

    if (parallel_num == 0)
        return {{}, ancilla};

    if (parallel_num == 1)
    {
        count = 0;
        auto inv_res_single = Inverison_Itoh_Tsujii_based(
            x1x2_list[0], n,
            slice(sqr, 0, n), slice(sqr, n, n), slice(sqr, 2 * n, n), slice(sqr, 3 * n, n),
            count, slice(ancilla, 0, num), slice(ancilla, ancilla.size() - num, num));
        std::vector<int> x1x2_inv = inv_res_single.res;

        size_t mid_ptr_inv = gm.current_pointer();
        CNOT_n(x1x2_inv, slice(inv_res, 0, n));
        gm.replay_reverse(start_ptr_inv, mid_ptr_inv);
    }
    else
    {
        count = 0;
        auto cm_res = cons_mul(x1x2_list, count, ancilla);
        std::vector<int> product_final = cm_res.product;
        std::vector<std::vector<std::vector<int>>> tree = cm_res.tree;

        count = 0;
        auto inv_res_root = Inverison_Itoh_Tsujii_based(
            product_final, n,
            slice(sqr, 0, n), slice(sqr, n, n), slice(sqr, 2 * n, n), slice(sqr, 3 * n, n),
            count, slice(ancilla, 0, num), slice(ancilla, ancilla.size() - num, num));
        std::vector<int> product_inv = inv_res_root.res;

        std::vector<std::vector<int>> inv_tree_level;
        inv_tree_level.push_back(product_inv);
        std::vector<std::vector<std::vector<int>>> inv_tree;
        inv_tree.push_back(inv_tree_level);

        for (int level = tree.size() - 2; level >= 0; --level)
        {
            std::vector<std::vector<int>> current_level = tree[level];
            std::vector<std::vector<int>> current_inv_level = inv_tree.back();
            std::vector<std::vector<int>> next_inv_level;

            int inv_idx = 0;
            int i = 0;
            int idx = 0;

            while (i < current_level.size())
            {
                if (i + 1 < current_level.size())
                {
                    std::vector<int> anc2_slice = slice(ancilla2, idx, n);
                    CNOT_n(current_inv_level[inv_idx], anc2_slice);
                    std::vector<int> copy = anc2_slice;

                    count = 0;
                    auto res_a = recursive_karatsuba(
                        current_inv_level[inv_idx], current_level[i + 1], n, count, slice(ancilla, num * i, num));
                    std::vector<int> a_inv = Reduction(res_a.res);

                    count = 0;
                    auto res_b = recursive_karatsuba(
                        copy, current_level[i], n, count, slice(ancilla, num * (i + 1), num));
                    std::vector<int> b_inv = Reduction(res_b.res);

                    CNOT_n(current_inv_level[inv_idx], copy);
                    next_inv_level.push_back(a_inv);
                    next_inv_level.push_back(b_inv);
                    idx += n;
                    i += 2;
                    inv_idx += 1;
                }
                else
                {
                    next_inv_level.push_back(current_inv_level[inv_idx]);
                    i += 1;
                    inv_idx += 1;
                }
            }
            inv_tree.push_back(next_inv_level);
        }

        std::vector<std::vector<int>> inverses = inv_tree.back();
        size_t mid_ptr_inv = gm.current_pointer();
        for (int i = 0; i < parallel_num; ++i)
        {
            CNOT_n(inverses[i], slice(inv_res, i * n, n));
        }
        gm.replay_reverse(start_ptr_inv, mid_ptr_inv);
    }

    this->cnt = 0;

    size_t start_ptr_mul = gm.current_pointer();
    std::vector<std::vector<int>> lambda_list_dirty;

    for (int i = 0; i < parallel_num; ++i)
    {
        count = 0;
        std::vector<int> current_inv = slice(inv_res, i * n, n);
        auto res_lambda = recursive_karatsuba(y1_list[i], current_inv, n, count, slice(ancilla, num * i, num));
        std::vector<int> lambda_i = Reduction(res_lambda.res);
        lambda_list_dirty.push_back(lambda_i);
    }

    size_t mid_ptr_mul = gm.current_pointer();
    std::vector<std::vector<int>> lambda_list_clean;
    for (int i = 0; i < parallel_num; ++i)
    {
        std::vector<int> target_slice = slice(target_copy, i * n, n);
        CNOT_n(lambda_list_dirty[i], target_slice);
        lambda_list_clean.push_back(target_slice);
    }

    gm.replay_reverse(start_ptr_mul, mid_ptr_mul);
    this->cnt = 0;

    return {lambda_list_clean, ancilla};
}

void QuantumContext::shor_level(
    int parallel_num,
    const std::vector<int> &q_list,
    std::vector<std::vector<int>> x1_list, std::vector<std::vector<int>> y1_list,
    const std::vector<int> &x2_vals, const std::vector<int> &y2_vals,  
    std::vector<int> ancilla, std::vector<std::vector<int>> x1x2_list,
    std::vector<int> ancilla2, std::vector<int> sqr, int level,
    std::vector<int> lambda_res, std::vector<int> mult_res, std::vector<int> inv_res,
    std::vector<int> cswap_ancilla)
{
    int a_const = 3; 
    int num = current_config.KARATSUBA_ANC_BASE;
    int n = current_config.N;

    std::vector<int> all_ones(n, 1);
    std::vector<int> x2_xor_a(n, 1);
    if (n > 0) x2_xor_a[0] = 0; 
    if (n > 1) x2_xor_a[1] = 0; 

    for (int i = 0; i < parallel_num; ++i)
    {
        CNOT_n(x1_list[i], x1x2_list[i]);
        
        CONST_ADD_n(y1_list[i], all_ones);    
        CONST_ADD_n(x1x2_list[i], all_ones);  
    }

    int count = 0;
    this->cnt = 0;
    auto mont_res = Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr, lambda_res, inv_res);
    std::vector<std::vector<int>> lamba = mont_res.lambda_list;
    this->cnt = 0;

    for (int i = 0; i < parallel_num; ++i)
    {
        CONST_ADD_n(y1_list[i], all_ones); 
        
        CONST_ADD_n(x1x2_list[i], x2_xor_a); 

        std::vector<int> combined = concat(lamba[i], x1x2_list[i]);
        auto sq_res = Squaring_plus_input(combined, n);
        lamba[i] = sq_res.r1;
        x1x2_list[i] = sq_res.r2;

        size_t start_mul = gm.current_pointer();
        count = 0;
        std::vector<int> sub_anc = slice(ancilla, num * i, num);
        auto res_mul = recursive_karatsuba(lamba[i], x1x2_list[i], n, count, sub_anc);
        std::vector<int> result2 = Reduction(res_mul.res);
        size_t end_mul = gm.current_pointer();

        std::vector<int> target_y = slice(mult_res, i * n, n);
        CNOT_n(result2, target_y);
        gm.replay_reverse(start_mul, end_mul);

        CONST_ADD_n(x1x2_list[i], all_ones); 
        CONST_ADD_n(target_y, all_ones);    
        CNOT_n(x1x2_list[i], target_y);

        int anc_idx = i * (2 * n - 1);
        std::vector<int> current_cswap_anc = slice(cswap_ancilla, anc_idx, 2 * n - 1);
        size_t copy_start = gm.current_pointer();
        std::vector<int> control_bits = copy_parallel(q_list[i], current_cswap_anc, 2 * n);
        size_t copy_end = gm.current_pointer();
        for (int b = 0; b < n; ++b) CSWAP(control_bits[b], x1_list[i][b], x1x2_list[i][b]);
        for (int b = 0; b < n; ++b) CSWAP(control_bits[n + b], y1_list[i][b], target_y[b]);
        gm.replay_reverse(copy_start, copy_end);
    }
}

ResourceStats QuantumContext::shor(int parallel_num, int level, int &length_out)
{
    gm.clear();
    cnt = 0;
    int n = current_config.N;
    long long nums1 = (long long)current_config.KARATSUBA_ANC_BASE * parallel_num;
    long long nums2 = (long long)current_config.TOFFOLI_BASE * (3 * (parallel_num - 1) + current_config.TOFFOLI_OFFSET);

    int idx = 0;
    std::vector<int> q_list = range_vec(idx, parallel_num);
    idx += parallel_num;
    std::vector<int> x1 = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> x1x2 = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> y1 = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;

    std::vector<int> x2_vals(parallel_num), y2_vals(parallel_num);
    for (int i = 0; i < parallel_num; ++i)
    {
        x2_vals[i] = -1; 
        y2_vals[i] = -1;
    }

    std::vector<int> ancilla = range_vec(idx, (int)nums1);
    idx += (int)nums1;
    Toffoli_qubits = range_vec(idx, (int)nums2);
    idx += (int)nums2;
    std::vector<int> sqr = range_vec(idx, 4 * n);
    idx += 4 * n;
    std::vector<int> ancilla2 = range_vec(idx, n * (parallel_num / 2));
    idx += n * (parallel_num / 2);

    std::vector<int> lambda_res = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> mult_res = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> inv_res = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;

    int cswap_len = (2 * n - 1) * parallel_num;
    std::vector<int> cswap_ancilla = range_vec(idx, cswap_len);
    idx += cswap_len;

    length_out = idx;

    std::vector<std::vector<int>> x1_list(parallel_num), y1_list(parallel_num), x1x2_list(parallel_num);
    for (int i = 0; i < parallel_num; ++i)
    {
        x1_list[i] = slice(x1, i * n, n);
        y1_list[i] = slice(y1, i * n, n);
        x1x2_list[i] = slice(x1x2, i * n, n);
    }

    shor_level(parallel_num, q_list, x1_list, y1_list, x2_vals, y2_vals, ancilla, x1x2_list, ancilla2, sqr, level, lambda_res, mult_res, inv_res, cswap_ancilla);
    auto [fd, cd, td, t_c, c_c] = gm.optimize_and_count(length_out);
    return ResourceStats{fd, cd, td, t_c, c_c};
}

// ==========================================
// Algorithm 2: Phase 3 (Reduction)
// ==========================================
void QuantumContext::shor_reduction_level(
    int parallel_num,
    std::vector<std::vector<int>> x1_list, std::vector<std::vector<int>> y1_list,
    std::vector<std::vector<int>> x2_list, std::vector<std::vector<int>> y2_list,
    std::vector<int> ancilla, std::vector<std::vector<int>> x1x2_list,
    std::vector<int> ancilla2, std::vector<int> sqr,
    std::vector<int> lambda_res, std::vector<int> mult_res, std::vector<int> inv_res)
{
    int a_const = 3;
    int num = current_config.KARATSUBA_ANC_BASE;
    int n = current_config.N;

    for (int i = 0; i < parallel_num; ++i)
    {
        CNOT_n(x1_list[i], x1x2_list[i]);
        CNOT_n(x2_list[i], x1x2_list[i]);
        CNOT_n(y2_list[i], y1_list[i]);
    }

    int count = 0;
    this->cnt = 0;
    auto mont_res = Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr, lambda_res, inv_res);
    std::vector<std::vector<int>> lamba = mont_res.lambda_list;
    this->cnt = 0;

    for (int i = 0; i < parallel_num; ++i)
    {
        CONST_ADD_n(x1x2_list[i], a_const);
        std::vector<int> combined = concat(lamba[i], x1x2_list[i]);
        auto sq_res = Squaring_plus_input(combined, n);
        lamba[i] = sq_res.r1;
        x1x2_list[i] = sq_res.r2;

        CNOT_n(x2_list[i], x1x2_list[i]); 

        size_t start_mul = gm.current_pointer();
        count = 0;
        std::vector<int> sub_anc = slice(ancilla, num * i, num);
        auto res_mul = recursive_karatsuba(x1x2_list[i], lamba[i], n, count, sub_anc);
        std::vector<int> result2 = Reduction(res_mul.res);
        size_t end_mul = gm.current_pointer();

        std::vector<int> target_y = slice(mult_res, i * n, n);
        CNOT_n(result2, target_y);
        gm.replay_reverse(start_mul, end_mul);

        CNOT_n(y2_list[i], target_y);
        CNOT_n(x2_list[i], x1x2_list[i]); 
        CNOT_n(x1x2_list[i], target_y);
    }
}

ResourceStats QuantumContext::shor_reduction(int parallel_num, int level, int &length_out)
{
    gm.clear();
    cnt = 0;
    int n = current_config.N;
    long long nums1 = (long long)current_config.KARATSUBA_ANC_BASE * parallel_num;
    long long nums2 = (long long)current_config.TOFFOLI_BASE * (3 * (parallel_num - 1) + current_config.TOFFOLI_OFFSET);

    int idx = 0;
    std::vector<int> x1 = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> x2 = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> x1x2 = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> y1 = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> y2 = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;

    std::vector<int> ancilla = range_vec(idx, (int)nums1);
    idx += (int)nums1;
    Toffoli_qubits = range_vec(idx, (int)nums2);
    idx += (int)nums2;
    std::vector<int> sqr = range_vec(idx, 4 * n);
    idx += 4 * n;
    std::vector<int> ancilla2 = range_vec(idx, n * (parallel_num / 2));
    idx += n * (parallel_num / 2);

    std::vector<int> lambda_res = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> mult_res = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;
    std::vector<int> inv_res = range_vec(idx, n * parallel_num);
    idx += n * parallel_num;

    length_out = idx;

    std::vector<std::vector<int>> x1_list(parallel_num), y1_list(parallel_num), x1x2_list(parallel_num), x2_list(parallel_num), y2_list(parallel_num);
    for (int i = 0; i < parallel_num; ++i)
    {
        x1_list[i] = slice(x1, i * n, n);
        y1_list[i] = slice(y1, i * n, n);
        x1x2_list[i] = slice(x1x2, i * n, n);
        x2_list[i] = slice(x2, i * n, n);
        y2_list[i] = slice(y2, i * n, n);
    }

    shor_reduction_level(parallel_num, x1_list, y1_list, x2_list, y2_list, ancilla, x1x2_list, ancilla2, sqr, lambda_res, mult_res, inv_res);
    auto [fd, cd, td, t_c, c_c] = gm.optimize_and_count(length_out);
    return ResourceStats{fd, cd, td, t_c, c_c};
}