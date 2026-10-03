#include "QuantumLib.h"
#include <iostream>
#include <map>
#include "config.h"
extern GFConfig current_config;

void QuantumContext::Modular_small(const std::vector<int>& input, const std::vector<int>& result, int size) {
    int n = current_config.N;
    
    if (n == 163) {
        // x^163 + x^7 + x^6 + x^3 + 1
        for(int i=0; i<size; ++i) CNOT(input[i], result[0+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[3+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[6+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[7+i]);
    } 
    else if (n == 233) {
        // x^233 + x^74 + 1
        for(int i=0; i<size; ++i) CNOT(input[i], result[0+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[74+i]);
    }
    else if (n == 283) {
        // x^283 + x^12 + x^7 + x^5 + 1
        for(int i=0; i<size; ++i) CNOT(input[i], result[0+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[5+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[7+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[12+i]);
    }
    else if (n == 571) {
        // x^571 + x^10 + x^5 + x^2 + 1
        for(int i=0; i<size; ++i) CNOT(input[i], result[0+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[2+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[5+i]);
        for(int i=0; i<size; ++i) CNOT(input[i], result[10+i]);
    }
}

std::vector<int> QuantumContext::Reduction(const std::vector<int>& res_in) {
    std::vector<int> result = res_in; 
    int n = current_config.N;

    if (n == 163) {
        for(int i=0; i<n-1; ++i) if(i < n) CNOT(result[i+n], result[i]);
        for(int i=0; i<n-1; ++i) if(i+3 < n) CNOT(result[i+n], result[i+3]);
        for(int i=0; i<n-1; ++i) if(i+6 < n) CNOT(result[i+n], result[i+6]);
        for(int i=0; i<n-1; ++i) if(i+7 < n) CNOT(result[i+n], result[i+7]);
        
        Modular_small(slice(result, 323, result.size()-323), result, 2);
        Modular_small(slice(result, 320, result.size()-320), result, 5);
        Modular_small(slice(result, 319, result.size()-319), result, 6);
    } 
    else if (n == 233) {
        for(int i=0; i<n-1; ++i) if(i < n) CNOT(result[i+n], result[i]);
        for(int i=0; i<n-1; ++i) if(i+74 < n) CNOT(result[i+n], result[i+74]);
        
        int start_idx = 392; 
        if (start_idx < result.size()) {
             Modular_small(slice(result, start_idx, result.size()-start_idx), result, 73);
        }
    }
    else if (n == 283) {
        for(int i=0; i<n-1; ++i) if(i < n) CNOT(result[i+n], result[i]);
        for(int i=0; i<n-1; ++i) if(i+5 < n) CNOT(result[i+n], result[i+5]);
        for(int i=0; i<n-1; ++i) if(i+7 < n) CNOT(result[i+n], result[i+7]);
        for(int i=0; i<n-1; ++i) if(i+12 < n) CNOT(result[i+n], result[i+12]);
        
        Modular_small(slice(result, 561, result.size()-561), result, 4);
        Modular_small(slice(result, 559, result.size()-559), result, 6);
        Modular_small(slice(result, 554, result.size()-554), result, 11);
    }
    else if (n == 571) {
        for(int i=0; i<n-1; ++i) if(i < n) CNOT(result[i+n], result[i]);
        for(int i=0; i<n-1; ++i) if(i+2 < n) CNOT(result[i+n], result[i+2]);
        for(int i=0; i<n-1; ++i) if(i+5 < n) CNOT(result[i+n], result[i+5]);
        for(int i=0; i<n-1; ++i) if(i+10 < n) CNOT(result[i+n], result[i+10]);

        Modular_small(slice(result, 1140, result.size()-1140), result, 1);
        Modular_small(slice(result, 1137, result.size()-1137), result, 4);
        Modular_small(slice(result, 1132, result.size()-1132), result, 9);
    }
    return slice(result, 0, n);
}

std::vector<int> QuantumContext::combine(const std::vector<int>& a, const std::vector<int>& b, const std::vector<int>& r, int n) {
    if (n % 2 != 0) {
        for(int i=0; i<n; ++i) CNOT(a[i], r[i]);
        for(int i=0; i<n-2; ++i) CNOT(b[i], r[i]);
        for(int i=0; i<n/2; ++i) CNOT(a[n/2+1+i], r[i]);
        for(int i=0; i<n/2; ++i) CNOT(b[i], r[n/2+1+i]);
        std::vector<int> out;
        for(int i=0; i<n/2+1; ++i) out.push_back(a[i]);
        for(int i=0; i<n; ++i) out.push_back(r[i]);
        int b_start = n/2;
        int count = (2*n - 1) - n/2 - 1 - n;
        for(int i=0; i<count; ++i) out.push_back(b[b_start+i]);
        return out;
    }
    int half_n = n/2;
    for(int i=0; i<n-1; ++i) { CNOT(a[i], r[i]); CNOT(b[i], r[i]); }
    for(int i=0; i<half_n-1; ++i) { CNOT(a[half_n+i], r[i]); CNOT(b[i], r[half_n+i]); }
    std::vector<int> result;
    for(int i=0; i<half_n; ++i) result.push_back(a[i]);
    for(int i=0; i<n-1; ++i) result.push_back(r[i]);
    for(int i=0; i<half_n; ++i) result.push_back(b[half_n-1+i]);
    return result;
}

QuantumContext::KaratsubaResult QuantumContext::recursive_karatsuba(std::vector<int> a, std::vector<int> b, int n, int count, std::vector<int> ancilla) {
    if (n == 1) {
        int c = Toffoli_qubits[cnt];
        Toffoli_gate(a[0], b[0], c);
        cnt++;
        return {{c}, count, ancilla};
    }
    int r_low = n/2;
    if (n % 2 != 0) r_low += 1;
    
    std::vector<int> r_a = slice(ancilla, count, r_low); count += r_low;
    std::vector<int> r_b = slice(ancilla, count, r_low); count += r_low;
    
    size_t start_ptr = gm.current_pointer();
    
    for(int i=0; i<r_low; ++i) CNOT(a[i], r_a[i]);
    for(int i=0; i<n/2; ++i) CNOT(a[r_low+i], r_a[i]);
    for(int i=0; i<r_low; ++i) CNOT(b[i], r_b[i]);
    for(int i=0; i<n/2; ++i) CNOT(b[r_low+i], r_b[i]);
    
    size_t end_ptr = gm.current_pointer();
    
    if (r_low == 1) {
        std::vector<int> c = slice(Toffoli_qubits, cnt, 3);
        cnt += 3;
        Toffoli_gate(a[0], b[0], c[0]);
        Toffoli_gate(a[1], b[1], c[2]);
        Toffoli_gate(r_a[0], r_b[0], c[1]);
        CNOT(c[0], c[1]);
        CNOT(c[2], c[1]);
        gm.replay_reverse(start_ptr, end_ptr);
        return {c, count, ancilla};
    }
    
    auto res_ca = recursive_karatsuba(slice(a, 0, r_low), slice(b, 0, r_low), r_low, count, ancilla);
    count = res_ca.count; ancilla = res_ca.ancilla;
    
    auto res_cb = recursive_karatsuba(slice(a, r_low, n-r_low), slice(b, r_low, n-r_low), n/2, count, ancilla);
    count = res_cb.count; ancilla = res_cb.ancilla;
    
    auto res_cr = recursive_karatsuba(slice(r_a, 0, r_low), slice(r_b, 0, r_low), r_low, count, ancilla);
    count = res_cr.count; ancilla = res_cr.ancilla;
    
    gm.replay_reverse(start_ptr, end_ptr);
    
    std::vector<int> final_res = combine(res_ca.res, res_cb.res, res_cr.res, n);
    return {final_res, count, ancilla};
}

QuantumContext::SquareResult QuantumContext::Squaring_plus_input(std::vector<int> x, int n) {
    std::string path = current_config.MATRIX_PATH + "square_Matrix_Squaring_plus.txt";
    MatrixData data = MatrixLoader::get_matrix_data(path, n);
    
    std::map<int, int> y;
    for (auto& op : data.ops) {
        if (op.control < x.size() && op.target < x.size()) {
            CNOT(x[op.control], x[op.target]);
            if (op.y_index != -1) y[op.y_index] = x[op.target];
        }
    }
    if (y.empty()) return {slice(x, 0, n), slice(x, 0, n)};
    
    int max_y = 0;
    if (!y.empty()) max_y = y.rbegin()->first;
    std::vector<int> y_list(max_y + 1, 0);
    
    for (auto const& [idx, val] : y) y_list[idx] = val;
    for (int i=0; i<y_list.size(); ++i) {
        if (y_list[i] == 0 && y.find(i) == y.end() && i < x.size()) {
            y_list[i] = x[i];
        }
    }
    while(y_list.size() < 2*n) y_list.push_back(0);
    return {slice(y_list, 0, n), slice(y_list, y_list.size()-n, n)};
}

QuantumContext::SquareResult QuantumContext::Square(std::vector<int> x, int n, int power) {
    std::string path = current_config.MATRIX_PATH + "square_Matrix_2_" + std::to_string(power) + ".txt";
    MatrixData data = MatrixLoader::get_matrix_data(path, n);
    
    std::map<int, int> y;
    for (auto& op : data.ops) {
        if (op.control < x.size() && op.target < x.size()) {
            CNOT(x[op.control], x[op.target]);
            if (op.y_index != -1) y[op.y_index] = x[op.target];
        }
    }
    if (y.empty()) return {slice(x, 0, n), slice(x, 0, n)};

    int max_y = 0;
    if (!y.empty()) max_y = y.rbegin()->first;
    std::vector<int> y_list(max_y + 1, 0);
    
    for (auto const& [idx, val] : y) y_list[idx] = val;
    for (int i=0; i<y_list.size(); ++i) {
        if (y_list[i] == 0 && y.find(i) == y.end() && i < x.size()) {
            y_list[i] = x[i];
        }
    }
    while(y_list.size() < 2*n) y_list.push_back(0);
    return {slice(y_list, 0, n), slice(y_list, y_list.size()-n, n)};
}