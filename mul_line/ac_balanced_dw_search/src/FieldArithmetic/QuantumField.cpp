#include "../BasicLogic/QuantumLib.h"
#include <map>
#include <stdexcept>

std::vector<int> QuantumContext::mul_matrix(const std::vector<int>& x, int n_rows, std::string file_path) {
    auto [ops, matrix_rows] = matrix_loader.get_matrix_data(file_path, n_rows);
    std::map<int, int> y;
    
    for(const auto& op : ops) {
        gm.add_CNOT(x[op.control], x[op.target]);
        if(op.y_index != -1) y[op.y_index] = x[op.target];
    }
    
    if(y.empty()) return {};
    
    std::vector<int> y_list(n_rows, 0);
    for(auto const& [key, val] : y) y_list[key] = val;
    
    for(std::size_t i=0; i<y_list.size(); ++i) {
        if(y_list[i] == 0 && y.find(i) == y.end()) {
             if(i < matrix_rows.size() && !matrix_rows[i].empty()) {
                 int j = matrix_rows[i][0];
                 if (j >= 0 && static_cast<std::size_t>(j) < x.size()) y_list[i] = x[static_cast<std::size_t>(j)];
             }
        }
    }
    return y_list;
}

std::pair<std::vector<int>, std::vector<int>> QuantumContext::Square(const std::vector<int>& x, int n_rows, int power) {
    std::string path = get_matrix_path("square_" + std::to_string(config.n), "square_Matrix_2_" + std::to_string(power) + ".txt");
    auto y_list = mul_matrix(x, n_rows, path);
    int half = n_rows / 2;
    return {std::vector<int>(y_list.begin(), y_list.begin() + half), 
            std::vector<int>(y_list.end() - half, y_list.end())};
}

std::pair<std::vector<int>, std::vector<int>> QuantumContext::Squaring_plus_input(const std::vector<int>& x, int n_rows) {
    std::string path = get_matrix_path("square_" + std::to_string(config.n), "square_Matrix_Squaring_plus.txt");
    auto y_list = mul_matrix(x, n_rows, path);
    int half = n_rows / 2;
    return {std::vector<int>(y_list.begin(), y_list.begin() + half), 
            std::vector<int>(y_list.end() - half, y_list.end())};
}

std::vector<int> QuantumContext::mul(const std::vector<int>& a0, const std::vector<int>& b0, const std::vector<int>& c0, int offset) {
    (void)offset;
    // In C++, offset identifies the start within Toffoli_qubits rather than passing a slice.
    // For a compact interface, slice c0 internally when the caller supplied a subvector.
    // Strengthened bounds check.
    const auto mul_size = static_cast<std::size_t>(config.mul_cost_size);
    if (c0.size() < mul_size) {
        throw std::invalid_argument("multiplication target register too short");
    }
    
    if (a0.size() < mul_size || b0.size() < mul_size) {
        throw std::invalid_argument("multiplication input register too short");
    }
    
    std::string path_td = get_matrix_path("CNOT_mul", "result_" + std::to_string(config.n) + "_TD_seq.txt");
    std::string path_a = get_matrix_path("CNOT_mul", "result_" + std::to_string(config.n) + "_A.txt");
    std::string path_c = get_matrix_path("CNOT_mul", "result_" + std::to_string(config.n) + "_C.txt");
    std::string path_inv = get_matrix_path("CNOT_mul", "result_" + std::to_string(config.n) + "_CT2D_inv_seq.txt");
    
    int slice1 = 2 * config.n - 1;
    int total_mul_size = config.mul_cost_size;
    
    size_t begin = gm.current_pointer();
    
    // a0[:slice1] construct
    std::vector<int> a0_sub(a0.begin(), a0.begin() + slice1);
    std::vector<int> b0_sub(b0.begin(), b0.begin() + slice1);
    
    auto a1 = mul_matrix(a0_sub, slice1, path_td);
    auto b1 = mul_matrix(b0_sub, slice1, path_td);
    
    // a1 + a0[slice1:total]
    std::vector<int> a_combined = a1;
    a_combined.insert(a_combined.end(), a0.begin() + slice1, a0.begin() + total_mul_size);
    std::vector<int> b_combined = b1;
    b_combined.insert(b_combined.end(), b0.begin() + slice1, b0.begin() + total_mul_size);
    
    auto a2 = mul_matrix(a_combined, total_mul_size, path_a);
    auto b2 = mul_matrix(b_combined, total_mul_size, path_a);
    
    size_t end_ptr = gm.current_pointer();
    
    for(int i=0; i<total_mul_size; ++i) {
        gm.add_Toffoli(a2[i], b2[i], c0[i]);
    }
    
    auto c1 = mul_matrix(std::vector<int>(c0.begin(), c0.begin() + total_mul_size), total_mul_size, path_c);
    auto c2 = mul_matrix(std::vector<int>(c1.begin(), c1.begin() + slice1), slice1, path_inv);
    
    gm.replay_reverse(begin, end_ptr);
    
    return std::vector<int>(c2.begin(), c2.begin() + config.n);
}
