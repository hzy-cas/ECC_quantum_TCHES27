// src/Infrastructure/MatrixLoader.h
#pragma once
#include <string>
#include <vector>
#include <map>
#include <tuple>

struct MatrixOp {
    int target;
    int control;
    int y_index; // -1 if None
};

class MatrixCache {
public:
    // Return (ops, matrix_rows).
    std::pair<std::vector<MatrixOp>, std::vector<std::vector<int>>> 
    get_matrix_data(const std::string& file_path, int n);

private:
    std::map<std::string, std::pair<std::vector<MatrixOp>, std::vector<std::vector<int>>>> cache;
};
