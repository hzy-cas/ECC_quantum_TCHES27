#ifndef MATRIXLOADER_H
#define MATRIXLOADER_H

#include <vector>
#include <string>
#include <map>

struct MatrixOp {
    int target;
    int control;
    int y_index; 
};

struct MatrixData {
    std::vector<MatrixOp> ops;
    std::map<int, int> y_map;
    bool valid;
};

class MatrixLoader {
public:
    static MatrixData get_matrix_data(const std::string& filename, int n);
};

#endif