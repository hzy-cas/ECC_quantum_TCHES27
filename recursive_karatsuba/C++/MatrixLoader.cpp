#include "MatrixLoader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>
#include <mutex> 

static std::map<std::string, MatrixData> cache;
static std::mutex cache_mutex;

MatrixData MatrixLoader::get_matrix_data(const std::string& filename, int n) {
    std::lock_guard<std::mutex> lock(cache_mutex);

    if (cache.count(filename)) return cache[filename];

    MatrixData data;
    data.valid = false;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        cache[filename] = data; 
        return data;
    }

    data.valid = true;
    std::string line;
    std::vector<std::string> lines;
    while(std::getline(file, line)) lines.push_back(line);

    int start_index = n;
    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].find("CNOT Depth") != std::string::npos) {
            start_index = i + 1;
            break;
        }
    }

    std::regex re(R"(x\[(\d+)\]\s*=\s*x\[\d+\]\s*\^\s*x\[(\d+)\](\s*y\[(\d+)\])?)");
    std::smatch match;

    for (size_t i = start_index; i < lines.size(); ++i) {
        if (lines[i].empty()) continue;
        if (std::regex_search(lines[i], match, re)) {
            MatrixOp op;
            op.target = std::stoi(match[1]);
            op.control = std::stoi(match[2]);
            op.y_index = -1;
            if (match[4].matched) {
                op.y_index = std::stoi(match[4]);
                data.y_map[op.y_index] = 1;
            }
            data.ops.push_back(op);
        }
    }
    
    cache[filename] = data;
    return data;
}