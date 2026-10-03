#include "matrix.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

std::vector<ROW> get_matrix(const std::string& filename, int& out_width)
{
    std::vector<ROW> m;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open matrix file: " << filename << std::endl;
        return m;
    }
    
    std::string line;
    int actual_size = 0;
    bool first_line = true;
    
    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        
        if (!line.empty() && line[line.length()-1] == '\r')
            line.erase(line.length()-1);

        if (first_line) {
            actual_size = line.length();
        
            out_width = actual_size; 

            if (actual_size > SIZE) {
                 std::cout << "Error: Matrix size " << actual_size 
                 << " exceeds defined SIZE " << SIZE << std::endl;
                 m.clear();
                 break;
            }
            first_line = false;
        }
        
        if (line.length() == actual_size) {
            ROW row(line); 
            m.push_back(row);
        } 
    }
    
    file.close();
    return m;
}

vector<ROW> get_matrix() {
    int dummy = 0;
    return get_matrix("a.txt", dummy); 
}