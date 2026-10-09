#include "MatrixLoader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>

// Helper: trim leading and trailing whitespace.
static std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::pair<std::vector<MatrixOp>, std::vector<std::vector<int>>>
MatrixCache::get_matrix_data(const std::string &file_path, int n)
{
    if (cache.count(file_path))
        return cache[file_path];

    std::vector<MatrixOp> ops;
    std::vector<std::vector<int>> matrix_rows;

    std::ifstream infile(file_path);
    if (!infile.is_open())
    {
        std::cerr << "[Error] Cannot open matrix file: " << file_path << std::endl;
        std::cerr << "Please ensure the 'data' folder is in the correct path." << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<std::string> lines;
    while (std::getline(infile, line))
        lines.push_back(line);
    infile.close();

    // =========================================================
    // Step 1: parse the leading binary matrix (the first n rows).
    // =========================================================
    int line_count = 0;
    for (const auto &l : lines)
    {
        if (line_count >= n) break;
        
        std::string stripped = trim(l);
        
        // Skip empty and marker lines.
        if (stripped.empty() || 
            stripped.find("Original") != std::string::npos || 
            stripped.find("Reduced") != std::string::npos)
            continue;

        // Extract the indices of all '1' entries.
        std::vector<int> ones;
        for (size_t j = 0; j < stripped.size(); ++j)
        {
            if (stripped[j] == '1')
                ones.push_back(static_cast<int>(j));
        }
        
        if (!ones.empty())
        {
            matrix_rows.push_back(ones);
            line_count++;
        }
    }

    // =========================================================
    // Step 2: find the "CNOT Depth" marker and locate the first operation.
    // =========================================================
    size_t start_index = n;
    for (size_t i = 0; i < lines.size(); ++i)
    {
        if (lines[i].find("CNOT Depth") != std::string::npos)
        {
            start_index = i + 1;
            break;
        }
    }

    // =========================================================
    // Step 3: parse CNOT operations: x[target]=x[target]^x[control] [y[index]].
    // =========================================================
    for (size_t i = start_index; i < lines.size(); ++i)
    {
        std::string l = trim(lines[i]);
        if (l.empty() || l.find('=') == std::string::npos)
            continue;

        try {
            // Split the left- and right-hand sides.
            size_t eq_pos = l.find('=');
            std::string left_part = l.substr(0, eq_pos);
            std::string right_part = l.substr(eq_pos + 1);

            // Parse target: x[target].
            size_t br1 = left_part.find('[');
            size_t br2 = left_part.find(']');
            
            if (br1 == std::string::npos || br2 == std::string::npos || br2 <= br1 + 1) {
                continue;
            }
            
            std::string target_str = trim(left_part.substr(br1 + 1, br2 - br1 - 1));
            if (target_str.empty()) continue;
            int target = std::stoi(target_str);

            // Parse control: x[control], following the ^ symbol.
            size_t xor_pos = right_part.find('^');
            if (xor_pos == std::string::npos) continue;

            std::string after_xor = right_part.substr(xor_pos + 1);
            size_t br3 = after_xor.find('[');
            size_t br4 = after_xor.find(']');
            
            if (br3 == std::string::npos || br4 == std::string::npos || br4 <= br3 + 1) {
                continue;
            }
            
            std::string control_str = trim(after_xor.substr(br3 + 1, br4 - br3 - 1));
            if (control_str.empty()) continue;
            int control = std::stoi(control_str);

            // Parse the optional y[index], when present.
            int y_idx = -1;
            size_t y_pos = l.find("y[");
            if (y_pos != std::string::npos)
            {
                size_t br_y = l.find(']', y_pos);
                if (br_y != std::string::npos && br_y > y_pos + 2) {
                    std::string y_str = trim(l.substr(y_pos + 2, br_y - (y_pos + 2)));
                    if (!y_str.empty()) {
                        y_idx = std::stoi(y_str);
                    }
                }
            }

            ops.push_back({target, control, y_idx});
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "[Warning] Invalid number in line: " << l << std::endl;
            std::cerr << "  Exception: " << e.what() << std::endl;
            continue;
        }
        catch (const std::out_of_range& e) {
            std::cerr << "[Warning] Number out of range in line: " << l << std::endl;
            std::cerr << "  Exception: " << e.what() << std::endl;
            continue;
        }
        catch (const std::exception& e) {
            std::cerr << "[Warning] Failed to parse line: " << l << std::endl;
            std::cerr << "  Exception: " << e.what() << std::endl;
            continue;
        }
    }

    auto res = std::make_pair(ops, matrix_rows);
    cache[file_path] = res;
    return res;
}
