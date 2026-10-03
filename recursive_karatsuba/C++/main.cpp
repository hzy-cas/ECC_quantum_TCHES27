// main.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>

#include "QuantumLib.h"
#include "config.h"

GFConfig current_config;

struct Result {
    int input_x;           
    long long toffoli;     
    long long cnot;        
    long long full_depth;     
    long long current_depth;
    long long toffoli_depth; 
    long long qubits;      
    double time_sec;       
};

Result run_accumulate_and_reduce(int total_inputs, int max_parallel_ops) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    long long total_Toffoli = 0;
    long long total_CNOT = 0;
    long long total_Full_Depth = 0;
    long long total_Current_Depth = 0;
    long long total_Toffoli_Depth = 0;
    long long peak_qubits = 0;
    long long clear = 0; 
    
    long long n = current_config.N; 
    
    int current_pool = total_inputs;
    
    QuantumContext ctx; 
    int length = 0;

    // ========================================
    // Phase 1: Initialization (First Layer)
    // ========================================
    int ops = max_parallel_ops;
    
    ResourceStats s1 = ctx.shor(ops, 0, length);
    
    long long fanout_depth = 2LL * static_cast<long long>(std::ceil(std::log2(2.0 * n)));
    
    long long phase1_f_depth = fanout_depth + 1LL + s1.full_depth;
    long long phase1_c_depth = fanout_depth + 1LL + s1.current_depth;
    long long phase1_t_depth = s1.toffoli_depth; 
    
    long long phase1_cnot = 2LL * ops * (2LL * n - 1LL) + 2LL * n * ops + s1.cnot;
    
    total_Toffoli += s1.toffoli;
    total_CNOT += phase1_cnot;
    total_Full_Depth += phase1_f_depth;
    total_Current_Depth += phase1_c_depth;
    total_Toffoli_Depth += phase1_t_depth;
    
    peak_qubits = length - ops;
    clear = length - 6LL * n * ops - ops;

    current_pool -= (2 * ops); 

    // ========================================
    // Phase 2: Accumulation Loop
    // ========================================
    while (current_pool > 0) {
        int ops_prime = std::min(max_parallel_ops, current_pool);
        
        ResourceStats s_batch = ctx.shor(ops_prime, 0, length);
        
        total_Toffoli += s_batch.toffoli;
        total_CNOT += s_batch.cnot;
        total_Full_Depth += s_batch.full_depth;
        total_Current_Depth += s_batch.current_depth;
        total_Toffoli_Depth += s_batch.toffoli_depth;
        
        long long demand = length - 2LL * n * ops_prime - ops_prime;
        
        long long needed_new = 0;
        if (demand > clear) {
            needed_new = demand - clear;
            peak_qubits += needed_new;
            clear = length - 6LL * n * ops_prime - ops_prime;
        } else {
            long long over = clear - demand;
            clear = length - 6LL * n * ops_prime - ops_prime + over;
        }
        
        current_pool -= ops_prime;
    }

    // ========================================
    // Phase 3: Tree Reduction
    // ========================================
    int current_items = max_parallel_ops;
    
    while (current_items > 1) {
        int ops_double_prime = current_items / 2;
        int leftover = current_items % 2;
        
        ResourceStats s_tree = ctx.shor_reduction(ops_double_prime, 0, length);
        
        total_Toffoli += s_tree.toffoli;
        total_CNOT += s_tree.cnot;
        total_Full_Depth += s_tree.full_depth;
        total_Current_Depth += s_tree.current_depth;
        total_Toffoli_Depth += s_tree.toffoli_depth;
        
        long long demand = length - 4LL * n * ops_double_prime;
        
        long long needed_new = 0;
        if (demand > clear) {
            needed_new = demand - clear;
            peak_qubits += needed_new;
            clear = length - 8LL * n * ops_double_prime;
        } else {
            long long over = clear - demand;
            clear = length - 8LL * n * ops_double_prime + over;
        }
        
        current_items = ops_double_prime + leftover;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;

    return {
        max_parallel_ops,                                      
        2LL * total_Toffoli,                                     
        2LL * total_CNOT + 2LL * n,   
        2LL * total_Full_Depth + 2LL,  
        2LL * total_Current_Depth + 2LL,
        2LL * total_Toffoli_Depth,
        peak_qubits + 4LL * n + 2LL,                              
        duration.count()                                       
    };
}

int main(int argc, char** argv) {
    try {
        int n = 163, w_start = 1, w_end = 1;
        std::string output, data;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "Usage: karatsuba_balanced_estimator [--n N] [--w W | --w-start W --w-end W] [--output CSV] [--data-root PATH]\n";
                return 0;
            }
            if (i + 1 >= argc) throw std::invalid_argument("missing value for " + arg);
            std::string value = argv[++i];
            if (arg == "--n") n = std::stoi(value);
            else if (arg == "--w") w_start = w_end = std::stoi(value);
            else if (arg == "--w-start") w_start = std::stoi(value);
            else if (arg == "--w-end") w_end = std::stoi(value);
            else if (arg == "--output") output = value;
            else if (arg == "--data-root") data = value;
            else throw std::invalid_argument("unknown option: " + arg);
        }
        setup_config(n);
        if (w_start < 1 || w_end < w_start || w_end > n + 1)
            throw std::invalid_argument("require 1 <= w-start <= w-end <= n+1");
        if (data.empty()) data = KARATSUBA_SOURCE_DIR;
        current_config.MATRIX_PATH = data + "/square_" + std::to_string(n) + "/";
        if (!std::filesystem::is_directory(current_config.MATRIX_PATH))
            throw std::invalid_argument("missing squaring data for requested field");
        std::ofstream file;
        if (!output.empty()) {
            file.open(output);
            if (!file) throw std::runtime_error("cannot open output CSV");
        }
        std::ostream& csv = output.empty() ? std::cout : file;
        csv << "Parallel_Ops,Toffoli,CNOT,Full_Depth,Current_Depth,Toffoli_Depth,Qubits,DW_Cost_Full,DW_Cost_Current,Time_Seconds\n";
        for (int w = w_start; w <= w_end; ++w) {
            const auto r = run_accumulate_and_reduce(2 * n + 2, w);
            csv << r.input_x << ',' << r.toffoli << ',' << r.cnot << ','
                << r.full_depth << ',' << r.current_depth << ',' << r.toffoli_depth
                << ',' << r.qubits << ',' << r.full_depth * r.qubits << ','
                << r.current_depth * r.qubits << ',' << r.time_sec << '\n';
            csv.flush();
        }
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
}
