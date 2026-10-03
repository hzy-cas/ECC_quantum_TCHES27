#pragma once
#include <string>
#include <iostream>

struct GFConfig {
    int n;
    int block_size;
    int nums2_const;
    int mul_cost_size;
    std::string data_root;

    
    static GFConfig setup(int target_n, std::string root_path) {
        GFConfig cfg;
        cfg.n = target_n;
        cfg.data_root = root_path;
        
        if (target_n == 163) {
            cfg.block_size = 743;
            cfg.nums2_const = 9;
            cfg.mul_cost_size = 906;
        } else if (target_n == 233) {
            cfg.block_size = 1108;
            cfg.nums2_const = 10;
            cfg.mul_cost_size = 1341;
        } else if (target_n == 283) {
            cfg.block_size = 1385;
            cfg.nums2_const = 11;
            cfg.mul_cost_size = 1668;
        } else if (target_n == 571) {
            cfg.block_size = 2998;
            cfg.nums2_const = 13;
            cfg.mul_cost_size = 3569;
        } else {
            std::cerr << "Unsupported N: " << target_n << std::endl;
            exit(1);
        }
        return cfg;
    }
};