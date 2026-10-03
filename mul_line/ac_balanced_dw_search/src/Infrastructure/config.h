#pragma once

#include <stdexcept>
#include <string>

// Fixed parameters for AC (algebraic-curve-based) field arithmetic.
// block_size = l - n, where l is the bilinear multiplication complexity; mul_cost_size = l.
struct GFConfig {
    int n{};
    int block_size{};
    int inversion_mul_blocks{};
    int mul_cost_size{};
    std::string data_root;

    int total_shor_inputs() const { return 2 * n + 2; }
    int scan_w_max() const { return n == 571 ? 256 : 128; }

    static GFConfig setup(int target_n, const std::string& root_path) {
        GFConfig cfg;
        cfg.n = target_n;
        cfg.data_root = root_path;

        switch (target_n) {
        case 163:
            cfg.block_size = 743;
            cfg.inversion_mul_blocks = 9;
            cfg.mul_cost_size = 906;
            break;
        case 233:
            cfg.block_size = 1108;
            cfg.inversion_mul_blocks = 10;
            cfg.mul_cost_size = 1341;
            break;
        case 283:
            cfg.block_size = 1385;
            cfg.inversion_mul_blocks = 11;
            cfg.mul_cost_size = 1668;
            break;
        case 571:
            cfg.block_size = 2998;
            cfg.inversion_mul_blocks = 13;
            cfg.mul_cost_size = 3569;
            break;
        default:
            throw std::invalid_argument("unsupported n; choose 163, 233, 283, or 571");
        }
        return cfg;
    }
};
