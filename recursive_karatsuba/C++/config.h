// config.h
#ifndef CONFIG_H
#define CONFIG_H

#include <cstdint>
#include <string>
#include <stdexcept>

constexpr int SHIFT_OP = 60;
constexpr int SHIFT_ARG = 30;
constexpr uint64_t MASK_ARG = (1ULL << SHIFT_ARG) - 1;

constexpr uint32_t OP_DATA = 0;
constexpr uint32_t OP_X = 1;
constexpr uint32_t OP_CNOT = 2;
constexpr uint32_t OP_TOFFOLI = 3;
constexpr uint32_t OP_NOP = 15;

struct GFConfig {
    int N;
    int KARATSUBA_ANC_BASE;
    int TOFFOLI_BASE;
    int TOFFOLI_OFFSET;
    std::string MATRIX_PATH;
};

extern GFConfig current_config;

inline void setup_config(int n) {
    if (n == 163) {
        current_config = {163, 8448, 4387, 9, "./square_163/"};
    } else if (n == 233) {
        current_config = {233, 12180, 6323, 10, "./square_233/"};
    } else if (n == 283) {
        current_config = {283, 19980, 10273, 11, "./square_283/"}; 
    } else if (n == 571) {
        current_config = {571, 61200, 31171, 13, "./square_571/"};
    } else {
        throw std::runtime_error("Unsupported GF(N). Options: 163, 233, 283, 571");
    }
}

#endif // CONFIG_H