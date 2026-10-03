#include "Infrastructure/config.h"
#include "BasicLogic/QuantumLib.h"
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        int n = 163;
        int inputs = 0;
        std::string data;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "Usage: ShorEstimator [--n 163|233|283|571] [--inputs N] [--data-root PATH]\n"
                          << "Default: N=2n+2. --inputs is for smaller schedule checks.\n";
                return 0;
            }
            if (i + 1 >= argc) throw std::invalid_argument("missing value for " + arg);
            std::string value = argv[++i];
            if (arg == "--n") n = std::stoi(value);
            else if (arg == "--inputs") inputs = std::stoi(value);
            else if (arg == "--data-root") data = value;
            else throw std::invalid_argument("unknown option: " + arg);
        }
        if (data.empty()) {
#ifdef AC_SOURCE_DATA_DIR
            data = AC_SOURCE_DATA_DIR;
#else
            data = "data";
#endif
        }
        if (inputs == 0) inputs = 2 * n + 2;
        if (inputs < 2) throw std::invalid_argument("--inputs must be at least 2");
        if (!std::filesystem::is_directory(data + "/quantum_" + std::to_string(n)))
            throw std::invalid_argument("missing arithmetic data for requested field");
        GFConfig cfg = GFConfig::setup(n, data);
        QuantumContext ctx(cfg);
        ctx.estimate_total_resources(inputs);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
}
