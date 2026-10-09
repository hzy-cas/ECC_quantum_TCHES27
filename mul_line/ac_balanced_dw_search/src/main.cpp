#include "AlgorithmFlow/BalancedResource.h"
#include "BasicLogic/QuantumLib.h"
#include "Infrastructure/config.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kModelVersion = "ac-balanced-clean-nct-v2";

struct Options {
    std::string command = "help";
    int n = 163;
    int w = 1;
    int w_start = 1;
    int w_end = 1;
    int g = 1;
    int lanes = 1;
    std::string mode = "accumulation";
    std::string kind = "multiplication";
    fs::path data_root;
    fs::path output_csv;
    fs::path layer_cache;
    bool force = false;
};

fs::path discover_data_root(const char* executable) {
    std::vector<fs::path> search_roots{fs::current_path()};
#ifdef AC_BALANCED_SOURCE_DIR
    search_roots.emplace_back(AC_BALANCED_SOURCE_DIR);
#endif
    std::error_code error;

    const auto proc_executable = fs::read_symlink("/proc/self/exe", error);
    if (!error) search_roots.push_back(proc_executable.parent_path());

    error.clear();
    const auto argument_path = fs::weakly_canonical(fs::absolute(executable), error);
    if (!error) search_roots.push_back(argument_path.parent_path());

    for (auto root : search_roots) {
        while (!root.empty()) {
            const auto candidate = root / "mul_line" / "C++" / "data";
            if (fs::is_directory(candidate)) return fs::weakly_canonical(candidate);
            const auto parent = root.parent_path();
            if (parent == root) break;
            root = parent;
        }
    }

    throw std::runtime_error(
        "cannot locate mul_line/C++/data; pass --data-root PATH explicitly");
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    return fields;
}

int parse_int(const std::string& value, const char* option) {
    std::size_t used = 0;
    int parsed = 0;
    try {
        parsed = std::stoi(value, &used);
    } catch (...) {
        throw std::invalid_argument(std::string("invalid integer for ") + option + ": " + value);
    }
    if (used != value.size()) {
        throw std::invalid_argument(std::string("invalid integer for ") + option + ": " + value);
    }
    return parsed;
}

Options parse_options(int argc, char** argv) {
    Options options;
    if (argc >= 2) options.command = argv[1];

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--force") {
            options.force = true;
            continue;
        }
        if (i + 1 >= argc) throw std::invalid_argument("missing value after " + arg);
        const std::string value = argv[++i];
        if (arg == "--n") options.n = parse_int(value, "--n");
        else if (arg == "--w") options.w = parse_int(value, "--w");
        else if (arg == "--w-start") options.w_start = parse_int(value, "--w-start");
        else if (arg == "--w-end") options.w_end = parse_int(value, "--w-end");
        else if (arg == "--g") options.g = parse_int(value, "--g");
        else if (arg == "--lanes") options.lanes = parse_int(value, "--lanes");
        else if (arg == "--mode") options.mode = value;
        else if (arg == "--kind") options.kind = value;
        else if (arg == "--data-root") options.data_root = value;
        else if (arg == "--output-csv") options.output_csv = value;
        else if (arg == "--layer-cache") options.layer_cache = value;
        else throw std::invalid_argument("unknown option: " + arg);
    }
    return options;
}

PointAddMode parse_mode(const std::string& mode) {
    if (mode == "accumulation") return PointAddMode::Accumulation;
    if (mode == "reduction") return PointAddMode::Reduction;
    throw std::invalid_argument("--mode must be accumulation or reduction");
}

const char* mode_name(PointAddMode mode) {
    return mode == PointAddMode::Accumulation ? "accumulation" : "reduction";
}

void ensure_parent(const fs::path& path) {
    if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
}

void print_stats_json(const LayerStats& stats) {
    std::cout << "{\"qubits\":" << stats.qubits
              << ",\"toffoli\":" << stats.toffoli
              << ",\"cnot\":" << stats.cnot
              << ",\"full_depth\":" << stats.full_depth
              << ",\"current_depth\":" << stats.current_depth
              << ",\"toffoli_depth\":" << stats.toffoli_depth << "}\n";
}

struct LayerKey {
    PointAddMode mode{};
    int lanes{};
    bool operator<(const LayerKey& other) const {
        if (mode != other.mode) return static_cast<int>(mode) < static_cast<int>(other.mode);
        return lanes < other.lanes;
    }
};

class LayerRepository {
public:
    LayerRepository(const GFConfig& config, fs::path cache_path)
        : config_(config), context_(config), cache_path_(std::move(cache_path)) {
        load();
    }

    LayerStats get(int lanes, PointAddMode mode) {
        const LayerKey key{mode, lanes};
        const auto found = cache_.find(key);
        if (found != cache_.end()) return found->second;

        std::cerr << "[layer] constructing n=" << config_.n
                  << " mode=" << mode_name(mode) << " lanes=" << lanes << std::endl;
        const auto start = std::chrono::steady_clock::now();
        const auto stats = context_.run_balanced_layer(lanes, mode);
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        cache_[key] = stats;
        append(key, stats, seconds);
        std::cerr << "[layer] completed lanes=" << lanes
                  << " in " << std::fixed << std::setprecision(3) << seconds << " s"
                  << ", peak-layer-qubits=" << stats.qubits << std::endl;
        return stats;
    }

private:
    GFConfig config_;
    QuantumContext context_;
    fs::path cache_path_;
    std::map<LayerKey, LayerStats> cache_;

    void load() {
        if (cache_path_.empty() || !fs::exists(cache_path_)) return;
        std::ifstream input(cache_path_);
        std::string line;
        std::getline(input, line); // header
        while (std::getline(input, line)) {
            if (line.empty()) continue;
            const auto f = split_csv(line);
            if (f.size() < 12 || f[0] != kModelVersion) {
                throw std::runtime_error("incompatible layer cache; use a new path for the shared NCT backend");
            }
            if (parse_int(f[1], "cache n") != config_.n) continue;
            LayerKey key{parse_mode(f[2]), parse_int(f[3], "cache lanes")};
            LayerStats stats;
            stats.qubits = std::stoull(f[4]);
            stats.toffoli = std::stoull(f[5]);
            stats.cnot = std::stoull(f[6]);
            stats.full_depth = std::stoull(f[7]);
            stats.current_depth = std::stoull(f[8]);
            stats.toffoli_depth = std::stoull(f[9]);
            cache_[key] = stats;
        }
    }

    void append(const LayerKey& key, const LayerStats& stats, double seconds) {
        if (cache_path_.empty()) return;
        ensure_parent(cache_path_);
        const bool needs_header = !fs::exists(cache_path_) || fs::file_size(cache_path_) == 0;
        std::ofstream output(cache_path_, std::ios::app);
        if (!output) throw std::runtime_error("cannot open layer cache: " + cache_path_.string());
        if (needs_header) {
            output << "model_version,n,mode,lanes,qubits,toffoli,cnot,full_depth,"
                      "current_depth,toffoli_depth,seconds,data_root\n";
        }
        output << kModelVersion << ',' << config_.n << ',' << mode_name(key.mode) << ','
               << key.lanes << ',' << stats.qubits << ',' << stats.toffoli << ','
               << stats.cnot << ',' << stats.full_depth << ',' << stats.current_depth << ','
               << stats.toffoli_depth << ',' << std::setprecision(12) << seconds << ','
               << fs::absolute(config_.data_root).string() << '\n';
        output.flush();
    }
};

std::set<int> completed_widths(const fs::path& csv, int n) {
    std::set<int> widths;
    if (csv.empty() || !fs::exists(csv)) return widths;
    std::ifstream input(csv);
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto f = split_csv(line);
        if (f.size() < 4 || f[0] != kModelVersion) continue;
        if (parse_int(f[1], "result n") == n) widths.insert(parse_int(f[3], "result w"));
    }
    return widths;
}

void append_result(const fs::path& csv, const BalancedResult& result, double seconds) {
    ensure_parent(csv);
    const bool needs_header = !fs::exists(csv) || fs::file_size(csv) == 0;
    std::ofstream output(csv, std::ios::app);
    if (!output) throw std::runtime_error("cannot open result CSV: " + csv.string());
    if (needs_header) {
        output << "model_version,n,total_inputs,w,accumulation_layers,reduction_layers,"
                  "toffoli,cnot,full_depth,current_depth,toffoli_depth,qubits,"
                  "dw_full,dw_current,tdw,log2_dw_full,log2_dw_current,log2_tdw,seconds\n";
    }
    output << kModelVersion << ',' << result.n << ',' << result.total_inputs << ',' << result.w
           << ',' << result.accumulation_layers << ',' << result.reduction_layers
           << ',' << result.toffoli << ',' << result.cnot << ',' << result.full_depth
           << ',' << result.current_depth << ',' << result.toffoli_depth << ',' << result.qubits
           << ',' << uint128_to_string(result.dw_full)
           << ',' << uint128_to_string(result.dw_current)
           << ',' << uint128_to_string(result.tdw)
           << ',' << std::setprecision(15) << static_cast<double>(uint128_log2(result.dw_full))
           << ',' << static_cast<double>(uint128_log2(result.dw_current))
           << ',' << static_cast<double>(uint128_log2(result.tdw))
           << ',' << std::setprecision(12) << seconds << '\n';
    output.flush();
}

void validate_data_root(const GFConfig& config) {
    const fs::path field_root = fs::path(config.data_root) /
                                fs::path("quantum_" + std::to_string(config.n));
    if (!fs::is_directory(field_root)) {
        throw std::runtime_error("data directory not found: " + field_root.string());
    }
}

void validate_w_range(const GFConfig& config, int start, int end) {
    if (start < 1 || end < start || end > config.scan_w_max()) {
        throw std::invalid_argument("w range must satisfy 1 <= start <= end <= " +
                                    std::to_string(config.scan_w_max()));
    }
    if (2 * end > config.total_shor_inputs()) {
        throw std::invalid_argument("w is too large for the Balanced initialization schedule");
    }
}

void print_help() {
    std::cout
        << "AC arithmetic + Balanced strategy resource estimator\n\n"
        << "Commands:\n"
        << "  scan      --n N --w-start A --w-end B --data-root PATH\n"
        << "            --output-csv FILE --layer-cache FILE [--g G] [--force]\n"
        << "  estimate  --n N --w W --data-root PATH [--layer-cache FILE]\n"
        << "  layer     --n N --lanes K --mode accumulation|reduction --data-root PATH\n"
        << "  primitive --n N --kind multiplication|inversion --data-root PATH\n\n"
        << "When --data-root is omitted, the executable searches parent directories\n"
        << "for the repository's mul_line/C++/data directory.\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        auto options = parse_options(argc, argv);
        if (options.command == "help" || options.command == "--help" || options.command == "-h") {
            print_help();
            return 0;
        }

        if (options.data_root.empty()) options.data_root = discover_data_root(argv[0]);

        auto config = GFConfig::setup(options.n, fs::absolute(options.data_root).string());
        validate_data_root(config);

        if (options.command == "primitive") {
            QuantumContext context(config);
            if (options.kind == "multiplication") print_stats_json(context.estimate_multiplication());
            else if (options.kind == "inversion") print_stats_json(context.estimate_inversion());
            else throw std::invalid_argument("--kind must be multiplication or inversion");
            return 0;
        }

        if (options.command == "layer") {
            QuantumContext context(config);
            print_stats_json(context.run_balanced_layer(options.lanes, parse_mode(options.mode)));
            return 0;
        }

        if (options.command == "estimate") {
            validate_w_range(config, options.w, options.w);
            LayerRepository layers(config, options.layer_cache);
            const auto start = std::chrono::steady_clock::now();
            const auto result = estimate_balanced(
                config, options.w,
                [&](int lanes, PointAddMode mode) { return layers.get(lanes, mode); });
            const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            std::cout << "{\"n\":" << result.n << ",\"w\":" << result.w
                      << ",\"toffoli\":" << result.toffoli
                      << ",\"cnot\":" << result.cnot
                      << ",\"full_depth\":" << result.full_depth
                      << ",\"current_depth\":" << result.current_depth
                      << ",\"toffoli_depth\":" << result.toffoli_depth
                      << ",\"qubits\":" << result.qubits
                      << ",\"dw_full\":\"" << uint128_to_string(result.dw_full)
                      << "\",\"tdw\":\"" << uint128_to_string(result.tdw)
                      << "\",\"seconds\":" << seconds << "}\n";
            return 0;
        }

        if (options.command == "scan") {
            validate_w_range(config, options.w_start, options.w_end);
            if (options.g < 1) {
                throw std::invalid_argument("--g must be at least 1");
            }
            if (options.output_csv.empty() || options.layer_cache.empty()) {
                throw std::invalid_argument("scan requires --output-csv and --layer-cache");
            }

            LayerRepository layers(config, options.layer_cache);
            auto completed = completed_widths(options.output_csv, config.n);
            for (int w = options.w_start;;) {
                if (!options.force && completed.count(w) != 0) {
                    std::cerr << "[scan] skip existing n=" << config.n << " w=" << w << std::endl;
                } else {
                    std::cerr << "[scan] estimate n=" << config.n << " w=" << w << std::endl;
                    const auto start = std::chrono::steady_clock::now();
                    const auto result = estimate_balanced(
                        config, w,
                        [&](int lanes, PointAddMode mode) { return layers.get(lanes, mode); });
                    const double seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - start).count();
                    append_result(options.output_csv, result, seconds);
                    std::cout << "RESULT n=" << config.n << " w=" << w
                              << " DW=" << uint128_to_string(result.dw_full)
                              << " TDW=" << uint128_to_string(result.tdw)
                              << " seconds=" << std::fixed << std::setprecision(3) << seconds
                              << std::endl;
                }
                // Avoid signed overflow for a very large g; include the endpoint only when aligned.
                if (options.g > options.w_end - w) break;
                w += options.g;
            }
            return 0;
        }

        throw std::invalid_argument("unknown command: " + options.command);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << std::endl;
        return 2;
    }
}
