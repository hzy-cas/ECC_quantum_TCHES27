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

constexpr const char* kModelVersion = "ac-balanced-clean-v1";

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
    std::string tables;
    std::string curve_a_hex;
    fs::path data_root = "data";
    fs::path output_csv;
    fs::path layer_cache;
    bool force = false;
};

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
        else if (arg == "--tables") options.tables = value;
        else if (arg == "--curve-a-hex") options.curve_a_hex = value;
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

std::vector<bool> parse_hex_bits(std::string value, int bits, const char* option) {
    if (value.rfind("0x", 0) == 0 || value.rfind("0X", 0) == 0) value.erase(0, 2);
    if (value.empty()) throw std::invalid_argument(std::string(option) + " is empty");
    std::vector<bool> result(static_cast<std::size_t>(bits), false);
    int output_bit = 0;
    for (auto it = value.rbegin(); it != value.rend(); ++it) {
        unsigned nibble = 0;
        if (*it >= '0' && *it <= '9') nibble = static_cast<unsigned>(*it - '0');
        else if (*it >= 'a' && *it <= 'f') nibble = 10U + static_cast<unsigned>(*it - 'a');
        else if (*it >= 'A' && *it <= 'F') nibble = 10U + static_cast<unsigned>(*it - 'A');
        else throw std::invalid_argument(std::string("invalid hex digit in ") + option);
        for (int j = 0; j < 4; ++j, ++output_bit) {
            const bool bit = ((nibble >> j) & 1U) != 0;
            if (output_bit < bits) result[static_cast<std::size_t>(output_bit)] = bit;
            else if (bit) throw std::invalid_argument(std::string(option) + " exceeds field width");
        }
    }
    return result;
}

PackedQromTable load_qrom_table(const fs::path& path, int address_bits, int n) {
    if (address_bits <= 0 || address_bits >= 31) {
        throw std::invalid_argument("QROM address bits outside 1..30");
    }
    const std::size_t coordinate_bytes = static_cast<std::size_t>((n + 7) / 8);
    const std::size_t entry_bytes = 2 * coordinate_bytes;
    const std::size_t entry_count = std::size_t{1} << address_bits;
    const std::size_t expected = entry_count * entry_bytes;
    if (!fs::is_regular_file(path) || fs::file_size(path) != expected) {
        throw std::invalid_argument("QROM file has wrong size: " + path.string());
    }
    std::ifstream input(path, std::ios::binary);
    std::vector<unsigned char> bytes(expected);
    input.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (!input) throw std::runtime_error("could not read QROM file: " + path.string());

    PackedQromTable table;
    table.address_bits = address_bits;
    table.word_bits = 2 * n;
    const std::size_t chunks = (static_cast<std::size_t>(table.word_bits) + 63U) / 64U;
    table.words.assign(entry_count * chunks, 0);
    for (std::size_t entry = 0; entry < entry_count; ++entry) {
        const std::size_t base = entry * entry_bytes;
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
            const std::size_t coordinate_base = base + static_cast<std::size_t>(coordinate) * coordinate_bytes;
            for (int bit = 0; bit < n; ++bit) {
                const std::size_t byte_index = coordinate_base + coordinate_bytes - 1U -
                                               static_cast<std::size_t>(bit / 8);
                if (((bytes[byte_index] >> (bit % 8)) & 1U) == 0) continue;
                const int packed_bit = coordinate * n + bit;
                table.words[entry * chunks + static_cast<std::size_t>(packed_bit) / 64U] |=
                    uint64_t{1} << (packed_bit % 64);
            }
        }
    }
    return table;
}

std::vector<PackedQromTable> parse_qrom_tables(const std::string& specs, int n) {
    if (specs.empty()) throw std::invalid_argument("window-layer requires --tables bits:path,...");
    std::vector<PackedQromTable> result;
    for (const auto& spec : split_csv(specs)) {
        const auto separator = spec.find(':');
        if (separator == std::string::npos) {
            throw std::invalid_argument("each --tables item must be bits:path");
        }
        const int bits = parse_int(spec.substr(0, separator), "table address bits");
        result.push_back(load_qrom_table(spec.substr(separator + 1), bits, n));
    }
    return result;
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
            if (f.size() < 12 || f[0] != kModelVersion || parse_int(f[1], "cache n") != config_.n) {
                continue;
            }
            if (fs::path(f[11]).lexically_normal() !=
                fs::path(config_.data_root).lexically_normal()) {
                continue;
            }
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
        << "  window-layer --n N --tables bits:path,... --curve-a-hex HEX --data-root PATH\n"
        << "  window-init  --n N --tables bits:path,... --data-root PATH\n"
        << "  qrom-pair    --n N --tables bits:path --data-root PATH\n"
        << "  early-tail   --n N --lanes K --curve-a-hex HEX --data-root PATH\n"
        << "  primitive --n N --kind multiplication|inversion --data-root PATH\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        const auto options = parse_options(argc, argv);
        if (options.command == "help" || options.command == "--help" || options.command == "-h") {
            print_help();
            return 0;
        }

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

        if (options.command == "window-layer") {
            if (options.curve_a_hex.empty()) {
                throw std::invalid_argument("window-layer requires --curve-a-hex");
            }
            QuantumContext context(config);
            const auto tables = parse_qrom_tables(options.tables, config.n);
            const auto curve_a = parse_hex_bits(options.curve_a_hex, config.n, "--curve-a-hex");
            print_stats_json(context.run_window_layer(tables, curve_a));
            return 0;
        }

        if (options.command == "window-init") {
            QuantumContext context(config);
            const auto tables = parse_qrom_tables(options.tables, config.n);
            print_stats_json(context.run_window_initialization(tables));
            return 0;
        }

        if (options.command == "qrom-pair") {
            QuantumContext context(config);
            const auto tables = parse_qrom_tables(options.tables, config.n);
            if (tables.size() != 1) {
                throw std::invalid_argument("qrom-pair requires exactly one table");
            }
            print_stats_json(context.run_qrom_pair(tables.front()));
            return 0;
        }

        if (options.command == "early-tail") {
            if (options.curve_a_hex.empty()) {
                throw std::invalid_argument("early-tail requires --curve-a-hex");
            }
            QuantumContext context(config);
            const auto curve_a = parse_hex_bits(options.curve_a_hex, config.n, "--curve-a-hex");
            print_stats_json(context.run_early_tail(options.lanes, curve_a));
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
