#include "src/BasicLogic/QuantumLib.h"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
#include <sys/resource.h>

namespace {
int integer(const std::string& value) {
    size_t used=0; int result=std::stoi(value,&used);
    if (used!=value.size()) throw std::invalid_argument("invalid integer");
    return result;
}
std::vector<bool> bits(std::string value,int n) {
    if (value.rfind("0x",0)==0) value.erase(0,2);
    std::vector<bool> result(n); size_t pos=0;
    for (auto i=value.rbegin();i!=value.rend();++i) {
        const std::string digits="0123456789abcdef";
        auto digit=digits.find(static_cast<char>(std::tolower(*i)));
        if (digit==std::string::npos) throw std::invalid_argument("invalid curve coefficient hex");
        for (int b=0;b<4;++b,++pos) {
            bool value=(digit>>b)&1;
            if (pos<result.size()) result[pos]=value;
            else if (value) throw std::invalid_argument("curve coefficient too wide");
        }
    }
    return result;
}
}

int main(int argc,char** argv) {
    try {
        if (argc==1 || std::string(argv[1])=="--help") {
            std::cout<<"retained_estimator --n N --lanes K --kind reduction|tail --data-root DIR "
                "--curve-a-hex L_BASIS_HEX [--architecture retained|balanced] [--verify-fixture FILE]\n"; return 0;
        }
        std::map<std::string,std::string> options;
        for (int i=1;i<argc;i+=2) {
            if (i+1==argc) throw std::invalid_argument("option missing value");
            std::string key=argv[i];
            if (key!="--n" && key!="--lanes" && key!="--kind" && key!="--data-root" &&
                key!="--curve-a-hex" && key!="--verify-fixture" && key!="--memory-limit-bytes" &&
                key!="--architecture") throw std::invalid_argument("unknown option "+key);
            if (!options.emplace(key,argv[i+1]).second) throw std::invalid_argument("duplicate option");
        }
        if (options.count("--memory-limit-bytes")) {
#ifdef __linux__
            auto bytes=std::stoull(options.at("--memory-limit-bytes"));
            rlimit limit{static_cast<rlim_t>(bytes),static_cast<rlim_t>(bytes)};
            if (setrlimit(RLIMIT_AS,&limit)!=0) throw std::runtime_error("cannot set job address-space limit");
#else
            throw std::invalid_argument("memory limit requires Linux");
#endif
        }
        int n=integer(options.at("--n")),k=integer(options.at("--lanes"));
        std::string kind=options.at("--kind");
        const std::string architecture=options.count("--architecture")?options.at("--architecture"):"retained";
        if (architecture!="retained" && architecture!="balanced") throw std::invalid_argument("invalid architecture");
        if (kind!="reduction" && kind!="tail") throw std::invalid_argument("invalid kind");
        if (k<1 || k>2*(n+1)) throw std::invalid_argument("invalid lanes");
        auto config=GFConfig::setup(n,options.at("--data-root"));
        if (!std::filesystem::is_directory(config.data_root)) throw std::invalid_argument("missing bundled data");
        auto start=std::chrono::steady_clock::now();
        QuantumContext context(config);
        const auto coefficient=bits(options.at("--curve-a-hex"),n);
        const auto fixture=options.count("--verify-fixture")?options.at("--verify-fixture"):"";
        const auto result=architecture=="retained"
            ?context.run_retained_layer(k,kind=="tail",coefficient,fixture)
            :context.run_balanced_window_layer(k,kind=="tail",coefficient,fixture);
        const auto& s=result.resources;
        rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
        uint64_t rss=usage.ru_maxrss;
#else
        uint64_t rss=static_cast<uint64_t>(usage.ru_maxrss)*1024;
#endif
        std::cout<<"{\"model_version\":\"window-"<<architecture<<"-clean-zv-v1\",\"n\":"<<n
            <<",\"lanes\":"<<k<<",\"kind\":\""<<kind<<"\",\"qubits\":"<<s.qubits
            <<",\"toffoli\":"<<s.toffoli<<",\"cnot\":"<<s.cnot
            <<",\"full_depth\":"<<s.full_depth<<",\"current_depth\":"<<s.current_depth
            <<",\"toffoli_depth\":"<<s.toffoli_depth<<",\"input_qubits\":"<<result.input_qubits
            <<",\"clean_qubits\":"<<result.clean_qubits<<",\"retained_new_qubits\":"<<result.retained_new_qubits
            <<",\"output_qubits\":"<<result.output_qubits
            <<",\"output_alias_input_qubits\":"<<result.output_alias_input_qubits
            <<",\"basis_verified\":"<<(result.basis_verified?"true":"false")
            <<",\"peak_rss_bytes\":"<<rss<<",\"seconds\":"
            <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
        return 0;
    } catch (const std::exception& e) { std::cerr<<"error: "<<e.what()<<'\n'; return 1; }
}
