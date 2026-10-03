#include "../BasicLogic/QuantumLib.h"
#include <fstream>
#include <limits>
#include <stdexcept>

namespace {
std::vector<int> slice(const std::vector<int>& v, size_t start, size_t length) {
    if (start > v.size() || length > v.size() - start) throw std::out_of_range("register slice");
    return {v.begin() + start, v.begin() + start + length};
}
std::vector<int> join(std::vector<int> a, const std::vector<int>& b) {
    a.insert(a.end(), b.begin(), b.end()); return a;
}
std::vector<int> allocate(uint64_t& cursor, uint64_t size) {
    if (cursor + size >= (1ULL << 30)) throw std::overflow_error("packed wire index limit");
    std::vector<int> result;
    result.reserve(size);
    for (uint64_t i = 0; i < size; ++i) result.push_back(static_cast<int>(cursor++));
    return result;
}
std::vector<std::vector<int>> split(const std::vector<int>& flat, int n, int k) {
    std::vector<std::vector<int>> result;
    for (int i = 0; i < k; ++i) result.push_back(slice(flat, static_cast<size_t>(i)*n, n));
    return result;
}
std::vector<uint8_t> hex_bits(std::string value, int n) {
    if (value.rfind("0x", 0) == 0) value.erase(0, 2);
    std::vector<uint8_t> out(n);
    int pos = 0;
    for (auto it = value.rbegin(); it != value.rend(); ++it) {
        int d = *it >= '0' && *it <= '9' ? *it-'0' :
                *it >= 'a' && *it <= 'f' ? *it-'a'+10 :
                *it >= 'A' && *it <= 'F' ? *it-'A'+10 : -1;
        if (d < 0) throw std::invalid_argument("invalid fixture hex");
        for (int b = 0; b < 4; ++b, ++pos) {
            if (pos < n) out[pos] = (d >> b) & 1;
            else if ((d >> b) & 1) throw std::invalid_argument("fixture value too wide");
        }
    }
    return out;
}
}

QuantumContext::RetainedMontgomery QuantumContext::montgomery_retained(
    const std::vector<std::vector<int>>& denominators,
    const std::vector<std::vector<int>>& numerators,
    const std::vector<int>& ancilla, const std::vector<int>& recovery_copy,
    const std::vector<int>& square_workspace) {
    const int k = static_cast<int>(denominators.size());
    const int n = config.n, block = config.block_size, m = config.mul_cost_size;
    int count = 0;
    auto main = slice(ancilla, 0, 2*block);
    auto side = slice(ancilla, ancilla.size()-2*block, 2*block);
    std::vector<std::vector<int>> inverses;
    if (k == 1) {
        inverses.push_back(inversion_logic(denominators[0], square_workspace, main, side, count));
    } else {
        auto product = cons_mul(denominators, count, ancilla);
        count = product.count;
        inverses.push_back(inversion_logic(product.product, square_workspace, main, side, count));
        for (int level = static_cast<int>(product.tree.size())-2; level >= 0; --level) {
            const auto& children = product.tree[level];
            std::vector<std::vector<int>> next;
            size_t copy_offset = 0;
            for (size_t i = 0; i < children.size(); i += 2) {
                const auto& parent = inverses[i/2];
                if (i+1 == children.size()) { next.push_back(parent); continue; }
                auto copy = slice(recovery_copy, copy_offset, n);
                copy_offset += n;
                // Only the linear copy is uncomputed; multiplication outputs remain live.
                CNOT_N(parent, copy);
                auto a = mul(join(parent, slice(ancilla, 2*block*i, block)),
                    join(children[i+1], slice(ancilla, 2*block*i+block, block)),
                    slice(Toffoli_qubits, count, m), 0);
                count += m;
                auto b = mul(join(copy, slice(ancilla, 2*block*(i+1), block)),
                    join(children[i], slice(ancilla, 2*block*(i+1)+block, block)),
                    slice(Toffoli_qubits, count, m), 0);
                count += m;
                CNOT_N(parent, copy);
                next.push_back(std::move(a)); next.push_back(std::move(b));
            }
            inverses = std::move(next);
        }
    }
    std::vector<std::vector<int>> lambdas;
    for (int i = 0; i < k; ++i) {
        lambdas.push_back(mul(join(numerators[i], slice(ancilla, 2*block*i, block)),
            join(inverses[i], slice(ancilla, 2*block*i+block, block)),
            slice(Toffoli_qubits, count, m), 0));
        count += m;
    }
    return {std::move(lambdas), count};
}

RetainedStats QuantumContext::run_retained_layer(int lanes, bool early_tail,
    const std::vector<bool>& curve_a, const std::string& fixture) {
    if (lanes < 1 || curve_a.size() != static_cast<size_t>(config.n))
        throw std::invalid_argument("invalid retained layer arguments");
    gm.clear();
    const uint64_t n=config.n, k=lanes, block=config.block_size, m=config.mul_cost_size;
    const uint64_t target_count=m*(3*(k-1)+config.inversion_mul_blocks+2*k);
    uint64_t cursor=0;
    auto rx=split(allocate(cursor,n*k),n,k), ry=split(allocate(cursor,n*k),n,k);
    std::vector<std::vector<int>> tx,ty;
    if (!early_tail) { tx=split(allocate(cursor,n*k),n,k); ty=split(allocate(cursor,n*k),n,k); }
    auto z=split(allocate(cursor,n*k),n,k), v=split(allocate(cursor,n*k),n,k);
    const auto z_input=z;
    auto ancilla=allocate(cursor,4*block*k);
    Toffoli_qubits=allocate(cursor,target_count);
    auto square=allocate(cursor,4*n);
    auto recovery=allocate(cursor,n*(k/2));
    const uint64_t width=cursor;

    if (!early_tail) for (int i=0;i<lanes;++i) {
        CNOT_N(rx[i],z[i]); CNOT_N(tx[i],z[i]);
        CNOT_N(ry[i],v[i]); CNOT_N(ty[i],v[i]);
    }
    auto mont=montgomery_retained(z,v,ancilla,recovery,square);
    std::vector<std::vector<int>> out_y;
    int count=mont.count;
    for (int i=0;i<lanes;++i) {
        CONST_ADD(z[i],curve_a);
        auto transformed=Squaring_plus_input(join(mont.lambdas[i],z[i]),2*n);
        mont.lambdas[i]=std::move(transformed.first); z[i]=std::move(transformed.second);
        CNOT_N(rx[i],z[i]);
        auto y=mul(join(z[i],slice(ancilla,2*block*i,block)),
            join(mont.lambdas[i],slice(ancilla,2*block*i+block,block)),
            slice(Toffoli_qubits,count,m),0);
        count+=m;
        CNOT_N(rx[i],z[i]); CNOT_N(z[i],y); CNOT_N(ry[i],y);
        out_y.push_back(std::move(y));
    }
    if (static_cast<uint64_t>(count)!=target_count) throw std::logic_error("target liveness mismatch");
    const auto counts=gm.optimize_and_count(width);
    bool verified=false;
    if (!fixture.empty()) {
        std::ifstream in(fixture);
        int fixture_n=0,fixture_k=0;
        if (!(in>>fixture_n>>fixture_k) || fixture_n!=config.n || fixture_k!=lanes)
            throw std::invalid_argument("fixture header mismatch");
        std::vector<uint8_t> state(width+2*n*k,0);
        std::vector<std::vector<uint8_t>> expected_x,expected_y;
        for (int i=0;i<lanes;++i) {
            std::string values[6];
            for (auto& value:values) if (!(in>>value)) throw std::invalid_argument("short fixture");
            auto a=hex_bits(values[0],n),b=hex_bits(values[1],n);
            auto c=hex_bits(values[2],n),d=hex_bits(values[3],n);
            expected_x.push_back(hex_bits(values[4],n)); expected_y.push_back(hex_bits(values[5],n));
            for (size_t j=0;j<n;++j) {
                state[rx[i][j]]=a[j]; state[ry[i][j]]=b[j];
                if (early_tail) { state[z_input[i][j]]=a[j]^c[j]; state[v[i][j]]=b[j]^d[j]; }
                else { state[tx[i][j]]=c[j]; state[ty[i][j]]=d[j]; }
            }
        }
        const auto initial=state;
        gm.apply_to_basis(state);
        for (int i=0;i<lanes;++i) for (size_t j=0;j<n;++j) {
            if (state[z[i][j]]!=expected_x[i][j] || state[out_y[i][j]]!=expected_y[i][j])
                throw std::runtime_error("affine point-add output mismatch");
            state[width+2*n*i+j]=state[z[i][j]];
            state[width+2*n*i+n+j]=state[out_y[i][j]];
        }
        for (auto q:ancilla) if (state[q]) throw std::runtime_error("arithmetic scratch is not clean");
        for (auto q:recovery) if (state[q]) throw std::runtime_error("recovery scratch is not clean");
        gm.apply_to_basis(state,true);
        for (size_t q=0;q<width;++q) if (state[q]!=initial[q])
            throw std::runtime_error("compute-copy-uncompute failed");
        for (int i=0;i<lanes;++i) for (size_t j=0;j<n;++j)
            if (state[width+2*n*i+j]!=expected_x[i][j] || state[width+2*n*i+n+j]!=expected_y[i][j])
                throw std::runtime_error("output copy failed");
        verified=true;
    }
    LayerStats resources{width,std::get<3>(counts),std::get<4>(counts),
        std::get<0>(counts),std::get<1>(counts),std::get<2>(counts)};
    const uint64_t inputs=4*n*k, clean=4*block*k+n*(k/2);
    gm.clear();
    return {resources,inputs,clean,width-inputs-clean,2*n*k,early_tail?n*k:0,verified};
}
