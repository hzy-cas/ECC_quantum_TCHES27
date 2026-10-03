#include "../BasicLogic/QuantumLib.h"
#include <fstream>
#include <stdexcept>

namespace {
std::vector<int> slice(const std::vector<int>& v, size_t start, size_t length) {
    if (start > v.size() || length > v.size()-start) throw std::out_of_range("balanced register slice");
    return {v.begin()+start,v.begin()+start+length};
}
std::vector<int> join(std::vector<int> a,const std::vector<int>& b) {
    a.insert(a.end(),b.begin(),b.end()); return a;
}
std::vector<int> allocate(uint64_t& cursor,uint64_t size) {
    if (cursor+size >= (1ULL<<30)) throw std::overflow_error("packed wire index limit");
    std::vector<int> result;
    result.reserve(size);
    for (uint64_t i=0;i<size;++i) result.push_back(static_cast<int>(cursor++));
    return result;
}
std::vector<std::vector<int>> split(const std::vector<int>& flat,int n,int k) {
    std::vector<std::vector<int>> result;
    for (int i=0;i<k;++i) result.push_back(slice(flat,static_cast<size_t>(i)*n,n));
    return result;
}
std::vector<uint8_t> hex_bits(std::string value,int n) {
    if (value.rfind("0x",0)==0) value.erase(0,2);
    std::vector<uint8_t> out(n);
    size_t pos=0;
    for (auto i=value.rbegin();i!=value.rend();++i) {
        int digit=*i>='0' && *i<='9'?*i-'0':
                  *i>='a' && *i<='f'?*i-'a'+10:
                  *i>='A' && *i<='F'?*i-'A'+10:-1;
        if (digit<0) throw std::invalid_argument("invalid fixture hex");
        for (int b=0;b<4;++b,++pos) {
            if (pos<out.size()) out[pos]=(digit>>b)&1;
            else if ((digit>>b)&1) throw std::invalid_argument("fixture value too wide");
        }
    }
    return out;
}
}

// The original Balanced Montgomery and coordinate compute-copy-uncompute
// construction. Tail inputs are the SAME clean z,v produced by the retained
// architecture's QROM pair. Reduction uses y1 ^= y2 as in reduction_point_add,
// saving an extra v register. No QROM gate or window partition changes here.
RetainedStats QuantumContext::run_balanced_window_layer(int lanes,bool early_tail,
    const std::vector<bool>& curve_a,const std::string& fixture) {
    if (lanes<1 || curve_a.size()!=static_cast<size_t>(config.n))
        throw std::invalid_argument("invalid balanced layer arguments");
    gm.clear();
    const uint64_t n=config.n,k=lanes,block=config.block_size,m=config.mul_cost_size;
    const uint64_t target_count=m*(3*(k-1)+config.inversion_mul_blocks);
    uint64_t cursor=0;
    auto rx=split(allocate(cursor,n*k),n,k),ry=split(allocate(cursor,n*k),n,k);
    std::vector<std::vector<int>> tx,ty;
    if (!early_tail) {
        tx=split(allocate(cursor,n*k),n,k); ty=split(allocate(cursor,n*k),n,k);
    }
    auto z=split(allocate(cursor,n*k),n,k);
    auto v=early_tail?split(allocate(cursor,n*k),n,k):ry;
    const auto z_input=z;
    auto ancilla=allocate(cursor,4*block*k);
    Toffoli_qubits=allocate(cursor,target_count);
    auto square=allocate(cursor,4*n);
    auto recovery=allocate(cursor,n*(k/2));
    auto lambda=allocate(cursor,n*k),out_y_flat=allocate(cursor,n*k),inverse=allocate(cursor,n*k);
    auto out_y=split(out_y_flat,n,k);
    const uint64_t width=cursor;

    if (!early_tail) for (int i=0;i<lanes;++i) {
        CNOT_N(rx[i],z[i]); CNOT_N(tx[i],z[i]); CNOT_N(ty[i],v[i]);
    }
    auto mont=montgomery_clean(z,v,ancilla,recovery,square,lambda,inverse);
    for (int i=0;i<lanes;++i) {
        CONST_ADD(z[i],curve_a);
        auto transformed=Squaring_plus_input(join(mont.lambdas[i],z[i]),2*n);
        mont.lambdas[i]=std::move(transformed.first); z[i]=std::move(transformed.second);
        const auto& px=early_tail?rx[i]:tx[i];
        const auto& py=early_tail?ry[i]:ty[i];
        CNOT_N(px,z[i]);
        const auto begin=gm.current_pointer();
        auto dirty_y=mul(join(z[i],slice(ancilla,2*block*i,block)),
            join(mont.lambdas[i],slice(ancilla,2*block*i+block,block)),
            slice(Toffoli_qubits,i*m,m),0);
        const auto end=gm.current_pointer();
        CNOT_N(dirty_y,out_y[i]);
        gm.replay_reverse(begin,end);
        CNOT_N(px,z[i]); CNOT_N(z[i],out_y[i]); CNOT_N(py,out_y[i]);
    }
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
                throw std::runtime_error("balanced affine point-add output mismatch");
            state[width+2*n*i+j]=state[z[i][j]];
            state[width+2*n*i+n+j]=state[out_y[i][j]];
        }
        // These and ONLY these blocks are declared reusable between layers.
        for (const auto* pool:{&ancilla,&Toffoli_qubits,&square,&recovery})
            for (auto q:*pool) if (state[q]) throw std::runtime_error("balanced reusable workspace is not clean");
        gm.apply_to_basis(state,true);
        for (size_t q=0;q<width;++q) if (state[q]!=initial[q])
            throw std::runtime_error("balanced compute-copy-uncompute failed");
        for (int i=0;i<lanes;++i) for (size_t j=0;j<n;++j)
            if (state[width+2*n*i+j]!=expected_x[i][j] || state[width+2*n*i+n+j]!=expected_y[i][j])
                throw std::runtime_error("balanced output copy failed");
        verified=true;
    }
    LayerStats resources{width,std::get<3>(counts),std::get<4>(counts),
        std::get<0>(counts),std::get<1>(counts),std::get<2>(counts)};
    const uint64_t inputs=4*n*k,clean=4*block*k+target_count+4*n+n*(k/2);
    gm.clear();
    return {resources,inputs,clean,width-inputs-clean,2*n*k,early_tail?n*k:0,verified};
}
