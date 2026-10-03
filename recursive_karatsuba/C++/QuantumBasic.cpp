#include "QuantumLib.h"
#include "config.h" 
#include <cmath>
#include <algorithm>
#include <limits>

extern GFConfig current_config;

std::vector<int> QuantumContext::X(int a) {
    gm.add_X(a);
    return {a};
}

std::vector<int> QuantumContext::Toffoli_gate(int a, int b, int c) {
    gm.add_Toffoli(a, b, c);
    return {a, b, c};
}

std::vector<int> QuantumContext::CNOT(int a, int b) {
    gm.add_CNOT(a, b);
    return {a, b};
}

void QuantumContext::CNOT_n(const std::vector<int>& a, const std::vector<int>& b) {
    int n = a.size();
    for(int i=0; i<n; ++i) CNOT(a[i], b[i]);
}

void QuantumContext::CONST_ADD_n(const std::vector<int>& target, int val) {
    int n = target.size();

    if (val == -1) {
        for(int i=0; i<n; ++i) {
            X(target[i]);
        }
        return;
    }

    // Read only the bits that are representable by int. Shifting an int by
    // its width (or more) is undefined behavior, which previously repeated
    // the low bits on some platforms for extension degrees n > 32.
    const unsigned int unsigned_val = static_cast<unsigned int>(val);
    const int value_bits = std::numeric_limits<unsigned int>::digits;
    for(int i=0; i<n && i<value_bits; ++i) {
        if ((unsigned_val >> i) & 1U) X(target[i]);
    }
}
void QuantumContext::CONST_ADD_n(const std::vector<int>& target, const std::vector<int>& val_bits) {
    int n = target.size();
    for(int i=0; i<n; ++i) {

        if (i < val_bits.size() && val_bits[i] == 1) {
            X(target[i]);
        }
    }
}
void QuantumContext::Round_constant_XOR(const std::vector<int>& k, int rc, int bit) {
    for(int i=0; i<bit; ++i) {
        if ((rc >> i) & 1) X(k[i]);
    }
}

void QuantumContext::CSWAP(int c, int t1, int t2) {
    // Fredkin Gate
    CNOT(t2, t1);
    Toffoli_gate(c, t1, t2);
    CNOT(t2, t1);
}

std::vector<int> QuantumContext::copy_parallel(int value, std::vector<int> ancillas, int n) {
    int divide = (int)std::log2(n);
    int last = n - (int)std::pow(2, divide);

    std::vector<int> copy_list;
    copy_list.push_back(value);
    for(int i=0; i<n-1; ++i) copy_list.push_back(ancillas[i]);

    for(int i=0; i<divide; ++i) {
        int limit = (int)std::pow(2, i);
        for(int j=0; j<limit; ++j) {
            if(i == 0) {
                CNOT(copy_list[0], copy_list[1]);
            } else {
                CNOT(copy_list[j], copy_list[(int)std::pow(2, i) + j]);
            }
        }
    }

    for(int i=0; i<last; ++i) {
        CNOT(copy_list[i], copy_list[(int)std::pow(2, divide) + i]);
    }
    return copy_list;
}


std::vector<int> range_vec(int start, int n) {
    std::vector<int> v(n);
    for(int i=0; i<n; ++i) v[i] = start + i;
    return v;
}

std::vector<int> concat(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> v = a;
    v.insert(v.end(), b.begin(), b.end());
    return v;
}

std::vector<int> slice(const std::vector<int>& v, int start, int len) {
    if (start >= v.size()) return {};
    int actual_len = std::min(len, (int)v.size() - start);
    return std::vector<int>(v.begin() + start, v.begin() + start + actual_len);
}
