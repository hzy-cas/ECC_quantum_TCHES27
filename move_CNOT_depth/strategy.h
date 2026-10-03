#ifndef STRATEGY_H
#define STRATEGY_H

#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <cstdint>
#include <random> 
#include <immintrin.h> 

using namespace std;

// Matrix configuration
#define SIZE 64*11
#define WORDS 11

// ==========================================
// High Performance Row Struct
// ==========================================
struct FastRow {
    uint64_t data[WORDS];

    FastRow() {
        memset(data, 0, sizeof(data));
    }

    FastRow(const std::string& s) {
        memset(data, 0, sizeof(data));
        int len = s.length();
        for (int i = 0; i < len; ++i) {
            if (s[len - 1 - i] == '1') {
                set(i);
            }
        }
    }

    inline void set(int index) {
        data[index / 64] |= (1ULL << (index % 64));
    }

    inline void reset() {
        memset(data, 0, sizeof(data));
    }

    inline bool test(int index) const {
        return (data[index / 64] >> (index % 64)) & 1ULL;
    }

    inline int count() const {
        int c = 0;
        for (int i = 0; i < WORDS; ++i) {
            c += __builtin_popcountll(data[i]);
        }
        return c;
    }

    inline void operator^=(const FastRow& other) {
        for (int i = 0; i < WORDS; ++i) {
            data[i] ^= other.data[i];
        }
    }
    
    inline int count_diff(const FastRow& other) const {
        int c = 0;
        const uint64_t* p1 = this->data;
        const uint64_t* p2 = other.data;
        
        #pragma GCC unroll 11
        for (int i = 0; i < WORDS; ++i) {
            c += __builtin_popcountll(p1[i] ^ p2[i]);
        }
        return c;
    }
    
    inline int count_diff_ptr(const uint64_t* p2) const {
        int c = 0;
        const uint64_t* p1 = this->data;
        
        #pragma GCC unroll 11
        for (int i = 0; i < WORDS; ++i) {
            c += __builtin_popcountll(p1[i] ^ p2[i]);
        }
        return c;
    }

    inline FastRow operator^(const FastRow& other) const {
        FastRow res;
        for (int i = 0; i < WORDS; ++i) {
            res.data[i] = this->data[i] ^ other.data[i];
        }
        return res;
    }

    bool operator==(const FastRow& other) const {
        for(int i=0; i<WORDS; ++i) if(data[i]!=other.data[i]) return false;
        return true;
    }
};

inline std::ostream& operator<<(std::ostream& os, const FastRow& r) {
    for (int i = SIZE - 1; i >= 0; i--) {
        os << (r.test(i) ? '1' : '0');
    }
    return os;
}

typedef FastRow ROW;

typedef struct{
    int src;
    int dst;
    bool flag; 
} xpair;

// Function declarations
vector<xpair> strgy1(vector<ROW> &m);
vector<xpair> strgy2(vector<ROW> &m);
vector<xpair> strgy3(vector<ROW> &m);

int get_ones(const vector<ROW> &m);
void get_trans_matrix(vector<ROW> &trans_m, const vector<ROW> &m);

int select_oper(const vector<ROW> &m, vector<xpair> &max_seq, int no_reduced, int opr_type, int sample_limit, std::mt19937 &gen);
// Overloaded version for compatibility if needed, or just update the main one
int select_oper(const vector<ROW> &m, vector<xpair> &max_seq, int no_reduced, int opr_type);

void build_table(const vector<ROW> &m, int tab[SIZE]);
vector<xpair> update_seq_str(const vector<xpair> &seq, int tab[SIZE]);

#endif