#include <memory.h>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <thread>
#include <iomanip>
#include "strategy.h"

using namespace std;

// random number generator
thread_local std::mt19937_64 rand_generator(std::random_device{}());

int get_ones(const vector<ROW> &m)
{
    int s = 0;
    for(int i = 0; i < m.size(); i++)
        s += m[i].count();
    return s;
}

void get_trans_matrix(vector<ROW> &trans_m, const vector<ROW> &m)
{
    int n = m.size();
    if (trans_m.size() != n) {
        trans_m.resize(n);
    }
    for(int i=0; i<n; ++i) trans_m[i].reset();

    const int BLOCK = 64; 
    for (int i = 0; i < n; i += BLOCK) {
        for (int j = 0; j < n; j += BLOCK) {
            int i_limit = std::min(i + BLOCK, n);
            int j_limit = std::min(j + BLOCK, n);
            for (int r = i; r < i_limit; ++r) {
                for (int c = j; c < j_limit; ++c) {
                    if (m[c].test(n - 1 - r)) {
                        trans_m[r].set(n - 1 - c);
                    }
                }
            }
        }
    }
}

int select_oper(const vector<ROW> &m, vector<xpair> &max_seq, int no_reduced, int opr_type)
{
    int rows = m.size();
    if (rows == 0) return no_reduced;

    static thread_local vector<int> row_weights;
    if (row_weights.size() < rows) row_weights.resize(rows);
    
    for (int i = 0; i < rows; ++i) {
        row_weights[i] = m[i].count();
    }

    for (int i = 0; i < rows; ++i) {
        int w_src = row_weights[i];
        if (w_src == 0) continue;

        const uint64_t* __restrict__ src_ptr = m[i].data;

        for (int j = 0; j < rows; ++j) {
            if (i == j) continue;

            
            const uint64_t* __restrict__ dst_ptr = m[j].data;
            int intersection_pop = 0;

            #pragma GCC unroll 18
            for (int k = 0; k < WORDS; ++k) {
                intersection_pop += __builtin_popcountll(src_ptr[k] & dst_ptr[k]);
            }

            int diff = (2 * intersection_pop) - w_src;

            if (diff > 0) {
                if (diff > no_reduced) {
                    no_reduced = diff;
                    max_seq.clear();
                    xpair new_ele;
                    new_ele.src = i;
                    new_ele.dst = j;
                    new_ele.flag = (opr_type == 1);
                    max_seq.push_back(new_ele);
                } else if (diff == no_reduced) {
                    xpair new_ele;
                    new_ele.src = i;
                    new_ele.dst = j;
                    new_ele.flag = (opr_type == 1);
                    max_seq.push_back(new_ele);
                }
            }
        }
    }
    return no_reduced;
}

vector<xpair> strgy1(vector<ROW> &m)
{
    vector<xpair> seq;
    int row_size = m.size();
    if (row_size == 0) return seq;
    int col_size = row_size; 

    int *mark = new int[row_size];
    memset(mark, 0, row_size * sizeof(int));
    xpair p;
    
    for (unsigned col = 1; col <= col_size; ++col)
    {
        unsigned r = 0;
        while (((r < row_size) && (!m[r].test(col_size - col))) || (mark[r] == 1))
        {
            ++r;
        }
        if (r >= row_size) continue;
        else mark[r] = 1;

        for (unsigned i = 0; i < row_size; ++i)
        {
            if (m[i].test(col_size - col) && (i != r))
            {
                m[i] ^= m[r];
                p.dst = i;
                p.src = r;
                p.flag = false;
                seq.push_back(p);
            }
        }
    }
    delete[] mark;
    return seq;
}
vector<xpair> strgy2(vector<ROW> &m)
{
    vector<ROW> trans_m;
    get_trans_matrix(trans_m, m);

    vector<xpair> seq = strgy1(trans_m);
    for(int i = 0; i < seq.size(); i++) seq[i].flag = true;
    
    get_trans_matrix(m, trans_m); 
    int tab[SIZE] = {0};
    build_table(m, tab);
    vector<xpair> final_seq = update_seq_str(seq, tab);
    return final_seq;
}

vector<xpair> strgy3(vector<ROW> &m)
{
    vector<xpair> tmp_seq;
    const int MAX_SAFETY_LOOPS = 50000; 
    int loops = 0;

    while(get_ones(m) != m.size() && loops < MAX_SAFETY_LOOPS)
    {
        loops++;
        vector<xpair> base_oper;
        base_oper.clear();
        int max_reduction = 0;

        max_reduction = select_oper(m, base_oper, max_reduction, 0);

        vector<ROW> trans_m;
        get_trans_matrix(trans_m, m);
        
        max_reduction = select_oper(trans_m, base_oper, max_reduction, 1);

        if (base_oper.size() >= 1)
        {
            int rand_num = rand_generator() % base_oper.size();
            xpair op = base_oper[rand_num];

            if (!op.flag) 
            {
                m[op.dst] ^= m[op.src];
                tmp_seq.push_back(op);
            }
            else 
            {
                trans_m[op.dst] ^= trans_m[op.src];
                get_trans_matrix(m, trans_m);
                tmp_seq.push_back(op);
            }
        }
        else
        {
            // No further optimization possible
            break;
        }
    }

    if (get_ones(m) != m.size())
    {
        int rnd = rand_generator() % 2;
        if(rnd == 0)
        {
            vector<xpair> seq_2(strgy1(m));
            tmp_seq.insert(tmp_seq.end(), seq_2.begin(), seq_2.end());
        }
        else
        {
            vector<xpair> seq_2(strgy2(m));
            tmp_seq.insert(tmp_seq.end(), seq_2.begin(), seq_2.end());
        }
    }

    int tab[SIZE] = {0};
    build_table(m, tab);
    vector<xpair> final_seq = update_seq_str(tmp_seq, tab);
    return final_seq;
}

void build_table(const vector<ROW> &m, int tab[SIZE])
{
    int ind = 0;
    for (int i = 0; i < m.size(); i++)
    {
        for (int j = 0; j < m.size(); j++)
        {
            if (m[i].test(j))
            {
                ind = m.size() - 1 - j;
                break;
            }
        }
        tab[ind] = i;
    }
}

vector<xpair> update_seq_str(const vector<xpair> &seq, int tab[SIZE])
{
    vector<xpair> tmp_seq;
    for(int i = 0; i < seq.size(); i++)
    {
        if(!seq[i].flag) tmp_seq.push_back(seq[i]);
    }
    for(int i = seq.size() - 1; i >= 0; i--)
    {
        if(seq[i].flag)
        {
            xpair new_ele;
            new_ele.src = tab[seq[i].dst];
            new_ele.dst = tab[seq[i].src];
            new_ele.flag = false;
            tmp_seq.push_back(new_ele);
        }
    }
    return tmp_seq;
}