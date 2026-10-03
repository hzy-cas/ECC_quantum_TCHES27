#ifndef REDUCE_H
#define REDUCE_H

#include "strategy.h"
#include "matrix.h"

#ifndef THREAD_NUM
#define THREAD_NUM 15 
#endif

vector<xpair> reduce(vector<ROW>& m);

bool reduce0(vector<xpair> &seq, const vector<vector<int>>& table); 

int reduce_step(vector<xpair> &seq);
int** get_table(vector<xpair> seq, int ** table, int osize, int nsize);

int get_ones(vector<ROW> m);

vector<ROW> get_reduced_matrix(const vector<xpair>& seq, const vector<ROW>& m);

int get_CNOT_Depth(vector<xpair>seq, int width);
int swap_reduce_CNOT_depth(vector<xpair>& seq, int width);
int move_reduce_CNOT_depth(const vector<xpair>& seq, vector<xpair>& seq_out);

#endif