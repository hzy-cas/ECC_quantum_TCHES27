#include "reduce.h"
#include <pthread.h>
#include <cstdlib>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>

vector<ROW> get_reduced_matrix(const vector<xpair>& seq, const vector<ROW>& m)
{
    vector<ROW> tmp_m = m; 
    for(int i = 0; i < seq.size(); i++)
        tmp_m[seq[i].dst] ^= tmp_m[seq[i].src];
    return tmp_m;
}


// ============================================================
// move_reduce_CNOT_depth
// ============================================================
int move_reduce_CNOT_depth(const vector<xpair>& seq, vector<xpair>& seq_out) {
    vector<vector<xpair>> layers;

    vector<ROW> layer_all; 
    
    vector<int> last_as_src(SIZE, -1);
    vector<int> last_as_dst(SIZE, -1);

    int estimated_layers = seq.size() / 2 + 10;
    layers.reserve(estimated_layers);
    layer_all.reserve(estimated_layers);

    for (const auto& gate : seq) {
        int c = gate.src;
        int t = gate.dst;
        int start_layer = 0;
        if (last_as_dst[c] > start_layer) start_layer = last_as_dst[c];
        if (last_as_src[t] > start_layer) start_layer = last_as_src[t];

        bool merged = false;

        for (int j = start_layer; j < layers.size(); ++j) {
            
            if (layer_all[j].test(c) || layer_all[j].test(t)) {
                continue; 
            }

            layers[j].push_back(gate);
            layer_all[j].set(c);
            layer_all[j].set(t);

            if (j > last_as_src[c]) last_as_src[c] = j;
            if (j > last_as_dst[t]) last_as_dst[t] = j;
            

            merged = true;
            break; 
        }

        if (!merged) {
            layers.push_back({gate});
            
            ROW mask_a;
            mask_a.set(c); 
            mask_a.set(t);
            layer_all.push_back(mask_a);
            
            int new_layer_idx = layers.size() - 1;
            
            if (new_layer_idx > last_as_src[c]) last_as_src[c] = new_layer_idx;
            if (new_layer_idx > last_as_dst[t]) last_as_dst[t] = new_layer_idx;
        }
    }

    seq_out.clear();
    seq_out.reserve(seq.size());
    for (const auto& layer : layers) {
        seq_out.insert(seq_out.end(), layer.begin(), layer.end());
    }

    return layers.size();
}
vector<xpair> reduce(vector<ROW>& m)
{
    vector<xpair> seq = strgy3(m);
    return seq;
}