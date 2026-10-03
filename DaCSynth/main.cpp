/**

* @file main.cpp

* @brief Depth Optimization Solver for Linear Reversible Quantum Circuits (DaCSynth Algorithm Implementation)

* * This program implements the DaCSynth algorithm from the paper "Reducing the Depth of Linear Reversible Quantum Circuits," specifically for the greedy strategy and graph decomposition method described in Section III-B. However, it superimposes zhu23's depth reduction algorithm at the end.

*/

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <bitset>
#include <cassert>
#include <map>
#include <queue>
#include <chrono>

// ================== Global Statistics Structure ==================
struct SolverStats {
    long long greedy_phases = 0;       // Total iterations in greedy phase
    long long greedy_ops = 0;          // Total row/column operations in greedy phase
    long long graph_decomp_tasks = 0;  // Number of tasks triggering graph decomposition
    long long graph_decomp_ops = 0;    // Total flip operations in graph decomposition phase
    long long sum_max_degree = 0;      // Sum of max degrees for all graph decomposition tasks (theoretical minimum total depth)

    void reset() {
        greedy_phases = 0; greedy_ops = 0;
        graph_decomp_tasks = 0; graph_decomp_ops = 0;
        sum_max_degree = 0;
    }

    void printReport() {
        std::cout << "\n\n===== Solver Statistics Report (Section III-B Full Version) =====" << std::endl;
        std::cout << "Greedy Phase (Row/Column Cost Minimization):" << std::endl;
        std::cout << "  Total Iterations: " << greedy_phases << std::endl;
        std::cout << "  Number of Gates Generated: " << greedy_ops << std::endl;
        std::cout << "Graph Decomposition Phase (Residual Elimination / Theorem 3.2):" << std::endl;
        std::cout << "  Number of Subproblems Solved: " << graph_decomp_tasks << std::endl;
        std::cout << "  Number of Flip Operations: " << graph_decomp_ops << std::endl;
        std::cout << "  Theoretical Minimum Total Depth (Based on Max Degree): " << sum_max_degree << std::endl;
        std::cout << "====================================================\n" << std::endl;
    }
};

SolverStats GLOBAL_STATS;

// ================== Basic type definitions ==================
using Row = std::vector<uint8_t>;
using Matrix = std::vector<Row>;
using Word = uint64_t; // For bit-packed matrix operations (64 bits)

struct CNOT { int control; int target; };

struct BlockCoord { int r, c; }; 

// ================== Matrix Operation Utility Functions ==================
// These functions efficiently handle Boolean matrix operations (GF(2))

std::vector<std::vector<Word>> packMatrix(const Matrix &M, int &cols) {
    int r = M.size();
    if (r == 0) return {};
    cols = M[0].size();
    int words = (cols + 63) / 64;
    std::vector<std::vector<Word>> packed(r, std::vector<Word>(words, 0));
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (M[i][j]) packed[i][j / 64] |= (1ULL << (j % 64));
        }
    }
    return packed;
}

Matrix unpackMatrix(const std::vector<std::vector<Word>> &packed, int rows, int cols) {
    Matrix M(rows, Row(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if ((packed[i][j / 64] >> (j % 64)) & 1ULL) M[i][j] = 1;
            else M[i][j] = 0;
        }
    }
    return M;
}

int getRank(const Matrix &m) {
    if (m.empty() || m[0].empty()) return 0;
    int cols;
    auto mat = packMatrix(m, cols); 
    int r = mat.size();
    int words = mat[0].size();
    int rank = 0;

    for (int j = 0; j < cols && rank < r; ++j) {
        int pivot = rank;
        int word_idx = j / 64;
        Word bit_mask = (1ULL << (j % 64));

        while (pivot < r && !(mat[pivot][word_idx] & bit_mask)) pivot++;

        if (pivot < r) {
            std::swap(mat[rank], mat[pivot]);
            for (int i = 0; i < r; ++i) {
                if (i != rank && (mat[i][word_idx] & bit_mask)) {
                    for (int w = 0; w < words; ++w) mat[i][w] ^= mat[rank][w];
                }
            }
            rank++;
        }
    }
    return rank;
}

Matrix matMul(const Matrix &A, const Matrix &B) {
    if (A.empty() || B.empty()) return {};
    int n = A.size();
    int m = B.size(); 
    int p_cols = B[0].size();
    
    int b_cols_dummy;
    auto packedB = packMatrix(B, b_cols_dummy);
    int words = packedB[0].size();
    std::vector<std::vector<Word>> packedC(n, std::vector<Word>(words, 0));

    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < m; ++k) {
            if (A[i][k]) {
                for (int w = 0; w < words; ++w) packedC[i][w] ^= packedB[k][w];
            }
        }
    }
    return unpackMatrix(packedC, n, p_cols);
}

Matrix matInverse(const Matrix &A) {
    int n = A.size();
    if (n == 0) return {};
    int cols;
    auto mat = packMatrix(A, cols);
    int words = mat[0].size();
    
    std::vector<std::vector<Word>> inv(n, std::vector<Word>(words, 0));
    for (int i = 0; i < n; ++i) inv[i][i / 64] |= (1ULL << (i % 64));

    for (int i = 0; i < n; ++i) {
        int word_idx = i / 64;
        Word bit_mask = (1ULL << (i % 64));
        int pivot = i;
        
        while (pivot < n && !(mat[pivot][word_idx] & bit_mask)) pivot++;
        if (pivot == n) return {};

        if (pivot != i) {
            std::swap(mat[i], mat[pivot]);
            std::swap(inv[i], inv[pivot]);
        }
        
        for (int j = 0; j < n; ++j) {
            if (i != j && (mat[j][word_idx] & bit_mask)) {
                for (int w = 0; w < words; ++w) {
                    mat[j][w] ^= mat[i][w]; 
                    inv[j][w] ^= inv[i][w]; 
                }
            }
        }
    }
    return unpackMatrix(inv, n, n);
}

// ================== Graph Decomposer (RegularizedDecomposition) ==================
class RegularizedDecomposition {
    struct Edge { int u, v, original_idx; bool active; };
    int n_left, n_right, max_degree;
    std::vector<Edge> edges;
    std::vector<std::vector<int>> adj;

public:
    std::vector<std::vector<BlockCoord>> solve(const std::vector<BlockCoord> &flip_tasks) {
        if (flip_tasks.empty()) return {};
        
        std::map<int, int> r_deg, c_deg;
        int max_r = 0, max_c = 0;
        for (const auto &b : flip_tasks) {
            r_deg[b.r]++; c_deg[b.c]++;
            max_r = std::max(max_r, b.r); max_c = std::max(max_c, b.c);
        }
        max_degree = 0;
        for (auto const &[k, v] : r_deg) max_degree = std::max(max_degree, v);
        for (auto const &[k, v] : c_deg) max_degree = std::max(max_degree, v);
        
        GLOBAL_STATS.graph_decomp_tasks++;
        GLOBAL_STATS.sum_max_degree += max_degree;

        int num_nodes = std::max(max_r, max_c) + 1;
        n_left = num_nodes; n_right = num_nodes;
        edges.clear();
        for (int i = 0; i < (int)flip_tasks.size(); ++i) 
            edges.push_back({flip_tasks[i].r, flip_tasks[i].c, i, true});
        
        regularizeGraph();
        
        std::vector<std::vector<BlockCoord>> layers;
        for (int d = 0; d < max_degree; ++d) {
            std::vector<int> matching = findPerfectMatching();
            std::vector<BlockCoord> current_layer;
            for (int edge_idx : matching) {
                if (edges[edge_idx].original_idx != -1) 
                    current_layer.push_back(flip_tasks[edges[edge_idx].original_idx]);
                edges[edge_idx].active = false;
            }
            if (!current_layer.empty()) layers.push_back(current_layer);
        }
        return layers;
    }

private:
    void regularizeGraph() {
        std::vector<int> l_deg(n_left, 0), r_deg(n_right, 0);
        for (const auto &e : edges) { l_deg[e.u]++; r_deg[e.v]++; }
        int l_ptr = 0, r_ptr = 0;
        while (l_ptr < n_left && r_ptr < n_right) {
            if (l_deg[l_ptr] == max_degree) { l_ptr++; continue; }
            if (r_deg[r_ptr] == max_degree) { r_ptr++; continue; }
            edges.push_back({l_ptr, r_ptr, -1, true}); 
            l_deg[l_ptr]++; r_deg[r_ptr]++;
        }
    }

    std::vector<int> findPerfectMatching() {
        adj.assign(n_left, {});
        for (int i = 0; i < (int)edges.size(); ++i) 
            if (edges[i].active) adj[edges[i].u].push_back(i);
            
        std::vector<int> pair_u(n_left, -1), pair_v(n_right, -1), dist(n_left + 1);
        
        while (bfs(pair_u, pair_v, dist)) {
            for (int u = 0; u < n_left; ++u) 
                if (pair_u[u] == -1) dfs(u, pair_u, pair_v, dist);
        }
        
        std::vector<int> res;
        for (int u = 0; u < n_left; ++u) {
            int v = pair_u[u];
            if (v != -1) {
                for (int edge_idx : adj[u]) 
                    if (edges[edge_idx].v == v && edges[edge_idx].active) { 
                        res.push_back(edge_idx); break; 
                    }
            }
        }
        return res;
    }

    bool bfs(const std::vector<int> &pair_u, const std::vector<int> &pair_v, std::vector<int> &dist) {
        std::queue<int> q;
        for (int u = 0; u < n_left; ++u) {
            if (pair_u[u] == -1) { dist[u] = 0; q.push(u); } else dist[u] = 1e9;
        }
        dist[n_left] = 1e9;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            if (dist[u] < dist[n_left]) {
                for (int edge_idx : adj[u]) {
                    int v = edges[edge_idx].v; 
                    int next_u = (pair_v[v] == -1) ? n_left : pair_v[v];
                    if (dist[next_u] == 1e9) { 
                        dist[next_u] = dist[u] + 1; q.push(next_u); 
                    }
                }
            }
        }
        return dist[n_left] != 1e9;
    }

    bool dfs(int u, std::vector<int> &pair_u, std::vector<int> &pair_v, std::vector<int> &dist) {
        if (u != n_left) {
            for (int edge_idx : adj[u]) {
                int v = edges[edge_idx].v; 
                int next_u = (pair_v[v] == -1) ? n_left : pair_v[v];
                if (dist[next_u] == dist[u] + 1) {
                    if (dfs(next_u, pair_u, pair_v, dist)) { 
                        pair_v[v] = u; pair_u[u] = v; return true; 
                    }
                }
            }
            dist[u] = 1e9; return false;
        }
        return true;
    }
};

// ================== Main Solver (DaCSynth Algorithm) ==================
using Circuit = std::vector<std::vector<CNOT>>;

class DaCSolver {
public:
    std::vector<std::vector<CNOT>> finalCircuit;
    std::vector<int> p_map;
    std::vector<int> qubitFreeTime;
    int N_CURRENT;

    DaCSolver() { N_CURRENT = 0; }

    void solve(Matrix inputM, Matrix &resultM, Circuit &outputCircuit) {
        finalCircuit.clear();
        N_CURRENT = inputM.size();

        if (getRank(inputM) != N_CURRENT) {
            throw std::runtime_error("Input matrix is singular (rank < N).");
        }

        qubitFreeTime.assign(N_CURRENT, 0);
        p_map.resize(N_CURRENT); 
        std::iota(p_map.begin(), p_map.end(), 0);
        resultM = inputM;
        try {
            recursiveDecompose(resultM, N_CURRENT, 0); 
            postProcessOptimization(outputCircuit);
        }
        catch (const std::exception &e) { std::cerr << "Solver Error: " << e.what() << std::endl; }
    }

private:

    std::vector<std::vector<CNOT>> optimizeSeq(const std::vector<CNOT>& seq, int num_qubits) {
        if (seq.empty()) return {};

        std::vector<std::vector<CNOT>> layers;
        std::vector<std::vector<bool>> layer_occupied; 

        std::vector<int> last_as_src(num_qubits, -1);
        std::vector<int> last_as_dst(num_qubits, -1);

        int estimated_layers = seq.size() / 2 + 10;
        layers.reserve(estimated_layers);
        layer_occupied.reserve(estimated_layers);

        for (const auto& gate : seq) {
            int c = gate.control;
            int t = gate.target;

            int start_layer = 0;
            if (last_as_dst[c] > start_layer) start_layer = last_as_dst[c];
            if (last_as_src[t] > start_layer) start_layer = last_as_src[t];
            if (last_as_dst[t] > start_layer) start_layer = last_as_dst[t];

            bool placed = false;

            for (int j = start_layer; j < (int)layers.size(); ++j) {
                if (layer_occupied[j][c] || layer_occupied[j][t]) {
                    continue;
                }

                layers[j].push_back(gate);
                layer_occupied[j][c] = true;
                layer_occupied[j][t] = true;

                if (j > last_as_src[c]) last_as_src[c] = j;
                if (j > last_as_dst[t]) last_as_dst[t] = j;

                placed = true;
                break;
            }

            if (!placed) {
                layers.push_back({gate});
                std::vector<bool> new_row(num_qubits, false);
                new_row[c] = true;
                new_row[t] = true;
                layer_occupied.push_back(new_row);

                int new_layer_idx = layers.size() - 1;

                if (new_layer_idx > last_as_src[c]) last_as_src[c] = new_layer_idx;
                if (new_layer_idx > last_as_dst[t]) last_as_dst[t] = new_layer_idx;
            }
        }
        return layers;
    }

    void postProcessOptimization(Circuit &outputCircuit) {
        std::vector<CNOT> flat_seq;
        for (const auto& layer : finalCircuit) {
            for (const auto& gate : layer) flat_seq.push_back(gate);
        }

        auto layers_1 = optimizeSeq(flat_seq, N_CURRENT);

        std::vector<CNOT> flat_seq_rev;
        for (const auto& layer : layers_1) {
            for (const auto& gate : layer) flat_seq_rev.push_back(gate);
        }
        std::reverse(flat_seq_rev.begin(), flat_seq_rev.end());

        auto layers_2 = optimizeSeq(flat_seq_rev, N_CURRENT);

        std::vector<CNOT> flat_seq_final;
        for (const auto& layer : layers_2) {
            for (const auto& gate : layer) flat_seq_final.push_back(gate);
        }
        std::reverse(flat_seq_final.begin(), flat_seq_final.end());

        outputCircuit = optimizeSeq(flat_seq_final, N_CURRENT);
        finalCircuit = outputCircuit;
    }

    void scheduleGate(int c, int t) {
        int avail = std::max(qubitFreeTime[c], qubitFreeTime[t]);
        if (avail >= (int)finalCircuit.size()) finalCircuit.resize(avail + 1);
        finalCircuit[avail].push_back({c, t});
        qubitFreeTime[c] = avail + 1; qubitFreeTime[t] = avail + 1;
    }

    void ensureInvertibleLogical(Matrix &M, int size, int offset) {
        int half = size / 2;
        std::vector<int> candidates; for(int i=0; i<size; ++i) candidates.push_back(offset+i);
        std::vector<int> selected;
        
        for(int k=0; k<half; ++k) {
            bool found = false;
            for(int i=0; i<(int)candidates.size(); ++i) {
                int cand = candidates[i];
                std::vector<int> rows = selected; rows.push_back(cand);
                
                Matrix sub(rows.size(), Row(half));
                for(size_t r=0; r<rows.size(); ++r) for(int c=0; c<half; ++c) sub[r][c] = M[p_map[rows[r]]][offset+c];
                
                if(getRank(sub) == (int)rows.size()) {
                    selected.push_back(cand); candidates.erase(candidates.begin()+i); found = true; break;
                }
            }
            if(!found) throw std::runtime_error("Singular block");
        }
        
        std::vector<int> new_map = p_map;
        for(int i=0; i<half; ++i) new_map[offset+i] = p_map[selected[i]];
        for(int i=0; i<(int)candidates.size(); ++i) new_map[offset+half+i] = p_map[candidates[i]];
        p_map = new_map;
    }

    void zeroBlockGreedy(Matrix &globalM, int r_t, int r_p, int c_s, int h, int w, int depth_level) {
        if(h<=0 || w<=0) return;

        Matrix Ap(w, Row(w)), At(h, Row(w));
        for(int i=0; i<w; ++i) for(int j=0; j<w; ++j) Ap[i][j] = globalM[p_map[r_p+i]][c_s+j];
        for(int i=0; i<h; ++i) for(int j=0; j<w; ++j) At[i][j] = globalM[p_map[r_t+i]][c_s+j];

        Matrix Ap_inv = matInverse(Ap);
        if (Ap_inv.empty()) {
            throw std::runtime_error("Singular Ap block in zeroBlockGreedy()");
        }

        Matrix B = matMul(At, Ap_inv); 

        struct Op { int src, dst, gain; };
        int iter = 0;
        int no_progress_count = 0; 
        int last_ones = -1;

        while (true) {
            iter++;
            int current_ones = 0;
            for(const auto& r : B) for(int v : r) current_ones += v;

            if (iter % 50 == 0 || iter == 1) {
                std::cout << "\r    [Greedy] Level " << depth_level << " | Size " << h << "x" << w 
                          << " | Iter " << std::setw(6) << iter 
                          << " | Residue Ones: " << std::setw(4) << current_ones << "    " << std::flush;
            }

            if (current_ones == last_ones) {
                no_progress_count++;
                if (no_progress_count > 5) break; 
            } else {
                no_progress_count = 0;
                last_ones = current_ones;
            }

            if (current_ones == 0) break;

            std::vector<Op> potential_row_ops;
            for (int i = 0; i < h; ++i) {
                bool empty_i = true; for(int k=0; k<w; ++k) if(B[i][k]) { empty_i = false; break; }
                if (empty_i) continue;

                for (int j = 0; j < h; ++j) {
                    if (i == j) continue;
                    int w_j = 0; for(int k=0; k<w; ++k) w_j += B[j][k];
                    int w_new = 0; for(int k=0; k<w; ++k) w_new += (B[j][k] ^ B[i][k]);
                    if (w_j - w_new > 0) potential_row_ops.push_back({i, j, w_j - w_new});
                }
            }
            std::sort(potential_row_ops.begin(), potential_row_ops.end(), [](const Op& a, const Op& b){ return a.gain > b.gain; });

            std::vector<Op> potential_col_ops;
            for (int i = 0; i < w; ++i) {
                bool empty_i = true; for(int k=0; k<h; ++k) if(B[k][i]) { empty_i = false; break; }
                if (empty_i) continue;

                for (int j = 0; j < w; ++j) {
                    if (i == j) continue;
                    int w_j_col = 0; for(int k=0; k<h; ++k) w_j_col += B[k][j];
                    int w_new_col = 0; for(int k=0; k<h; ++k) w_new_col += (B[k][j] ^ B[k][i]);
                    if (w_j_col - w_new_col > 0) potential_col_ops.push_back({i, j, w_j_col - w_new_col});
                }
            }
            std::sort(potential_col_ops.begin(), potential_col_ops.end(), [](const Op& a, const Op& b){ return a.gain > b.gain; });

            long long total_row_gain = 0;
            std::vector<Op> selected_row_ops;
            std::vector<bool> r_busy(h, false);
            for(auto &op : potential_row_ops) {
                if (!r_busy[op.src] && !r_busy[op.dst]) {
                    selected_row_ops.push_back(op); 
                    r_busy[op.src] = true; r_busy[op.dst] = true;
                    total_row_gain += op.gain;
                }
            }

            long long total_col_gain = 0;
            std::vector<Op> selected_col_ops;
            std::vector<bool> c_busy(w, false);
            for(auto &op : potential_col_ops) {
                if (!c_busy[op.src] && !c_busy[op.dst]) {
                    selected_col_ops.push_back(op); 
                    c_busy[op.src] = true; c_busy[op.dst] = true;
                    total_col_gain += op.gain;
                }
            }

            if (selected_row_ops.empty() && selected_col_ops.empty()) break; 

            GLOBAL_STATS.greedy_phases++;

            if (total_row_gain >= total_col_gain) {
                for(auto &op : selected_row_ops) {
                    for(int k=0; k<w; ++k) B[op.dst][k] ^= B[op.src][k];
                    int pc = p_map[r_t + op.src]; int pt = p_map[r_t + op.dst];
                    scheduleGate(pc, pt);
                    for(int k=0; k<N_CURRENT; ++k) globalM[pt][k] ^= globalM[pc][k];
                    GLOBAL_STATS.greedy_ops++;
                }
            } else {
                for(auto &op : selected_col_ops) {
                    for(int k=0; k<h; ++k) B[k][op.dst] ^= B[k][op.src];
                    int pc = p_map[r_p + op.dst]; int pt = p_map[r_p + op.src];
                    scheduleGate(pc, pt);
                    for(int k=0; k<N_CURRENT; ++k) globalM[pt][k] ^= globalM[pc][k];
                    GLOBAL_STATS.greedy_ops++;
                }
            }
        }

        std::vector<BlockCoord> flip_tasks;
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                if (B[i][j]) flip_tasks.push_back({i, j});
            }
        }

        if (!flip_tasks.empty()) {
            RegularizedDecomposition graph_solver;
            auto layers = graph_solver.solve(flip_tasks);

            for (const auto& layer : layers) {
                for (const auto& coord : layer) {
                    int i = coord.r; int j = coord.c;
                    int pc = p_map[r_p + j]; 
                    int pt = p_map[r_t + i]; 
                    
                    scheduleGate(pc, pt);
                    for(int k=0; k<N_CURRENT; ++k) globalM[pt][k] ^= globalM[pc][k];
                    B[i][j] = 0;
                    GLOBAL_STATS.graph_decomp_ops++;
                }
            }
        }
        std::cout << "\r" << std::string(80, ' ') << "\r"; 
    }

    void recursiveDecompose(Matrix &M, int size, int offset, int depth_level = 0) {
        if(size <= 1) return;
        int half = size/2, rest = size-half;
        ensureInvertibleLogical(M, size, offset);
        zeroBlockGreedy(M, offset+half, offset, offset, rest, half, depth_level);
        zeroBlockGreedy(M, offset, offset+half, offset+half, half, rest, depth_level);
        recursiveDecompose(M, half, offset, depth_level + 1);
        recursiveDecompose(M, rest, offset+half, depth_level + 1);
    }
};
// ================== IO and File Handling ==================
Matrix readMatrixFromFile(const std::string &fn) {
    std::ifstream f(fn); 
    std::vector<std::string> lines; std::string l;
    while(std::getline(f, l)) {
        std::string cl=""; for(char c:l) if(c=='0'||c=='1') cl+=c;
        if(!cl.empty()) lines.push_back(cl);
    }
    int n = lines.size(); if (n == 0) return {};
    Matrix m; m.reserve(n);
    for(const auto& s : lines) {
        Row r; r.reserve(n);
        for(int i = 0; i < (int)s.size() && i < n; ++i) r.push_back(s[i] - '0');
        if ((int)r.size() < n) r.resize(n, 0); 
        m.push_back(r);
    }
    return m;
}

std::string insert_seq_before_ext(const std::string &fn) {
    size_t p = fn.rfind(".txt"); return (p!=std::string::npos ? fn.substr(0, p) : fn) + "_seq.txt";
}

void print_seq_for_python_reverse(const Matrix &m_orig, const Circuit &circuit, 
                                  const std::vector<int>& solver_p_map, const std::string &filename) {
    int n = m_orig.size(); 
    std::vector<int> real_to_file(n, -1);
    for(int i=0; i<n; ++i) if(i < (int)solver_p_map.size()) real_to_file[solver_p_map[i]] = i;
    int fill = n; for(int i=0; i<n; ++i) if(real_to_file[i] == -1) real_to_file[i] = fill++; 

    std::vector<CNOT> seq;
    for (const auto &layer : circuit) for (const auto &gate : layer) seq.push_back(gate);
    std::reverse(seq.begin(), seq.end()); 

    std::vector<int> wire_free_time(n, 0); 
    int true_depth = 0;
    for(const auto& g : seq) {
        if (g.control >= n || g.target >= n) continue; 
        int file_c = real_to_file[g.control]; int file_t = real_to_file[g.target];
        int gate_time = std::max(wire_free_time[file_c], wire_free_time[file_t]) + 1;
        wire_free_time[file_c] = gate_time; wire_free_time[file_t] = gate_time;
        if(gate_time > true_depth) true_depth = gate_time;
    }

    std::ofstream f(filename);
    for (const auto &row : m_orig) { for (int v : row) f << (v ? '1' : '0'); f << std::endl; }
    f << "CNOT Depth = " << true_depth << std::endl; 
    f << "CNOT Count = " << seq.size() << std::endl; 

    std::vector<bool> y_printed(n, false);
    std::vector<int> last_op_index(n + 5000, -1);
    for(size_t i=0; i<seq.size(); ++i) {
        if (seq[i].target < n) {
            int file_t = real_to_file[seq[i].target];
            if(file_t < (int)last_op_index.size()) last_op_index[file_t] = i;
        }
    }

    for(size_t i=0; i<seq.size(); ++i) {
        int real_c = seq[i].control; int real_t = seq[i].target;
        if (real_c >= n || real_t >= n) continue; 
        int file_c = real_to_file[real_c]; int file_t = real_to_file[real_t];
        f << "x[" << file_t << "] = x[" << file_t << "] ^ x[" << file_c << "]";
        if (i == (size_t)last_op_index[file_t]) {
            int k = real_t; if (k < n) { f << "   y[" << k << "]"; y_printed[k] = true; }
        }
        f << std::endl;
    }
    for(int k=0; k<n; ++k) {
        if(!y_printed[k]) {
            int file_w = real_to_file[k];
            if (n > 1) {
                f << "x[" << file_w << "] = x[" << file_w << "] ^ x[" << (file_w==0?1:0) << "]\n";
                f << "x[" << file_w << "] = x[" << file_w << "] ^ x[" << (file_w==0?1:0) << "]   y[" << k << "]\n";
            } else f << "x[" << file_w << "]   y[" << k << "]\n";
        }
    }
    std::cout << "[Output] Written to " << filename << " (N=" << n << ", Depth: " << true_depth << ")" << std::endl;
}

// ================== Main Function ==================
int main() {
    std::string inputDir = "input", outputDir = "result";
    std::filesystem::create_directories(outputDir);
    GLOBAL_STATS.reset();

    for (const auto &entry : std::filesystem::directory_iterator(inputDir)) {
        if (!entry.is_regular_file()) continue;
        std::string fn = entry.path().filename().string();
        std::cout << "\n[Processing] " << fn << std::endl;

        try {
            Matrix M = readMatrixFromFile(entry.path().string());
            int n_dynamic = M.size();
            if (n_dynamic == 0) continue;
            std::cout << "  > Detected matrix size N = " << n_dynamic << std::endl;

            Matrix M_dummy; Circuit circ; DaCSolver solver;
            
            auto start_time = std::chrono::high_resolution_clock::now();
            
            solver.solve(M, M_dummy, circ); 
            
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end_time - start_time;
            
            std::cout << "  > Time elapsed: " << diff.count() << " seconds" << std::endl;
            std::cout << "  > Number of gates: " << circ.size() << " (uncompressed layers)" << std::endl;
            
            print_seq_for_python_reverse(M, circ, solver.p_map, insert_seq_before_ext(outputDir + "/result_" + fn));
        }
        catch (const std::exception &e) { std::cerr << "  [Error] " << e.what() << std::endl; }
    }
    GLOBAL_STATS.printReport();
    return 0;
}