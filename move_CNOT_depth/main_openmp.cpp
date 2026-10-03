#include <fstream>
#include <sys/time.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <omp.h>
#include <vector>
#include <atomic>
#include <climits>
#include "strategy.h"
#include "matrix.h" 
#include "reduce.h"

using namespace std;

void print_seq(vector<ROW> m, vector<ROW> tmp_m, vector<xpair> seq, const std::string& output_filename, int final_depth, int original_width);

volatile int global_min_depth = INT_MAX;
volatile int global_counter = INT_MAX;
vector<xpair> best_seq;
vector<ROW> best_tmp_m;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <matrix_file_path>" << endl;
        return 1;
    }
    string input_filename = argv[1];
    string output_filename = "result_" + input_filename.substr(input_filename.find_last_of("/\\") + 1);

    int bound_cnot = INT_MAX;
    //Number of rounds
    const int TOTAL_LOOPS = 200;
    struct timeval start, end;
    
    gettimeofday(&start, NULL);

    cout << "Using OpenMP with up to " << omp_get_max_threads() << " threads." << endl;

    int original_width = 0;
    vector<ROW> original_m = get_matrix(input_filename, original_width);

    if (original_m.empty())
    {
        cout << "Matrix read failed" << endl;
        return -1;
    }

    cout << "Detected Matrix Size: " << original_m.size() << " rows x " << original_width << " cols." << endl;
    cout << "Internal SIZE: " << SIZE << endl;

    std::atomic<int> progress_cnt(0);

    #pragma omp parallel
    {
        vector<ROW> m; 
        m.reserve(original_m.size()); 
        vector<xpair> seq; 
        vector<xpair> seq_out;

        #pragma omp for schedule(dynamic, 1)
        for (int say = 0; say < TOTAL_LOOPS; say++)
        {
            m = original_m;
            seq = reduce(m);
            if (seq.size() > bound_cnot) { 
                progress_cnt++;
                continue; 
            }

            vector<ROW> tmp_m = get_reduced_matrix(seq, original_m);
            seq_out.clear();
            int depth = move_reduce_CNOT_depth(seq, seq_out);
            reverse(seq_out.begin(), seq_out.end());
            seq.clear();
            depth = move_reduce_CNOT_depth(seq_out, seq);
            reverse(seq.begin(), seq.end());


            bool potentially_better = false;
            if (depth < global_min_depth) potentially_better = true;
            else if (depth == global_min_depth && seq.size() < global_counter) potentially_better = true;

            if (potentially_better) {
                #pragma omp critical (update_best_solution)
                {
                    if ((depth < global_min_depth) || (depth == global_min_depth && seq.size() < global_counter)) 
                    {
                        global_min_depth = depth;
                        global_counter = seq.size();
                        best_seq = seq;
                        best_tmp_m = tmp_m;
                        
                        cout << "\n[Thread " << omp_get_thread_num() << "] >>> New Record Found! <<<" << endl;
                        cout << "  Depth: " << depth << " CNOT: " << seq.size() << endl;
                        print_seq(original_m, best_tmp_m, best_seq, output_filename, global_min_depth, original_width);
                    }
                }
            }
            
            // progress printing
            int current_val = ++progress_cnt;
            if (current_val % 5 == 0 || current_val == TOTAL_LOOPS) {
                #pragma omp critical (print_prog)
                {
                    cout << "Progress: " << current_val << " / " << TOTAL_LOOPS << " completed." << endl;
                }
            }
        }
    }

    gettimeofday(&end, NULL);
    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    cout << "Total Running Time: " << elapsed_time << " s" << endl;

    return 0;
}

void print_seq(vector<ROW> m, vector<ROW> tmp_m, vector<xpair> seq, const std::string& output_filename, int final_depth, int original_width)
{
    ofstream f;
    f.open(output_filename, std::ios::out | std::ios::trunc);
    if (!f.is_open()) return;

    f << "Original Matrix:" << endl;
    for (size_t i = 0; i < m.size(); i++) {
        for (int k = original_width - 1; k >= 0; k--) {
            f << (m[i].test(k) ? '1' : '0');
        }
        f << endl;
    }

    f << endl << endl;
    f << "Reduced Matrix:" << endl;
    for (size_t i = 0; i < tmp_m.size(); i++) {
        for (int k = original_width - 1; k >= 0; k--) {
            f << (tmp_m[i].test(k) ? '1' : '0');
        }
        f << endl;
    }

    f << endl << endl;
    f << "CNOT Count = " << seq.size() << endl;
    f << "CNOT Depth = " << final_depth << endl; 

    int tab[SIZE] = {0};
    for (size_t i = 0; i < tmp_m.size(); i++)
    {
        for (size_t j = 0; j < tmp_m.size(); j++)
        {
            if (tmp_m[i].test(tmp_m.size() - 1 - j))
            {
                tab[i] = j;
                break;
            }
        }
    }

    for (int i = seq.size() - 1; i >= 0; i--)
    {
        f << "x[" << tab[seq[i].dst] << "] = x[" << tab[seq[i].dst] << "] ^ x[" << tab[seq[i].src] << "]";
        tmp_m[seq[i].dst] ^= tmp_m[seq[i].src];
        bool flag = false;
        for (size_t j = 0; j < m.size(); j++)
        {
            if (tmp_m[seq[i].dst] == m[j])
            {
                flag = true;
                f << "    y[" << j << "]" << endl;
                break;
            }
        }

        if (!flag) f << endl;
    }
    f.close();
}