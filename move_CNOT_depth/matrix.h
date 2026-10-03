#ifndef MATRIX_H
#define MATRIX_H

#include "strategy.h" 
#include <vector>
#include <string>
#include <iostream>

using namespace std;

typedef struct
{
    vector<xpair> seq;
    int gap;
    int start;
    int len;
} thread_data;

std::vector<ROW> get_matrix(const std::string& filename, int& out_width);
std::vector<ROW> get_matrix();

#endif