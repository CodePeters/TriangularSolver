/**
 * @author:   George Petrou georgepetrou08@gmail.com
 * @description: This file contains the declarations of functions, varaibles, constants
 * and library includes.
 */

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <vector>
#include <iostream>
#include <fstream>
#include "sprsmatrix.h"
#include "sprsmatrix.cpp"
#include <string>
#include <queue>
#include <math.h> 
#include <omp.h>
using namespace std;

template<typename T> SprsMatrix<T> read_MM_format(string filename, vector<int>*& graph);
template<typename T> void read_MM_format_b(string filename, T*& b, vector<int>& nz_b_indexes);
template<typename T> void serial_lsolve (int n, SprsMatrix<T> &matrix, T* x);
template<typename T> void serial_lsolve2 (int n, SprsMatrix<T> &L , T* x, vector<int> st);
template<typename T> void parallel_lsolve (int n, SprsMatrix<T> &L , T* x, vector<int> st);
template<typename T> void parallel_lsolve2 (int n, SprsMatrix<T> &L , T* x, vector<vector<int>>& lev);

#endif