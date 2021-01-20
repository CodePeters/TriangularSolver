/**
 * @author:   George Petrou georgepetrou08@gmail.com
 * @description: This file contains the declarations of functions, varaibles, constants
 * and library includes.
 */

#include "utils.h"

#define EPSILON 1e-9

void topological_sort(vector<int>* graph, int v, bool* visited, vector<int>& st);
void reachSet_dfs(int n, vector<int>* graph, vector<int>& st, vector<int> nz_b_indexes);
void level_order_Khan(int n, vector<int>* graph, vector<int> st, vector<vector<int>>& lev_order);
bool are_equal(double a, double b);