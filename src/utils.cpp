/**
 * @author:   George Petrou georgepetrou08@gmail.com
 * @description: This file contains the functions for topological oredring
 * and level ordering of the DAG graph. It also contains functions for
 * checking corectness of results.
 */

#include "utils.h"

/**
 * Reads the sparse matrix in Market Matrix format and stores it in 
 * CSC format. While reading the input it also build the column
 * dependency DAG.
 *
 * @param graph: adjacency list that repsresnts the DAG graph
 * @param filename: the name of the file containing the sparse matrix in CSC format
 * @return_param graph: the dependency graph between columns 
 * @returns: the sparse matrix object
 */
template<typename T> 
SprsMatrix<T> read_MM_format(string filename, vector<int>*& graph)
{
	ifstream file(filename);
	int num_row, num_col, num_lines;

	// Ignore comments headers
	while (file.peek() == '%') file.ignore(2048, '\n');

	// Read number of rows and columns
	file >> num_row>> num_col >> num_lines;

	SprsMatrix<T> matrix(num_col);

	// initialize graph in form of adjacency list
	graph = new vector<int>[num_row];

	// to keep track when to change column in CSC array of pointers
	int prev_col = 0; 
	// we need to keep track of how element belong in lower trinagular
	// part to set column pointers correctly
	int lower_trinagular_cnt = 0;
	// the position in column pointers array that we set
	int pos = 0;

	// fill the matrix and graph with data
	for (int l = 0; l < num_lines; l++)
	{
	    T data;
	    int row, col;
	    file >> row >> col >> data;

	    //skip upper triangular data
	    if (row<col) {
	    	continue;
	    }

	    // add value to sparse matrix values vector
	    matrix.addVal(data);
	    // similarly add the corresonding row index
	    matrix.addIndex(row-1);

	    // we changed column
	    if (prev_col < col){
	    	pos = col-1;
	    	// set the pos position of column pointers array to
	    	// the value lower_trinagular_cnt i.e the elements
	    	// that we have read until now
	    	matrix.addPtr(pos, lower_trinagular_cnt);
	    }

	    //add edge in graph
	    if (col != row) {
	    	graph[col-1].push_back(row-1);
	    }
	    prev_col = col;

	    // update the number of element we have read
	    lower_trinagular_cnt++;
	}
	// (n+1)th element of column pointers
	matrix.addPtr(pos+1, lower_trinagular_cnt);
	file.close();
	return matrix;
}


/**
 * Reads and stores the b vector 
 *
 * @param filename: the name of the file containing the b elements 
 * @return_param b: the array that reprsents the b vector of Lx=b system
 * @return_param nz_b_indexes: the indexes of non zero elements of b
 */
template<typename T> 
void read_MM_format_b(string filename, T*& b, vector<int>& nz_b_indexes)
{
	ifstream file(filename);
	int num_row, num_col, num_lines;
	// Ignore comments headers
	while (file.peek() == '%') file.ignore(2048, '\n');

	// Read number of rows and columns
	file >> num_row>> num_col >> num_lines;
	b = new T[num_row];

	for (int l = 0; l < num_lines; l++)
	{
	    T data;
	    int row, col;
	    file >> row >> col >> data;
	    b[row-1] = data;
	    // add index of non zero elements of b
	    nz_b_indexes.push_back(row-1);
	}
	file.close();
}


/**
 * Solves Lower trinagular system: Lx=b. This is the
 * naive serial implementation
 *
 * @param n: the dimension of square sparse matrix 
 * @param L: the square sparse array L
 * @return_param x: the matrix initially containg b values. In the parameter
 * th results are stored.
 */
template<typename T> 
void serial_lsolve (int n, SprsMatrix<T> &L , T* x) {
	double start_time = omp_get_wtime();

	for (int j = 0 ; j < n ; ++j) {
		// skip zero X elements
		if (abs(x[j])<1e-9) continue;
		x[j] /= L.getVal(L.getPtr(j));
		for (int p = L.getPtr(j)+1 ; p < L.getPtr(j+1) ; ++p) {
			x[L.getIndex(p)] -= L.getVal(p) * x[j] ;
		}	
	}
	double elapsed_time = omp_get_wtime() - start_time;
	cout << "Serial time: " << elapsed_time << " μsec" << endl; 
}


/**
 * Solves Lower trinagular system: Lx=b. This is a better
 * serial implementation where we only iterate in non zero x 
 * elements in a topological sorted (in the DAG) order. 
 *
 * @param n: the dimension of square sparse matrix 
 * @param L: the square sparse array L
 * @return_param x: the matrix initially containg b values. In the parameter
 *   th results are stored.
 * param st: the vector-stack with the indexes of non zero x, returnd 
 *   by the topological sorting.
 */
template<typename T> 
void serial_lsolve2 (int n, SprsMatrix<T> &L, T* x, vector<int> st) {
	double start_time = omp_get_wtime();

    for (auto i = st.rbegin(); i != st.rend(); ++i ) { 
    	int j = *i;
		x[j] /= L.getVal(L.getPtr(j));
		for (int p = L.getPtr(j)+1 ; p < L.getPtr(j+1) ; ++p) {
			x[L.getIndex(p)] -= L.getVal(p)*x[j];
		}
	}
	double elapsed_time = omp_get_wtime() - start_time;
	cout << "Serial time: " << elapsed_time << " μsec" << endl; 
}


/**
 * Solves Lower trinagular system: Lx=b. This is a naive
 * parallel implementation where we only iterate in non zero x 
 * elements in a topological sorted (in the DAG) order. The parallel
 * part is in subtracting the Lij*x[i] value from corresponding rows
 * in parallel.
 *
 * @param n: the dimension of square sparse matrix 
 * @param L: the square sparse array L
 * @return_param x: the matrix initially containg b values. In the parameter
 *   th results are stored.
 * param st: the vector-stack with the indexes of non zero x, returnd 
 *   by the topological sorting.
 */
template<typename T> 
void parallel_lsolve (int n, SprsMatrix<T> &L , T* x, vector<int> st) {
	int p;
	omp_set_num_threads(2);
	double start_time = omp_get_wtime();

	for (auto i = st.rbegin(); i != st.rend(); ++i ) {
		int j = *i;
		x[j] /= L.getVal(L.getPtr(j));

	    #pragma omp parallel for private(p) schedule(static)
			for ( p = L.getPtr(j)+1 ; p < L.getPtr(j+1) ; ++p) {
				x[L.getIndex(p)] -= L.getVal(p)*x[j];
			}
	}
	double elapsed_time = omp_get_wtime() - start_time;
	cout << "Parallel time: " << elapsed_time << " μsec" << endl;
}


/**
 * Solves Lower trinagular system: Lx=b. This is a more advanced
 * parallel implementation where we iterate in non zero x 
 * elements of topological sort and also we parallelize the 
 * procedure for element that are in the same level and thus don't depend
 * on each other. 
 *
 * @param n: the dimension of square sparse matrix 
 * @param L: the square sparse array L
 * @return_param x: the matrix initially containg b values. In the parameter
 *   th results are stored.
 * param lev: the vector of vectors containg elements in same level, e.g lev[0]
 *   is a vector with elemnts at level 0.
 */
template<typename T> 
void parallel_lsolve2 (int n, SprsMatrix<T> &L , T* x, vector<vector<int>>& lev) {
	int p,j;
    omp_set_num_threads(2);
    double start_time = omp_get_wtime();

    for (auto i = lev.begin(); i != lev.end(); ++i) {
	    #pragma omp for  private(j,p) schedule(static)
	    for (auto w= (*i).begin(); w != (*i).end(); ++w) {  	
			j = *w;
			x[j] /= L.getVal(L.getPtr(j));

			for ( p = L.getPtr(j)+1 ; p < L.getPtr(j+1) ; ++p) {
				#pragma omp atomic
				x[L.getIndex(p)] -= L.getVal(p)*x[j];
			}
			
		}
		
	}
	double elapsed_time = omp_get_wtime() - start_time;
	cout << "Parallel time: " << elapsed_time << " μsec" << endl;
}
