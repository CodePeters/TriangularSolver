/**
 * @author:   George Petrou georgepetrou08@gmail.com
 * @description: This is the mentry point where xecution starts.
 *
 */

#include "utils.h"
#include "utils.cpp"
#include "ordering.h"

int main(int argc, char **argv) {

  // check for input format
  if (argc != 3){
        cout << "Usage " << argv[0] <<" L_Matix B_Vector" << endl;
        exit(0);
  }

  //efficient io
  cout.sync_with_stdio(false); 

  vector<int> *graph;            //the DAG representing column dependancies of sparse Matrix
  vector<int>  nz_b_indexes;     // indexes of non zero elements of b
  vector<vector<int>> lev_order; // level order of final (prunned) DAG
  vector<int> st;                // stack for topological ordering of reach set of b vector  
  double *b;                     // the array with the final solution of serial implementation
  double *b2;                    // the array with the final solution of parallel implementation
  double *b3;                    // the array with the final solution of // implementation
  int n;                         // the dimension of Sparse matrix 

  // read b vector
  read_MM_format_b<double>(argv[2], b, nz_b_indexes);

  // read sparse matrix
  SprsMatrix<double> M = read_MM_format<double>(argv[1], graph);

  n = M.getDim();

  // find the reach set of non zero elements of b
  // and topologicaly orders this (prunned)sub DAG of initial dag 
  reachSet_dfs(n, graph, st, nz_b_indexes);

  // takes the topological elements in stack (st)
  // and produces a level order
  level_order_Khan(n, graph, st, lev_order);

  // two new b vectors for checking
  // corectnes of two parallel algorithms
  b2 = new double[n];
  copy(b, b+n, b2);
  b3 = new double[n];
  copy(b, b+n, b3);

  //serial_lsolve(n, M, b);
  serial_lsolve2(n, M, b, st);
  parallel_lsolve(n, M, b2, st);
  parallel_lsolve2(n, M, b3, lev_order);

  // checking for corectness of results
  for (int i = 0; i < n; ++i) {
  		if (are_equal(b[i], b2[i]) == false) {
        cout << "Different results for serail and parallel code !!" << endl;
        exit(0);
  		}
      if (are_equal(b[i], b3[i]) == false) {
        cout << "Different results for serail and 2nd parallel code !!" << endl;
        exit(0);
      }
  }

  //free heap allocated memory
  delete[] b;
  delete[] b2;
  delete[] b3;
  delete[] graph;

  return 0;
}
