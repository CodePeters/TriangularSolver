/**
 * @author:   George Petrou georgepetrou08@gmail.com
 * @description: This file contains the functions for topological oredring
 * and level ordering of the DAG graph. It also contains functions for
 * checking corectness of results.
 */

#include "ordering.h"

/**
 * Implements topological sort of the DAG graph using DFS variation algorithm.
 *
 * @param graph: adjacency list that repsresnts the DAG graph
 * @param v: the current node
 * @param visited: 
 * @return_param st: the vector that represents a stack with the results.
 */
void topological_sort(vector<int>* graph, int v, bool* visited, vector<int>& st) { 
    // Mark the current node as visited. 
    visited[v] = true; 

    // Recur for all neighbour vertices 
    for (auto i = graph[v].begin(); i != graph[v].end(); ++i) {
        if (!visited[*i]) {
            topological_sort(graph, *i, visited, st); 
        }
    }
  
    // Push current vertex to stack which stores result 
    st.push_back(v); 
} 


/**
 * Finds the Reach set of non zero elements of b in a topological
 * ordering using DFS. Function reachSet_dfs only starts from non visited
 * non zero B elements and calls topological_sort for dfs traversal.
 *
 * @param n: the number of nodes in graph
 * @param graph: adjacency list that repsresnts the DAG graph
 * @param nz_b_indexes: non zero elements of b vector
 * @return_param st: the vector that represents a stack with the results.
 */
void reachSet_dfs(int n, vector<int>* graph, vector<int>& st, vector<int> nz_b_indexes) {
  // arrays keeping visited nodes
	bool* visited = new bool[n];
	fill(visited, visited + n, false);

  //recur all non zero elements of b and call topological_sort for dfs
	for (auto i = nz_b_indexes.begin(); i != nz_b_indexes.end(); ++i) {
      if (visited[*i] == false) { 
            topological_sort(graph, *i, visited, st);
      }
  }

  delete[] visited;
}


/**
 * Modified Khan's algorithm for topological ordering. Since Khan's
 * algorithm finds a topological ordering but also orders
 * the DAG by levels, we use it for level ordering of prunned graph (i.e the
 * graph that contains the nodes discovered by dfs).
 *
 * @param n: the number of nodes in graph
 * @param graph: adjacency list that repsresnts the DAG graph
 * @param st: the vector-stack with the results from topological order.
 * @return_param lev_order: vector that contains vectors with nodes at each level
 */
void level_order_Khan(int n, vector<int>* graph, vector<int> st, vector<vector<int>>& lev_order) { 
  vector<int> in_degree(n, 0);
  queue<int> q; 

  for (auto i = st.begin(); i != st.end(); ++i) {
      for (auto j = graph[*i].begin(); j != graph[*i].end(); ++j) {
        in_degree[*j]++;
      }
  }

  for (auto i = st.begin(); i != st.end(); ++i) {
      if (in_degree[*i] == 0) 
            q.push(*i); 
  }

  while (!q.empty()) { 
        // Extract front of queue 
        // (or perform dequeue) 
        // and add it to topological order 
        int s = q.size();
        vector<int> local_level;
        while (s) {

          int u = q.front(); 
          q.pop(); 
          s--;
          local_level.push_back(u); 
  
          // Iterate through all its 
          // neighbouring nodes 
          // of dequeued node u and 
          // decrease their in-degree 
          // by 1 
          for (auto itr = graph[u].begin(); itr != graph[u].end(); itr++) {
              // If in-degree becomes zero, 
              // add it to queue 
              if (--in_degree[*itr] == 0) 
                  q.push(*itr); 
          }
        }
        lev_order.push_back(local_level);
  } 
}


/**
 * Checks if two doubles are essentially equal, see here:
 * https://stackoverflow.com/a/253874/6671342
 * Epsilon used is defined in header file
 *
 * param a: first double
 * param b: second double to compare
 */
bool are_equal(double a, double b) {
    return std::abs(a - b) <= ( (std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * EPSILON);
}
