/**
 * @author:   George Petrou georgepetrou08@gmail.com
 * @description: This file contains the declarations of sparse matrix class
 */

#ifndef SPRSMATRIX_H_INCLUDED
#define SPRSMATRIX_H_INCLUDED

#include <vector>
#include <iostream>
using namespace std;

template<typename T> class SprsMatrix {

	public:
		// In the scope of this assignement
		// we only support square sparse matrices
		// but it would be easy to be modified
		SprsMatrix(int n);
		int getDim(void) const;

		// Get/Set functions to add and 
		// access elements of Sparse array
		void addVal(T val);
		T getVal(int pos);
		void addIndex(int index);
		int getIndex(int pos);
		void addPtr(int pos, int ptrValue);
		int getPtr(int pos);

		// functions to print contents
		// of sparse matrix mainly for debugging
		void printVals();
		void printIndexes();
		void printPtrs();
		void printMatrix();

	protected:		
		int n;
		vector<T> vals;
		// we don't name e.g column_pointers but we keep
		// naming general to support other formats in future
		vector<int> pointers, indexes;
};

#endif 
