/**
 * @author:   George Petrou georgepetrou08@gmail.com
 * @description: This file contains the definitions of all member functions
 * of sparse matrix class.
 */

#include "sprsmatrix.h"

// constructor for square matrices only
template<typename T>SprsMatrix<T>::SprsMatrix(int n)
{
	this->n = n;
	this->vals = vector<T>();
	this->indexes = vector<int>();
	this->pointers = vector<int>(n + 1);
}

template<typename T>int SprsMatrix<T>::getDim(void) const
{
	return this->n;
}

template<typename T>void SprsMatrix<T>::addVal(T val)
{
	this->vals.push_back(val);
}

template<typename T>T SprsMatrix<T>::getVal(int pos)
{
	return this->vals[pos];
}

template<typename T>void SprsMatrix<T>::addIndex(int index)
{
	this->indexes.push_back(index);
}

template<typename T>int SprsMatrix<T>::getIndex(int pos)
{
	return this->indexes[pos];
}

template<typename T>void SprsMatrix<T>::addPtr(int pos, int ptrValue)
{
	this->pointers[pos] = ptrValue;
}

template<typename T>int SprsMatrix<T>::getPtr(int pos)
{
	return this->pointers[pos];
}

template<typename T>void SprsMatrix<T>::printVals(void)
{
	for (auto i = this->vals.begin(); i != this->vals.end(); ++i)
	    cout << *i << ' ';
	cout << endl;
}

template<typename T>void SprsMatrix<T>::printIndexes(void)
{
	for (auto i = this->indexes.begin(); i != this->indexes.end(); ++i)
	    cout << *i << ' ';
	cout << endl;
}

template<typename T>void SprsMatrix<T>::printPtrs(void)
{
	for (auto i = this->pointers.begin(); i != this->pointers.end(); ++i)
	    cout << *i << ' ';
	cout << endl;
}

template<typename T>void SprsMatrix<T>::printMatrix(void)
{
	printVals();
	printIndexes();
	printPtrs();
}

