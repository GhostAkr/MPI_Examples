#include "../include/Matr.h"

using namespace std;

// Creation & destroying

double** createMatr(int _rows, int _cols) {
	double** result = new double*[_rows];
	for (int i = 0; i < _rows; ++i) {
		result[i] = new double[_cols];
		for (int j = 0; j < _cols; ++j) {
			result[i][j] = 0.0;
		}
	}
	return result;
}

void deleteMatr(double** _matr, int _rows) {
	for (int i = 0; i < _rows; ++i) {
		delete[] _matr[i];
	}
	delete[] _matr;
}

// Operations

double** matrixMult(double** _source1, double** _source2, int m, int n, int s) {
	double** result = new double*[m];
	for (int i = 0; i < m; ++i) {
		result[i] = new double[s];
		for (int j = 0; j < s; ++j) {
			result[i][j] = 0;
		}
	}
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < s; j++) {
				result[i][j] += _source1[i][k] * _source2[k][j];
			}
		}
	}
	return result;
}

// Staff things

void printMatr(double** _mesh, int _rows, int _cols) {
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			cout << _mesh[i][j] << " ";
		}
		cout << endl;
	}
}
