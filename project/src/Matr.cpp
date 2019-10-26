#include "../include/Matr.h"

using namespace std;

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

void printMatr(double** _mesh, int _rows, int _cols) {
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			cout << _mesh[i][j] << " ";
		}
		cout << endl;
	}
}
