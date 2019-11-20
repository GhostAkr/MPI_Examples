#include <iostream>
#include "../include/mpiProduction.h"

using namespace std;

// Getting lines and cols for calculations

double** ALine(double** _A, int _nOfLines, int _nOfCols, int _i, int _height) {
	int rows = _height;
	int cols = _nOfCols;
	double** result = new double*[rows];
	for (int i = 0; i < rows; ++i) {
		result[i] = new double[cols];
		for (int j = 0; j < cols; ++j) {
			result[i][j] = _A[(_i - 1) * _height + i][j];
		}
	}
	return result;
}

double** BLine(double** _B, int _nOfCols, int _nOfLines, int _j, int _width) {
	int rows = _nOfLines;
	int cols = _width;
	double** result = new double*[rows];
	for (int i = 0; i < rows; ++i) {
		result[i] = new double[cols];
		for (int j = 0; j < cols; ++j) {
			result[i][j] = _B[i][(_j - 1) * _width + j];
		}
	}
	return result;
}

// Block production

double** matrixBlockMult(double** _A, double** _B, int _m, int _n, int _s) {
	int nOfCores = 4;  // In parallel version it's an argument
	int nOfLines = nOfCores / 2;
	int nOfCols = 2;
	int height = _m / nOfLines;
	int width = _s / nOfCols;
	double** result = new double*[_m];
	for (int i = 0; i < _m; ++i) {
		result[i] = new double[_s];
	}
	// Calculating
	for (int i = 1; i <= nOfLines; ++i) {
		for (int j = 1; j <= nOfCols; ++j) {
			double** block = matrixMult(ALine(_A, _m, _n, i, height), BLine(_B, _s, _n, j, width), height, _n, width);  // Calculating block [i][j]
			cout << "i = " << i << "; j = " << j << "; block is" << endl;
			printMatr(block, height, width);
			// Filling block [i][j]
			for (int k = 0; k < height; ++k) {
				for (int l = 0; l < width; ++l) {
					result[(i - 1) * height + k][(j - 1) * width + l] = block[k][l];
				}
			}
		}
	}
	return result;
}

// Additional methods

bool compareMatrices(double** _source1, int m1, int n1, double** _source2, int m2, int n2) {
	double epsNull = 1e-5;
	if (m1 != m2 || n1 != n2) {
		return false;
	}
	for (int i = 0; i < m1; ++i) {
		for (int j = 0; j < n1; ++j) {
			if (abs(_source1[i][j] - _source2[i][j]) > epsNull) {
				return false;
			}
		}
	}
	return true;
}
