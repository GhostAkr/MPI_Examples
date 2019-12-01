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

double* ALine(double* _A, int _nOfLines, int _nOfCols, int _i, int _height) {
	int rows = _height;
	int cols = _nOfCols;
	double* result = new double[rows * cols];
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			result[i * cols + j] = _A[((_i - 1) * _height + i) * _nOfCols + j];
		}
	}
	return result;
}

double* BLine(double* _B, int _nOfCols, int _nOfLines, int _j, int _width) {
	int rows = _nOfLines;
	int cols = _width;
	double* result = new double[rows * cols];
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			result[i * cols + j] = _B[i * _nOfCols + ((_j - 1) * _width + j)];
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

double* matrixBlockMultParal(double* _A, double* _B, int _m, int _n, int _s, const int nOfCores){
	int nOfCols = 2;
	int nOfLines = nOfCores / nOfCols;
	const int nOfCoresArg = 2;  // Need to set number of cores manually
	/*
		B matrix is transposed from the beginning
	*/
	int height = _m / nOfLines;   
	int width = _n / nOfCols;
	int targetBlockSize = height * width;
	int ALineSize = height * _n;
	int BLineSize = width * _n;
	cout << "ALineSize = " << ALineSize << endl;
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double* ABlock = new double[ALineSize];
	double* BBlock = new double[BLineSize];
	double* result = new double[_m * _s];
	MPI_Scatter(_A, ALineSize, MPI_DOUBLE, ABlock, ALineSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(_B, BLineSize, MPI_DOUBLE, BBlock, BLineSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	double* buffBBlock = transpose(BBlock, width, _n);
	double* block = matrixMult(ABlock, buffBBlock, height, _n, width);
	/*if (rank == 0) {
		cout << "Now on " << rank << " rank" << endl;
		cout << "Calculated block is" << endl;
		printMatr(block, height, width);
		cout << "ABlock is" << endl;
		printMatr(ABlock, height, _n);
		cout << "BBlock is" << endl;
		printMatr(BBlock, width, _n);
	}*/
	/*if (rank == 1) {
		cout << "Now on " << rank << " rank" << endl;
		cout << "Calculated block is" << endl;
		printMatr(block, height, width);
		cout << "ABlock is" << endl;
		printMatr(ABlock, height, _n);
		cout << "BBlock is" << endl;
		printMatr(BBlock, width, _n);
	}*/
	if (rank == 2) {
		cout << "Now on " << rank << " rank" << endl;
		cout << "Calculated block is" << endl;
		printMatr(block, height, width);
		cout << "ABlock is" << endl;
		printMatr(ABlock, height, _n);
		cout << "BBlock is" << endl;
		printMatr(BBlock, width, _n);
	}
	/*if (rank == 3) {
		cout << "Now on " << rank << " rank" << endl;
		cout << "Calculated block is" << endl;
		printMatr(block, height, width);
		cout << "ABlock is" << endl;
		printMatr(ABlock, height, _n);
		cout << "BBlock is" << endl;
		printMatr(BBlock, width, _n);
	}*/
	delete[] buffBBlock;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(block, targetBlockSize, MPI_DOUBLE, result, targetBlockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete[] block;
	delete[] ABlock;
	delete[] BBlock;
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
