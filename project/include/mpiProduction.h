#pragma once
#include "mpi.h"
#include "../include/Matr.h"

// Getting lines and cols for calculations

/*
	A * B = C

	   A|* * * * * *|  *   B|* x|  =  C|* *|  [* � x --- �����]
		|* * * * * *|		|* x|	   |* *|
		|* * * * * *|		|* x|	   |* *|	
		|x x x x x x|		|* x|      |* x|
		|* * * * * *|		|* x|      |* *|
		|* * * * * *|		|* x|      |* *|

	[nOfCores --- ����� ����]
	[A ������� �� nOfCores / 2 �����, B ������� �� 2 �������. ����� ������ ���� ������� ������������ ����� �������.]
	[height --- ������ ������� � A, width --- ������ ������� � B]
*/
double** ALine(double** _A, int _nOfLines, int _nOfCols, int _i, int _height);  // _i --- number of line (count from 1); _A[_nOfLines x _nOfCols]
double** BLine(double** _B, int _nOfCols, int _nOfLines, int _j, int _width);  // _j --- number of column (count from 1); _B[_nOfLines x _nOfCols]

// Block production

double** matrixBlockMult(double** _A, double** _B, int _m, int _n, int _s);  // _A[m x n] * _B[n x s]
double** matrixBlockMultParal(double** _A, double** _B, int _m, int _n, int _s, int nOfCores);

// Additional methods

bool compareMatrices(double** _source1, int m1, int n1, double** _source2, int m2, int n2);
