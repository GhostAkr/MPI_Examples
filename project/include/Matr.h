#pragma once

#include <iostream>

using namespace std;

// Creation & destroying

double** createMatr(int _rows, int _cols);
double* createMatrLine(int _rows, int _cols);
void deleteMatr(double** _matr, int _rows);

// Operations

double** matrixMult(double** _source1, double** _source2, int m, int n, int s);  // A[m x n] * B[n x s]
double* matrixMult(double* _source1, double* _source2, int m, int n, int s);  // For line-view matrix
double* transpose(double* _source, int m, int n);  // Transpose from A[m x n] to A^T[n x m]
// TODO: Make transpose function for double**

// Staff things

void printMatr(double** _mesh, int _rows, int _cols);
void printMatr(double* _mesh, int _rows, int _cols);

