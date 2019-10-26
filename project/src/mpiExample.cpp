#include <iostream>
#include "../include/mpiExample.h"

using namespace std;

double* matrixMult(double* _source1, double* _source2, int m, int n, int s) {
	double* result = new double[m * s];
	for (int i = 0; i < m * s; ++i) {
		result[i] = 0.0;
	}
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < s; j++) {
				result[i * s + j] += _source1[i * n + k] * _source2[k * s + j];
			}
		}
	}
	return result;
}