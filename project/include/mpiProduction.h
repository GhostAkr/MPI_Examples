#pragma once
#include "mpi.h"


double* matrixMult(double* _source1, double* _source2, int m, int n, int s);  // Usual matrix production. Sizes are [m x n] and [n x s]