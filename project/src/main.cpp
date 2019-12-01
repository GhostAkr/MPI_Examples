#include <iostream>
#include "../include/mpiProduction.h"
#include "../include/Matr.h"
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int nOfCores;
	MPI_Comm_size(MPI_COMM_WORLD, &nOfCores);
	if (rank == 0) {
		int m = 2, n = 2;  // A[m x n]
		int s = 2;  // B[n x s]
		double* A = createMatrLine(m, n);
		double* B = createMatrLine(n, s);
		double* C = matrixBlockMultParal(A, B, m, n, s, nOfCores);
		cout << "C is " << endl;
		printMatr(C, m, s);
	}
	MPI_Finalize();
	system("pause");
}