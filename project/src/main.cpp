#include <iostream>
#include "../include/mpiProduction.h"
#include "../include/Matr.h"
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int nOfCores = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &nOfCores);
	// Source matrices
	int m = 4, n = 4;  // A[m x n]
	int s = 4;  // B[n x s]
	double* A = NULL;
	double* B = NULL;
	// Initialization
	if (rank == 0) {
		A = createMatrLine(m, n);
		cout << "A is" << endl;
		printMatr(A, m, n);
		B = createMatrLine(n, s);
		cout << "B is" << endl;
		printMatr(A, n, s);
		cout << "Correct answer is " << endl;
		double* buffB = transpose(B, n, s);
		printMatr(matrixMult(A, buffB, m, n, s), m, s);
		delete[] buffB;
	}
	double* C = matrixBlockMultParal(A, B, m, n, s, nOfCores);
	if (rank == 0) {
		cout << "Calculated answer is " << endl;
		printMatr(C, m, s);
	}
	delete[] A;
	delete[] B;
	MPI_Finalize();
	//system("pause");
}