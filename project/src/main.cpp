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
	int m = 8, n = 5;  // A[m x n]
	int s = 4;  // B[n x s]
	double* A = NULL;
	double* B = NULL;
	// B is common for all of processes
	B = createMatrLine(n, s);
	// Initialization
	if (rank == 0) {
		cout << "B is" << endl;
		printMatr(B, n, s);
		// A is initialized only in root process
		A = createMatrLine(m, n);
		cout << "A is" << endl;
		printMatr(A, m, n);
		cout << "Correct answer is " << endl;
		printMatr(matrixMult(A, B, m, n, s), m, s);
	}
	double* C = matrixBlockMultParal(A, B, m, n, s, nOfCores);
	if (rank == 0) {
		cout << "Calculated answer is " << endl;
		printMatr(C, m, s);
	}
	delete[] A;
	delete[] B;
	delete[] C;
	MPI_Finalize();
	//system("pause");
}