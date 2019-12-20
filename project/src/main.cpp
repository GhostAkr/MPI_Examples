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
	int m = 2064, n = 2064;  // A[m x n]
	int s = 2064;  // B[n x s]
	double* A = NULL;
	double* B = NULL;
	double* correctC = NULL;
	double tInit = 0., tFinal = 0.;
	// B is common for all of processes
	B = createMatrLine(n, s);
	// Initialization
	if (rank == 0) {
		// A is initialized only in root process
		A = createMatrLine(m, n);
		//tInit = MPI_Wtime();
		//correctC = matrixMult(A, B, m, n, s);
		//tFinal = MPI_Wtime();
		//cout << "Time spent with sequental production: " << tFinal - tInit << endl;
	}
	tInit = MPI_Wtime();
	double* C = matrixBlockMultParal(A, B, m, n, s, nOfCores);
	tFinal = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		/*if (compareMatrices(correctC, m, s, C, m, s)) {
			cout << "Answer is correct" << endl;
		}
		else {
			cout << "Answer is NOT correct" << endl;
		}*/
		cout << "Time spent with parallel production:" << tFinal - tInit << endl;
	}
	delete[] A;
	delete[] B;
	delete[] C;
	MPI_Finalize();
	//system("pause");
}