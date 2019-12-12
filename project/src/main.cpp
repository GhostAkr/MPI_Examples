#include <iostream>
#include "../include/mpiProduction.h"
#include "../include/Matr.h"
#include "../include/mpiGelmgolts.h"
#include "mpi.h"

#define PRODUCTION 0
#define GELMGOLTS 1

using namespace std;
 
int main(int argc, char* argv[]) {
	if (PRODUCTION) {
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
	if (GELMGOLTS) {
		MPI_Init(&argc, &argv);
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		// Init block
		double rows = 0.0;
		double cols = 0.0;
		double step = 0.1;
		double* mesh1 = NULL;
		double* mesh2 = NULL;
		double t1 = 0.0, t2 = 0.0;
		int nOfIters = 15000;
		// Assigning block
		mesh1 = createMesh(1, 1, step, &rows, &cols);
		//cout << "Rows = " << rows << endl;
		//cout << "Cols = " << cols << endl;
		// Jacobi
		
		//cout << "Jacobi" << endl;
		JacobiParall(mesh1, rows, cols, 0, step, nOfIters,2);

		// Zeidel
		//cout << "Zeidel" << endl;
		//Zeidel(mesh2, rows, cols, 1, step, nOfIters);

		// Cleanup
		delete[] mesh1;
		delete[] mesh2;
		//system("pause");
		MPI_Finalize();
	}

}