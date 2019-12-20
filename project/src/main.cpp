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
		int nproc = 0;
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		double tInit = 0.;
		double tFinal = 0.;

		// Init block
		double rows = 0.0;
		double cols = 0.0;
		double step = 0.1;
		int nOfNodes = 100;
		double* mesh1 = NULL;
		double* mesh2 = NULL;
		double t1 = 0.0, t2 = 0.0;
		int nOfIters = 15000;

		// Assigning block
		//mesh1 = createMesh(1, 1, step, &rows, &cols);
		mesh1 = createMeshFromNodes(1, 1, &step, &rows, &cols, nOfNodes);

		// Calculations block

		tInit = MPI_Wtime();
		// Jacobi
		//Jacobi(mesh1, rows, cols, 150, step, nOfIters);
		JacobiParall(mesh1, rows, cols, 150, step, nOfIters,nproc);
		//printMatr(mesh1, rows, cols);

		// Zeidel
		//Zeidel(mesh1, rows, cols, 1, step, nOfIters);
		//ZeidelParall(mesh1, rows, cols, 150, step, nOfIters, nproc);
		//printMatr(mesh1, rows, cols);
		tFinal = MPI_Wtime();
		if (rank == 0) {
			double deltaT = tFinal - tInit;
			cout << "Total time: " << deltaT << endl;
		}

		// Cleanup
		delete[] mesh1;
		//system("pause");  // For sequental methods
		MPI_Finalize();
	}

}