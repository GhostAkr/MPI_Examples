#include <iostream>
#include "../include/mpiProduction.h"
#include "../include/Matr.h"
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	double t1 = MPI_Wtime();
	cout << "Hello, world!" << endl;
	int size = 0;
	int rank = 0;
	//double* M = new double[3];
	double** M = new double*[3];
	for (int i = 0; i < 3; i++)
		M[i] = new double[3];
	/*for (int j = 0; j < 3; j++) {
		M[j] = j + 1;
	}*/
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			M[i][j] = j + 3 * i;
	double** A = new double*[3];
	for (int i = 0; i < 3; i++)
		A[i] = new double[3];
	//double* A = new double[3];
	double D[3][3] = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	double C[3][3];
	int disp[4] = { 0, 1, 0, 1 };
	int counts[4] = { 2, 2, 2, 2 };
	int q = 0;
	MPI_Request send, rec;
	MPI_Status stat;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//MPI_Type_vector
	MPI_Scatterv(M, counts, disp, MPI_DOUBLE, A, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		cout << "size = " << size << endl;
		//for (int i = 0; i < 4; i++)
		//	for (int j = 0; j < 4; j++)
		//		M[i][j] = j + 4 * i;
		//for (int i = 0; i < 4; i++) {
		//	for (int j = 0; j < 4; j++)
		//		cout << M[i][j];
		//	cout << endl;
		//}
		//MPI_Isend(D, 4, MPI_DOUBLE, 1, 8, MPI_COMM_WORLD, &send);

	}
	if (rank == 1) {
		cout << "qergqerge" << endl;
		//MPI_Irecv(C, 4, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, &rec);
		//MPI_Test(&rec, &q, &stat);
		//MPI_Wait(&rec, &stat);
		//cout << A[0][0] << endl;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++)
			cout << A[i][j];
			//cout << A[j] << " ";
			cout << endl;
		}
	}
	//cout << C[0][0] << endl;
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++)
	//		cout << A[i][j];
	//	cout << endl;
	//}
	double t2 = MPI_Wtime();
	cout << "Time: " << t2 - t1 << endl;
	MPI_Finalize();
	system("pause");
}