#include <iostream>
#include "../include/mpiProduction.h"
#include "../include/Matr.h"

using namespace std;

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	double** m1 = new double* [2];
	double** m2 = new double* [2];
	for (int i = 0; i < 2; ++i) {
		m1[i] = new double[2];
		m2[i] = new double[2];
	}
	m1[0][0] = 1; m1[0][1] = 2; m1[1][0] = 3; m1[1][1] = 4;
	m2[0][0] = 5; m2[0][1] = 6; m2[1][0] = 7; m2[1][1] = 8;
	cout << "m1" << endl;
	printMatr(m1, 2, 2);
	cout << "m2" << endl;
	printMatr(m2, 2, 2);
	cout << "m1 * m2" << endl;
	double** res = matrixMult(m1, m2, 2, 2, 2);
	printMatr(res, 2, 2);
	MPI_Finalize();
	system("pause");
}