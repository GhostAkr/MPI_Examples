#include <iostream>
#include "../include/mpiProduction.h"

using namespace std;

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	double t1 = MPI_Wtime();
	cout << "Hello, world!" << endl;
	double t2 = MPI_Wtime();
	cout << "Time: " << t2 - t1 << endl;
	MPI_Finalize();
	//system("pause");
}