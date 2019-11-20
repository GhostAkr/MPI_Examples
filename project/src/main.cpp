#include <iostream>
#include "../include/mpiProduction.h"
#include "../include/Matr.h"

using namespace std;

int main(int argc, char* argv[]) {
	double** tMatr1 = createMatr(4, 4);
	double** tMatr2 = createMatr(4, 4);
	double** res1 = matrixMult(tMatr1, tMatr2, 4, 4, 4);
	cout << "res 1" << endl;
	printMatr(res1, 4, 4);
	double** res2 = matrixBlockMult(tMatr1, tMatr2, 4, 4, 4);
	cout << "res 2" << endl;
	printMatr(res2, 4, 4);
	if (compareMatrices(res1, 4, 4, res2, 4, 4)) {
		cout << "Answer is correct" << endl;
	}
	else {
		cout << "Answer is INcorrect" << endl;
	}
	system("pause");
}