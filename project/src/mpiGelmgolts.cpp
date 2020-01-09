#include "../include/mpiGelmgolts.h"
#include "../include/Matr.h"

void leftBoundary(double* _mesh, int _rows, int _cols, double _boundValue) {
	for (int i = 0; i < _rows; ++i) {
		_mesh[i * _cols + 0] = _boundValue;
	}
}
void rightBoundary(double* _mesh, int _rows, int _cols, double _boundValue) {
	for (int i = 0; i < _rows; ++i) {
		_mesh[i * _cols + (_cols - 1)] = _boundValue;
	}
}
void topBoundary(double* _mesh, int _rows, int _cols, double _boundValue) {
	for (int i = 0; i < _cols; ++i) {
		_mesh[(_rows - 1) * _cols + i] = _boundValue;
	}
}
void bottomBoundary(double* _mesh, int _rows, int _cols, double _boundValue) {
	for (int i = 0; i < _cols; ++i) {
		_mesh[0 * _cols + i] = _boundValue;
	}
}

double f(double _x, double _y, double _k) {
	return 2.0 * sin(Pi * _y) + _k * _k * (1 - _x) * _x * sin(Pi * _y) + \
		Pi * Pi * (1 - _x) * _x * sin(Pi * _y);
}

double* rightPart(double _step, int _rows, int _cols, double _k) {
	double* result = createMatrLine(_rows, _cols);
	double _step2 = _step * _step;
	for (int i = 1; i < _rows - 1; ++i) {
		for (int j = 1; j < _cols - 1; ++j) {
			result[i * _cols + j] = f(i * _step, j * _step, _k) * _step2;
		}
	}
	double boundaryValue = 0.0;
	leftBoundary(result, _rows, _cols, boundaryValue);
	rightBoundary(result, _rows, _cols, boundaryValue);
	topBoundary(result, _rows, _cols, boundaryValue);
	bottomBoundary(result, _rows, _cols, boundaryValue);
	return result;
}

double* copyMesh(double* _mesh, int _rows, int _cols) {
	double* result = createMatrLine(_rows, _cols);
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			result[i * _cols + j] = _mesh[i * _cols + j];
		}
	}
	return result;
}

bool checkResult(double* _result, int _rows, int _cols, double _step) {
	cout << "Step while checking = " << _step << endl;
	double epsNull = 1e-4;
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			if (fabs(_result[i * _cols + j] - exactSolution(i * _step, j * _step)) > epsNull) {
				cout << "j * step = " << j * _step << endl;
				cout << "Result[" << i << "][" << j << "] = " << _result[i * _cols + j] << "; Exact = " << exactSolution(i * _step, j * _step) << endl;
				return false;
			}
		}
	}
	return true;
}

double exactSolution(double _x, double _y) {
	return (1 - _x) * _x * sin(Pi * _y);
}

void zeroLayer(double* _mesh, int _rows, int _cols) {
	double zeroValue = 0.0;
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			_mesh[i * _cols + j] = zeroValue;
		}
	}
	double boundaryValue = 0.0;
	leftBoundary(_mesh, _rows, _cols, boundaryValue);
	rightBoundary(_mesh, _rows, _cols, boundaryValue);
	topBoundary(_mesh, _rows, _cols, boundaryValue);
	bottomBoundary(_mesh, _rows, _cols, boundaryValue);
}

double* createMesh(double _xBorder, double _yBorder, double _step, double* _rows, double* _cols) {
	double zeroPoint = 0.0;
	// Step must be valid. Here is no exception for it!
	int rows = (_yBorder - zeroPoint) / _step;
	int cols = (_xBorder - zeroPoint) / _step;
	cout << "Rows = " << rows << endl;
	cout << "Cols = " << cols << endl;
	*_rows = rows;
	*_cols = cols;
	double* resMesh = createMatrLineGelm(rows, cols);
	return resMesh;
}

double* createMeshFromNodes(double _xBorder, double _yBorder, double* _step, double* _rows, double* _cols, int _nOfNodes) {
	double zeroPoint = 0.0;
	*_step = (_yBorder - zeroPoint) / (_nOfNodes - 1);
	*_rows = _nOfNodes;
	*_cols = _nOfNodes;
	double* resMesh = createMatrLineGelm(*_rows, *_cols);
	return resMesh;
}

double* createMatrLineGelm(int _rows, int _cols) {
	double* result = new double[_rows * _cols];
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			result[i * _cols + j] = 0.0;
		}
	}
	return result;
}

void Jacobi(double* _mesh, int _rows, int _cols, double _k, double _step, int ITERAT) {
	double* rPart = rightPart(_step, _rows, _cols, _k);
	double c = 1 / (4 + _k * _k *_step*_step);
	double* previousLayer = copyMesh(_mesh, _rows, _cols);
	for (int s = 0; s < ITERAT; ++s) {
		if (s % 2 == 0) {
			for (int i = 1; i < _rows - 1; ++i) {
				for (int j = 1; j < _cols - 1; ++j) {
					_mesh[i * _cols + j] = c * (previousLayer[(i - 1) * _cols + j] + previousLayer[(i + 1) * _cols + j] + previousLayer[i * _cols + (j - 1)] + \
						previousLayer[i * _cols + (j + 1)] + rPart[i * _cols + j]);
				}
			}
		}
		/*if (s == 0)
			printMatr(_mesh, _rows, _cols);*/
		else {
			for (int i = 1; i < _rows - 1; ++i) {
				for (int j = 1; j < _cols - 1; ++j) {
					previousLayer[i * _cols + j] = c * (_mesh[(i - 1) * _cols + j] + _mesh[(i + 1) * _cols + j] + _mesh[i * _cols + (j - 1)] + \
						_mesh[i * _cols + (j + 1)] + rPart[i * _cols + j]);
				}
			}
		}
		/*if (checkResult(_mesh, _rows, _cols, _step)) {
		cout << "Accuracy was reached on " << s << " iteration" << endl;
		break;
		}*/
	}
	if (checkResult(_mesh, _rows, _cols, _step)) {
		cout << "Answer is correct" << endl;
	}
	else {
		cout << "Answer is INcorrect" << endl;
	}
	delete[] previousLayer;
	delete[] rPart;
}

void Zeidel(double* _mesh, int _rows, int _cols, double _k, double _step, int ITERAT) {
	double* rPart = rightPart(_step, _rows, _cols, _k);
	double c = 1.0 / (4.0 + _k * _k * _step * _step);
	double* previousLayer = copyMesh(_mesh, _rows, _cols);
	double* buff = nullptr;
	for (int s = 0; s <= ITERAT; ++s) {
		buff = previousLayer;
		previousLayer = _mesh;
		_mesh = buff;
		for (int i = 1; i < _rows - 1; i += 1) {
			for (int j = 2 - (i % 2); j < _cols - 1; j += 2) {
				_mesh[i * _cols + j] = c * (previousLayer[(i - 1) * _cols + j] + previousLayer[(i + 1) * _cols + j] + previousLayer[i * _cols + (j - 1)] + \
					previousLayer[i * _cols + (j + 1)] + rPart[i * _cols + j]);
			}
		}
		for (int i = 1; i < _rows - 1; i += 1) {
			for (int j = (i % 2) + 1; j < _cols - 1; j += 2) {
				_mesh[i * _cols + j] = c * (_mesh[(i - 1) * _cols + j] + _mesh[(i + 1) * _cols + j] + _mesh[i * _cols + (j - 1)] + \
					_mesh[i * _cols + (j + 1)] + rPart[i * _cols + j]);
			}
		}
		/*if (checkResult(_mesh, _rows, _cols, _step)) {
		cout << "Accuracy was reached on " << s << " iteration" << endl;
		break;
		}*/
	}
	if (checkResult(_mesh, _rows, _cols, _step)) {
		cout << "Answer is correct" << endl;
	}
	else {
		cout << "Answer is INcorrect" << endl;
	}
	// Buff points to _mesh or to previousLayer
	delete[] buff;
	delete[] rPart;
}



void JacobiParall(double* _mesh, int _rows, int _cols, double _k, double _step, int ITERAT, const int nOfCores) {
	double* rPart = rightPart(_step, _rows, _cols, _k);
	double c = 1. / (4. + _k * _k *_step*_step);
	MPI_Status stat;
	MPI_Request sendreq;
	int ireq;
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int ps = _rows / nOfCores;
	double tInit1 = 0., tInit2 = 0., tFinal1 = 0., tFinal2 = 0.;
	double tCalc = 0.;
	// Results from one core
	double* res = new double[_cols * (ps + 2)];
	for (int i = 0; i < _cols * (ps + 2); ++i) {
		res[i] = 0.;
	}
	// Buffer layer for one core
	double* BufLayer = new double[_cols * (ps + 2)];
	for (int i = 0; i < _cols * (ps + 2); ++i) {
		BufLayer[i] = 0.;
	}
	double* previousLayer = copyMesh(_mesh, _rows, _cols);
	// Tags for sending
	MPI_Request* reqsS = new MPI_Request[nOfCores];
	for (int i = 0; i < nOfCores; ++i) {
		reqsS[i] = MPI_REQUEST_NULL;
	}
	// Tags for recieving
	MPI_Request* reqsR = new MPI_Request[nOfCores];
	for (int i = 0; i < nOfCores; ++i) {
		reqsR[i] = MPI_REQUEST_NULL;
	}
	if (rank == 0) {
		MPI_Send_init(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &reqsS[rank]);
		MPI_Recv_init(res + ps * _cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &reqsR[rank]);
	}
	else if (rank < nOfCores - 1) {
		MPI_Send_init(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &reqsS[rank]);
		MPI_Recv_init(res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &reqsR[rank]);

		MPI_Send_init(res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &reqsS[rank]);
		MPI_Recv_init(res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &reqsR[rank]);
	}
	else {
		MPI_Send_init(res + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &reqsS[rank]);
		MPI_Recv_init(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &reqsR[rank]);
	}
	for (int s = 0; s < ITERAT; ++s) {
		// Even iterations
		if (s % 2 == 0) {
			tInit1 = MPI_Wtime();
			if (rank == 0) {  // First part
				for (int i = 1; i < ps; ++i) {
					for (int j = 1; j < _cols - 1; ++j) {
						res[i * _cols + j] = c * (BufLayer[(i - 1) * _cols + j] + BufLayer[(i + 1) * _cols + j] + BufLayer[i * _cols + (j - 1)] + \
							BufLayer[i * _cols + (j + 1)] + rPart[i * _cols + j]);
					}
				}
			}
			else if (rank == nOfCores - 1) {  // Last part
				for (int i = 2; i < ps + 1; ++i) {
					for (int j = 1; j < _cols - 1; ++j) {
						res[i * _cols + j] = c * (BufLayer[(i - 1) * _cols + j] + BufLayer[(i + 1) * _cols + j] + BufLayer[i * _cols + (j - 1)] + \
							BufLayer[i * _cols + (j + 1)] + rPart[((i - 2) + rank * ps) * _cols + j]);
					}
				}
			}
			else {
				for (int i = 1; i < ps + 1; ++i) {  // Other parts
					for (int j = 1; j < _cols - 1; ++j) {
						res[i * _cols + j] = c * (BufLayer[(i - 1) * _cols + j] + BufLayer[(i + 1) * _cols + j] + BufLayer[i * _cols + (j - 1)] + \
							BufLayer[i * _cols + (j + 1)] + rPart[((i - 1) + rank * ps) * _cols + j]);
					}
				}
			}
			tFinal1 = MPI_Wtime();
			tCalc += tFinal1 - tInit1;
			// Send-Recieve block
			//MPI_Barrier(MPI_COMM_WORLD);
			//if (rank == 06) {
			//	MPI_Sendrecv(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			//}
			//else if (rank == nOfCores - 1){
			//	MPI_Sendrecv(res + 2 *_cols, _cols, MPI_DOUBLE, rank - 1, 42, res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			//}
			//else {
			//    MPI_Sendrecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			//	MPI_Sendrecv(res + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			//}
			//MPI_Barrier(MPI_COMM_WORLD);
			//MPI_Startall(nOfCores, reqs);
			//MPI_Start(&reqs[rank]);
			/*MPI_Start(&reqs[rank]);
			MPI_Wait(&reqs[rank], &stat);*/
			if (rank == 0) {
				cout << "000" << endl;
				MPI_Start(&req1);
				//MPI_Wait(&req2, &stat);
				//MPI_Start(&req2);

			}
			else {
				cout << "111" << endl;
				MPI_Start(&req2);
				//MPI_Wait(&req1, &stat);
				//MPI_Start(&req1);
			}
			if (rank == 0){
				MPI_Request_free(&req1);
				MPI_Request_free(&req2);
				//MPI_Request_free(&req3);
				//MPI_Request_free(&req4);
			}
			if (rank == 1){
				//MPI_Request_free(&req1);
				//MPI_Request_free(&req2);
				MPI_Request_free(&req3);
				MPI_Request_free(&req4);
			}
			//if (rank == 0 && s == 0) {
			//	cout << "Before" << &reqsS[0] << endl;
			//}
			//cout << "Test2" << endl;
			//if (rank == 0) {
			//	MPI_Start(&reqsS[rank]);
			//	cout << "req0 = " << &reqsS[rank + 1] << endl;
			//	MPI_Wait(&reqsS[rank + 1], &stat);
			//	MPI_Start(&reqsR[rank]);
			//	for (int i = 0; i < nOfCores; ++i) {
			//		//reqsS[i] = MPI_REQUEST_NULL;
			//		//reqsR[i] = MPI_REQUEST_NULL;
			//		//MPI_Request_free(&reqsS[i]);
			//		//MPI_Request_free(&reqsR[i]);
			//	}
			//}
			//else if (rank < nOfCores - 1) {
			//	MPI_Start(&reqsS[rank]);
			//	MPI_Wait(&reqsS[rank - 1], &stat);
			//	MPI_Start(&reqsR[rank]);

			//	MPI_Start(&reqsS[rank]);
			//	MPI_Wait(&reqsS[rank + 1], &stat);
			//	MPI_Start(&reqsR[rank]);
			//	for (int i = 0; i < nOfCores; ++i) {
			//		//reqsS[i] = MPI_REQUEST_NULL;
			//		//reqsR[i] = MPI_REQUEST_NULL;
			//		//MPI_Request_free(&reqsS[i]);
			//		//MPI_Request_free(&reqsR[i]);
			//	}
			//}
			//else {
			//	MPI_Start(&reqsS[rank]);
			//	cout << "req1 = " << &reqsS[rank + 1] << endl;
			//	MPI_Wait(&reqsS[rank - 1], &stat);
			//	MPI_Start(&reqsR[rank]);
			//	for (int i = 0; i < nOfCores; ++i) {
			//		//reqsS[i] = MPI_REQUEST_NULL;
			//		//reqsR[i] = MPI_REQUEST_NULL;
			//		//MPI_Request_free(&reqsS[i]);
			//		//MPI_Request_free(&reqsR[i]);
			//	}
			//}
			//if (rank == 0 && s == 0) {
			//	cout << "After" << &reqsS[0] << endl;
			//}




			//MPI_Barrier(MPI_COMM_WORLD);
			//if (rank == 0 && s == 0) {
			//	//printMatr(res, ps + 2, _cols);
			//}

			//if (rank == 0) {
			//	MPI_Isend(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &sendreq);
			//	MPI_Irecv(res + ps * _cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &ireq);
			//}
			//else if (rank < nOfCores - 1) {
			//	MPI_Isend(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &sendreq);
			//	MPI_Irecv(res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &ireq);

			//	MPI_Isend(res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &sendreq);
			//	MPI_Irecv(res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &ireq);
			//}
			//else {
			//	MPI_Isend(res + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &sendreq);
			//	MPI_Irecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &ireq);
			//}
			//MPI_Barrier(MPI_COMM_WORLD);
		}
		// Odd iterations
		else {
			tInit2 = MPI_Wtime();
			if (rank == 0) {  // First part
				for (int i = 1; i < ps; ++i) {
					for (int j = 1; j < _cols - 1; ++j) {
						BufLayer[i * _cols + j] = c * (res[(i - 1) * _cols + j] + res[(i + 1) * _cols + j] + res[i * _cols + (j - 1)] + \
							res[i * _cols + (j + 1)] + rPart[i * _cols + j]);
					}
				}
			}
			else if (rank == nOfCores - 1) {  // Last part
				for (int i = 2; i < ps + 1; ++i) {
					for (int j = 1; j < _cols - 1; ++j) {
						BufLayer[i * _cols + j] = c * (res[(i - 1) * _cols + j] + res[(i + 1) * _cols + j] + res[i * _cols + (j - 1)] + \
							res[i * _cols + (j + 1)] + rPart[((i - 2) + rank * ps) * _cols + j]);
					}
				}
			}
			else {  // Other parts
				for (int i = 1; i < ps + 1; ++i) {
					for (int j = 1; j < _cols - 1; ++j) {
						BufLayer[i * _cols + j] = c * (res[(i - 1) * _cols + j] + res[(i + 1) * _cols + j] + res[i * _cols + (j - 1)] + \
							res[i * _cols + (j + 1)] + rPart[((i - 1) + rank * ps) * _cols + j]);
					}
				}
			}
			tFinal2 = MPI_Wtime();
			tCalc += tFinal2 - tInit2;
			// Send-Recieve block
			//MPI_Barrier(MPI_COMM_WORLD);
			//if (rank == 0) {
			//	MPI_Sendrecv(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			//}
			//else if (rank == nOfCores - 1){
			//	MPI_Sendrecv(res + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			//}
			//else {
			//	MPI_Sendrecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			//	MPI_Sendrecv(res + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			//}
			//MPI_Barrier(MPI_COMM_WORLD);
			//cout << "Test" << endl;
			if (rank == 0) {
				MPI_Isend(BufLayer + (ps - 1) *_cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(BufLayer + ps * _cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &ireq);
			}
			else if (rank < nOfCores - 1) {
				MPI_Isend(BufLayer + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(BufLayer, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &ireq);
				//MPI_Barrier(MPI_COMM_WORLD);
				MPI_Isend(BufLayer + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(BufLayer + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &ireq);
			}
			else {
				MPI_Isend(BufLayer + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(BufLayer + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &ireq);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	// Constructing result
	double* resBuf = NULL;
	if (rank == 0) {
		resBuf = BufLayer;
	}
	else if (rank == nOfCores - 1) {
		resBuf = BufLayer + 2 * _cols;
	}
	else {
		resBuf = BufLayer + _cols;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(resBuf, ps * _cols, MPI_DOUBLE, previousLayer, ps * _cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0) {  // Check answer
		//printMatr(previousLayer, _rows, _cols);
		if (checkResult(previousLayer, _rows, _cols, _step)) {
			cout << "Answer is correct" << endl;
		}
		else {
			cout << "Answer is INcorrect" << endl;
		}
		cout << "Error = " << error(previousLayer, _rows, _cols, _step) << endl;
		cout << "Calculation time: " << tCalc << endl;
	}
	if (rank == 0) {
		cout << "Writing to file" << endl;
		ofstream file("C:\\Users\\ighos\\Documents\\Programming\\Univercity\\MPI_Examples\\fileJac.dat");
		if (!file) {
			cout << "Error while opening file" << endl;
		}
		for (int i = 1; i < _rows; ++i) {
			for (int j = 1; j < _cols - 1; ++j) {
				file << previousLayer[i * _cols + j] << " ";
			}
			file << endl;
		}
		file.close();
	}
	delete[] previousLayer;
	delete[] rPart;
}

void ZeidelParall(double* _mesh, int _rows, int _cols, double _k, double _step, int ITERAT, const int nOfCores) {
	double* rPart = rightPart(_step, _rows, _cols, _k);
	double c = 1 / (4 + _k * _k *_step*_step);
	MPI_Status stat;
	MPI_Request sendreq;
	int ireq;
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int ps = _rows / nOfCores;
	double tInit1 = 0., tInit2 = 0., tFinal1 = 0., tFinal2 = 0.;
	double tCalc = 0.;
	// Results from one core
	double* res = new double[_cols * (ps + 2)];
	for (int i = 0; i < _cols * (ps + 2); ++i) {
		res[i] = 0.;
	}
	// Buffer layer for one core
	double* BufLayer = new double[_cols * (ps + 2)];
	for (int i = 0; i < _cols * (ps + 2); ++i) {
		BufLayer[i] = 0.;
	}
	double* previousLayer = copyMesh(_mesh, _rows, _cols);
	for (int s = 0; s < ITERAT; ++s) {
		// Even iterations
		if (s % 2 == 0) {
			tInit1 = MPI_Wtime();
			if (rank == 0) {  // First part
				for (int i = 1; i < ps; ++i) {
					for (int j = 1; j < _cols - 1; j += 2) {
						res[i * _cols + j] = c * (BufLayer[(i - 1) * _cols + j] + BufLayer[(i + 1) * _cols + j] + BufLayer[i * _cols + (j - 1)] + \
							BufLayer[i * _cols + (j + 1)] + rPart[i * _cols + j]);
					}
				}
			}
			else if (rank == nOfCores - 1) {  // Last part
				for (int i = 2; i < ps + 1; ++i) {
					for (int j = 1; j < _cols - 1; j += 2) {
						res[i * _cols + j] = c * (BufLayer[(i - 1) * _cols + j] + BufLayer[(i + 1) * _cols + j] + BufLayer[i * _cols + (j - 1)] + \
							BufLayer[i * _cols + (j + 1)] + rPart[((i - 2) + rank * ps) * _cols + j]);
					}
				}
			}
			else {
				for (int i = 1; i < ps + 1; ++i) {  // Other parts
					for (int j = 1; j < _cols - 1; j += 2) {
						res[i * _cols + j] = c * (BufLayer[(i - 1) * _cols + j] + BufLayer[(i + 1) * _cols + j] + BufLayer[i * _cols + (j - 1)] + \
							BufLayer[i * _cols + (j + 1)] + rPart[((i - 1) + rank * ps) * _cols + j]);
					}
				}
			}
			tFinal1 = MPI_Wtime();
			tCalc += tFinal1 - tInit1;
			BufLayer = res;
			//BufLayer = copyMesh(res, ps + 2, _cols);
			// Send-Recieve block
			//MPI_Barrier(MPI_COMM_WORLD);
			//if (rank == 0) {
			//	MPI_Sendrecv(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			//}
			//else if (rank == nOfCores - 1){
			//	MPI_Sendrecv(res + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			//}
			//else {
			//	MPI_Sendrecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			//	MPI_Sendrecv(res + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row

			//}
			if (rank == 0) {
				MPI_Isend(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(res + ps * _cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &ireq);
			}
			else if (rank < nOfCores - 1){
				MPI_Isend(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &ireq);

				MPI_Isend(res + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &ireq);
			}
			else {
				MPI_Isend(res + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &ireq);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		// Odd iterations
		else {
			//MPI_Barrier(MPI_COMM_WORLD);
			tInit2 = MPI_Wtime();
			if (rank == 0) {  // First part
				for (int i = 1; i < ps; ++i) {
					for (int j = 2; j < _cols - 1; j += 2) {
						BufLayer[i * _cols + j] = c * (res[(i - 1) * _cols + j] + res[(i + 1) * _cols + j] + res[i * _cols + (j - 1)] + \
							res[i * _cols + (j + 1)] + rPart[i * _cols + j]);
					}
				}
			}
			else if (rank == nOfCores - 1) {  // Last part
				for (int i = 2; i < ps + 1; ++i) {
					for (int j = 2; j < _cols - 1; j += 2) {
						BufLayer[i * _cols + j] = c * (res[(i - 1) * _cols + j] + res[(i + 1) * _cols + j] + res[i * _cols + (j - 1)] + \
							res[i * _cols + (j + 1)] + rPart[((i - 2) + rank * ps) * _cols + j]);
					}
				}
			}
			else {
				for (int i = 1; i < ps + 1; ++i) {  // Other parts
					for (int j = 2; j < _cols - 1; j += 2) {
						BufLayer[i * _cols + j] = c * (res[(i - 1) * _cols + j] + res[(i + 1) * _cols + j] + res[i * _cols + (j - 1)] + \
							res[i * _cols + (j + 1)] + rPart[((i - 1) + rank * ps) * _cols + j]);
					}
				}
			}
			tFinal2 = MPI_Wtime();
			tCalc += tFinal2 - tInit2;
			// Send-Recieve block
			//MPI_Barrier(MPI_COMM_WORLD);
			//if (rank == 0) {
			//	MPI_Sendrecv(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			//}
			//else if (rank == nOfCores - 1){
			//	MPI_Sendrecv(res + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			//}
			//else {
			//	MPI_Sendrecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			//	MPI_Sendrecv(res + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			//}
			if (rank == 0) {
				MPI_Isend(BufLayer + (ps - 1) *_cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(BufLayer + ps * _cols, _cols, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, &ireq);
			}
			else if (rank < nOfCores - 1) {
				MPI_Isend(BufLayer + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(BufLayer, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &ireq);
				//MPI_Barrier(MPI_COMM_WORLD);
				MPI_Isend(BufLayer + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(BufLayer + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &ireq);
			}
			else {
				MPI_Isend(BufLayer + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &sendreq);
				MPI_Irecv(BufLayer + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &ireq);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	//cout << "Test" << endl;
	// Constructing result
	double* resBuf = NULL;
	if (rank == 0) {
		resBuf = BufLayer;
	}
	else if (rank == nOfCores - 1) {
		resBuf = BufLayer + 2 * _cols;
	}
	else  {
		resBuf = BufLayer + _cols;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(resBuf, ps * _cols, MPI_DOUBLE, previousLayer, ps * _cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		//printMatr(previousLayer, _rows, _cols);
		if (checkResult(previousLayer, _rows, _cols, _step)) {
			cout << "Answer is correct" << endl;
		}
		else {
			cout << "Answer is INcorrect" << endl;
		}
		cout << "Error = " << error(previousLayer, _rows, _cols, _step) << endl;
		cout << "Calculation time: " << tCalc << endl;
	}
	if (rank == 0) {
		cout << "Writing to file" << endl;
		ofstream file("C:\\Users\\ighos\\Documents\\Programming\\Univercity\\MPI_Examples\\fileZeid.dat");
		if (!file) {
			cout << "Error while opening file" << endl;
		}
		for (int i = 1; i < _rows; ++i) {
			for (int j = 1; j < _cols - 1; ++j) {
				file << previousLayer[i * _cols + j] << " ";
			}
			file << endl;
		}
		file.close();
	}
	delete[] previousLayer;
	delete[] rPart;
}

double error(double* _result, int _rows, int _cols, double _step) {
	double max = 0.;
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			double diff = fabs(_result[i * _cols + j] - exactSolution(i * _step, j * _step));
			if (diff > max) {
				max = diff;
			}
		}
	}
	return max;
}
