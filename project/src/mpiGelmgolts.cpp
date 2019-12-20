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
	double epsNull = 1e-2;
	for (int i = 0; i < _rows; ++i) {
		for (int j = 0; j < _cols; ++j) {
			if (fabs(_result[i * _cols + j] - exactSolution(i * _step, j * _step)) > epsNull) {
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
	double c = 1 / (4 + _k * _k *_step*_step);
	MPI_Status stat;
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int ps = _rows / nOfCores;
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
			// Send-Recieve block
			MPI_Barrier(MPI_COMM_WORLD);
			if (rank == 0) {
				MPI_Sendrecv(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			}
			else if (rank == nOfCores - 1){
				MPI_Sendrecv(res + 2 *_cols, _cols, MPI_DOUBLE, rank - 1, 42, res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			}
			else {
				MPI_Sendrecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
				MPI_Sendrecv(res + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
				
			}
			MPI_Barrier(MPI_COMM_WORLD);
		} 
		// Odd iterations
		else {
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
			// Send-Recieve block
			MPI_Barrier(MPI_COMM_WORLD);
			if (rank == 0) {
				MPI_Sendrecv(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			}
			else if (rank == nOfCores - 1){
				MPI_Sendrecv(res + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			}
			else {
				MPI_Sendrecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
				MPI_Sendrecv(res + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
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
	else  {
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
	}
	delete[] previousLayer;
	delete[] rPart;
}

void ZeidelParall(double* _mesh, int _rows, int _cols, double _k, double _step, int ITERAT, const int nOfCores) {
	double* rPart = rightPart(_step, _rows, _cols, _k);
	double c = 1 / (4 + _k * _k *_step*_step);
	MPI_Status stat;
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int ps = _rows / nOfCores;
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
			BufLayer = copyMesh(res, ps + 2, _cols);
			// Send-Recieve block
			MPI_Barrier(MPI_COMM_WORLD);
			if (rank == 0) {
				MPI_Sendrecv(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			}
			else if (rank == nOfCores - 1){
				MPI_Sendrecv(res + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			}
			else {
				MPI_Sendrecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
				MPI_Sendrecv(res + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row

			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		// Odd iterations
		else {
			MPI_Barrier(MPI_COMM_WORLD);
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
			// Send-Recieve block
			MPI_Barrier(MPI_COMM_WORLD);
			if (rank == 0) {
				MPI_Sendrecv(res + (ps - 1) *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + ps * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
			}
			else if (rank == nOfCores - 1){
				MPI_Sendrecv(res + 2 * _cols, _cols, MPI_DOUBLE, rank - 1, 42, res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
			}
			else {
				MPI_Sendrecv(res + _cols, _cols, MPI_DOUBLE, rank - 1, 42, res, _cols, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &stat);  // First row
				MPI_Sendrecv(res + ps *_cols, _cols, MPI_DOUBLE, rank + 1, 42, res + (ps + 1) * _cols, _cols, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD, &stat);  // Last row
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
	}
	delete[] previousLayer;
	delete[] rPart;
}
