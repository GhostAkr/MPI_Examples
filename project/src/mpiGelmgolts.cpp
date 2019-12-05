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
	int rows = (_yBorder - zeroPoint) / _step + 1;
	int cols = (_xBorder - zeroPoint) / _step + 1;
	*_rows = rows;
	*_cols = cols;
	double* resMesh = createMatrLineGelm(rows, cols);
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
