#ifndef MPIGELMGOLTS_H
#define MPIGELMGOLTS_H

#include <cmath>
#define Pi 3.14159265358979323846264338327950288419716939937510582
//#define Pi 3.14
#include "mpi.h"
#include <iostream>
#include <fstream>

using namespace std;

// Additional

double* rightPart(double _step, int _rows, int _cols, double _k);
double f(double _x, double _y, double _k);
double* copyMesh(double* _mesh, int _rows, int _cols);
bool checkResult(double* _result, int _rows, int _cols, double _step);
double exactSolution(double _x, double _y);
void zeroLayer(double* _mesh, int _rows, int _cols);
double* createMesh(double _xBorder, double _yBorder, double _step, double* _rows, double* _cols);  // Step must be valid
double* createMeshFromNodes(double _xBorder, double _yBorder, double* _step, double* _rows, double* _cols, int _nOfNodes);
double* createMatrLineGelm(int _rows, int _cols);
double error(double* _result, int _rows, int _cols, double _step);

// Boundaries

void leftBoundary(double* _mesh, int _rows, int _cols, double _boundValue);
void rightBoundary(double* _mesh, int _rows, int _cols, double _boundValue);
void topBoundary(double* _mesh, int _rows, int _cols, double _boundValue);
void bottomBoundary(double* _mesh, int _rows, int _cols, double _boundValue);

// Sequental realizations

void Zeidel(double* _mesh, int _rows, int _cols, double _k, double _step, int ITERAT);
void Jacobi(double* _mesh, int _rows, int _cols, double _k, double _step, int ITERAT);

// Parallel realizations

void JacobiParall(double* _mesh, int _rows, int _cols, double _k, double _step, int ITERAT, const int nOfCores);
void ZeidelParall(double* _mesh, int _rows, int _cols, double _k, double _step, int ITERAT, const int nOfCores);

#endif  // MPIGELMGOLTS_H