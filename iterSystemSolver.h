#ifndef GAUSS_ITERSYSTEMSOLVER_H
#define GAUSS_ITERSYSTEMSOLVER_H

#include <iostream>
#include <vector>
#include <string>
#include "SubFunctions.h"

using std::vector;

typedef vector<vector<double>> Matrix;

void solveSystem(const Matrix& LU, vector<double>& B, const vector<int>& T, vector<double>& X);

#endif //GAUSS_ITERSYSTEMSOLVER_H
