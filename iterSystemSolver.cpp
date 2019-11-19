//
// Created by ovv on 26.02.17.
//

#include "iterSystemSolver.h"


void solveSystem(const Matrix& LU, vector<double>& B, const vector<int>& T, vector<double>& X) {
    // умножаем T на B
    size_t n = B.size() - 1;
    vector<double> B_temp = B;
    for (size_t i = 0; i < T.size(); i++) {
        B[i] = B_temp[T[i]];
    }


    vector<double> y(B.size(), 0);
    y[0] = B[0];
    for (int i = 1; i <= n; i++) {
        double sum = 0;
        for (int j = 0; j <= i - 1; j++) {
            sum += LU[i][j] * y[j];
        }
        y[i] = B[i] - sum;
    }

    X[n] = y[n] / LU[n][n];
    for (int i = n - 1; i >= 0 ; i--) {
        double sum = 0;
        for (int j = i + 1; j <= n; j++) {
            sum += LU[i][j] * X[j];
        }
        X[i] = (y[i] - sum) / LU[i][i];
    }
}
