//
// Created by ovv on 15.03.17.
//

#include "RK4.h"

#include <iostream>

void RK4(
    void (*function)(
        int dimension,
        double x0,
        vector<double> y,
        vector<double>& dy,
        int& exit_signal,
        unsigned int step_number),
    unsigned int dimension,
    int steps,
    double x0,
    double dx,
    vector<double> y,
    int& exit_signal,
    int& step_counter,
    Matrix& result)
{
    step_counter = 0;
    double x = x0;
    exit_signal = 0;

    vector<double> ak1(dimension, 0);
    vector<double> ak2(dimension, 0);
    vector<double> ak3(dimension, 0);
    vector<double> ak4(dimension, 0);
    vector<double> y1(dimension, 0);

    result = Matrix(steps + 1, vector<double>(dimension + 5, 0));

    result[0][0] = x;

    for (int i = 0; i < dimension; i++) {
        result[0][i + 1] = y[i];
    }

    vector<double> dy(dimension, 0);

    for (int i = 0; i < steps; i++) {
        function(dimension, x, y, dy, exit_signal, step_counter);
        if (exit_signal == 1) {
            break;
        }
        for (int j = 0; j < dimension; j++) {
            ak1[j] = dx * dy[j];
            y1[j] = y[j] + ak1[j] / 2.0;
        }

        double x1 = x + dx / 2.0;

        function(dimension, x1, y1, dy, exit_signal, step_counter);
        for (int j = 0; j < dimension; j++) {
            ak2[j] = dx * dy[j];
            y1[j] = y[j] + ak2[j] / 2.0;
        }

        function(dimension, x1, y1, dy, exit_signal, step_counter);
        for (int j = 0; j < dimension; j++) {
            ak3[j] = dx * dy[j];
            y1[j] = y[j] + ak3[j];
        }

        function(dimension, x + dx, y1, dy, exit_signal, step_counter);

        result[i + 1][0] = x + dx;

        for (int j = 0; j < dimension; j++) {
            ak4[j] = dx * dy[j];
            y[j] += (ak1[j] + 2.0 * ak2[j] + 2.0 * ak3[j] + ak4[j]) / 6.0;
            result[i + 1][j + 1] = y[j];
        }

        x += dx;
        step_counter++;
    }
}
