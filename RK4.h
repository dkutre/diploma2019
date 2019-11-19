//
// Created by ovv on 15.03.17.
//

#ifndef GAUSS_RK4_H
#define GAUSS_RK4_H

#include <vector>
#include "gauss.h"


void RK4(
    void (*function)(
        int dimension,
        double x0,
        vector<double> y,
        vector<double>& dy,
        int& exit_signal,
        unsigned int step_number),
    unsigned int dimension, int steps,
    double x0, double dx, vector<double> y,
    int& exit_signal, int& step_counter, Matrix& result);

#endif //GAUSS_RK4_H
