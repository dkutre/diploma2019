
#ifndef GAUSS_CLASSIC_NEWTON_H
#define GAUSS_CLASSIC_NEWTON_H

#include "SubFunctions.h"
#include "gauss.h"
#include "iterSystemSolver.h"
#include <algorithm>



void newtonClassic(
    void (*getResidual) (
        const vector<double>&,
        vector<double>&),
    unsigned int dimension,
    vector<double>& x,
    vector<double> residuals,
    int iteration);

#endif //GAUSS_CLASSIC_NEWTON_H
