#include "Classic_Newton.h"

void newtonClassic(
    void (*getResidual)(
        const vector<double> &,
        vector<double> &),
    unsigned int dimension,
    vector<double> & x,
    vector<double> residuals, int iteration) 
{

  //  std::cout << "NEWTON_START\n";
    double E = 1e-06;
    int max_iterations = 100; //M1
    iteration = 1; // NS1

    for (; iteration < max_iterations;) {
        getResidual(x, residuals);

        /* Проверка сходимости */
        int R1 = 0;
        for (int i = 0; i < dimension; i++) {
            if (fabs(residuals[i]) > E) {
                R1 = 1;
            }
        }
        if (R1 == 0) {
            return;
        }
        /* конец проверки сходимости */

        // что-то странное
        vector<double> B = residuals;
        for (int i = 0; i < B.size(); i++) {
            B[i] = -B[i];
        }
        //printContainer(B);
        //


        /* Вычисление матрицы Якоби */
        Matrix A(dimension, vector<double>(dimension, 0));
        vector<double> x_temp = x;
        for (int i = 0; i < dimension; i++) {
            double jStep = E * fabs(x_temp[i]); // шаг по i-ой переменной
            if (fabs(jStep) < 1e-10) {
                jStep = E;
            }
            x[i] = x_temp[i] + jStep;

            getResidual(x, residuals);

            // цикл по вычислению производных по невязкам
            for (int j = 0; j < dimension; j++) {
                A[j][i] = (residuals[j] + B[j]) / jStep;
            }
            x[i] = x_temp[i];
        }

        /* конец вычисления м. Якоби*/

        iteration++;
        // проверка сходимости по кол-ву итераций

        // решение СЛАУ
        Matrix LU = A;
        vector<int> T;
        string error = "";
        LUP(A, LU, T, error);
        if (!error.empty()) {}

        vector<double> dx(B.size(), 0);
        solveSystem(LU, B, T, dx);

        R1 = 0;
       // printContainer(x);
        for (int i = 0; i < dimension; i++) {
            x[i] += dx[i];
        }
    //    std::cout << "NEWTON_ITERATION " << iteration;
    //    printContainer(x);
    }
  //  std::cout << "NEWTON_END\n";
}
