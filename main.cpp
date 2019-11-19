#include "gauss.h"
#include "SubFunctions.h"
#include "iterSystemSolver.h"
#include "RK4.h"
#include "Classic_Newton.h"
#include <cmath>
#include <algorithm>

double H = 2000.0;
double V = 100.0;
double am = 20000.0;
double g = 9.81;
double S = 10.75;
double P1 =  -0.15;

double hkm = H / 1000;
double ro = 1.22235514 - 0.11422776 * hkm + 0.00368329 * hkm * hkm - 0.0000404 * hkm * hkm * hkm;
double Xaer = 0.8 * ro * V * V * S / 2.0;

double P2 = P1 * V / (g - Xaer / am);
double P3 = -0.7;
double T = 10.0;
double u_max = 305.0;
double beta = 2766.0;
double U = 0.0;
double U_switch = 0;


int rkStepsFirst  = 0; // N2
int rkStepsSecond = 0; // N3
int rungeKuttaSteps = 1000;

double Hpontr = 0;
double Fper = 0;
double tabs = 0;

Matrix RESULT(1000, vector<double>(12, 0));
Matrix OPTIMAL(1000, vector<double>(12, 0));

void FUN2 (int dimension, double x0, vector<double> y, vector<double>& dy, int& exit_signal, unsigned int step_number)
{
    if (rkStepsFirst == 0) {
        U = 0;
    } else if (rkStepsFirst != 0 && rkStepsSecond == rungeKuttaSteps - step_number - rkStepsFirst) {
        U = U_switch;
    } else if (rkStepsFirst != 0 && rkStepsSecond != rungeKuttaSteps - step_number - rkStepsFirst) {
        U = u_max;
    }

    hkm = y[0] / 1000;
    ro = 1.22235514 - 0.11422776 * hkm + 0.00368329 * hkm * hkm - 0.0000404 * hkm * hkm * hkm;
    Xaer = 0.8 * ro * y[1] * y[1] * S / 2.0;


    double drodH = (-0.11422776 + 2.0 * 0.00368329 * hkm - 3.0 * 0.0000404 * hkm * hkm) * 0.001;
    dy[0] = -y[1] * T ;
    dy[1] = (-(beta * U + Xaer) / y[2] + g) * T;
    dy[2] = - U * T;
    dy[3] = (y[4] * 0.8 * y[1] * y[1] * S * drodH * T) / (2.0 * y[2]);
    dy[4] = (y[3] + y[4] * Xaer * 2.0 / y[1] / y[2]) * T;
    dy[5] = -y[4] * (beta * U + Xaer) * T / (y[2] * y[2]);

    if (fabs(1.0 / rungeKuttaSteps * (step_number + rkStepsFirst) - x0) < 1e-10) {
        Hpontr = (-y[3] * y[1] + y[4] * (g - (beta * U + Xaer) / y[2]) - y[5] * U - U) * T;
        Fper = -y[4] * beta / y[2] - y[5] - 1.0;
        tabs = x0 * T;

        RESULT[rkStepsFirst + step_number][dimension + 1] = Hpontr;// ! запись дополнительных переменных
        RESULT[rkStepsFirst + step_number][dimension + 2] = Fper;// ! запись дополнительных переменных
        RESULT[rkStepsFirst + step_number][dimension + 3] = U;// ! запись дополнительных переменных
        RESULT[rkStepsFirst + step_number][dimension + 4] = tabs;// ! запись дополнительных переменных
        RESULT[rkStepsFirst + step_number][dimension + 5] = Xaer;// ! запись дополнительных переменных


        if (rkStepsSecond == 0 && Fper > 0) {
            exit_signal = 1;
        }
        if (step_number >= rkStepsSecond && rkStepsSecond != 0) {
            exit_signal = 1;
        }
    }
}

void resetGlobalVariables()
{
    RESULT.clear();
    RESULT = Matrix(1001, vector<double>(12, 0));
    rkStepsFirst = 0;
    rkStepsSecond = 0;

}

void getResidual (const vector<double> & x, vector<double>& residual)
{
    resetGlobalVariables();

    int dimension = 6;
    rkStepsFirst  = 0; // N2
    rkStepsSecond = 0; // N3

    vector<double> y(6, 0);
    Matrix result(rungeKuttaSteps + 1, vector<double>(7, 0));
    int exit_signal;
    residual = vector<double>(x.size(), 0);


    hkm = H / 1000;
    ro = 1.22235514- 0.11422776 * hkm + 0.00368329 * hkm * hkm - 0.0000404 * hkm * hkm * hkm;
    Xaer = 0.8 * ro * V * V * S / 2.0;


    double x0 = 0.0;
    double dx = 1.0 / rungeKuttaSteps;
    P1 = x[0];
    P3 = x[1];
    //!P2=(P1*V+P3*U+U)/(g-(b*U+Xaer)/am)
    P2 = (P1 * V) / (g - Xaer / am);
    //! if (XX(3).Lt.5.d0) XX(3)=5.d0
    T = x[2];

    y[0] = H;
    y[1] = V;
    y[2] = am;
    y[3] = P1;
    y[4] = P2;
    y[5] = P3;

    int temp_step = 0;
    RK4(FUN2, dimension, rungeKuttaSteps, x0, dx, y, exit_signal, temp_step, result);
    rkStepsFirst = temp_step;

    for (int i = 0; i < rkStepsFirst + 1; i++) {
        for (int j = 0; j < dimension + 1; j++) {
            RESULT[i][j] = result[i][j];
        }
    }


    double t_switch = result[rkStepsFirst - 1][0]
                      - RESULT[rkStepsFirst - 1][dimension + 2]
                        * (result[rkStepsFirst][0] - result[rkStepsFirst - 1][0])
                        / (RESULT[rkStepsFirst][dimension + 2] - RESULT[rkStepsFirst - 1][dimension + 2]);

    double ksi = (t_switch - result[rkStepsFirst - 1][0]) / (result[rkStepsFirst][0] - result[rkStepsFirst - 1][0]);
    U_switch = u_max * (1.0 - ksi);

    y[0] = result[rkStepsFirst - 1][1];
    y[1] = result[rkStepsFirst - 1][2];
    y[2] = result[rkStepsFirst - 1][3];
    y[3] = result[rkStepsFirst - 1][4];
    y[4] = result[rkStepsFirst - 1][5];
    y[5] = result[rkStepsFirst - 1][6];

    rkStepsSecond = rungeKuttaSteps - rkStepsFirst;
    x0 = result[rkStepsFirst - 1][0];

    // второй вызов РК
    RK4(FUN2, dimension, rungeKuttaSteps, x0, dx, y, exit_signal, temp_step, result);

    for (int i = rkStepsFirst - 1; i < rungeKuttaSteps + 1; i++) {
        for (int j = 0; j < dimension + 1; j++) {
            RESULT[i][j] = result[i - rkStepsFirst + 1][j];
        }
    }


    residual[0] = RESULT[rungeKuttaSteps][1];
    residual[1] = RESULT[rungeKuttaSteps][2];
    residual[2] = RESULT[rungeKuttaSteps][6];
}


void setInitalValue()
{
    Hpontr = (-P1 * V + P2 * (g - (beta * U + Xaer) / am) - P3 * U - U) * T;
    Fper = -P2 * beta / am - P3 - 1.0;

    RESULT[0][7] = Hpontr;
    RESULT[0][8] = Fper;
    RESULT[0][9] = U;
    RESULT[0][10] = 0; //tabs
    RESULT[0][11] = Xaer;
}

double findSwitchMoment()
{
    for (auto elem : RESULT) {
        if (elem[8] > 0) {
            return elem[10];
        }
    }
}

void massAnalyz()
{
    Matrix correlation = Matrix(1000, vector<double>(3, 0)); // зависимость масса время посадки и расход топлива

    int newton_steps = 0;
    vector<double> residuals(3, 1);
    vector<double> args = {P1, P3, T};
    newtonClassic(&getResidual, 3, args, residuals, newton_steps);
    correlation[0][0] = RESULT[0][3]; // записываем массу
    correlation[0][1] = RESULT[1000][10]; // время посадки
    correlation[0][2] = RESULT[0][3] - RESULT[1000][3]; // затраченное топливо

    printMatrix(RESULT);
    for (size_t i = 1; i < correlation.size(); i++) {
        H += 4;
        setInitalValue();
        newton_steps = 0;
        newtonClassic(&getResidual, 3, args, {1, 1, 1}, newton_steps);
        correlation[i][0] = RESULT[0][1]; // записываем массу
        correlation[i][1] = findSwitchMoment(); // момент переключения
        correlation[i][2] = RESULT[0][3] - RESULT[1000][3]; // затраченное топливо
        if (i % 100 == 0) {
            std::cout << "ke" << std::endl;
        }
        //printMatrix(RESULT);
    }

    printMatrix(correlation);
}


void FUN2_2 (int dimension, double x0, vector<double> y, vector<double>& dy, int& exit_signal, unsigned int step_number)
{
    hkm = y[0] / 1000;
    ro = 1.22235514 - 0.11422776 * hkm + 0.00368329 * hkm * hkm - 0.0000404 * hkm * hkm * hkm;
    Xaer = 0.8 * ro * y[1] * y[1] * S / 2.0;

    double a = -0.8 * ro * S / 1.0 / y[2];
    double b = -beta / y[2];
    double c = a * y[1];

//    std::cout << y[0] << "," << y[1] << "," << y[2] << " ";
//    std::cout << a << "_" << b << "_" << c << std::endl;

    double P12 = 1.0 / b;
    double P22 = 2.0 * c / b / b + sqrt((4.0 * c * c + b * b - 2.0 * P12 * b * b) / b / b / b / b);
    double dx1 = y[0] - OPTIMAL[step_number][1];
    double dx2 = y[1] - OPTIMAL[step_number][2];
    double dU = -b * (P12 * dx1 + P22 * dx2);

    RESULT[rkStepsFirst + step_number][dimension + 4] = dU * T;// ! запись дополнительных переменных
    RESULT[rkStepsFirst + step_number][dimension + 5] = dx1;// ! запись дополнительных переменых
    RESULT[rkStepsFirst + step_number][dimension + 6] = dx2;// ! запись дополнительных переменых
    RESULT[rkStepsFirst + step_number][dimension + 7] = y[3];// ! запись дополнительных переменных
    //std::cout <<step_number << " " << -b << " " << " " << P12 << " " << P22 << " " << dx1 << " " << dx2 <<  " " << dU * T << std::endl;

    dy[3] = dU * T;
    if (dU * T > 0)
        y[3] = std::min(y[3] + dy[3], u_max);
    dy[0] = -y[1] * T;
    dy[1] = (-(beta * y[3] + Xaer) / y[2] + g) * T;
    dy[2] = - y[3] * T;

    if (fabs(1.0 / rungeKuttaSteps * (step_number + rkStepsFirst) - x0) < 1e-10) {
        tabs = x0 * T;
        //RESULT[rkStepsFirst + step_number][dimension + 1] = y[3];// ! запись дополнительных переменных
        RESULT[rkStepsFirst + step_number][dimension + 2] = tabs;// ! запись дополнительных переменных
        RESULT[rkStepsFirst + step_number][dimension + 3] = Xaer;// ! запись дополнительных переменных
        if (step_number >= rkStepsSecond && rkStepsSecond != 0) {
            exit_signal = 1;
        }
    }
}

void execRK (const vector<double> & x)
{
    resetGlobalVariables();

    int dimension = 4;
    rkStepsFirst  = 0; // N2
    rkStepsSecond = 0; // N3

    vector<double> y(4, 0);
    Matrix result(rungeKuttaSteps + 1, vector<double>(7, 0));
    int exit_signal;

    hkm = H / 1000;
    ro = 1.22235514- 0.11422776 * hkm + 0.00368329 * hkm * hkm - 0.0000404 * hkm * hkm * hkm;
    Xaer = 0.8 * ro * V * V * S / 2.0;
    U = 0;

    double x0 = 0.0;
    double dx = 1.0 / rungeKuttaSteps;
    P1 = x[0];
    P3 = x[1];
    //!P2=(P1*V+P3*U+U)/(g-(b*U+Xaer)/am)
    P2 = (P1 * V) / (g - Xaer / am);
    //! if (XX(3).Lt.5.d0) XX(3)=5.d0
    T = x[2];

    y[0] = H;
    y[1] = V;
    y[2] = am;

    int temp_step = 0;
    RK4(FUN2_2, dimension, rungeKuttaSteps, x0, dx, y, exit_signal, temp_step, result);
    rkStepsFirst = temp_step;

    for (int i = 0; i < rkStepsFirst + 1; i++) {
        for (int j = 0; j < dimension + 1; j++) {
            RESULT[i][j] = result[i][j];
        }
    }
}

void synthesis()
{
    std::cout << "synthesis()";
    int newton_steps = 0;
    vector<double> residuals(3, 1);
    vector<double> args = {P1, P3, T};
    newtonClassic(&getResidual, 3, args, residuals, newton_steps);
    OPTIMAL = RESULT;
    execRK(args);
    printMatrix(RESULT);
    //printMatrix(OPTIMAL);

}

int main()
{
    synthesis();
}
