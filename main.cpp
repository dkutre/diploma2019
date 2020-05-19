#include "gauss.h"
#include "SubFunctions.h"
#include "iterSystemSolver.h"
#include "RK4.h"
#include "Classic_Newton.h"
#include <cmath>
#include <algorithm>

#include <random>
#include <cmath>

double H = 2000.0;
double V = 100.0;
double am = 20000.0;
double g = 9.81;
double S = 10.75;
double P1 =  -0.15;

double hkm = H / 1000;
double ro = 1.22235514 - 0.11422776 * hkm + 0.00368329 * hkm * hkm - 0.0000404 * hkm * hkm * hkm;
double cx = 0.8;
double Xaer = cx* ro * V * V * S / 2.0;

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

Matrix RESULT(rungeKuttaSteps, vector<double>(12, 0));
Matrix OPTIMAL(rungeKuttaSteps, vector<double>(12, 0));
Matrix INTERMEDIATE_RESULT(rungeKuttaSteps + 1, vector<double>(7, 0));

double get_ro(double h_km)
{
    static std::random_device rd{};
    static std::mt19937 gen{rd()};
    static std::normal_distribution<> du{1, 0.01};
    double random = du(gen);
    return random * (1.22235514 - 0.11422776 * hkm + 0.00368329 * hkm * hkm - 0.0000404 * hkm * hkm * hkm);
}

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
    Xaer = cx * ro * y[1] * y[1] * S / 2.0;


    dy[0] = -y[1] * T ;
    dy[1] = (-(beta * U + Xaer) / y[2] + g) * T;
    dy[2] = - U * T;
    double drodH = (-0.11422776 + 2.0 * 0.00368329 * hkm - 3.0 * 0.0000404 * hkm * hkm) * 0.001;
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
    Xaer = cx * ro * V * V * S / 2.0;


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
    //    correlation[i][1] = findSwitchMoment(); // момент переключения
        correlation[i][2] = RESULT[0][3] - RESULT[1000][3]; // затраченное топливо
        if (i % 100 == 0) {
            std::cout << "ke" << std::endl;
        }
        //printMatrix(RESULT);
    }

    printMatrix(correlation);
}

std::size_t findSwitchIndex(const Matrix & matrix)
{
    for (size_t i = 0; i < matrix.size(); i++)
        if (matrix[i][9] > 0) // тяга стала больше 0
            return i - 1;
    throw std::runtime_error("нет момента переключения");
}


void free_fall (int dimension, double x0, vector<double> y, vector<double>& dy, int& exit_signal, unsigned int step_number)
{
    hkm = y[0] / 1000;
    ro = get_ro(hkm);
    Xaer = cx * ro * y[1] * y[1] * S / 2.0;
    /*
    double a = -0.8 * ro * S / 1.0 / y[2];
    double b = -beta / y[2];
    double c = a * y[1];

    double P12 = 1.0 / b;
    double P22 = 2.0 * c / b / b + sqrt((4.0 * c * c + b * b - 2.0 * P12 * b * b) / b / b / b / b);
    double dx1 = INTERMEDIATE_RESULT[step_number][1] - OPTIMAL[step_number][1];
    double dx2 = INTERMEDIATE_RESULT[step_number][2] - OPTIMAL[step_number][2];
    double dU = -b * (P12 * dx1 + P22 * dx2);
    U = OPTIMAL[step_number + 1][9] + dU;
    U = std::max(0.0, std::min(U, 305.0));
    т.к методе Ньютора когда достигаем момента переключения делаем шаг назад и заново просчитываем шаг уже с учетом U_switch

    RESULT[step_number][dimension + 4] = dU;// ! запись дополнительных переменных
    RESULT[step_number][dimension + 5] = dx1;// ! запись дополнительных переменых
    RESULT[step_number][dimension + 6] = dx2;// ! запись дополнительных переменых
    RESULT[step_number][dimension + 7] = U;//y[3];// ! запись дополнительных переменных
    */
    //std::cout <<step_number << " " << -b << " " << " " << P12 << " " << P22 << " " << dx1 << " " << dx2 <<  " " << dU * T << std::endl;

    U = 0;
    dy[0] = -y[1] * T;
    dy[1] = (-(beta * U + Xaer) / y[2] + g) * T;
    dy[2] = - U * T;

    static auto local_step_number = step_number;
    if (local_step_number == step_number)
        return;
    //std::cout << step_number;
    local_step_number = step_number;
    auto switch_index = findSwitchIndex(OPTIMAL);
    if (y[1] > OPTIMAL[switch_index][2] || y[0] < OPTIMAL[switch_index][1])
    {
        //std::cout << "\nexit_signal = 1, step" << step_number << std::endl;
        //std::cout << y[1] << " " << OPTIMAL[switch_index][2] << " " << y[0] << " " << OPTIMAL[switch_index][1];
        exit_signal = 1;
    }
}

double linear_interpolation(double x, double x0, double y0, double x1, double y1)
{
    if (x1 == x0)
        throw std::runtime_error("Can't do linear interpolation, x0 == x1");
    //std::cout << "x: " << x << ", x0: " << x0 << ", y0: " << y0 << ", x1:" << x1 << ", y1:" << y1 << std::endl;
    const double y = y0 + ((y1 - y0) / (x1 - x0)) * (x - x0);
    //std::cout << "y: " << y << std::endl;
    return y;
}

double getVByH(const Matrix& stats, double H)
{
    for (size_t i = 1; i < stats.size(); i++) // должна быть точка слева
    {
        if (stats[i][1] < H)
        {
           const std::size_t left = i - 1;
           const std::size_t right = i;
           double res = linear_interpolation(H, stats[left][1], stats[left][2], stats[right][1], stats[right][2]);
           return res;
        }
    }
    return 0;
}

void engine_on (int dimension, double x0, vector<double> y, vector<double>& dy, int& exit_signal, unsigned int step_number)
{
    hkm = y[0] / 1000;
    ro = get_ro(hkm);
    Xaer = cx * ro * y[1] * y[1] * S / 2.0;
    /*
    double a = -0.8 * ro * S / 1.0 / y[2];
    double b = -beta / y[2];
    double c = a * y[1];

    double P12 = 1.0 / b;
    double P22 = 2.0 * c / b / b + sqrt((4.0 * c * c + b * b - 2.0 * P12 * b * b) / b / b / b / b);
    double dx1 = INTERMEDIATE_RESULT[step_number][1] - OPTIMAL[step_number][1];
    double dx2 = INTERMEDIATE_RESULT[step_number][2] - OPTIMAL[step_number][2];
    double dU = -b * (P12 * dx1 + P22 * dx2);
    U = OPTIMAL[step_number + 1][9] + dU;
    U = std::max(0.0, std::min(U, 305.0));
    т.к методе Ньютора когда достигаем момента переключения делаем шаг назад и заново просчитываем шаг уже с учетом U_switch

//   RESULT[step_number][dimension + 4] = dU;// ! запись дополнительных переменных
//    RESULT[step_number][dimension + 5] = dx1;// ! запись дополнительных переменых
//    RESULT[step_number][dimension + 6] = dx2;// ! запись дополнительных переменых
//    RESULT[step_number][dimension + 7] = U;//y[3];// ! запись дополнительных переменных
//    */
    //std::cout <<step_number << " " << -b << " " << " " << P12 << " " << P22 << " " << dx1 << " " << dx2 <<  " " << dU * T << std::endl;

    static auto local_step_number = -1;
    static double akp = 25;
    static double delta_v = 0;
    static double delta_u = akp * delta_v;
    //std::cout << "y[1]: " << y[1] << ", step: " << step_number << std::endl;
    if (local_step_number != step_number)
    {
        delta_v = y[1] - getVByH(OPTIMAL, y[0]);
        delta_u = akp * delta_v;
    //    std::cout << "speed: " << y[1] << " " << getVByH(OPTIMAL, y[0]) << "\n";
    //    std::cout << "delta_v: " << delta_v << std::endl;
    //    std::cout << "delta_u: " << delta_u << std::endl;
        local_step_number = step_number;
    }

    U = u_max + delta_u;
    dy[0] = -y[1] * T;
    dy[1] = (-(beta * U + Xaer) / y[2] + g) * T;
    dy[2] = - U * T;

    if (y[1] < OPTIMAL.back()[2] || y[0] < OPTIMAL.back()[1])
    {
    //    std::cout << "\nexit_signal = 1, step" << step_number << std::endl;
        exit_signal = 1;
    }
    //std::cout << "=======" << std::endl;
}

void end_free_fall (int dimension, double x0, vector<double> y, vector<double>& dy, int& exit_signal, unsigned int step_number)
{
    hkm = y[0] / 1000;
    ro = get_ro(hkm);
    Xaer = cx * ro * y[1] * y[1] * S / 2.0;
    /*
    double a = -0.8 * ro * S / 1.0 / y[2];
    double b = -beta / y[2];
    double c = a * y[1];

    double P12 = 1.0 / b;
    double P22 = 2.0 * c / b / b + sqrt((4.0 * c * c + b * b - 2.0 * P12 * b * b) / b / b / b / b);
    double dx1 = INTERMEDIATE_RESULT[step_number][1] - OPTIMAL[step_number][1];
    double dx2 = INTERMEDIATE_RESULT[step_number][2] - OPTIMAL[step_number][2];
    double dU = -b * (P12 * dx1 + P22 * dx2);
    U = OPTIMAL[step_number + 1][9] + dU;
    U = std::max(0.0, std::min(U, 305.0));
    т.к методе Ньютора когда достигаем момента переключения делаем шаг назад и заново просчитываем шаг уже с учетом U_switch

    RESULT[step_number][dimension + 4] = dU;// ! запись дополнительных переменных
    RESULT[step_number][dimension + 5] = dx1;// ! запись дополнительных переменых
    RESULT[step_number][dimension + 6] = dx2;// ! запись дополнительных переменых
    RESULT[step_number][dimension + 7] = U;//y[3];// ! запись дополнительных переменных
    */
    //std::cout <<step_number << " " << -b << " " << " " << P12 << " " << P22 << " " << dx1 << " " << dx2 <<  " " << dU * T << std::endl;

    U = 0;
    dy[0] = -y[1] * T;
    dy[1] = (-(beta * U + Xaer) / y[2] + g) * T;
    dy[2] = - U * T;

    static auto local_step_number = step_number;
    if (local_step_number == step_number)
        return;
    //std::cout << step_number;
    local_step_number = step_number;
    if (y[0] < 0)
    {
        std::cout << "\nexit_signal = 1, step" << step_number << std::endl;
        exit_signal = 1;
    }
}

void append(const Matrix& source, const size_t source_size, Matrix& dest)
{
    for (size_t i = 0; i <= source_size && i < source.size(); i++)
        dest.push_back(source[i]);
}

void appendWithoutFirstRow(const Matrix& source, const size_t source_size, Matrix& dest)
{
    for (size_t i = 1; i <= source_size && i < source.size(); i++)
        dest.push_back(source[i]);
}

// для последней строчки считаем значения при 0 высоте
void calcLastRow(Matrix& matrix)
{
    //std::cout << "calcLastRow()" << std::endl;
    if (matrix.size() < 2)
        throw std::runtime_error("Should be more then one row");
    const auto second_from_end_row = matrix[matrix.size() - 2]; // вторая с конца строка
    //std::cout << "second_from_end_row";
    printContainer(second_from_end_row);
    const auto last_row = matrix.back();
    //std::cout << "last_row";
    printContainer(last_row);

    auto new_last_row = last_row;
    for (std::size_t i = 2; i < last_row.size(); i++)
    {
       if (second_from_end_row[i] != last_row[i])
           new_last_row[i] = linear_interpolation(
               0,
               second_from_end_row[1],
               second_from_end_row[i],
               last_row[1],
               last_row[i]
               );
    }
    new_last_row[1] = 0;
    //std::cout << "new_last_row";
    printContainer(new_last_row);

    matrix.back() = new_last_row;
}

void execRK (const vector<double> & x, const vector<double> init_y, Matrix & all_result)
{
    //std::cout << "execRK()" << std::endl;
    resetGlobalVariables();

    int dimension = 3;
    rkStepsFirst  = 0; // N2
    rkStepsSecond = 0; // N3

    vector<double> y = init_y;
    int exit_signal;

    hkm = y[0] / 1000;
    ro = 1.32235514 - 0.11422776 * hkm + 0.00368329 * hkm * hkm - 0.0000404 * hkm * hkm * hkm;
    Xaer = cx * ro * y[1] * y[1] * S / 2.0;
    U = 0;

    double x0 = 0.0;
    double dx = 1.0 / rungeKuttaSteps;
    P1 = x[0];
    P3 = x[1];
    //!P2=(P1*V+P3*U+U)/(g-(b*U+Xaer)/am)
    P2 = (P1 * V) / (g - Xaer / am);
    //! if (XX(3).Lt.5.d0) XX(3)=5.d0
    T = x[2];

    //y[0] = H;
    //y[1] = V;
    //y[2] = am;

    int temp_step = 0;
    std::cout << "RK4(), free_fall" << std::endl;
    RK4(free_fall, dimension, rungeKuttaSteps, x0, dx, y, exit_signal, temp_step, INTERMEDIATE_RESULT);
    std::cout << "/RK4()" << std::endl;
    std::cout << temp_step;
    append(INTERMEDIATE_RESULT, temp_step, all_result);

    INTERMEDIATE_RESULT = Matrix(rungeKuttaSteps + 1, vector<double>(7, 0));
    std::cout << "RK4(), engine_on" << std::endl;
    y[0] = all_result.back()[1];
    y[1] = all_result.back()[2];
    y[2] = all_result.back()[3];
    RK4(engine_on, dimension, rungeKuttaSteps, x0, dx, y, exit_signal, temp_step, INTERMEDIATE_RESULT);
    std::cout << "/RK4()" << std::endl;
    appendWithoutFirstRow(INTERMEDIATE_RESULT, temp_step, all_result);

    // недолетели, т.е отрицательная скорость на ненулевой высоте (летим вверх)
    if (all_result.back()[1] > OPTIMAL.back()[1])
    {
        INTERMEDIATE_RESULT = Matrix(rungeKuttaSteps + 1, vector<double>(7, 0));
        std::cout << "RK4(), free fall after engine off" << std::endl;
        y[0] = all_result.back()[1];
        y[1] = all_result.back()[2];
        y[2] = all_result.back()[3];
        RK4(end_free_fall, dimension, rungeKuttaSteps, x0, dx, y, exit_signal, temp_step, INTERMEDIATE_RESULT);
        std::cout << "/RK4()" << std::endl;
        appendWithoutFirstRow(INTERMEDIATE_RESULT, temp_step, all_result);
    }

    calcLastRow(all_result); // высчитать значения при нулевой высоте

    RESULT = all_result;
    std::cout << "/execRK()" << std::endl;
}

Matrix h_distribution()
{
    vector<double> args = {P1, P3, T};

    double h_variance = H * 0.03;
    std::random_device h_rd{};
    std::mt19937 h_gen{h_rd()};
    std::normal_distribution<> dh{H, h_variance};

    Matrix results_with_random(0, vector<double>(12, 0));
    for (std::size_t i = 0; i < 1000; i++)
    {
        Matrix result(0, vector<double>(12, 0));
        H = dh(h_gen);
        vector<double> init_y = {H, V, am};
        execRK(args, init_y, result);
        if (result.back()[2] < 0)
            continue;
        results_with_random.push_back(result.back());
        results_with_random.back().back() = H;
    }
    return results_with_random;
}

Matrix v_distribution()
{
    vector<double> args = {P1, P3, T};

    double variance = V * 0.05;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{V, variance};

    Matrix results_with_random(0, vector<double>(12, 0));
    for (std::size_t i = 0; i < 10000; i++)
    {
        Matrix result(0, vector<double>(12, 0));
        double local_V = d(gen);
        vector<double> init_y = {H, local_V, am};
        execRK(args, init_y, result);
        if (result.back()[2] < 0)
            continue;
        results_with_random.push_back(result.back());
        results_with_random.back().back() = local_V;
    }
    return results_with_random;
}

Matrix m_distribution()
{
    vector<double> args = {P1, P3, T};

    double variance = am * 0.01;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{am, variance};

    Matrix results_with_random(0, vector<double>(12, 0));
    for (std::size_t i = 0; i < 1000; i++)
    {
        Matrix result(0, vector<double>(12, 0));
        double local_am = d(gen);
        vector<double> init_y = {H, V, local_am};
        execRK(args, init_y, result);
        if (result.back()[2] < 0)
            continue;
        results_with_random.push_back(result.back());
        results_with_random.back().back() = local_am;
    }
    return results_with_random;
}

Matrix cx_distribution()
{
    vector<double> args = {P1, P3, T};

    double variance = cx * 0.03;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{cx, variance};

    Matrix results_with_random(0, vector<double>(12, 0));
    for (std::size_t i = 0; i < 10000; i++)
    {
        Matrix result(0, vector<double>(12, 0));
        cx = d(gen);
        vector<double> init_y = {H, V, am};
        execRK(args, init_y, result);
        if (result.back()[2] < 0)
            continue;
        results_with_random.push_back(result.back());
        results_with_random.back().back() = cx;
    }
    return results_with_random;
}

Matrix all_distribution()
{
    vector<double> args = {P1, P3, T};

    double h_variance = H * 0.0075;
    std::random_device h_rd{};
    std::mt19937 h_gen{h_rd()};
    std::normal_distribution<> dh{H, h_variance};

    double v_variance = V * 0.0075;
    std::random_device v_rd{};
    std::mt19937 v_gen{v_rd()};
    std::normal_distribution<> dv{V, v_variance};

    double cx_variance = cx * 0.0075;
    std::random_device cx_rd{};
    std::mt19937 cx_gen{cx_rd()};
    std::normal_distribution<> dcx{cx, cx_variance};

    double m_variance = am * 0.0075;
    std::random_device m_rd{};
    std::mt19937 m_gen{m_rd()};
    std::normal_distribution<> dm{am, m_variance};

    double u_variance = u_max * 0.0075;
    std::random_device u_rd{};
    std::mt19937 u_gen{u_rd()};
    std::normal_distribution<> du{u_max, u_variance};

    Matrix results_with_random(0, vector<double>(12, 0));
    for (std::size_t i = 0; i < 10000; i++)
    {
        Matrix result(0, vector<double>(12, 0));

        cx = dcx(cx_gen);
        double local_am = dm(m_gen);
        double local_v = dv(v_gen);
        double local_h = dh(h_gen);
        u_max = du(u_gen);

        vector<double> init_y = {local_h, local_v, local_am};
        execRK(args, init_y, result);

        if (result.back()[2] < 0)
            continue;
        results_with_random.push_back(result.back());
        results_with_random.back().back() = cx;
    }
    return results_with_random;
}

void synthesis()
{
    std::cout << "synthesis()" << std::endl;
    int newton_steps = 0;
    vector<double> residuals(3, 1);
    vector<double> args = {P1, P3, T};
    std::cout << "newtonClassic()" << std::endl;
    newtonClassic(&getResidual, 3, args, residuals, newton_steps);
    std::cout << "/newtonClassic()" << std::endl;
    OPTIMAL = RESULT;


    //Matrix result(0, vector<double>(12, 0));
    //vector<double> init_y = {H, 1998.422, am};
    //execRK(args, init_y, result);
    //printMatrix(result);
    //const auto results_to_print = h_distribution();
    //const auto results_to_print = v_distribution();
    //const auto results_to_print = m_distribution();
    //const auto results_to_print = cx_distribution();
    const auto results_to_print = all_distribution();
    printMatrix(results_to_print);
    //printMatrix(OPTIMAL);
    std::cout << "/synthesis()" << std::endl;
}

int main()
{
    synthesis();
}
