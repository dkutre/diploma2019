#include "gauss.h"


/**
 * \file
 * \brief файл с реализацией метода LUP разложения
 */

void LUP(const Matrix &A, Matrix &LU, vector<int> &T, string& error) {
    //n - размерность исходной матрицы
    const int n = A.size();

    LU = A;


    //загружаем в матрицу T начальный вектор перестановок
    T = vector<int>(A.size(), 0);
    for (int i = 0; i < n; i++) {
        T[i] = i;
    }

    for (int i = 0; i < n; i++) {
        //поиск опорного элемента
        double pivotValue = 0;
        int pivot = -1;
        for (int row = i; row < n; row++) {
            if (fabs(LU[row][i]) > pivotValue) {
                pivotValue = fabs(LU[row][i]);
                pivot = row;
            }
        }
        if (pivotValue != 0) {
            //меняем местами i-ю строку и строку с опорным элементом
            std::swap(T[pivot], T[i]);
            LU[pivot].swap(LU[i]);
            for (int j = i + 1; j < n; j++) {
                LU[j][i] /= LU[i][i];
                for (int k = i + 1; k < n; k++) {
                    LU[j][k] -= LU[j][i] * LU[i][k];
                }
            }
        } else {
            error = "Не удалось выбрать ведущий элемент";
            return;
        }
    }
}


