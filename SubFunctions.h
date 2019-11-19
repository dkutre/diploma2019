/*!
 * \file
 * \brief Файл с общими фунциями для остальных модулей
 *
*/

#ifndef GAUSS_SUBFUNCTIONS_H
#define GAUSS_SUBFUNCTIONS_H

#include <vector>
#include <iostream>
#include <string>

using std::vector;

typedef vector<vector<double>> Matrix;

void printMatrix(const Matrix& m);


template <class Container>
void printContainer(const Container& container) {
    std::cout << "\n";
    for (typename Container::const_iterator it = container.cbegin(); it != container.cend(); it++) {
        std::cout << *it << " ";
    }
    std::cout << "\n";
};

#endif //GAUSS_SUBFUNCTIONS_H
