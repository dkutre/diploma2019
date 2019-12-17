#include "SubFunctions.h"

void printMatrix(const Matrix& m) {
    std::cout << "\n";
    for (size_t i = 0; i < m.size(); i++) {
        for (size_t j = 0; j < m[i].size(); j++) {
            std::cout.width(10);
            std::cout.precision(7);
            std::cout.fill(' ');
            std::cout << std::defaultfloat << m[i][j] << "; ";
        }
        std::cout << "\n";
    }
}
