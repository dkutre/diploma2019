#include "SubFunctions.h"

void printMatrix(const Matrix& m) {
    std::cout << "\n";
    for (size_t i = 0; i < m.size(); i++) {
        for (size_t j = 0; j < m[i].size(); j++) {
            std::cout.width(9);
            std::cout << m[i][j] << "; ";
        }
        std::cout << "\n";
    }
}
