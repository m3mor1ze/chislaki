#include "gauss.h"
#include <cmath>
#include <algorithm>
#include <cstring>
#include <iostream>

void gauss(double** a, double* b, double* solutions, int size){
    for (int i = 0; i < size; ++i) {
        for (int j = i + 1; j < size; ++j) {
            if (std::abs(a[j][i]) > std::abs(a[i][i])) {
                std::swap(a[i], a[j]);
                std::swap(b[i], b[j]);
            }
        }
        double amain = a[i][i];
        if (amain == 0) {
            std::cout << "Ошибка" << std::endl;
            exit(1);
        }
        for (int j = i; j < size; ++j) {
            a[i][j] /= amain;
        }
        b[i] /= amain;
        for (int j = i + 1; j < size; ++j) {
            double s = a[j][i];
            for (int k = i; k < size; ++k) {
                a[j][k] -= s * a[i][k];
            }
            b[j] -= s * b[i];
        }
    }
    for (int i = size - 1; i >= 0; --i) {
        solutions[i] = b[i];
        for (int j = size - 1; j > i; --j) {
            solutions[i] -= a[i][j] * solutions[j];
        }
    }
}
