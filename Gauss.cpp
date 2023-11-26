#include "Gauss.h"

bool gauss(std::vector<std::vector<double>>& a, std::vector<double>& b, int n, std::vector<double>& x) {
    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            if (abs(a[i][k]) > abs(a[k][k])) {
                std::swap(a[i], a[k]);
                std::swap(b[i], b[k]);
            }
        }
        double A_Main = a[k][k];
        if (A_Main == 0) {
            return false;
        }
        for (int i = k; i < n; i++) {
            a[k][i] /= A_Main;
        }
        b[k] /= A_Main;
        for (int i = k + 1; i < n; i++) {
            double s = a[i][k];
            for (int j = k; j < n; j++) {
                a[i][j] -= s * a[k][j];
            }
            b[i] -= s * b[k];
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        x[k] = b[k];
        for (int i = n - 1; i > k; i--) {
            x[k] -= a[k][i] * x[i];
        }
    }
    return true;
}
