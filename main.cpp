#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

// Функции, представляющие систему уравнений
double f1(double x1, double x2) {
    return cos(0.4 * x2 + x1 * x1) + x2 * x2 + x1 * x1 - 1.6;
}

double f2(double x1, double x2) {
    return 1.5 * x1 * x1 - x2 * x2 / 0.36 - 1;
}

// Производные функций по аргументам
double df1_dx1(double x1, double x2) {
    return -2 * x1 * sin(0.4 * x2 + x1 * x1) + 3 * x1;
}

double df1_dx2(double x1, double x2) {
    return -0.4 * sin(0.4 * x2 + x1 * x1) + 2 * x2;
}

double df2_dx1(double x1, double x2) {
    return 3 * x1;
}

double df2_dx2(double x1, double x2) {
    return -2 * x2 / 0.36;
}

int main() {
    double x1 = 1.0; // Начальное приближение для x1
    double x2 = -1.0; // Начальное приближение для x2
    const double e1 = 1e-9;
    const double e2 = 1e-9;
    const int NIT = 100;
    int k = 1;

    while (k <= NIT) {
        double F1 = f1(x1, x2);
        double F2 = f2(x1, x2);

        double J11 = df1_dx1(x1, x2);
        double J12 = df1_dx2(x1, x2);
        double J21 = df2_dx1(x1, x2);
        double J22 = df2_dx2(x1, x2);

        // Решение системы линейных уравнений
        double det = J11 * J22 - J12 * J21;
        double dx1 = (J22 * F1 - J12 * F2) / det;
        double dx2 = (-J21 * F1 + J11 * F2) / det;

        // Уточнение решения
        x1 -= dx1;
        x2 -= dx2;

        // Вычисление b1 и b2
        double b1 = max(abs(F1), abs(F2));
        double b2 = max(abs(dx1), abs(dx2));
        b2 = max(b2, abs(dx1 / x1));
        b2 = max(b2, abs(dx2 / x2));

        // Вывод текущих значений
        cout << "Iteration " << k << ": x1 = " << x1 << ", x2 = " << x2 << ", b1 = " << b1 << ", b2 = " << b2 << endl;

        // Проверка критерия завершения
        if (b1 <= e1 && b2 <= e2) {
            cout << "Converged to the desired precision." << endl;
            break;
        }

        // Проверка условия k >= NIT
        if (k >= NIT) {
            cout << "Iteration limit reached. IER = 2" << endl;
            break;
        }

        k++;
    }

    return 0;
}
