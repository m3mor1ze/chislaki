#include "trapezoid.h"
#include <cmath>

double trapezoid(double a, double b, int n, double(*F)(double)) {
    const double width = (b - a) / n;
    double trapezoidal_integral = F(a);
    for (int i = 1; i < n; i++) {
        double xx = a + width * i;
        trapezoidal_integral += 2 * F(xx);
    }
    trapezoidal_integral += F(b);
    trapezoidal_integral *= width / 2;
    return trapezoidal_integral;
}

double trapezoid_wrapped(double a, double b, int n, double(*F)(double)) {
    double base_n = n;
    double n_2 = base_n * 2;
    double eps = 1e-6;
    double trapezoid1 = trapezoid(a, b, base_n, F);
    double trapezoid2 = trapezoid(a, b, n_2, F);
    while (!(std::abs(trapezoid1 - trapezoid2) <= 3 * eps)) {
        base_n = n_2;
        n_2 *= 2;
        trapezoid1 = trapezoid(a, b, base_n, F);
        trapezoid2 = trapezoid(a, b, n_2, F);
    }
    return trapezoid2;
}

