#include "simpson.h"
#include <cmath>

double simpson(double a, double b, int m, double (*F)(double)) {
    int n = 2 * m;
    const double width = (b - a) / n;
    double simpson_integral = F(a);
    double sum1 = 0;
    for (int i = 1; i <= m; i++) {
        double xx = a + width * (2 * i - 1);
        sum1 += 4 * F(xx);
    }
    double sum2 = 0;
    for (int i = 1; i < m; i++) {
        double xx = a + width * 2 * i;
        sum2 += 2 * F(xx);
    }
    simpson_integral += sum1;
    simpson_integral += sum2;
    double xlast = a + width * n;
    simpson_integral += F(xlast);
    simpson_integral *= width / 3;
    return simpson_integral;
}

double simpson_wrapped(double a, double b, int m, double (*F)(double)) {
    double base_subintervals = m;
    double intervals = base_subintervals * 2;
    double simpson1 = simpson(a, b, base_subintervals, F);
    double simpson2 = simpson(a, b, intervals, F);
    while (!(std::abs(simpson1 - simpson2) <= 1e-6)) {
        base_subintervals = intervals;
        intervals *= 2;
        simpson1 = simpson(a, b, base_subintervals, F);
        simpson2 = simpson(a, b, intervals, F);
    }
    return simpson2;
}

double cubical_simpson(double a, double b, double c, double d, int n, int m, double (*F)(double, double)){
    const double width_x = (b-a)/(2*n);
    const double width_y = (d-c)/(2*m);
    double integral = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            double x = a + width_x*2*i;
            double y = c + width_y*2*j;
            double f1 = F(x, y);

            x = a + width_x*(2*i+1);
            y = c + width_y*(2*j);
            double f2 = 4*F(x, y);

            x = a + width_x*(2*i+2);
            y = c + width_y*(2*j);
            double f3 = F(x, y);

            x = a + width_x*(2*i);
            y = c + width_y*(j+1);
            double f4 = 4*F(x, y);

            x = a + width_x*(2*i+1);
            y = c + width_y*(2*j+1);
            double f5 = 16*F(x, y);

            x = a + width_x*(2*i+2);
            y = c + width_y*(2*j+1);
            double f6 = 4*F(x, y);

            x = a + width_x*(2*i);
            y = c + width_y*(2*j+2);
            double f7 = F(x, y);

            x = a + width_x*(2*i+1);
            y = c + width_y*(2*j+2);
            double f8 = 4*F(x, y);

            x = a + width_x*(2*i+2);
            y = c + width_y*(2*j+2);
            double f9 = F(x, y);

            double inner_result = f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9;
            integral += inner_result;
        }
    }
    integral *= width_x*width_y/9;
    return integral;
}
