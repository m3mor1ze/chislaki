#include "trapezoid.h"
#include "simpson.h"
#include <iostream>
#include <cmath>
using namespace std;

double f1(double x){

    return (pow(x, 2) + 1)/(pow(x, 3) + 1);
}

double f2(double x, double y){
    return 4 - pow(x, 2) - pow(y, 2);
}

int main(){
    double A = 3.0, B = 4.254, N = 50, A2d = -1.0, B2d = 1.0, C2d = -1.0, D2d = 1.0, N2d = 100, M2d = 100;
    cout << "F(x) = (pow(x, 2) + 1)/(pow(x, 3) + 1)" << endl;
    cout << "TRAPEZOID: " << trapezoid_wrapped(A, B, N, f1) << endl;
    cout << "SIMPSON: " << simpson_wrapped(A, B, N, f1) << endl;
    cout << "CUBICAL SIMPSON: " << cubical_simpson(A2d, B2d, C2d, D2d, N2d, M2d, f2) << endl;
    return 0;
}
