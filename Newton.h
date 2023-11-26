#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

void printNumerically(double M);

void printAnalytically();

void computeJacobianNumerically(double x1, double x2, std::vector<std::vector <double>>& J, double M);

void computeJacobianAnalytically(double x1, double x2, std::vector<std::vector <double>>& J);

double f1(double x1, double x2);

double f2(double x1, double x2);
