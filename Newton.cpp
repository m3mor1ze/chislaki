#include "Newton.h"
#include "Gauss.h"

using namespace std;

double f1(double x1, double x2) {
    return cos(0.4 * x2 + x1 * x1) + x2 * x2 + x1 * x1 - 1.6;
}

double f2(double x1, double x2) {
    return 1.5 * x1 * x1 - x2 * x2 / 0.36 - 1;
}

void computeJacobianNumerically(double x1, double x2, vector<vector <double>>& J, double M) {
    double f1_val = f1(x1, x2);
    double f2_val = f2(x1, x2);
    double df1_dx1 = (f1(x1 + M, x2) - f1_val) / M;
    double df2_dx1 = (f2(x1 + M, x2) - f2_val) / M;
    double df1_dx2 = (f1(x1, M + x2) - f1_val) / M;
    double df2_dx2 = (f2(x1, M + x2) - f2_val) / M;
    J[0][0] = df1_dx1;
    J[0][1] = df1_dx2;
    J[1][0] = df2_dx1;
    J[1][1] = df2_dx2;
}

void computeJacobianAnalytically(double x1, double x2, vector<vector <double>>& J) {
    double df1_dx1 = -2 * x1 * sin(0.4 * x2 + x1 * x1) + 2 * x1;
    double df1_dx2 = -0.4 * sin(0.4 * x2 + x1 * x1) + 2 * x2;
    double df2_dx1 = 3 * x1;
    double df2_dx2 = -2 * x2 / 0.36;
    J[0][0] = df1_dx1;
    J[0][1] = df1_dx2;
    J[1][0] = df2_dx1;
    J[1][1] = df2_dx2;
}

void printNumerically(double M){
    const double e1 = 1e-9;
    const double e2 = 1e-9;
    const int NIT = 100;
    double x1 = 1.0;
    double x2 = -1.0;
    int k = 1;
    cout<<"##############################################################################"<<endl;
    cout<<"M = "<<M<<endl;
    cout<<"Iteration |" << setw(12) << "x1" << setw(6)<<"|"<< setw(12) << "x2" << setw(6)<<"|"<< setw(12) << "b1" << setw(6)<<"|" << setw(12)<< "b2" << endl;
    while (k <= NIT) {
        vector<vector<double>> J(2, vector<double>(2));
        computeJacobianNumerically(x1, x2, J, M);
        vector<double> F = {f1(x1, x2), f2(x1, x2)};
        vector<double> increments(2);
        if(!gauss(J, F, 2, increments))
            break;
        x1-=increments[0];
        x2-=increments[1];
        double b1 = max(abs(F[0]), abs(F[1]));
        double b2 = max(abs(increments[0]), abs(increments[1]));
        b2 = max(b2, abs(increments[0] / x1));
        b2 = max(b2, abs(increments[1] / x2)); 
        cout << k ;
        if(k<10)
            cout<<setw(10);
        else
            cout<<setw(9);
        cout<<"|"<<setw(12)<< x1 << setw(6)<<"|"<< setw(12) << x2 << setw(6)<<"|"<< setw(12) << b1 << setw(6)<<"|" << setw(12)<< b2 << endl;        
        if (b1 <= e1 && b2 <= e2) {
            cout << "Converged to the desired precision." << endl;            
            break;
        }
        if (k >= NIT) {
            cout << "Iteration limit reached. IER = 2" << endl;
            break;
        }
    k++;
    }
}

void printAnalytically(){
    const double e1 = 1e-9;
    const double e2 = 1e-9;
    const int NIT = 100;
    double x1 = 1.0;
    double x2 = -1.0;
    int k = 1;
    cout<<"##############################################################################"<<endl;
    cout<<"Jacobian Analytically:"<<endl;
    cout<<"Iteration |" << setw(12) << "x1" << setw(6)<<"|"<< setw(12) << "x2" << setw(6)<<"|"<< setw(12) << "b1" << setw(6)<<"|" << setw(12)<< "b2" << endl;
    while (k <= NIT) {
        vector<vector<double>> J(2, vector<double>(2));
        computeJacobianAnalytically(x1, x2, J);
        vector<double> F = {f1(x1, x2), f2(x1, x2)};
        vector<double> increments(2);
        if(!gauss(J, F, 2, increments))
            break;
        x1-=increments[0];
        x2-=increments[1];
        double b1 = max(abs(F[0]), abs(F[1]));
        double b2 = max(abs(increments[0]), abs(increments[1]));
        b2 = max(b2, abs(increments[0] / x1));
        b2 = max(b2, abs(increments[1] / x2)); 
        cout << k ;
        if(k<10)
            cout<<setw(10);
        else
            cout<<setw(9);
        cout<<"|"<<setw(12)<< x1 << setw(6)<<"|"<< setw(12) << x2 << setw(6)<<"|"<< setw(12) << b1 << setw(6)<<"|" << setw(12)<< b2 << endl;        
        if (b1 <= e1 && b2 <= e2) {
            cout << "Converged to the desired precision." << endl;            
            break;
        }
        if (k >= NIT) {
            cout << "Iteration limit reached. IER = 2" << endl;
            break;
        }
    k++;
    }
}
