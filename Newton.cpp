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
    J[0][0] = (f1(x1 + x1 * M, x2) - f1(x1, x2)) / (x1 * M);
    J[0][1] = (f1(x1, x2 + x2 * M) - f1(x1, x2)) / (x2 * M);
    J[1][0] = (f2(x1 + x1 * M, x2) - f2(x1, x2)) / (x1 * M);
    J[1][1] = (f2(x1, x2 + x2 * M) - f2(x1, x2)) / (x2 * M);
}

void computeJacobianAnalytically(double x1, double x2, vector<vector <double>>& J) {
    J[0][0] = -2 * x1 * sin(0.4 * x2 + x1 * x1) + 2 * x1;
    J[0][1] = -0.4 * sin(0.4 * x2 + x1 * x1) + 2 * x2;
    J[1][0] = 3 * x1;
    J[1][1] = -2 * x2 / 0.36;
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
        x1 -= increments[0];
        x2 -= increments[1];

        double b1 = max(abs(F[0]), abs(F[1]));
        double b2;
        for (int i = 0; i < 2; i++) {
            if (abs(F[i]) < 1)
                b2 = max(b2, abs(F[i] - (F[i] - increments[i])));
            else
                b2 = max(b2, abs((F[i] - (F[i] - increments[i])) / F[i]));
        }
        k++;
        cout << k ;
        if(k<10)
            cout<<setw(10);
        else
            cout<<setw(9);
        cout<<"|"<<setw(12)<< x1 << setw(6)<<"|"<< setw(12) << x2 << setw(6)<<"|"<< setw(12) << b1 << setw(6)<<"|" << setw(12)<< b2 << endl;        
        if (b1 <= e1 && b2 <= e2) {
            cout << "Converged to the desired precision." << endl;     
            cout << "f1 = " << f1(x1, x2)<< endl << "f2 = " << f2(x1, x2)<<endl;       
            break;
        }
        if (k >= NIT) {
            cout << "f1 = " << f1(x1, x2)<< endl << "f2 = " << f2(x1, x2)<<endl;
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
        x1 -= increments[0];
        x2 -= increments[1];
        double b1 = max(abs(F[0]), abs(F[1]));
        double b2;
        for (int i = 0; i < 2; i++) {
            if (abs(F[i]) < 1)
                b2 = max(b2, abs(F[i] - (F[i] - increments[i])));
            else
                b2 = max(b2, abs((F[i] - (F[i] - increments[i])) / F[i]));
        } 
        cout << k ;
        if(k<10)
            cout<<setw(10);
        else
            cout<<setw(9);
        cout<<"|"<<setw(12)<< x1 << setw(6)<<"|"<< setw(12) << x2 << setw(6)<<"|"<< setw(12) << b1 << setw(6)<<"|" << setw(12)<< b2 << endl; 
        if (b1 <= e1 && b2 <= e2) {
            cout << "Converged to the desired precision." << endl;            
            cout << "f1 = " << f1(x1, x2)<< endl << "f2 = " << f2(x1, x2)<<endl;
            break;
        }
        if (k >= NIT) {
            cout << "f1 = " << f1(x1, x2)<< endl << "f2 = " << f2(x1, x2)<<endl;
            cout << "Iteration limit reached. IER = 2" << endl;
            break;
        }
    k++;
    }
}
