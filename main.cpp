#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

void gauss(vector<vector<double>>& a, vector<double>& b, int n, vector<double>& x) {
    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {//свап строк
            if (abs(a[i][k]) > abs(a[k][k])) {
                swap(a[i], a[k]);
                swap(b[i], b[k]);
            }
        }
        //выбор ведущего элемента
        double A_Main = a[k][k];
        if (A_Main == 0) {
            cout << "Ошибка" << endl;
            return;
        }
        //делим на ведущий элемент строку
        for (int i = k; i < n; i++) {
            a[k][i] /= A_Main;
        }
        b[k] /= A_Main;
        //исключение переменной
        for (int i = k + 1; i < n; i++) {
            double s = a[i][k];
            for (int j = k; j < n; j++) {
                a[i][j] -= s * a[k][j];
            }
            b[i] -= s * b[k];
        }
    }
    //обратный ход
    for (int k = n - 1; k >= 0; k--) {
        x[k] = b[k];
        for (int i = n - 1; i > k; i--) {
            x[k] -= a[k][i] * x[i];
        }
    }
}


vector<double> vectornev(vector<vector<double>>& a, vector<double>& b, vector<double>& x, int n) {
    vector<double> f(n);
    for (int i = 0; i < n; i++) {
        f[i] = -b[i];
        for (int j = 0; j < n; j++) {
            f[i] += a[i][j] * x[j];
        }
    }
    return f;
}


double Norma(vector<double>& f, int n) {
    double norma = abs(f[0]);
    for (int i = 0; i < n; i++) {
        norma = max(f[i], norma);
    }
    return norma;
}

int main() {
    int n = 3;
    vector<vector<double>> a = {{2.75, 1.78, 1.11}, {3.28, 0.71, 1.15}, {1.15, 2.7, 3.58}};
    vector<double> b = {15.71, 43.78, 37.11};
    vector<double> x(n);
    gauss(a, b, n, x);

    for (int i = 0; i < n; i++) {
        cout << x[i] << " ";
    }
    cout << endl;

    vector<double> f = vectornev(a, b, x, n);

    cout << "Норма: " << Norma(f, n) << endl;

    return 0;
}

