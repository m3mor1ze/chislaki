#include "approximation.h"
#include <cmath>

#define M 2

double* approx(double* F, double* v, int n){
    double** sums = new double*[M];
    for(int i = 0; i < M; i++){
        sums[i] = new double[M];
        for(int j = 0; j < M; j++){
            double sum0 = 0;
            for(int k = 0; k < n; k++){
                sum0 += pow(F[k], i+j);
            }
            sums[i][j] = sum0;
        }
    }
    double* sums_y = new double[M];
    for(int i = 0; i < M; i++){
        double sumxy = 0;
        for(int j = 0; j < n; j++){
            sumxy += v[j]*pow(F[j], i);
        }
        sums_y[i] = sumxy;
    }
    double *solutions = new double[M];
    gauss(sums, sums_y, solutions, M);
    return solutions;
}

double* approx_interlayer(double* F, double* v, int n){
    double* tmpv = new double[n];
    double* tmpF = new double[n];
    for(int i = 0; i < n; i++){
        tmpv[i] = log(v[i]);
        tmpF[i] = log(F[i]);
    }
    double* solutions = approx(tmpF, tmpv, n);
    return solutions;
}
