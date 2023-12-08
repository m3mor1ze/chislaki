#include "approximation.h"
#include "qcustomplot.h"
#include <Qt>
#include <QPen>
#include <QBrush>
#include <QVector>
#include <QPushButton>
#include <QApplication>
#include <iostream>
#include <fstream>
#include <cmath>

#define WINDOW_W 800
#define WINDOW_H 600

using namespace std;

int main(int argc, char **argv){
    ifstream datafile = ifstream("../src/data.txt");
    int len;
    datafile >> len;
    double* F = new double[len];
    for(int i = 0; i < len; i++){
        datafile >> F[i];
    }
    double* v = new double[len];
    for(int i = 0; i < len; i++){
        datafile >> v[i];
    }

    double* odds = approx_interlayer(F, v, len);
    double c = exp(odds[0]);
    double e = -1.0 /odds[1];
    cout << "c = "<< c << endl<< "e = " << e << endl;
    cout<<"v = c * F^(-1/e)"<<endl;
    /*    
    for(int i = 0; i < len; i++){
        cout << v[i] << " = " << c * pow(F[i], -1.0/e)<<endl;
    }
    */

    QApplication app(argc, argv);    
    QCustomPlot* customPlot = new QCustomPlot();
    customPlot->addGraph();
    customPlot->graph(0)->setPen(QPen(Qt::red));
    customPlot->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20)));
    QVector<double> x(13), y(13);
    x[0] = 0;
    y[0] = 0;
    for(int i = 1; i < 13; i++){
        x[i] = F[i];    
        y[i] = v[i];    
    }
    customPlot->xAxis2->setVisible(true);
    customPlot->xAxis2->setTickLabels(false);
    customPlot->yAxis2->setVisible(true);
    customPlot->yAxis2->setTickLabels(false);
    customPlot->graph(0)->setData(x, y);
    customPlot->graph(0)->rescaleAxes();
    customPlot->resize(WINDOW_W, WINDOW_H);
    customPlot->show();
    delete[] F;
    delete[] v;
    delete[] odds;
    return app.exec();
}
