#pragma once
#ifndef CANNYDETECTER_H
#define CANNYDETECTER_H

#include <Eigen/Dense>
#include <QString>
#include <QDebug>
#include "mainwindow.h"
#include <QPalette>
#include <QColor>
#include <QProgressBar>
using namespace Eigen;

class QProgressBar;
class CannyDetecter
{
public:
    CannyDetecter();

    QProgressBar* parent;
    //QImage Raw;
    //QImage *Gray;
    MatrixXi GrayPic;

    MatrixXi GauPic;
    //QImage *Gau;

    MatrixXi EI;//Edge Intensity
    MatrixXi EI_x;
    MatrixXi EI_y;
    MatrixXi EdgeDir;//Edge Intensity rotation

    MatrixXi Edges;


    //QImage *EImap;

    int Height;
    int Width;
    void LoadGray(const QImage &);
    void Load(const MatrixXi &);
    void GaussFilter();
    void ApplySobel();
    bool hasAppliedNonMaxDe;
    bool hasAppliedDoubleThreshold;
    void ApplyNonMaxDe();
    void ApplyDoubleThreshold();

    void ApplyDoubleThreshold(float,float);
    void ApplySingleThreshold(float);
    static Matrix3i GauKernel;
    static int GauKernelSum;
    static Matrix3i SobelKernelX;
    static Matrix3i SobelKernelY;
    void getThreshold(int&,int&,int&);
    int OTSU(float&,float&);
    void ApplyOTSUThreshold();
private:
    char step;//0->Nothing,1->exported,2->Guassianed
};


#endif // CANNYDETECTER_H
