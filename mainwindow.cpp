#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QString>
#include <Eigen/Dense>
#include "cannydetecter.h"
#include <iostream>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    Canny =new CannyDetecter;
    Canny->parent=this;

    Canny->GauKernel<<57,94,57,
                                  94,155,94,
                                  57,94,57;
    Canny->GauKernelSum=Canny->GauKernel.array().sum();
    Canny->SobelKernelX<<-1,0,1,
                                        -2,0,2,
                                        -1,0,1;
    Canny->SobelKernelY<<-1,-2,-1,
                                        0,0,0,
                                         1,2,1;

    qDebug("构造函数运行完毕");
}

MainWindow::~MainWindow()
{
    delete Canny;
    delete ui;
}


void MainWindow::on_LoadPic_clicked()
{
    QString PicPath=QFileDialog::getOpenFileName(this,"Select an image","","Image Files(*.png *.jpg)");
    if (PicPath.isEmpty())
        return;

    QImage Raw(PicPath);
    ui->ShowRaw->setPixmap(QPixmap::fromImage(Raw));
    Gray=Raw.convertToFormat(QImage::Format_Grayscale8);
    ui->ShowGray->setPixmap(QPixmap::fromImage(Gray));
    //ui->ShowGray->setAutoFillBackground(true);
    Canny->step=1;
    Canny->LoadGray(Gray);
}

void MainWindow::on_doGaussian_clicked()
{
    if (Canny->step<1)return;
    Canny->GaussFilter();
    Gau=QImage(Gray);
    qDebug("成功初始化Gau图片");
    for(int r=0;r<Canny->Height;r++)
    {
        unsigned char*sL=Gau.scanLine(r);
        for(int c=0;c<Canny->Width;c++)
            sL[c]=Canny->GauPic(r+1,c+1);
    }
    ui->ShowGaussian->setPixmap(QPixmap::fromImage(Gau));
    //ui->ShowGaussian->setAutoFillBackground(true);
    Canny->step=2;
}



void MainWindow::on_EdgeIntensity_clicked()
{
    if(Canny->step<2)return;
    Canny->ApplySobel();
    EImap=QImage(Canny->Height,Canny->Width,QImage::Format_Grayscale8);
    int maxI=Canny->EI.array().maxCoeff();
    for(int r=0;r<Canny->Height;r++)
    {
        unsigned char*sL=EImap.scanLine(r);
        for(int c=0;c<Canny->Width;c++)
            sL[c]=255-(unsigned char)(Canny->EI(r+1,c+1)*255/maxI);
    }
    ui->ShowEdgeIntensity->setPixmap(QPixmap::fromImage(EImap));
    Canny->step=3;
}

void MainWindow::on_NonMax_clicked()
{
    if (Canny->step<3)return;
    //std::cout<<Canny->EdgeDir;
    Canny->ApplyNonMaxDe();

    NonMaxDeMap=QImage(Canny->Height,Canny->Width,QImage::Format_Grayscale8);
    int maxI=Canny->EI.array().maxCoeff();
    for(int r=0;r<Canny->Height;r++)
    {
        unsigned char*sL=NonMaxDeMap.scanLine(r);
        for(int c=0;c<Canny->Width;c++)
            sL[c]=(unsigned char)(Canny->EI(r+1,c+1)*255/maxI);
    }
    ui->ShowNonMax->setPixmap(QPixmap::fromImage(NonMaxDeMap));
    Canny->step=4;
}

void MainWindow::on_DoubleT_clicked()
{
    if(Canny->step<4)return;
    Canny->ApplyDoubleThreshold();
    NonMaxDeMap=QImage(Canny->Height,Canny->Width,QImage::Format_Grayscale8);
    //int maxI=Canny->EI.array().maxCoeff();
    qDebug("开始生成边缘识别结果");
    for(int r=0;r<Canny->Height;r++)
    {
        unsigned char*sL=NonMaxDeMap.scanLine(r);
        for(int c=0;c<Canny->Width;c++)
        {
            //qDebug()<<"r="<<r<<";  c="<<c<<";";
            sL[c]=(Canny->Edges(r+1,c+1)>0)?255:0;
        }

    }
    ui->ShowDoubeT->setPixmap(QPixmap::fromImage(NonMaxDeMap));
    Canny->step=5;
}
