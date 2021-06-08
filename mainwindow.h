#pragma once
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <Eigen/Dense>


//class CannyDetecter;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

#include "cannydetecter.h"

class CannyDetecter;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    CannyDetecter *Canny;
    QImage Gray;
    QImage Gau;
    QImage EImap;
    QImage NonMaxDeMap;
    QImage DoubleThresholdMap;
private slots:
    void on_LoadPic_clicked();

    void on_doGaussian_clicked();

    void on_EdgeIntensity_clicked();

    void on_NonMax_clicked();

    void on_DoubleT_clicked();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
