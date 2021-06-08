#include "cannydetecter.h"
#include <cmath>
#include <iostream>
Matrix3i CannyDetecter::GauKernel=Matrix3i::Zero(3,3);
int CannyDetecter::GauKernelSum=0;
Matrix3i CannyDetecter::SobelKernelX=Matrix3i::Zero(3,3);
Matrix3i CannyDetecter::SobelKernelY=Matrix3i::Zero(3,3);

CannyDetecter::CannyDetecter()
{
//Gray=QImage();
//Raw=QImage();
step=0;
parent=NULL;
}

void CannyDetecter::LoadGray(const QImage& Src)
{
    Height=Src.height();
    Width=Src.width();
    GrayPic.setZero(Height+2,Width+2);
    for(int r=0;r<Height;r++)
    {
        unsigned char*Current=(unsigned char*)Src.scanLine(r);
        for(int c=0;c<Width;c++)
            GrayPic(r+1,c+1)=Current[c];
    }

    GrayPic.block(0,1,1,Width)=GrayPic.block(1,1,1,Width);
    GrayPic.block(Height+1,1,1,Width)=GrayPic.block(Height,1,1,Width);
    GrayPic.block(1,0,Height,1)=GrayPic.block(1,1,Height,1);
    GrayPic.block(1,Width+1,Height,1)=GrayPic.block(1,Width,Height,1);
    GrayPic(0,0)=GrayPic(1,1);
    GrayPic(Height+1,0)=GrayPic(Height,1);
    GrayPic(0,Width+1)=GrayPic(1,Width);
    GrayPic(Height+1,Width+1)=GrayPic(Height,Width);
    qDebug("成功加载灰度图为Eigen矩阵并处理边界");
}

void CannyDetecter::GaussFilter()
{
    if(step<1)return;
    GauPic.setZero(Height+2,Width+2);
    for(int r=1;r<=Height;r++)
        for(int c=1;c<=Width;c++)
        {
            auto Toki=GrayPic.block(r-1,c-1,3,3).array()*GauKernel.array();
            GauPic(r,c)=Toki.sum()/GauKernelSum;
        }
    GauPic.block(0,1,1,Width)=GauPic.block(1,1,1,Width);
    GauPic.block(Height+1,1,1,Width)=GauPic.block(Height,1,1,Width);
    GauPic.block(1,0,Height,1)=GauPic.block(1,1,Height,1);
    GauPic.block(1,Width+1,Height,1)=GauPic.block(1,Width,Height,1);
    GauPic(0,0)=GauPic(1,1);
    GauPic(Height+1,0)=GauPic(Height,1);
    GauPic(0,Width+1)=GauPic(1,Width);
    GauPic(Height+1,Width+1)=GauPic(Height,Width);
    qDebug("高斯滤波完成");
}

void CannyDetecter::ApplySobel()
{
    if(step<2)return;
    EI.setZero(Height+2,Width+2);
    EdgeDir.setZero(Height,Width);
    EI_x.setZero(Height+2,Width+2);
    EI_y.setZero(Height+2,Width+2);
    float Gx;
    float Gy;
    for(int r=1;r<=Height;r++)
        for(int c=1;c<=Width;c++)
        {
            Gx=((GauPic.block(r-1,c-1,3,3).array()*SobelKernelX.array()).sum());
            Gy=((GauPic.block(r-1,c-1,3,3).array()*SobelKernelY.array()).sum());
            EI_x(r,c)=Gx;
            EI_y(r,c)=Gy;
            EI(r,c)=round(sqrt(float(Gx*Gx+Gy*Gy)));
            EdgeDir(r-1,c-1)=round((atan2(Gx,Gy)/M_PI)*180.0f);
        }


    EdgeDir=(EdgeDir.array()>=0).select(EdgeDir,EdgeDir.array()+180);

    EI.block(0,1,1,Width)=EI.block(1,1,1,Width);
    EI.block(Height+1,1,1,Width)=EI.block(Height,1,1,Width);
    EI.block(1,0,Height,1)=EI.block(1,1,Height,1);
    EI.block(1,Width+1,Height,1)=EI.block(1,Width,Height,1);
    EI(0,0)=EI(1,1);
    EI(Height+1,0)=EI(Height,1);
    EI(0,Width+1)=EI(1,Width);
    EI(Height+1,Width+1)=EI(Height,Width);
    qDebug("计算梯度完成");
    hasAppliedNonMaxDe=false;
    hasAppliedDoubleThreshold=false;
    qDebug()<<"梯度强度最大值="<<EI.maxCoeff()<<"梯度强度最小值="<<EI.minCoeff();
}

void CannyDetecter::ApplyNonMaxDe()
{
    if(hasAppliedNonMaxDe||hasAppliedDoubleThreshold)return;

    /*int tempMax=EI.array().maxCoeff();
    EI*=255;
    EI/=tempMax;*/

    MatrixXi treated(Height,Width);
    treated.setZero();

    MatrixXf weight=(EdgeDir.cast<float>()*M_PI/180.0f).array().sin();

    MatrixXi indexOffset=EdgeDir/45;
    indexOffset=(indexOffset.array()<4).select(indexOffset,0);
    //std::cout<<rotation<<std::endl;
    int BegOffsetR[8]={0,-1,-1,0,0,1,1,0};
    int BegOffsetC[8]={1,0,0,-1,-1,0,0,1};
    int EndOffsetR[8]={-1,-1,-1,-1,1,1,1,1};
    int EndOffsetC[8]={1,1,-1,-1,-1,-1,1,1};

    int BegVal=0,EndVal=0;
    float dT1=0,dT2=0;
    for(int r=0;r<Height;r++)
        for(int c=0;c<Width;c++)
        {
            if(!EI(r+1,c+1))continue;
            BegVal=EI(1+r+BegOffsetR[indexOffset(r,c)],1+c+BegOffsetC[indexOffset(r,c)]);
            EndVal=EI(1+r+EndOffsetR[indexOffset(r,c)],1+c+EndOffsetC[indexOffset(r,c)]);
            dT1=weight(r,c)*BegVal+(1.0f-weight(r,c))*EndVal;
            BegVal=EI(1+r+BegOffsetR[4+indexOffset(r,c)],1+c+BegOffsetC[4+indexOffset(r,c)]);
            EndVal=EI(1+r+EndOffsetR[4+indexOffset(r,c)],1+c+EndOffsetC[4+indexOffset(r,c)]);
            dT2=weight(r,c)*BegVal+(1.0f-weight(r,c))*EndVal;

            if(EI(r+1,c+1)<dT1||EI(r+1,c+1)<dT2)
                treated(r,c)=0;
            else
                treated(r,c)=EI(r+1,c+1);
        }

    EI.block(1,1,Height,Width)=treated;

    hasAppliedNonMaxDe=true;
    hasAppliedDoubleThreshold=false;
    qDebug("非极大值抑制完成");
    qDebug()<<"梯度强度最大值="<<EI.maxCoeff()<<"梯度强度最小值="<<EI.minCoeff();
}

void CannyDetecter::getThreshold(int&Low,int&High,int&Max)
{
    if(hasAppliedDoubleThreshold)return;
    if(!hasAppliedNonMaxDe)return;
    if(step<4)return;
//qDebug("开始计算自适应双阈值");
    auto e_XY=EI_x.array().max(EI_y.array());

    High=(e_XY*EI.array()).sum()/(e_XY.sum());
    Low=High/2;
    Max=EI.maxCoeff();
qDebug("计算自适应双阈值完毕");
}

void CannyDetecter::ApplyDoubleThreshold()
{
    if(hasAppliedDoubleThreshold)return;
    if(!hasAppliedNonMaxDe)return;
    if(step<4)return;

    int Low=-1,High=-1,Max=-1;
    getThreshold(Low,High,Max);
    if(Low>=0&&High>Low&&Max>=High)
    {
        float low=float(Low)/Max;
        float high=float(High)/Max;
        qDebug("将使用自适应双阈值处理图像");
        ApplyDoubleThreshold(low,high);
    }
    return;
}

void CannyDetecter::ApplyDoubleThreshold(float low, float high)
{
    if(hasAppliedDoubleThreshold)return;
    if(!hasAppliedNonMaxDe)return;
    if(step<4)return;
    if(low<0.0f||high<low||low>1.0f||high>1.0f)
    {
        ApplyDoubleThreshold();
        return;
    }

    int Max=EI.maxCoeff();
    int Low=low*Max;
    int High=high*Max;
    qDebug()<<"低阈值："<<Low<<"；高阈值："<<High;
    //Edges;
    Edges.setZero(Height+2,Width+2);
qDebug("rue");
    Edges=(EI.array()>Low).select(EI,0);
    Edges=(EI.array()<High).select(Edges,Max);
/*
qDebug("将要拓宽EI矩阵");
    EI<<MatrixXi::Zero(1,Width+2),
          MatrixXi::Zero(Height,1),EI,MatrixXi::Zero(Height,1),
          MatrixXi::Zero(1,Width+2);
qDebug("成功拓宽EI矩阵");*/
qDebug("rua!");
qDebug()<<"EI尺寸："<<EI.rows()<<"×"<<EI.cols();
qDebug()<<"Edges尺寸："<<Edges.rows()<<"×"<<Edges.cols();

    for(int r=0;r<Height;r++)
        for(int c=0;c<Width;c++)
        {
            if(!Edges(r+1,c+1)||Edges(r+1,c+1)>=Max)
                continue;
            //qDebug()<<"r="<<r<<";  c="<<c<<";";
            Edges(r+1,c+1)=(EI.block(r,c,3,3).array()>=Max).any()*Max;
        }

    //EI=Edges;
    qDebug("双阈值检测连接边缘已完成");
    hasAppliedDoubleThreshold=true;
    return;
}
