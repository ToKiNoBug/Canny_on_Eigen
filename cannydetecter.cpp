#include "cannydetecter.h"
#include <cmath>
#include <iostream>
#include <vector>
Matrix3i CannyDetecter::GauKernel=Matrix3i::Zero(3,3);
int CannyDetecter::GauKernelSum=0;
Matrix3i CannyDetecter::SobelKernelX=Matrix3i::Zero(3,3);
Matrix3i CannyDetecter::SobelKernelY=Matrix3i::Zero(3,3);

using namespace std;

CannyDetecter::CannyDetecter()
{
//Gray=QImage();
//Raw=QImage();
step=0;
parent=NULL;
}

void CannyDetecter::Load(const MatrixXi & Src)
{
    Height=Src.rows();
    Width=Src.cols();
    GrayPic.setZero(Height+2,Width+2);
    GrayPic.block(1,1,Height,Width)=Src;
    GrayPic.block(0,1,1,Width)=GrayPic.block(1,1,1,Width);
    GrayPic.block(Height+1,1,1,Width)=GrayPic.block(Height,1,1,Width);
    GrayPic.block(1,0,Height,1)=GrayPic.block(1,1,Height,1);
    GrayPic.block(1,Width+1,Height,1)=GrayPic.block(1,Width,Height,1);
    GrayPic(0,0)=GrayPic(1,1);
    GrayPic(Height+1,0)=GrayPic(Height,1);
    GrayPic(0,Width+1)=GrayPic(1,Width);
    GrayPic(Height+1,Width+1)=GrayPic(Height,Width);
    qDebug("成功加载灰度图为Eigen矩阵并处理边界");
    step=1;
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
    step=1;
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
    step=2;
}

void CannyDetecter::ApplySobel()
{
    if(step<2)return;
    EI.setZero(Height+2,Width+2);
    EdgeDir.setZero(Height+2,Width+2);
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
            EdgeDir(r,c)=round((atan2(Gy,Gx)/M_PI)*180.0f);
        }

    EdgeDir=(EdgeDir.array()>=0).select(EdgeDir,EdgeDir.array()+180);
    EdgeDir=(EdgeDir.array()!=180).select(EdgeDir,0);
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
    step=3;
}

void CannyDetecter::ApplyNonMaxDe()
{
    if(step<3)return;
    if(hasAppliedNonMaxDe||hasAppliedDoubleThreshold)return;
    qDebug()<<"minDir:"<<EdgeDir.minCoeff()<<"; maxDir:"<<EdgeDir.maxCoeff();


    //weight.setZero(Height+2,Width+2);
    auto EdgeDirMod45=EdgeDir.array()-(EdgeDir.array()/45)*45;
    MatrixXf weight=(EdgeDirMod45.cast<float>().array()*M_PI/180.0f).tan();
    qDebug()<<"minWeight:"<<weight.minCoeff()<<"; maxWeight:"<<weight.maxCoeff();
    MatrixXi treated=EI;treated.setZero();
    int dT1=0,dT2=0;
    int BegI,EndI;
    for(int r=1;r<=Height;r++)
        for(int c=1;c<Width;c++)
        {
            if(EI(r,c)<=1)
            {
                treated(r,c)=0;
                continue;
            }
            if(EdgeDir(r,c)>=0&&EdgeDir(r,c)<45)
            {
                BegI=EI(r,c+1);EndI=EI(r-1,c+1);
                dT1=BegI*weight(r,c)+EndI*(1-weight(r,c));
                BegI=EI(r,c-1);EndI=EI(r+1,c-1);
                dT2=BegI*weight(r,c)+EndI*(1-weight(r,c));

                treated(r,c)=(EI(r,c)<dT1||EI(r,c)<dT2)?0:EI(r,c);
                continue;
            }
            if(EdgeDir(r,c)>=45&&EdgeDir(r,c)<90)
            {
                BegI=EI(r-1,c+1);EndI=EI(r-1,c);
                dT1=BegI*weight(r,c)+EndI*(1-weight(r,c));
                BegI=EI(r+1,c-1);EndI=EI(r+1,c);
                dT2=BegI*weight(r,c)+EndI*(1-weight(r,c));

                treated(r,c)=(EI(r,c)<dT1||EI(r,c)<dT2)?0:EI(r,c);
                continue;
            }
            if(EdgeDir(r,c)>=90&&EdgeDir(r,c)<135)
            {
                BegI=EI(r-1,c);EndI=EI(r-1,c-1);
                dT1=BegI*weight(r,c)+EndI*(1-weight(r,c));
                BegI=EI(r+1,c);EndI=EI(r+1,c+1);
                dT2=BegI*weight(r,c)+EndI*(1-weight(r,c));

                treated(r,c)=(EI(r,c)<dT1||EI(r,c)<dT2)?0:EI(r,c);
                continue;
            }
            //if(EdgeDir(r,c)>=135&&EdgeDir(r,c)<=180)
            //{
                BegI=EI(r-1,c-1);EndI=EI(r,c-1);
                dT1=BegI*weight(r,c)+EndI*(1-weight(r,c));
                BegI=EI(r+1,c+1);EndI=EI(r,c+1);
                dT2=BegI*weight(r,c)+EndI*(1-weight(r,c));

                treated(r,c)=(EI(r,c)<dT1||EI(r,c)<dT2)?0:EI(r,c);
                //continue;
            //}
        }
    EI=treated;
    /*MatrixXi treated(Height,Width);
    treated.setZero();

    //MatrixXf weight=(EdgeDir.cast<float>()*M_PI/180.0f).array().sin();
    MatrixXf weight=(EdgeDir.cast<float>()*M_PI/180.0f).array().sin();
            //EI_y.cast<float>().array().abs()/(EI_x.cast<float>().array().abs()+1e-10);

    qDebug()<<"minWeight:"<<weight.minCoeff()<<"; maxWeight:"<<weight.maxCoeff();

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
            if(EI(r+1,c+1)<=1)
            {
                treated(r,c)=0;
                continue;
            }
            BegVal=EI(1+r+BegOffsetR[indexOffset(r,c)],1+c+BegOffsetC[indexOffset(r,c)]);
            EndVal=EI(1+r+EndOffsetR[indexOffset(r,c)],1+c+EndOffsetC[indexOffset(r,c)]);
            dT1=weight(r,c)*BegVal+(1.0f-weight(r,c))*EndVal;
            BegVal=EI(1+r+BegOffsetR[4+indexOffset(r,c)],1+c+BegOffsetC[4+indexOffset(r,c)]);
            EndVal=EI(1+r+EndOffsetR[4+indexOffset(r,c)],1+c+EndOffsetC[4+indexOffset(r,c)]);
            dT2=weight(r,c)*BegVal+(1.0f-weight(r,c))*EndVal;



            if(EI(r+1,c+1)>dT1&&EI(r+1,c+1)>dT2)
                treated(r,c)=EI(r+1,c+1);
            else
                treated(r,c)=0;
        }

    EI.block(1,1,Height,Width)=treated;*/

    hasAppliedNonMaxDe=true;
    hasAppliedDoubleThreshold=false;
    qDebug("非极大值抑制完成");
    qDebug()<<"梯度强度最大值="<<EI.maxCoeff()<<"梯度强度最小值="<<EI.minCoeff();
    step=4;
}

int CannyDetecter::OTSU(float &lowMean, float &lowVariance)
{
    //get Histogram
    VectorXf HistFrequency,HistValue;
    VectorXf &p_i=HistFrequency;
    int Max=EI.block(1,1,Height,Width).maxCoeff();
    HistFrequency.setZero(Max);
    HistValue.setLinSpaced(Max,1,Max);
    for(int r=1;r<=Height;r++)
        for(int c=1;c<=Width;c++)
        {
            if(EI(r,c)>0)
                HistFrequency(EI(r,c)-1)++;
        }
    HistFrequency/=HistFrequency.sum();
    VectorXf p1;p1.setZero(Max);p1(0)=HistFrequency(0);
    for(int val=2;val<=Max;val++)
        p1(val-1)=HistFrequency(val-1)+p1(val-2);
    p1/=p1.maxCoeff();
    VectorXf m_k;m_k.setZero(Max);m_k(0)=p_i(0)*HistValue(0);
    for(int val=2;val<=Max;val++)
        m_k(val-1)=HistValue(val-1)*p_i(val-1)+m_k(val-2);

    float m_G=m_k.maxCoeff();

    auto sigma2=(m_G*p1.array()-m_k.array()).square()/(p1.array()*(1.0f-p1.array())+1e-10f);

    int OtsuVal=0;
    sigma2.maxCoeff(&OtsuVal);
    OtsuVal++;
    qDebug()<<"OTSU="<<OtsuVal;
    //Low area index=0->Otsu-1

    lowMean=(p_i.segment(0,OtsuVal-1).array()*HistValue.segment(0,OtsuVal-1).array()).sum()/p_i.segment(0,OtsuVal-1).sum();
    qDebug()<<"低区段均值："<<lowMean;
    auto diffSquare=(HistValue.segment(0,OtsuVal-1).array()-lowMean).square();

    lowVariance=(diffSquare*p_i.segment(0,OtsuVal-1).array()).sum()/p_i.segment(0,OtsuVal-1).sum();



    qDebug()<<"低区段方差："<<lowVariance;
    return OtsuVal;
}

void CannyDetecter::getThreshold(int&Low,int&High,int&Max)
{
    //if(hasAppliedDoubleThreshold)return;
    if(!hasAppliedNonMaxDe)return;
    if(step<4)return;
    Max=EI.maxCoeff();
//qDebug("开始计算自适应双阈值");
    /*auto e_XY=EI_x.array().max(EI_y.array());

    High=(e_XY*EI.array()).sum()/(e_XY.sum());
    Low=High/2;*/
    /*
    VectorXi Hist;
    VectorXf accumAvg,accumCount;//不统计0
    Hist.setZero(Max);
    accumAvg.setZero(Max);
    accumCount.setZero(Max);//val:1->Max; index:0->Max-1
    for(int r=1;r<=Height;r++)
        for(int c=1;c<=Width;c++)
        {
            if(EI(r,c)>0)
                Hist(EI(r,c)-1)++;
        }
    accumAvg(0)=1*Hist(0);
    accumCount(0)=Hist(0);
    for(int val=2;val<=Max;val++)
    {
            accumAvg(val-1)=val*Hist(val-1)+accumAvg(val-2);
            accumCount(val-1)=Hist(val-1)+accumCount(val-2);
    }
    accumAvg.array()/=accumCount.array()+1e-10f;
    int globalAvg=accumAvg(Max-1);
    int OTSU=0;
    auto p1=accumCount/accumCount.sum();

    auto sigma2=(globalAvg*p1.array()-accumAvg.array())/(p1.array()*(1.0f-p1.array())+1e-4f);
    sigma2.maxCoeff(&OTSU);
    OTSU++;
qDebug()<<"OTSU="<<OTSU;
    VectorXf ValVec;
    ValVec.setLinSpaced(Max,1,Max);

    float miu=(ValVec.array()*Hist.cast<float>().array()).segment(0,OTSU-1).mean();
    float sigmaSquare=((ValVec.segment(0,OTSU-1).array()-miu).square()*Hist.segment(0,OTSU-1).cast<float>().array()).sum()/Hist.segment(0,OTSU-1).cast<float>().sum();
qDebug()<<"miu="<<miu;
qDebug()<<"sigma^2="<<sigmaSquare;*/

    float lowMean,lowVariance;
    OTSU(lowMean,lowVariance);
    Low=lowMean-sqrt(lowVariance);
    High=lowMean+sqrt(lowVariance);
    qDebug()<<"低阈值="<<Low;
    qDebug()<<"高阈值="<<High;

qDebug("计算自适应双阈值完毕");
}

void CannyDetecter::ApplyDoubleThreshold()
{
    //if(hasAppliedDoubleThreshold)return;
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
    //if(hasAppliedDoubleThreshold)return;
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
    qDebug()<<"低阈值："<<100.0f*Low/Max<<"%；高阈值："<<100.0f*High/Max<<"%";
    //Edges;
    Edges=EI;
qDebug("rue");
    Edges=(Edges.array()>Low).select(Edges,0);
    Edges=(Edges.array()<High).select(Edges,Max);

//qDebug("rua!");
//qDebug()<<"EI尺寸："<<EI.rows()<<"×"<<EI.cols();
//qDebug()<<"Edges尺寸："<<Edges.rows()<<"×"<<Edges.cols();

    for(int r=0;r<Height;r++)
        for(int c=0;c<Width;c++)
        {
            if(!Edges(r+1,c+1)||Edges(r+1,c+1)>=High)
                continue;
            //qDebug()<<"r="<<r<<";  c="<<c<<";";
            Edges(r+1,c+1)=(EI.block(r,c,3,3).array()>=High).any()?Max:0;
        }

    //EI=Edges;
    qDebug("双阈值检测连接边缘已完成");
    hasAppliedDoubleThreshold=true;
    return;
}

void CannyDetecter::ApplySingleThreshold(float thre)
{
    if(!hasAppliedNonMaxDe)return;
    if(step<4)return;

    int Max=EI.maxCoeff();
    int T=thre*Max;
    Edges=(EI.array()>=T).select(MatrixXi::Ones(Height+2,Width+2),0);
    step=5;
}

void CannyDetecter::ApplyOTSUThreshold()
{
    if(!hasAppliedNonMaxDe)return;
    if(step<4)return;
    float a,b;
    ApplySingleThreshold(float(OTSU(a,b))/EI.maxCoeff());
}
