/*
 * =====================================================================================
 *
 *       Filename:  DefClass.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/17/2015 03:25:49 PM
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Dr. Fritz Mehner (fgm), mehner.fritz@fh-swf.de
 *   Organization:  FH Südwestfalen, Iserlohn
 *
 * =====================================================================================
 */

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex>
#include <fstream>
#include "DefClass.h"
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_vsl.h"

RK4::RK4()
{
}

void RK4::initial()
{
    double *x=new double [N_];
    int i,j,Nw;
    
    Nw=10;
    for (i=0;i<N_;i++)
    {
        *(x+i)=(-1.0*Nw+2.0*Nw/(N_-1)*i)*ws_;
    }
    dx_=*(x+1)-*(x+0);

    V0_=new double [N_];
    V1_=new double [N_];
    VI_=new double [N_];

    for (i=0;i<N_;i++)
    {
        *(V0_+i)=0.0;
        *(V1_+i)=0.0;
        *(VI_+i)=0.0;
    }

    for (i=0;i<N_;i++)
    {
        for (j=0;j<Nw;j++)
        {
            *(V0_+i)+=-p_*(exp(-pow((*(x+i)+(2*j+1)*ws_/2.0)/wx_,6.0))+exp(-pow((*(x+i)-(2*j+1)*ws_/2.0)/wx_,6.0)));
            *(V1_+i)+=-p_*mu_*pow(-1,j)*(exp(-pow((*(x+i)+(2*j+1)*ws_/2.0)/wx_,6.0))-exp(-pow((*(x+i)-(2*j+1)*ws_/2.0)/wx_,6.0)));
            *(VI_+i)+=-p_*alpha_*pow(-1,j)*(exp(-pow((*(x+i)+(2*j+1)*ws_/2.0)/wx_,6.0))-exp(-pow((*(x+i)-(2*j+1)*ws_/2.0)/wx_,6.0)));
        }
//        std::cout << *(x+i) << " " << *(V0_+i) << " " << *(V1_+i) << " " << *(VI_+i) << std::endl;
    }

    char *uplo="U";
    int *n=new int [1];
    *n=N_;
//    int kla=1;
//    double *a=new double [N_*2];
//    int *lda=new int [1];
//    *lda=2;
    double epsout;
    int loop;
    double *emin=new double [1];
    *emin=-3.0;
    double *emax=new double [1];
    *emax=*emin+1.0;
    int *m0=new int [1];
    *m0=N_;
    double *e=new double[*m0];
    double *outx=new double [*m0*N_];
    int *m=new int [1];
    double *res=new double[*m0];
    int *info=new int [1];
    I_=std::complex<double> (0,1);

    dV_=-1/(2.0*k_)/pow(dx_,2.0);

    double *a;
    int *ia,*ja;

    if (bMode_==1)
    {
        a=new double [2*N_-1];
        ia=new int [N_+1];
        ja=new int [2*N_-1];

        for (i=0;i<N_-1;i++)
        {
            *(a+2*i)=V0_[i]-2.0*dV_;
            *(a+2*i+1)=dV_;
        }
        *(a+2*N_-2)=V0_[N_-1]-2.0*dV_;

        for (i=0;i<N_;i++)
        {
            *(ia+i)=2*i+1;
        }
        *(ia+N_)=2*N_;

        for (i=0;i<N_-1;i++)
        {
            *(ja+2*i)=i+1;
            *(ja+2*i+1)=i+2;
        }
        *(ja+2*N_-2)=N_;
    }
    else if (bMode_==2)
    {
        a=new double [2*N_];
        ia=new int [N_+1];
        ja=new int [2*N_];

        for (i=1;i<N_-1;i++)
        {
            *(a+2*i+1)=V0_[i]-2.0*dV_;
            *(a+2*i+1+1)=dV_;
        }
        *a=V0_[0]-2.0*dV_;
        *(a+1)=dV_;
        *(a+2)=dV_;
        *(a+2*N_-1)=V0_[N_-1]-2.0*dV_;

        for (i=1;i<N_;i++)
        {
            *(ia+i)=2*i+1+1;
        }
        *ia=1;
        *(ia+N_)=2*N_+1;

        for (i=1;i<N_-1;i++)
        {
            *(ja+2*i+1)=i+1;
            *(ja+2*i+1+1)=i+2;
        }
        *ja=1;
        *(ja+1)=2;
        *(ja+2)=N_;
        *(ja+2*N_-1)=N_;
    }

    //    for (i=0;i<N_;i++)
    //    {
    //        a[i*2]=dV_;
    //        a[i*2+1]=(V0_[i]+1/k_/pow(dx_,2.0));
//    }

    int *fpm=new int [128];
    feastinit(fpm);

//    dfeast_sbev(uplo,n,&kla,a,lda,fpm,&epsout,&loop,emin,emax,m0,e,outx,m,res,info);
    dfeast_scsrev(uplo,n,a,ia,ja,fpm,&epsout,&loop,emin,emax,m0,e,outx,m,res,info);
    std::cout << "-----eigenvalue-----" << std::endl;
    for (i=0;i<*m;i++)
    {
        std::cout << *(e+i) << std::endl;
    }

    if (1==2)
    {

        std::cout<< *info << std::endl;
        std::cout<< *m << std::endl;

        std::cout << "-----eigenvalue-----" << std::endl;
        for (i=0;i<*m;i++)
        {
            std::cout << *(e+i) << std::endl;
        }
        std::cout << "-----eigenvector-----" << std::endl;

        for (i=0;i<2*N_;i++)
        {
            std::cout << *(outx+i) << std::endl;
        }
    }

    E_=new std::complex<double> [N_];
    E=new std::complex<double> [N_];
    for (i=0;i<N_;i++)
    {
        *(E_+i)=(-*(outx+i)-*(outx+i+N_))/sqrt(2.0);
    }

    V_=new std::complex<double> [N_];

    for (i=0;i<4;i++)
    {
        K_[i]=new std::complex<double> [N_];
    }

}

void RK4::dE(std::complex<double> *E,std::complex<double> *k,double t)
{
    int i;
    getV(t);

    if (bMode_==1)
    {
        *k=(V_[0]**E+dV_**(E+1))/I_;
        *(k+N_-1)=(V_[N_-1]**(E+N_-1)+dV_**(E+N_-2))/I_;
    }
    else if (bMode_==2)
    {
        *k=(V_[0]**E+dV_*(*(E+1)+*(E+N_-1)))/I_;
        *(k+N_-1)=(V_[N_-1]**(E+N_-1)+dV_*(*(E+N_-2)+*E))/I_;
    }
    for (i=1;i<N_-1;i++)
    {
        *(k+i)=(V_[i]**(E+i)+dV_*(*(E+i-1)+*(E+i+1)))/I_;
    }
}

void RK4::getV(double t)
{
    int i;
    for (i=0;i<N_;i++)
    {
        V_[i]=V0_[i]+V1_[i]*(sin(w_*t)+f_*sin(2.0*w_*t+phi_))+I_*VI_[i];
    }
}

void RK4::onestep()
{
    dE(E_,K_[0],t_);
    int i;
    for (i=0;i<N_;i++)
    {
        *(E+i)=*(E_+i)+h_/2.0**(K_[0]+i);
    }
    dE(E,K_[1],t_+h_/2.0);
    for (i=0;i<N_;i++)
    {
        *(E+i)=*(E_+i)+h_/2.0**(K_[1]+i);
    }
    dE(E,K_[2],t_+h_/2.0);
    for (i=0;i<N_;i++)
    {
        *(E+i)=*(E_+i)+h_**(K_[2]+i);
    }
    dE(E,K_[3],t_+h_);

    for (i=0;i<N_;i++)
    {
        *(E_+i)=*(E_+i)+h_/6.0*(*(K_[0]+i)+2.0**(K_[1]+i)+2.0**(K_[2]+i)+*(K_[3]+i));
    }
}

void RNHQS::go()
{
    int i;
    for (i=0;i<numdt_;i++)
    {
        t_=h_*i;
        RK4::onestep();
        if (i%1000==0)
        {
            record();
        }
    }
}

void RNHQS::record()
{
    int i;
    for (i=0;i<N_;i++)
    {
        outE_ << std::norm(*(E_+i)) << std::endl;
    }
}

void RNHQS::disp()
{
}

RNHQS::RNHQS(int tN, int bMode, double h, int N, double wx, double phi, int numdt)
{
    t_=0.0;
    bMode_=bMode;
    h_=h;
    numdt_=numdt;
    N_=N;
    k_=1.0;ws_=3.2;wx_=wx;p_=3.0;mu_=0.07;alpha_=0.014;f_=0.25;w_=0.2168;
    phi_=phi;
    RK4::initial();

    char filename[20];

    sprintf(filename,"EF%d.txt",tN);

    outE_.open(filename,std::ostream::out);
}
