#ifndef DEFCLASS_H
#define DEFCLASS_H
#include <fstream>
#include <complex>
#include "mkl.h"

class RK4
{
    public:
        RK4();
    protected:
        std::complex<double> *E_,*E,*V_,I_,*K_[4],meanx_,meanp_;
        void initial();
        void dE(std::complex<double> *,std::complex<double> *,double);
        void getV(double);
        void getxp();
        double getDisorder();
        void onestep();
        double k_,ws_,wx_,p_,mu_,alpha_,f_,w_,phi_,normE_,W_;
        double dx_,t_,h_;
        double *V0_,*V1_,*VI_,dV_,*xp_,*disW_;
        int N_,bMode_,oMode_,nI_,nW_,tN_;

        void LP_initital();
        void LP_onestep();
        double *LPE11_,*LPE12_,*LPE21_,*LPE22_,*LPV1_,*LPV2_,*LPtemp_,*LPtemp2_,*LPV2nD_;

        VSLStreamStatePtr stream_;
};

class RNHQS : protected RK4
{
    public:
        RNHQS(int,int,int,double,int,double,double,double,double,double,int,int,double,double);
        void go();
        void record();
        void disp();
    private:
        std::ofstream outE_,outxp_;
        int numdt_;
};

class RNHQS_mean
{
    public:
        RNHQS_mean(int,int,int,double,double,double,double,double,double,double,int,double);
        void go();
        void record();
        void disp();
    private:
        std::ofstream out_;
        void onestep();
        double getDis();
        void getW_();
        double dC(int,double,double*);
        
        int tN_,bMode_,oMode_,numdt_,i,nW_;
        double h_,gama_,v_,A_,w_,f_,phi_,t_,W_,Wr_,Wi_;
        double *c_,*disW_;

        VSLStreamStatePtr stream_;
};
#endif
