/*
 * =====================================================================================
 *
 *       Filename:  DefClass.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/17/2015 03:08:28 PM
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Dr. Fritz Mehner (fgm), mehner.fritz@fh-swf.de
 *   Organization:  FH Südwestfalen, Iserlohn
 *
 * =====================================================================================
 */

#ifndef DEFCLASS_H
#define DEFCLASS_H
#include <fstream>
#include <complex>

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
        void onestep();
        double k_,ws_,wx_,p_,mu_,alpha_,f_,w_,phi_,normE_;
        double dx_,t_,h_;
        double *V0_,*V1_,*VI_,dV_,*xp_;
        int N_,bMode_,oMode_;

        void LP_initital();
        void LP_onestep();
        double *LPE11_,*LPE12_,*LPE21_,*LPE22_,*LPV1_,*LPV2_,*LPtemp_,*LPtemp2_;
};

class RNHQS : protected RK4
{
    public:
        RNHQS(int,int,int,double,int,double,double,double,double,int);
        void go();
        void record();
        void disp();
    private:
        std::ofstream outE_,outxp_;
        int numdt_;
};
#endif
