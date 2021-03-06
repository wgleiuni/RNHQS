#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <string>
#include <fstream>
#include "DefClass.h"
#include "mkl.h"

int main (int argc, char *argv[])
{
    int tN,bMode,oMode,N,numdt,nI;
    double h,mu,alpha,wx,f,phi,W,w,sigma,ri;

    if (argc>1)
    {
        if (argc==15)
        {
            tN=(int)atof(argv[1]);
            bMode=(int)atof(argv[2]);
            oMode=(int)atof(argv[3]);
            h=(double)atof(argv[4]);
            N=(int)atof(argv[5]);
            mu=(double)atof(argv[6]);
            alpha=(double)atof(argv[7]);
            wx=(double)atof(argv[8]);
            f=(double)atof(argv[9]);
            phi=(double)atof(argv[10]);
            numdt=(int)atof(argv[11]);
            nI=(int)atof(argv[12]);
            W=(double)atof(argv[13]);
            w=(double)atof(argv[14]);
            RNHQS one(tN,bMode,oMode,h,N,mu,alpha,wx,f,phi,numdt,nI,W,w);         /* bMode=2: Period Boundary */
            //    RNHQS one(4,1,0.0005,1001,1);
            one.go();
        }
        else if (argc==13)
        {
            tN=(int)atof(argv[1]);
            bMode=(int)atof(argv[2]);
            oMode=(int)atof(argv[3]);
            h=(double)atof(argv[4]);
            wx=(double)atof(argv[5]);
            mu=(double)atof(argv[6]);
            alpha=(double)atof(argv[7]);
            w=(double)atof(argv[8]);
            f=(double)atof(argv[9]);
            phi=(double)atof(argv[10]);
            numdt=(int)atof(argv[11]);
            W=(double)atof(argv[12]);
            RNHQS_mean two(tN,bMode,oMode,h,wx,mu,alpha,w,f,phi,numdt,W);
            two.go();
        }
        else if (argc==8)
        {
            tN=(int)atof(argv[1]);
            h=(double)atof(argv[2]);
            sigma=(double)atof(argv[3]);
            ri=(double)atof(argv[4]);
            w=(double)atof(argv[5]);
            phi=(double)atof(argv[6]);
            numdt=(int)atof(argv[7]);
            RNHQS_ZB three(tN,h,sigma,ri,w,phi,numdt);
            three.go();
        }
    }

    return 0;
}
