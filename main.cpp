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
    double h,alpha,wx,f,phi;

    if (argc>1)
    {
        tN=(int)atof(argv[1]);
        bMode=(int)atof(argv[2]);
        oMode=(int)atof(argv[3]);
        h=(double)atof(argv[4]);
        N=(int)atof(argv[5]);
        alpha=(double)atof(argv[6]);
        wx=(double)atof(argv[7]);
        f=(double)atof(argv[8]);
        phi=(double)atof(argv[9]);
        numdt=(int)atof(argv[10]);
        nI=(int)atof(argv[11]);
    }

    RNHQS one(tN,bMode,oMode,h,N,alpha,wx,f,phi,numdt,nI);         /* bMode=2: Period Boundary */
//    RNHQS one(4,1,0.0005,1001,1);
    one.go();

    return 0;
}
