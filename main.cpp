/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/18/2015 07:17:24 AM
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
#include <complex>
#include <string>
#include <fstream>
#include "DefClass.h"
#include "mkl.h"

int main (int argc, char *argv[])
{
    int tN,bMode,N,numdt;
    double h,alpha,wx,phi;

    if (argc>1)
    {
        tN=(int)atof(argv[1]);
        bMode=(int)atof(argv[2]);
        h=(double)atof(argv[3]);
        N=(int)atof(argv[4]);
        alpha=(double)atof(argv[5]);
        wx=(double)atof(argv[6]);
        phi=(double)atof(argv[7]);
        numdt=(int)atof(argv[8]);
    }

    RNHQS one(tN,bMode,h,N,alpha,wx,phi,numdt);         /* bMode=2: Period Boundary */
//    RNHQS one(4,1,0.0005,1001,1);
    one.go();

    return 0;
}
