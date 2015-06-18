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
    RNHQS one(3,0.00001,401,1);
//    RNHQS one(3,0.0001,101,1);
    one.go();

    return 0;
}
