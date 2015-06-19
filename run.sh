#!/bin/bash

tN=6
bMode=2
# 1: closed; 2: period 
h=0.005
N=1001
wx=1.0
phi=1.570796326794897
numdt=5000000
./a.out $tN $bMode $h $N $wx $phi $numdt

if [ -e "EF$I.mat" ]
then
    rm EF$I.mat
fi

txt2mat EF$I.txt EF$I EF$I.mat
rm EF$I.txt
