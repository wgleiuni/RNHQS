#!/bin/bash

tN=9
bMode=2
# 1: closed; 2: period 
h=0.0005
N=1001
alpha=0.001
wx=1.0
phi=0
numdt=5000000
./a.out $tN $bMode $h $N $alpha $wx $phi $numdt

if [ -e "EF$tN.mat" ]
then
    rm EF$tN.mat
fi

txt2mat EF$tN.txt EF$tN EF$tN.mat
rm EF$tN.txt
