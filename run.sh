#!/bin/bash

tN=12
bMode=2
# 1: closed; 2: period 
oMode=3
# 1: E; 2: x & p; 3: E & x & p
h=0.002
N=1000
mu=0.07
alpha=0.0001
wx=0.3
f=0.25;
phi=0
numdt=5000000
nI=4
./a.out $tN $bMode $oMode $h $N $mu $alpha $wx $f $phi $numdt $nI

if [ -e "EF$tN.mat" ]
then
    rm EF$tN.mat
fi

if [ -e "XP$tN.mat" ]
then
    rm XP$tN.mat
fi

#if [$oMode == 1]
#then
#    txt2mat EF$tN.txt EF$tN EF$tN.mat
#    rm EF$tN.txt
#fi

#if [$oMode == 2]
#then
#    txt2mat XP$tN.txt XP$tN XP$tN.mat
#    rm XP$tN.txt
#fi

#if [$oMode == 3]
#then
#    txt2mat EF$tN.txt EF$tN EF$tN.mat
#    txt2mat XP$tN.txt XP$tN XP$tN.mat
#    rm EF$tN.txt XP$tN.txt
#fi
