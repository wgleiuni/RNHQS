#!/bin/bash

tN=12
bMode=1
# 1: closed; 2: period 
oMode=3
# 1: E; 2: x & p; 3: E & x & p
h=0.002
N=101
mu=0.07
alpha=0.014
wx=0.3
f=0.0
phi=0
numdt=100000
nI=2
W=0.0
w=0.2168
./a.out $tN $bMode $oMode $h $N $mu $alpha $wx $f $phi $numdt $nI $W $w

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
