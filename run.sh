#!/bin/bash

tN=1
bMode=2
# 1: closed; 2: period 
oMode=3
# 1: E; 2: x & p; 3: E & x & p
h=0.002
N=1000
mu=0.07
alpha=0.00005
wx=0.3
f=0.0
phi=0
numdt=10000
nI=2
W=0.0
w=0.2168
#./a.out $tN $bMode $oMode $h $N $mu $alpha $wx $f $phi $numdt $nI $W $w
sigma=2.1
ri=0.5
w=0.0
h=0.0063
phi=0.5

./a.out $tN $h $sigma $ri $w $phi $numdt

if [ -e "EF$tN.mat" ]
then
    rm EF$tN.mat
fi

if [ -e "XP$tN.mat" ]
then
    rm XP$tN.mat
fi

if [ -e "ZB$tN.mat" ]
then
    rm ZB$tN.mat
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
