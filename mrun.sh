#!/bin/bash

tN=12
bMode=3
# 1: closed; 2: period; 3: mean
oMode=1
h=0.002
gama=0.1
v=0.0501
A=1
omega=10
f=0.25
phi=0.5
numdt=20000000
W=0.0
./a.out $tN $bMode $oMode $h $gama $v $A $omega $f $phi $numdt $W

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
