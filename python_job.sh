#!/bin/bash
#Generate qsub file with arguments from Python.
tN=$1
h=$2
sigma=$3
ri=$4
omega=$5
phi=$6
numdt=$7
echo "#PBS -l nodes=1:ppn=1 
#PBS -j oe
#PBS -o log
#PBS -l walltime=240:00:00
cd /home/glwang/RNHQS/

./a.out $tN $h $sigma $ri $omega $phi $numdt
txt2mat ZB$I.txt ZB$I ZB$I.mat
rm ZB$I.txt" > single.sh
#qsub single.sh
