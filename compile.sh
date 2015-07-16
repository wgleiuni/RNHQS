#!/bin/bash

icpc main.cpp DefClass.cpp DefClass.h -mkl -O3
#icpc -g -shared-intel -debug -O0 main.cpp DefClass.cpp DefClass.h -mkl
