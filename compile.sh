#!/bin/bash

icpc main.cpp DefClass.cpp DefClass.h -mkl -O2
#icpc -g -shared-intel -debug -O0 main.cpp DefClass.cpp DefClass.h -mkl
