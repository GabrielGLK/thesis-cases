#!/bin/sh

export OMP_NUM_THREADS=4
qcc -Wall -O3 -fopenmp rising.c -o rising  -L$BASILISK/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
./rising>>log 2>>out

