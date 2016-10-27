#!/bin/bash
module load libraries/atlas-3.10.1-gcc-4.4.6-nehalem 
gcc main.c -o dtrsm_nehal -L/opt/tools/libraries/atlas/3.10.1-nehalem-gcc-4.4.6/lib -lcblas -latlas -lm -std=c99 -D NEHALEM
for i in {1..5..1}
do
./dtrsm_nehal <size$i.txt
done

