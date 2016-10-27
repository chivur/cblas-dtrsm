#!/bin/bash
module load libraries/atlas-3.10.1-gcc-4.4.6-opteron
gcc main.c -o dtrsm_opteron -L/opt/tools/libraries/atlas/3.10.1-opteron-gcc-4.4.6/lib -lcblas -latlas -lm -std=c99 -D OPTERON
for i in {1..5..1}
do
./dtrsm_opteron <size$i.txt
done

