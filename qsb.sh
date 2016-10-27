#!/bin/bash
qsub -cwd -q ibm-nehalem.q nehal.sh
qsub -cwd -q ibm-opteron.q opt.sh
