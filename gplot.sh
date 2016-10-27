#!/bin/bash
gnuplot -persist -e "set title \"Nehalem vs Opteron\";
set ylabel \"size (MxN)\";
set xlabel \"time (s) \"; 
plot \"time_normal.txt\" with lines title \"normal_neh\" ,
\"time_opt.txt\" with lines title \"opt_neh\", 
\"time_blas.txt\" with lines title \"Blas_neh\",
\"time_normal_o.txt\" with lines title \"Normal_opteron \",
\"time_opt_o.txt\" with lines title \" Opt_opteron\",
\"time_blas_o.txt\" with lines title \"Blas_opteron\";
"

#rm time_normal.txt
#rm time_normal_o.txt
#rm time_opt.txt
#rm time_opt_o.txt
#rm time_blas.txt
#rm time_blas_o.txt
