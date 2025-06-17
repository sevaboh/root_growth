#!/bin/bash
ulimit -c unlimited
#OMP_THREAD_LIMIT=8
#export OMP_THREAD_LIMIT

#simulate root system growing from the beginning
if [ $1 != "best" ]; then
time ./mt2d_run 0.3 1 1 50 1_beginning_$1 6400000 3500 70 0 1 "fl 2468 fi_file fit_input_beginning.txt pv_file pv_best_fit10.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 25 max_tau_0 200 rpf 3 \
rg_r $2 rg_la $3 rg_lb $4 rg_ln $5 rg_N $6 rg_angle $7 rg_max_len $8 rg_w1 $9 root_systems_out_file root_systems_$1.txt"
fi
# simulate actual measurements with generated root system
time ./mt2d_run 0.3 1 1 50 1_$1${11} 1400000 3500 70 0 1 "fl 2468 fi_file fit_input_${10}.txt pv_file pv_best_fit10.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 200 max_tau_0 200 rpf 3 root_systems_file root_systems_$1.txt no_growth 1 \
rg_r $2 rg_la $3 rg_lb $4 rg_ln $5 rg_N $6 rg_angle $7 rg_max_len $8 rg_w1 $9"
grep "#checking err" log_fit1_$1.txt
