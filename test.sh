#!/bin/sh
ulimit -c unlimited
#OMP_THREAD_LIMIT=8
#export OMP_THREAD_LIMIT
#g++ -O3 -fopenmp --param max-inline-insns-auto=4000 mt2d_root_growing_hg.cpp -o mt2d
#exit
#fl 2468
#fl 2342 - pic gen
#to_fit 496
# simulate actual measurements

time ./mt2d_run 0.3 1 1 50 clean 1400000 3500 70 0 1 "fl 2468 fi_file fit_input_.txt pv_file pv_clean.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400"
time ./mt2d_run 0.3 1 1 50 a0 1400000 3500 70 0 1 "fl 2468 fi_file fit_input_.txt pv_file pv_best_fit00.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400"
time ./mt2d_run 0.3 1 1 50 a1 1400000 3500 70 0 1 "fl 2468 fi_file fit_input_.txt pv_file pv_best_fit10.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400"
time ./mt2d_run 0.3 1 1 50 a2 1400000 3500 70 0 1 "fl 2468 fi_file fit_input_.txt pv_file pv_best_fit20.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400"

#cat res_*.js > res_1.js1
#rm res_*.js
