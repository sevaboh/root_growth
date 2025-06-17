#!/bin/sh
ulimit -c unlimited
#OMP_THREAD_LIMIT=8
#export OMP_THREAD_LIMIT
#g++ -O3 -fopenmp --param max-inline-insns-auto=4000 mt2d_root_growing_hg.cpp -o mt2d
#exit
#fl 2468
#fl 2342 - pic gen
#to_fit 496
# calibrate the model
time ./mt2d_run 0.3 1 1 50 0 800000 3500 70 0 0 "Fit 1 to_fit 464 bounds_11_0 10 bounds_11_1 10 fl 2468 fi_file fit_input_.txt pv_file pv_.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400 ls_min_tau 100"
time ./mt2d_run 0.3 1 1 50 1 800000 3500 70 0 0 "Fit 1 to_fit 496 bounds_11_0 10 bounds_11_1 10 fl 2468 fi_file fit_input_.txt pv_file pv_.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400 ls_min_tau 100"
time ./mt2d_run 0.3 1 1 50 2 800000 3500 70 0 0 "Fit 1 to_fit 499 bounds_11_0 10 bounds_11_1 10 fl 2468 fi_file fit_input_for_fr.txt pv_file pv_.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400 ls_min_tau 100"

time ./mt2d_run 0.3 1 1 50 1a 800000 3500 70 0 0 "fit_func 1 Fit 1 to_fit 496 bounds_11_0 10 bounds_11_1 10 fl 2468 fi_file fit_input_.txt pv_file pv_.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400 ls_min_tau 100"
time ./mt2d_run 0.3 1 1 50 1b1 800000 3500 70 0 0 "fit_vgm 1 Fit 1 to_fit 496 bounds_11_0 10 bounds_11_1 10 fl 2468 fi_file fit_input_.txt pv_file pv_.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400 ls_min_tau 100"
time ./mt2d_run 0.3 1 1 50 1b2 800000 3500 70 0 0 "fit_vgm 2 Fit 1 to_fit 496 bounds_11_0 10 bounds_11_1 10 fl 2468 fi_file fit_input_.txt pv_file pv_.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400 ls_min_tau 100"
time ./mt2d_run 0.3 1 1 50 1b3 800000 3500 70 0 0 "fit_vgm 3 Fit 1 to_fit 496 bounds_11_0 10 bounds_11_1 10 fl 2468 fi_file fit_input_.txt pv_file pv_.txt rl_file Rl.txt vgm vgm_makariv_2023_20.txt ls_eps 1e-13 ls_max_iter 500 max_tau_irr 400 max_tau_0 400 ls_min_tau 100"
