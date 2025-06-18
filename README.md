# root_growth
Fractional-order moisture transport solver with root growth simulation

Moisture transport and root growth solver with PSO for parameters identification:

mt2d_root_growing_hg.cpp - main file and moisture transport model + PSO for parameters identification

root_growing.h - root growing model

test2.sh - a script to run different parameters identification scenarios

test.sh - a script to run checks of solution quality with parameter values obtained by test2.sh


Root system optimization:

rg_genetic.py - a genetic algorithm for the selection of parameters for the root system growth model

run_once.sh - a script to run 1 simulation of root growth (used by rg_genetic.py)

mt2d_run - a script to run the solver (used by run_once.sh, test.sh, test2.sh)

rg_params.json - a file with genetic algorithm parameters (used by rg_genetic.py)

run_genetic.sh - a script to run genetic algorithm and check the quality of solution


Input files:

vgm_makariv_2023_20.txt - parameters of van Genuchten and Mualem models for a soil in a field in Makariv district, Kyiv region, Ukraine

fit_input_.txt - water heads, precipitation and irrigation rates, ET measured in 2024 in the field in Makariv district mid-season

fit_input_beginning.txt - ET measured in 2024 in the field in Makariv district from the beginning of observations

Rl.txt - default root system depths

pv_.txt - default values of moisture transport model's parameters (used for further parameters identification)

pv_clean.txt - values of moisture transport model's parameters without calibration (used as the `base case` for comparison)


Results of simulations:

pv_best_fit0.txt - integer-order moisture transport model's parameters without fitting filtration coefficients

pv_best_fit1.txt - integer-order moisture transport model's parameters fitting filtration coefficients

pv_best_fit2.txt - fractional-order moisture transport model's parameters fitting filtration coefficients

pv_best_fit1a.txt - integer-order moisture transport model's parameters fitting filtration coefficients (relative error fitness function)

pv_best_fit1b1.txt - integer-order moisture transport model's parameters fitting filtration coefficients (heterogeneous case)

pv_best_fit1b2.txt - integer-order moisture transport model's parameters fitting filtration coefficients (heterogeneous case + theta_r)

root_systems_best1.txt - optimized root system form for the integer-order moisture transport model
root_systems_best2.txt - optimized root system form for the fractional-order moisture transport model

