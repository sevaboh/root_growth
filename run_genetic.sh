#!/bin/sh
#g++ -O3 -fopenmp --param max-inline-insns-auto=4000 mt2d_root_growing_hg.cpp -o mt2d

./rg_genetic.py > genetic_log1.txt

./run_once.sh best `tail -n 2 genetic_log1.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`

./run_once.sh 0 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
./run_once.sh 1 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
./run_once.sh 2 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
./run_once.sh 3 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
./run_once.sh 4 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
./run_once.sh 5 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
./run_once.sh 6 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
./run_once.sh 7 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
./run_once.sh 8 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
./run_once.sh 9 `tail -n 2 genetic_log.txt | head -n 1 | sed 's/\,//g' | sed 's/\[//g' | sed 's/\]//g'`
