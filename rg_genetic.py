#!/usr/bin/python3
import subprocess
import numpy.random as rnd
import sys
import json
import os
# on the base of the code from https://machinelearningmastery.com/simple-genetic-algorithm-from-scratch-in-python/

# tournament selection
def selection(pop, scores, k=3):
	# first random selection
	selection_ix = rnd.randint(len(pop))
	for ix in rnd.randint(0, len(pop), k-1):
		# check if better (e.g. perform a tournament)
		if scores[ix] < scores[selection_ix]:
			selection_ix = ix
	return pop[selection_ix]

# crossover two parents to create two children
def crossover(p1, p2, r_cross):
	# children are copies of parents by default
	c1, c2 = p1.copy(), p2.copy()
	# check for recombination
	if rnd.rand() < r_cross:
		# linear combination
		r = rnd.rand()
		for i in range(len(c1)):
			c1[i] = r*p1[i] + (1.0 - r)*p2[i]
			c2[i] = r*p2[i] + (1.0 - r)*p1[i]
	return [c1, c2]

# mutation operator
def mutation(p1, r_mut, mins, maxs):
	for i in range(len(p1)):
		# check for a mutation
		if rnd.rand() < r_mut:
			# set to random value
			p1[i] = mins[i]+rnd.rand()*(maxs[i]-mins[i])

# genetic algorithm
def genetic_algorithm(objective_start, objective_collect,n_bits, n_iter, n_pop, r_cross, r_mut, mins, maxs):
	# initial population
	pop = []
	procs = []
	scores = []
	for i in range(n_pop):
		pop.append([])
		procs.append(0)
		scores.append(0)
		for j in range(n_bits):
			if i==0:
				if j==0:
					pop[i].append(mins[j])
				else:
					pop[i].append(maxs[j])
			else:
				if i==1:
					if j==0:
						pop[i].append(maxs[j])
					else:
						pop[i].append(mins[j])
				else:
					pop[i].append(mins[j]+rnd.rand()*(maxs[j]-mins[j]))
	# enumerate generations
	for gen in range(n_iter):
		# evaluate all candidates in the population
		for i in range(len(pop)):
			procs[i]=objective_start(pop[i],i)
			if params["gen_parallel"]!=2:
				scores[i] = objective_collect(pop[i],procs[i],i)[0]
		if params["gen_parallel"]==2:
			for i in range(len(pop)):
				scores[i] = objective_collect(pop[i],procs[i],i)[0]
		if gen == 0:
			best = 0
			best_eval = scores[0]
		# check for new best solution
		for i in range(n_pop):
			if scores[i] < best_eval:
				best, best_eval, best_i = pop[i], scores[i],i
				os.replace("root_systems_"+str(best_i)+".txt","root_systems_best.txt");
		print("%d best f(%s) = %.3f" % (gen,  best, best_eval))
		# select parents
		selected = [selection(pop, scores) for _ in range(n_pop)]
		# create the next generation
		children = list()
		for i in range(0, n_pop, 2):
			# get selected parents in pairs
			p1, p2 = selected[i], selected[i+1]
			# crossover and mutation
			for c in crossover(p1, p2, r_cross):
				# mutation
				mutation(c, r_mut, mins, maxs)
				# store for next generation
				children.append(c)
		# replace population
		pop = children
	return [best, best_eval]

# objective function - run simulation
def objective_start(p,i):
	proc=[]
	if params["gen_parallel"]!=0:
		print("starting "+str(p)+" "+str(i))
		proc.append(subprocess.Popen(["./run_once.sh",str(i),str(float(p[0])),str(float(p[1])),str(float(p[2])),str(float(p[3])),str(float(p[4])),str(float(p[5])),str(float(p[6])),str(float(p[7]))],stdout=subprocess.PIPE))
	return proc

# objective function - get simulation results
def objective_collect_one(p,proc,i):
	out,err = proc.communicate()
	strs=out.decode('utf-8').split("\n")
	# volume per 1 pipeline per 1 m in the last string of output
	total_err=float((strs[-2].split(" "))[2])
	print(str(i)+":("+str(p)+")->"+str(total_err))
	return [total_err]

def objective_collect(p,proc,i):
	s=0
	if params["gen_parallel"]==0:
		print("starting "+str(p)+" "+str(i))
		pr=subprocess.Popen(["./run_once.sh",str(i),str(float(p[0])),str(float(p[1])),str(float(p[2])),str(float(p[3])),str(float(p[4])),str(float(p[5])),str(float(p[6])),str(float(p[7]))],stdout=subprocess.PIPE)
		of=objective_collect_one(p,pr,i)
		s=s+of[0]
	else:
		for i in range(len(proc)):
			of=objective_collect_one(p,proc[i],i)
			s=s+of[0]
	return [s] 

# simulation parameters and ranges (in mt2d_params.json)
#mins, maxs - variable 1 - depth of drip pipeline, variable 2 - number of pipelines per 10 m

params_fname="rg_params.json"
if len(sys.argv)>1:
	params_fname=str(sys.argv[1])
# read parameters json
print("reading params from "+params_fname)
jsf=open(params_fname)
params=json.load(jsf)
jsf.close()
# parse argv for <param_name> <param_value> pairs
if len(sys.argv)>2:
	for i in range(int((len(sys.argv)-2)/2)):
		s=str(sys.argv[2*i+3])
		if s[0].isalpha():
			s="\""+s+"\""
		params[str(sys.argv[2*i+2])]=json.loads(s)
		print("set "+str(sys.argv[2*i+2])+"="+s)
s=json.dumps(params,indent=4)
print(s)
# run optimization
best,best_eval=genetic_algorithm(objective_start, objective_collect,8,params["gen_niter"],params["gen_npop"],params["gen_rcross"],params["gen_rmut"],params["mins"],params["maxs"])
print(best)
print(best_eval)