#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 16:17:57 2021

@author: jakub
"""

from tqdm import tqdm
from random import *
import IsoSpecPy
import IsoSpecPy.Distributions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_proteins(file):
   
    with open(file) as f:
        
        proteins_in_clusters = {}
        
        for line in f:
            protein1, protein2 = line.split(",")[2:4]
            cluster_n = line.strip().split(",")[-1]
            if cluster_n in proteins_in_clusters.keys():
                proteins_in_clusters[cluster_n].add(protein1)
                proteins_in_clusters[cluster_n].add(protein2)
            else:
                proteins_in_clusters[cluster_n] = set([protein1, protein2])
    f.close()
    
    clusters = []
    
    for cluster in proteins_in_clusters.values():
        proteins_in_this_cluster = []
        for protein in cluster:
            formula = protein
            mass = mass = IsoSpecPy.Iso(formula).getTheoreticalAverageMass()
            proteins_in_this_cluster.append((mass, formula))
        clusters.append(proteins_in_this_cluster)
    return clusters

# Gaussian noise

def gaussian(clusters, stdev, cluster_n):
    
    with open('../data/' + str(cluster_n) + "gaussians" + str(stdev) + ".txt", "w+") as outfile:
        
        noise = IsoSpecPy.Distributions.Gaussian(stdev, 0.001, 0.99)
        
        for cluster in [clusters[cluster_n]]:
            print(cluster_n, "gaussian")
            the_protein = cluster[int(len(cluster) / 2 )]
            the_protein_formula = the_protein[1]
            the_protein_mass = the_protein[0]
            the_protein_spectre = IsoSpecPy.IsoTotalProb(0.99, the_protein_formula)
            spectre_binned = the_protein_spectre.binned(0.01)
            spectre_binned_unsettled = spectre_binned * noise
            outfile.write("<" + str(the_protein_spectre.wassersteinDistance(spectre_binned_unsettled)) + ", " + the_protein_formula + ", " + the_protein_formula + "\n")
            for protein in cluster:
                if protein == the_protein:
                    continue
                else:
                    protein_formula = protein[1]
                    protein_mass = protein[0]
                    protein_spectre = IsoSpecPy.IsoTotalProb(0.99, protein_formula)
                    outfile.write(str(protein_spectre.wassersteinDistance(spectre_binned_unsettled)) + ", " + protein_formula + ", " + the_protein_formula + "\n")

klastry = get_proteins('../data/clusters.txt')

# 0
gaussian(klastry, 1e-06, 0)
gaussian(klastry, 1.778279410038923e-06, 0)
gaussian(klastry, 3.162277660168379e-06, 0)
gaussian(klastry, 5.623413251903491e-06, 0)
gaussian(klastry, 1e-05, 0)
gaussian(klastry, 1.778279410038923e-05, 0)
gaussian(klastry, 3.1622776601683795e-05, 0)
gaussian(klastry, 5.623413251903491e-05, 0)
gaussian(klastry, 0.0001, 0)
gaussian(klastry, 0.00017782794100389227, 0)
gaussian(klastry, 0.00031622776601683794, 0)
gaussian(klastry, 0.0005623413251903491, 0)
gaussian(klastry, 0.001, 0)
gaussian(klastry, 0.0017782794100389228, 0)
gaussian(klastry, 0.0031622776601683794, 0)
gaussian(klastry, 0.005623413251903491, 0)
gaussian(klastry, 0.01, 0)
gaussian(klastry, 0.01778279410038923, 0)
gaussian(klastry, 0.03162277660168379, 0)
gaussian(klastry, 0.05623413251903491, 0)
gaussian(klastry, 0.1, 0)
gaussian(klastry, 0.1778279410038923, 0)
gaussian(klastry, 0.31622776601683794, 0)
gaussian(klastry, 0.5623413251903491, 0)
gaussian(klastry, 1.0, 0)
gaussian(klastry, 1.7782794100389228, 0)
gaussian(klastry, 3.1622776601683795, 0)
gaussian(klastry, 5.623413251903491, 0)
gaussian(klastry, 10.0, 0)

# 1
gaussian(klastry, 1e-06, 1)
gaussian(klastry, 1.778279410038923e-06, 1)
gaussian(klastry, 3.162277660168379e-06, 1)
gaussian(klastry, 5.623413251903491e-06, 1)
gaussian(klastry, 1e-05, 1)
gaussian(klastry, 1.778279410038923e-05, 1)
gaussian(klastry, 3.1622776601683795e-05, 1)
gaussian(klastry, 5.623413251903491e-05, 1)
gaussian(klastry, 0.0001, 1)
gaussian(klastry, 0.00017782794100389227, 1)
gaussian(klastry, 0.00031622776601683794, 1)
gaussian(klastry, 0.0005623413251903491, 1)
gaussian(klastry, 0.001, 1)
gaussian(klastry, 0.0017782794100389228, 1)
gaussian(klastry, 0.0031622776601683794, 1)
gaussian(klastry, 0.005623413251903491, 1)
gaussian(klastry, 0.01, 1)
gaussian(klastry, 0.01778279410038923, 1)
gaussian(klastry, 0.03162277660168379, 1)
gaussian(klastry, 0.05623413251903491, 1)
gaussian(klastry, 0.1, 1)
gaussian(klastry, 0.1778279410038923, 1)
gaussian(klastry, 0.31622776601683794, 1)
gaussian(klastry, 0.5623413251903491, 1)
gaussian(klastry, 1.0, 1)
gaussian(klastry, 1.7782794100389228, 1)
gaussian(klastry, 3.1622776601683795, 1)
gaussian(klastry, 5.623413251903491, 1)
gaussian(klastry, 10.0, 1)

# 2
gaussian(klastry, 1e-06, 2)
gaussian(klastry, 1.778279410038923e-06, 2)
gaussian(klastry, 3.162277660168379e-06, 2)
gaussian(klastry, 5.623413251903491e-06, 2)
gaussian(klastry, 1e-05, 2)
gaussian(klastry, 1.778279410038923e-05, 2)
gaussian(klastry, 3.1622776601683795e-05, 2)
gaussian(klastry, 5.623413251903491e-05, 2)
gaussian(klastry, 0.0001, 2)
gaussian(klastry, 0.00017782794100389227, 2)
gaussian(klastry, 0.00031622776601683794, 2)
gaussian(klastry, 0.0005623413251903491, 2)
gaussian(klastry, 0.001, 2)
gaussian(klastry, 0.0017782794100389228, 2)
gaussian(klastry, 0.0031622776601683794, 2)
gaussian(klastry, 0.005623413251903491, 2)
gaussian(klastry, 0.01, 2)
gaussian(klastry, 0.01778279410038923, 2)
gaussian(klastry, 0.03162277660168379, 2)
gaussian(klastry, 0.05623413251903491, 2)
gaussian(klastry, 0.1, 2)
gaussian(klastry, 0.1778279410038923, 2)
gaussian(klastry, 0.31622776601683794, 2)
gaussian(klastry, 0.5623413251903491, 2)
gaussian(klastry, 1.0, 2)
gaussian(klastry, 1.7782794100389228, 2)
gaussian(klastry, 3.1622776601683795, 2)
gaussian(klastry, 5.623413251903491, 2)
gaussian(klastry, 10.0, 2)

# 3
gaussian(klastry, 1e-06, 3)
gaussian(klastry, 1.778279410038923e-06, 3)
gaussian(klastry, 3.162277660168379e-06, 3)
gaussian(klastry, 5.623413251903491e-06, 3)
gaussian(klastry, 1e-05, 3)
gaussian(klastry, 1.778279410038923e-05, 3)
gaussian(klastry, 3.1622776601683795e-05, 3)
gaussian(klastry, 5.623413251903491e-05, 3)
gaussian(klastry, 0.0001, 3)
gaussian(klastry, 0.00017782794100389227, 3)
gaussian(klastry, 0.00031622776601683794, 3)
gaussian(klastry, 0.0005623413251903491, 3)
gaussian(klastry, 0.001, 3)
gaussian(klastry, 0.0017782794100389228, 3)
gaussian(klastry, 0.0031622776601683794, 3)
gaussian(klastry, 0.005623413251903491, 3)
gaussian(klastry, 0.01, 3)
gaussian(klastry, 0.01778279410038923, 3)
gaussian(klastry, 0.03162277660168379, 3)
gaussian(klastry, 0.05623413251903491, 3)
gaussian(klastry, 0.1, 3)
gaussian(klastry, 0.1778279410038923, 3)
gaussian(klastry, 0.31622776601683794, 3)
gaussian(klastry, 0.5623413251903491, 3)
gaussian(klastry, 1.0, 3)
gaussian(klastry, 1.7782794100389228, 3)
gaussian(klastry, 3.1622776601683795, 3)
gaussian(klastry, 5.623413251903491, 3)
gaussian(klastry, 10.0, 3)

# 4
gaussian(klastry, 1e-06, 4)
gaussian(klastry, 1.778279410038923e-06, 4)
gaussian(klastry, 3.162277660168379e-06, 4)
gaussian(klastry, 5.623413251903491e-06, 4)
gaussian(klastry, 1e-05, 4)
gaussian(klastry, 1.778279410038923e-05, 4)
gaussian(klastry, 3.1622776601683795e-05, 4)
gaussian(klastry, 5.623413251903491e-05, 4)
gaussian(klastry, 0.0001, 4)
gaussian(klastry, 0.00017782794100389227, 4)
gaussian(klastry, 0.00031622776601683794, 4)
gaussian(klastry, 0.0005623413251903491, 4)
gaussian(klastry, 0.001, 4)
gaussian(klastry, 0.0017782794100389228, 4)
gaussian(klastry, 0.0031622776601683794, 4)
gaussian(klastry, 0.005623413251903491, 4)
gaussian(klastry, 0.01, 4)
gaussian(klastry, 0.01778279410038923, 4)
gaussian(klastry, 0.03162277660168379, 4)
gaussian(klastry, 0.05623413251903491, 4)
gaussian(klastry, 0.1, 4)
gaussian(klastry, 0.1778279410038923, 4)
gaussian(klastry, 0.31622776601683794, 4)
gaussian(klastry, 0.5623413251903491, 4)
gaussian(klastry, 1.0, 4)
gaussian(klastry, 1.7782794100389228, 4)
gaussian(klastry, 3.1622776601683795, 4)
gaussian(klastry, 5.623413251903491, 4)
gaussian(klastry, 10.0, 4)

# 5
gaussian(klastry, 1e-06, 5)
gaussian(klastry, 1.778279410038923e-06, 5)
gaussian(klastry, 3.162277660168379e-06, 5)
gaussian(klastry, 5.623413251903491e-06, 5)
gaussian(klastry, 1e-05, 5)
gaussian(klastry, 1.778279410038923e-05, 5)
gaussian(klastry, 3.1622776601683795e-05, 5)
gaussian(klastry, 5.623413251903491e-05, 5)
gaussian(klastry, 0.0001, 5)
gaussian(klastry, 0.00017782794100389227, 5)
gaussian(klastry, 0.00031622776601683794, 5)
gaussian(klastry, 0.0005623413251903491, 5)
gaussian(klastry, 0.001, 5)
gaussian(klastry, 0.0017782794100389228, 5)
gaussian(klastry, 0.0031622776601683794, 5)
gaussian(klastry, 0.005623413251903491, 5)
gaussian(klastry, 0.01, 5)
gaussian(klastry, 0.01778279410038923, 5)
gaussian(klastry, 0.03162277660168379, 5)
gaussian(klastry, 0.05623413251903491, 5)
gaussian(klastry, 0.1, 5)
gaussian(klastry, 0.1778279410038923, 5)
gaussian(klastry, 0.31622776601683794, 5)
gaussian(klastry, 0.5623413251903491, 5)
gaussian(klastry, 1.0, 5)
gaussian(klastry, 1.7782794100389228, 5)
gaussian(klastry, 3.1622776601683795, 5)
gaussian(klastry, 5.623413251903491, 5)
gaussian(klastry, 10.0, 5)

# Gather the data and plot it

def get_gaussian_data(files):
	
	stdevs = []
	ys = []
	
	for file in files:
		stdev0 = file[10:]
		stdev = stdev0.split("txt")[0].strip(".")
		stdevs.append(stdev)

	nr_klastra = files[0][0]

	for i in range(len(stdevs)):
		if "-" in stdevs[i]:
			if "." not in stdevs[i]:
				stdevs[i] = stdevs[i][0] + ".0" + stdevs[i][1:]
			number = stdevs[i]
			lead, power = number.split("e-")
			a, b = lead.split(".")
			number = "0."+"0"*(int(power) - 1) + a + b
			stdevs[i] = float(number)
		else:
			stdevs[i] = float(stdevs[i])

	for i in range(len(files)):
		with open("../data/"+files[i]) as f:
			y = []
			for line in f:
				if "<" in line:
					y.append(float(line.split(",")[0].strip("<")))
				else:
					y.append(float(line.split(",")[0]))
		
		f.close()
		ys.append(y)

	df = pd.DataFrame(ys, index = stdevs)
	df_main = df[0]
	df = df.transpose()
	df = df.sort_values(by = stdevs[1])

	fig, ax = plt.subplots()
	
	plt.plot(df.transpose()[int(0.1*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.2*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.3*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.4*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.5*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.6*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.7*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.8*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.9*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(1.0*len(df)) - 1], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df_main.transpose(), color = 'red', linewidth=0.75)
	plt.xscale('log')
	plt.title(f'Cluster no {nr_klastra} - Gaussian noise')
	plt.xlabel('Standard deviation (Da)')
	plt.ylabel('Wasserstein distance')
	plt.savefig('../plots/' + f'gauss{nr_klastra}.png', dpi = 200)

# Plot the gaussians
for i in range(0, 6):
	get_gaussian_data([f"{i}gaussians1e-06.txt", f"{i}gaussians1.778279410038923e-06.txt", f"{i}gaussians3.162277660168379e-06.txt", f"{i}gaussians5.623413251903491e-06.txt", f"{i}gaussians1e-05.txt", f"{i}gaussians1.778279410038923e-05.txt", f"{i}gaussians3.1622776601683795e-05.txt", f"{i}gaussians5.623413251903491e-05.txt", f"{i}gaussians0.0001.txt", f"{i}gaussians0.00017782794100389227.txt", f"{i}gaussians0.00031622776601683794.txt", f"{i}gaussians0.0005623413251903491.txt", f"{i}gaussians0.001.txt", f"{i}gaussians0.0017782794100389228.txt", f"{i}gaussians0.0031622776601683794.txt", f"{i}gaussians0.005623413251903491.txt", f"{i}gaussians0.01.txt", f"{i}gaussians0.01778279410038923.txt", f"{i}gaussians0.03162277660168379.txt", f"{i}gaussians0.05623413251903491.txt", f"{i}gaussians0.1.txt", f"{i}gaussians0.1778279410038923.txt", f"{i}gaussians0.31622776601683794.txt", f"{i}gaussians0.5623413251903491.txt", f"{i}gaussians1.0.txt", f"{i}gaussians1.7782794100389228.txt", f"{i}gaussians3.1622776601683795.txt", f"{i}gaussians5.623413251903491.txt", f"{i}gaussians10.0.txt"])	

# Electronic noise

def electronic(clusters, signal, noise, cluster_n, iteration):

	with open('../data/' + str(iteration) + "electronic" + str(signal) + ".txt", "w+") as outfile:
		for cluster in [clusters[cluster_n]]:
			the_protein = cluster[int(len(cluster) / 2)]
			the_protein_formula = the_protein[1]
			print((cluster_n), "electronic")
			the_protein_spectre = IsoSpecPy.IsoTotalProb(0.99, the_protein_formula)
			the_protein_spectre_binned = the_protein_spectre.binned(0.1)
			M = []
			for mass, prob in the_protein_spectre_binned:
				M.append((mass, prob))
			M.sort()
			min_mass = M[0][0]
			max_mass = M[-1][0]
			masses = list(np.arange(min_mass - 0.1, max_mass + 0.1, 0.1))
			probs = [np.random.poisson(100.0) for _ in masses]
			S = IsoSpecPy.IsoDistribution(masses = masses, probs = probs)
			S.normalize()
			the_protein_spectre_unsettled = IsoSpecPy.IsoDistribution.LinearCombination([the_protein_spectre, S], [signal, noise])
			the_protein_spectre_unsettled_binned = the_protein_spectre_unsettled.binned(0.1)
			the_protein_spectre_unsettled_binned.normalize()
			the_protein_spectre.normalize()
			outfile.write("<" + str(the_protein_spectre.wassersteinDistance(the_protein_spectre_unsettled_binned)) + ", " + the_protein_formula + ", " + the_protein_formula + "\n")
			for protein in cluster:
				if protein == the_protein:
					continue
				protein_formula = protein[1]
				spectre = IsoSpecPy.IsoTotalProb(0.99, protein_formula)
				spectre.normalize()
				outfile.write(str(spectre.wassersteinDistance(the_protein_spectre_unsettled_binned)) + ", " + protein_formula + "; " + the_protein_formula + "\n")
			break        
	outfile.close()

for i in range(6):
    electronic(klastry, 1.0, 0.0, i, i +1)
    electronic(klastry, 0.95, 0.05, i, i +1)
    electronic(klastry, 0.9, 0.1, i, i +1)
    electronic(klastry, 0.85, 0.15, i, i +1)
    electronic(klastry, 0.80, 0.20, i, i +1)
    electronic(klastry, 0.75, 0.25, i, i +1)
    electronic(klastry, 0.70, 0.30, i, i +1)
    electronic(klastry, 0.65, 0.35, i, i +1)
    electronic(klastry, 0.60, 0.40, i, i +1)
    electronic(klastry, 0.55, 0.45, i, i +1)
    electronic(klastry, 0.50, 0.50, i, i +1)


def get_electronic_data(files):
	
	noises = []
	ys = []

	for file in files:
		nois0 = file[11:]
		nois = nois0.split("txt")[0].strip(".")
		noises.append(1.0 - float(nois))

	nr_klastra = files[0][0]

	for i in range(len(files)):
		with open('../data/'+files[i]) as f:
			y = []
			for line in f:
				if "<" in line:
					y.append(float(line.split(",")[0].strip("<")))
				else:
					y.append(float(line.split(",")[0]))
		
		f.close()
		ys.append(y)

	df = pd.DataFrame(ys, index = noises)
	df_main = df[0]
	df = df.transpose()
	df = df.sort_values(by = noises[1])

	fig, ax = plt.subplots()

	plt.plot(df.transpose()[int(0.1*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.2*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.3*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.4*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.5*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.6*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.7*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.8*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.9*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(1.0*len(df)) - 1], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df_main.transpose(), color = 'red', linewidth=0.75)
	plt.title(f'Cluster no {nr_klastra} - electronic noise')
	plt.xlabel('Noise level')
	plt.ylabel('Wasserstein distance')
	plt.savefig('../plots/' + f'electronic{nr_klastra}.png', dpi = 200)

#electronic
for i in range(1, 7):
	get_electronic_data([f"{i}electronic1.0.txt" ,f"{i}electronic0.95.txt" ,f"{i}electronic0.9.txt", f"{i}electronic0.85.txt", f"{i}electronic0.8.txt", f"{i}electronic0.75.txt", f"{i}electronic0.7.txt", f"{i}electronic0.65.txt", f"{i}electronic0.6.txt", f"{i}electronic0.55.txt", f"{i}electronic0.5.txt"])










