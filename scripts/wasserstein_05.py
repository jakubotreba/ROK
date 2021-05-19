#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 12:58:09 2021

@author: jakub
"""

from random import *
from tqdm import tqdm
import linecache
import csv
import IsoSpecPy
import argparse

parser = argparse.ArgumentParser()
    
parser.add_argument("-K", help = 'Enter the sorted mass and formula file')

args = parser.parse_args()

def count_totalprob(mass, formula, prob):
    
    s = IsoSpecPy.IsoTotalProb(prob, formula)
    s.normalize()
    
    return (mass, formula, s)

def count_wasserstein(spect1, spect2):
    
    wasserstein_distance = spect1.wassersteinDistance(spect2)
    
    return wasserstein_distance

def cluster_getter(window):
    
    with open(args.K) as infile:
        with open('data/wasserstein_results_05.csv', 'w') as outfile:
            
            writer = csv.writer(outfile)
            considered = []
            spectres = {}
            
            
            for i, line in tqdm(enumerate(infile), total = 2001070):
                
                if i > 1 and i < 2001069:
                    
                    mass, formula = line.split()
                    mass = float(mass[:-1])
                    considered.append((mass, formula))
                    counter = i
                    
                    # move backward
                    
                    while True:
                        
                        line_ = linecache.getline(args.K, counter - 1)
                        mass_, formula_ = line_.split()
                        mass_ = float(mass_[:-1])
                        if mass_ + window >= mass and (mass_, formula_) not in considered:
                            considered.append((mass_, formula_))
                            counter -= 1
                        else:
                            break
                        
                    # move forward
                        
                    while True:
                        
                        line_ = linecache.getline(args.K, counter + 1)
                        mass_, formula_ = line_.split()
                        mass_ = float(mass_[:-1])
                        if mass_ - window <= mass and (mass_, formula_) not in considered:
                            considered.append((mass_, formula_))
                            counter += 1
                        else:
                            break
                        
                    if len(considered) > 1:
                        for record in considered:
                            isototal = count_totalprob(record[0], record[1], 0.99)
                            # mass | spectre | formula
                            spectres[record[1]] = (isototal[0], isototal[2], isototal[1])
                            
                        for j in range(1, len(considered)):
                            wass = count_wasserstein(spectres[considered[j][1]][1], spectres[considered[0][1]][1])
                            mass_difference = abs(spectres[considered[j][1]][0] - spectres[considered[0][1]][0])
                            formula1, formula2 = spectres[considered[j][1]][2], spectres[considered[0][1]][2]
                            writer.writerow([mass_difference, wass, formula1, formula2])
                    
                    considered = []
                    spectres = {}
                
                else:
                    continue
                    
        outfile.close()
    infile.close()
                
cluster_getter(0.5)


















