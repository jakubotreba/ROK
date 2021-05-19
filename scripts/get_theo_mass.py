#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 12:32:24 2021

@author: jakub
"""


import IsoSpecPy
from tqdm import tqdm
import argparse


parser = argparse.ArgumentParser()
    
parser.add_argument("-I", help = 'Enter the formulas file')

args = parser.parse_args()


# Function to get theoretical average mass of each formula
def get_theo_mass():
    
    list_of_mass_formulas = []
    
    # Set of visited proteins
    seen = set()
    
    with open(args.I) as file:
        for line in tqdm(file, total = 2228232):
            if ">" in line:
                continue
            elif line in seen:
                continue
            
            seen.add(line)
            
            formula = line.strip()
            mass = IsoSpecPy.Iso(formula).getTheoreticalAverageMass()
            list_of_mass_formulas.append((mass, formula))

    list_of_mass_formulas.sort()
    
    with open('data/sorted_formulas_masses.txt', 'w') as outfile:
        for item in list_of_mass_formulas:
            outfile.write(str(item[0]) + ', ' + str(item[1]) + '\n')
    outfile.close()
    
get_theo_mass()



















