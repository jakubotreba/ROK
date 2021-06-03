#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 11:11:29 2021

@author: jakub
"""

import matplotlib.pyplot as plt
import csv
from tqdm import tqdm
import IsoSpecPy
import numpy as np

#sth = ['nucleons', 'C', 'H', 'O', 'N', 'S']

def makeplot(sth):
    
    x = []
    y = []
    z = []

    with open("../data/wasserstein_results_05.csv", 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=",")
        for row in tqdm(plots):
            #x.append(float(row[0]))
            #y.append(float(row[1]))
            conf1 = next(IsoSpecPy.IsoOrderedGenerator(row[2], get_confs = True).__iter__())[2]
            conf2 = next(IsoSpecPy.IsoOrderedGenerator(row[3], get_confs = True).__iter__())[2]
            #C H N O S
            if sth == 'nucleons':
                nukleons1 = [0,0,0]
                nukleons2 = [0,0,0]
                i = 0
                nukleons1[0]+=conf1[i][0]*12
                nukleons2[0]+=conf2[i][0]*12
                nukleons1[0]+=conf1[i][1]*13
                nukleons2[0]+=conf2[i][1]*13
                i = 1
                nukleons1[0]+=conf1[i][0]*1
                nukleons2[0]+=conf2[i][0]*1
                nukleons1[0]+=conf1[i][1]*2
                nukleons2[0]+=conf2[i][1]*2
                i = 2
                nukleons1[0]+=conf1[i][0]*14
                nukleons2[0]+=conf2[i][0]*14
                nukleons1[0]+=conf1[i][1]*15
                nukleons2[0]+=conf2[i][1]*15
                i = 3
                nukleons1[0]+=conf1[i][0]*16
                nukleons2[0]+=conf2[i][0]*16
                nukleons1[0]+=conf1[i][1]*17
                nukleons2[0]+=conf2[i][1]*17
                nukleons1[0]+=conf1[i][2]*18
                nukleons2[0]+=conf2[i][2]*18
                i = 4
                nukleons1[0]+=conf1[i][0]*32
                nukleons2[0]+=conf2[i][0]*32
                nukleons1[0]+=conf1[i][1]*33
                nukleons2[0]+=conf2[i][1]*33
                nukleons1[0]+=conf1[i][2]*34
                nukleons2[0]+=conf2[i][2]*34
                nukleons1[0]+=conf1[i][3]*36
                nukleons2[0]+=conf2[i][3]*36
                
                diff = abs(nukleons1[0] - nukleons2[0])
                if diff != 0:
                    z.append(diff)
                    x.append(float(row[0]))
                    y.append(float(row[1]))
                else:
                    continue

            elif sth == 'C':
                cs1 = 0
                cs2 = 0
                cs1 += conf1[0][0] + conf1[0][1]
                cs2 += conf2[0][0] + conf2[0][1]
                diff = abs(cs1 - cs2)
                if diff != 0:
                    z.append(diff)
                    x.append(float(row[0]))
                    y.append(float(row[1]))
                else:
                    continue

            elif sth == 'H':
                hs1 = 0
                hs2 = 0
                hs1 += conf1[1][0] + conf1[1][1]
                hs2 += conf2[1][0] + conf2[1][1]
                diff = abs(hs1 - hs2)
                if diff != 0:
                    z.append(diff)
                    x.append(float(row[0]))
                    y.append(float(row[1]))

            elif sth == 'N':
                ns1 = 0
                ns2 = 0
                ns1 += conf1[2][0] + conf1[2][1]
                ns2 += conf2[2][0] + conf2[2][1]
                diff = abs(ns1 - ns2)
                if diff != 0:
                    z.append(diff)
                    x.append(float(row[0]))
                    y.append(float(row[1]))

            elif sth == 'O':
                os1 = 0
                os2 = 0
                os1 += conf1[3][0] + conf1[3][1] + conf1[3][2]
                os2 += conf2[3][0] + conf2[3][1] + conf2[3][2]
                diff = abs(os1 - os2)
                if diff != 0:
                    z.append(diff)
                    x.append(float(row[0]))
                    y.append(float(row[1]))

            elif sth == 'S':
                ss1 = 0
                ss2 = 0
                ss1 += conf1[4][0] + conf1[4][1] + conf1[4][2] + conf1[4][3]
                ss2 += conf2[4][0] + conf2[4][1] + conf2[4][2] + conf2[4][3]
                diff = abs(ss1 - ss2)
                if diff != 0:
                    z.append(diff)
                    x.append(float(row[0]))
                    y.append(float(row[1]))

    fig, ax = plt.subplots()
    plot = ax.scatter(x, y, c = z, alpha = 0.3, edgecolors = None, cmap = 'YlGnBu', s = 5.0)
    legend = ax.legend(*plot.legend_elements(alpha = 1.0), loc = "lower right", bbox_to_anchor = (0.94, 0.0), fontsize = 'xx-small', ncol = 3, title=('%s' % (sth)) + ' ' + "number difference")
    ax.add_artist(legend)
    plt.xlabel('Mean mass difference')
    plt.ylabel('Wasserstein distance')
    plt.title('Mean mass difference, Wasserstein distance and' + ' ' + ('%s' % (sth)) + ' ' + 'difference plot ')
    #plt.gray()
    plt.legend()
    #plt.show()
    plt.tight_layout()
    plt.savefig("../plots/Plot" + ('%s' % (sth)) + ".png", dpi = 200)

makeplot('C')
makeplot('H')
makeplot('O')
makeplot('N')
makeplot('S')
makeplot('nucleons')





