#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 12:06:27 2021

@author: jakub
"""

import argparse
import IsoSpecPy
from Bio import SeqIO


# Function to obtain chemical formulas
# from a fasta file

def fasta_to_formula():

	parser = argparse.ArgumentParser()
	parser.add_argument("-U", help = 'Uniprot fasta file')

	args = parser.parse_args()

	results = open('data/fasta_to_formula_results.txt',"a")

	with open(args.U, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			#C H N O S		
			formula = [0,0,0,0,0]
			for j in range(len(record.seq)):
				if record.seq[j]=='A':
					formula[0]+=3
					formula[1]+=7
					formula[2]+=1
					formula[3]+=2
				elif record.seq[j]=='C':
					formula[0]+=3
					formula[1]+=7
					formula[2]+=1
					formula[3]+=2
					formula[4]+=1
				elif record.seq[j]=='D':
					formula[0]+=4
					formula[1]+=7
					formula[2]+=1
					formula[3]+=4
				elif record.seq[j]=='E':
					formula[0]+=5
					formula[1]+=9
					formula[2]+=1
					formula[3]+=4
				elif record.seq[j]=='F':
					formula[0]+=9
					formula[1]+=11
					formula[2]+=1
					formula[3]+=2
				elif record.seq[j]=='G':
					formula[0]+=2
					formula[1]+=5
					formula[2]+=1
					formula[3]+=2
				elif record.seq[j]=='H':
					formula[0]+=6
					formula[1]+=9
					formula[2]+=3
					formula[3]+=2
				elif record.seq[j]=='I':
					formula[0]+=6
					formula[1]+=13
					formula[2]+=1
					formula[3]+=2
				elif record.seq[j]=='K':
					formula[0]+=6
					formula[1]+=14
					formula[2]+=2
					formula[3]+=2
				elif record.seq[j]=='L':
					formula[0]+=6
					formula[1]+=13
					formula[2]+=1
					formula[3]+=2
				elif record.seq[j]=='M':
					formula[0]+=5
					formula[1]+=11
					formula[2]+=1
					formula[3]+=2
					formula[4]+=1
				elif record.seq[j]=='N':
					formula[0]+=4
					formula[1]+=8
					formula[2]+=2
					formula[3]+=3
				elif record.seq[j]=='P':
					formula[0]+=5
					formula[1]+=9
					formula[2]+=1
					formula[3]+=2
				elif record.seq[j]=='Q':
					formula[0]+=5
					formula[1]+=10
					formula[2]+=2
					formula[3]+=3
				elif record.seq[j]=='R':
					formula[0]+=6
					formula[1]+=14
					formula[2]+=4
					formula[3]+=2
				elif record.seq[j]=='S':
					formula[0]+=3
					formula[1]+=7
					formula[2]+=1
					formula[3]+=3
				elif record.seq[j]=='T':
					formula[0]+=4
					formula[1]+=9
					formula[2]+=1
					formula[3]+=3
				elif record.seq[j]=='V':
					formula[0]+=5
					formula[1]+=11
					formula[2]+=1
					formula[3]+=2
				elif record.seq[j]=='W':
					formula[0]+=11
					formula[1]+=12
					formula[2]+=2
					formula[3]+=2
				elif record.seq[j]=='Y':
					formula[0]+=9
					formula[1]+=11
					formula[2]+=1
					formula[3]+=3

			formula[0]='C'+str(formula[0])
			formula[1]='H'+str(formula[1])
			formula[2]='N'+str(formula[2])
			formula[3]='O'+str(formula[3])
			formula[4]='S'+str(formula[4])
			formula=''.join(formula)

			results.write(">"+str(record.id)+"\n"+formula+'\n')

	results.close()

fasta_to_formula()


























