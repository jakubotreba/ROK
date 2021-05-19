# ROK - Protein recognizability project
Project made during Architecture of Large Bioinformatics Projects class at University of Warsaw.

# Members:
Jakub Otręba, Ewa Kizling, Julia Rymuza
# Description:
We want to use IsoSpecPy, which is a fine structure isotopic calculator. 
It can reveal for you a given fraction of the isotopic distribution of mass for a given molecule 
based only on its molecular composition. Mass is not unique due to the presence of isotopes that occur 
randomly in Nature, albeit with well studied frequencies. We are going to use this calculator
to model protein spectra from fasta files downloaded from UniProt for Homo Sapiens (circa 50 000 sequences), 
such as in mass spectrometry, and then try to compare
protein recognizability by a classical approach, which is basicaly just looking at
a mass differences with an approach using Wasserstein metric.
