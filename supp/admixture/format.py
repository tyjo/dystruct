#!/usr/bin/env python

# Converts sample files to .geno format for ADMIXTURE

import os
import numpy as np
from sys import argv

folder_name = argv[1][8:]
os.makedirs("formatted_data/" + folder_name)
outfile = open("formatted_data/" + folder_name + "/samples.ped", "w")

genotypes = np.loadtxt("../data/" + folder_name + "/samples", delimiter="\t")[1:,]

for row in genotypes.T:
    outfile.write("FAM001 0 0 0 0 0 ")
    genotypes = []
    for g in row:
        if g == 0:
            genotypes.append(1)
            genotypes.append(1)
        elif g == 1:
            genotypes.append(1)
            genotypes.append(2)
        else:
            genotypes.append(2)
            genotypes.append(2)
    outfile.write(" ".join([str(g) for g in genotypes]))
    outfile.write("\n")

