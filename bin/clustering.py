
from time import time
import numpy as np
import pylab as pl
import argparse

from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn.datasets import load_digits

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
argparser.add_argument('-i', help='Input file (freuquencies of mutations to be clustered)')
argparser.add_argument('-o', help='Output file location (Cluster centers go to this file)')
argparser.add_argument('-c', help='Number of clusters', type=int)

argparser.parse_args();


arguments = argparser.parse_args()
num_clusters = arguments.c
input_file_location = arguments.i
output_file_location = arguments.o

file = open(input_file_location,'r')
for i in range(3):
    words=file.readline().split()
    if   i==0:
        num_mutations=int(words[1])
    elif i==1:
        num_samples=int(words[1])
    elif i==2:
        error_rate=float(words[1])


f = [[0 for i in range(num_samples)] for j in range(num_mutations)]

for i in range(num_mutations):

    line=''
    while len(line.split())==0:
        line=file.readline()

    words=line.split()
    if len(words):
        for s in range(num_samples):
            f[i][s]=float(words[s])


data = KMeans(init='random', n_clusters=num_clusters, n_init=1000).fit(f)
cluster_centers = data.cluster_centers_
labels  = data.labels_
inertia = data.inertia_

output = open(output_file_location, 'w')
output.write("Num_mutations " + num_mutations.__str__() + " \n")
output.write("Num_clusters " + num_clusters.__str__() + " \n")
output.write("Num_samples " + num_samples.__str__() + " \n")
output.write("Inertia " + inertia.__str__() + " \n")
output.write("\n \n")

for cluster in cluster_centers:
    for s in cluster:
        output.write('{:10.5f}'.format(s) + ' ')
    output.write("\n")

output.write("\n")

for s in labels:
    output.write('{:6}'.format(s) + ' ')

output.write("\n\n\n")

for mutation_i in f:
	for sample_j in mutation_i:
		output.write('{:6}'.format(sample_j) + ' ')
	output.write("\n")

