
import csv
import sys
import logging
import os
import glob
import ConfigParser
import itertools
import argparse
import string
import ConfigParser
import shutil
from collections import *
import numpy as np
import pandas as pd
from sklearn import mixture
import pypeliner

def read_frequencies(freq_filename):
    with open(freq_filename, 'r') as f:
        metadata = dict()
        data = list()
        for line in f:
            if ':' in line:
                entry = line.rstrip().split()
                metadata[entry[0]] = entry[1]
            else:
                for entry in line.rstrip().split():
                    data.append(entry)
        num_mutations = int(metadata['Num_mutations:'])
        num_samples = int(metadata['Num_samples:'])
        data = list(reversed(data))
        freq = np.zeros([num_mutations, num_samples])
        for s in range(num_samples):
            for i in range(num_mutations):
                freq[i,s] = float(data.pop())
        return freq

def read_results(results_filename):
    with open(results_filename, 'r') as f:
        metadata = dict()
        for line in f:
            entry = line.rstrip().split()
            if len(entry) == 2:
                metadata[entry[0]] = entry[1]
        results_info = dict()
        results_info['objective_value'] = float(metadata['Objective_value'])
        results_info['cplex_status'] = metadata['Cplex_status']
        results_info['cplex_time'] = float(metadata['Cplex_hours'])
        results_info['num_intermediate_clusters'] = int(metadata['Num_clusters'])
        results_info['tree_index'] = int(metadata['Tree_index'])
        results_info['num_nodes'] = int(metadata['Num_nodes'])
        results_info['num_samples'] = int(metadata['Num_samples'])
        results_info['num_mutations'] = int(metadata['Num_mutations'])
        return results_info
    
def estimate_error_rate(freq_filename):
    freq = read_frequencies(freq_filename)
    gmm_bics = list()
    for n in range(1, min(21, freq.shape[0]+1)):
        gmm = mixture.GMM(n_components=n, covariance_type='spherical')
        gmm.fit(freq)
        gmm_bics.append((gmm.bic(freq), n))
    n = sorted(gmm_bics, key=lambda a: a[0])[0][1]
    gmm = mixture.GMM(n_components=n, covariance_type='spherical')
    gmm.fit(freq)
    return np.sqrt(np.mean(gmm.covars_))
    
def compute_bic(error_rate, objective_value, num_samples, num_nodes, num_mutations):
    likelihood = objective_value / (2.0 * error_rate * error_rate)
    bic = likelihood + num_samples * (num_nodes - 1.0) * np.log(num_mutations)
    return bic

 
if __name__ == '__main__':
        
    options = sys.argv[1:]
    working_directory = options[0]
    freq_filename = working_directory + "/Frequencies.txt"
    error_rate = estimate_error_rate(freq_filename)
    best_score = -1.0
    stats = dict()
    matching_pattern = working_directory + "/Results_*_*.txt"
    files = glob.glob(matching_pattern)
    for results_filename in files: 
        res_info = read_results(results_filename)
        bic = compute_bic(error_rate, res_info['objective_value'], res_info['num_samples'], res_info['num_nodes'], res_info['num_mutations'])
        if best_score < 0 or bic < best_score:
            best_score = bic
        stats[results_filename] = bic

    print "Best tree(s) is: "  
    for results_filename in files:
        if stats[results_filename] == best_score:
            print results_filename, best_score
