import csv
import sys
import logging
import os
import ConfigParser
import itertools
import argparse
import string
import gzip
import ConfigParser
import shutil
from collections import *
import numpy as np
import pandas as pd
from sklearn import mixture

import pypeliner

import citup
import utils.BeyerHedetmieni
import utils.treenode

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pypeliner.easypypeliner.add_arguments(argparser)
    argparser.add_argument('input_freqs', help='Input Mutation Frequencies')
    argparser.add_argument('output_solution', help='Output Solution Tree')
    argparser.add_argument('output_all_trees', help='Output For All Trees')

    cfg = pypeliner.easypypeliner.Config(vars(argparser.parse_args()))
    pyp = pypeliner.easypypeliner.EasyPypeliner([citup], cfg)

    citup_bin_directory = os.path.join(os.path.dirname(citup.__file__))
    clustering_script = os.path.join(citup_bin_directory, 'clustering.py')
    citup_cplex_tool = os.path.join(citup_bin_directory, 'CplexApplication.exe')
    gamma_directory = os.path.join(citup_bin_directory, os.pardir, 'GammaAdjMatrices')

    lowmem = {'mem':1,'ncpu':1}
    mltcpu = {'mem':8,'ncpu':cfg.cplex_threads}
    
    pyp.sch.transform('create_trees', (), lowmem,
        citup.create_trees,
        pyp.sch.oobj('trees', ('tree',)),
        int(cfg.min_nodes),
        int(cfg.max_nodes),
        int(cfg.max_children_per_node))
    
    pyp.sch.commandline('cluster_freq', ('tree',), lowmem, 
        sys.executable, clustering_script,
        '-i', pyp.sch.input(cfg.input_freqs),
        '-o', pyp.sch.ofile('freqclusters', ('tree',)),
        '-c', pyp.sch.iobj('trees', ('tree',)).prop('num_nodes'))

    pyp.sch.commandline('run_citup', ('tree',), mltcpu, 
        citup_cplex_tool, 
        '-o', pyp.sch.ofile('results', ('tree',)),
        '-i', pyp.sch.ifile('freqclusters', ('tree',)),
        '-g', gamma_directory,
        '-r', cfg.cplex_threads,
        '-d', cfg.cplex_time_limit,
        '-t', pyp.sch.iobj('trees', ('tree',)).prop('tree_index'),
        '-n', pyp.sch.iobj('trees', ('tree',)).prop('num_nodes'),
        '-q')

    pyp.sch.transform('select_optimal_tree', (), lowmem, 
        citup.select_optimal_tree,
        None,
        pyp.sch.input(cfg.input_freqs),
        pyp.sch.iobj('trees', ('tree',)),
        pyp.sch.ifile('results', ('tree',)),
        pyp.sch.output(cfg.output_solution), 
        pyp.sch.output(cfg.output_all_trees))

    pyp.run()

else:

    class TreeInfo(object):
        def __init__(self, num_nodes, tree_index, tree):
            self.num_nodes = num_nodes
            self.tree_index = tree_index
            self.tree = tree
        @property
        def unlabeled_tree_string(self):
            return self.tree.create_unlabled_tree_string()

    def generate_trees(min_nodes, max_nodes, max_children_per_node):
        for num_nodes in xrange(min_nodes, max_nodes + 1):
            for tree_index, parent_array in enumerate(utils.BeyerHedetmieni.getParentArrays(num_nodes, max_children_per_node)):
                yield TreeInfo(num_nodes, tree_index, utils.treenode.create_from_parent_array(parent_array))
                tree_string = utils.treenode.create_from_parent_array(parent_array).create_unlabled_tree_string()

    def create_trees(min_nodes, max_nodes, max_children_per_node):
        return dict(enumerate(generate_trees(min_nodes, max_nodes, max_children_per_node)))

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

    def read_results(tree, results_filename):
        with open(results_filename, 'r') as f:
            results_info = dict()
            results_info['num_nodes'] = tree.num_nodes
            results_info['tree_index'] = tree.tree_index
            results_info['tree_string'] = tree.unlabeled_tree_string
            metadata = dict()
            for line in f:
                entry = line.rstrip().split()
                if len(entry) == 2:
                    metadata[entry[0].lower()] = entry[1]
            results_info['objective_value'] = float(metadata['objective_value'])
            results_info['num_nodes'] = int(metadata['num_nodes'])
            results_info['num_samples'] = int(metadata['num_samples'])
            results_info['num_mutations'] = int(metadata['num_mutations'])
            try:
                results_info['cplex_status'] = metadata['cplex_status']
                results_info['cplex_time'] = float(metadata['cplex_hours'])
            except KeyError:
                pass
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
    
    def select_optimal_tree(freq_filename, trees, results_filenames, optimal_filename, all_trees_filename):
        error_rate = estimate_error_rate(freq_filename)
        results_table = list()
        for tree_id in trees.keys():
            results_table.append(pd.DataFrame(read_results(trees[tree_id], results_filenames[tree_id]), index=[tree_id]))
        results_table = pd.concat(results_table)
        results_table['error_rate'] = error_rate
        results_table['likelihood'] = results_table['objective_value'] / (2.0 * error_rate * error_rate)
        results_table['bic'] = 2.0 * results_table['likelihood'] + results_table['num_samples'] * (results_table['num_nodes'] - 1.0) * np.log(results_table['num_mutations'])
        results_table.sort('bic', inplace=True)
        results_table['optimal'] = False
        results_table['optimal'].iloc[0] = True
        results_table.to_csv(all_trees_filename, sep='\t', index=False)
        optimal_index = results_table.index[0]
        shutil.copyfile(results_filenames[optimal_index], optimal_filename)
    
