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

import pypeliner

import citup

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pypeliner.easypypeliner.add_arguments(argparser)
    argparser.add_argument('input_freqs', help='Input Mutation Frequencies')
    argparser.add_argument('output_solution', help='Output Solution Tree')
    argparser.add_argument('output_all_trees', help='Output For All Trees')

    cfg = pypeliner.easypypeliner.Config(vars(argparser.parse_args()))
    pyp = pypeliner.easypypeliner.EasyPypeliner([citup], cfg)

    citup_bin_directory = os.path.join(os.path.dirname(citup.__file__))
    citup_iter_tool = os.path.join(citup_bin_directory, 'CitupIter.exe')

    lowmem = {'mem':1,'ncpu':1}

    pyp.sch.transform('create_trees', (), lowmem,
        citup.create_trees,
        pyp.sch.oobj('trees', ('tree',)),
        int(cfg.min_nodes),
        int(cfg.max_nodes),
        int(cfg.max_children_per_node))
    
    pyp.sch.commandline('run_citup_iter', ('tree',), lowmem, 
        citup_iter_tool,
        '-t', pyp.sch.iobj('trees', ('tree',)).prop('unlabeled_tree_string'),
        '-f', pyp.sch.input(cfg.input_freqs),
        '-s', pyp.sch.ofile('results', ('tree',)))

    pyp.sch.transform('select_optimal_tree', (), lowmem,
        citup.select_optimal_tree,
        None,
        pyp.sch.input(cfg.input_freqs),
        pyp.sch.iobj('trees', ('tree',)),
        pyp.sch.ifile('results', ('tree',)),
        pyp.sch.output(cfg.output_solution), 
        pyp.sch.output(cfg.output_all_trees))

    pyp.run()

