#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
import common
import yaml

# import matplotlib.pyplot as plt
# from scipy import stats
# from sklearn.decomposition import PCA
# import seaborn as sns
# from scipy import stats
# from glob import glob
# import os
# from itertools import product
# import matplotlib.cm as cm
# from matplotlib.patches import Ellipse
# import matplotlib.patheffects as path_effects
# import rpy2

np.random.seed(0)

parser = ArgumentParser(description='Calculates standard and FL-based score for diagnosed sample'
                                    'and predicts diagnosis, i.e. positive, false positive, uninformative, negative')

parser.add_argument('counts', help='TSV file with the corrected number or reads '
                                   'for each assumed fragment length (50-220, organised in columns) '
                                   'and autosomes (chr1..chr22, organised in rows)')

parser.add_argument('param_file', help='YAML file with trained parameters (output of the train.py script)')

args = parser.parse_args()

with open(args.param_file) as param_f:
    params = yaml.load(param_f)

counts = common.load_counts(args.counts)
counts_per_chromosome = counts.sum(axis=1)


def calc_zscore(diagnosed_chromosome):

    global params
    global counts
    global counts_per_chromosome

    mapped_reads = counts_per_chromosome.sum()
    refs = counts_per_chromosome[common.REFERENCE_CHROMOSOMES[diagnosed_chromosome]].sum() / mapped_reads
    tris = counts_per_chromosome[diagnosed_chromosome] / mapped_reads
    ratio = tris / refs

    m = params[common.METHOD_NCV][diagnosed_chromosome][common.ATTR_MEAN]
    s = params[common.METHOD_NCV][diagnosed_chromosome][common.ATTR_STD]

    return (ratio - m) / s


for chromosome in common.DIAGNOSED_CHROMOSOMES:
    print(chromosome, calc_zscore(chromosome))

