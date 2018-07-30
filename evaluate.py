#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
import common

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

args = parser.parse_args()


def calc_zscore(chr_tris, chromosome_counts):

    chr_refs, p1, p2, c = {
        12: ([2, 3, 6], 0.036426625, 0.200130699, 0.000584512),
        17: ([3, 6, 7, 8, 9, 15], 0.02915150251875600000, 0.29958905234496900000, 0.00026512091729272300),
        20: ([0, 3, 7, 9, 18, 19, 21], 0.012822216, 0.311503396, 0.000176532)
    }[chr_tris]

    mapped_reads = chromosome_counts.sum()

    refs = chromosome_counts[chr_refs].sum() / mapped_reads
    tris = chromosome_counts[chr_tris] / mapped_reads
    ps = tris / refs
    m = p1 / p2

    # Since we do not have the set of train samples to estimate standard deviation directly,
    # we use for testing purpose its interpolation,
    # see https://link.springer.com/article/10.1186/s40488-018-0083-x
    ss = np.sqrt((1 / mapped_reads) * ((p1 / p2) ** 2) * (1 / p1 + 1 / p2) + c ** 2)

    return (ps - m) / ss


counts = common.load_counts(args.counts)
counts_per_chromosome = counts.sum(axis=1)
for diagnosed_chromosome in common.DIAGNOSED_CHROMOSOMES:
    print(diagnosed_chromosome, calc_zscore(diagnosed_chromosome, counts_per_chromosome))

