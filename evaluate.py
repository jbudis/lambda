#!/usr/bin/env python3

from argparse import ArgumentParser
import pandas as pd
# import numpy as np
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
# np.random.seed(0)

parser = ArgumentParser(description='Calculates standard and FL-based score for diagnosed sample'
                                    'and predicts diagnosis, i.e. positive, false positive, uninformative, negative')

parser.add_argument('counts', help='TSV file with the corrected number or reads '
                                   'for each assumed fragment length (50-220, organised in columns) '
                                   'and autosomes (chr1..chr22, organised in rows)')

args = parser.parse_args()

AUTOSOMES = ['chr{}'.format(i) for i in range(22)]
MIN_FL, MAX_FL = 50, 220
FLS = list(range(MIN_FL, MAX_FL))

counts = pd.read_csv(args.counts, sep='\t', header=0, index_col=0).values
n_rows, n_cols = counts.shape

assert n_rows == len(AUTOSOMES), 'Number or rows of input count file ' \
                                 'does not correspond to {} autosomes'.format(len(AUTOSOMES))

assert n_cols == len(FLS), 'Number of columns does not correspond to ' \
                           'number of used fragment lengths, from {} to {}, ' \
                           'i.e. {} columns'.format(MIN_FL, MAX_FL, MAX_FL - MIN_FL)

