#!/usr/bin/env python3

from argparse import ArgumentParser
# import pandas as pd
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

parser.add_argument('counts', help='TSV file with the number or reads '
                                   'for each assumed fragment length (50-220, organised in columns) '
                                   'and autosomes (chr1..chr22, organised in rows)')

args = parser.parse_args()

print(args)
