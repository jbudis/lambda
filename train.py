#!/usr/bin/env python3

from argparse import ArgumentParser
import pandas as pd
import numpy as np
from glob import glob
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

parser = ArgumentParser(description='Trains parameters for FL a NCV models.'
                                    'Parameters are then utilized in the evaluate.py script to infer predictions '
                                    'for aneuploidy of diagnosed fetus')

parser.add_argument('train_dir', help='Directory with TSV count files. All files with ".tsv" suffix in the directory '
                                      'would be used for training. '
                                      'Count file contains the corrected number or reads '
                                      'for each assumed fragment length (50-220, organised in columns) '
                                      'and autosomes (chr1..chr22, organised in rows)')

parser.add_argument('param_file', help='Output YAML file with trained parameters')

args = parser.parse_args()

count_files = glob('{}/*.tsv'.format(args.train_dir))
assert len(count_files) > 0, 'No *.tsv file in the provided train directory'

