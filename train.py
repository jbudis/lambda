#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
from glob import glob
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

parser = ArgumentParser(description='Trains parameters for FL a NCV models.'
                                    'Parameters are then utilized in the evaluate.py script to infer predictions '
                                    'for aneuploidy of diagnosed fetus')

parser.add_argument('train_dir', help='Directory with TSV count files for sample with healthy euploid fetus. '
                                      'All files with ".tsv" suffix in the directory would be used for training. '
                                      'Count file contains the corrected number or reads '
                                      'for each assumed fragment length (50-220, organised in columns) '
                                      'and autosomes (chr1..chr22, organised in rows)')

parser.add_argument('param_file', help='Output YAML file with trained parameters')

args = parser.parse_args()

count_files = glob('{}/*.tsv'.format(args.train_dir))
assert len(count_files) > 0, 'No *.tsv file in the provided train directory'

afls = np.array([common.load_counts(cf) for cf in count_files])


def train_sehnert(ref, fls):
    mapped_reads = fls.sum(axis=(1, 2))
    refs = fls[:, common.REFERENCE_CHROMOSOMES[ref], :].sum(axis=(1, 2)) / mapped_reads
    tris = fls[:, ref, :].sum(axis=1) / mapped_reads
    ratios = tris / refs
    return np.mean(ratios), np.std(ratios)


m13, s13 = train_sehnert(common.CHR13, afls)
m18, s18 = train_sehnert(common.CHR18, afls)
m21, s21 = train_sehnert(common.CHR21, afls)

with open(args.param_file, 'w') as out:
    yaml.dump({
        'NCV': {
            common.CHR13: {
                common.ATTR_MEAN: float(m13),
                common.ATTR_STD: float(s13)
            },
            common.CHR18: {
                common.ATTR_MEAN: float(m18),
                common.ATTR_STD: float(s18)
            },
            common.CHR21: {
                common.ATTR_MEAN: float(m21),
                common.ATTR_STD: float(s21)
            }
        }
    }, out, default_flow_style=False)
