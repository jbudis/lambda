#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
from glob import glob
import common
import yaml

np.random.seed(0)

# Parsing arguments

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


# Loading input files

count_files = glob('{}/*.tsv'.format(args.train_dir))
assert len(count_files) > 0, 'No *.tsv file in the provided train directory'

afls = np.array([common.load_counts(cf) for cf in count_files])


# NCV Score

def train_ncv(diagnosed_chromosome, fls):
    mapped_reads = fls.sum(axis=(1, 2))
    refs = fls[:, common.REFERENCE_CHROMOSOMES[diagnosed_chromosome], :].sum(axis=(1, 2)) / mapped_reads
    tris = fls[:, diagnosed_chromosome, :].sum(axis=1) / mapped_reads
    ratios = tris / refs
    return float(np.mean(ratios)), float(np.std(ratios))


# FL score

def train_fl(diagnosed_chromosome, fls):
    max_lambdas = np.array([common.calc_max_lambda_score(counts, diagnosed_chromosome) for counts in fls])
    return float(np.mean(max_lambdas)), float(np.std(max_lambdas))


# Error ellipse

def train_ellipse_variance(diagnosed_chromosome, fls, methods_params):

    size_counts = fls[:, :, :common.SZ_FL - common.MIN_FL]
    size_scores = np.array([common.calc_score_ncv(diagnosed_chromosome, counts,
                                                  methods_params[common.METHOD_SZ]) for counts in size_counts])

    fl_scores = np.array([common.calc_score_fl(diagnosed_chromosome, counts,
                                               methods_params[common.METHOD_FL]) for counts in fls])

    scores = np.array([size_scores, fl_scores])
    eigenvalues, eigenvectors = np.linalg.eig(np.corrcoef(scores))

    vx = float(max(eigenvalues[0], eigenvalues[1]))
    vy = float(min(eigenvalues[0], eigenvalues[1]))

    return vx, vy


# Calculating and storing parameters

params = {
    common.METHOD_NCV: {chromosome: {} for chromosome in common.DIAGNOSED_CHROMOSOMES},
    common.METHOD_SZ: {chromosome: {} for chromosome in common.DIAGNOSED_CHROMOSOMES},
    common.METHOD_FL: {chromosome: {} for chromosome in common.DIAGNOSED_CHROMOSOMES},
    common.METHOD_SZ_FL: {chromosome: {} for chromosome in common.DIAGNOSED_CHROMOSOMES}
}


for chromosome in common.DIAGNOSED_CHROMOSOMES:

    # NCV parameters
    ncv_mean, ncv_std = train_ncv(chromosome, afls)
    params[common.METHOD_NCV][chromosome][common.ATTR_MEAN] = ncv_mean
    params[common.METHOD_NCV][chromosome][common.ATTR_STD] = ncv_std

    # Size selection parameters
    size_fls = afls[:, :, :common.SZ_FL - common.MIN_FL]
    size_mean, size_std = train_ncv(chromosome, size_fls)
    params[common.METHOD_SZ][chromosome][common.ATTR_MEAN] = size_mean
    params[common.METHOD_SZ][chromosome][common.ATTR_STD] = size_std

    # FL parameters
    fl_mean, fl_std = train_fl(chromosome, afls)
    params[common.METHOD_FL][chromosome][common.ATTR_MEAN] = fl_mean
    params[common.METHOD_FL][chromosome][common.ATTR_STD] = fl_std

    # Size + FL parameters
    vx, vy = train_ellipse_variance(chromosome, afls, params)
    params[common.METHOD_SZ_FL][chromosome][common.ATTR_VX] = vx
    params[common.METHOD_SZ_FL][chromosome][common.ATTR_VY] = vy

with open(args.param_file, 'w') as out:
    yaml.dump(params, out, default_flow_style=False)
