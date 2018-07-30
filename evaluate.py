#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
import common
import yaml
import pandas as pd
import rpy2.robjects
import rpy2.robjects.numpy2ri
import rpy2.robjects.pandas2ri

np.random.seed(0)

# Parsing arguments

parser = ArgumentParser(description='Calculates standard and FL-based score for diagnosed sample'
                                    'and predicts diagnosis, i.e. positive, false positive, uninformative, negative')

parser.add_argument('counts', help='TSV file with the corrected number or reads '
                                   'for each assumed fragment length (50-220, organised in columns) '
                                   'and autosomes (chr1..chr22, organised in rows)')

parser.add_argument('param_file', help='YAML file with trained parameters (output of the train.py script)')

args = parser.parse_args()


# Loading input files

with open(args.param_file) as param_f:
    params = yaml.load(param_f)

all_counts = common.load_counts(args.counts)
size_counts = all_counts[:, :common.SZ_FL-common.MIN_FL]


# NCV Score

def calc_score_ncv(diagnosed_chromosome, counts, method_params):

    counts_per_chromosome = counts.sum(axis=1)

    mapped_reads = counts_per_chromosome.sum()
    refs = counts_per_chromosome[common.REFERENCE_CHROMOSOMES[diagnosed_chromosome]].sum() / mapped_reads
    tris = counts_per_chromosome[diagnosed_chromosome] / mapped_reads
    ratio = tris / refs

    m = method_params[diagnosed_chromosome][common.ATTR_MEAN]
    s = method_params[diagnosed_chromosome][common.ATTR_STD]

    return (ratio - m) / s


# FL score

def calc_score_fl(diagnosed_chromosome, counts, method_params):

    max_lambda = common.calc_max_lambda_score(counts, diagnosed_chromosome)

    m = method_params[diagnosed_chromosome][common.ATTR_MEAN]
    s = method_params[diagnosed_chromosome][common.ATTR_STD]
    return (max_lambda - m) / s


# Combination of scores

def r_chi_square(z1_np, z2_np):
    command = 'newZ = function(z1, z2) qnorm(pchisq(z1**2 +z2**2, df=2, lower.tail=FALSE, log.p=TRUE), lower.tail=FALSE, log.p=TRUE)'
    rpy2.robjects.r(command)
    z1 = rpy2.robjects.numpy2ri.py2ri(np.array([z1_np]))
    z2 = rpy2.robjects.numpy2ri.py2ri(np.array([z2_np]))
    rz = rpy2.robjects.r.newZ(z1, z2)
    return rpy2.robjects.numpy2ri.ri2py(rz)[0]


# Storing results

METHOD, CHROMOSOME, SCORE = 'Method', 'Chromosome', 'Score'

zscore_items = []
for chromosome in common.DIAGNOSED_CHROMOSOMES:

    # NCV score
    score_ncv = calc_score_ncv(chromosome, all_counts, params[common.METHOD_NCV])
    zscore_items.append({
        METHOD: common.METHOD_NCV,
        CHROMOSOME: chromosome+1,
        SCORE: score_ncv
    })

    # Size score
    score_sz = calc_score_ncv(chromosome, size_counts, params[common.METHOD_SZ])
    zscore_items.append({
        METHOD: common.METHOD_SZ,
        CHROMOSOME: chromosome+1,
        SCORE: score_sz
    })

    # FL score
    score_fl = calc_score_fl(chromosome, all_counts, params[common.METHOD_FL])
    zscore_items.append({
        METHOD: common.METHOD_FL,
        CHROMOSOME: chromosome+1,
        SCORE: score_fl
    })

    # NCV + FL score
    zscore_items.append({
        METHOD: common.METHOD_NCV_FL,
        CHROMOSOME: chromosome+1,
        SCORE: r_chi_square(score_ncv, score_fl)
    })


zscores = pd.DataFrame(zscore_items)
print(zscores.pivot(index=METHOD, columns=CHROMOSOME, values=SCORE))
