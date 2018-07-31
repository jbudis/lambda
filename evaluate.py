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


# Combination of scores


def r_chi_square_distance(distance):
    command = '''
        newZS = function(distance)  
            qnorm(pchisq(distance, df=2, lower.tail=FALSE, log.p=TRUE), lower.tail=FALSE, log.p=TRUE)
    '''
    rpy2.robjects.r(command)
    distance = rpy2.robjects.numpy2ri.py2ri(np.array([distance]))
    rz = rpy2.robjects.r.newZS(distance)
    return rpy2.robjects.numpy2ri.ri2py(rz)


def r_chi_square_scores(z1_np, z2_np):
    distance = z1_np**2 + z2_np**2
    return r_chi_square_distance(distance)[0]


# Ellipse correction

def calc_ellipse_dist(diagnosed_chromosome, zx, zy, method_params):
    vx = method_params[diagnosed_chromosome][common.ATTR_VX]
    vy = method_params[diagnosed_chromosome][common.ATTR_VY]
    alpha = np.radians(45)
    distance = (zx*np.cos(alpha) + zy*np.sin(alpha))**2 / vx + \
               (zx*np.sin(alpha) - zy*np.cos(alpha))**2 / vy
    return r_chi_square_distance(distance)[0]


# Storing results

METHOD, CHROMOSOME, SCORE = 'Method', 'Chromosome', 'Score'

zscore_items = []
for chromosome in common.DIAGNOSED_CHROMOSOMES:

    # NCV score
    score_ncv = common.calc_score_ncv(chromosome, all_counts, params[common.METHOD_NCV])
    zscore_items.append({
        METHOD: common.METHOD_NCV,
        CHROMOSOME: chromosome+1,
        SCORE: score_ncv
    })

    # Size score
    score_sz = common.calc_score_ncv(chromosome, size_counts, params[common.METHOD_SZ])
    zscore_items.append({
        METHOD: common.METHOD_SZ,
        CHROMOSOME: chromosome+1,
        SCORE: score_sz
    })

    # FL score
    score_fl = common.calc_score_fl(chromosome, all_counts, params[common.METHOD_FL])
    zscore_items.append({
        METHOD: common.METHOD_FL,
        CHROMOSOME: chromosome+1,
        SCORE: score_fl
    })

    # NCV + FL score
    zscore_items.append({
        METHOD: common.METHOD_NCV_FL,
        CHROMOSOME: chromosome+1,
        SCORE: r_chi_square_scores(score_ncv, score_fl)
    })

    # SIZE + FL score
    zscore_items.append({
        METHOD: common.METHOD_SZ_FL,
        CHROMOSOME: chromosome+1,
        SCORE: calc_ellipse_dist(chromosome, score_sz, score_fl, params[common.METHOD_SZ_FL])
    })


zscores = pd.DataFrame(zscore_items)
print(zscores.pivot(index=METHOD, columns=CHROMOSOME, values=SCORE))
