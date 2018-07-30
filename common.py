import pandas as pd
import numpy as np

AUTOSOMES = ['chr{}'.format(i) for i in range(22)]
MIN_FL, MAX_FL, SZ_FL = 50, 220, 150

FLS = list(range(MIN_FL, MAX_FL))
CHR13, CHR18, CHR21 = 12, 17, 20
DIAGNOSED_CHROMOSOMES = [CHR13, CHR18, CHR21]
REFERENCE_CHROMOSOMES = {
    CHR13: [2, 3, 6],
    CHR18: [3, 6, 7, 8, 9, 15],
    CHR21: [0, 3, 7, 9, 18, 19, 21]
}
ATTR_MEAN = 'mean'
ATTR_STD = 'std'
ATTR_VX = 'vx'
ATTR_VY = 'vy'

METHOD_NCV = 'NCV'
METHOD_SZ = 'SZ'
METHOD_FL = 'FL'
METHOD_NCV_FL = 'NCV+FL'
METHOD_SZ_FL = 'SZ+FL'
METHODS = [METHOD_NCV, METHOD_SZ, METHOD_FL, METHOD_NCV_FL, METHOD_SZ_FL]


def load_counts(count_file):

    counts = pd.read_csv(count_file, sep='\t', header=0, index_col=0).values
    n_rows, n_cols = counts.shape

    assert n_rows == len(AUTOSOMES), 'Number or rows of input count file ' \
                                     'does not correspond to {} autosomes'.format(len(AUTOSOMES))

    assert n_cols == len(FLS), 'Number of columns does not correspond to ' \
                               'number of used fragment lengths, from {} to {}, ' \
                               'i.e. {} columns. ' \
                               'Loading file {}'.format(MIN_FL, MAX_FL, MAX_FL - MIN_FL, count_file)
    return counts


MIN_STEP, MAX_STEP = 75, 96


# NCV Score

def calc_score_ncv(diagnosed_chromosome, counts, method_params):

    counts_per_chromosome = counts.sum(axis=1)

    mapped_reads = counts_per_chromosome.sum()
    refs = counts_per_chromosome[REFERENCE_CHROMOSOMES[diagnosed_chromosome]].sum() / mapped_reads
    tris = counts_per_chromosome[diagnosed_chromosome] / mapped_reads
    ratio = tris / refs

    m = method_params[diagnosed_chromosome][ATTR_MEAN]
    s = method_params[diagnosed_chromosome][ATTR_STD]

    return (ratio - m) / s


# FL score

def calc_lambda_score(fls, diagnosed_chromosome, p, step):
    chr_l = fls[diagnosed_chromosome, :step+1].sum()
    nl = fls[:, :step+1].sum()
    return (chr_l - nl*p) / np.sqrt(nl*p*(1-p))


def calc_max_lambda_score(fls, diagnosed_chromosome):
    p = fls[diagnosed_chromosome, :].sum() / fls.sum()
    return max(calc_lambda_score(fls, diagnosed_chromosome, p, step) for step in range(MIN_STEP, MAX_STEP))


def calc_score_fl(diagnosed_chromosome, counts, method_params):

    max_lambda = calc_max_lambda_score(counts, diagnosed_chromosome)

    m = method_params[diagnosed_chromosome][ATTR_MEAN]
    s = method_params[diagnosed_chromosome][ATTR_STD]
    return (max_lambda - m) / s
