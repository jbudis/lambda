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
METHOD_NCV = 'NCV'
METHOD_SZ = 'SZ'
METHOD_FL = 'FL'
METHODS = [METHOD_NCV, METHOD_SZ, METHOD_FL]


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


def calc_lambda_score(fls, diagnosed_chromosome, p, step):
    chr_l = fls[diagnosed_chromosome, :step+1].sum()
    nl = fls[:, :step+1].sum()
    return (chr_l - nl*p) / np.sqrt(nl*p*(1-p))


def calc_max_lambda_score(fls, diagnosed_chromosome):
    p = fls[diagnosed_chromosome, :].sum() / fls.sum()
    print(max(calc_lambda_score(fls, diagnosed_chromosome, p, step) for step in range(MIN_STEP, MAX_STEP)))
    return max(calc_lambda_score(fls, diagnosed_chromosome, p, step) for step in range(MIN_STEP, MAX_STEP))
