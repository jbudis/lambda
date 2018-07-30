import pandas as pd

AUTOSOMES = ['chr{}'.format(i) for i in range(22)]
MIN_FL, MAX_FL = 50, 220
FLS = list(range(MIN_FL, MAX_FL))
CHR13, CHR18, CHR21 = 12, 17, 20
DIAGNOSED_CHROMOSOMES = [CHR13, CHR18, CHR21]


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