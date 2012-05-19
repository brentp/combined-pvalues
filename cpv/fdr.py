"""
perform Benjamini-Hochberg FDR correction on a BED file with p-values.
"""
import argparse
from _common import bediter
from itertools import izip
from _common import get_col_num
import numpy as np

def run(args):
    # get rid of N, just keep the correlation.
    col_num = get_col_num(args.c)
    if args.null is not None:
       col_null = get_col_num(args.null)
       # TODO: call qvality.

    for bh, l in fdr(args.bed_file, col_num):
        print "%s\t%.4g" % (l.rstrip("\r\n"), bh)

def fdr(fbed_file, col_num):
    from scikits.statsmodels.sandbox.stats.multicomp import fdrcorrection0
    pvals = np.array([b["p"] for b in bediter(fbed_file, col_num)],
                        dtype=np.float64)
    bh_pvals = fdrcorrection0(pvals, method='indep')[1]
    fh = open(fbed_file)
    line = fh.readline()
    # drop header
    if not (line[0] == "#" or line.split()[0] == "chrom"):
        fh.seek(0)
    for bh, l in izip(bh_pvals, fh):
        yield bh, l

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-c", dest="c", help="column number of the pvalues",
                   type=int, default=-1)
    p.add_argument("--null", type=int,
        help="""(optional) column number of the pvalues under the null, e.g.
             for shuffled data. if given, qvality (which must be on the path)
             is used to do the correction. Otherwise, Benjamini-Hochberg is
             used""")
    p.add_argument('bed_file', help='bed file to correct')
    args = p.parse_args()
    return run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
