"""
   calculate the autocorrelation of a *sorted* bed file with a set
   of *distance* lags.
"""
import argparse
import os
from _common import bediter
from itertools import izip
import sys
from scikits.statsmodels.sandbox.stats.multicomp import fdrcorrection0

def run(args):
    # get rid of N, just keep the correlation.
    col_num = args.c if args.c < 0 else (args.c - 1)
    pvals = [b["p"] for b in bediter(args.bed_file, col_num)]
    bh_pvals = fdrcorrection0(pvals, alpha=args.alpha,
            method='indep')[1]
    for bh, l in izip(bh_pvals, open(args.bed_file)):
        print "%s\t%.4g" % (l.rstrip("\r\n"), bh)

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-c", dest="c", help="column number of the pvalues",
                   type=int, default=-1)
    p.add_argument("--alpha", dest="alpha", default=0.05, type=float, help="cutoff"
            " for significance after benjamini hochberg FDR correction")
    p.add_argument('bed_file', help='bed file to correct')
    args = p.parse_args()
    return run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
