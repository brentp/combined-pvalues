"""
"""
import argparse
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from _common import read_acf, bediter, get_col_num
from itertools import groupby, combinations
from stouffer_liptak import stouffer_liptak


def get_corr(dist, acfs):
    """
    extract the correlation from the acf sigma matrix
    given a distance.
    """
    # it's very close. just give it the next up.
    # TODO: should probably not do this. force them to start at 0.
    # acfs[0] is like (lag_min, lag_max), corr
    # so this is checking if it's < smallest lag...
    if dist < acfs[0][0][0]:
        return acfs[0][1]
    for (lag_min, lag_max), corr in acfs:
        if lag_min <= dist <= lag_max:
            return corr
    return 0

def walk(chromlist, lag_max):
    """
    for each item in chromlist, yield the item and its neighborhood
    within lag-max. These yielded values are then used to generate
    the sigma autocorrelation matrix.
    """
    L = list(chromlist)
    N = len(L)
    imin = imax = 0
    for ithis, xbed in enumerate(L):
        # move up the bottom of the interval
        while xbed["start"] - L[imin]["end"] > lag_max:
            imin += 1
        if imax == N: imax -= 1
        while L[imax]["start"] - xbed["end"] < lag_max:
            imax += 1
            if imax == N: break
        assert imin <= ithis <= imax
        # dont need to add 1 to imax because we got outside of the range above.
        yield xbed, L[imin: imax]

def gen_sigma_matrix(group, acfs, cached={}):
    a = np.eye(len(group))
    group = enumerate(group)
    for (i, ibed), (j, jbed) in combinations(group, 2):
        # j is always right of i. but could overlap
        dist = jbed["start"] - ibed["end"]
        if dist < 0: dist = 0
        # symmetric.
        # cached speeds things up a bit...
        if not dist in cached:
            cached[dist] = get_corr(dist, acfs)
        a[j, i] = a[i, j] = cached[dist]
    return a

def slk_chrom(chromlist, lag_max, acfs):
    """
    calculate the slk for a given chromosome
    """
    for xbed, xneighbors in walk(chromlist, lag_max):

        sigma = gen_sigma_matrix(xneighbors, acfs)
        pvals = [g['p'] for g in xneighbors]
        r = stouffer_liptak(pvals, sigma)

        yield (xbed["chrom"], xbed["start"], xbed["end"], xbed["p"],
                r["p"])
        if not r["OK"]:
            print >>sys.stderr, "# non-invertible %s\t%i\t%i" % \
                    (xbed["chrom"], xbed["start"], xbed["end"])

def _slk_chrom(args):
    return list(slk_chrom(*args))

def adjust_pvals(fnames, col_num0, acfs):
    lag_max = acfs[-1][0][1]

    # parallelize if multiprocesing is installed.
    try:
        from multiprocessing import Pool
        pool = Pool()
        imap = pool.imap
    except ImportError:
        import itertools
        imap = itertools.imap

    for fname in fnames:
        arg_iter = ((list(chromlist), lag_max, acfs) for key, chromlist in
                        groupby(bediter(fname, col_num0), lambda a: a["chrom"]))

        for results in imap(_slk_chrom, arg_iter):
            for r in results:
                yield r


def run(args):
    acf_vals = read_acf(args.acf)
    col_num = get_col_num(args.c)
    for row in adjust_pvals(args.files, col_num, acf_vals):
        sys.stdout.write("%s\t%i\t%i\t%.3g\t%.3g\n" % row)

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--acf", dest="acf", help="acf file containing the lagged "
                   "correlations")
    p.add_argument("-c", dest="c", help="column number that has the value to take the"
            " acf", type=int, default=4)
    p.add_argument('files', nargs='+', help='files to process')
    args = p.parse_args()
    if (len(args.files) == 0):
        sys.exit(not p.print_help())
    return run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
