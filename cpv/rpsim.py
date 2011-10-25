"""
   calculate the autocorrelation of a *sorted* bed file with a set
   of *distance* lags.
"""
import argparse
from toolshed import reader
from _common import bediter, get_col_num
import numpy as np
from itertools import chain
from slk import gen_sigma_matrix

from scipy.stats import norm
from numpy.linalg import cholesky as chol
qnorm = norm.ppf
pnorm = norm.cdf

def gen_correlated(sigma, n):
    """
    generate autocorrelated data according to the matrix
    sigma
    """
    sigma = np.asmatrix(sigma)
    X = np.random.uniform(0, 1, size=(sigma.shape[0], n))
    return pnorm(chol(sigma) * qnorm(X))

def calc_w(ps, truncate_at):
    # product of ps that are less than truncate_at
    # avoid underflow by taking log.
    return np.exp(np.sum(np.log(ps[ps <= truncate_at])))

def sim(sigma, ps, nsims, truncate):
    # see: https://gist.github.com/1306786#file_zaykin_truncated.py
    Y = gen_correlated(sigma, nsims)
    w0 = calc_w(ps, truncate)
    return sum(calc_w(row, truncate) <= w0 for row in Y.T) / float(nsims)

def run(args):
    col_num = get_col_num(args.c)
    #all_ps = np.array([b["p"] for b in bediter(args.pvals, col_num)])
    rpsim(args.pvals, args.regions, col_num, args.N, args.tau)

def rpsim(fpvals, fregions, col_num, nsims, tau):
    piter = chain(bediter(fpvals, col_num), [None])
    prow = piter.next()
    # just use 2 for col_num, but dont need the p from regions.
    no = -1
    for nr, toks in enumerate((l.rstrip("\r\n").split("\t")
                                   for l in open(fregions))):
        rchrom = toks[0]
        rstart, rend = map(int, toks[1:3])
        ps = []
        # grab the p-values in the bed file that are within the current region
        while (prow["chrom"], prow["start"]) < (rchrom, rstart):
            prow = piter.next()
            if prow is None: break
        while (rchrom, rend) >= (prow["chrom"], prow["end"]):
            ps.append(prow)
            prow = piter.next()
            if prow is None: break
        if ps:
            # TODO: calc ACF out to len(p). maybe collect ps for all
            # regions first and then create ACF for the longest one.
            sigma = gen_sigma_matrix(ps, acfs)
            psim = sim(sigma, ps, nsims)
            print "%s\t%.4g" % ("\t".join(toks), psim)
            no += 1
    assert nr == no, (nr, no)

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", dest="pvals", help="BED containing all the p values"
                  " used to generate `regions`")
    p.add_argument("-r", dest="regions", help="BED containing all the regions")
    p.add_argument("-t", dest="tau", help="tau cutoff", type=float,
                  default=0.05)
    p.add_argument("-N", dest="N", help="number of simulations to perform",
                   type=int, default=1000)
    p.add_argument("-c", dest="c", help="column number containing the p-value"
                   " of interest", type=int, default=1)
    args = p.parse_args()
    if not (args.regions and args.pvals):
        import sys
        sys.exit(not p.print_help())
    return run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
