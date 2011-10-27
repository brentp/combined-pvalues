"""
   simulate a p-value of a region using the truncated product method described
   in Zaykin et al. 2002, "Truncated Product Method for Combining p-values".
   Genet Epidemiol.
"""
import argparse
import sys
import numpy as np

from _common import bediter, get_col_num
from itertools import chain
from slk import gen_sigma_matrix
from acf import acf

from scipy.stats import norm
from numpy.linalg import cholesky as chol
qnorm = norm.ppf
pnorm = norm.cdf

def gen_correlated(sigma, n, X=None):
    """
    generate autocorrelated data according to the matrix
    sigma. if X is None, then data will be sampled from
    the uniform distibution. Otherwise, it will be sampled
    from X. Where X is then *all* observed
    p-values.
    """
    sigma = np.asmatrix(sigma)
    if X is None:
        X = np.random.uniform(0, 1, size=(sigma.shape[0], n))
    else:
        idxs = np.random.random_integers(0, len(X) - 1,
                                         size=sigma.shape[0] * n)
        X = X[idxs].reshape((sigma.shape[0], n))

    return pnorm(chol(sigma) * qnorm(X))

def calc_w(ps, truncate_at):
    # product of ps that are less than truncate_at
    # avoid underflow by taking log.
    return np.exp(np.sum(np.log(ps[ps <= truncate_at])))

def sim(sigma, ps, nsims, truncate, sample_distribution=None):
    # see: https://gist.github.com/1306786#file_zaykin_truncated.py
    assert isinstance(ps[0], (int, float, long))
    B = 0.
    w0 = calc_w(ps, truncate)
    for i in range(10):
        Y = gen_correlated(sigma, nsims/10, sample_distribution)
        B += sum(calc_w(row, truncate) <= w0 for row in Y.T)
    return B / nsims

def run(args):
    col_num = get_col_num(args.c)
    rpsim(args.pvals, args.regions, col_num, args.N, args.tau, args.step,
            args.random)

def _gen_acf(region_info, fpvals, col_num, step):
    # calculate the ACF as far out as needed...
    # [1] is the length of each region.
    max_len = max(r[1] for r in region_info)
    print >>sys.stderr, "# calculating ACF out to: %i" % max_len

    lags = range(1, max_len, step)
    if lags[-1] < max_len: lags.append(lags[-1] + step)
    print >>sys.stderr, "# with lags: %s" % lags
    if len(lags) > 100:
        print >>sys.stderr, "# !! this could take a looong time"
        print >>sys.stderr, "# !!!! consider using a larger step size (-s)"
    acfs = acf(fpvals, lags, col_num, simple=True)
    print >>sys.stderr, "# Done with one-time ACF calculation"
    return acfs

def rpsim(fpvals, fregions, col_num, nsims, tau, step, random=False):
    piter = chain(bediter(fpvals, col_num), [None])
    prow = piter.next()
    # just use 2 for col_num, but dont need the p from regions.
    region_info = []

    if(sum(1 for _ in open(fregions) if _[0] != "#") == 0):
        print >>sys.stderr, "no regions in %s" % (fregions, )
        sys.exit()

    for nr, region_line in enumerate((l.rstrip("\r\n")
                                   for l in open(fregions))):
        toks = region_line.split("\t")
        rchrom = toks[0]
        rstart, rend = map(int, toks[1:3])
        prows = []
        # grab the p-values in the bed file that are within the current region
        while (prow["chrom"] != rchrom or prow["start"] < rstart):
            prow = piter.next()
            if prow is None: break
        while (rchrom, rend) >= (prow["chrom"], prow["end"]):
            prows.append(prow)
            prow = piter.next()
            if prow is None: break
        assert prows, (region_line)
        region_len = rend - rstart + 1
        region_info.append((region_line, region_len, prows))
        del prows
    assert nr + 1 == len(region_info), (nr, len(region_info))

    acfs = _gen_acf(region_info, (fpvals,), col_num, step)
    # regions first and then create ACF for the longest one.
    if not random:
        sample_distribution = np.array([b["p"] for b in bediter(fpvals,
                                                                col_num)])
    else:
        sample_distribution = None
    for region_line, region_len, prows in region_info:
        # gen_sigma expects a list of bed dicts.
        sigma = gen_sigma_matrix(prows, acfs)
        ps = np.array([prow["p"] for prow in prows])
        psim = sim(sigma, ps, nsims, tau, sample_distribution)
        print "%s\t%.4g" % (region_line, psim)

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", dest="pvals", help="BED containing all the p values"
                  " used to generate `regions`")
    p.add_argument("-r", dest="regions", help="BED containing all the regions")
    p.add_argument("-t", dest="tau", help="tau cutoff", type=float,
                  default=0.05)
    p.add_argument("--random", default=False, action="store_true",
            help="for simulations, the default is to sample from the p-values"
            " from -p argument; if this flag is set, it will sample from"
            " the uniform distribution instead")
    p.add_argument("-s", dest="step", type=int, default=50,
            help="step size for acf calculation. should be the same "
            " value as the step sent to -d arg for acf")
    p.add_argument("-N", dest="N", help="number of simulations to perform",
                   type=int, default=2000)
    p.add_argument("-c", dest="c", help="column number containing the p-value"
                   " of interest", type=int, default=-1)
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
