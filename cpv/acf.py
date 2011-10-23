"""
   calculate the autocorrelation of a *sorted* bed file with a set
   of *distance* lags.
"""
import argparse
from chart import chart
import sys
import os
import numpy as np
from itertools import groupby, combinations
from _common import bediter, pairwise

try:
    import scikits.statsmodels.api as sm
    from scikits.statsmodels.stats.stattools import durbin_watson

    def dw(xs, ys):
        model = sm.OLS(xs, ys)
        r = model.fit()
        return durbin_watson(r.resid)
    HAS_DW = False # leave for now, not working.
except ImportError:
    def dw(xs, ys): return None
    HAS_DW = False

def acf(fnames, lags, col_num0):
    acfs = {}
    for lag_min, lag_max in pairwise(lags):
        xs, ys = [], [] # these hold the lagged values.
        for fname in fnames:
            # groupby chromosome.
            for key, chromlist in groupby(bediter(fname, col_num0), lambda a: a["chrom"]):
                chromlist = list(chromlist)
                for ix, xbed in enumerate(chromlist):
                    # find all lines within lag of xbed.
                    for iy in xrange(ix + 1, len(chromlist)):
                        ybed = chromlist[iy]
                        # y is always > x so dist calc is simplified.
                        # too close:
                        if ybed['start'] - xbed['end'] < lag_min: continue
                        # too far.
                        if ybed['start'] - xbed['end'] > lag_max: break

                        xs.append(xbed['p'])
                        ys.append(ybed['p'])
        xs, ys = np.array(xs), np.array(ys)
        acfs[(lag_min, lag_max)] = (np.corrcoef(xs, ys)[0, 1], len(xs),
                                    dw(xs, ys))
    return sorted(acfs.items())

def run(args):
    d = map(int, args.d.split(":"))
    d[1] += 1 # adjust for non-inclusive end-points...
    assert len(d) == 3
    lags = range(*d)
    acf_vals = acf(args.files, lags, args.c - 1)
    values = [float(v[0]) for k, v in acf_vals]
    xlabels = "|".join("%s-%s" % k for k, v in acf_vals)
    print "#", chart(values, xlabels)
    print "lag_min\tlag_max\tcorrelation\tN" + \
                      ("\tdurbin-watson" if HAS_DW else "")
    for k,v in sorted(acf_vals):
        line = "%i\t%i\t%.4g\t%i" % (k[0], k[1], v[0], v[1])
        if HAS_DW: line += "\t%.2f" % v[2]
        print line

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", dest="d", help="start:stop:stepsize of distance. e.g."
            " %(default)s means check acf at distances of:"
            "[15, 65, 115, 165, 215, 265, 315, 365, 415, 465]",
            type=str, default="15:500:50")
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
