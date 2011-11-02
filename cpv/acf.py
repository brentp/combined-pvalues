"""
   calculate the autocorrelation of a *sorted* bed file with a set
   of *distance* lags.
"""
import argparse
from array import array
from chart import chart
import sys
import numpy as np
from itertools import groupby
from _common import bediter, pairwise, get_col_num

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

def acf(fnames, lags, col_num0, partial=True, simple=False):
    """
    calculate the correlation of the numbers in `col_num0` from the bed files
    in `fnames` at various lags. The lags are specified by distance. Partial
    autocorrelation may be calculated as well.

    Since the bed files may be very large, this attempts to be as memory
    efficient as possible while still being very fast for a pure python
    implementation.
    """
    acfs = []
    for lag_min, lag_max in pairwise(lags):
        acfs.append((lag_min, lag_max,
            # array uses less memory than list.
            {"x": array("f"), "y": array("f") }))
    # reversing allows optimization below.
    acfs = acfs[::-1]

    max_lag = max(a[1] for a in acfs)
    for fname in fnames:
        # groupby chromosome.
        for key, chromlist in groupby(bediter(fname, col_num0), lambda a: a["chrom"]):
            chromlist = list(chromlist)
            for ix, xbed in enumerate(chromlist):
                # find all lines within lag of xbed.
                for iy in xrange(ix + 1, len(chromlist)):
                    ybed = chromlist[iy]
                    # y is always > x so dist calc is simplified.
                    dist = ybed['start'] - xbed['end']
                    if dist > max_lag: break

                    for lag_min, lag_max, xys in acfs:
                        # can break out of loop because we reverse-sorted acfs
                        # above. this is partial, but we merge below if needed.
                        if lag_min <= dist < lag_max:
                            xys["x"].append(xbed['p'])
                            xys["y"].append(ybed['p'])
                        elif dist > lag_max:
                            break
    acf_res = {}
    xs = np.array([], dtype='f')
    ys = np.array([], dtype='f')
    # iterate over it backwards and remove to reduce memory.
    while len(acfs):
        lmin, lmax, xys = acfs.pop()
        if partial:
            xs, ys = np.array(xys["x"]), np.array(xys["y"])
        else:
            # add the inner layers as we move out.
            xs = np.hstack((xs, xys["x"]))
            ys = np.hstack((ys, xys["y"]))
        if len(xs) == 0:
            print >>sys.stderr, "no values found at lag: %i-%i. skipping" \
                    % (lmin, lmax)
            continue
        if simple:
            acf_res[(lmin, lmax)] = np.corrcoef(xs, ys)[0, 1]
        else:
            acf_res[(lmin, lmax)] = (np.corrcoef(xs, ys)[0, 1], len(xs),
                                    dw(xs, ys))
    return sorted(acf_res.items())

def run(args):
    """
    general function that takes an args object (from argparse)
    with the necessary options and calls acf()
    """
    d = map(int, args.d.split(":"))
    d[1] += 1 # adjust for non-inclusive end-points...
    assert len(d) == 3
    lags = range(*d)
    acf_vals = acf(args.files, lags, get_col_num(args.c), partial=(not
                                                            args.full))
    write_acf(acf_vals, sys.stdout)

def write_acf(acf_vals, out):
    # write acf to a file and return only [((lag_min, lag_max), corr)...]
    simple_acf = []
    values = [float(v[0]) for k, v in acf_vals]
    xlabels = "|".join("%s-%s" % k for k, v in acf_vals)
    print >>out, "#", chart(values, xlabels)
    print >> out, "#lag_min\tlag_max\tcorrelation\tN" + \
                      ("\tdurbin-watson" if HAS_DW else "")
    for k,v in sorted(acf_vals):
        line = "%i\t%i\t%.4g\t%i" % (k[0], k[1], v[0], v[1])
        if HAS_DW: line += "\t%.2f" % v[2]
        print >> out, line
        simple_acf.append((k, v[0]))
    return simple_acf

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", dest="d", help="start:stop:stepsize of distance. e.g."
            " %(default)s means check acf at distances of:"
            "[15, 65, 115, 165, 215, 265, 315, 365, 415, 465]",
            type=str, default="15:500:50")
    p.add_argument("-c", dest="c", help="column number that has the value to"
                   "take the  acf", type=int, default=4)
    p.add_argument("--full", dest="full", action="store_true",
                   default=False, help="do full autocorrelation (default"
                   " is partial")
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
