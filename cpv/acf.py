"""
   calculate the autocorrelation of a *sorted* bed file with a set
   of *distance* lags.
"""
import argparse
from array import array
import sys
import numpy as np
import scipy.stats as ss
from itertools import groupby, izip, chain
from _common import bediter, pairwise, get_col_num, get_map


def create_acf_list(lags):
    acfs = []
    for lag_min, lag_max in pairwise(lags):
        acfs.append((lag_min, lag_max,
            # array uses less memory than list.
            {"x": array("f"), "y": array("f")}))
    acfs.reverse()
    return acfs

def _acf_by_chrom(args):
    """
    calculate the ACF for a single chromosome
    chromlist is the data for a single chromsome
    """
    chromlist, lags = args
    acfs = create_acf_list(lags)
    if not isinstance(chromlist, list):
        chromlist = list(chromlist)
    max_lag = max(a[1] for a in acfs)
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
    return acfs


def merge_acfs(unmerged):
    """
    utitlity function to merge the chromosomes after
    they've been calculated, and before the correlation
    is calculated.
    """
    merged = unmerged.pop()
    for um in unmerged:
        # have to merge at each lag.
        for (glag_min, glag_max, gxys), (ulag_min, ulag_max, uxys) in \
                                                izip(merged, um):
            assert glag_min == ulag_min and glag_max == ulag_max
            gxys["x"].extend(uxys["x"])
            gxys["y"].extend(uxys["y"])
            # reduce copies in memory.
            uxys = {}
    return merged


def acf(fnames, lags, col_num0, partial=True, simple=False, mlog=True):
    """
    calculate the correlation of the numbers in `col_num0` from the bed files
    in `fnames` at various lags. The lags are specified by distance. Partial
    autocorrelation may be calculated as well.

    Since the bed files may be very large, this attempts to be as memory
    efficient as possible while still being very fast for a pure python
    implementation.
    """
    # reversing allows optimization below.
    imap = get_map()

    arg_list = [] # chaining
    for fname in fnames:
        # groupby chromosome.
        arg_list = chain(arg_list, ((list(chromlist), lags) for chrom, \
                    chromlist in \
                    groupby(bediter(fname, col_num0), lambda a: a["chrom"])))

    unmerged_acfs = [] # separated by chrom. need to merge later.
    for chrom_acf in imap(_acf_by_chrom, arg_list):
        unmerged_acfs.append(chrom_acf)

    acfs = merge_acfs(unmerged_acfs)
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
        if mlog:
            xs[xs == 0] = 1e-12
            ys[ys == 0] = 1e-12
            xs, ys = -np.log10(xs), -np.log10(ys)
        #slope, intercept, corr, p_val, stderr = ss.linregress(xs, ys)
        # NOTE: using pearson correlation, which assumes normality.
        # could switch to spearman as below.
        corr, p_val = ss.spearmanr(xs, ys)
        if simple:
            acf_res[(lmin, lmax)] = corr
        else:
            acf_res[(lmin, lmax)] = (corr, len(xs), p_val)
    return sorted(acf_res.items())

def run(args):
    """
    general function that takes an args object (from argparse)
    with the necessary options and calls acf()
    """
    d = map(int, args.d.split(":"))
    assert len(d) == 3, ("-d argument must in in the format start:end:step")
    d[1] += 1 # adjust for non-inclusive end-points...
    lags = range(*d)

    acf_vals = acf(args.files, lags, get_col_num(args.c), partial=(not
                                                            args.full))
    write_acf(acf_vals, sys.stdout)

def write_acf(acf_vals, out):
    # write acf to a file and return only [((lag_min, lag_max), corr)...]
    simple_acf = []
    values = [float(v[0]) for k, v in acf_vals]
    xlabels = "|".join("%s-%s" % k for k, v in acf_vals)
    print >> out, "#lag_min\tlag_max\tcorrelation\tN\tp"
    for k, v in sorted(acf_vals):
        print >> out, "%i\t%i\t%.4g\t%i\t%.4g" % (k[0], k[1], v[0], v[1], v[2])
        simple_acf.append((k, v[0]))
    return simple_acf

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", dest="d", help="start:stop:stepsize of distance. e.g."
            " %(default)s means check acf at distances of:"
            "[15, 65, 115, 165, 215, 265, 315, 365, 415, 465]",
            type=str, default="15:500:50")
    p.add_argument("-c", dest="c", help="column number with p-values for acf "
                   "calculations", type=int, default=4)
    p.add_argument("--full", dest="full", action="store_true",
                   default=False, help="do full autocorrelation (default"
                   " is partial)")
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
