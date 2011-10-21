"""
   calculate the autocorrelation of a *sorted* bed file with a set
   of *distance* lags.
"""
import argparse
from toolshed import reader
from chart import chart
import sys
from tempfile import mktemp
import os
import numpy as np
from itertools import groupby, tee, izip, combinations

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def bediter(fname, col_num):
    """
    iterate over a bed file. turn col_num into a float
    and the start, stop column into an int and yield a dict
    for each row.
    """
    for l in reader(fname, header=False):
        if l[0][0] == "#": continue
        p = float(l[col_num])
        if p == 1: p-= 1e-10 # the stouffer correction doesnt like values == 1
        yield  {"chrom": l[0], "start": int(l[1]), "end": int(l[2]),
                "p": p} # "stuff": l[3:][:]}

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

        acfs[(lag_min, lag_max)] = (np.corrcoef(xs, ys)[0, 1], len(xs))
    return sorted(acfs.items())

def get_corr(dist, acfs):
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
        # a is always left of b
        dist = jbed["start"] - ibed["end"]
        # symmetric.
        # cached speeds things up a bit...
        if not dist in cached:
            cached[dist] = get_corr(dist, acfs)
        a[j, i] = a[i, j] = cached[dist]
    return a

def adjust_pvals(fnames, col_num0, acfs):
    from stouffer_liptak import stouffer_liptak
    lag_max = acfs[-1][0][1]
    for fname in fnames:
        for key, chromlist in groupby(bediter(fname, col_num0), lambda a: a["chrom"]):
            for xbed, xneighbors in walk(chromlist, lag_max):

                sigma = gen_sigma_matrix(xneighbors, acfs)
                pvals = [g['p'] for g in xneighbors]
                adjusted = stouffer_liptak(pvals, sigma)["p"]
                yield (xbed["chrom"], xbed["start"], xbed["end"], xbed["p"],
                        adjusted)

def run(args):
    d = map(int, args.d.split(":"))
    d[1] += 1 # adjust for non-inclusive end-points...
    assert len(d) == 3
    lags = range(*d)
    acf_vals = acf(args.files, lags, args.c - 1)
    print >>sys.stderr, "#", chart([(k, v[0]) for k, v in acf_vals])
    print >>sys.stderr, "lag_min-lag_max\tcorrelation\tN"
    for k,v in sorted(acf_vals):
        print >>sys.stderr, "%i-%i\t%.4g\t%i" % (k[0], k[1], v[0], v[1])

    # get rid of N, just keep the correlation.
    acf_vals = [(k, v[0]) for k, v in acf_vals]
    # use a temp file because we don't want to keep all these in memory
    # during the FDR calc.
    tmp_fh = open(mktemp(), "w")
    adjusted = []
    for row in adjust_pvals(args.files, args.c - 1, acf_vals):
        adjusted.append(row[-1])
        tmp_fh.write("%s\t%i\t%i\t%.3g\t%.3g\n" % row)
    tmp_fh.close()
    from scikits.statsmodels.sandbox.stats.multicomp import fdrcorrection0
    rejected, bh_pvals = fdrcorrection0(adjusted, alpha=args.alpha,
            method='indep')[:2]

    write_file(tmp_fh.name, rejected, bh_pvals)

def write_file(fname, rejected, bh_pvals):
    tmp_fh = open(fname, "r")
    print "#chrom\tstart\tend\tp_orig\tp_stouffer\trejected\tp_bh"
    for line, rej, bhp in izip(tmp_fh, rejected, bh_pvals):
        sys.stdout.write("%s\t%s\t%.4g\n" % (line.rstrip("\r\n"),
                                       "T" if rej else "F",
                                       bhp))
    tmp_fh.close()
    os.unlink(tmp_fh.name)

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", dest="d", help="start:stop:stepsize of distance. e.g."
            " %default means check acf at distances of:"
            "[15, 65, 115, 165, 215, 265, 315, 365, 415, 465]",
            type=str, default="15:500:50")
    p.add_argument("-c", dest="c", help="column number that has the value to take the"
            " acf", type=int, default=4)
    p.add_argument("-a", dest="alpha", default=0.05, type=float, help="cutoff"
            " for significance after benjamini hochberg FDR correction")
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
