"""
   calculate the autocorrelation of a *sorted* bed file with a set
   of *distance* lags.
"""
import argparse
from toolshed import reader
import sys
import numpy as np
from itertools import groupby, tee, izip

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def bediter(fname, col_num):
    for l in reader(fname, header=False):
        if l[0][0] == "#": continue
        yield  {"chrom": l[0], "start": int(l[1]), "end": int(l[2]),
                "p": float(l[col_num])} # "stuff": l[3:][:]}


def acf(fnames, lags, col_num0):
    acfs = {}
    for lag_min, lag_max in pairwise(lags):
        # groupby chromosome.
        xs, ys = [], []
        for fname in fnames:
            for key, chromlist in groupby(bediter(fname, col_num0), lambda a: a["chrom"]):
                chromlist = list(chromlist)
                for ix, xbed in enumerate(chromlist):
                    for iy in xrange(ix + 1, len(chromlist)):
                        ybed = chromlist[iy]
                        # y is always > x
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
    if dist < acfs[0][0][0]:
        return acfs[0][1]
    for (lag_min, lag_max), corr in acfs:
        if lag_min <= dist <= lag_max:
            return corr
    return 0

def adjust_pvals(fnames, col_num0, acfs):

    lag_max = acfs[-1][0][1]
    D = len(acfs)
    for fname in fnames:
        for key, chromlist in groupby(bediter(fname, col_num0), lambda a: a["chrom"]):
            chromlist = list(chromlist)
            for ix, xbed in enumerate(chromlist):
                sigma = [1.0]
                pvals = [xbed['p']]
                dists = [0]
                # have to go backward,
                for iy in xrange(max(ix - 5 * D, 0), len(chromlist)):
                    if ix == iy: continue
                    ybed = chromlist[iy]
                    dist = abs(ybed['start'] - xbed['end']) if ix < iy else abs(xbed['start'] - ybed['end'])
                    if dist > lag_max: break
                    dists.append(dist)
                    pvals.append(ybed['p'])
                    sigma.append(get_corr(dist, acfs))
                # convert sigma to a matrix
                sigma, pvals, dists = map(np.array, (sigma, pvals, dists))
                idxs = dists.argsort()
                print "i", idxs

                print "s", sigma[idxs]
                print "p", pvals[idxs]
                print "d", dists[idxs]
                1/0




def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", dest="d", help="start:stop:step of distance. e.g."
            " %default means check acf at distances of:"
            "[15, 65, 115, 165, 215, 265, 315, 365, 415, 465]",
            type=str, default="15:500:50")
    p.add_argument("-c", dest="c", help="column number that has the value to take the"
            " acf", type=int)
    p.add_argument("--adjust", dest="adjust", default=False,
        action="store_true", help="after the acf, adjust the p-values")
    p.add_argument('files', nargs='+', help='files to process')
    args = p.parse_args()
    if (args.d is None or len(args.files) == 0):
        sys.exit(not p.print_help())

    d = map(int, args.d.split(":"))
    assert len(d) == 3
    lags = range(*d)
    acf_vals = acf(args.files, lags, args.c - 1)
    # get rid of N, just keep the correlation.
    acf_vals = [(k, v[0]) for k, v in acf_vals]

    adjust_pvals(args.files, args.c - 1, acf_vals)
    for k,v in sorted(acf_vals):
        print k, v

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
