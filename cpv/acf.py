"""
   calculate the autocorrelation of a *sorted* bed file with a set
   of *distance* lags
"""
import argparse
from toolshed import reader
import sys
import numpy as np
from itertools import groupby, tee, izip

def dist(beda, bedb):
    # not on same chrom
    if beda[0] != bedb[0]: return None
    # beda left-of bedb
    if beda[1] < bedb[0]:
        return bedb[0] - beda[1]
    # bedb left-of beda
    if bedb[1] < beda[0]:
        return beda[0] - bedb[1]
    # overlapping
    return 0

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def bediter(fname, col_num):
    for l in reader(fname, header=False):
        if l[0][0] == "#": continue
        yield  {"chrom": l[0], "start": int(l[1]), "end": int(l[2]),
                "p": float(l[col_num]), "stuff": l[3:][:]}

def acf(fnames, lags, col_num):
    acfs = {}
    for lag_min, lag_max in pairwise(lags):
        # groupby chromosome.
        xs, ys = [], []
        for fname in fnames:
            for key, chromlist in groupby(bediter(fname, col_num), lambda a: a["chrom"]):
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
    return acfs



def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", dest="d", help="start:stop:step of distance. e.g."
            " %default means check acf at distances of:"
            "[35, 85, 135, 185, 235, 285, 335, 385, 435, 485]",
            type=str, default="15:500:50")
    p.add_argument("-c", dest="c", help="column number that has the value to take the"
            " acf", type=int)
    p.add_argument('files', nargs='+', help='files to process')
    args = p.parse_args()
    if (args.d is None or len(args.files) == 0):
        sys.exit(not p.print_help())

    d = map(int, args.d.split(":"))
    assert len(d) == 3
    lags = range(*d)
    acf_vals = acf(args.files, lags, args.c - 1)
    for k,v in sorted(acf_vals.items()):
        print k, v

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
