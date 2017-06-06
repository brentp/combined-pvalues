"""
   calculate the step-size that should be used for the ACF calculations.
   The step-size is calculated as::
       median(distance-between-adjacent-starts)
   This heuristic seems to work well for creating bins with equal amounts of
   records for the ACF.
"""
import argparse
from _common import get_col_num, bediter
from operator import itemgetter
from itertools import groupby
import numpy as np

def stepsize(bed_files, col):

    D1 = []
    for bed_file in bed_files:
        for _, chromlist in groupby(bediter(bed_file, col), itemgetter('chrom')):
            L = list(chromlist)
            if len(L) < 2: continue

            last_start = 0
            for i, ibed in enumerate(L):
                assert ibed['start'] >= last_start
                # look around ibed. nearest could be up or down-stream
                if i + 2 == len(L): break
                D1.append(L[i + 1]['start'] - ibed['start'])
        # round up to the nearest 10
    return int(round(np.median(D1) + 5, -1))

def run(args):
    col = get_col_num(args.c)
    print stepsize((args.bed_file,), col)

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-c", dest="c", help="column number that has the value to"
                   "take the  acf", type=int, default=-1)

    p.add_argument('bed_file', help='bed file to process')
    run(p.parse_args())

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
