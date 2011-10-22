"""
   calculate the autocorrelation of a *sorted* bed file with a set
   of *distance* lags.
"""
import argparse
from toolshed import reader
from _common import bediter
import numpy as np
from itertools import chain

def sim(all_ps, this_ps, N):
    """
    N times: take len(this_ps) random p-values
             from all_ps.
    return the number of times those simulated
    values are less than the values in this_ps
    by rank, and by actual value
    """
    all_rank = np.arange(len(all_ps))

def run(args):
    col_num = args.c if args.c < 0 else (args.c - 1)
    all_ps = np.array([b["p"] for b in bediter(args.pvals, col_num)])

    piter = chain(bediter(args.pvals, col_num), [None])
    prow = piter.next()
    # just use 2 for col_num, but dont need the p from regions.
    for toks in (l.rstrip("\r\n").split("\t") for l in open(args.regions)):
        rchrom = toks[0]
        rstart, rend = map(int, toks[1:3])
        ps = []
        while (prow["chrom"], prow["start"]) < (rchrom, rstart):
            prow = piter.next()
            if prow is None: break
        while (rchrom, rend) >= (prow["chrom"], prow["end"]):
            ps.append(prow)
            prow = piter.next()
            if prow is None: break
        if ps:
            sim(all_ps, ps, args.N)






def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", dest="pvals", help="BED containing all the p values")
    p.add_argument("-r", dest="regions", help="BED containing all the regions")

    p.add_argument("-N", dest="N", help="number of simulations to perform",
                   type=int, default=1000)
    p.add_argument("-c", dest="c", help="column number of PVALS containing the p-value"
                   " of interest", type=int, default=1)
    args = p.parse_args()
    if not (args.regions and args.pvals):
        sys.exit(not p.print_help())
    return run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
