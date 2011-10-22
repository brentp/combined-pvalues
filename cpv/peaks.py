"""
find peaks or troughs in sorted bed files

for a bedgraph file with pvalues in the 4th column. usage would be:

    $ python peaks.py --dist 100 --seed 0.01 some.bed > some.regions.bed

where regions.bed contains the start and end of the region and (currently) the
lowest p-value in that region.
"""
from itertools import groupby
import operator
from toolshed import reader
import argparse
import sys

def bediter(fname, col_num):
    """
    iterate over a bed file. turn col_num into a float
    and the start, stop column into an int and yield a dict
    for each row.
    """
    for l in reader(fname, header=False):
        if l[0][0] == "#": continue
        yield  {"chrom": l[0], "start": int(l[1]), "end": int(l[2]),
                "p": float(l[col_num])} # "stuff": l[3:][:]}

# use class to keep track of written peaks.
class _write_peaks(object):
    def __init__(self):
        self.peak_count = 0
    def __call__(self, peaks, seed, out, scmp):
        # could have a list with only those passing the threshold.
        if not any(scmp(p["p"], seed) for p in peaks): return None
        if len(peaks) == 0: return None
        # peak_count unused...
        self.peak_count += 1
        peak_start = peaks[0]["start"]
        peak_end = peaks[-1]["end"]
        peak_count = len(peaks)
        # TODO: something better than keep best p-value ?
        pbest = peaks[0]["p"]
        for p in peaks:
            if scmp(p["p"], pbest): pbest = p["p"]
        out.write("%s\t%i\t%i\t%.4g\t%i\n" % (
            peaks[0]["chrom"], peak_start, peak_end, pbest, peak_count))

write_peaks = _write_peaks()

def walk(chromiter, thresh, seed, dist, out=None, scmp=operator.le):
    assert(scmp(seed, thresh))
    for key, bedlist in groupby(chromiter, lambda c: c["chrom"]):
        last_start = -1
        peaks = []
        for b in bedlist:
            assert last_start <= b["start"], ("enforce sorted")
            last_start = b["start"]
            # this comparison gets both thresh and seed.
            if scmp(b["p"], thresh):
                # we have to add this to peaks.
                # first check distance.
                # if distance is too great, we create a new peak
                if peaks != [] and b["start"] - peaks[-1]["end"] > dist:
                    if out is None:
                       if any(scmp(p, seed) for p in peaks):
                           for p in peaks: yield p
                    else:
                        write_peaks(peaks, seed, out, scmp)
                    peaks = []

                #add new peak regardless
                peaks.append(b)

        if out is None:
           if any(scmp(p, seed) for p in peaks):
               for p in peaks: yield p
        else:
            write_peaks(peaks, seed, out, scmp)


if __name__ == "__main__":
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--dist", dest="dist", help="Maximum dist to skip before "
             " finding a seed/thresh value. If this number is exceeded, the"
             " region is ended.", type=int)

    p.add_argument("--seed", dest="seed", help="A value must be at least this"
                 " large/small in order to seed a region.", type=float)
    p.add_argument("--threshold", dest="threshold", help="After seeding, a value"
                 " of at least this number can extend a region. ",
                 type=float, default=3.0)
    p.add_argument("--invert", dest="invert", action="store_true",
            help="by default, the test is for a value less-than seed or"
            " thresh--e.g. for p-values. If this flag is specified, the test"
            " is for greater-than--e.g. for scores or -log10(p-values)")
    p.add_argument("-c", type=int, help="column number containing the value "
                  "for which to find peaks.", default=4)
    p.add_argument("bed_file")
    try:
        args = p.parse_args()
    except TypeError:
        sys.exit(not p.print_help())

    if not(args.dist and args.seed):
        sys.exit(not p.print_help())

    if args.threshold is None:
        args.threshold = args.seed
        print >>sys.stderr, "setting threshold == seed"


    chromiter = bediter(args.bed_file, args.c - 1)
    scmp = operator.ge if args.invert else operator.le
    assert scmp(args.seed, args.threshold)
    # call list because the walk function is an iterator.
    list(walk(chromiter, args.threshold, args.seed, args.dist, sys.stdout))

