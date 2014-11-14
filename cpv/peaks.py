"""
find peaks or troughs in sorted bed files

for a bedgraph file with pvalues in the 4th column. usage would be:

    $ python peaks.py --dist 100 --seed 0.01 some.bed > some.regions.bed

where some.regions.bed contains the start and end of the region and (currently)
the lowest p-value in that region.
"""
from itertools import groupby
import operator
from toolshed import reader
import argparse
import sys

# don't use _common because I want this to be stand-alone.
def bediter(fname, col_num):
    for i, l in enumerate(reader(fname, header=False)):
        if l[0][0] == "#": continue
        try:
            yield  {"chrom": l[0], "start": int(l[1]), "end": int(l[2]),
                "p": float(l[col_num])} # "stuff": l[3:][:]}
        except:
            print >>sys.stderr, l
            if i != 0:
                raise

# use class to keep track of written peaks.
def write_peaks(peaks, seed, out, scmp):
    # could have a list with only those passing the threshold.
    if len(peaks) == 0: return None
    # peak_count unused...
    peak_start = peaks[0]["start"]
    # dont konw the length of the regions and they are only sorted
    # by start.
    peak_end = max(p["end"] for p in peaks)
    peak_count = len(peaks)
    # TODO: something better than keep best p-value ? truncated product?
    pbest = peaks[0]["p"]
    for p in peaks:
        if scmp(p["p"], pbest): pbest = p["p"]
    out.write("%s\t%i\t%i\t%.4g\t%i\n" % (
        peaks[0]["chrom"], peak_start, peak_end, pbest, peak_count))

def trim_peaks(peaks, seed, thresh, scmp):
    """
    if thresh was greater than seed, we trim the region
    so the ends are < seed, but middle values can be seed < p < thresh
    """
    if seed == thresh: return peaks
    try:
        i_start = next(i for i, p in enumerate(peaks) if scmp(p['p'], seed))
    except StopIteration:
        return []
    i_end = len(peaks) - next(i for i, p in enumerate(reversed(peaks)) if scmp(p['p'], seed))
    return peaks[i_start:i_end]

def walk(chromiter, thresh, seed, dist, out=None, scmp=operator.le):
    assert(scmp(seed, thresh))
    for key, bedlist in groupby(chromiter, lambda c: c["chrom"]):
        last_start = -1
        # have to track max end because intervals are sorted only by start.
        max_end, peaks = 0, []
        for b in bedlist:
            assert last_start <= b["start"], ("enforce sorted", last_start, b)
            last_start = b["start"]
            # this comparison gets both thresh and seed.
            if scmp(b["p"], thresh):
                # we have to add this to peaks.
                # first check distance.
                # if distance is too great, we create a new peak
                if peaks != [] and b["start"] - max_end > dist:

                    peaks = trim_peaks(peaks, seed, thresh, scmp)
                    if out is None:
                        for p in peaks: yield p
                    else:
                        write_peaks(peaks, seed, out, scmp)
                    peaks = []
                    max_end = 0

                #add new peak regardless
                peaks.append(b)
                max_end = max(b['end'], max_end)

        if out is None:
           if any(scmp(p['p'], seed) for p in peaks):
               for p in peaks: yield p
        else:
            write_peaks(peaks, seed, out, scmp)

def peaks(fbedfile, col_num, threshold, seed, dist, fout, scmp):
    chromiter = bediter(fbedfile, col_num)
    # TODO: make this yield...
    for _ in walk(chromiter, threshold, seed, dist, fout, scmp):
        yield _

def run(args):
    col_num = args.c if args.c < 0 else args.c - 1
    scmp = operator.ge if args.invert else operator.le
    assert scmp(args.seed, args.threshold)
    # call list because the walk function is an iterator.
    for peak in peaks(args.bed_file, col_num, args.threshold, args.seed, args.dist,
            sys.stdout, scmp):
        yield peak

def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--dist", dest="dist", help="Maximum dist to skip before "
             " finding a seed/thresh value. If this number is exceeded, the"
             " region is ended.", type=int)

    p.add_argument("--seed", dest="seed", help="A value must be at least this"
                 " large/small in order to seed a region.", type=float)
    p.add_argument("--threshold", dest="threshold", help="After seeding, a value"
                 " of at least this number can extend a region. ",
                 type=float)
    p.add_argument("--invert", dest="invert", action="store_true",
            help="by default, the test is for a value less-than seed or"
            " thresh--e.g. for p-values. If this flag is specified, the test"
            " is for greater-than--e.g. for scores or -log10(p-values)")
    p.add_argument("-c", type=int, help="column number containing the value "
                  "for which to find peaks.", default=-1)
    p.add_argument("bed_file")

    args = p.parse_args()

    if not(args.dist and args.seed):
        sys.exit(not p.print_help())

    if args.threshold is None:
        args.threshold = args.seed
    for peak in run(args):
        print peak

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()

