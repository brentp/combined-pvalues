"""
draw a histogram of the distribution of a given column
and check for uniformity with the chisq test.
"""
import argparse
from _common import get_col_nums
import toolshed as ts
from itertools import cycle

def run(args):
    col_nums = get_col_nums(args.c)

    # if they didn't specify labels, just the the column numbers
    labels = [x.strip() for x in args.labels.strip().split(",")] \
            if args.labels else map(str, [c + 1 for c in col_nums])

    ps = dict((c, []) for c in col_nums)

    file_iter =  (l.rstrip("\r\n").split("\t")
                  for l in ts.nopen(args.file) if l[0] != "#")

    chrom, start_end = args.region.split(":")
    start, end = map(int, start_end.split("-"))

    for b in file_iter:
        if b[0] != chrom: continue
        bstart, bend = int(b[1]), int(b[2])
        if bstart < start: continue
        if bend > end: break
        for c in col_nums:
            ps[c].append((bstart, bend, float(b[c])))

    plot(col_nums, ps, labels, args.lines, chrom)

def plot(col_nums, ps, labels, lines, chrom, colors=cycle("brgkyc")):
    from matplotlib import pyplot as plt
    f = plt.figure(figsize=(8, 3))
    ax = f.add_axes((0.08, 0.03, 0.9, 0.94))

    # choose the symbol based on the number of points.
    sym = '' if lines else \
                    ('o' if len(ps[col_nums[0]]) < 600 else \
                    ('.' if len(ps[col_nums[0]]) < 1200 else ","))
    lw = 2.5
    for col, label in zip(col_nums, labels):
        xs, xends, ys = zip(*ps[col])
        if len(xs) < 600:
            plt.hlines(ys, xs, xends, colors=colors.next(), label=label,
                    linewidth=lw, alpha=0.8)
            if lw >= 1: lw -= 0.6
        else:
            plt.plot(xs, ys, colors.next() + sym, label=label)
    rng = (max(xs) - min(xs))
    rng = ("%.2fKB" % (rng / 1000.)) if rng > 2000 else str(rng) + "bp"
    txts = (chrom, min(xs), max(xs), rng, len(xs))
    plt.text(0.02, 0.92, '%s:%i-%i (%s, %i probes)' % txts,
            transform=ax.transAxes)
    plt.legend()
    plt.xticks([])
    plt.ylabel("p-value")
    plt.ylim(ymin=-0.02, ymax=1.02)
    plt.show()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-c", dest="c", help="column number(s) to plot."
                   "specify multiple columns seperated by a ','.",
                   type=str, default='-1')
    p.add_argument("--labels", help="optional labels for the columns",
                   type=str, default=None)
    p.add_argument("-r", help="region to plot, e.g. chrY:1000-2000",
                   type=str, default=None, dest="region")
    p.add_argument("--lines", help="plot as lines (default is points)",
                   action="store_true", default=False)
    p.add_argument('file', help='bed file to correct')
    args = p.parse_args()
    return run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
