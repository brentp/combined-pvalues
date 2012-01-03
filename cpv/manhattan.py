"""
    %prog [options] files

plot a manhattan plot of the input file(s).
"""

import argparse
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from itertools import groupby, cycle
from operator import itemgetter
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from cpv._common import bediter, get_col_num

def chr_cmp(a, b):
    a, b = a[0], b[0]
    a = a.lower().replace("_", ""); b = b.lower().replace("_", "")
    achr = a[3:] if a.startswith("chr") else a
    bchr = b[3:] if b.startswith("chr") else b

    try:
        return cmp(int(achr), int(bchr))
    except ValueError:
        if achr.isdigit() and not bchr.isdigit(): return -1
        if bchr.isdigit() and not achr.isdigit(): return 1
        # X Y
        return cmp(achr, bchr)

def manhattan(fname, col_num, image_path, no_log, colors, title, lines, ymax,
             bonferonni=False, regions=None):
    """
    regions is keyed by chromosome with [(start, stop), ...] extents of
    the regions to highlight
    """
    xs, ys, cs = [], [], []
    region_xys = [] # highlight certain regions.
    colors = cycle(colors)
    chrom_centers = []

    last_x = 0
    nrows = 0
    giter = [(seqid, list(rlist)) for seqid, rlist \
        in groupby(bediter(fname, col_num), key=itemgetter('chrom'))]

    region_xs, region_ys = [], []
    for seqid, rlist in sorted(giter, cmp=chr_cmp):
        color = colors.next()
        nrows += len(rlist)
        # since chroms are on the same plot. add this chrom to the end of the
        # last chrom

        region_xs = [last_x + r['start'] for r in rlist]
        xs.extend(region_xs)
        ys.extend([r['p'] for r in rlist])
        cs.extend([color] * len(rlist))

        if regions and seqid in regions:
            regions_bounds = regions[seqid]
            region_xys.extend([(last_x + r['start'], r['p']) for r in rlist \
                  if any((s <= r['start'] <= e) for s, e in regions_bounds)])

        # save the middle of the region to place the label
        chrom_centers.append((seqid, (region_xs[0] + region_xs[-1]) / 2))
        # keep track so that chrs don't overlap.
        last_x = xs[-1]

    xs = np.array(xs)
    ys = np.array(ys) if no_log else -np.log10(ys)

    plt.close()
    f = plt.figure()
    ax = f.add_axes((0.1, 0.09, 0.88, 0.85))

    bonferonni_p = 0.05 / nrows

    if title is not None:
        plt.title(title)

    ax.set_ylabel('' if no_log else '-log10(p)')
    if lines:
        ax.vlines(xs, 0, ys, colors=cs, alpha=0.5)
    else:
        ax.scatter(xs, ys, s=1.2, c=cs, edgecolors='none')

    if regions:
        rxs, rys = zip(*region_xys)
        if not no_log: rys = -np.log10(rys)
        ax.scatter(rxs, rys, s=1.5, c='#ffff00', alpha=1.0, edgecolors='none')

    # plot 0.05 line after multiple testing. always nlog10'ed since
    # that's the space we're plotting in.
    if bonferonni:
        ax.axhline(y=-np.log10(bonferonni_p), color='0.5', linewidth=2)
    plt.axis('tight')
    plt.xlim(0, xs[-1])
    plt.ylim(ymin=0)
    if ymax is not None: plt.ylim(ymax=ymax)
    plt.xticks([c[1] for c in chrom_centers], [c[0] for c in chrom_centers], rotation=-90, size=8.5)
    #plt.show()
    print >>sys.stderr, "Bonferonni-corrected p-value for %i rows: %.3g" \
            % (nrows, 0.05 / nrows)
    print >>sys.stderr, "values less than Bonferonni-corrected p-value: %i " \
            % (ys > -np.log10(bonferonni_p)).sum()

    if False:
        ax_qq = f.add_axes((0.74, 0.12, 0.22, 0.22), alpha=0.2)

        pys = np.sort(10**-ys) # convert back to actual p-values
        qqplot(ys, ax_qq)

        ax_hist = f.add_axes((0.12, 0.12, 0.22, 0.22), frameon=True, alpha=0.6)
        hist(pys, ax_hist)

    print >>sys.stderr, "saving to: %s" % image_path
    plt.savefig(image_path)

    return image_path

def hist(pys, ax_hist):
    ax_hist.hist(pys, bins=40, color='0.75')
    ax_hist.set_xticks([])
    ax_hist.set_yticks([])

def qqplot(lpys, ax_qq):
    lunif = -np.log10(np.arange(1, len(lpys) + 1) / float(len(lpys)))[::-1]
    ax_qq.plot(lunif, np.sort(lpys), marker=',', linestyle='none')
    ax_qq.set_xticks([])
    ax_qq.plot(lunif, lunif, 'r--')
    ax_qq.set_xlabel('')
    ax_qq.set_ylabel('')
    ax_qq.set_yticks([])
    ax_qq.axis('tight')
    ax_qq.axes.set_frame_on(True)

def read_regions(fregions):
    if not fregions: return None
    regions = {}
    for toks in (l.split("\t") for l in open(fregions) if l[0] != "#"):
        if not toks[0] in regions: regions[toks[0]] = []
        regions[toks[0]].append((int(toks[1]), int(toks[2])))
    return regions

def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--no-log", help="the p-value is already -log10'd, don't "
                "re -log10", action='store_true', default=False)
    p.add_argument("-b", dest="bonferonni", 
            help="plot a line for the bonferonni of 0.05")
    p.add_argument("--col", dest="col", help="index of the column containing"
                   " the the p-value", default=-1, type=int)
    p.add_argument("--colors", dest="colors", help="cycle through these colors",
                default="bk")
    p.add_argument("--image", dest="image", help="save the image to this file."
                   " e.g. %(default)s", default="manhattan.png")
    p.add_argument("--title", help="title for the image.", default=None,
                   dest="title")
    p.add_argument("--ymax", help="max (-log) y-value for plot", dest="ymax",
                   type=float)
    p.add_argument("--lines", default=False, dest="lines", action="store_true",
        help="plot the p-values as lines extending from the x-axis rather than"
             " points in space. plotting will take longer with this option.")
    p.add_argument("--regions", dest="regions",
                   help="points in these bed regions are colored differently")

    p.add_argument('bed_file', help="bed-file to plot")

    args = p.parse_args()
    if (not args.bed_file):
        sys.exit(not p.print_help())

    column = get_col_num(args.col)
    regions = read_regions(args.regions)

    manhattan(args.bed_file, column, args.image, args.no_log, args.colors,
              args.title, args.lines, args.ymax, regions=regions,
             bonferonni=args.bonferonni)

if __name__ == "__main__":
    main()
