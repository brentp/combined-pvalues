#!/usr/bin/env python
import sys
__actions = ("acf", "slk", "fdr", "peaks", "region_p", "hist", "splot",
             "manhattan")

def main():
    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        print >>sys.stderr,"""
Tools for viewing and adjusting p-values in BED files.

   Contact: Brent Pedersen - bpederse@gmail.com
   License: BSD

To run, indicate one of:

   acf       - calculate autocorrelation within BED file
   slk       - Stouffer-Liptak-Kechris correction of correlated p-values
   fdr       - Benjamini-Hochberg correction of p-values
   peaks     - find peaks in a BED file.
   region_p  - generate SLK p-values for a region (of p-values)
   hist      - plot a histogram of a column and check for uniformity.
   splot     - a scatter plot of column(s) in a bed file for a given region.
   manhattan - a manhattan plot of values in a BED file.

NOTE: most of these assume *sorted* BED files.
SEE: https://github.com/brentp/combined-pvalues for documentation
    """
        sys.exit()
    if not sys.argv[1] in __actions:
        sys.argv.pop(1)
        _pipeline()
    else:
        module = __import__(sys.argv[1])
        # remove the action
        sys.argv.pop(1)
        module.main()

def _pipeline():
    import sys
    import os.path as op
    sys.path.insert(0, op.join(op.dirname(__file__), ".."))
    from cpv import acf, slk, fdr, peaks, region_p
    from _common import get_col_num
    import argparse
    import operator

    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-c", dest="c", help="column number that has the value to"
                   "take the  acf", type=int, default=4)
    p.add_argument("--tau", help="tau for the truncated product",
                   type=float, default=0.1)
    p.add_argument("--dist", dest="dist", help="Maximum dist to extend the"
             " ACF calculation", type=int)
    p.add_argument("--step", dest="step", help="step size for the"
             " ACF calculation", type=int)

    p.add_argument("--seed", dest="seed", help="A value must be at least this"
                 " large/small in order to seed a region.", type=float,
                 default=0.1)
    p.add_argument("--threshold", dest="threshold", help="After seeding, a value"
                 " of at least this number can extend a region. ",
                 type=float)
    p.add_argument("--random", action="store_true", default=False,
            help="If specified, sample from the uniform distribution rather "
            "than the observed p-values. The latter deflates significance.")

    p.add_argument("-p", dest="prefix", help="prefix for output files",
                   default=None)

    p.add_argument('files', nargs='+', help='files to process')

    args = p.parse_args()

    if not (args.prefix and args.step):
        sys.exit(p.print_help())

    if not args.threshold:
        args.threshold = args.seed

    lags = range(1, args.dist, args.step)
    lags.append(lags[-1] + args.step)
    col_num = get_col_num(args.c)

    acf_vals = acf.acf(args.files, lags, col_num, simple=False)

    with open(args.prefix + ".acf.txt", "w") as fh:
        acf_vals = acf.write_acf(acf_vals, fh)

    print >>sys.stderr, "ACF:\n", open(args.prefix + ".acf.txt").read()
    with open(args.prefix + ".slk.bed", "w") as fh:
        for row in slk.adjust_pvals(args.files, col_num, acf_vals):
            fh.write("%s\t%i\t%i\t%.4g\t%.4g\n" % row)
        print >>sys.stderr, "wrote: %s" % fh.name

    with open(args.prefix + ".fdr.bed", "w") as fh:
        for bh, l in fdr.fdr(args.prefix + ".slk.bed", -1, 0.05):
            fh.write("%s\t%.4g\n" % (l.rstrip("\r\n"), bh))
        print >>sys.stderr, "wrote: %s" % fh.name

    fregions = args.prefix + ".regions.bed"
    with open(fregions, "w") as fh:
        peaks.peaks(args.prefix + ".fdr.bed", -1, args.threshold, args.seed,
            args.step, fh, operator.le)
        print >>sys.stderr, "wrote: %s" % fh.name
    print >>sys.stderr, "%i regions" % (sum(1 for _ in open(fregions)))

    with open(args.prefix + ".sim-p.regions.bed", "w") as fh:
        fh.write("#chrom\tstart\tend\tmin-p\tn-probes\tslk-p\tslk-sidak-p\n")
        # use -2 for original, uncorrected p-values in slk.bed
        for region_line, slk_p, slk_sidak_p in region_p.region_p(
                               args.prefix + ".slk.bed",
                               args.prefix + ".regions.bed", -2,
                               5000, args.tau, args.step,
                               random=args.random):
            fh.write("%s\t%.4g\t%.4g\n" % (region_line, slk_p, slk_sidak_p))
        print >>sys.stderr, "wrote: %s" % fh.name

if __name__ == "__main__":
    main()
