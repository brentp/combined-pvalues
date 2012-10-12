import sys
import os.path as op

def main():
    import argparse
    from _common import get_col_num

    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-c", dest="c", help="column number that has the value to"
                   "take the  acf", default='4')
    p.add_argument("--dist", dest="dist", help="Maximum dist to extend the"
             " ACF calculation", type=int)
    p.add_argument("--step", dest="step", help="step size for bins in the"
             " ACF calculation", type=int)
    p.add_argument("--seed", dest="seed", help="A value must be at least this"
                 " large/small in order to seed a region.", type=float,
                 default=0.1)
    p.add_argument("--threshold", dest="threshold", help="After seeding, a value"
                 " of at least this number can extend a region. ",
                 type=float)
    p.add_argument("-p", "--prefix", dest="prefix",
            help="prefix for output files", default=None)

    p.add_argument("--mlog", dest="mlog", action="store_true",
                   default=False, help="do the correlation on the -log10 of"
                   "the p-values. Default is to do it on the raw values")

    p.add_argument("--region-filter-p", help="max adjusted region-level p-value to be reported"
                 "in final output", type=float, default=1)

    p.add_argument("--region-filter-n", help="require at least this many probes"
                 "for a region to be reported in final output", type=int, default=1)

    p.add_argument('bed_files', nargs='+', help='sorted bed file to process')

    args = p.parse_args()


    if not (args.prefix):
        sys.exit(p.print_help())

    if not args.threshold:
        args.threshold = args.seed
    assert op.exists(args.bed_files[0])

    col_num = get_col_num(args.c, args.bed_files[0])
    return pipeline(col_num, args.step, args.dist, args.prefix,
            args.threshold, args.seed,
            args.bed_files, mlog=args.mlog,
		region_filter_p=args.region_filter_p,
		region_filter_n=args.region_filter_n)

def pipeline(col_num, step, dist, prefix, threshold, seed, bed_files, mlog=False,
    region_filter_p=1, region_filter_n=1):
    sys.path.insert(0, op.join(op.dirname(__file__), ".."))
    from cpv import acf, slk, fdr, peaks, region_p, stepsize, filter
    import operator

    if step is None:
        step = stepsize.stepsize(bed_files, col_num)
        print >>sys.stderr, "calculated stepsize as: %i" % step

    lags = range(1, dist, step)
    lags.append(lags[-1] + step)

    # go out to max requested distance but stop once an autocorrelation
    # < 0.05 is added.

    putative_acf_vals = acf.acf(bed_files, lags, col_num, simple=False,
                                mlog=mlog)
    acf_vals = []
    for a in putative_acf_vals:
        # a is ((lmin, lmax), (corr, N))
        # this heuristic seems to work. stop just above the 0.08 correlation
        # lag.
        if a[1][0] < 0.04 and len(acf_vals) > 2: break
        acf_vals.append(a)
        if a[1][0] < 0.04 and len(acf_vals): break

    # save the arguments that this was called with.
    with open(prefix + ".args.txt", "w") as fh:
        print >>fh, " ".join(sys.argv[1:]) + "\n"
        import datetime
        print >>fh, "date: %s" % datetime.datetime.today()

    with open(prefix + ".acf.txt", "w") as fh:
        acf_vals = acf.write_acf(acf_vals, fh)
        print >>sys.stderr, "wrote: %s" % fh.name
    print >>sys.stderr, "ACF:\n", open(prefix + ".acf.txt").read()

    with open(prefix + ".slk.bed", "w") as fh:
        for row in slk.adjust_pvals(bed_files, col_num, acf_vals):
            fh.write("%s\t%i\t%i\t%.4g\t%.4g\n" % row)
        print >>sys.stderr, "wrote: %s" % fh.name

    with open(prefix + ".fdr.bed", "w") as fh:
        for bh, l in fdr.fdr(prefix + ".slk.bed", -1):
            fh.write("%s\t%.4g\n" % (l.rstrip("\r\n"), bh))
        print >>sys.stderr, "wrote: %s" % fh.name

    fregions = prefix + ".regions.bed"
    with open(fregions, "w") as fh:
        peaks.peaks(prefix + ".fdr.bed", -1, threshold, seed,
            step, fh, operator.le)
    n_regions = sum(1 for _ in open(fregions))
    print >>sys.stderr, "wrote: %s (%i regions)" % (fregions, n_regions)

    with open(prefix + ".regions-p.bed", "w") as fh:
        N = 0
        fh.write("#chrom\tstart\tend\tmin_p\tn_probes\tslk_p\tslk_sidak_p\n")
        # use -2 for original, uncorrected p-values in slk.bed
        for region_line, slk_p, slk_sidak_p, sim_p in region_p.region_p(
                               prefix + ".slk.bed",
                               prefix + ".regions.bed", -2,
                               0, step):
            fh.write("%s\t%.4g\t%.4g\n" % (region_line, slk_p, slk_sidak_p))

            fh.flush()
            N += int(slk_sidak_p < 0.05)
        print >>sys.stderr, "wrote: %s, (regions with corrected-p < 0.05: %i)" \
                % (fh.name, N)

    regions_bed = fh.name
    with open(prefix + ".regions-t.bed", "w") as fh:
        N = 0
        for i, toks in enumerate(filter.filter(bed_files[0], regions_bed)):
	    if i == 0: toks[0] = "#" + toks[0]
            else:
		if float(toks[6]) > region_filter_p: continue
		if int(toks[4]) < region_filter_n: continue
                N += 1
            print >>fh, "\t".join(toks)
        print >>sys.stderr, "wrote: %s, (regions with region-p < %.3f and n-probes >= %i: %i)" \
                % (fh.name, region_filter_p, region_filter_n, N)
