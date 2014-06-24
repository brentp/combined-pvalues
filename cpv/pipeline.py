import sys
import os.path as op
import gzip

def main():
    import argparse
    from _common import get_col_num

    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-c", dest="c", help="column number that has the value to"
                   "take the  acf", default='4')
    p.add_argument("--dist", "--distance", dest="dist", help="Maximum dist to extend the"
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
    p.add_argument("-z", "--z-score", action="store_true", default=False,
            help="use z-score correction instead of liptak")

    p.add_argument("--genomic-control", dest="genomic_control",
            help="perform the genomic control correction on the input"
            " pvalues", action="store_true", default=False)

    p.add_argument("--mlog", "--nlog", dest="mlog", action="store_true",
                   default=False, help="do the correlation on the -log10 of"
                   "the p-values. Default is to do it on the raw values")

    p.add_argument("--region-filter-p", help="max adjusted region-level p-value"
                 " to be reported "
                 "in final output. this requires the input bed file to have"
                 " chrom, start, end, 't' columns", type=float, default=1)

    p.add_argument("--region-filter-n", help="require at least this many probes"
                 "for a region to be reported in final output. "
                 " this requires the input bed file to have chrom, start, "
                 "end, 't' columns", type=int, default=None)
    p.add_argument("--annotate", help="annotate with refGen from this db" \
            "in UCSC (e.g. hg19) requires cruzdb", default=None)

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
            region_filter_n=args.region_filter_n,
            genome_control=args.genomic_control,
            db=args.annotate,
            z=args.z_score)

def pipeline(col_num, step, dist, prefix, threshold, seed, bed_files, mlog=False,
    region_filter_p=1, region_filter_n=None, genome_control=False, db=None,
    z=False):
    sys.path.insert(0, op.join(op.dirname(__file__), ".."))
    from cpv import acf, slk, fdr, peaks, region_p, stepsize, filter
    from cpv._common import genome_control_adjust, genomic_control, bediter
    import operator


    if step is None:
        step = stepsize.stepsize(bed_files, col_num)
        print >>sys.stderr, "calculated stepsize as: %i" % step

    lags = range(1, dist, step)
    lags.append(lags[-1] + step)

    prefix = prefix.rstrip(".")
    putative_acf_vals = acf.acf(bed_files, lags, col_num, simple=False,
                                mlog=mlog)
    acf_vals = []
    # go out to max requested distance but stop once an autocorrelation
    # < 0.05 is added.
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
        from .__init__ import __version__
        print >>fh, "version:", __version__

    with open(prefix + ".acf.txt", "w") as fh:
        acf_vals = acf.write_acf(acf_vals, fh)
        print >>sys.stderr, "wrote: %s" % fh.name
    print >>sys.stderr, "ACF:\n", open(prefix + ".acf.txt").read()

    spvals, opvals = [], []
    with open(prefix + ".slk.bed", "w") as fhslk:

        for row in slk.adjust_pvals(bed_files, col_num, acf_vals, z=z):
            fhslk.write("%s\t%i\t%i\t%.4g\t%.4g\n" % row)
            opvals.append(row[-2])
            spvals.append(row[-1])

    print >>sys.stderr, "# original lambda: %.2f" % genomic_control(opvals)
    del opvals

    gc_lambda = genomic_control(spvals)
    print >>sys.stderr, "wrote: %s with lambda: %.2f" % (fhslk.name, gc_lambda)

    if genome_control:
        fhslk = open(prefix + ".slk.gc.bed", "w")
        adj = genome_control_adjust([d['p'] for d in bediter(prefix + ".slk.bed", -1)])
        for i, line in enumerate(open(prefix + ".slk.bed")):
            print >>fhslk, "%s\t%.5g" % (line.rstrip("\r\n"), adj[i])

        fhslk.close()
        print >>sys.stderr, "wrote: %s" % fhslk.name

    with open(prefix + ".fdr.bed", "w") as fh:
        for bh, l in fdr.fdr(fhslk.name, -1):
            fh.write("%s\t%.4g\n" % (l.rstrip("\r\n"), bh))
        print >>sys.stderr, "wrote: %s" % fh.name

    fregions = prefix + ".regions.bed"
    with open(fregions, "w") as fh:
        list(peaks.peaks(prefix + ".fdr.bed", -1, threshold, seed,
            step, fh, operator.le))
    n_regions = sum(1 for _ in open(fregions))
    print >>sys.stderr, "wrote: %s (%i regions)" % (fregions, n_regions)

    with open(prefix + ".regions-p.bed", "w") as fh:
        N = 0
        fh.write("#chrom\tstart\tend\tmin_p\tn_probes\t{correction}_p\t{correction}_sidak_p\n".format(correction=('z'
            if z else 'slk')))
        # use -2 for original, uncorrected p-values in slk.bed
        for region_line, slk_p, slk_sidak_p, sim_p in region_p.region_p(
                               prefix + ".slk.bed",
                               prefix + ".regions.bed", -2,
                               0, step, mlog=mlog, z=z):
            fh.write("%s\t%.4g\t%.4g\n" % (region_line, slk_p, slk_sidak_p))
            fh.flush()
            N += int(slk_sidak_p < 0.05)
        print >>sys.stderr, "wrote: %s, (regions with corrected-p < 0.05: %i)" \
                % (fh.name, N)

    regions_bed = fh.name
    header = (gzip.open(bed_files[0]) if bed_files[0].endswith(".gz")
            else open(bed_files[0])).next().split("\t")
    #if all(h in header for h in ('t', 'start', 'end')):
    if region_filter_p != 1 or region_filter_n != None:
        if region_filter_n is None: region_filter_n = 0
        with open(prefix + ".regions-t.bed", "w") as fh:
            N = 0
            for i, toks in enumerate(filter.filter(bed_files[0], regions_bed,
                p_col_name=col_num)):
                if i == 0: toks[0] = "#" + toks[0]
                else:
                    if float(toks[6]) > region_filter_p: continue
                    if int(toks[4]) < region_filter_n: continue
                    N += 1
                print >>fh, "\t".join(toks)
            print >>sys.stderr, ("wrote: %s, (regions with region-p"
                                "< %.3f and n-probes >= %i: %i)") \
                    % (fh.name, region_filter_p, region_filter_n, N)

    try:
        from cpv import manhattan
        regions = manhattan.read_regions(fh.name)

        manhattan.manhattan(prefix + ".slk.bed", 3, prefix.rstrip(".") + ".manhattan.png",
                         False, ['#959899', '#484B4C'], "", False, None,
                         regions=regions, bonferonni=True)
    except ImportError:
        pass # they dont have matplotlib


    if db is not None:
        from cruzdb import Genome
        g = Genome(db)
        lastf = fh.name
        with open(prefix + ".anno.%s.bed" % db, "w") as fh:
            fh.write('#')
            g.annotate(lastf, ("refGene", "cpgIslandExt"), out=fh,
                    feature_strand=True, parallel=len(spvals) > 500)
        print >>sys.stderr, "wrote: %s annotated with %s" % (fh.name, db)
