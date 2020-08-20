from __future__ import print_function
import sys
import array
import os.path as op
import toolshed as ts

def main():
    import argparse
    from ._common import get_col_num

    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-c", dest="c", help="column number that has the value to"
                   "take the  acf", default='4')
    p.add_argument("--dist", "--distance" "--peak-dist", dest="dist", help="Maximum dist to "
            " search for adjacent peaks.", type=int, required=True)
    p.add_argument("--acf-dist", help="distance/window-size to use for "
            " smoothing. Defaults to 1/3 * peak-dist ", type=int, default=None)

    p.add_argument("--step", dest="step", help="step size for bins in the"
             " ACF calculation", type=int)
    p.add_argument("--seed", dest="seed", help="A value must be at least this"
                 " large/small in order to seed a region.", type=float,
                 default=0.05)
    p.add_argument("--threshold", dest="threshold", help="After seeding, a value"
                 " of at least this number can extend a region. ",
                 type=float)
    p.add_argument("--no-fdr", dest="no_fdr", help="Don't use FDR-corrected p-values "
            "for finding peaks (either way, we still do multiple-testing correction "
            "on the p-values for the regions).", action='store_true',
            default=False)
    p.add_argument("-p", "--prefix", dest="prefix",
            help="prefix for output files", default=None)

    p.add_argument("--genomic-control", dest="genomic_control",
            help="perform the genomic control correction on the input"
            " pvalues", action="store_true", default=False)

    p.add_argument("--region-filter-p", help="max adjusted region-level p-value"
                 " to be reported "
                 "in final output. this requires the input bed file to have"
                 " chrom, start, end, 't' columns", type=float, default=1)

    p.add_argument("--region-filter-n", help="require at least this many probes"
                 "for a region to be reported in final output. "
                 " this requires the input bed file to have chrom, start, "
                 "end, 't' columns", type=int, default=None)
    p.add_argument("--annotate", help="annotate with a gene table from this db " \
            "in UCSC (e.g. hg19) requires cruzdb", default=None)
    p.add_argument("--table", help="annotate with this gene table from a db " \
            "in UCSC (default is refGene) requires cruzdb", default='refGene')

    p.add_argument('bed_files', nargs='+', help='sorted bed file to process')

    args = p.parse_args()

    if not (args.prefix):
        sys.exit(p.print_help())

    if not args.threshold:
        args.threshold = args.seed
    assert op.exists(args.bed_files[0])

    if args.acf_dist is None:
        args.acf_dist = int(round(0.33333 * args.dist, -1))
        sys.stderr.write("setting --acf-dist to 0.33 * --dist == %i\n" %
                args.acf_dist)

    col_num = get_col_num(args.c, args.bed_files[0])
    return pipeline(col_num, args.step,
            args.dist, args.acf_dist, args.prefix,
            args.threshold, args.seed, args.table,
            args.bed_files,
            region_filter_p=args.region_filter_p,
            region_filter_n=args.region_filter_n,
            genome_control=args.genomic_control,
            db=args.annotate,
            use_fdr=not args.no_fdr)

def pipeline(col_num, step, dist, acf_dist, prefix, threshold, seed, table,
        bed_files, mlog=True, region_filter_p=1, region_filter_n=None,
        genome_control=False, db=None, use_fdr=True):
    sys.path.insert(0, op.join(op.dirname(__file__), ".."))
    from cpv import acf, slk, fdr, peaks, region_p, stepsize, filter
    from cpv._common import genome_control_adjust, genomic_control, bediter
    import operator


    if step is None:
        step = min(acf_dist, stepsize.stepsize(bed_files, col_num))
        print("calculated stepsize as: %i" % step, file=sys.stderr)

    lags = list(range(1, acf_dist, step))
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
        print(" ".join(sys.argv[1:]) + "\n", file=fh)
        import datetime
        print("date: %s" % datetime.datetime.today(), file=fh)
        from .__init__ import __version__
        print("version:", __version__, file=fh)

    with open(prefix + ".acf.txt", "w") as fh:
        acf_vals = acf.write_acf(acf_vals, fh)
        print("wrote: %s" % fh.name, file=fh)
    print("ACF:\n", open(prefix + ".acf.txt").read(), file=sys.stderr)

    spvals, opvals = array.array('f'), array.array('f')
    with ts.nopen(prefix + ".slk.bed.gz", "w") as fhslk:
        fhslk.write('#chrom\tstart\tend\tp\tregion-p\n')
        for chrom, results in slk.adjust_pvals(bed_files, col_num, acf_vals):
            fmt = chrom + "\t%i\t%i\t%.4g\t%.4g\n"
            for row in results:
                row = tuple(row)
                fhslk.write(fmt % row)
                opvals.append(row[-2])
                spvals.append(row[-1])

    print("# original lambda: %.2f" % genomic_control(opvals), file=sys.stderr)
    del opvals

    gc_lambda = genomic_control(spvals)
    print("wrote: %s with lambda: %.2f" % (fhslk.name, gc_lambda),
            file=sys.stderr)

    if genome_control:
        fhslk = ts.nopen(prefix + ".slk.gc.bed.gz", "w")
        adj = genome_control_adjust([d['p'] for d in bediter(prefix + ".slk.bed.gz", -1)])
        for i, line in enumerate(ts.nopen(prefix + ".slk.bed.gz")):
            print("%s\t%.5g" % (line.rstrip("\r\n"), adj[i]), file=fhslk)

        fhslk.close()
        print("wrote: %s" % fhslk.name, file=sys.stderr)

    with ts.nopen(prefix + ".fdr.bed.gz", "w") as fh:
        fh.write('#chrom\tstart\tend\tp\tregion-p\tregion-q\n')
        for bh, l in fdr.fdr(fhslk.name, -1):
            fh.write("%s\t%.4g\n" % (l.rstrip("\r\n"), bh))
        print("wrote: %s" % fh.name, file=sys.stderr)
    fregions = prefix + ".regions.bed.gz"
    with ts.nopen(fregions, "w") as fh:
        list(peaks.peaks(prefix + ".fdr.bed.gz", -1 if use_fdr else -2, threshold, seed,
            dist, fh, operator.le))
    n_regions = sum(1 for _ in ts.nopen(fregions))
    print("wrote: %s (%i regions)" % (fregions, n_regions), file=sys.stderr)
    if n_regions == 0:
        sys.exit()

    with ts.nopen(prefix + ".regions-p.bed.gz", "w") as fh:
        N = 0
        fh.write("#chrom\tstart\tend\tmin_p\tn_probes\tz_p\tz_sidak_p\n")
        # use -2 for original, uncorrected p-values in slk.bed
        for region_line, slk_p, slk_sidak_p, sim_p in region_p.region_p(
                               prefix + ".slk.bed.gz",
                               prefix + ".regions.bed.gz", -2,
                               step):
            fh.write("%s\t%.4g\t%.4g\n" % (region_line, slk_p, slk_sidak_p))
            fh.flush()
            N += int(slk_sidak_p < 0.05)
        print("wrote: %s, (regions with corrected-p < 0.05: %i)" \
                % (fh.name, N), file=sys.stderr)

    regions_bed = fh.name
    #if all(h in header for h in ('t', 'start', 'end')):
    if region_filter_n is None: region_filter_n = 0
    with ts.nopen(prefix + ".regions-t.bed", "w") as fh:
        N = 0
        for i, toks in enumerate(filter.filter(bed_files[0],
            regions_bed, p_col_name=col_num)):
            if i == 0: toks[0] = "#" + toks[0]
            else:
                if float(toks[6]) > region_filter_p: continue
                if int(toks[4]) < region_filter_n: continue
                #if region_filter_t and "/" in toks[7]:
                #    # t-pos/t-neg. if the lower one is > region_filter_t?
                #    vals = map(int, toks[7].split("/"))
                #    if min(vals) > region_filter_t: continue

                N += 1
            print("\t".join(toks), file=fh)
        print(("wrote: %s, (regions with region-p "
                            "< %.3f and n-probes >= %i: %i)") \
                % (fh.name, region_filter_p, region_filter_n, N),
                file=sys.stderr)

    try:
        from cpv import manhattan
        regions = manhattan.read_regions(fh.name)

        manhattan.manhattan(prefix + ".slk.bed.gz", 3, prefix.rstrip(".") + ".manhattan.png",
                         False, ['#959899', '#484B4C'], "", False, None,
                         regions=regions, bonferonni=False)
    except ImportError:
        pass # they dont have matplotlib


    if db is not None:
        from cruzdb import Genome
        g = Genome(db)
        lastf = fh.name
        with open(prefix + ".anno.%s.bed" % db, "w") as fh:
            fh.write('#')
            g.annotate(lastf, (table, "cpgIslandExt"), out=fh,
                    feature_strand=True, parallel=len(spvals) > 500)
        print("wrote: %s annotated with %s %s" % (fh.name, db, table), file=sys.stderr)

