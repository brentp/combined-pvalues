"""
   simulate a p-value of a region using:

    + The Stouffer-Liptak method.

    + the truncated product method described
      in Zaykin et al. 2002, "Truncated Product Method for Combining p-values".
      Genet Epidemiol.
"""
import argparse
import sys
import numpy as np

from _common import bediter, get_col_num
from itertools import chain, groupby
from operator import itemgetter
from slk import gen_sigma_matrix
from acf import acf
from stouffer_liptak import stouffer_liptak

from scipy.stats import norm
from numpy.linalg import cholesky as chol
qnorm = norm.ppf
pnorm = norm.cdf

def gen_correlated(sigma, n, observed=None):
    """
    generate autocorrelated data according to the matrix
    sigma. if X is None, then data will be sampled from
    the uniform distibution. Otherwise, it will be sampled
    from X. Where X is then *all* observed
    p-values.
    """
    C = np.matrix(chol(sigma))
    if observed is None:
        X = np.random.uniform(0, 1, size=(n, sigma.shape[0]))
    else:
        assert n * sigma.shape[0] < observed.shape[0]
        idxs = np.random.random_integers(0, len(observed) - 1,
                                         size=sigma.shape[0] * n)
        X = observed[idxs].reshape((n, sigma.shape[0]))

    Q = np.matrix(qnorm(X))
    for row in  np.array(1 - norm.sf((Q * C).T)).T:
        yield row

def sl_sim(sigma, ps, nsims, sample_distribution=None):
    N = 0
    print >>sys.stderr, "nsims:", nsims
    w0 = stouffer_liptak(ps, sigma)["p"]
    # TODO parallelize here.
    for i in range(10):
        for prow in gen_correlated(sigma, nsims/10, sample_distribution):
            s = stouffer_liptak(prow, sigma)
            if not s["OK"]: 1/0
            if s["p"] <= w0: N += 1

    return N / float(nsims)

def run(args):
    col_num = get_col_num(args.c)
    # order in results is slk, uniform, sample
    for region_line, slk, slk_sidak, sim_p in region_p(args.pvals, args.regions,
            col_num, args.N, args.step, args.mlog):
        #if sim_p != "NA":
        #    sim_p = "%.4g" % (sim_p)
        print "%s\t%.4g\t%.4g" % (region_line, slk, slk_sidak)

def _gen_acf(region_info, fpvals, col_num, step, mlog):
    # calculate the ACF as far out as needed...
    # [1] is the length of each region.
    max_len = max(r[1] for r in region_info)
    print >>sys.stderr, "# calculating ACF out to: %i" % max_len

    lags = range(1, max_len, step)
    if lags[-1] < max_len: lags.append(lags[-1] + step)
    if len(lags) > 20:
        repr_lags = "[" + ", ".join(map(str, lags[1:4])) + \
                    " ... " + \
                    ", ".join(map(str, lags[-5:])) + "]"
    else:
        repr_lags = str(lags)
    print >>sys.stderr, "#           with %-2i lags: %s" \
            % (len(lags), repr_lags)

    if len(lags) > 100:
        print >>sys.stderr, "# !! this could take a looong time"
        print >>sys.stderr, "# !!!! consider using a larger step size (-s)"
    acfs = acf(fpvals, lags, col_num, simple=True, mlog=mlog)
    print >>sys.stderr, "# Done with one-time ACF calculation"
    return acfs

def get_total_coverage(fpvals, col_num, out_val):
    """
    Calculate total bases of coverage in `fpvals`.
    Used for the sidak correction
    """
    total_coverage = 0
    for key, chrom_iter in groupby(bediter(fpvals, col_num),
            itemgetter('chrom')):
        bases = set([])
        for feat in chrom_iter:
            bases.update(range(feat['start'], feat['end']))
        total_coverage += len(bases)
    out_val.value = total_coverage

def _get_total_coverage(fpvals, col_num):
    from multiprocessing import Process, Value
    val = Value('f')
    p = Process(target=get_total_coverage, args=(fpvals, col_num, val))
    p.start()
    return p, val

def sidak(p, region_length, total_coverage):
    """
    see: https://github.com/brentp/combined-pvalues/issues/2
    """
    k = total_coverage / float(region_length)
    p_sidak = 1 - (1 - p)**k
    # print "bonferroni:", min(p * k, 1)
    return min(p_sidak, 1)

def gen_regions(fregions):
    for region_line in (l.rstrip("\r\n")
                                   for l in open(fregions) if l[0] != "#"):
        toks = region_line.split("\t")
        rchrom = toks[0]
        rstart, rend = map(int, toks[1:3])
        yield rchrom, rstart, rend, region_line

def _get_ps_in_regions(fregions, fpvals, col_num):
    """
    find the pvalues associated with each region
    """
    region_info = []
    piter = chain(bediter(fpvals, col_num), [None])
    prow = piter.next()
    nr = 0
    for rchrom, rstart, rend, region_line in sorted(gen_regions(fregions),
                                                key=itemgetter(0, 1)):
        prows = []
        nr += 1
        # grab the p-values in the bed file that are within the current region
        while (prow["chrom"] != rchrom or prow["start"] < rstart):
        #while (prow["chrom"] != rchrom or prow["end"] < rstart):
            prow = piter.next()
            if prow is None: break

        #while prow is not None and (rchrom, rend) > (prow["chrom"], prow["start"]):
        while prow is not None and (rchrom, rend) >= (prow["chrom"], prow["end"]):
            prows.append(prow)
            prow = piter.next()
            if prow is None: break
        if not prows:
            print >>sys.stderr, "missed,:", prows, (region_line)
        region_len = rend - rstart
        region_info.append((region_line, region_len, prows[:]))
        del prows
    assert nr == len(region_info), (nr, len(region_info))
    return region_info

def region_p(fpvals, fregions, col_num, nsims, step, mlog=False):
    # just use 2 for col_num, but dont need the p from regions.

    if(sum(1 for _ in open(fregions) if _[0] != "#") == 0):
        print >>sys.stderr, "no regions in %s" % (fregions, )
        sys.exit()

    process, total_coverage_sync = _get_total_coverage(fpvals, col_num)
    region_info = _get_ps_in_regions(fregions, fpvals, col_num)

    acfs = _gen_acf(region_info, (fpvals,), col_num, step, mlog=mlog)
    process.join()
    total_coverage = total_coverage_sync.value

    # regions first and then create ACF for the longest one.
    print >>sys.stderr, "%i bases used as coverage for sidak correction" % \
                                (total_coverage)
    sample_distribution = np.array([b["p"] for b in bediter(fpvals,
                                                                col_num)])
    for region_line, region_len, prows in region_info:
        # gen_sigma expects a list of bed dicts.
        sigma = gen_sigma_matrix(prows, acfs)
        ps = np.array([prow["p"] for prow in prows])
        if ps.shape[0] == 0:
            print >>sys.stderr,("bad region", region_line)
            continue
        #assert ps.shape[0] != 0, ("bad region", region_line)
        #assert ps.shape[0] == sigma.shape[0], ("bad_region", region_line)

        # calculate the SLK for the region.

        region_slk = stouffer_liptak(ps, sigma)

        assert region_slk["OK"] is True
        slk_p = region_slk["p"]

        sidak_slk_p = sidak(slk_p, region_len, total_coverage)

        result = [region_line, slk_p, sidak_slk_p]

        # corroborate those with p-values < 0.1 by simulation
        #"""
        if nsims > 0:

            # adjust nsims so it's an adjusted p-value.
            q_nsims = int(0.5 + total_coverage / float(region_len))
            assert sample_distribution is not None
            # trim sigma because we may have trimmed the ps above.
            sim_p = sl_sim(sigma, ps, q_nsims, sample_distribution)
            result.append(sim_p)
        else:
            result.append("NA")
        #"""
        #result.append("NA")
        yield result

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", dest="pvals", help="BED containing all the p values"
                  " used to generate `regions`")
    p.add_argument("-r", dest="regions", help="BED containing all the regions")
    p.add_argument("-s", dest="step", type=int, default=50,
            help="step size for acf calculation. should be the same "
            " value as the step sent to -d arg for acf")
    p.add_argument("--mlog", dest="mlog", action="store_true",
                   default=False, help="do the correlation on the -log10 of"
                   "the p-values. Default is to do it on the raw values")
    p.add_argument("-N", dest="N", help="number of simulations to perform",
                   type=int, default=0)
    p.add_argument("-c", dest="c", help="column number containing the p-value"
                   " of interest", type=int, default=-1)
    args = p.parse_args()
    if not (args.regions and args.pvals):
        import sys
        sys.exit(not p.print_help())
    from toolshed import nopen
    header = nopen(args.regions).next()
    if header.startswith("#") or (not header.split("\t")[2].isdigit()):
        print "%s\tslk_p\tslk_sidak_p" % (header.rstrip("\r\n"),)
    return run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
