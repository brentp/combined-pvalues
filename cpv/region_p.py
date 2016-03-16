"""
   calculate a p-value of a region using the Stouffer-Liptak method or the
   z-score method.
"""
import argparse
import sys
import numpy as np
import toolshed as ts
from collections import defaultdict

from interlap import InterLap

from _common import bediter, get_col_num
from itertools import chain, groupby
from operator import itemgetter
from slk import gen_sigma_matrix
from acf import acf
from stouffer_liptak import stouffer_liptak, z_score_combine

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
            col_num, args.step, z=True):
        #if sim_p != "NA":
        #    sim_p = "%.4g" % (sim_p)
        print "%s\t%.4g\t%.4g" % (region_line, slk, slk_sidak)

def _gen_acf(region_info, fpvals, col_num, step):
    # calculate the ACF as far out as needed...
    # keys of region_info are (chrom, start, end)
    max_len = max(int(r[2]) - int(r[1]) for r in region_info)
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
    acfs = acf(fpvals, lags, col_num, simple=True)
    print >>sys.stderr, "# Done with one-time ACF calculation"
    return acfs

def get_total_coverage(fpvals, col_num, step, out_val):
    """
    Calculate total bases of coverage in `fpvals`.
    Used for the sidak correction
    """
    total_coverage = 0
    for key, chrom_iter in groupby(bediter(fpvals, col_num),
            itemgetter('chrom')):
        bases = set([])
        for feat in chrom_iter:
            s, e = feat['start'], feat['end']
            if s == e: e += 1
            #e = max(e, s + step)
            bases.update(range(s, e))
        total_coverage += len(bases)
    out_val.value = total_coverage

def _get_total_coverage(fpvals, col_num, step):
    from multiprocessing import Process, Value
    val = Value('f')
    p = Process(target=get_total_coverage, args=(fpvals, col_num, step, val))
    p.start()
    return p, val

def sidak(p, region_length, total_coverage, message=[False]):
    """
    see: https://github.com/brentp/combined-pvalues/issues/2
    """
    if region_length == 0:
        region_length = 1
        if not message[0]:
            message[0] = True
            sys.stderr.write(""""warning: 0-length region found.
does input have 0-length intervals? using length of 1 and not reporting
further 0-length intervals""")
    # use 1.1 as heuristic to account for limit in available regions
    # of a given size as the region_length increases
    # TODO: base that on the actual number of regiosn of this length
    # that could be seen based on the distance constraint.
    k = total_coverage / (np.float64(region_length)**1.0)
    if k < 1: k = total_coverage
    p_sidak = 1 - (1 - p)**k
    if p_sidak == 0:
        assert p < 1e-16, (p, k, total_coverage, region_length)
        p_sidak = (1 - (1 - 1e-16)**k) / (p / 1e-16)
        p_sidak = min(p_sidak, p * k)

    # print "bonferroni:", min(p * k, 1)
    return min(p_sidak, 1)

def _get_ps_in_regions(tree, fpvals, col_num):
    """
    find the pvalues associated with each region
    """
    region_info = defaultdict(list)
    for row in bediter(fpvals, col_num):
        for region in tree[row['chrom']].find((row['start'], row['end'])):
            region_len = max(1, region[1] - region[0])
            region_tup = tuple(region[-1])
            region_info[region_tup].append(row)
    assert sum(len(v) for v in tree.values()) >= len(region_info)
    if sum(len(v) for v in tree.values()) > len(region_info):
        sys.stderr.write("# note: not all regions contained measurements\n")
    return region_info

def read_regions(fregions):
    tree = defaultdict(InterLap)
    for i, toks in enumerate(ts.reader(fregions, header=False)):
        if i == 0 and not (toks[1] + toks[2]).isdigit(): continue
        tree[toks[0]].add((int(toks[1]), int(toks[2]), toks))
    sys.stderr.write("# read %i regions from %s\n" \
            % (sum(len(v) for v in tree.values()), fregions))
    return tree

def region_p(fpvals, fregions, col_num, step, z=True):
    # just use 2 for col_num, but dont need the p from regions.

    tree = read_regions(fregions)
    process, total_coverage_sync = _get_total_coverage(fpvals, col_num, step)

    region_info = _get_ps_in_regions(tree, fpvals, col_num)

    acfs = _gen_acf(region_info, (fpvals,), col_num, step)
    process.join()
    total_coverage = total_coverage_sync.value

    # regions first and then create ACF for the longest one.
    print >>sys.stderr, "%i bases used as coverage for sidak correction" % \
                                (total_coverage)
    sample_distribution = np.array([b["p"] for b in bediter(fpvals,
                                                                col_num)])

    combine = z_score_combine if z else stouffer_liptak
    for region, prows in region_info.iteritems():
        # gen_sigma expects a list of bed dicts.
        sigma = gen_sigma_matrix(prows, acfs)
        ps = np.array([prow["p"] for prow in prows])
        if ps.shape[0] == 0:
            print >>sys.stderr,("bad region", region)
            continue

        # calculate the SLK for the region.
        region_slk = combine(ps, sigma)
        if not region_slk["OK"]:
            print >>sys.stderr, "problem with:", region_slk, ps

        slk_p = region_slk["p"]

        sidak_slk_p = sidak(slk_p, int(region[2]) - int(region[1]), total_coverage)

        result = ["\t".join(region), slk_p, sidak_slk_p, "NA"]
        yield result

def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", dest="pvals", help="BED containing all the p values"
                  " used to generate `regions`")
    p.add_argument("-r", dest="regions", help="BED containing all the regions")
    p.add_argument("-s", "--step", dest="step", type=int, default=50,
            help="step size for acf calculation. should be the same "
            " value as the step sent to -d arg for acf")
    p.add_argument("-c", dest="c", help="column number containing the p-value"
                   " of interest", type=str, default=-1)
    p.add_argument("-z", dest="z", help="use z-score correction",
                    action="store_true")
    args = p.parse_args()
    if not (args.regions and args.pvals):
        import sys
        sys.exit(not p.print_help())
    header = ts.nopen(args.regions).next()
    if header.startswith("#") or (not header.split("\t")[2].isdigit()):
        print "%s\tslk_p\tslk_sidak_p" % (header.rstrip("\r\n"),)

    header = ts.header(args.pvals)
    if args.c in header:
        args.c = header.index(args.c) + 1
    else:
        args.c = int(args.c)
    return run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
