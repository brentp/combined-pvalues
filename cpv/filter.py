"""
count the number of switches in sign in the regions. Since the region
calculation is based on the p-value only, it could be that a region is
discovered that has both high and low t-scores.
This script will output the original region_bed intervals, along with
sum of positive t-scores and the sum of negative t-scores.
"""
import argparse
from toolshed import reader, header as get_header
from operator import itemgetter
from itertools import groupby
from tempfile import mktemp
from math import exp

def ilogit(v):
    return 1 / (1 + exp(-v))

def fix_header(fname):
    r = reader(fname, header=False)
    h = r.next()
    if not h[0].startswith(("#", )) and (h[1] + h[2]).isdigit():
        return fname
    tname = mktemp()
    fh = open(tname, "w")
    print >>fh, "#" + "\t".join(h)
    for toks in r:
        print >>fh, "\t".join(toks)
    fh.close()
    return tname


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", dest="p", help="p-value column from region_bed",
            default="P.Value")
    p.add_argument("-t", dest="t", default="t",
            help="t-statistic or directionality from p_bed")
    p.add_argument("--filter", help="don't print row if there's a swith in \
            t-scores", action="store_true", default=False)
    p.add_argument("--max-p", help="filter regions with any p-value > this value", type=float)
    p.add_argument('region_bed', help='file containing the regions')
    p.add_argument('p_bed', help='file containing the raw p-values')
    args = p.parse_args()

    for row in filter(args.p_bed, args.region_bed, args.max_p):
        print "\t".join(row)

def filter(p_bed, region_bed, max_p=None, p_col_name="P.Value"):
    ph = ['p' + h for h in get_header(p_bed)]
    rh = get_header(region_bed)
    if isinstance(p_col_name, (int, long)):
        p_col_name = ph[p_col_name][1:]

    a = dict(p_bed=p_bed, region_bed=region_bed)
    a['p_bed'] = fix_header(a['p_bed'])

    yield rh + ["t.pos/t.neg", "t.sum", "n.gt.p05", 'avg.diff', 'ilogit.diff']
    for group, plist in groupby(reader('|bedtools intersect -b %(p_bed)s -a %(region_bed)s -wo' % a,
            header=rh + ph), itemgetter('chrom','start','end')):
        plist = list(plist)
        plist = [x for x in plist if (int(x['start']) <= int(x['pstart']) <= int(x['pend'])) and ((int(x['start']) <= int(x['pend']) <= int(x['end'])))]
        tscores = [float(row['pt']) for row in plist if 'pt' in row]

        if max_p:
            if any(float(row['p' + p_col_name]) > max_p for row in plist):
                continue

        ngt05  = sum(1 for row in plist if float(row['p' + p_col_name]) > 0.05)
        tpos = sum(1 for ts in tscores if ts > 0)
        tneg = sum(1 for ts in tscores if ts < 0)
        tpn = "%i/%i" % (tpos, tneg)
        tsum = sum(ts for ts in tscores)
        coef = (sum(float(row['plogFC']) for row in plist) /
                                    len(plist))
        # since we probably had the data logit transformed, here we
        # do the inverse and subtract 0.5 since ilogit(0) == 0.5
        icoef = (sum(ilogit(float(row['plogFC'])) for row in plist) /
                                    len(plist)) - 0.5
        frow = [plist[0][h] for h in rh] + [tpn, str(tsum),
                str(ngt05), '%.3f' % coef, '%.3f' % icoef]
        yield frow

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
