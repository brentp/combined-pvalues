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

    ph = ['p' + h for h in get_header(args.p_bed)]
    rh = get_header(args.region_bed)

    a = dict(p_bed=args.p_bed, region_bed=args.region_bed)
    print "#" + "\t".join(rh + ["t-pos", "t-neg"])
    for group, plist in groupby(reader('|bedtools intersect -b %(p_bed)s -a %(region_bed)s -wo' % a,
            header=rh + ph), itemgetter('chrom','start','end')):
        plist = list(plist)
        tscores = [float(row['p' + args.t]) for row in plist]

        tpos = sum(1 for ts in tscores if ts > 0)
        tneg = sum(1 for ts in tscores if ts < 0)
        if args.filter and (tpos * tneg != 0):
            continue
        if args.max_p:
            if any(float(row['p' + args.p]) > args.max_p for row in plist):
                continue

        frow = [plist[0][h] for h in rh] + [str(tpos), str(tneg)]
        print "\t".join(frow)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
