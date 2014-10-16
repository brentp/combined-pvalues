"""
count the number of switches in sign in the regions. Since the region
calculation is based on the p-value only, it could be that a region is
discovered that has both high and low t-scores.
This script will output the original region_bed intervals, along with
sum of positive t-scores and the sum of negative t-scores.
"""
import argparse
from operator import itemgetter
from itertools import groupby
from tempfile import mktemp
from math import exp
import atexit
import os

import toolshed as ts

def ilogit(v):
    return 1 / (1 + exp(-v))

def fix_bed(fname):
    """
    a lot of bed files will have no header or have e.g.
    8e6 instead of 8000000 for start/end. this just fixes that
    so we can send to bedtools
    """
    r = ts.reader(fname, header=False)
    h = next(r)
    assert not (h[1] + h[2]).isdigit(), "must have header for filtering"
    tname = mktemp()
    fh = ts.nopen(tname, "w")
    fh.write("#" + "\t".join(h) + "\n")
    for toks in r:
        toks[1:3] = map(str, (int(float(t)) for t in toks[1:3]))
        fh.write("%s\n" % "\t".join(toks))
    fh.close()
    atexit.register(os.unlink, tname)
    return tname


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", dest="p", help="p-value column from `p_bed`",
            default="P.Value")
    p.add_argument("-t", dest="t", default="t",
            help="t-statistic or directionality from p_bed")
    p.add_argument("--coef", default="logFC",
            help="name of coefficient column in BED")
    p.add_argument("--filter", help="don't print row if there's a swith in \
            t-scores", action="store_true", default=False)
    p.add_argument("--max-p", help="filter regions with any p-value > this value", type=float)
    p.add_argument("--region-p", help="filter regions with combined p-value > this value", type=float)
    p.add_argument('region_bed', help='file containing the regions')
    p.add_argument('p_bed', help='file containing the raw p-values')
    args = p.parse_args()

    for row in filter(args.p_bed, args.region_bed, args.max_p,
            region_p=args.region_p,
            p_col_name=args.p,
            coef_col_name=args.coef):
        print "\t".join(row)

def filter(p_bed, region_bed, max_p=None, region_p=None, p_col_name="P.Value",
                    coef_col_name="logFC"):

    ph = ts.header(p_bed)
    if (ph[1] + ph[2]).isdigit():
        raise Exception('need header in p-value file to run filter')
    assert ph[1] == 'start' and ph[2] == 'end' and ph[0] == 'chrom', \
            ('must have chrom, start, end header for', p_bed)
    ph = ['p' + h for h in ph]

    rh = ts.header(region_bed)
    header = not (rh[1] + rh[2]).isdigit()

    if isinstance(p_col_name, str) and p_col_name.isdigit():
        p_col_name = int(p_col_name) - 1

    if isinstance(p_col_name, (int, long)):
        p_col_name = ph[p_col_name][1:]

    a = dict(p_bed=p_bed, region_bed=region_bed)
    a['p_bed'] = fix_bed(a['p_bed'])
    a['header'] = ""

    j = 0
    for group, plist in groupby(
            ts.reader('|bedtools intersect -b %(p_bed)s \
                         -a %(region_bed)s -wo %(header)s' % a,
            header=rh + ph), itemgetter('chrom','start','end')):
        plist = list(plist)

        if region_p:
            r = plist[0] # first cols are all the same
            region_p_key = 'slk_sidak_p' if 'slk_sidak_p' in r \
                                         else 'z_sidak_p' if 'z_sidak_p' in r \
                                         else None
            if region_p_key is None: raise Exception
            if float(r[region_p_key]) > region_p:
                continue

        try:
            plist = [x for x in plist if (int(x['start']) <= int(x['pstart']) <= int(x['pend'])) and ((int(x['start']) <= int(x['pend']) <= int(x['end'])))]
        except:
            print plist
            raise
        tscores = [float(row['pt']) for row in plist if 'pt' in row]

        if max_p:
            if any(float(row['p' + p_col_name]) > max_p for row in plist):
                continue

        ngt05  = sum(1 for row in plist if float(row['p' + p_col_name]) > 0.05)

        # logic to try to find t and coef headers and skip if not found
        extra_header = []
        extra = []
        if tscores:
            tpos = sum(1 for ts in tscores if ts > 0)
            tneg = sum(1 for ts in tscores if ts < 0)
            tpn = "%i/%i" % (tpos, tneg)

            tsum = str(sum(ts for ts in tscores))
            extra_header += ["t.pos/t.neg", "t.sum"]
            extra += [tpn, tsum]
        else:
            tsum = tpn = "NA"
        if 'p' + coef_col_name not in plist[0] and 'pcoefficient' in plist[0]:
            coef_col_name = 'coefficient'
        if 'p' + coef_col_name in plist[0]:
            coef = (sum(float(row['p' + coef_col_name]) for row in plist) /
                                    len(plist))

            # since we probably had the data logit transformed, here we
            # do the inverse and subtract 0.5 since ilogit(0) == 0.5
            icoef = (sum(ilogit(float(row['p' + coef_col_name])) for row in plist) /
                                    len(plist)) - 0.5
            extra_header += ["avg.diff", "ilogit.diff"]
            extra += ["%.3f" % coef, "%.3f" % icoef]
        else:
            coef = icoef = float('nan')

        frow = [plist[0][h] for h in rh] + extra
        if j == 0:
            yield rh + extra_header
            j = 1
        yield frow

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
