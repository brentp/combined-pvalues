import toolshed as ts
from itertools import tee, izip
import signal
import sys

def get_col_nums(c):
    """
    >>> get_col_nums(4)
    [3]
    >>> get_col_nums("4,5")
    [3, 4]
    >>> get_col_nums("4,-1")
    [3, -1]
    """
    if isinstance(c, int): return [get_col_num(c)]
    return [get_col_num(int(n)) for n in c.rstrip().split(",")]

def get_col_num(c, bed_file=None):
    """
    adjust the col number so it does intutive stuff
    for command-line interface
    >>> get_col_num(4)
    3
    >>> get_col_num(-1)
    -1
    """
    if isinstance(c, basestring) and c.isdigit():
        c = int(c)
    if isinstance(c, (int, long)):
        return c if c < 0 else (c - 1)
    header = ts.header(bed_file)
    assert c in header
    return header.index(c)


def bediter(fnames, col_num, delta=None):
    """
    iterate over a bed file. turn col_num into a float
    and the start, stop column into an int and yield a dict
    for each row.
    """
    last_chrom = chr(0)
    last_start = -1
    if isinstance(fnames, basestring):
        fnames = [fnames]
    for fname in fnames:
        for i, l in enumerate(ts.reader(fname, header=False)):
            if l[0][0] == "#": continue
            if i == 0: # allow skipping header
                try:
                    float(l[col_num])
                except ValueError:
                    continue
            chrom = l[0]
            start = int(float(l[1]))
            if chrom == last_chrom:
                assert start >= last_start, ("error at line: %i, %s"
                        % (i, "\t".join(l)), "file is not sorted")
            else:
                assert last_chrom < chrom, ("error at line: %i, %s "
                        " with file: %s" % (i, "\t".join(l), fname),
                        "chromosomes must be sorted as characters",
                        last_chrom, "is not < ", chrom)
                last_chrom = chrom

            last_start = start

            p = float(l[col_num])
            if not delta is None:
                if p > 1 - delta: p-= delta # the stouffer correction doesnt like values == 1
                if p < delta: p = delta # the stouffer correction doesnt like values == 0

            v = {"chrom": l[0], "start": start, "end": int(float(l[2])),
                 "p": p} # "stuff": l[3:][:]}
            if v['end'] - v['start'] > 100000:
                print("warning! large interval at %s will increase memory use." % v)
            yield v

def genomic_control(pvals):
    """
    calculate genomic control factor, lambda
    >>> genomic_control([0.25, 0.5, 0.75])
    1.0000800684096998
    >>> genomic_control([0.025, 0.005, 0.0075])
    15.715846578113579
    """
    from scipy import stats
    import numpy as np
    pvals = np.asarray(pvals)
    return np.median(stats.chi2.ppf(1 - pvals, 1)) / 0.4549

def genome_control_adjust(pvals):
    """
    adjust p-values by the genomic control factor, lambda
    >>> genome_control_adjust([0.4, 0.01, 0.02])
    array([ 0.8072264 ,  0.45518836,  0.50001716])
    """
    import numpy as np
    from scipy import stats
    pvals = np.asarray(pvals)
    qchi = stats.chi2.ppf(1 - pvals, 1)
    gc = np.median(qchi) / 0.4549
    return 1 - stats.chi2.cdf(qchi / gc, 1)

def genome_control_adjust_bed(bedfiles, colnum, outfh):
    c = colnum
    adj = genome_control_adjust([d['p'] for d in bediter(bedfiles, colnum)])

    diff = 0
    if len(bedfiles) > 1:
        print >>sys.stderr, "can't do genomic control adjustment with more than 1 bed file"
        sys.exit(4)
    for j, bedfile in enumerate(bedfiles):
        for i, toks in enumerate(line.rstrip("\r\n").split("\t") \
                for line in ts.nopen(bedfile)):
            try:
                float(toks[c])
            except ValueError: # headder
                if i == 0 == j:
                    print >>outfh, "\t".join(toks)
                    diff = 1
                    continue
                elif i == 0:
                    continue
                else:
                    raise
            toks[c] = "%.5g" % adj[i - diff]
            print >>outfh, "\t".join(toks)

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def read_acf(acf_file):
    acf_vals = {}
    for row in ts.nopen(acf_file):
        if row[0] == "#": continue
        row = row.split("\t")
        if row[0] == "lag_min": continue
        acf_vals[(int(row[0]), int(row[1]))] = float(row[2])
    return sorted(acf_vals.items())

def pool_sig():
    return signal.signal(signal.SIGINT, signal.SIG_IGN)

# from aljunberg:  https://gist.github.com/aljungberg/626518 
from multiprocessing.pool import IMapIterator
def wrapper(func):
    def wrap(self, timeout=None):
        return func(self, timeout=timeout or 1e100)
    return wrap
IMapIterator.next = wrapper(IMapIterator.next)

def get_map():
    try:
        from multiprocessing import Pool
        p = Pool(None, pool_sig)
        imap = p.imap
    except ImportError:
        from itertools import imap
    return imap

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
