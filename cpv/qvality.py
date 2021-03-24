"""
qvality must be on the $PATH.
run on a file with pvalues in column 4:

    $ python qvality.py --column 4 input.bed > output.bed

output file will have extra columns for PEP (posterior error prob)
and q-value.

if input.bed has a header, it is retained.
TODO: handle null via --null [column_number]
"""
from __future__ import print_function
import tempfile
import subprocess as sp
import sys
import bisect
try:
    from itertools import izip
except ImportError:
    izip = zip
from array import array

def qvality(pvals, null=None, **kwargs):
    f = open(tempfile.mktemp(suffix=".qvality"), "w")
    f.write("\n".join(map(str, pvals)))
    f.close()
    if null:
        fnull = open(tempfile.mktemp(suffix=".qvality"), "w")
        fnull.write("\n".join(map(str, null)))
        fnull.close()
    cmd = ['qvality']
    for k, v in kwargs.iteritems():
        flag = ("-%s" if len(k) == 1 else "--%s") % k
        if v is not None: flag += " %s" % v
        cmd.append(flag)

    cmd.append(f.name)
    if null:
        cmd.append(fnull.name)

    print(cmd, file=sys.stderr)

    p = sp.Popen(cmd, stderr=sp.PIPE,
                         stdout=sp.PIPE)
    ps, peps, qs = [array('d') for _ in range(3)]
    for pmax, pep, q in (map(float, l.split("\t"))
                         for i, l in enumerate(p.stdout) if i > 0):
        ps.append(pmax)
        peps.append(pep)
        qs.append(q)
    # make sure they are sorted.
    for i, pset in enumerate((ps, qs)):
        assert all(a >= b for a, b in zip(pset[1:], pset[:-1])), (i, pset[:10], pset[-10:], len(pset))
    p.wait()
    print(p.stderr.read(), file=sys.stderr)
    if p.returncode != 255: # ?
        print(p.stderr.read(), p.returncode, file=sys.stderr)

    for p in pvals:
        # find the index in the ps
        idx = bisect.bisect_left(ps, p)
        # and return the posterior prob, q
        try:
            yield p, peps[idx], qs[idx]
        except IndexError:
            yield p, peps[idx - 1], qs[idx - 1]

def main(f_pvals, column_no):
    # accept negative indexes and adjust to 0-based indexing if positive.
    col = column_no - int(column_no > 0)
    pvals = []
    has_header = False
    for i, toks in enumerate(l.rstrip("\r\n").split("\t")
                                         for l in open(f_pvals)):
        if i == 0:
            try:
                float(toks[col])
            except ValueError:
                has_header = True
                continue
        pvals.append(float(toks[col]))

    assert all(0 <= p <= 1 for p in pvals)

    fh = open(f_pvals)
    if has_header:
        print( "%s\tPEP\tq-value" % (fh.readline().rstrip("\r\n"), ))

    for (p, pep, q), line in  izip(qvality(pvals), fh):
        toks = line.rstrip("\r\n").split("\t")
        assert float(toks[col]) == p, (line, p)
        print("%s\t%.6g\t%.6g" % (line.rstrip("\r\n"), pep, q))



if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--column", default=-1,
                   help="column number the file containing p-values. default \
                   is last column", type=int)
    p.add_argument("pvals")

    args = p.parse_args()

    main(args.pvals, args.column)
