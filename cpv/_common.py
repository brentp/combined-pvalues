from toolshed import reader
from itertools import tee, izip

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

def get_col_num(c):
    """
    adjust the col number so it does intutive stuff
    for command-line interface
    >>> get_col_num(4)
    3
    >>> get_col_num(-1)
    -1
    """
    return c if c < 0 else (c - 1)

def bediter(fname, col_num):
    """
    iterate over a bed file. turn col_num into a float
    and the start, stop column into an int and yield a dict
    for each row.
    """
    for l in reader(fname, header=False):
        if l[0][0] == "#": continue
        p = float(l[col_num])
        if p == 1: p-= 1e-8 # the stouffer correction doesnt like values == 1
        yield  {"chrom": l[0], "start": int(l[1]), "end": int(l[2]),
                "p": p} # "stuff": l[3:][:]}

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def read_acf(acf_file):
    acf_vals = {}
    for row in open(acf_file):
        if row[0] == "#": continue
        row = row.split("\t")
        if row[0] == "lag_min": continue
        acf_vals[(int(row[0]), int(row[1]))] = float(row[2])
    return sorted(acf_vals.items())

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
