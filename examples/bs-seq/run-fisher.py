import sys
import fisher
from pybedtools import BedTool
"""
chr1	108	109	16	0	CG+
chr1	109	110	56	4	CG-
chr1	114	115	19	1	CG+
chr1	115	116	57	0	CG-
chr1	160	161	7	11	CG+
chr1	161	162	14	7	CG-
chr1	309	310	16	2	CG+
chr1	310	311	7	1	CG-
chr1	499	500	22	7	CG+
chr1	500	501	8	0	CG-
"""

for b in BedTool(sys.argv[1]).intersect(sys.argv[2], stream=True, wo=True):
    p = fisher.pvalue(int(b[3]), int(b[4]), int(b[9]), int(b[10]))
    print "\t".join((b.chrom, b[1], b[2], ("%.4g" % p.two_tail)))
