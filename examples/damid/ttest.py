import sys
from scipy.stats import ttest_ind as ttest
a = [float(x) for x in open(sys.argv[1])]
b = [float(x) for x in open(sys.argv[2])]

print ttest(a, b)
