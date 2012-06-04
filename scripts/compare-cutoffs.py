"""
To calculate a true false-discovery rate for a given set
of data, one can run the entire pipeline on the original
data, and on a set of data on which the clinical data has
been shuffled.

This script accepts the regions-p.bed files from the original
and shuffled data respectively and prints the `n_probes` and
`slk_sidak_p` cutoffs that minimize the false discovery rate.

These cutoffs can be used to filter the putative DMRs for those
in down-stream analyses.
"""
from toolshed import reader


prange = [1e-5, 2e-5, 1e-4, 2e-4, 1e-3, 2e-3, 1e-2, 0.05, 0.1]
nrange = range(2, 9)

actual = [(float(d['slk_sidak_p']), int(d['n_probes']))
                                        for d in reader(1)]

shuff = [(float(d['slk_sidak_p']), int(d['n_probes']))
                                        for d in reader(2)]

best_rates = []
for p_cutoff in prange:
    for n_cutoff in nrange:
        n_actual = sum(1 for a in actual if a[0] <= p_cutoff and a[1] >= n_cutoff)
        n_shuff = sum(1 for s in shuff if s[0] <= p_cutoff and s[1] >= n_cutoff)
        best_rates.append((float(n_shuff) / n_actual, p_cutoff, n_cutoff, n_actual, n_shuff))

best_rates.sort()
n = 0
print "#fdr\tslk_sidak_p\tn_probes\tn_actual\tn_shuff"
for i, b in enumerate(best_rates[:-1]):
    bnext = best_rates[i + 1]
    if b[0] == bnext[0] and b[3] < bnext[3]: continue
    if b[0] == bnext[0] and b[2] < bnext[2]: continue
    n += 1
    print "\t".join(map(str, b))
    if n > 10: break
