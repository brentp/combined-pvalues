Here we show a dataset with some regions of low p-values for which a
traditional approach would find no p-values.
For example, after Benjamini-Hochberg, there are no regions with a q-value
less that 0.05:

```Shell
$ comb-p fdr -c 4 fake/pvalues.bed | awk '$5 < 0.05' | wc -l
0
```

However, using the `comb-p` pipeline, we can find 34 significant regions.

```Shell
comb-p pipeline --seed 0.01 --dist 300 --step 300 pvalues.bed -p fake
wrote: fake.acf.txt
ACF:
# https://chart.googleapis.com/chart?cht=bvs&&chd=t:0.059&chs=74x250&chco=224499&chxt=x,y&chxl=0:|1-301&chxr=1,0,0.08&chm=N,000000,0,-1,12&chbh=42,6,12&chds=a
#lag_min    lag_max correlation N   p
1   301 0.05881 544300  0

wrote: fake.slk.bed
wrote: fake.fdr.bed
wrote: fake.regions.bed (75 regions)
# calculating ACF out to: 497
#           with 3  lags: [1, 301, 601]
# Done with one-time ACF calculation
21789882 bases used as coverage for sidak correction
wrote: fake.regions-p.bed, (regions with corrected-p < 0.05: 34)
```
