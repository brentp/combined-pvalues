Here we present an example using CHARM (Comprehensive high-throughput arrays
for relative methylation) for measuring methylatin using a microarray.

We use data from Irizarry **et al.** Nature Genetics 41, 178 - 186 (2009).
It contains data for different tissues and for tumor vs normal in colon
tissue. Where the original authors use a t-test approach in the R CHARM
package (http://bioconductor.org/packages/release/bioc/html/charm.html),
to do pairwise tests, we use fit a linear model of the form:

```R
model = lm(methp ~ DiseaseState + TissueType, data=df)
```

in the R programming language before sending to comb-p.
It's also possible to find tissue-specific DMR's more specific
than DMR's that vary among the tissues.

```R
model = lm(methp ~ DiseaseState + colon + frontalcortex + liver + spleen)
```

However, we do not report these per-tissue DMR's.

We fit the first model at each of the 2.1 million probes of the CHARM array.
Here, for the purpose of validation and demonstration, we find disease-specific 
DMRs (from DiseaseState) and compare to Irizarry et al's c-dmrs and we find
tissue-specific DMR's and compare to Irizarry et al's t-dmrs.

The full pipeline can be found in:
https://github.com/brentp/combined-pvalues/blob/master/examples/charm/run.sh
with the R code in:
https://github.com/brentp/combined-pvalues/blob/master/examples/charm/scripts/fit.lm.R

Disease-Specific c-dmrs
-----------------------

To compare to Irizarry, we filter `comb-p` candidate DMR's to those
with a p-values less than 0.1 and with more than 8 probes in the DMR.

```Shell
awk '$7 < 0.1 || $5 > 8' data/quantile/p.disease/p.disease.regions-p.bed \
    > comb-p.filtered.bed
```
This gives 3105 DMR's compared to Irizarry's 2707. Of those 2289 overlap,
meaning 85% of Irizarry's DMR's are represented in the `comb-p` set.

```Shell
 $ intersectBed -a comb-p.filtered.bed -b data/c-dmrs.txt | wc -l
 2289
```
The differences are likely due to the difference in methods where Irizarry use
a t-test directly and we pull the t-score from the linear model.

Tissue-Specific t-dmrs
----------------------

We again filter our tissue-specific DMR's to get 11,061 to compare to
Irizarry's 16,379.
This time we filter to a corrected p-value of 0.05 and require at least
8 probes in the DMR.

```Shell
 $ awk '$7 < 0.05 && $5 > 7' data/quantile/p.tissue/p.tissue.regions-p.bed \
     > tissue.filtered.bed

 $ grep -cv "#" tissue.filtered.bed 
 11061
```
Of the `comb-p` DMRs, 8275 or 75% appear in Irizarry's

```Shell
$ intersectBed -a tissue.filtered.bed -b data/t-dmrs.txt  | wc -l
8275
```

For tissue-specific DMR's Irizarry reports an FDR q-value. We can
compare this between his DMRs that overlap our DMR's and those that do
not.

```Shell
$ intersectBed -a data/t-dmrs.txt -b tissue.filtered.bed -v \
      | cut -f 4 > overlap.txt                    
$ intersectBed -a data/t-dmrs.txt -b tissue.filtered.bed  \
      | cut -f 4 > no-overlap.txt

$ python ../damid/ttest.py overlap.txt no-overlap.txt 
(25.948028554040935, 1.2957717193202179e-145)
```
So the p-value of 1.296e-145 indicates that we can enrich for
lower q-values from the Irizarry data by intersecting with our
DMR's. So the Irizarry DMR's that overlap with ours are less
likely to be false discoveries.

Further Exploration
-------------------

The included scripts provide the basis for further exploration of
tissue-specific DMR's and tumor or disease-specific DMR's.
