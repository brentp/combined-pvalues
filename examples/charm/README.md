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
with Sidak corrected region p-values less than 0.1 and with more than 8 probes in the DMR.

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

We can use the manhattan plot to get a view of the p-values with the command:

```Shell
comb-p manhattan data/pvalues.bed -c 5 --image data/manhattan.tissue.png
```
to create this plot:
![Manhattan Plot of P-values](https://github.com/brentp/combined-pvalues/raw/master/examples/charm/data/manhattan.tissue.png "Manhattan Plot") 

We can see that there are *many* significant single probes.
We find tissue-specific DMRs and filter to get 11,061 to compare to
Irizarry's 16,379.
This time we filter to a corrected region p-value of 0.05 and require at least
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
      | cut -f 4 > no-overlap.txt                    
$ intersectBed -a data/t-dmrs.txt -b tissue.filtered.bed  \
      | cut -f 4 > overlap.txt

$ python ../damid/ttest.py overlap.txt no-overlap.txt 
(25.948028554040935, 1.2957717193202179e-145)
```
So the p-value of 1.296e-145 indicates that we can enrich for
lower q-values from the Irizarry data by intersecting with our
DMR's. So the Irizarry DMR's that overlap with ours are less
likely to be false discoveries.


Shuffling
---------

We can shuffle the clinical data (DiseaseStatus and TissueType relative
to the lab data and check what the pipeline will call.
In fact, when running comb-p on shuffled data, it most often finds *no*
DMR's. On 20 repeats of randomly shuffled data, we find a total of only
37 false-positive DMR's; 34 of those come from 2 of the random sets,
indicating that they happened to coincide fairly well with the case-control
status of the original data.
This indicates a false discovery rate of essentially 0.0016.
We did this with a command like

```Shell
$ comb-p pipeline -c 6 -s --seed 0.0005 --dist 80 --step 40 -p $PREFIX $SHUFFLED_PVALS
```
The shuffled p-values were generated with the script:
https://github.com/brentp/combined-pvalues/blob/master/examples/charm/scripts/fit.lm.R

We then count only the regions that are significant after the SLK-Sidak
correction.

Power
-----
Conversely, we can show that we find more total probes that are differentially
methylated when using our SLK method than using a simple FDR correction.

```Shell
$ python ../../cpv/fdr.py -c 4 /tmp/pvalues.bed | awk '$8 < 0.0005' | wc -l
1238
```
shows that there are 1238 probes with a p-value less than our chosen cutoff.
However, if we count the total number of probes in regions with a corrected
p-value less then 0.05 (note the distinction between probe and region
p-values):

```Shell
$ awk 'BEGIN{sum=0}(NR > 1 && $7 < 0.05){ sum += $5; }END{ print sum }' /tmp/p.disease.regions-p.bed 
31514
```
We have a much larger number *31,514* vs *1,238*. (The number is higher if we
don't limit to significant regions).

Of course, this is necessary in studies where we have smaller differences and the original
p-values at each probe are not low enough to survive multiple-testing correction.
In fact, in most cases, we find 0 probes with a Benjamini-Hochberg q-value <
0.05.

Further Exploration
-------------------

The included scripts provide the basis for further exploration of
tissue-specific DMR's and tumor or disease-specific DMR's.
