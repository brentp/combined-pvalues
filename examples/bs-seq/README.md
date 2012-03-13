BS-Seq
======

We use data from *Arabidopsis thaliana* embryo and endosperm bisulfite-treated
sequences and perform a fisher's exact test on the contingency table produced
from the number of converted and un-converted (methylated) cytosines in each
tissue.
See: https://github.com/brentp/combined-pvalues/blob/master/examples/bs-seq/run.sh
For the entire pipeline including downloading the reads and reference and
mapping the reads using Methylcoder and GSNAP.

The comb-p command used is

```Shell
comb-p pipeline -c 4 --dist 300 \
        --step 60 --seed 0.01 \
        -p data/cpv \
        data/embryo-endo.fisher.bed

```

where `embryo-endo.fisher.bed` contains the p-values from the comparison
of methylation between endopserm and embryo.
This comparison requires a p-value of 0.01 to start a region and looks up
to 300 bases to find another value of 0.01 to extend the region.

The file `data/cpv.regions-p.bed` the contains the final SLK and Sidak adjusted
p-values for each putative region.
