A library to calculate adjusted p-values from spatially autocorrelated tests.
This is useful for ChIP-Seq probes and Tiling arrays.


See
===

    Kechris et al. 2010:
    Generalizing Moving Averages for Tiling
    Arrays Using Combined P-Value Statistics

note
====
::

   This changes that implementation by allowing lags by *distance* (presumably)
   in bases, rather than by an index offset as is generally done with ACF.
   This makes the implementation quite a bit slower but provides more
   flexibility for probes/p-values that are not evenly spaced.

Example
=======

with data in this repository::

    python cpv/acf.py -d 15:500:50 -c 5 data/pvals.bed -a 0.1 >  data/pvalues.adjusted.bed

This takes `this BED file <https://github.com/brentp/combined-pvalues/blob/master/data/pvals.bed>`_ with the p-values in column 5, finds the autocorrelation
at lags starting at *15* and going up to *500* in steps of *50*.
Even if the start is 0, it will not include the self in the autocorrelation.

BED Output
----------
It then:

 + adjusts the p-values by stouffer-liptak with moving average.
 + performs a Benjamini-Hochberg FDR test, outputs adjusted values and
   indicates significance.and 
 + outputs a new BED file with columns:

*chr*, *start*, *end*, *pval*, *stouffer-pval*, *bh_pval*

ACF
---

The ACF step outputs someting like::

    # {link}
    lag_min-lag_max correlation N
    15-65   0.06377 4374
    65-115  0.03168 2638
    115-165 0.04158 3689
    165-215 0.0176  2679
    215-265 0.01616 3064
    265-315 0.01906 2880
    315-365 0.05305 2491
    365-415 0.04555 3097
    415-465 0.03252 2460

To STDERR, where the first column indicates the lag-bin, the second is the
autocorrelation at that lag, and the last is the number of pairs used in
calculating the autocorrelation.
If that number is too small, the correlation values may be un-reliable.
We expect the correlation to decrease with increase lag (unless there is some
periodicity).

The first line of the output is a link to an image of the ACF data represented
in the table. It looks something like with parameter (-d 1:500:60):

.. image:: http://bit.ly/qovwTN

Or, with more bins (-d 1:500:30)

.. image:: http://goo.gl/4fI5V

Regions
-------
We are often interested in entire regions. After running the above example, we
can find the extent of any regions using::

    $ python cpv/peaks.py --dist 150 \
                          --seed 0.1 \
                          -c 6 \
                          data/pvalues.adjusted.bed > data/pvalues.regions.bed

Where the seed allows a region to start, *-c* indicates where to find the
p-values to merge, and `--dist` tells the program to merge peaks (in this case
troughs) within 150 bases of the other.
The output file is a BED file with each region and the lowest (currently)
p-value in the region.

The cpv/peaks.py script is quite flexibile. Run it without arguments for
further usage.

TODO
====

1. Separate into component steps. acf.py does too much. We may just want to check
   the ACF for various step sizes before continuing on.

   The pipeline shoud look like::

    $ python acf.py -d 1:500:50 -c 5 data/pvals.bed > data/acf.txt
    $ python combine.py --acf data/acf.txt -c 5 data/pvals.bed > data/pvals.acf.bed

    # for fdr_correct and peaks.py -c could default to last_column.
    $ python fdr.py -a 0.05 -c 6 data/pvals.acf.bed > data/pvals.adjusted.bed
    $ python peaks.py --dist 500 --seed 0.05 -c 7 \
                     data/pvals.adjusted.bed > data/pvals.regions.bed

   so with an executable comb-p, all steps together look like::

    comb-p -d 1:500:50 -c 5 -a 0.05 data/pvals.bed -o data/prefix

   this will run set --seed = -a and --dist == 500 (though these can also be
   specified explicitly) and will create::

    data/prefix.acf.txt # the acf correlations.
    data/prefix.acf.bed # the acf corrected bed
    data/prefix.adj.bed # the acf + fdr corrected bed
    data/prefix.regions.bed # the regions that have been run.

   each individual step and be run as::

    comb-p acf
    comb-p combine
    comb-p fdr
    comb-p peaks

2. **Rigorous p-values for regions**.
   Since we have the stouffer-liptak for combined p-values, it should be used
   to do a correction for all p-values in a peak-region.
   This will require calculating the ACF on the input so it should be optional.
   Probably go out a given distance and then fit with a function so dont have
   to actually calculate the ACF for the full set of lags (can have very large
   regions).
   This will require keep the non-significant p-values for a region as well.
   Maybe this should be a seperate step.::

    comb-p region-correct --peaks data/prefix.regions.bed \
                          --pvals data/prefix.adj.bed \
                          -c 6 \
                          -d 1:500:50 > data/prefix.regions.pvals.ped

   Where --pvals is the file used to generated --peaks. But, if comb-p peaks
   (optionally) output all p-values in a region, we wouldn't need --pvals
   Then could have --acf and mirror comb-p combine...
