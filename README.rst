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

*chr*, *start*, *end*, *pval*, *stouffer-pval*, *reject_null_after_bh*, *bh_pval*

ACF Output
----------

It also outputs someting like::

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

Regions
-------
We are often interested in entire regions. After running the above example, we
can find the extent of any regions using::

    $ python cpv/peaks.py --dist 150 \
                          --seed 0.1 \
                          -c 7 \
                          data/pvalues.adjusted.bed > data/pvalues.regions.bed

Where the seed allows a region to start, *-c* indicates where to find the
p-values to merge, and `--dist` tells the program to merge peaks (in this case
troughs) within 150 bases of the other.
The output file is a BED file with each region and the lowest (currently)
p-value in the region.

The cpv/peaks.py script is quite flexibile. Run it without arguments for
further usage.
