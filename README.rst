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

.. image:: https://chart.googleapis.com/chart?cht=bvs&&chd=t:0.165,0.136,0.116,0.092,0.076,0.065,0.056,0.049&chs=424x250&chco=224499&chxt=x,y&chxl=0:|1-61|61-121|121-181|181-241|241-301|301-361|361-421|421-481&chxr=1,0,0.18&chm=N,000000,0,-1,12&chbh=42,6,12&chds=a

Or, with more bins (-d 1:500:30)

.. image:: https://chart.googleapis.com/chart?cht=bvs&&chd=t:0.175,0.154,0.139,0.130,0.123,0.109,0.097,0.087,0.079,0.072,0.068,0.063,0.059,0.053,0.051,0.048&chs=824x250&chco=224499&chxt=x,y&chxl=0:|1-31|31-61|61-91|91-121|121-151|151-181|181-211|211-241|241-271|271-301|301-331|331-361|361-391|391-421|421-451|451-481&chxr=1,0,0.20&chm=N,000000,0,-1,12&chbh=42,6,12&chds=a

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
