A library to combine, analyze, group and correct p-values in BED files.
Unique tools involve correction for spatial autocorrelation tests.
This is useful for ChIP-Seq probes and Tiling arrays, or any data with spatial
correlation.

About
=====

The Bioinformatics Applications Note manuscript is available here:
    http://bioinformatics.oxfordjournals.org/content/28/22/2986.full

free link:
    http://bioinformatics.oxfordjournals.org/cgi/reprint/bts545?ijkey=ZTTOnczUJYLfKgw&keytype=ref

It includes an explanation of 3 examples in the examples directory
of this repository.

The software is distributed under the MIT license.

QuickStart
==========

If your data is a sorted BED (first columns are chrom, start, stop) with a column for
p-value in the 4th column, you can find DMRs as::

    comb-p pipeline \
        -c 4 \
        --mlog \
        --step 100 \
        --dist 200
        --seed 1e-3 \
        -p $OUT_PREFIX \
        --region-filter-p 0.1 \
        --anno mm9 \
        $PVALS

Where the `seed` indicates the maximal p-value that can seed a region and distance is
how far to extend a region without seeing another p-value that low. In the end, regions
with a corrected (Sidak) p-value of 0.1 will be annotated with genome version mm9 (any
version from UCSC database is supported).
The output will look like:

    https://github.com/brentp/combined-pvalues/blob/master/manuscript/anno.tsv

Commands below give finer control over each step.

Installation
============

run::

    sudo python setup.py install

to have `comb-p` installed on your path.
Otherwise, you can use the python scripts in the cpv subdirectory.
E.g.
::

    python cpv/peaks.py

corresponds to the command::

    comb-p peaks


Invocation
==========
The program is run with::

   $ comb-p

This message is displayed::

    To run, indicate one of:

       acf   - calculate autocorrelation within BED file
       slk   - Stouffer-Liptak-Kechris correction of spatially correlated p-values
       fdr   - Benjamini-Hochberg correction of p-values
       peaks - find peaks in a BED file.
       region_p  - generate SLK p-values for a region (of p-values)
       hist      - plot a histogram of a column and check for uniformity.
       splot     - a scatter plot of column(s) in a bed file for a given region.
       manhattan - a manhattan plot of values in a BED file.


    NOTE: most of these assume *sorted* BED files.


Where each of the listed modules indicates an available program.
Running any of the above will result in a more detailed help message. e.g.::

    $ comb-p acf -h

Gives::

    usage: comb-p [-h] [-d D] [-c C] files [files ...]

       calculate the autocorrelation of a *sorted* bed file with a set
       of *distance* lags.

    positional arguments:
      files       files to process

    optional arguments:
      -h, --help  show this help message and exit
      -d D        start:stop:stepsize of distance. e.g. 15:500:50 means check acf
                  at distances of:[15, 65, 115, 165, 215, 265, 315, 365, 415, 465]
      -c C        column number with p-values for acf calculations


Indicating that it can be run as::

    $ .comb-p acf -d 1:500:50 -c 5 data/pvals.bed > data/acf.txt

Each module is described in detail below.

Example
=======

Find and merge peaks/troughs within a bed file
----------------------------------------------
::

     python cpv/peaks.py --seed 0.05 --dist 1000 data/pvals.bed > data/pvals.peaks.bed

This will seed peaks with values < 0.05 and merge any adjacent values
within 1KB. The output is a BED file containing the extent of the troughs.
If the argument `--invert` is specified, the program will find look for
values larger than the seed.

Pipeline
========

The default steps are to:

 1) calculate the ACF
 2) use the ACF to do the Stouffer-Liptak correction
 3) do the Benjamini-Hochberg FDR correction
 4) find regions from the adjusted p-values.

Inputs and outputs to each step are BED files.

Note that any of these steps can be run independently, e.g. to do multiple
testing correction on a BED file with p-values, just call the fdr.py script.

ACF
---
To calclulate autocorrelation from 1 to 500 bases with a stepsize of 50
on the p-values in column 5, the command would look something like:

    $ python cpv/acf.py -d 1:500:50 -c 5 data/pvals.bed > data/acf.txt

The ACF will look something like::

    # {link}
    lag_min lag_max correlation N
    1   51  0.06853 2982
    51  101 0.04583 4182
    101 151 0.02719 2623
    151 201 0.0365  3564
    201 251 0.0005302   2676
    251 301 0.02595 3066
    301 351 0.04935 2773
    351 401 0.04592 2505
    401 451 0.03923 2972

Where the first and second columns indicate the lag-bin, the third is the
autocorrelation at that lag, and the last is the number of pairs used in
calculating the autocorrelation.
If that number is too small, the correlation values may be unreliable.
We expect the correlation to decrease with increase lag (unless there is some
periodicity).

The first line of the output is a link to an image of the ACF data represented
in the table. For parameter -d 1:500:60 it looks like:

.. image:: https://raw.github.com/brentp/combined-pvalues/master/data/1_500_60.png

Or, with more bins -d 1:500:30:

.. image:: https://raw.github.com/brentp/combined-pvalues/master/data/1_500_30.png

That output should be directed to a file for use in later steps.

Combine P-values with Stouffer-Liptak-Kechris correction
--------------------------------------------------------

See
+++

    Kechris et al. 2010:
    Generalizing Moving Averages for Tiling
    Arrays Using Combined P-Value Statistics

    This changes that implementation by allowing lags by *distance* (presumably)
    in bases, rather than by an index offset as is generally done with ACF.
    This makes the implementation quite a bit slower but provides more
    flexibility for probes/p-values that are not evenly spaced.

Usage
+++++

The ACF output is then used to do the Stouffer-Liptak-Kechris correction.
A call like::

    $ python cpv/slk.py --acf data/acf.txt -c 5 data/pvals.bed > data/pvals.acf.bed

 + adjusts the p-values by stouffer-liptak with values from the autocorrelation
   in the step above.
 + outputs a new BED file with columns:

*chr*, *start*, *end*, *pval*, *combined-pval*

Benjamini-Hochberg Correction
-----------------------------

This performs BH FDR correction on the pvalues. A call looks like::

    $ python cpv/fdr.py data/pvals.acf.bed > data/pvals.adjusted.bed

where the new file has one additional column, the corrected p-value. By
default, it uses the last column as the p-value input, but another column can
be used by specifying *-c*.

Regions
-------
We are often interested in entire regions. After running the above example, we
can find the extent of any regions using::

    $ python cpv/peaks.py --dist 500 --seed 0.1 \
                     data/pvals.adjusted.bed > data/pvals.regions.bed

where the seed inidicates a minimum p-value to start a region.
Again, *-c* can be used to indicate the column containing the p-values
(defaults to last column)`--dist` tells the program to merge peaks (in this case
troughs) within 500 bases of the other.
The output file is a BED file with each region and the lowest (currently)
p-value in the region.

The cpv/peaks.py script is quite flexible. Run it without arguments for
further usage.

ScatterPlot (splot)
-------------------

The command::

    comb-p splot -c 5,6 data/pvals.adjusted.bed \
                                -r chrY:2717613-2728613 \
                                --labels original,adjusted

will plot columns 5 and 6 from the region `-r`, resulting in

.. image:: https://raw.github.com/brentp/combined-pvalues/master/data/scatter.png

larger regions will automatically be plotted as points.
You may specify any number of columns to plot.


Region P-values (region_p)
--------------------------

Currently, the reported p-value is a Stouffer-Liptak *p-value* for the entire
region. This is done by taking a file of regions, and the original,
uncorrected p-values, calculating the ACF out to the length of the longest
region, and then using that ACF to perform the Stouffer-Liptak correction on
each region based on the original p-values.
The 1-step Sidak correction for multiple testing is performed on the p-value
for the region. Because the original p-values are sent in, we know the
coverage of the input. The Sidak correction is then based on the number of
possible regions of the current size that could be created from the total
coverage. The extra columns added to the output file are the Stouffer-Liptak
p-value of the region and the Sidak correction of that p-value.


An invocation::

   $ comb-p region_p -p data/pvals.bed \
                     -r data/regions.bed \
                     -s 50 \
                     -c 5 > data/regions.sig.bed

Will extract p-values from column 5 of pvals.bed for lines within regions in
regions.bed. It will set tau to (-t) 0.1, use a step-size of 50 for the ACF
calculation.

Frequently Asked Questions
==========================

See the Wiki `F.A.Q.`_

.. _`F.A.Q.`: https://github.com/brentp/combined-pvalues/wiki/F.A.Q.
