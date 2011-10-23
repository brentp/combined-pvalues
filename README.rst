A library to calculate adjusted p-values from spatially autocorrelated tests.
This is useful for ChIP-Seq probes and Tiling arrays.

Invocation
==========
When the program is run with::

   $ ./cpv/comb-p.py

This message is displayed::

    To run, indicate one of:

       acf   - calculate autocorrelation within BED file
       slk   - Stouffer-Liptak-Kechris correction of spatially correlated p-values
       fdr   - Benjamini-Hochberg correction of p-values
       peaks - find peaks in a BED file.
       rpsim - generate p-values for a region (of p-values) by simulation.**
       hist  - plot a histogram of a column and check for uniformity.

    NOTE: most of these assume *sorted* BED files.


Where each of the listed modules indicates an available program.
Running any of the above will result in a more detailed help message. e.g.::

    $ cpv/comb-p.py acf -h

Gives::

    usage: comb-p.py [-h] [-d D] [-c C] files [files ...]

       calculate the autocorrelation of a *sorted* bed file with a set
       of *distance* lags.

    positional arguments:
      files       files to process

    optional arguments:
      -h, --help  show this help message and exit
      -d D        start:stop:stepsize of distance. e.g. 15:500:50 means check acf
                  at distances of:[15, 65, 115, 165, 215, 265, 315, 365, 415, 465]
      -c C        column number that has the value to take the acf


Indicating that it can be run as::

    $ ./cpv/comb-p.py acf -d 1:500:50 -c 5 data/pvals.bed > data/acf.txt

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

The ACF step something like::

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

Where the first column indicates the lag-bin, the second is the
autocorrelation at that lag, and the last is the number of pairs used in
calculating the autocorrelation.
If that number is too small, the correlation values may be un-reliable.
We expect the correlation to decrease with increase lag (unless there is some
periodicity).

The first line of the output is a link to an image of the ACF data represented
in the table. It looks something like with parameter (-d 1:500:60):

.. image:: https://raw.github.com/brentp/combined-pvalues/master/data/1_500_60.png

Or, with more bins (-d 1:500:30)

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

*chr*, *start*, *end*, *pval*, *stouffer-pval*

Benjamini-Hochberg Correction
-----------------------------

This performas BH FDR correction on the pvalues. A call looks like::

    $ python cpv/fdr.py --alpha 0.05 data/pvals.acf.bed > data/pvals.adjusted.bed

where the new file has one a additional column that is the corrected p-value.
By default, it uses the last column as the p-value, but another column can
be used by specifying *-c*.

Regions
-------
We are often interested in entire regions. After running the above example, we
can find the extent of any regions using::

    $ python cpv/peaks.py --dist 500 --seed 0.1 \
                     data/pvals.adjusted.bed > data/pvals.regions.bed

where the seed inidicates a minimum value that must be see to start a region.
Again, *-c* can be used to indicate the column containing the p-values
(defaults to last column)`--dist` tells the program to merge peaks (in this case
troughs) within 150 bases of the other.
The output file is a BED file with each region and the lowest (currently)
p-value in the region.

The cpv/peaks.py script is quite flexible. Run it without arguments for
further usage.

Region Sims
-----------

Given a region we want to generate, by simulation a *p-value* for the entire
region. To do this, we do the following for each region:

 + sort the p-values
 + set S = the sum of their ranks (or the sum of the actual p-vals).
 + set RL = number of probes in the region
 + set P = number of probes in the original input file.
 + N times, do:

   - generate RL numbers between 1 and P
   - store the sum of those RL values (and the sum of the p-values at those indexes)
 + p = [report the number of times the simulated ranks or sums is < RL)/ N

This adds a column for a p-value based on the ranks and a column for a p-value
based on actual p's to the BED file. (NOT DONE).::

   $ python cpv/rpsim.py -r data/pvals.regions.bed \
                         -p data/pvals.bed \
                         -c 4 \
                         -N 10000 > data/region-ps.bed


TODO
====

1. meta script to run steps with sensible defaults.

   create executable comb-p, all steps together look like::

    comb-p -d 1:500:50 -c 5 -a 0.05 data/pvals.bed -o data/prefix

   this will run set --seed = -a and --dist == 500 (though these can also be
   specified explicitly) and will create::

    data/prefix.acf.txt # the acf correlations.
    data/prefix.acf.bed # the acf corrected bed
    data/prefix.adj.bed # the acf + fdr corrected bed
    data/prefix.regions.bed # the regions that have been run.

   each individual step and be run as::

    comb-p acf
    comb-p slk
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
   Then could have --acf as an argument. and mirror comb-p combine...
