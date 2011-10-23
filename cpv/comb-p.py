#!/usr/bin/env python
import sys
__actions = ("acf", "slk", "fdr", "peaks", "rpsim", "hist")

def main():
    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        print >>sys.stderr,"""
Tools for viewing and adjusting p-values in BED files.

   Contact: Brent Pedersen - bpederse@gmail.com
   License: BSD

To run, indicate one of:

   acf   - calculate autocorrelation within BED file
   slk   - Stouffer-Liptak-Kechris correction of spatially correlated p-values
   fdr   - Benjamini-Hochberg correction of p-values
   peaks - find peaks in a BED file.
   rpsim - generate p-values for a region (of p-values) by simulation.**
   hist  - plot a histogram of a column and check for uniformity.

NOTE: most of these assume *sorted* BED files.
    """
        sys.exit()
    if not sys.argv[1] in __actions:
        _pipeline()
    else:
        module = __import__(sys.argv[1])
        # remove the action
        sys.argv.pop(1)
        module.main()

def _pipeline():
    # TODO
    pass

if __name__ == "__main__":
    main()
