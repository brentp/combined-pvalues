#!/usr/bin/env python
import sys
__actions = ["acf", "combine", "fdr", "peaks"]
def main():
    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        print >>sys.stderr,"""
Tools for adjusting p-values in BED files.
        Contact: Brent Pedersen - bpederse@gmail.com
        License: BSD


To run, indicate one of:

    acf     - calculate autocorrelation within BED file
    combine - Stouffer-Liptak correction of p-values
    fdr     - Benjamini-Hochberg correction of p-values
    peaks   - find peaks in a BED file.

NOTE: all of these assume *sorted* BED files.
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
