#!/bin/bash

# doctests
python cpv/stouffer_liptak.py
python cpv/_common.py

##################################################################
python cpv/peaks.py --dist 50 --seed 0.02 data/close_peaks.bed > t

d=$(wc -l t | awk '{ print $1}')
test $d -ne 1 && echo "ERROR" $d

##################################################################
python cpv/peaks.py --dist 1 --seed 0.02 data/close_peaks.bed > t

d=$(wc -l t | awk '{ print $1}')
test $d -ne 2 && echo "ERROR" $d

##########################################################################

for dist in 1 1000 10000000; do
    python cpv/peaks.py --dist $dist --seed 0.05 data/overlapping_peaks.bed > t
    d=$(wc -l t | awk '{ print $1}')
    test $d -ne 1 && echo "ERROR" $d
done

./cpv/comb-p.py hist -c 5 cpv/tests/data/pvals.bed > t

# unittests.
python cpv/tests/test.py

# the sum of the partials should add up the the final number in
# the non-partial
psum=$(python cpv/acf.py cpv/tests/data/pvals.bed -c 5 -d 1:500:50 \
    | awk '(NR > 2){ s += $4; }END{ print s}')
npsum=$(python cpv/acf.py cpv/tests/data/pvals.bed -c 5 -d 1:500:50 --full \
    | awk '{ s=$4 }END{ print s }')

test $psum -ne $npsum && echo "ERROR in ACF"

rm t
