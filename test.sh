#!/bin/bash

python cpv/peaks.py --dist 50 --seed 0.02 data/close_peaks.bed > t

d=$(wc -l t | awk '{ print $1}')
test $d -ne 1 && echo "ERROR" $d

python cpv/peaks.py --dist 1 --seed 0.02 data/close_peaks.bed > t

d=$(wc -l t | awk '{ print $1}')
test $d -ne 2 && echo "ERROR" $d

rm t
