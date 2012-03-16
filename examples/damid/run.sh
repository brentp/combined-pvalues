#BSUB -J combp
#BSUB -n 12
#BSUB -e comb.err
#BSUB -o comb.out

# get the probe p-values from Kechris data
# wget http://amc-sandbox.ucdenver.edu/~brentp/2012/comb-p/data/probes.sorted.bed
comb-p pipeline \
    --seed 0.1 \
    -p cnew \
    -c 4 \
    --dist 1600 \
    probes.sorted.bed

# the regions from Kechris' paper to compare to
wget http://amc-sandbox.ucdenver.edu/~brentp/2012/comb-p/data/kechris.regions.bed
