The original Kechris paper finds enrichment regions from the Dam-Id
technology. Here, we attempt to recreate their results.

Their p-values ::

    wget http://amc-sandbox.ucdenver.edu/~brentp/2012/comb-p/data/probes.sorted.bed


We run `comb-p` on the p-values with a seed of 0.1 and a maximum distance
of 1600. ::

    comb-p pipeline \
        --seed 0.1 \
        -p cnew \
        -c 4 \
        --dist 1600 \
        probes.sorted.bed

We then compare the regions from Kechris' paper to compare to the ones founds
with `comb-p`

::

    wget http://amc-sandbox.ucdenver.edu/~brentp/2012/comb-p/data/kechris.regions.be
