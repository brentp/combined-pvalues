# wget http://rafalab.jhsph.edu/data/shores/natgen2009.csv
# wget http://rafalab.jhsph.edu/data/shores/natgen2009.tgz

<<NORMALIZE
R --slave --args data/natgen2009.csv data/methp.txt \
    SampleID quantile < scripts/charm.normalize.R
sed -i "1s/^/ID\t/" data/methp.txt
NORMALIZE

NCOL=2162407
STEP=8000
for i in $(awk -v cols=$NCOL -v step=$STEP 'BEGIN{for(i=2;i<cols;i+=step){print i }}'); do
    start=$i
    end=$(($i + $STEP - 1))
    nf=data/fit/$start-$end.split
    bedp=data/fit/$start-$end.bed
    cut -f 1,$start-$end data/methp.txt > $nf
    echo "R --slave < scripts/fit.lm.R --args $nf > $bedp" \
        | bsub -J "$start-$end" \
        -e logs/$start-$end.err \
        -o logs/$start-$end.out \
        -R "rusage[mem=7]"
done
     
