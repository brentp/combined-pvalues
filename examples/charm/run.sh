#BSUB -J combp[4-25]
#BSUB -R rusage[mem=6]
#BSUB -e logs/combp.%I.%J.err
#BSUB -o logs/combp.%I.%J.out
#BSUB -n 12
#BSUB -R span[hosts=1]

mkdir -p logs/ 
mkdir -p data
set -e

#H=/home/brentp/src/combined-pvalues/examples/charm
#cd $H

SECTION=$1
echo $SECTION

if [ "$SECTION" = "GET" ]; then
cd data
# remove extra space and save.
wget -O - http://rafalab.jhsph.edu/data/shores/natgen2009.csv \
    | perl -pe 's/,\s/,/' > natgen2009.csv

mkdir -p xys/ 
cd xys/
wget http://rafalab.jhsph.edu/data/shores/natgen2009.tgz
tar xzvf natgen2009.tgz
exit;
fi

if [ "$SECTION" = "NORMALIZE" ]; then

R --slave --args data/natgen2009.csv data/methp.txt \
    SampleID quantile < scripts/charm.normalize.R
sed -i "1s/^/ID\t/" data/methp.txt
exit;

fi


if [ "$SECTION" = "FIT" ]; then
mkdir -p data/fit
NCOL=2162407
STEP=8000

# get rid of errant space that messes up matching.
for i in `awk -v cols=$NCOL -v step=$STEP 'BEGIN{for(i=2;i<cols;i+=step){print i }}'`; do
    start=$i
    end=$(($i + $STEP - 1))
    nf=data/fit/$start-$end.split
    bedp=data/fit/$start-$end.bed
    #cut -f 1,$start-$end data/methp.txt > $nf
    echo "~/local/bin/R --slave < scripts/fit.lm.R --args $nf > $bedp" \
        | bsub  \
        -e logs/$start-$end.err \
        -o logs/$start-$end.out \
        -R "rusage[mem=6]"
done
exit;
fi


if [ "$SECTION" = "MERGE" ]; then
cat data/fit/*.bed | awk 'NR == 1 || $1 != "#chrom"' \
    | sort -k1,1 -k2,2n > data/pvalues.bed
exit;
fi


if [ "$SECTION" = "COMB" ]; then

COLS=(fake chrom start end p.disease p.tissue)

COL=${COLS[$LSB_JOBINDEX]}

if [[ $LSB_JOBINDEX -gt 5 ]]; then
    COL="p.disease-$LSB_JOBINDEX"
fi

if [[ $LSB_JOBINDEX -gt 25 ]]; then
    COL="p.tissue-$LSB_JOBINDEX"
fi

PRE=data/quantile/$COL/$COL
mkdir -p data/quantile/$COL

comb-p pipeline \
    -c $LSB_JOBINDEX \
    -s \
    --seed 0.0005 \
    --dist 80 --step 40 \
    -p $PRE \
    data/pvalues.bed

    awk 'NR == 1 || ($5 > 5 && $7 < 0.001)' ${PRE}.regions-p.bed \
            > ${PRE}.sig.regions.bed
exit;

fi

if [ "$SECTION" = "SPLIT_FIT" ]; then

mkdir -p data/split_fit
mkdir -p logs/split_fit

for start in $(seq 2 8000 2058592); do

    f=data/split_fit/f$start.bed
    end=$((start + 8000))
    awk -vstart=$start -v end=$end 'NR == 1 || (NR > start && NR < end)' data/pvalues.bed | cut -f 1-4 > $f
    echo "comb-p pipeline -c 4 -s --seed 0.0005 --dist 80 --step 40 -p data/split_fit/s$(basename $f .bed) $f" | bsub -J $(basename $f) \
            -e logs/split_fit/$(basename $f .bed).err \
            -o logs/split_fit/$(basename $f .bed).out \
            -R "rusage[mem=4]"
    echo $f
done
exit;

fi


echo "send in a section to run, e.g.: bash run.sh FIT"
echo "sections should be run in order:
   GET
   NORMALIZE
   FIT
   MERGE
   COMB"



