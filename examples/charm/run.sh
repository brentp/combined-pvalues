#BSUB -J combp[4-9]
#BSUB -R rusage[mem=6]
#BSUB -e logs/combp.%I.%J.err
#BSUB -o logs/combp.%I.%J.out
#BSUB -n 12

mkdir -p logs/
set -e

SECTION=$1
echo $SECTION

if [ "$SECTION" = "GET" ]; then
cd data
wget http://rafalab.jhsph.edu/data/shores/natgen2009.csv
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
for i in `awk -v cols=$NCOL -v step=$STEP 'BEGIN{for(i=2;i<cols;i+=step){print i }}'`; do
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
exit;
fi


if [ "$SECTION" = "MERGE" ]; then
cat data/fit/*.bed | awk 'NR == 1 || $1 != "#chrom"' \
    | sort -k1,1 -k2,2n > data/pvalues.bed
exit;
fi


if [ "$SECTION" = "COMB" ]; then

COLS=(fake chrom start end p.disease p.tissue p.colon p.frontalcortex p.liver p.spleen)
COL=${COLS[$LSB_JOBINDEX]}

PRE=data/quantile/$COL/$COL
mkdir -p data/quantile/$COL

comb-p pipeline \
    -c $LSB_JOBINDEX \
    -s \
    --seed 0.0005 \
    --dist 80 --step 40 \
    -p $PRE \
    data/pvalues.bed

    exit;
fi


if [ "$SECTION" = "FILTER" ]; then
    awk 'NR == 1 || ($5 > 5 && $7 < 0.001)' ${PRE}.regions-p.bed \
            > ${PRE}.sig.regions.bed
    exit;
fi


echo "send in a section to run, e.g.: bash run.sh FIT"
echo "sections should be run in order:
   GET
   NORMALIZE
   FIT
   MERGE
   COMB
   FILTER"



