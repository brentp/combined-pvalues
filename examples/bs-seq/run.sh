#BSUB -J gsnap
#BSUB -n 12
#BSUB -e logs/methylcoder.err
#BSUB -o logs/methylcoder.out

READSA=reads/2_endosperm
READSB=reads/2_embryo
<<DONE

rm -f $READSA*
rm -f $READSB*
#mkdir -p reads/

for i in 1 2 3 4 5; do
    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/WT_endosperm_BS_seq_raw_batch-${i}.1.fastq >> ${READSA}_1.fastq &
    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/WT_endosperm_BS_seq_raw_batch-${i}.2.fastq >> ${READSA}_2.fastq
    wait

    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/Embryo_BS_seq_raw_batch-${i}.1.fastq >> ${READSB}_1.fastq &
    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/Embryo_BS_seq_raw_batch-${i}.2.fastq >> ${READSB}_2.fastq
    wait

done
exit;
DONE
<<REF
mkdir -p reference && cd reference
rm -f thaliana_v10.fasta
for i in `seq 1 5` C M
do 
        wget -O - ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr${i}.fas >> thaliana_v10.fasta
done
perl -pi -e "s/^>([^\s]+).*/>\1/;tr/C/c/" thaliana_v10.fasta

REF

mkdir -p data/$(basename $READSA)
python gsnap-meth.py \
    -r reference/thaliana_v10.fasta \
    -t 12 --out-dir data/$(basename $READSA) \
    --extra-args "--npaths 1 --quiet-if-excessive" \
    ${READSA}_1.fastq ${READSA}_2.fastq
    > data/$(basename ${READSA}).meth.txt

mkdir -p data/$(basename $READSB)
python gsnap-meth.py \
    -r reference/thaliana_v10.fasta \
    -t 12 --out-dir data/$(basename $READSB) \
    --extra-args "--npaths 1 --quiet-if-excessive" \
    ${READSB}_1.fastq ${READSB}_2.fastq
    > data/$(basename ${READSB}).meth.txt
