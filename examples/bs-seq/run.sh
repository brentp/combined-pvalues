#BSUB -J gsnap
#BSUB -n 12
#BSUB -e logs/gsnap.err
#BSUB -o logs/gsnap.out

READSA=reads/WT_endosperm
READSB=reads/WT_embryo
<<DONE
rm -rf reads/*
mkdir -p reads/

#1 2 3 4 5
for i in 2 
do
    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/WT_endosperm_BS_seq_raw_batch-${i}.1.fastq >> ${READSA}_1.fastq &
    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/WT_endosperm_BS_seq_raw_batch-${i}.2.fastq >> ${READSA}_2.fastq
    wait

    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/Embryo_BS_seq_raw_batch-${i}.1.fastq >> ${READSB}_1.fastq &
    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/Embryo_BS_seq_raw_batch-${i}.2.fastq >> ${READSB}_2.fastq
    wait

done
exit;
mkdir -p reference && cd reference
rm -f thaliana_v10.fasta
for i in `seq 1 5` C M
do 
        wget -O - ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr${i}.fas >> thaliana_v10.fasta
done
perl -pi -e "s/^>([^\s]+).*/>\1/;tr/C/c/" thaliana_v10.fasta

DONE

python gsnap-meth.py \
    -r reference/thaliana_v10.fasta \
    -t 12 --out-dir data/ \
    --extra-args "--npaths 1 --quiet-if-excessive" \
        reads/WT_endosperm_1.fastq reads/WT_endosperm_2.fastq
