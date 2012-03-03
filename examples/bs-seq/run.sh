#BSUB -J combp
#BSUB -n 12
#BSUB -e logs/combp.err
#BSUB -o logs/combp.out

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

<<GSNAP
mkdir -p data/\$(basename $READSA)
echo "python gsnap-meth.py \
    -r reference/thaliana_v10.fasta \
    -t 12 --out-dir data/$(basename $READSA) \
    --extra-args '--npaths 1 --quiet-if-excessive' \
    ${READSA}_1.fastq \
    > data/\$(basename ${READSA}).meth.txt" \
    | bsub -n 12 -J endosperm -e logs/endo.err -o logs/endo.out
exit;
mkdir -p data/\$(basename $READSB)
echo "python gsnap-meth.py \
    -r reference/thaliana_v10.fasta \
    -t 12 --out-dir data/\$(basename $READSB) \
    --extra-args '--npaths 1 --quiet-if-excessive' \
    ${READSB}_1.fastq \
    > data/\$(basename ${READSB}).meth.txt" \
    | bsub -n 12 -J embryo -e logs/embryo.err -o logs/embryo.out
GSNAP
<<FISHER

awk 'BEGIN{OFS=FS="\t"}($5 ~ /^CG/){ print $1,$2-1,$2,$3,$4,$5 }' \
    data/2_endosperm.meth.txt > data/endosperm.meth.bed

awk 'BEGIN{OFS=FS="\t"}($5 ~ /^CG/){ print $1,$2-1,$2,$3,$4,$5 }' \
    data/2_embryo.meth.txt > data/embryo.meth.bed
echo "python run-fisher.py data/embryo.meth.bed  data/endosperm.meth.bed \
    | sort -k1,1 -k2,2n > data/embryo-endo.fisher.bed" \
    | bsub -o logs/fisher.log -e logs/fisher.err

FISHER

comb-p pipeline -c 4 --dist 300 \
        --step 60 --seed 0.01 \
        -p data/cpv \
        data/embryo-endo.fisher.bed

awk '$7 < 0.05 && $5 > 3' data/cpv.regions-p.bed > data/cpv.regions-p.sig.bed
