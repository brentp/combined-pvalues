<<GETDATA
mkdir -p data/sra
cd data/sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP007%2FSRP007563/SRR315573/SRR315573.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP007%2FSRP007563/SRR315574/SRR315574.sra
/vol1/home/kenjones/bin/illumina-dump --table-path SRR315573.sra --outdir .  -qseq 1
/vol1/home/kenjones/bin/illumina-dump --table-path SRR315574.sra --outdir .  -qseq 1
wget https://bitbucket.org/lullis/galaxy-tools-dev/raw/c3f1e7cb83a6/src/rgenetics/qseq2fastq/qseq2fastq.pl
# modify qseq2fastq.pl to print to stdout...

rm -f ../SRR315573.fastq; for f in SRR315573_*_qseq.txt ; do perl qseq2fastq.pl $f >> ../SRR315573.fastq; done
rm -f ../SRR315574.fastq; for f in SRR315574_*_qseq.txt ; do perl qseq2fastq.pl $f >> ../SRR315574.fastq; done
GETDATA


#!/bin/bash
#BSUB -J bwa.map
#BSUB -o logs/%J.bwa.out
#BSUB -e logs/%J.bwa.err
#BSUB -R "rusage[mem=13500]"
#BSUB -n 1
PATH=$PATH:/vol2/home/brentp/src/bwa/

set -e

FASTA=~/data/hg19.fa
THREADS=12
#FQ=data/SRR315573.fastq
FQ=data/SRR315574.fastq
GROUP=$(basename $FQ .fastq)
<<DONE
#bwa index -a bwtsw $FASTA
bwa aln -q 15 -t $THREADS $FASTA $FQ > data/$GROUP.sai

# bwa doesn't trim the quals after trimming the reads use awk to correct.
bwa samse $FASTA data/$GROUP.sai $FQ | \
    awk 'BEGIN{FS=OFS="\t"}
         ($1 ~ /^@/){ print $0} 
         ($1 !~ /^@/){ $11 = substr($11, 0, length($10)); print $0}' | \
    samtools view -bSF 4 - > data/$GROUP.unsorted.bam

samtools sort -m 3500000000 data/$GROUP.unsorted.bam data/$GROUP.q
samtools index data/$GROUP.q.bam
DONE

<<FIND_SHIFT
echo "shift	count" > data/$GROUP.counts.txt
for i in $(seq 175 1 185); do
    bamToBed -i data/$GROUP.q.bam \
        | awk -f ~/src/bioawk/chipfrag.awk -v shift=$i -v pow=6 -v window=1 \
        >> data/$GROUP.counts.txt
    echo $i
done

python scripts/plot-peaks.py data/$GROUP.counts.txt
sort -k2,2n data/$GROUP.counts.txt | tail -n 5
FIND_SHIFT

<<CORRECT
bamToBed -i data/$GROUP.q.bam \
    | awk -f scripts/correct-chip-offset.awk -v half=90 \
    | sort -k1,1 -k2,2n > data/$GROUP.shifted.bed
exit;
CORRECT

<<CHECK_NEW_PEAKS
for i in $(seq -20 3 60); do
    awk -f ~/src/bioawk/chipfrag.awk -v shift=$i -v pow=6 -v window=1 data/$GROUP.shifted.bed \
        >> data/$GROUP.corrected.counts.txt
    echo $i
done
python scripts/plot-peaks.py data/$GROUP.corrected.counts.txt
CHECK_NEW_PEAKS

WINDOW=100
<<MAKE_WINDOW
# overlap is w - s
bedtools makewindows -g ~/data/hg19.genome -w $WINDOW -s 80 \
    | bedtools intersect -wa -a - -b data/$GROUP.shifted.bed -c -sorted \
    | awk '($4 != 0)' > data/$GROUP.$WINDOW.counts.bed
exit
MAKE_WINDOW

<<UNION
bedtools unionbedg -header -names A3 A4 -i data/SRR31557{3,4}.$WINDOW.counts.bed \
    | awk 'BEGIN{OFS=FS="\t"}
      (NR == 1){ print "A3","A4" } # r expects 1 fewer cols in the header.
      (NR > 1) { print $1":"$2"-"$3,$4,$5 }' \
    > data/both.counts.$WINDOW..bed
UNION

echo "R --slave < scripts/run-deseq.R --args \
        ~/tinker/bentley-chipseq/data/both.counts.bed \
        > data/deseq.pvals.$WINDOW.txt" \
    | bsub -R "rusage[mem=25000]" -J deseq -e deseq.err -o deseq.out
