#BSUB -J chip[1-4]
#BSUB -e logs/chip.%I.%J.err
#BSUB -o logs/chip.%I.%J.out

GROUP=Rep2Gm12878ControlV2 # 50
GROUP=Rep1Gm12878H3k4me3V2 # 53
GROUP=Rep2Gm12878H3k4me3V2 # 51
GROUP=Rep1Gm12878ControlV2 # 50

CGROUPS=(Rep1Gm12878ControlV2 Rep2Gm12878ControlV2 Rep1Gm12878H3k4me3V2 Rep2Gm12878H3k4me3V2 Rep1Gm12878ControlV2)
GROUP=${CGROUPS[$LSB_JOBINDEX]}
WINDOW=40

set -e
set -o pipefail

<<FIND_SHIFT
echo "shift count" > data/$GROUP.counts.txt
for i in $(seq -15 2 15); do
        awk -f ~/src/ialbert-bioawk/chipfrag.awk -v shift=$i -v pow=6 -v window=1 \
            data/$GROUP.shifted.bed \
        >> data/$GROUP.shift.counts.txt
    echo $i
done

python ../chip-seq/scripts/plot-peaks.py data/$GROUP.shift.counts.txt
sort -k2,2n data/$GROUP.counts.txt | tail -n 5
exit
FIND_SHIFT

<<CORRECT_SHIFT
awk -f scripts/correct-chip-offset.awk -v half=25 data/$GROUP.bed \
           | sort -k1,1 -k2,2n > data/$GROUP.shifted.bed

CORRECT_SHIFT

#<<MAKE_WINDOW
# overlap is w - s
bedtools makewindows -g ~/data/hg19.genome -w $WINDOW -s 10 \
        | bedtools intersect -wa -a - -b data/$GROUP.shifted.bed -c -sorted \
        | awk '($4 != 0)' > data/$GROUP.$WINDOW.counts.bed
exit

MAKE_WINDOW
