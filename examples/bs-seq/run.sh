#PBS -o logs/map-sperm.log.out
#PBS -e logs/map-sperm.log.err
#PBS -N sperm

set -e
HERE=/home/brentp/src/methylcode/fisher-project
H=$HERE
<<DONE
python $HERE/../gsnap-meth.py -r $HERE/../reference/thaliana_v10.fasta --out-dir $HERE/data/sperm --extra-args "--quality-protocol sanger" $HERE/data/sperm-new.fastq > $HERE/data/sperm/sperm.per-base.txt
DONE

<<DONE
for grp in sperm veg; do
    for ctx in CG CHG CHH; do
        awk -v ctx=$ctx 'BEGIN{OFS=FS="\t"}($5 ~ ctx){ print $1,$2-1,$2,$3,$4,$5 }' data/$grp/$grp.per-base.txt  > data/$grp/$grp.$ctx.bed
    done
done

for ctx in CG CHG CHH; do
    echo "python $H/run-fisher.py $H/data/veg/veg.$ctx.bed $H/data/sperm/sperm.$ctx.bed > $H/data/p.$ctx.fisher.bed" | qsub
done


for ctx in CHG CHH; do
    python ~/src/combined-pvalues/cpv/comb-p.py pipeline -c 4 --dist 300 \
        --step 60 --seed 0.01  data/p.$ctx.fisher.bed -p data/cpv.$ctx
    <<DONE
    awk 'BEGIN{OFS="\t"}($7 < 0.01 && $5 > 6 || NR == 1){ print $1,$2,$3,$5,$7 }' data/cpv.CG.regions-p.bed > data/cpv.CG.regions-sig.bed
    DONE
    awk '$7 < 0.05 && $5 > 6 || NR == 1'  data/cpv.$ctx.regions-p.bed \
                | cut -f 1-3,5,7 > data/$ctx.regions.bed

    python intersect.py data/$ctx.regions.bed > data/$ctx.genes.bed
done

echo	"#chrom	start	end	n-probes	slk-sidak-p	nearest	gene	cs	ts" > $H/data/cpv.CG.sperm.annotated.bed
echo	"#chrom	start	end	n-probes	slk-sidak-p	nearest	gene	cs	ts" > $H/data/cpv.CG.veg.annotated.bed
for chrom in `seq 1 5`; do
  grep "^chr${chrom}" $H/data/sperm/sperm.CG.bed | bedtools intersect -header -a $H/data/cpv.CG.annotated.bed -b - -wo > tmp
  bedtools groupby -g 1,2,3,4,5,6,7 -c 11,12 -o sum,sum -i tmp >> $H/data/cpv.CG.sperm.annotated.bed
  grep "^chr${chrom}" $H/data/veg/veg.CG.bed | bedtools intersect -header -a $H/data/cpv.CG.annotated.bed -b - -wo > tmp
  bedtools groupby -g 1,2,3,4,5,6,7 -c 11,12 -o sum,sum -i tmp >> $H/data/cpv.CG.veg.annotated.bed
done;

DONE

paste $H/data/cpv.CG.sperm.annotated.bed $H/data/cpv.CG.veg.annotated.bed \
    | awk 'BEGIN{OFS=FS="\t"}
(NR == 1) { print $1,$2,$3,$4,$5,$6,$7,"sperm-meth","veg-meth" }
(NR > 1) { print $1,$2,$3,$4,$5,$6,$7,$8/($8 + $9),$17/($17 + $18) }' > $H/data/methp.CG.annotated.bed
