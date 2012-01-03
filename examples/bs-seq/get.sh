<<DONE
mkdir -p reference && cd reference
rm -f thaliana_v10.fasta
for i in `seq 1 5` C M
do
    wget -O - ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr${i}.fas >> thaliana_v10.fasta
done
perl -pi -e "s/^>([^\s]+).*/>\1/;tr/C/c/" thaliana_v10.fasta
exit;
DONE

H=/home/brentp/src/methylcode/fisher-project
#<<DONE2
#rm -f $H/data/sperm2.fastq
cd $H/data/
for i in `seq 1 9`; do
    wget http://dzlab.pmb.berkeley.edu:8080/work/raw/110627_HS2B/Project_Fischer_R/Sample_TFH_45/TFH_45_NoIndex_L001_R1_00${i}.fastq.gz
    echo $i
done

for i in `seq 10 20`; do
    wget http://dzlab.pmb.berkeley.edu:8080/work/raw/110627_HS2B/Project_Fischer_R/Sample_TFH_45/TFH_45_NoIndex_L001_R1_0${i}.fastq.gz
    echo $i
done
exit;
DONE2

rm -f $H/data/veg.fastq
for i in `seq 1 9`; do
    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/raw/110627_HS2B/Project_Fischer_R/Sample_TFH_47/TFH_47_NoIndex_L003_R1_00${i}.fastq.gz | gunzip -c - >> $H/data/veg.fastq
    echo $i
done

for i in `seq 10 16`; do
    wget -O - http://dzlab.pmb.berkeley.edu:8080/work/raw/110627_HS2B/Project_Fischer_R/Sample_TFH_47/TFH_47_NoIndex_L003_R1_0${i}.fastq.gz | gunzip -c - >> $H/data/veg.fastq
    echo $i
done
