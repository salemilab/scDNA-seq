
ml samtools

mkdir SIV

bam2=$(ls *mapped.bam_sorted.bam)

samtools view -b ${bam2} KU892415.1 -o SIV_sorted.bam

for barcode in $(cat *cells_full_genome.txt);

do
sbatch 02_splitSIV.sh ${barcode};

done
