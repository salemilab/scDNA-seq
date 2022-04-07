
ml samtools

bam1=$(ls *mapped.bam)

samtools sort ${bam1} -o ${bam1}_sorted.bam

bam2=$(ls *mapped.bam_sorted.bam)

samtools index ${bam2}
