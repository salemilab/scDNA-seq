#!/bin/bash
#SBATCH --job-name=split_SIV   #Job name	
#SBATCH --mail-type=NONE   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=droskel@ufl.edu   # Where to send mail	
#SBATCH --time=36:00:00   # Walltime (hours-min-sec)
#SBATCH --output=bwa_mem.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi
#SBTACH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb


ml samtools

mkdir SIV

bam2=$(ls *mapped.bam_sorted.bam)

samtools view -b ${bam2} KU892415.1 -o SIV_sorted.bam

for barcode in $(cat *cells_full_genome.txt);

do
sbatch 02_splitSIV.sh ${barcode};

done
