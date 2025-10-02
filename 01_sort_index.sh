#!/bin/bash
#SBATCH --job-name=sort_index   #Job name	
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=brittany.magalis@louisville.edu   # Where to send mail	
#SBATCH --time=04:00:00   # Walltime (hours-min-sec)
#SBATCH --output=bwa_mem.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi-b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb


ml samtools

bam1=$(ls *mapped.bam)

samtools sort ${bam1} -o ${bam1%.bam}_sorted.bam

bam2=$(ls *mapped_sorted.bam)

samtools index ${bam2}

date
