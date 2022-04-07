#!/bin/bash
#SBATCH --job-name=sort_index   #Job name	
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=droskel@ufl.edu   # Where to send mail	
#SBATCH --time=02:00:00   # Walltime (hours-min-sec)
#SBATCH --output=sort_index.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi
#SBTACH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=3gb


ml samtools

bam1=$(ls *mapped.bam)

samtools sort ${bam1} -o ${bam1}_sorted.bam

bam2=$(ls *mapped.bam_sorted.bam)

samtools index ${bam2}

