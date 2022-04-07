#!/bin/bash
#SBATCH --job-name=split_SIV   #Job name	
#SBATCH --mail-type=NONE   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=droskel@ufl.edu   # Where to send mail	
#SBATCH --time=36:00:00   # Walltime (hours-min-sec)
#SBATCH --output=bwa_mem.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi-b
#SBTACH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb


ml samtools
samtools view -r ${1} SIV_sorted.bam -b > ./SIV/${1}_SIV.bam
