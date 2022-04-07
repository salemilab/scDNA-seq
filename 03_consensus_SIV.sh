#!/bin/bash
#SBATCH --job-name=consensus   #Job name	
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=droskel@ufl.edu   # Where to send mail	
#SBATCH --time=1:00:00   # Walltime (hours-min-sec)
#SBATCH --output=consensus.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi
#SBTACH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1gb
#SBATCH --array=1-99

file1=$(ls *-1*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)

module purge

ml samtools
ml ivar

file2=$(ls /blue/salemi/share/tapestri_pipelines/SIV_reference/*fa)
   	
samtools mpileup -d 600000 -A -Q 0 -f ${file2} ${file1} | ivar consensus -p ${file1%.bam}_consensus -q 20 -t 0.5