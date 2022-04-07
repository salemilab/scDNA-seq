#!/bin/bash
#SBATCH --job-name=split_TRIM5   #Job name	
#SBATCH --mail-type=NONE   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=droskel@ufl.edu   # Where to send mail	
#SBATCH --time=36:00:00   # Walltime (hours-min-sec)
#SBATCH --output=bwa_mem.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi
#SBTACH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4gb
#SBATCH --array=1-3000


## For TRIM5a:

ml samtools

mkdir TRIM5a

bam2=$(ls *mapped.bam_sorted.bam)

samtools view -b ${bam2} NC_041767 -o TRIM5a_sorted.bam

split -l 3000 -d *barcodes_all.txt array

for ARRAY in $(ls array*);

	do
		barcode=$(cat ${ARRAY} | sed -n ${SLURM_ARRAY_TASK_ID}p);

		samtools view -r ${barcode} TRIM5a_sorted.bam -b > ./TRIM5a/${barcode}_TRIM5a.bam;

	done;