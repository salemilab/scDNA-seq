#!/bin/bash
#SBATCH --job-name=split_TRIM5   #Job name	
#SBATCH --mail-type=NONE   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=brittany.rife@ufl.edu   # Where to send mail	
#SBATCH --time=36:00:00   # Walltime (hours-min-sec)
#SBATCH --output=bwa_mem.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi-b
#SBTACH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb

## For TRIM5a:
date;hostname;pwd

ml samtools

## Generate local list of barcodes
if [[ ! -f "barcodes.txt" ]]; then
    cp ../log/*cells.txt barcodes.txt
fi


## Create a new directory for individual bam files
if [[ ! -d "TRIM5a" ]]; then
    mkdir TRIM5a
fi

# Filter and index all reads that map to the TRIM5a reference sequence
    bam2=$(ls *mapped_sorted.bam)
    samtools view -b ${bam2} NC_041767.1 -o TRIM5a_sorted.bam # For newer runs
    samtools index TRIM5a_sorted.bam


 
# Start indexing sorted bam file for individual barcodes 


cd TRIM5a

ml dibig_tools

BAM=$1
BARCODES=$2
NBC=$3

if [[ -z $BAM ]];
then
    BAM="../TRIM5a_sorted.bam"
fi
if [[ -z $BARCODES ]];
then
    BARCODES="../barcodes.txt"
fi
if [[ -z $NBC ]];
then
    NBC=1000
fi


echo $(date): splitting barcodes

split -l $NBC $BARCODES bclist

ns=$(ls bclist* | wc -l)
echo $(date): starting $ns jobs

nj=1
for bc in bclist*; do
    outdir="_split${nj}"
    mkdir -p $outdir
    submit -W -o --mem=80G generic.qsub /blue/salemi/share/bin/splitbam.py \
	   -p "${outdir}/" \
	   -r splitbam-report-${nj}.txt \
	   $bc $BAM &
    nj=$((nj+1))
done
wait

# Old version (200 mb)
# samtools view -r $1 SIV_sorted.bam -b > $2

rm bclist*

# Move all BAM files back to top-level dir
for od in _split*/; do
    mv ${od}*.bam .
    \rm -rf $od
done


date





