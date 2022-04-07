
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
