
file1=$(ls *-1*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)

module purge

ml samtools
ml ivar

file2=$(ls /blue/salemi/share/tapestri_pipelines/SIV_reference/*fa)
   	
samtools mpileup -d 600000 -A -Q 0 -f ${file2} ${file1} | ivar consensus -p ${file1%.bam}_consensus -q 20 -t 0.5
