#!/bin/bash
#SBATCH --job-name=consensus   #Job name	
#SBATCH --time=10:00:00   # Walltime (hours-min-sec)
#SBATCH --output=consensus.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi-b
#SBTACH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=800mb
##SBATCH --array=1-3000

date;hostname;pwd

module purge

module load samtools dibig_tools

echo "$(date): removing empty BAM files and generating consensus sequences"

# for bam in *bam; do
#     nreads=$(samtools view -c $bam)
#     if [[ $nreads -ge 10 ]]; then
#         echo $bam >> goodbams.txt
#     fi
# done &


function bam_fa {
    bam=$1
    if [[ ! -f ${bam%.bam}_consensus.fa ]]; 
    then
        nreads=$(samtools view -c $bam)
        if [[ $nreads -ge 10 ]];
        then
            submit -o --account=salemi,--qos=salemi-b /blue/salemi/share/tapestri/consensus.qsub $bam
        else
            rm $bam
        fi
    fi
}

for bam in *.bam; 
do
    n_jobs=$(qmine | wc -l)
    while [ $n_jobs -ge 1000 ]; 
    do
        sleep 0.5
        n_jobs=$(qmine | wc -l)
    done
    bam_fa $bam &
done


wait

\rm -rf empty_fastas
for fa in *fa;
do
    seq_line=$(tail -n 1 $fa)
    if echo $seq_line | grep -qvE "[ACGTN]";
    then
        bc=$(echo ${fa%_consensus.fa})
        echo $bc >> empty_fastas
    fi
done

\rm -rf qos_errors
for qsub in consensus.qsub*;
do
    if grep -q "error" $qsub;
    then
        bc=$(grep "-1" $qsub)
        echo $bc >> qos_errors
    fi
done


echo "$(date): done."
