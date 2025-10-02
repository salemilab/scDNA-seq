#!/bin/bash
#SBATCH --job-name=05_map_deletions   #Job name	
#SBATCH --mail-type=NONE   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=brittany.magalis@louisville.edu   # Where to send mail	
#SBATCH --time=03:00:00   # Walltime (hours-min-sec)
#SBATCH --output=05_map_deletions.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi-b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb

ml R/4.2

bn=$(basename "$PWD")
for i in */; do
    cd "${i}"
    if [ -f hypermut_seqs.txt ]; then
        Rscript /blue/salemi/share/tapestri/05_map_deletions.R
    fi 
    cd ../
done
