#!/bin/bash
#SBATCH --job-name=rename   #Job name	
#SBATCH --mail-type=FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=brittany.magalis@louisville.edu   # Where to send mail	
#SBATCH --time=04:00:00   # Walltime (hours-min-sec)
#SBATCH --output=04_rename.%j.out   # Name output file 
#SBATCH --account=salemi
#SBATCH --qos=salemi-b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=el9
#SBATCH --mem=2gb


ml R/4.2

Rscript /blue/salemi/share/tapestri/04_seq_rename.R -m $1 -d $2 -t $3

date
