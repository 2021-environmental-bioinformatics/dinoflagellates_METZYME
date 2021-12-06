#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=mothur_16S               # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=2                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=mothur_16S%j.log           # Standard output/error

export OMP_NUM_THREADS=12

module load anaconda/5.1
source activate mothur_1.46.1

mothur mothur_16S_plasmids.txt
#This activates mothur then runs the classify.seqs mothur_16S_plasmids.txt script written
