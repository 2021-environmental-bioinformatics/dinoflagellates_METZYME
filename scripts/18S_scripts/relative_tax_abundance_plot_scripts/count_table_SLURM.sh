#!/bin/bash
#SBATCH --partition=compute                  # Queue selection
#SBATCH --job-name=18S_count_table              # Job name
#SBATCH --mail-type=ALL                      # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu       # Where to send mail
#SBATCH --ntasks=1                           # Run a single task
#SBATCH --cpus-per-task=36                   # Number of CPU cores per task
#SBATCH --mem=100gb                          # Job memory request
#SBATCH --time=24:00:00                      # Time limit hrs:min:sec
#SBATCH --output=18S_count_table.log            # Standard output/error
#export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate R_environment

Rscript 18S_rRNA-tax-loop.R
