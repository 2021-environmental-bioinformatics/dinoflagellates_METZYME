#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=rnaspade_coassembly   # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=rnaspades_coassembly.log# Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate spades

/vortexfs1/home/selkassas/.conda/envs/megahit/bin/rnaspades.py -t 36 --s1 30B8Z_S11_001_40m_mRNA.fasta --s2 30B90_S12_001_70m_mRNA.fasta --s3 30B91_S28_001_380m_mRNA.fasta -o metaT_assembly

#Flags 
#The -t flag specifies how many threads to use
#The -s# flags are the already merged metaT reads, read in as single-end
#The -o flag specifies the name of the output assembly
