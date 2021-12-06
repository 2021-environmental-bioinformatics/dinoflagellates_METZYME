#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_sortmerna                # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=elks@mit.edu                     # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=180gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_sortmerna%j.log            # Standard output/error

export OMP_NUM_THREADS=16

#tells slurm to enter your specific sortmerna conda environment

module load anaconda/5.1
source activate sortmerna

#this tells poseidon to go into the folder where all samples are located

cd /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired

#this removes the run directory each time a sample is processed. If you do not do this, sortmerna will quit with an error. 

rm -r /vortexfs1/home/selkassas/sortmerna/run

#40m
sortmerna --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/rfam-5.8s-database-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/rfam-5s-database-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-arc-16s-id95.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-arc-23s-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-bac-16s-id90.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-bac-23s-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-euk-18s-id95.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-euk-28s-id98.fasta --reads 30B8Z_S11_001.paired_trimmed.fasta --aligned 30B8Z_S11_001_40m_rRNA --fastx --other 30B8Z_S11_001_40m_mRNA -v

rm -r /vortexfs1/home/selkassas/sortmerna/run

#70m
sortmerna --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/rfam-5.8s-database-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/rfam-5s-database-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-arc-16s-id95.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-arc-23s-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-bac-16s-id90.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-bac-23s-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-euk-18s-id95.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-euk-28s-id98.fasta --reads 30B90_S12_001.paired_trimmed.fasta --aligned 30B90_S12_001_70m_rRNA --fastx --other 30B90_S12_001_70m_mRNA -v

rm -r /vortexfs1/home/selkassas/sortmerna/run

#380m
sortmerna --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/rfam-5.8s-database-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/rfam-5s-database-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-arc-16s-id95.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-arc-23s-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-bac-16s-id90.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-bac-23s-id98.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-euk-18s-id95.fasta --ref /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/rRNA_databases/silva-euk-28s-id98.fasta --reads 30B91_S28_001.paired_trimmed.fasta --aligned 30B91_S28_001_380m_rRNA --fastx --other 30B91_S28_001_380m_mRNA -v

#Flags
#The --ref flag tells sortmerna where the formatted databases are located. There can be any amount of references, as long as each is preceded by a —-ref flag. 
#The --reads flag is your input file (metaT file)
#The —-aligned specifies reads that aligned to the LSU and SSU db, i.e. the rRNA 
#The —-other flag specifies all reads that did not align to the databases, i.e. the mRNA
#The --fastx flag outputs in fasta format
#The -v flag gives verbose output in the log file

