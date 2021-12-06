#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=hmmsearch_pf3             # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=hmmsearch_pf3%j.log         # Standard output/error
#export OMP_NUM_THREADS=16

module load anaconda/5.1
source activate hmmer

hmmbuild PF01036_uniprot_proteorhodopsin.hmm PF01036_uniprot_proteorhodopsin.txt 
hmmbuild PF03713_iron_storage_protein_uniprot.hmm  PF03713_iron_storage_protein_uniprot.txt
hmmbuild PF07692_low_iron_inducible_periplasmic_protein_uniprot.hmm  PF07692_low_iron_inducible_periplasmic_protein_uniprot.txt

#run with coassembly ONLY

#PF01036_uniprot_proteorhodopsin
hmmsearch -o dino_metzyme_annotated_coassembly_proteorhodopsin_out PF01036_uniprot_proteorhodopsin.hmm dino_metzyme_annotated_coassembly.faa

#PF03713_iron_storage_protein
hmmsearch -o dino_metzyme_annotated_coassembly_iron_storage_protein_out PF03713_iron_storage_protein_uniprot.hmm dino_metzyme_annotated_coassembly.faa

#PF07692_low_iron_inducible_periplasmic_protein
hmmsearch -o dino_metzyme_annotated_coassembly_low_iron_inducible_periplasmic_protein_out PF07692_low_iron_inducible_periplasmic_protein_uniprot.hmm dino_metzyme_annotated_coassembly.faa

#Flags
#The -o flag is the name of the output. You then input the built hmm file and lastly enter the .faa file to be searched. 
