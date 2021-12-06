srun -p scavenger --time=04:00:00 --ntasks-per-node 2 --pty bash 
conda activate qc-trim

#files must be unzipped first
for file in *.fq.gz
do 
gunzip ${file}
done

#Proceed to fq2fa
#usage: fq2fa --merge --filter read_1.fq read_2.fq read.fa
#40m 
fq2fa --merge --filter WHOI040_16S_TCGACGAG_R1_val_1.fq WHOI040_16S_TCGACGAG_R2_val_2.fq WHOI040_16S_TCGACGAG_40m_merged.fasta
#70m 
fq2fa --merge --filter WHOI020_16S_TCGACGAG_R1_val_1.fq WHOI020_16S_TCGACGAG_R2_val_2.fq WHOI020_16S_TCGACGAG_70m_merged.fasta
#380m
fq2fa --merge --filter WHOI010_16S_AGAGTCAC_R1_val_1.fq WHOI010_16S_AGAGTCAC_R2_val_2.fq WHOI010_16S_AGAGTCAC_380m_merged.fasta

#Flags
#The --merge flag merges the two paired reads inputted
#The --filter flag converts the inputted fq reads to fasta

