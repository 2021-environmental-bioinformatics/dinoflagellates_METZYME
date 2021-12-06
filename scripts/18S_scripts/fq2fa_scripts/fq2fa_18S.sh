srun -p scavenger --time=04:00:00 --ntasks-per-node 2 --pty bash 
#40m
fq2fa --merge --filter WHOI021-V9_S105_L001_R1_001_val_1.fq WHOI021-V9_S105_L001_R2_001_val_2.fq WHOI021-V9_S105_L001_40m_merged.fasta
#70
fq2fa --merge --filter WHOI040-V9_S124_L001_R1_001_val_1.fq WHOI040-V9_S124_L001_R2_001_val_2.fq WHOI040-V9_S124_L001_70m_merged.fasta
#380
fq2fa --merge --filter WHOI010-V9_S94_L001_R1_001_val_1.fq WHOI010-V9_S94_L001_R2_001_val_2.fq WHOI010-V9_S94_L001_380m_merged.fasta

#Flags
#The --merge flag merges the two paired reads inputted
#The —-filter flag converts the inputted fq reads to fasta
