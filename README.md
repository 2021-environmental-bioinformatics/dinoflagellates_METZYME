# dinoflagellates_METZYME

We chose the paper “Dinoflagellates alter their carbon and nutrient metabolic strategies across environmental gradients in the central Pacific Ocean” by Cohen et al. We propose to reanalyze the meta- genomic and transcriptomic data to regenerate the bar plots of relative community abundance shown in Figures 1b and 1c, “Protistan community composition across the METZYME transect”. We will also re-annotate the metatranscriptomic data to recreate the TPM-normalized gene expression heat map shown in Figure 3a, “Distinct dinoflagellate functional metabolism between the euphotic and mesopelagic zones of the central Pacific". The metatranscriptomic and 16S and 18S rRNA metabarcoding data generated for this paper was from 42 seawater samples collected for biomass at 7 sites and between 3 and 13 depths per site, with the majority between 40 and 600m. To simplify the analysis, we will focus on Station 9, which has samples from three different depths (40m, 70m, and 380m). This dataset is multi-pronged and stored in many different repositories: the Proteome Xchange Consortium, Ocean Protein Portal, NSF’s Biological and Chemical Oceanography Data Management Office, and NCBI, though we had to get the data from one of the author. The dataset is approximately 15Gb. 

Cohen, Natalie R., Matthew R. McIlvin, Dawn M. Moran, Noelle A. Held, Jaclyn K. Saunders, Nicholas J. Hawco, Michael Brosnahan et al. "Dinoflagellates alter their carbon and nutrient metabolic strategies across environmental gradients in the central Pacific Ocean." Nature Microbiology 6, no. 2 (2021): 173-186.

Data availability
——————————————————————————————————————————————
The data is available at the following links:

Metatranscriptomic reads: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA555787
16S reads: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA555787
located under Biosample accession numbers: SAMN12331629–SAMN12331670 
18S reads:https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA555787
located under Biosample accession numbers:  SAMN12332710–SAMN12332751

Metatranscriptome processing

Merging forward and reverse MetaT reads then Removing rRNA reads from MetaT with SortMeRNA 4.2.0
————————————————————————————————————————————————————————————————————————————————————————————————
Programs: SortMeRNA 4.2.0

Metatranscriptome QC - did not QC this data; all are fasta (.fna) files

#Create the conda environment
```
conda create --name sortmerna
conda activate sortmerna
conda install -c bioconda sortmerna
```


SortMeRNA 4.2.0
```
1. Create databases used for rRNA screening. Based on riboPicker’s program defaults (what they used), I will use SILVA Small subunit reference database (16S/18S) & SILVA Large subunit reference database (23S/28S). No version was given, so I will use SILVA 111 (also used for the 16S and 18S - This WAS specified in the paper.)

Source to cite:
Pruesse, E., C. Quast, K. Knittel, B. Fuchs, W. Ludwig, J. Peplies, and F.
O. Glockner. SILVA: a comprehensive online resource for quality checked and
aligned ribosomal RNA sequence data compatible with ARB. 
Nuc. Acids Res. 2007; Vol. 35, No. 21, p. 7188-7196
```

#This will give you the newest version. It is 4.2.0. Any 'indexdb_rna' commands will not work. There is only one command now, 'sortmerna'. Everything, including indexing your databases will be done with this. 

2. Clone the GitHub repo and its contents:

git clone https://github.com/biocore/sortmerna.git

3. Now there will be a directory called "sortmerna" in conda path. 
cd ~/.conda/envs/sortmerna
cd sortmerna
#the sortmerna directory is a subdirectory of the main sortmerna path created by conda

4. Reformat SSU and LSU databases using bash code below
Obtain SSU and LSU by going to the SILVA archives; we used v.132 —> https://www.arb-silva.de/download/arb-files/
```
#!/bin/sh
#
#   script: (1-2) This script takes as input an rRNA database in the format 
#           of a fasta file and an xml file, both exported from the 
#           ARB package for some complete rRNA database (ex. SILVA).
#           We use the xml description tags to filter out rRNA
#           which are not 16S or 23S rRNA. 
#           ARB source http://www.arb-home.de/
#
#			(3) The resulting (only 16S or 23S)
#           fasta file is further filtered by PRINSEQ to remove
#           long sequences and sequences with more than 1% occurrences
#           of the ambiguous letter 'N'. Long sequences are defined to be
#           >1600 nucleotides for 16S rRNA and >5000 nucleotides for 
#           23S rRNA. This aids to remove sequences labeled as `complete 
#           genome', which can include a complete rRNA molecule along
#           with other coding or non-coding nucleic acids appended to 
#           either end of the molecule (mRNA, tRNA ..). Since only rRNA
#           sequences are masked out from the NCBI complete bacterial 
#           genomes, the resulting non-rRNA reads can match to an rRNA
#           database if sequences labeled as `complete genome' are not 
#           filtered out. 
#           PRINSEQ source http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.17.4.tar.gz
#
#			(4) Finally, we apply UCLUST to the filtered database to create
#           a representative rRNA database with x% id.
#	        UCLUST source http://www.drive5.com/usearch/
#
#   contact: evguenia.kopylova@lifl.fr
#
#   date: July 2012
#


# path to ARB-exported rRNA sequences xml file
SILVAxml=""
# path to ARB-exported rRNA sequences fasta file
SILVAfasta=""
# path to PRINSEQ
PRINSEQ="prinseq-lite.pl"
# path to UCLUST
UCLUST="usearch5.2.32_i86linux32"
# UCLUST -id value
ID=98

# Step (1) find non-23s sequences

cat $SILVAxml | grep "5S" -B 45 | sed -n 's/<species name="\(.*\)">/\1/p' > rrna5s.txt
cat $SILVAxml | grep "5.8S" -B 45 | sed -n 's/<species name="\(.*\)">/\1/p' > rrna5.8s.txt
cat $SILVAxml | grep "16S" -B 45 | sed -n 's/<species name="\(.*\)">/\1/p' > rrna16s.txt
cat $SILVAxml | grep "18S" -B 45 | sed -n 's/<species name="\(.*\)">/\1/p' > rrna18s.txt
cat $SILVAxml | grep "26S" -B 45 | sed -n 's/<species name="\(.*\)">/\1/p' > rrna26s.txt
cat $SILVAxml | grep "28S" -B 45 | sed -n 's/<species name="\(.*\)">/\1/p' > rrna28s.txt

cat rrna*.txt > totalnon23srrna_temp.txt
sort totalnon23srrna_temp.txt | uniq -c | awk '{print $2;} ' > totalnon23srrna.txt
rm rrna*.txt
rm totalnon23srrna_temp.txt

echo "Found all non-23S rRNA."
echo `wc -l totalnon23srrna.txt`

# Step (2) remove non-23s sequences

# put fasta file into 1 line per entry
awk '/^>/ {printf("\n%s\t",$0);next;} {printf("%s",$0);} END {printf("\n");}' < $SILVAfasta \
	| egrep -v '^$' > ${SILVAfasta%.fasta}-1line.fasta

while read not23s
do
	sed -i "/$not23s /d" ${SILVAfasta%.fasta}-1line.fasta
done < "totalnon23srrna.txt"
rm totalnon23srrna.txt

# put back into 2 lines per entry
tr "\t" "\n" < ${SILVAfasta%.fasta}-1line.fasta > ${SILVAfasta%.fasta}-no23s.fasta

echo "Removed all non-23S rRNA."

# Step (3) PRINSEQ (maximum length 5000, maximum ambigious N's 1%)
$PRINSEQ -fasta ${SILVAfasta%.fasta}-no23s.fasta \
   -max_len 5000 -ns_max_p 1 -verbose -out_good ${SILVAfasta%.fasta}-no23s-prinseq -out_bad null

# Step (4) UCLUST
$UCLUST -sort ${SILVAfasta%.fasta}-no23s-prinseq.fasta -output ${SILVAfasta%.fasta}-no23s-prinseq-sorted.fasta
$UCLUST -cluster ${SILVAfasta%.fasta}-no23s-prinseq-sorted.fasta -seedsout ${SILVAfasta%.fasta}-database-id$ID.fasta \
	-id 0.$ID

echo "Results:"

echo `grep -c '>' ${SILVAfasta%.fasta}-database-id$ID.fasta`

echo "Done."
```

#SortmeRNA code - merged, paired reads 
Script name in GitHub: sortmerna_metaT.sh
```
	#general format/usage: 

sortmerna --ref <insert path to reference database 1> --ref <insert path to reference database 2> --ref <insert path to reference database 3> --reads <insert merged reads file from Flash (will end in .extendedFrags.fastq> --aligned <insert file name_rRNA> --fastx --other <insert file name_mRNA> --otu_map -v

#Script

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_sortmerna                # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=email                            # Where to send mail
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
```

Source: Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611.

Assembling MetaT reads using SPAdes with RNASPAdes feature 
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Programs used: SPAdes v.3.15.3

Coassembled using RNASPAdes
Script name in GitHub: RNASPAdes_coassembly.sh
```
#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=rnaspade_coassembly   # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=email               # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=rnaspades_coassembly.log# Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate spades
 
/vortexfs1/home/selkassas/.conda/envs/megahit/bin/rnaspades.py -t 36 —s1 30B8Z_S11_001_40m_mRNA.fasta —s2 
30B90_S12_001_70m_mRNA.fasta —s3 30B91_S28_001_380m_mRNA.fasta -o metaT_assembly

#Flags 
#The -t flag specifies how many threads to use
#The -s# flags are the already merged metaT reads, read in as single-end
#The -o flag specifies the name of the output assembly
```
 

Annotating metaT assembly using FragGeneScan v.1.31
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Programs used: FragGeneScan v.1.31

Description: FragGeneScan is an application for finding (fragmented) genes in short reads. It can also be applied to predict prokaryotic genes in incomplete assemblies or complete genomes.

documentation on how to run: https://github.com/gaberoo/FragGeneScan

They use version 1.16 in the paper, but it was not available via conda, so I just used 1.31. 

Setting up the conda environment: 
```
conda create --name FragGeneScan
conda activate FragGeneScan
conda install -c bioconda fraggenescan
```

Code to run FragGeneScan:
Script name in GitHub: fraggenescan_script_coassembly_and_mRNA.sh
```
USAGE: ./FragGeneScan -s [seq_file_name] -o [output_file_name] -w [1 or 0] -t [train_file_name] (-p [thread_num])

       Mandatory parameters
       [seq_file_name]:    sequence file name including the full path
       [output_file_name]: output file name including the full path
       [1 or 0]:           1 if the sequence file has complete genomic sequences
                           0 if the sequence file has short sequence reads
       [train_file_name]:  file name that contains model parameters; this file should be in the "train" directory
                           Note that four files containing model parameters already exist in the "train" directory
                           [complete] for complete genomic sequences or short sequence reads without sequencing error
                           [sanger_5] for Sanger sequencing reads with about 0.5% error rate
                           [sanger_10] for Sanger sequencing reads with about 1% error rate
                           [454_5] for 454 pyrosequencing reads with about 0.5% error rate
                           [454_10] for 454 pyrosequencing reads with about 1% error rate
                           [454_30] for 454 pyrosequencing reads with about 3% error rate
                           [illumina_5] for Illumina sequencing reads with about 0.5% error rate
                           [illumina_10] for Illumina sequencing reads with about 1% error rate

       Optional parameter
       [thread_num]:       the number of threads used by FragGeneScan; default is 1 thread.

#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=FragGeneScan_annot    # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=12               # Number of CPU cores per task
#SBATCH --mem=80gb                       # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=FragGeneScan_annot.log  # Standard output/error
#export OMP_NUM_THREADS=12

module load anaconda/5.1
source activate FragGeneScan

#Coassembly
FragGeneScan -s /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/metaT_assembly/transcripts.fasta -o /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/dino_metzyme_annotated_coassembly -w 0 -t illumina_5 -p 8

#mRNA file for each depth. This will be necessary for the diamond step later on in the pipeline. 

FragGeneScan -s /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B8Z_S11_001_40m_mRNA.fasta -o /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B8Z_S11_001_40m_mRNA_annotated -w 0 -t illumina_5 -p 8

FragGeneScan -s /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B90_S12_001_70m_mRNA.fasta -o /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B90_S12_001_70m_mRNA_annotated -w 0 -t illumina_5 -p 8

FragGeneScan -s /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B91_S28_001_380m_mRNA.fasta -o /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B91_S28_001_380m_mRNA_annotated -w 0 -t illumina_5 -p 8

#Flags
#The -s flag is the assembly path. Please note that the output of RNASPAdes that you want to use as the input for FragGeneScan is the transcripts.fasta file in the main assembly folder. 
#The -o flag is the output directory 
#The -w flag is the read type; 0 if the sequence file has short sequence reads and 1 if the sequence file has complete genomic sequences
#The -t flag is the train file that FragGeneScan should use based on the sequencer of the reads
#The -p flag is the number of threads to use. default is 1. 
```
 
Mapping reads to ORFs using the Burrows–Wheeler Aligner-MEM (0.7.17) and Samtools (1.14)————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Programs used: Burrows-Wheeler Aligner 0.7.17 and Samtools 1.14
—-From what I understand, the mRNA reads need to be mapped to annotated coassembly, so we will have 3 align files and sam and bam files to go with each. FragGeneScan generated the following outputs: 
	•	The .out file contains a summary of the number of features annotated.
	•	The .faa file contains the protein sequences of the genes annotated.
	•	The .ffn file contains the nucleotide sequences of the genes annotated

I already had this conda environment created for HW3. 
#Slurm script for submitting bwa mem code on mRNA NUCLEOTIDE files: 
Script name in GitHub: bwa_mem_mapping_samtools_NUCLEOTIDE.sh
```
#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=bwa_mem_samtools      # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=bwa_mem_samtools.log    # Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate bwa_samtools 

#did this before running — bwa index dino_metzyme_annotated_coassembly.ffn

#bwa code to map to coassembly  

#40m
bwa mem -t 8 dino_metzyme_annotated_coassembly.ffn 30B8Z_S11_001_40m_mRNA.fasta  -o 30B8Z_S11_001_40m_mRNA_mapped.sam
#70m
bwa mem -t 8 dino_metzyme_annotated_coassembly.ffn 30B90_S12_001_70m_mRNA.fasta -o 30B90_S12_001_70m_mRNA_mapped.sam
#380m
bwa mem -t 8 dino_metzyme_annotated_coassembly.ffn 30B91_S28_001_380m_mRNA.fasta -o 30B91_S28_001_380m_mRNA_mapped.sam

#Flags 
#The -t flag specifies the number of threads to use
#The -o flag specifies the sam output name

#use samtools to convert sam to bam
#40m
samtools view -S -b 30B8Z_S11_001_40m_mRNA_mapped.sam > 30B8Z_S11_001_40m_mRNA_mapped.bam
#70m
samtools view -S -b 30B90_S12_001_70m_mRNA_mapped.sam > 30B90_S12_001_70m_mRNA_mapped.bam
#380m
samtools view -S -b 30B91_S28_001_380m_mRNA_mapped.sam > 30B91_S28_001_380m_mRNA_mapped.bam

#Flags
#The -S flag indicates you have a sam input
#The -b flag indicates you want to have a bam output

#use samtools to sort bam files
#40m
samtools sort 30B8Z_S11_001_40m_mRNA_mapped.bam -o 30B8Z_S11_001_40m_mRNA_mapped_sorted.bam
#70m
samtools sort 30B90_S12_001_70m_mRNA_mapped.bam -o 30B90_S12_001_70m_mRNA_mapped_sorted.bam
#380m
samtools sort 30B91_S28_001_380m_mRNA_mapped.bam -o 30B91_S28_001_380m_mRNA_mapped_sorted.bam

#Flags
#The -o flag specifies output file
```

#Slurm script for submitting bwa mem code on mRNA AMINO ACID files: Ultimately, we did end up using these output to calculate TPM for figure 3a.
Script name in GitHub: bwa_mem_mapping_samtool_AMINOACID.sh

```
#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=bwa_mem_samtools      # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=bwa_mem_samtools.log    # Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate bwa_samtools 

#did this before running: bwa index dino_metzyme_annotated_coassembly.faa

#bwa code to map to translated coassembly  

#40m
bwa mem -t 8 dino_metzyme_annotated_coassembly.faa 30B8Z_S11_001_40m_mRNA_annotated.faa -o 30B8Z_S11_001_40m_mRNA_translated_mapped.sam
#70m
bwa mem -t 8 dino_metzyme_annotated_coassembly.faa 30B90_S12_001_70m_mRNA_annotated.faa -o 30B90_S12_001_70m_mRNA_translated_mapped.sam
#380m
bwa mem -t 8 dino_metzyme_annotated_coassembly.faa 30B91_S28_001_380m_mRNA_annotated.faa -o 30B91_S28_001_380m_mRNA_translated_mapped.sam

#Flags 
#The -t flag specifies the number of threads to use
#The -o flag specifies the sam output name

#use samtools to convert sam to bam
#40m
samtools view -S -b 30B8Z_S11_001_40m_mRNA_translated_mapped.sam > 30B8Z_S11_001_40m_mRNA_translated_mapped.bam
#70m
samtools view -S -b 30B90_S12_001_70m_mRNA_translated_mapped.sam > 30B90_S12_001_70m_mRNA_translated_mapped.bam
#380m
samtools view -S -b 30B91_S28_001_380m_mRNA_translated_mapped.sam > 30B91_S28_001_380m_mRNA_translated_mapped.bam

#Flags
#The -S flag indicates you have a sam input
#The -b flag indicates you want to have a bam output

#use samtools to sort bam files in order to input into diamond
#40m
samtools sort 30B8Z_S11_001_40m_mRNA_translated_mapped.bam -o 30B8Z_S11_001_40m_mRNA_translated_mapped_sorted.bam
#70m
samtools sort 30B90_S12_001_70m_mRNA_translated_mapped.bam -o 30B90_S12_001_70m_mRNA_translated_mapped_sorted.bam
#380m
samtools sort 30B91_S28_001_380m_mRNA_translated_mapped.bam -o 30B91_S28_001_380m_mRNA_translated_mapped_sorted.bam

#Flags
#The -o flag specifies output file
```


Assigning taxonomy and function to ORFs using Diamond v2.0.13 with BlastP function and PhyloDB v.1.075————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Programs used: diamond v.2.0.13 PhyloDB v.1.075

Creating the conda environment: 
```
conda create —-name diamond
conda activate diamond
conda install -c bioconda diamond
```

Great resource for how to run this code: https://github.com/bb

The next step is to setup a binary DIAMOND database file that can be used for subsequent searches against the database:
```
diamond makedb --in /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/PhyloDB/phylodb_1.076.pep.fa.gz -d PhyloDB
—> —in must be a fasta file 
```

A binary file (PhyloDB.dmnd) containing the database sequences will be created in the current working directory. We will now conduct a search against this database using the translated mRNA mapped and sorted binary files as queries:
Script name in GitHub: diamond_faa_mRNA.sh
```
#!/bin/bash
#SBATCH --partition=compute            # Queue selection
#SBATCH --job-name=diamond_faa_k1    # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=12               # Number of CPU cores per task
#SBATCH --mem=50gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=diamond_faa_k1%j.log  # Standard output/error
#export OMP_NUM_THREADS=12

module load anaconda/5.1
source activate diamond

#Did this before running this script! 
diamond makedb --in /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/PhyloDB/phylodb_1.076.pep.fa.gz -d PhyloDB

#Flags
#The --in flag indicates the input database file. MUST BE A FASTA FILE
#The -d flag give that indexed database a name to refer to in the blastp scripts below

#40m
diamond blastp -q 30B8Z_S11_001_40m_mRNA_annotated.faa -d PhyloDB --evalue 0.001 -k 1 -o 30B8Z_S11_001_40m_mRNA_diamond_faa_out.tsv 
#70m
diamond blastp -q 30B90_S12_001_70m_mRNA_annotated.faa -d PhyloDB --evalue 0.001 -k 1 -o 30B90_S12_001_70m_mRNA_diamond_faa_out.tsv 
#380m
diamond blastp -q 30B91_S28_001_380m_mRNA_annotated.faa -d PhyloDB --evalue 0.001 -k 1 -o 30B91_S28_001_380m_mRNA_diamond_faa_out.tsv

#Flags
#The -q	flag indicates the query sequence
#The -d	flag points to the diamond-formatted indexed database, made using the diamond makedb function
#The --evalue specifies	the false positive rate
#The -k	flag specifies the number of hits to show. We chose one, since we only wanted the top-scoring hit. 
#The -o	flag specifies the name of the output file
```

Diamond code for coassembly against PhyloDB:
Script name in GitHub: diamond_faa_coassembly.sh 
```
#!/bin/bash
#SBATCH --partition=scavenger                   # Queue selection
#SBATCH --job-name=diamond_coassembly_faa_k1    # Job name
#SBATCH --mail-type=ALL                         # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kabbott@whoi.edu   		  # Where to send mail
#SBATCH --ntasks=1                              # Run a single task
#SBATCH --cpus-per-task=36                      # Number of CPU cores per task
#SBATCH --mem=180gb                             # Job memory request
#SBATCH --time=24:00:00                         # Time limit hrs:min:sec
#SBATCH --output=diamond_coassembly_faa_k1%j.log # Standard output/error
#export OMP_NUM_THREADS=12

module load anaconda/5.1
source activate diamond

dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/dino_metzyme_annotated_coassembly

diamond blastp -q ${dir}/dino_metzyme_annotated_coassembly.faa -d PhyloDB --evalue 0.001 -k 1 -o dino_metzyme_annotated_coassembly_diamond_out.tsv

#Flags
#The -q	flag indicates the query sequence
#The -d	flag points to the diamond-formatted indexed database, made using the diamond makedb function
#The --evalue specifies	the false positive rate
#The -k	flag specifies the number of hits to show. We chose one, since we only wanted the top-scoring hit. 
#The -o	flag specifies the name of the output file
```

Getting taxonomy from diamond output 
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
(modifying from this repository: https://github.com/Lswhiteh/phylodbannotation)

Script name in GitHub: diamond_taxonomy.py
```
srun -p compute --ntasks-per-node=2 --mem=100G --time=05:00:00 --pty bash

conda activate python_jupyter
python diamond_taxonomy.py

#contents of diamond_taxonomy.py
#!/usr/bin/env python_jupyter

import pandas as pd
import os
import glob

os.chdir('/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA')

#Adjust this depending on which files you want to attach taxonomy to
diamond_files=['dino_metzyme_annotated_coassembly_diamond_out']

full_col_names = ['Query ID', 'Subject ID', 'Percentage of identical matches', 'Alignment length', 'Number of mismatches', 'Number of gap openings', 'Start of alignment in query', 'End of alignment in query', 'Start of alignment in subject', 'End of alignment in subject','Expected value', 'Bit score']

short_col_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch','gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

### Modified from Logan Whitehouse's lab, 
# https://github.com/Lswhiteh/phylodbannotation/blob/master/fastaannotation.py

taxonomy_file = "/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/PhyloDB/phylodb_1.076.taxonomy.txt"
gene_file = "/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/PhyloDB/phylodb_1.076.annotations.txt"

tax_dict = {}
gene_dict = {}

with open(taxonomy_file) as taxfile:
    for line in taxfile:
        row = line.strip().split("\t")
        tax_dict[row[0]] = row[1:]

with open(gene_file) as genefile:
    for line in genefile:
        row = line.strip().split("\t")
        gene_dict[row[0]] = row[1:]

for dmd_file in diamond_files:
        print(f"Processing {dmd_file}")
        df_sample = pd.read_csv(dmd_file+'.tsv', sep='\t', names=short_col_names)
        gene_mapping = df_sample['sseqid'].map(gene_dict)
        ## Massage dataframe to get just the second column (organism classification)
        organism_df = pd.DataFrame(gene_mapping.values.tolist())[1]
        taxonomy_mapping = organism_df.map(tax_dict)
        tax_df = pd.DataFrame(taxonomy_mapping.values.tolist())[1]
        #splitting up taxonomy into individual groupings
        tax_df = pd.DataFrame(tax_df.str.split(';').tolist())
        tax_df.columns = ['Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Strain_name']
        df_sample_tax = pd.concat([df_sample, tax_df], axis=1)
        df_dinophyta = df_sample_tax[df_sample_tax['Phylum'] == 'Dinophyta']
        df_sample_tax.to_csv(dmd_file + '_taxonomy.tsv', sep='\t')
        df_dinophyta.to_csv(dmd_file + '_dinophyta.tsv', sep='\t')
        print(f"Done processing {dmd_file}")
~
```

Annotate genes using KEGG via GhostKOALA
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
"To maximize functional assignments, sequences were searched against the KEGG and EuKaryotic Orthologous Groups (KOG) tools and conserved protein domain families were identified using HMMER v.3.1b2 (ref. 121) with Pfam”

Based on this information, we chose to use the online web server KEGG tools. Note that GhostKOALA was chosen because it had slightly higher accuracy than KofamKOALA, and works better on large sequences than BlastKOALA (which only accepts up to 10,000 sequences for a single run).

Need to split up amino acid fasta file into two files because the limit is 300MB for an upload (and we are at ~330MB). 

#Use this code to split files up
```
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1500000==0){file=sprintf("dino_metzyme_annotated_coassembly_%d.faa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < dino_metzyme_annotated_coassembly.faa
```

Submit by adding .faa file to the online web form at https://www.kegg.jp/ghostkoala/ and enter email address

Once completed (took up to 12 hours), download 3 files (user_ko.txt, user_ko_definition.txt, user.out.top.gz) to local computer and scp back to Vortex.

Script name in GitHub: generate_KEGG_orthology.txt
```
#!/bin/bash

# Taken from Meren Lab's handy tutorial
# https://merenlab.org/2018/01/17/importing-ghostkoala-annotations/

# turns hierarchical .keg file into tab-delimited file where each row is a gene with different columsn describing layers of classification

cd /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/KEGG

wget 'https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=' -O ko00001.keg

kegfile="ko00001.keg"

while read -r prefix content
do
    case "$prefix" in A) col1="$content";; \
                      B) col2="$content" ;; \
                      C) col3="$content";; \
                      D) echo -e "$col1\t$col2\t$col3\t$content";;
    esac 
done < <(sed '/^[#!+]/d;s/<[^>]*>//g;s/^./& /' < "$kegfile") > KO_Orthology_ko00001.txt
```
 
Recombine into one file with commands similar to the following:

```
cat dino_metzyme_annotated_coassembly_0_KO.txt >> dino_metzyme_annotated_coassembly_KEGG_annotation.txt
cat dino_metzyme_annotated_coassembly_150000_KO.txt >> dino_metzyme_annotated_coassembly_KEGG_annotation.txt
```

Note that naming convention now is dino_metzyme_annotated_coassembly_KEGG_annotation.txt: KEGG annotations for each gene, empty if it could not classify. dino_metzyme_annotated_coassembly_KEGG_annotation_full.txt: KEGG annotations for each gene, including brief written description, as well as GhostKOALA score and second-best classification. 

Finding conserved protein domain families using HMMER with PFAM db and KOG db
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Programs used: HMMER 3.3.2, PFAM 35.1

Creating the conda environment:
```
conda create -n hmmer
conda activate hmmer
conda install -c bioconda hmmer
```

Download HMM profile of Pfam
```
#Download the profile and put it in: /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
```

HMMsearch just for the 3 protein families of interest: PF01036, PF03713, and PF07692
You can access these files by querying on http://pfam.xfam.org/search#tabview=tab1 by keyword. Then click on “alignments” on the lefthand of the page, scroll down to “format an alignment”, fill in the UniProt bubble, then click “Generate” to download the SELEX file for that protein. 
Script name in GitHub: hmmsearch_3_proteins.sh
```
#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=hmmsearch_pf3         # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=hmmsearch_pf3%j.log     # Standard output/error
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
```

There is a single KOG gene included in figure 3, so we also ran HMMER against that: KOG2348
Script name in GitHub: KOG_hmmsearch.sh
```
#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=hmmsearch_kog             # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kaabbott@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=hmmsearch_kog3%j.log         # Standard output/error
#export OMP_NUM_THREADS=16

module load anaconda/5.1
source activate hmmer

kog_dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/KOG
data_dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA

cd $data_dir

echo $data_dir


hmmsearch -o hmmer_output/dino_metzyme_annotated_coassembly_KOG2348_out -E 0.001 ${kog_dir}/KOG2348.hmm dino_metzyme_annotated_coassembly.faa
```

TPM normalization
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Need to index bam files before we can use pysam to read them into Python
```
conda activate bwa_samtools 
samtools index 30B8Z_S11_001_40m_mRNA_mapped_sorted.bam
samtools index 30B90_S12_001_70m_mRNA_mapped_sorted.bam
samtools index 30B91_S28_001_380m_mRNA_mapped_sorted.bam
```

TPM mapping/calculations are all in the TPM_normalization.ipynb Jupyter Notebook
```
srun -p compute --ntasks-per-node=1 --mem=10G --time=05:00:00 --pty bash
jupyter notebook --no-browser —port=8899
#Open TPM_normalization from the home directory, dinoflagellates_METZYME  
```

TPM figure (figure 3a in the paper)
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Done in python. See directions below
```
Run `conda activate python_jupyter`

Allocate an interactive job on the HPC with 
`srun -p compute --ntasks-per-node=1 --time=4:00:00 --mem=50G --pty bash` and set up a Jupyter Notebook session with `jupyter notebook --no-browser --port=8888`

- `TPM_normalization_and_plotting.ipynb` links the KEGG, KOG and Pfam annotations to the taxonomic classification TSV for dinoflagellates. It additionally uses the number of reads mapped to each ORF (BWA alignment) to calculate the transcripts per million (TPM). The z-scores for each gene of interest are calculated, and the results are plotted as a heat map for the top 50 KEGG-annotated genes identified in the Cohen et al. paper. 
```

Supergroup transcript abundance (figure 1c)
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Done in python. See directions below
```
Run `conda activate python_jupyter`

Allocate an interactive job on the HPC with 
`srun -p compute --ntasks-per-node=1 --time=4:00:00 --mem=50G --pty bash` and set up a Jupyter Notebook session with `jupyter notebook --no-browser --port=8888`

- `metaT_taxonomy.ipynb` links the output of the Diamond taxonomic classification with the PhyloDB taxonomy, so that it can be used for creating plots of relative community abundance. The outputs are two TSVs: One with the full taxonomic classification for all ORFs, and one filtered by just dinoflagellate ORFs. These outputs are used in all other Jupyter Notebooks and should be run first.

- `Reclassify_transcripts_for_relative_abundance.ipynb` splits up the full taxonomic tables into individual TSVs for different taxa used by Cohen et al. It can be modified depending on what level of taxonomic classification the reader is interested in. The outputs of this notebook are used for `community-abundance_from_transcripts.ipynb` to plot relative community abundance for both the eukaryotic and whole communities. 

- `community-abundance_from_transcripts.ipynb` takes the output of the reclassify_transcripts notebook and the output of the BWA alignment to calculate the total number of transcripts that map to each taxa of interest and normalizes it by the total number of transcripts to get relative community abundance. This relative abundance is plotted for both the eukaryotic and whole communities.
```

16S rRNA processing
Quality Controlling the data
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Programs = trim_galore v.0.6.7 + fastqc v.0.11.9

16S Trimming and FastQC on default settings
#files are only a maximum of 13MB, so ran this on srun
```
srun -p scavenger --time=04:00:00 --ntasks-per-node 2 --pty bash 
conda activate qc-trim

#Station 9 40m
trim_galore --paired --fastqc WHOI040_16S_TCGACGAG_R1.fastq.gz WHOI040_16S_TCGACGAG_R2.fastq.gz
#Station 9 70m
trim_galore --paired --fastqc WHOI020_16S_TCGACGAG_R1.fastq.gz WHOI020_16S_TCGACGAG_R2.fastq.gz
#Station 9 380m
trim_galore --paired --fastqc WHOI010_16S_AGAGTCAC_R1.fastq.gz  WHOI010_16S_AGAGTCAC_R2.fastq.gz

#The --paired flag indicates the input is paired-end reads
#The --fastqc flag tells trimgalore! to run a quality check on the reads post-trimming
```

How trimmed 16s compares in size to untrimmed:

raw total 14M
-rwx------ 1 selkassas sg-envbio-mgr 2.9M Oct 20 21:18 WHOI010_16S_AGAGTCAC_R2.fastq.gz
-rwx------ 1 selkassas sg-envbio-mgr 2.2M Oct 20 21:18 WHOI020_16S_TCGACGAG_R1.fastq.gz
-rwx------ 1 selkassas sg-envbio-mgr 2.2M Oct 20 21:18 WHOI040_16S_TCGACGAG_R2.fastq.gz
-rwx------ 1 selkassas sg-envbio-mgr 3.2M Oct 20 21:18 WHOI020_16S_TCGACGAG_R2.fastq.gz
-rwx------ 1 selkassas sg-envbio-mgr 1.9M Oct 20 21:18 WHOI010_16S_AGAGTCAC_R1.fastq.gz
-rwx------ 1 selkassas sg-envbio-mgr 1.5M Oct 20 21:18 WHOI040_16S_TCGACGAG_R1.fastq.gz

trimmed total 13M
-rwx------ 1 selkassas sg-envbio-mgr 2.0M Oct 20 21:32 WHOI040_16S_TCGACGAG_R2_val_2.fq.gz
-rwx------ 1 selkassas sg-envbio-mgr 1.4M Oct 20 21:32 WHOI040_16S_TCGACGAG_R1_val_1.fq.gz
-rwx------ 1 selkassas sg-envbio-mgr 2.9M Oct 20 21:32 WHOI020_16S_TCGACGAG_R2_val_2.fq.gz
-rwx------ 1 selkassas sg-envbio-mgr 2.1M Oct 20 21:32 WHOI020_16S_TCGACGAG_R1_val_1.fq.gz
-rwx------ 1 selkassas sg-envbio-mgr 2.6M Oct 20 21:31 WHOI010_16S_AGAGTCAC_R2_val_2.fq.gz
-rwx------ 1 selkassas sg-envbio-mgr 1.9M Oct 20 21:31 WHOI010_16S_AGAGTCAC_R1_val_1.fq.gz

Classifying 16S using KRAKEN2
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Programs used: fq2fa 1.1.3, kraken2 2.1.2 

Setting up the conda environment:
```
conda create —-name kraken2 
conda activate kraken2
conda install -c bioconda kraken2
```

Change .fastq files to .fasta files using fq2fa (IDBA_UD extension). We will need this as an input for later programs.
Script name in Github: fq2fa_16S.sh
```
#get fq2fa
git clone git@github.com:loneknightpy/idba.git

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
```

KRAKEN2 Code
Had to get the kraken formatted silva132 db from the kraken website first then unzip the tar file:
```
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/16S_Silva132_20200326.tgz
#unzip tar file
tar -xvzf 16S_Silva132_20200326.tgz
```

Code to run Kraken2: 
```
#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=kraken2               # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kaabbott@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=kraken2_paired.log             # Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate kraken2

#40m, SSU
kraken2 --db /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db/kraken_silva --output WHOI040_16S_TCGACGAG_40m_taxonomy --report 40m.kreport --paired /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI040_16S_TCGACGAG_R1_val_1.fq /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI040_16S_TCGACGAG_R2_val_2.fq

#70m, SSU
kraken2 --db /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db/kraken_silva --output WHOI020_16S_TCGACGAG_70m_taxonomy --report 70m.kreport --paired /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI020_16S_TCGACGAG_R1_val_1.fq /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI020_16S_TCGACGAG_R2_val_2.fq

#380m, SSU
kraken2 --db /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db/kraken_silva --output WHOI010_16S_AGAGTCAC_380m_taxonomy --report 380m.kreport --paired /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI010_16S_AGAGTCAC_R1_val_1.fq /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI010_16S_AGAGTCAC_R2_val_2.fq

#Flags
#The --db flag specifies the database to check reads against to id taxonomy
#The --output flag names the output file
#The --report flag should generate a kreport, though more processing is necessary to generate this.
#The —-paired flag indicated the input are two paired forward and reverse reads
```

Classifying 16S plasmids using mothur classify.seqs
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
mothur 1.46.1 

Mapping the 16S to the Protist Ribosomal Database v.4.14.0 will label plastid sequences. This is a recent update. So, we classified the 16S again using PR2, but only kept plastid sequences, denoted by :plas.

Prep the classify.seqs file 
```
nano mothur_16S_plasmids.txt
##Plasmids - 16S 

#40m
classify.seqs(fasta=WHOI040_16S_TCGACGAG_40m_merged.fasta, count=WHOI040_16S_TCGACGAG_40m_final.count_table, reference=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.fasta, taxonomy=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.tax)

#70m
classify.seqs(fasta=WHOI020_16S_TCGACGAG_70m_merged.fasta, count=WHOI020_16S_TCGACGAG_70m_final.count_table, reference=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.fasta, taxonomy=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.tax)

#380m
classify.seqs(fasta=WHOI010_16S_AGAGTCAC_380m_merged.fasta, count=WHOI010_16S_AGAGTCAC_380m_final.count_table, reference=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.fasta, taxonomy=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.tax)
```

SlURM script
```
nano mothur_16S_SLURM.sh

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_mothur_16S               # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=2                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_mothur_16S%j.log           # Standard output/error

export OMP_NUM_THREADS=12

module load anaconda/5.1
source activate mothur_1.46.1

mothur mothur_16S_plasmids.txt

#Submit to slurm: 

sbatch mothur_16S_SLURM.sh
#activates mothur then runs the commands in mothur_18S.txt
```


Subsampling 144 plastids and 2500 18S using shuf (in bash)
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

filtering out plastids from mothur plastid taxonomy
```
grep :plas WHOI040_16S_TCGACGAG_40m_merged.0_SSU_mothur.wang.taxonomy > WHOI040_16S_TCGACGAG_40m_plastids_only
grep :plas WHOI020_16S_TCGACGAG_70m_merged.0_SSU_mothur.wang.taxonomy > WHOI020_16S_TCGACGAG_70m_plastids_only
grep :plas WHOI010_16S_AGAGTCAC_380m_merged.0_SSU_mothur.wang.taxonomy > WHOI010_16S_AGAGTCAC_380m_plastids_only
```

subsampling 144 plastids
Note: 70m only has 142 plastid sequences, so did not subsample this file
40m has 1454 and 380m has 1484 
```
shuf WHOI040_16S_TCGACGAG_40m_plastids_only | head -n 144 > WHOI040_16S_TCGACGAG_40m_plastids_only_144_subsampled
shuf WHOI010_16S_AGAGTCAC_380m_plastids_only | head -n 144 > WHOI010_16S_AGAGTCAC_380m_plastids_only_144_subsampled
```

Removing plastids from 16S using grep, subsampling 6000 16S using shuf (in bash)
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

First, select first column (containing the accession numbers/sample names) from the plastid taxonomy and make it into a new file. We are going to search the full 16S taxonomy and select only sequences that DO NOT match the plastid sequences. 

```
#select first column
#40m
cut -f 1 WHOI040_16S_TCGACGAG_40m_plastids_only > WHOI040_16S_TCGACGAG_40m_plastids_only_accnos.txt
#70m
cut -f 1 WHOI020_16S_TCGACGAG_70m_plastids_only > WHOI020_16S_TCGACGAG_70m_plastids_only_accnos.txt
#380m
cut -f 1 WHOI010_16S_AGAGTCAC_380m_plastids_only > WHOI010_16S_AGAGTCAC_380m_plastids_only_accnos.txt
```

```
#search full 16S files for everything not in the plastids file
#40m
grep -vf WHOI040_16S_TCGACGAG_40m_plastids_only_accnos.txt WHOI040_16S_TCGACGAG_40m_taxonomy > WHOI040_16S_TCGACGAG_40m_no_plastids
#70m
grep -vf WHOI020_16S_TCGACGAG_70m_plastids_only_accnos.txt WHOI020_16S_TCGACGAG_70m_taxonomy > WHOI020_16S_TCGACGAG_70m_no_plastids
#380m
grep -vf WHOI010_16S_AGAGTCAC_380m_plastids_only_accnos.txt WHOI010_16S_AGAGTCAC_380m_taxonomy > WHOI010_16S_AGAGTCAC_380m_no_plastids
```

```
#shuffle sequences then subsample 6000
#40m
shuf WHOI040_16S_TCGACGAG_40m_no_plastids | head -n 6000 > WHOI040_16S_TCGACGAG_40m_no_plastids_6000_subsampled
#70m
shuf WHOI020_16S_TCGACGAG_70m_no_plastids | head -n 6000 > WHOI020_16S_TCGACGAG_70m_no_plastids_6000_subsampled
#380m
shuf WHOI010_16S_AGAGTCAC_380m_no_plastids | head -n 6000 >WHOI010_16S_AGAGTCAC_380m_no_plastids_6000_subsampled
```

Generate a parsable kraken report with taxonomic information on 16S. This will help to make the count table then the relative abundance plot (fig 1c)
Script name in GitHub: kraken_taxonomy.sh
```
#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=kraken_taxon               # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kaabbott@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=10gb                      # Job memory request
#SBATCH --time=02:00:00                  # Time limit hrs:min:sec
#SBATCH --output=kraken_taxon.log             # Standard output/error
#export OMP_NUM_THREADS=1

module load anaconda/5.1
source activate kraken2

#set directories
script_dir=/vortexfs1/home/kaabbott/.conda/pkgs/krakentools-1.2-pyh5e36f6f_0/python-scripts
database_dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db/kraken_silva
data_dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/kraken_taxonomy

cd $data_dir

for taxonomy in WHOI*_no_plastids_6000_subsampled
do
        python ${script_dir}/make_kreport.py -i $taxonomy -t ${database_dir}/16S_SILVA132_taxonomy_file -o ${taxonomy}_report
done

#rename
mv WHOI040_16S_TCGACGAG_40m_no_plastids_6000_subsampled_report 40m_16S_no_plastids_subsample_6000.kreport
mv WHOI020_16S_TCGACGAG_70m_no_plastids_6000_subsampled_report 70m_16S_no_plastids_subsample_6000.kreport
mv WHOI010_16S_AGAGTCAC_380m_no_plastids_6000_subsampled_report 380m_16S_no_plastids_subsample_6000.kreport
```

Preparing 16S taxonomy for stacked barplot figure creation
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Selecting only the taxonomy shown in the figure in the paper: 
```
#40m 
cut -f 2,6 40m_16S_no_plastids_subsample_6000.kreport | grep -E 'Bacteroidota|Betaproteobacteria|Gammaproteobacteria|Prochlorococcus|Alphaproteobacteria|Desulfobacteria|Archaea|Myxococcota|SAR324|Bdellvibrionota|Synechococcus' > 40m_16S_count_table.csv

other bacteria = 1060
other proteobacteria = 2
other cyanobacteria = 1254

#70m 
cut -f 2,6 70m_16S_no_plastids_subsample_6000.kreport | grep -E 'Bacteroidota|Betaproteobacteria|Gammaproteobacteria|Prochlorococcus|Alphaproteobacteria|Desulfobacteria|Archaea|Myxococcota|SAR324|Bdellvibrionota|Synechococcus' > 70m_16S_count_table.csv

other bacteria = 3111  
other proteobacteria = 7
other cyanobacteria =  186

#380m
cut -f 2,6 380m_16S_no_plastids_subsample_6000.kreport | grep -E 'Bacteroidota|Burkholderiales|Gammaproteobacteria|Prochlorococcus|Alphaproteobacteria|Desulfobacteria|Archaea||Myxococcota|SAR324|Bdellvibrionota|Synechococcus' > 380m_16S_count_table.csv

other bacteria = 567 
other proteobacteria = 0
other cyanobacteria = 675

#column 2 is the counts; column 6 is the taxon name
```

Note: selected Myxococcota, Desulfobacterota, SAR324, and Bdellvibrionota, since they are the new designations for deltaproteobacteria. Kraken2 also grouped all betaproteobacteria under Burkholderiales. 

To make the plot: 
Script name in Github: 16S_plot.R 
```
library(tidyverse)
library(reshape2) 
library(RColorBrewer)
library(ggplot2)
library(dplyr)

#PLOTTING
##With NEW designations
station_9_16S_tax_table <- read.delim("/Users/sabrinaelkassas/Desktop/16S/all_16S_count_table_normalized.txt", header = TRUE, sep = "\t")

#This selects the colors for the graphs
n <- 14
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


#ggplot code - 16S
ggplot(station_9_16S_tax_table, aes(x = fct_relevel(Sample, "9_40m", "9_70m", "9_380"), y = Percent, fill = Class)) + 
  geom_bar(stat = "identity", width = 0.5) + theme(axis.text.x.bottom = element_text(angle = 45)) + scale_fill_manual(values=col_vector)

#PLOTTING
##With OLD designation of deltaproteobacteria
station_9_16S_delta_tax_table <- read.delim("/Users/sabrinaelkassas/Desktop/16S/all_16S_count_table_deltaproteobacteria_grouped_normalized.txt", header = TRUE, sep = "\t")

#This selects the colors for the graphs
n <- 14
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

#ggplot code - 16S
ggplot(station_9_16S_delta_tax_table, aes(x = fct_relevel(Sample, "9_40m", "9_70m", "9_380"), y = Percent, fill = Class)) + 
  geom_bar(stat = "identity", width = 0.5) + theme(axis.text.x.bottom = element_text(angle = 45)) + scale_fill_manual(values=col_vector)
```



18S rRNA processing
Quality Controlling the data
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Programs = trim_galore v.0.6.7 + fastqc v.0.11.9

18S Trimming and FastQC on default settings 
#files are only a maximum of 32MB, so ran this on srun
Script name in GitHub: 18S_trim.sh
```
srun -p scavenger --time=04:00:00 --ntasks-per-node 2 --pty bash 
conda activate qc-trim

#Station 9 40m
trim_galore --paired --fastqc WHOI021-V9_S105_L001_R1_001.fastq WHOI021-V9_S105_L001_R2_001.fastq
#Station 9 70m
trim_galore --paired --fastqc WHOI040-V9_S124_L001_R1_001.fastq WHOI040-V9_S124_L001_R2_001.fastq
#Station 9 380m
trim_galore --paired fastqc WHOI010-V9_S94_L001_R1_001.fastq WHOI010-V9_S94_L001_R2_001.fastq

#The --paired flag indicates the input is paired-end reads
#The --fastqc flag tells trimgalore! to run a quality check on the reads post-trimming
```

How trimmed 18s compares in size to untrimmed:

raw total 11M
-rwx------ 1 kaabbott sg-envbio-mgr 2.0M Oct 14 15:37 WHOI021-V9_S105_L001_R1_001.fastq
-rwx------ 1 kaabbott sg-envbio-mgr 1.9M Oct 14 15:37 WHOI040-V9_S124_L001_R2_001.fastq
-rwx------ 1 kaabbott sg-envbio-mgr 1.2M Oct 14 15:37 WHOI010-V9_S94_L001_R1_001.fastq
-rwx------ 1 kaabbott sg-envbio-mgr 2.0M Oct 14 15:37 WHOI021-V9_S105_L001_R2_001.fastq
-rwx------ 1 kaabbott sg-envbio-mgr 1.9M Oct 14 15:37 WHOI040-V9_S124_L001_R1_001.fastq
-rwx------ 1 kaabbott sg-envbio-mgr 1.2M Oct 14 15:37 WHOI010-V9_S94_L001_R2_001.fastq

total 11M
-rwx------ 1 selkassas sg-envbio-mgr 2.0M Oct 20 21:51 WHOI021-V9_S105_L001_R2_001_val_2.fq
-rwx------ 1 selkassas sg-envbio-mgr 2.0M Oct 20 21:51 WHOI021-V9_S105_L001_R1_001_val_1.fq
-rwx------ 1 selkassas sg-envbio-mgr 1.9M Oct 20 21:51 WHOI040-V9_S124_L001_R2_001_val_2.fq
-rwx------ 1 selkassas sg-envbio-mgr 1.9M Oct 20 21:51 WHOI040-V9_S124_L001_R1_001_val_1.fq
-rwx------ 1 selkassas sg-envbio-mgr 1.2M Oct 20 21:51 WHOI010-V9_S94_L001_R2_001_val_2.fq
-rwx------ 1 selkassas sg-envbio-mgr 1.2M Oct 20 21:51 WHOI010-V9_S94_L001_R1_001_val_1.fq

-quality stats for 18s are very bad. Running through trim-galore again with a -q score of 25. 
NOTE: I tried with different length and quality cutoffs (—length 60 and -length 80 and with -q scores of 25 and 30) but did not see much of an improvement with the data. Left on default and continued with the pipeline

Classifying 18S by first merged forward and reverse reads using fq2fa then using mothur classify.seqs
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Programs used: fq2fa 1.1.3, mothur 1.46.1

Got the Protist Ribosomal Reference db v 4.14.0 from https://github.com/pr2database/pr2database/releases/tag/4.14.0
```
wget https://github.com/pr2database/pr2database/releases/download/4.14.0/pr2_version_4.11.0_merged.tsv.gz
wget https://github.com/pr2database/pr2database/releases/download/4.14.1/pr2_version_4.11.1_metadata.tsv.gz
wget https://github.com/pr2database/pr2database/releases/download/4.11.1/pr2_version_4.11.1_taxo_long.fasta.gz

NOTE: the dinoREF db is already attached to the Protist Ribosomal Reference db as of version 4.9.9.
``` 

Just in case, can get dinoREF db from https://figshare.com/ndownloader/articles/5568454/versions/2 anyways
```
wget https://figshare.com/ndownloader/articles/5568454/versions/2
```

Creating the conda environment: 
```
conda create —-name mothur_1.46.1
conda activate mothur_1.46.1
conda install -c bioconda mothur

#add fq2fa to this environment while it is activated: 
git clone git@github.com:loneknightpy/idba.git
```

First, merge forward and reverse reads and convert to fasta using fq2fa 
Script name in GitHub: fq2fa_18S.sh
```
srun -p scavenger --time=04:00:00 --ntasks-per-node 2 --pty bash 
#40m
fq2fa --merge --filter WHOI021-V9_S105_L001_R1_001_val_1.fq WHOI021-V9_S105_L001_R2_001_val_2.fq WHOI021-V9_S105_L001_40m_merged.fasta
#70m
fq2fa --merge --filter WHOI040-V9_S124_L001_R1_001_val_1.fq WHOI040-V9_S124_L001_R2_001_val_2.fq WHOI040-V9_S124_L001_70m_merged.fasta
#380m
fq2fa --merge --filter WHOI010-V9_S94_L001_R1_001_val_1.fq WHOI010-V9_S94_L001_R2_001_val_2.fq WHOI010-V9_S94_L001_380m_merged.fasta
#Flags
#The --merge flag merges the two paired reads inputted
#The —-filter flag converts the inputted fq reads to fasta
```

Creating the conda environment: 
```
conda create —-name mothur_1.46.1
conda activate mothur_1.46.1
conda install -c bioconda mothur
```

Prepping the classify.seqs script:
Script name in GitHub: mothur_18S.txt
```
nano mothur_18S.txt
#paste what’s below

##18S

#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=mothur_18S            # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=100gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=mothur_18S%j.log         # Standard output/error
#export OMP_NUM_THREADS=16

#40m
classify.seqs(fasta=WHOI021_V9_S105_L001_40m_merged.fasta, count=WHOI021_V9_S105_L001_40m_final.count_table, reference=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.fasta, taxonomy=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.tax)

#70m
classify.seqs(fasta=WHOI040_V9_S124_L001_70m_merged.fasta, count=WHOI040_V9_S124_L001_70m_final.count_table, reference=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.fasta, taxonomy=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.tax)

#380m
classify.seqs(fasta=WHOI010_V9_S94_L001_380m_merged.fasta, count=WHOI010_V9_S94_L001_380m_final.count_table, reference=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.fasta, taxonomy=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.14.0_SSU_mothur.tax)

#Flags
#The fasta input is the file you need classified taxonomically
#The count should generate a count table (it usually doesn’t thought)
#The reference is the pr2 database formatted for mothur
#The taxonomy is the pr2 tax file formatted for mothur
```

Slurm script
Script name in GitHub: mothur_18S_SLURM.sh
```
nano mothur_18S_SLURM.sh

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_mothur_18S               # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=2                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_mothur_18S%j.log           # Standard output/error

export OMP_NUM_THREADS=12

module load anaconda/5.1
source activate mothur_1.46.1

mothur mothur_18S.txt

#activates mothur then runs the commands in mothur_18S.txt

#Submit to slurm: 

sbatch mothur.sh
```

Subsampling 2,500 18S sequences using shuf in bash then head
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
subsampling 2,500 18S sequences
```
shuf WHOI021_V9_S105_L001_40m_merged.0_SSU_mothur.wang.taxonomy | head -n 2500 > WHOI021_V9_S105_L001_40m_18S_taxonomy_2500_subsampled
shuf WHOI040_V9_S124_L001_70m_merged.0_SSU_mothur.wang.taxonomy | head -n 2500 > WHOI040_V9_S124_L001_70m_18S_taxonomy_2500_subsampled
shuf WHOI010_V9_S94_L001_380m_merged.0_SSU_mothur.wang.taxonomy | head -n 2500 > WHOI010_V9_S94_L001_380m_18S_taxonomy_2500_subsampled
```


18S relative abundance stacked bar plot figure creation
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Processing subsampled tax files to prep for plotting
Script name in GitHub: 18S_processing_and_plot.R

```
library(tidyverse)
library(reshape2) 
library(RColorBrewer)
library(ggplot2)

#Read in each 18S .csv file, but will have to make formatting changes, hence the variable name "bad"
bad40m <- read.csv("/Users/sabrinaelkassas/Desktop/WHOI021_V9_S105_L001_40m_18S_taxonomy_2500_subsampled.csv", header = FALSE)
bad70m <- read.csv("/Users/sabrinaelkassas/Desktop/WHOI040_V9_S124_L001_70m_18S_taxonomy_2500_subsampled.csv", header = FALSE)
bad380m <- read.csv("/Users/sabrinaelkassas/Desktop/WHOI010_V9_S94_L001_380m_18S_taxonomy_2500_subsampled.csv", header = FALSE)

#This then reformats the columns correctly for the tax_table code
good_40m <- bad40m %>% filter(!is.na(V1)) %>% 
  select(x = V1)
good_70m <- bad70m %>% filter(!is.na(V1)) %>% 
  select(x = V1)
good_380m <- bad380m %>% filter(!is.na(V1)) %>% 
  select(x = V1)

write.csv(good_40m, "/Users/sabrinaelkassas/Desktop/WHOI021_V9_S105_L001_40m_18S_taxonomy_2500_subsampled.NEW.csv")
write.csv(good_70m, "/Users/sabrinaelkassas/Desktop/WHOI040_V9_S124_L001_70m_18S_taxonomy_2500_subsampled.NEW.csv")
write.csv(good_380m, "/Users/sabrinaelkassas/Desktop/WHOI010_V9_S94_L001_380m_18S_taxonomy_2500_subsampled.NEW.csv")
```

Generating the tax_table - Did this on Poseidon
Script name in GitHub: 18S_rRNA-tax-loop.R
#for loop to import rRNA read taxonomy assignments —> copied this to an .R file and then ran on the hpc with the slurm parameters below it. 
#Credits: Sarah Hu
##12/01/21 18S rRNA for dino_METZYME Project 

#install packages
library(tidyverse)

#Set up files
#select correct files for HPC
taxa_raw <- list.files(path = "/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/18S_reads/trimmed_18S/mothur_18S_taxonomy/subsampled_18S", pattern = "subsampled.NEW.csv", full.names = FALSE)

#path to files for HPC
path_data <- "/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/18S_reads/trimmed_18S/mothur_18S_taxonomy/subsampled_18S/"

#For loop

for(a in taxa_raw){
  #import files, use paste to string together path and file names
  imported_tax <- read.csv(paste(path_data, a, sep = ""))
  #Extract sample name from "a", and split at ".wang"
  sample_names <- unlist(strsplit(a, ".wang"))
  #Modify imported data
  output_tmp <- imported_tax %>%
    #Adding in the sample name
    mutate(SAMPLE = sample_names[1]) %>%
    #filter out unknowns
    filter(!(grepl("unknown_unclassified", x))) %>%
    select(useful = x) %>%
    separate(useful, into = c("ACCESSION_NUMBER", "taxonomy"), sep = "\t") %>%
    # use regex to modify taxonomy column
    mutate(new_tax = str_replace_all(taxonomy, pattern = "\\(\\d+\\)", replacement = "")) %>%
    # parse taxonomy lineage name by semicolon
    separate(new_tax, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
    # add artifical count column
    add_column(COUNT = 1) %>%
    # add sample information
    group_by(sample_names[1], Phylum, Class, Order, Family, Genus, Species) %>%
    summarise(SUM = sum(COUNT))
  cat("Processing...", sample_names[1], "/n/n")
  # if else statement to facilitate row bind
  if (!exists("tax_table")){
    tax_table <- output_tmp
  } else {
    tax_table <- bind_rows(tax_table, output_tmp)
  }
  rm(output_tmp)
}
#rm(output_tmp)
#run the rm(tax_table) every time you run the for-loop so it doesn't add files twice.
#rm(tax_table)

write.csv(tax_table, file  = "output-tax-table_18S.csv")
```

To submit the script to the cluster to generate the count table, used the following SLURM script: 
Script name in GitHub: count_table_SLURM.sh
```
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
```

To be done in RStudio: 
Script name in GitHub: 18S_processing_and_plot.R
```
#PLOTTING
station_9_tax_table <- read.delim("/Users/sabrinaelkassas/Desktop/output-tax-table_18S_normalized.txt", header = TRUE, sep = "\t")

#This selects the colors for the graphs
n <- 13
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

#ggplot code - 18S
ggplot(station_9_tax_table, aes(x = Sample, y = Percent, fill = Class)) + 
  geom_bar(stat = "identity", width = 0.5) + theme(axis.text.x.bottom = element_text(angle = 45)) + scale_fill_manual(values=col_vector)
```

Preparing eukaryotic taxonomy for eukaryotic transcript counts figure
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
```
grep Chlorophyta dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_chlorophyta.tsv
grep Cryptophyta dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_crytophyta.tsv
grep Haptophyta dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_haptophyta.tsv
#Diatoms and other stramenopiles:
grep Stramenopiles dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_stramenopiles.tsv
grep Bacillariophyta transcripts_stramenopiles.tsv > transcripts_diatoms.tsv
cat transcripts_diatoms.tsv | sed 's/Bacillariophyta/Diatom/g' > transcripts_diatoms.tsv 
grep -v Bacillariophyta transcripts_stramenopiles.tsv > transcripts_other_stramenopiles.tsv
grep Ciliophora dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_Ciliophora.tsv
grep Metazoa dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_metazoa.tsv
grep Dinophyta dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_dinophyta.tsv
grep Rhizaria dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_Rhizaria.tsv
grep Excavata dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_Excavata.tsv
grep Fungi dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_fungi.tsv
grep Amoebozoa dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > transcripts_amoebozoa.tsv
nano Other_eukaryota.txt
grep -vf Other_eukaryota.txt dino_metzyme_annotated_coassembly_diamond_out_taxonomy.tsv > Other_eukaryota.tsv
grep -vE ‘Bacteria|Archaea|Virus Other_eukaryota.tsv > transcripts_other_eukaryota
#Flags 
#The -E flag is regex and allow multiple grep queries at once
#The -v flag selects everything but the query
```

Eukaryotic transcript counts figure code
————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
This was done in Python. See directions below. 
```
Run `conda activate python_jupyter`

Allocate an interactive job on the HPC with 
`srun -p compute --ntasks-per-node=1 --time=4:00:00 --mem=50G --pty bash` and set up a Jupyter Notebook session with `jupyter notebook --no-browser --port=8888`

- `metaT_taxonomy.ipynb` links the output of the Diamond taxonomic classification with the PhyloDB taxonomy, so that it can be used for creating plots of relative community abundance. The outputs are two TSVs: One with the full taxonomic classification for all ORFs, and one filtered by just dinoflagellate ORFs. These outputs are used in all other Jupyter Notebooks and should be run first.

- `Reclassify_transcripts_for_relative_abundance.ipynb` splits up the full taxonomic tables into individual TSVs for different taxa used by Cohen et al. It can be modified depending on what level of taxonomic classification the reader is interested in. The outputs of this notebook are used for `community-abundance_from_transcripts.ipynb` to plot relative community abundance for both the eukaryotic and whole communities. 

- `community-abundance_from_transcripts.ipynb` takes the output of the reclassify_transcripts notebook and the output of the BWA alignment to calculate the total number of transcripts that map to each taxa of interest and normalizes it by the total number of transcripts to get relative community abundance. This relative abundance is plotted for both the eukaryotic and whole communities.
```












We chose the paper “Dinoflagellates alter their carbon and nutrient metabolic strategies across environmental gradients in the central Pacific Ocean” by Cohen et al. We propose to reanalyze the meta- genomic and transcriptomic data to regenerate the bar plots of relative community abundance shown in Figures 1b and 1c, “Protistan community composition across the METZYME transect”. We will also re-annotate the metatranscriptomic data to recreate the TPM-normalized gene expression heat map shown in Figure 3a, “Distinct dinoflagellate functional metabolism between the euphotic and mesopelagic zones of the central Pacific''. The metatranscriptomic and 16S and 18S rRNA metabarcoding data generated for this paper was from 42 seawater samples collected for biomass at 7 sites and between 3 and 13 depths per site, with the majority between 40 and 600m. To simplify the analysis, we will focus on Station 9, which has samples from three different depths. This dataset is multi-pronged and stored in many different repositories: the Proteome Xchange Consortium, Ocean Protein Portal, NSF’s Biological and Chemical Oceanography Data Management Office, and NCBI, resulting in approximately 350Gb of data. 

Cohen, Natalie R., Matthew R. McIlvin, Dawn M. Moran, Noelle A. Held, Jaclyn K. Saunders, Nicholas J. Hawco, Michael Brosnahan et al. "Dinoflagellates alter their carbon and nutrient metabolic strategies across environmental gradients in the central Pacific Ocean." Nature Microbiology 6, no. 2 (2021): 173-186.

