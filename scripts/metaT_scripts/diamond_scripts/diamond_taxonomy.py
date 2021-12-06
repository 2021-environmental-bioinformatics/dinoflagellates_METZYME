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
