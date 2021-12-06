#for loop to import rRNA read taxonomy assignments
##Input path to all data files
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
