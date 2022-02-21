#!/usr/bin/env Rscript

## Collapse GTEx annotations to gene-individual pair level
## The annotations and their transformations across multiple rare variants
## - af (min)
## - num_rare_variants (sum)
## TODO
## - donor_ss (max)
## - acceptor_ss (max)
## - ppt_region (max)
## - donor_ss_window (max)
## - acceptor_ss_winow (max)

library(data.table)
library(dplyr)

rv.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.txt'
rv = fread(rv.file)

# Rename columns
rv = rv %>% rename(GeneName = Gene, SubjectID = Ind)

# Gene-Ind level transformation
rv = rv %>% group_by(SubjectID, GeneName) %>%
  summarise(af=min(AF), num_rare_variants=n()) %>%
  select(c("SubjectID", "GeneName", "af", "num_rare_variants"))

outfile = paste0('/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.gtex_anno.collapse.tsv')
write.table(rv, outfile, quote = F, sep = '\t', row.names = F, col.names = T)
