#!/usr/bin/env Rscript

## Collapse PhyloP (UCSC conservation scores) to gene-individual pair level
## Take maximum across multiple rare variants

library(data.table)
library(dplyr)

rv.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.txt'
phylop.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.phyloP100way.bed'

rv = fread(rv.file)
phylop = fread(phylop.file, col.names=c("Chrom", "Start", "End", "phylop"))

# Match rare variants with their annotations
phylop.collapse = left_join(rv, phylop, by = c("Chrom", "Start"))

# Rename columns
phylop.collapse = phylop.collapse %>% rename(GeneName = Gene, SubjectID = Ind)

# Gene-Ind level transformation
phylop.collapse = phylop.collapse %>% group_by(GeneName, SubjectID) %>%
  summarise(phylop=max(phylop)) %>%
  relocate("SubjectID")

outfile = paste0(tools::file_path_sans_ext(phylop.file), '.collapse.tsv')
write.table(phylop.collapse, outfile, quote = F, sep = '\t', row.names = F, col.names = T)
