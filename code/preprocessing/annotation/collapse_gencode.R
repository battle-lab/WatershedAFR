#!/usr/bin/env Rscript

## Collapse gencode annotations to gene-individual pair level
## Take minimum across multiple rare variants

library(data.table)
library(dplyr)

gencode.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.gencode.txt'

gencode = fread(gencode.file)

gencode.collapse = gencode %>% group_by(Gene, Ind) %>% 
  summarise(distTSS=min(distTSS), distTES=min(distTES)) %>%
  select("Ind", "Gene", "distTSS", "distTES")

colnames(gencode.collapse) = c("SubjectID", "GeneName", "distTSS", "distTES")

outfile = paste0(tools::file_path_sans_ext(gencode.file), '.collapse.tsv')
write.table(gencode.collapse, outfile, quote = F, sep = '\t', row.names = F, col.names = T)
