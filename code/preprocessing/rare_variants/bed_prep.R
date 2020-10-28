#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

# Read command line arguments
outliers_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/outlier_calling/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt'
pop_list_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/gtex_v8_wgs_individuals_EUR.txt'
regions_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed'
indiv_per_rv_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants/gtex_EUR_rare.QC.indiv.txt'
outdir = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants/rv_bed_EUR'


# open files
outliers = read.table(outliers_file, header = TRUE, check.names = FALSE)  # multi-tissue outliers
pop_list = as.character(read.table(pop_list_file, header = FALSE)$V1) # list of samples from a population (population list)
regions = read.table(regions_file, col.names = c("CHROM","START","END", "Gene", "SCORE", "STRAND"), check.names = FALSE) # bed file of 10kb with TSS of protein coding and lincRNA coding genes (regions)
indiv_per_rv = read.table(indiv_per_rv_file, header = TRUE, check.names = FALSE) # list of samples with each rare variant (indiv per rare variant)

# filter multi-tissue outliers for individuals from the population list
outliers = filter(outliers, (Ind %in% pop_list))

# filter multi-tissue outliers for genes that have at least 1 outlier individual from the population list
gene_list = unique(filter(outliers, Y == "outlier")$Gene)
outliers = filter(outliers, (Gene %in% gene_list))

# combine outlier dataframe and regions dataframe by gene
df = left_join(x = outliers, y = regions, by = "Gene")

# split dataframes (df and indiv_per_rv) by individual
df = arrange(df, Ind, CHROM, START)
indiv_per_rv = arrange(indiv_per_rv, SAMPLE)

indiv_list = unique(df$Ind)
write.table(indiv_list, paste0(outdir,'/indiv_list'), col.names = FALSE, row.names = FALSE, quote = FALSE)

for (indiv in indiv_list) {
  indiv_df = filter(df, Ind==indiv)
  indiv_rv = filter(indiv_per_rv, SAMPLE==indiv)
  indiv_rv$END = indiv_rv$POS + 1
  
  #gene-indiv bed
  gene.indiv.bed = indiv_df[,c('CHROM','START','END','Gene','Ind')]
  gene.indiv.bed_file = paste0(outdir,'/',indiv,'.gene-indiv.bed')
  write.table(gene.indiv.bed, gene.indiv.bed_file, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  # rv bed
  rv.bed = indiv_rv[,c('CHROM', 'POS', 'END', 'REF', 'ALT', 'SAMPLE')]
  rv.bed_file = paste0(outdir,'/',indiv,'.rv.bed')
  write.table(rv.bed, rv.bed_file, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}
