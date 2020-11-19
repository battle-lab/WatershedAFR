#!/usr/bin/env Rscript

## Script to create two bed files for each individual:
## 1. 10kb +/- window around the gene in the gene-individual pair
## 2. all rare variants belonging to the individual


library(dplyr)
library(data.table)
library(optparse)

# Read command line arguments
option_list = list(make_option(c('--outliers'), type = 'character', default = NULL, help = 'multi-tissue outliers with Zscore data'),
                   make_option(c('--population'), type = 'character', default = NULL, help = 'list of samples from a population'),
                   make_option(c('--regions'), type = 'character', default = NULL, help = 'bed file of 10kb +/- window around the protein coding and lincRNA coding genes'),
                   make_option(c('--indiv_at_rv'), type = 'character', default = NULL, help = 'list of individuals with each rare variant'),
                   make_option(c('--outdir'), type = 'character', default = NULL, help = 'directory to save output bed files'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

outliers_file = opt$outliers
pop_list_file = opt$population
regions_file = opt$regions
indiv_at_rv_file = opt$indiv_at_rv
outdir = opt$outdir

# open files
print('Reading files...')
outliers = read.table(outliers_file, header = TRUE, check.names = FALSE)  # multi-tissue outliers
pop_list = as.character(read.table(pop_list_file, header = FALSE)$V1) # list of individuals from a population (population list)
regions = read.table(regions_file, col.names = c("CHROM","START","END", "Gene", "SCORE"), check.names = FALSE) # bed file of 10kb +/- window around the protein coding and lincRNA coding genes
indiv_at_rv = read.table(indiv_at_rv_file, header = TRUE, check.names = FALSE) # list of individuals with each rare variant (indiv per rare variant)
print('Done')

# filter multi-tissue outliers for individuals from the population list
outliers = filter(outliers, (Ind %in% pop_list))

# filter multi-tissue outliers for genes that have at least 1 outlier individual from the population list
gene_list = unique(filter(outliers, Y == "outlier")$Gene)
outliers = filter(outliers, (Gene %in% gene_list))

# combine outlier dataframe and regions dataframe by gene
df = left_join(x = outliers, y = regions, by = "Gene")

# split dataframes (df and indiv_at_rv) by individual
df = arrange(df, Ind, CHROM, START)
indiv_at_rv = arrange(indiv_at_rv, SAMPLE)

# create output directory if it doesn't already exist
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

print('Writing list of individuals from gene-individual pairs')
indiv_list = unique(df$Ind)
indiv_list_file = paste0(outdir,'/indiv_list')
write.table(indiv_list, indiv_list_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
print(indiv_list_file)

print('Writing bed files for individual:')
for (indiv in indiv_list) {
  print(indiv)
  indiv_df = filter(df, Ind==indiv)
  indiv_rv = filter(indiv_at_rv, SAMPLE==indiv)
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
