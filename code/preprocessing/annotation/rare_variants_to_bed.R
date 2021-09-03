# converts rare variants file to bed input for querying PhyloP 100way conservation scores from UCSC

library(data.table)
library(dplyr)
library(stringr)
library(optparse)

# Read command line arguments
option_list = list(make_option(c('--RV'), type = 'character', default = NULL, help = 'Rare variant file'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

rv.file = opt$RV
#rv.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.txt'

# open rare variants
rv = fread(rv.file)

# select relevant columns
rv = rv %>% select(Chrom, Start, End)

# convert chromosome to number only
chrom_num = str_remove(rv$Chrom, 'chr')
rv$Chrom = as.numeric(chrom_num)

# sort and remove duplicate rows
rv = rv %>% arrange(Chrom, Start) %>% distinct()

# add the 'chr' prefix back to chromosome column
rv = rv %>% mutate(Chrom = paste0('chr',Chrom))

# write to same folder
outfile = paste0(tools::file_path_sans_ext(rv.file), '.bed')
write.table(rv, outfile, quote = F, sep = '\t', row.names = F, col.names = F)