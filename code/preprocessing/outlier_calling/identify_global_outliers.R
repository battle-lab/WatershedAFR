#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
# library(ggplot2)
library(RColorBrewer)
library(optparse)

option_list = list(make_option(c('--OUTLIERS'), type = 'character', default = NULL, help = 'path to the Z-score data'),
                   make_option(c('--METHOD'), type = 'character', default = 'proportion', help = 'indicate whether to determine global outlier threshold based on proportion of outliers per person or number, options = [proportion, number]'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

outlier_file = opt$OUTLIERS
global_method = opt$METHOD

medz_data = fread(outlier_file, data.table=F)
medz_outliers = filter(medz_data, Y == 'outlier')

# Look at count of outliers per individual
medzCount = as.data.frame(table(medz_outliers$Ind))

# Get count of genes per individual
outliers_per_ind = as.data.frame(table(medz_data$Ind))
medzCount = merge(medzCount,outliers_per_ind, by='Var1') %>% mutate(PropOut = Freq.x/Freq.y)

# Remove individuals with more outliers than Q3 + 1.5*IQ 
if (global_method == 'proportion') {
  q3 = quantile(medzCount$PropOut)[4] 
  q1 = quantile(medzCount$PropOut)[2]
  qthres = q3 + 1.5*(q3-q1)
  indsRemove = unique(filter(medzCount, PropOut > qthres)$Var1)
} else {
  q3 = quantile(medzCount$Freq.x)[4] 
  q1 = quantile(medzCount$Freq.x)[2]
  qthres = q3 + 1.5*(q3-q1)
  indsRemove = unique(filter(medzCount, Freq.x > qthres)$Var1)
}

medz_data = filter(medz_data, !(Ind %in% indsRemove))

out_file = paste0(strsplit(basename(outlier_file), '.txt')[[1]][1], '_globalOutliersRemoved.txt')
write.table(medz_data, file=paste0(dirname(outlier_file), '/', out_file), sep='\t', quote=F, row.names=F)




