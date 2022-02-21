# add pvalues using zscore column provided

library(data.table)
library(dplyr)
library(stringr)
library(optparse)

# Read command line arguments
option_list = list(make_option(c('--ZFILE'), type = 'character', default = NULL, help = 'File with Z-scores'), 
                   make_option(c('--ZCOL'), type = 'character', default = "MedZ", help = 'Name of column holding zscores'),
                   make_option(c('--OUTFILE'), type = 'character', default = NULL, help = 'Name of desired outputfile'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

input.file = opt$ZFILE
#input.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt'
output.file = opt$OUTFILE
output.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved_PVAL.txt'

# open input file
input = fread(input.file)

# select relevant columns
input$PVAL = pnorm(input$MedZ, lower.tail=FALSE)

fwrite(input,file=output.file, sep="\t")

