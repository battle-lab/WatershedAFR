# standardize all of the columns except the ones called GeneName and SubjectID

library(data.table)
library(dplyr)
library(stringr)
library(optparse)

# Read command line arguments
option_list = list(make_option(c('--ANNOT'), type = 'character', default = NULL, help = 'File with annotation columns'), 
                   make_option(c('--OUTFILE'), type = 'character', default = NULL, help = 'Name of desired outputfile'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

input.file = opt$ANNOT
if (is.na(opt$OUTFILE)){
  output.file = paste0(tools::file_path_sans_ext(input.file), '_std',tools::file_ext(input.file))
} else{
output.file = opt$OUTFILE
}

# open input file
input = fread(input.file)

# select relevant columns
input2 <- input %>%
  mutate_at(vars(-SubjectID, -GeneName), ~(scale(., center=TRUE, ) %>% as.vector))

fwrite(input,file=output.file, sep="\t", col.names = TRUE)

