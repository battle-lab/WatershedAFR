#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)
library(optparse)

#--------------- MAIN

## Read command line arguments
option_list = list(make_option(c('--MAP'), type = 'character', default = NULL, help = 'path to sample-to-tissue map'),
                   make_option(c('--POP.LIST'), type = 'character', default = NULL, help = 'path to list of individuals from a population'),
                   make_option(c('--OUT'), type = 'character', default = NULL, help = 'output file name'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

map_file = opt$MAP
pop_list_file = opt$POP.LIST
outfile = opt$OUT


# read files
map = read.table(map_file, header = F, stringsAsFactors = F, col.names = c('Sample', 'Tissue'))

pop_list = read.table(pop_list_file, header = F, stringsAsFactors = F)[,1]

# subset by population
map = map %>% separate(col = Sample, into = c('GTEX', 'ID'), sep = '-', extra = 'drop', remove = F) %>%
  unite(col = 'GTEX-ID', GTEX:ID, sep = '-', remove = T) %>% filter(`GTEX-ID` %in% pop_list) %>%
  select(c('Sample', 'Tissue'))

# save
write.table(map, outfile, quote = F, sep = '\t', row.names = F, col.names = F)

