library(data.table)
library(dplyr)
library(stringr)
library(optparse)
library(tidyr)

# Read command line arguments
option_list = list(
  make_option(c('--cadd'), type = 'character', default = NULL, help = 'Annotation output file from CADD-scripts'),
  make_option(c('--out'), type = 'character', default = NULL, help = 'Annotation output after munging. it should be a tsv.gz'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

annot.file = opt$cadd
out.name = opt$out

# "GC"
# "CpG"
# "SIFTval"
# "PolyPhenVal"
# "priPhCons"
# "mamPhCons"
# "verPhCons"
# "priPhyloP"
# "mamPhyloP"
# "verPhyloP"
# "bStatistic"
# "GerpN"
# "GerpS"
# "RawScore"
# "PHRED"

desired_cols = c( "Consequence", "ConsScore", "ConsDetail", "GC", "CpG", "SIFTcat", "SIFTval", "PolyPhenCat", "PolyPhenVal",
                 "priPhCons", "mamPhCons", "verPhCons", "priPhyloP", "mamPhyloP", "verPhyloP", "bStatistic", "GerpRS",
                 "GerpN", "GerpS", "PHRED", "RawScore")
# annot.file="/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.CADD.chr12.tsv.gz"

# open output from cadd-scripts
annot <-fread(annot.file, skip = 1, header = TRUE) 
annot <- annot[ , ..desired_cols]

# combine possibly_damaging and probably_damaging categories in PolyPhenCat
annot$PolyPhenCat[annot$PolyPhenCat == "probably_damaging"] <- "possibly_damaging"

out.annot <- annot %>% 
  tidyr::pivot_wider( names_from=PolyPhenCat, 
                      values_from=PolyPhenCat, 
                      names_prefix="PolyPhenCat_", 
                      values_fn=function(x){1}, 
                      values_fill=0) %>%
  tidyr::pivot_wider( names_from=SIFTcat, 
                      values_from=SIFTcat, 
                      names_prefix="SIFTcat_", 
                      values_fn=function(x){1}, 
                      values_fill=0)

fwrite(out.annot, file=out.name, sep='\t', compress="gzip")
