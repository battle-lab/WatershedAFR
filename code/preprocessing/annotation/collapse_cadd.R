#!/usr/bin/env Rscript

## Collapse CADD annotations to gene-individual pair level
## Take maximum across multiple rare variants

library(data.table)
library(dplyr)

rv.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.txt'
out.file ='/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.CADD.collapse.tsv'

rv = fread(rv.file)
header = TRUE
for (chrom in 1:22){
  cadd.file = paste0('/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.CADD.chr', chrom, '.tsv.gz')
  cadd = fread(cadd.file) %>%
    rename("Chrom" = "#Chrom") %>%
    mutate(Chrom = paste0("chr",Chrom))
  
  # Convert genomic coordinates in vep.loftee to be consistent with rare variants (0-based)
  cadd$Pos = cadd$Pos - 1
  # Match rare variants with their annotations and rename columns
  cadd.collapse = left_join(rv, cadd, by = c("Chrom" = "Chrom", "Start" = "Pos")) %>%
    rename(GeneName = Gene, SubjectID = Ind)


# Gene-Ind level transformation
cadd.collapse = cadd.collapse %>% group_by(SubjectID, GeneName) %>%
  summarise(SIFTcat_deleterious=max(SIFTcat_deleterious),
            SIFTcat_tolerated=max(SIFTcat_tolerated), 
            SIFTval=max(SIFTval),
            PolyPhenCat_benign=max(PolyPhenCat_benign),
            PolyPhenCat_possibly_damaging=max(PolyPhenCat_possibly_damaging),
            PolyPhenCat_unknown=max(PolyPhenCat_unknown),
            PolyPhenVal=max(PolyPhenVal),
            bStatistic=max(bStatistic),
            priPhCons=max(priPhCons),
            mamPhCons=max(mamPhCons),
            verPhCons=max(verPhCons),
            priPhyloP=max(priPhyloP),
            mamPhyloP=max(mamPhyloP),
            verPhyloP=max(verPhyloP),
            GerpN=max(GerpN),
            GerpS=max(GerpS),
            RawScore=max(RawScore)) %>%
  select(SubjectID, everything(.))
  # relocate("SubjectID")

#outfile = paste0(tools::file_path_sans_ext(vep.loftee.file), '.collapse.tsv')
fwrite(cadd.collapse, out.file, quote = F, sep = '\t', row.names = F, col.names = header, append=!header)
header = FALSE
}