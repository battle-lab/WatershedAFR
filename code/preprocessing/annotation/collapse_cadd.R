#!/usr/bin/env Rscript

## Collapse CADD annotations to gene-individual pair level
## Take maximum across multiple rare variants

library(data.table)
library(dplyr)
library(tidyr)

rv.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.txt'
out.file ='/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.CADD.collapse.tsv'

rv = fread(rv.file)
header = TRUE
for (chrom in 1:22){
  cadd.file = paste0('/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.CADD.chr', chrom, '.tsv.gz')
  cadd = fread(cadd.file) %>%
    rename("Chrom" = "#Chrom") %>%
    mutate(Chrom = paste0("chr",Chrom))
  
  # Convert genomic coordinates in CADD annotations to be consistent with rare variants (0-based)
  cadd$Pos = cadd$Pos - 1
  # Match rare variants with their annotations and rename columns
  cadd.collapse = filter(rv, Chrom == paste0("chr",chrom)) %>%
    left_join(., cadd, by = c("Chrom" = "Chrom", "Start" = "Pos")) %>%
    rename(GeneName = Gene, SubjectID = Ind)


# Gene-Ind level transformation
cadd.collapse = cadd.collapse %>% group_by(SubjectID, GeneName) %>%
  summarise(GC=max(GC),
            CpG=max(CpG),
            SIFTcat_deleterious=max(SIFTcat_deleterious),
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
            PHRED=max(PHRED)) %>%
  select(SubjectID, everything(.))

# Impute missing values
cadd.collapse = cadd.collapse %>%
  mutate(GC=replace_na(GC, 0.418),
         CpG=replace_na(CpG, 0.024),
         SIFTcat_deleterious=replace_na(SIFTcat_deleterious,0),
         SIFTcat_tolerated=replace_na(SIFTcat_tolerated,0),
         SIFTval=replace_na(SIFTval, 0),
         PolyPhenCat_benign=replace_na(PolyPhenCat_benign, 0),
         PolyPhenCat_possibly_damaging=replace_na(PolyPhenCat_possibly_damaging, 0),
         PolyPhenCat_unknown=replace_na(PolyPhenCat_unknown, 0),
         PolyPhenVal=replace_na(PolyPhenVal,0),
         bStatistic=replace_na(bStatistic, 800.261),
         priPhCons=replace_na(priPhCons, 0.115),
         mamPhCons=replace_na(mamPhCons, 0.079),
         verPhCons=replace_na(verPhCons, 0.094),
         priPhyloP=replace_na(priPhyloP, -0.033),
         mamPhyloP=replace_na(mamPhyloP, -0.038),
         verPhyloP=replace_na(verPhyloP, 0.017),
         GerpN=replace_na(GerpN, 1.909),
         GerpS=replace_na(GerpS, -0.2),
         PHRED=replace_na(PHRED, 0))
         

fwrite(cadd.collapse, out.file, quote = F, sep = '\t', row.names = F, col.names = header, append=!header)
header = FALSE
}