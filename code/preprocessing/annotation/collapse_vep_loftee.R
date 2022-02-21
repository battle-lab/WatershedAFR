#!/usr/bin/env Rscript

## Collapse VEP and LOFTEE annotations to gene-individual pair level
## Take maximum across multiple rare variants

library(data.table)
library(dplyr)

rv.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.txt'
vep.loftee.file = '/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.vep.loftee.snv.tsv'

rv = fread(rv.file)
vep.loftee = fread(vep.loftee.file)

# Convert genomic coordinates in vep.loftee to be consistent with rare variants (0-based)
vep.loftee$Pos = vep.loftee$Pos - 1

# Match rare variants with their annotations
vep.loftee.collapse = left_join(rv, vep.loftee, by = c("Chrom" = "Chrom", "Start" = "Pos"))

# Rename columns
vep.loftee.collapse = vep.loftee.collapse %>% rename(GeneName = Gene, SubjectID = Ind)

# Gene-Ind level transformation
vep.loftee.collapse = vep.loftee.collapse %>% group_by(SubjectID, GeneName) %>%
  summarise(`3_prime_UTR_variant`=max(`3_prime_UTR_variant`), `5_prime_UTR_variant`=max(`5_prime_UTR_variant`), 
            TF_binding_site_variant=max(TF_binding_site_variant), downstream_gene_variant=max(downstream_gene_variant),
            intergenic_variant=max(intergenic_variant), intron_variant=max(intron_variant), missense_variant=max(missense_variant),
            non_coding_transcript_exon_variant=max(non_coding_transcript_exon_variant),
            regulatory_region_variant=max(regulatory_region_variant), splice_acceptor_variant=max(splice_acceptor_variant),
            splice_donor_variant=max(splice_donor_variant), splice_region_variant=max(splice_region_variant),
            stop_gained=max(stop_gained), synonymous_variant=max(synonymous_variant),
            upstream_gene_variant=max(upstream_gene_variant), LoF_HC=max(LoF_HC), LoF_LC=max(LoF_LC)) %>%
  relocate("SubjectID")

outfile = paste0(tools::file_path_sans_ext(vep.loftee.file), '.collapse.tsv')
write.table(vep.loftee.collapse, outfile, quote = F, sep = '\t', row.names = F, col.names = T)
