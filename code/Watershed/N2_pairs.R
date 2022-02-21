library(data.table)
library(dplyr)

rv_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.txt'

rv = fread(rv_file)
rv = rv %>% rename(GeneName = Gene, SubjectID = Ind)

set.seed(12)
rv = rv[sample(nrow(rv)),]

by_gene_rv = rv %>% group_by(GeneName, Chrom, Start, Ref, Alt) 

# Get N2 pairs (2+ individuals sharing rare variant nearby a gene)
n2_pair = by_gene_rv %>% filter(n() >= 2)
n2_pair_id = n2_pair %>% group_indices()
n2_pair = n2_pair %>% ungroup %>% mutate(N2pair = n2_pair_id)

by_gene_ind = rv %>% group_by(SubjectID, GeneName) %>% group_keys() 
out = left_join(by_gene_ind, n2_pair, by = c("SubjectID", "GeneName"))
out_not_n2 = out %>% filter(is.na(N2pair))
out_n2 = out %>% filter(!is.na(N2pair))

# In a given N2 pair, select exactly two individuals and ignore all other individuals
out_n2 = out_n2 %>% group_by(N2pair) %>% slice(1:2) %>% ungroup()

out = bind_rows(out_n2, out_not_n2) %>% select(c("SubjectID", "GeneName", "N2pair"))

outfile = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.N2pairs.tsv'
write.table(out, outfile, quote = F, sep = '\t', row.names = F, col.names = T)
