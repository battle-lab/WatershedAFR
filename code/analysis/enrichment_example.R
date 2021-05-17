source('./enrichment.R')

# List of individuals used in analysis. This example is list of all GTEx Europeans that have WGS data
pop_subset_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/gtex_v8_wgs_individuals_EUR.txt'

# Multi-tissue outliers with columns:
# Ind - Name of individual
# Gene - Gene name
# N - Number of individuals with this gene (not neccessary for enrichment)
# Df - Number of tissues that have gene expression for this gene-individual pair (not neccessary for enrichment)
# MedZ - Median gene expression across tissues tissues
# Y - outlier status (not neccessary for enrichment)
z_scores_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/outlier_calling/WATERSHED.gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt'

# table of gene-individual pairs with their rare variants
rare_var_pairs_file = '/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/all_rare_variants_SNPs_10kb_genebody.txt' 
#column names for rare_var_pairs_file, not needed if rare_var_pairs_file already has column names
col.names <- c("Ind", "Gene", "Chrom", "Start", "End") 

# Outlier thresholds to compute enrichment at
thresholds = seq(1,5,0.5)

# Compute enrichment
enrich.output <- enrichmentTipTail(pop_subset_file = pop_subset_file, z_scores_file = z_scores_file, 
                                   rare_var_pairs_file = rare_var_pairs_file, thresholds = thresholds, 
                                   output_dir = "./", title.mod.pop = "Watershed European", 
                                   col.names = col.names,
                                   save_thresh_files = FALSE)


enrich.output
