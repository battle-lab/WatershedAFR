source('./enrichment.R')

thresholds <- sort(unique(c(seq(1,2,.2), seq(1,5,.5))))
rootdir <- '/scratch/groups/abattle4/victor/WatershedAFR'
datadir <- paste0(rootdir,'/','data/')

# African




# files needed for african variants
afr_z_scores_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/outlier_calling/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt'
afr_rare_var_pairs_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants/gene-AFR-rv.txt' # table of gene-individual pairs with their rare variants
afr_pop_list_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/gtex_v8_wgs_individuals_AFR.txt'


afr.output <- enrichmentTipTail(pop_subset_file = afr_pop_list_file,
                                z_scores_file = afr_z_scores_file,
                                rare_var_pairs_file = afr_rare_var_pairs_file,
                                output_dir = datadir,
                                title.mod.pop = "African",
                                thresholds = thresholds,
                                outlier_bin_width = 0.5)



saveRDS(afr.output, file=paste0(datadir,'african_1to5_binw05.rds'))

