source('./enrichment.R')

thresholds <- sort(unique(c(seq(.2,2,.2), seq(1,5,.5))))
rootdir <- '/scratch/groups/abattle4/victor/WatershedAFR'
datadir <- paste0(rootdir,'/','data/enrichment/')
outlier_tag <- ''

# African


# files needed for african variants
afr_z_scores_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/outlier_calling/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt'
afr_rare_var_pairs_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants/gene-AFR-rv.txt' # table of gene-individual pairs with their rare variants
afr_pop_list_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/gtex_v8_wgs_individuals_AFR.txt'


output_dir <-file.path(datadir, basename(tools::file_path_sans_ext(afr_pop_list_file)))
dir.create(output_dir,showWarnings = FALSE)



afr.output <- enrichmentTipTail(pop_subset_file = afr_pop_list_file,
                                z_scores_file = afr_z_scores_file,
                                rare_var_pairs_file = afr_rare_var_pairs_file,
                                output_dir = output_dir,
                                title.mod.pop = "African",
                                thresholds = thresholds)



saveRDS(afr.output, file=paste0(datadir,'/african_1to5.rds'))

