source('./enrichment.R')

thresholds <- sort(unique(c(seq(1,2,.2), seq(1,5,.5))))

rootdir <- '/scratch/groups/abattle4/victor/WatershedAFR'
datadir <- file.path(rootdir,'data','enrichment')
outlier_tag <- ''

# European

# files needed for european variants
euro_z_scores = '/scratch/groups/abattle4/victor/WatershedAFR/data/outlier_calling/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt'
# table of gene-individual pairs with their rare variants
euro_rare_var_pairs =
  # '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants/old/gene-EUR-rv.txt' # table of gene-individual pairs with their rare variants
  '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants/gene-EUR-rv.txt'

euro_pop_list = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/gtex_v8_wgs_individuals_EUR.txt'  #list of samples from the desired population

output_dir <-file.path(datadir, paste0(basename(tools::file_path_sans_ext(euro_pop_list)),outlier_tag))
dir.create(output_dir,showWarnings = FALSE)

euro.output <- enrichmentTipTail(pop_subset_file = euro_pop_list, 
                                 z_scores_file = euro_z_scores, 
                                 rare_var_pairs_file = euro_rare_var_pairs , 
                                 output_dir = datadir, title.mod.pop = "European",
                                 thresholds = thresholds)

saveRDS(euro.output, file=paste0(datadir,'european',outlier_tag, '_1to5.rds'))
