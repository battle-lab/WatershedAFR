source('./enrichment.R')
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied bin width value", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  outlier_tag = ""
  args[2] = "NULL"
}
if (args[2] == 'outlier'){
  outlier_tag = '_outliers'
}
binw <- as.numeric(args[1])

thresholds <- sort(unique(c(seq(.2,2,.2), seq(.5,5,.5))))
rootdir <- '/scratch/groups/abattle4/victor/WatershedAFR'
datadir <- file.path(rootdir,'data','enrichment')
outlier_tag <- ''


# African


# files needed for african variants
if (outlier_tag == '_outliers'){
  afr_z_scores_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/outlier_calling/test/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt'
} else{
afr_z_scores_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/outlier_calling/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt'
}
afr_rare_var_pairs_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants/gene-AFR-rv.txt' # table of gene-individual pairs with their rare variants
afr_pop_list_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/gtex_v8_wgs_individuals_AFR.txt'

output_dir <-file.path(datadir, paste0(basename(tools::file_path_sans_ext(afr_pop_list_file)),outlier_tag))
dir.create(output_dir,showWarnings = FALSE)

afr.output <- enrichmentTipTail(pop_subset_file = afr_pop_list_file,
                                z_scores_file = afr_z_scores_file,
                                rare_var_pairs_file = afr_rare_var_pairs_file,
                                output_dir = output_dir,
                                title.mod.pop = "African",
                                thresholds = thresholds,
                                outlier_bin_width = binw)



saveRDS(afr.output, file=file.path(datadir,paste0('african_1to5_binw', binw, outlier_tag,'.rds')))

