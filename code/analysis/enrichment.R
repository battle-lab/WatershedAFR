library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)


summarizeRareInlierOutlier <- function(pop_subset_file, z_scores_file, rare_var_pairs_file, 
                                       thresholds=seq(1,5,0.5), outlier_bin_width = NULL, 
                                       outlier_bin_percentile=NULL,
                                       output_dir = "./", output_prefix=NULL, col.names = NULL,
                                       save_thresh_files = TRUE){
  names(thresholds) <- thresholds
  
  if (is.null(output_prefix )){
    output_prefix <- basename(tools::file_path_sans_ext(pop_subset_file))
  }
  if (! is.null(outlier_bin_width)){
    output_prefix <- paste0(output_prefix,"_binw",outlier_bin_width)
  }
  
  
  pop_list <- as.character(fread(pop_subset_file, header = FALSE)$V1)
  
  if (! is.null(col.names)){
    rare_var_pairs <- fread(rare_var_pairs_file, header=FALSE) 
    colnames(rare_var_pairs) <- col.names
  }
  else {
    rare_var_pairs <- fread(rare_var_pairs_file, header=TRUE) 
  }
  
  rare_var_pairs <- rare_var_pairs %>%
    mutate(easykey = paste0(Gene,Ind))
  # remove "duplicate" gene-individual pairs with rare variants NOTE: Only first occurance in table is kept
  rare_var_pairs <- rare_var_pairs[!duplicated(rare_var_pairs$easykey),] 
  
  # create a subtable that will be more useful ## DO NOT CHANGE ITS ORDER AFTER THIS STEP
  pairs.df.pop <- fread(z_scores_file, header=TRUE) %>%
    filter(Gene %in% rare_var_pairs$Gene) %>%
    mutate(easykey = paste0(Gene,Ind)) %>% # create an id column for gene-indiv pairs
    select(easykey,Gene,Ind,MedZ) %>%
    filter(Ind %in% pop_list) %>% #limit pairs.df.pop to individuals within the population
    mutate(has_rare=ifelse(easykey %in% rare_var_pairs$easykey,1,0))  %>% # add column to pairs.df.pop indicating if the given gene-individual pair has at least 1 rare variant within 10kb of gene TSS
    ungroup()
  gc()
  
  maximum_z <- max(abs(pairs.df.pop$MedZ))
  #if we don't want to use an outlier bin width, set it to the maximum zscore
  if (is.null(outlier_bin_width )){
    thresh_upper_bound = maximum_z
  }
  
  if (! is.null(outlier_bin_percentile)){
    qthresh <- quantile(pairs.df.pop$MedZ[abs(pairs.df.pop$MedZ) > as.numeric(t)], outlier_bin_percentile)
  }
  
  # frequency table to count inliers and outliers with or without rare variants
  counts.list <- lapply(thresholds, function(t){
    print(t)
    if (! is.null(outlier_bin_width )){
      thresh_upper_bound = t + outlier_bin_width
    }
    if (! is.null(outlier_bin_percentile)){
      thresh_upper_bound = max(thresh_upper_bound,qthresh)
    }
    
    gktmp <- 
      pairs.df.pop %>%
      mutate(outlier = ifelse( abs(pairs.df.pop$MedZ) > as.numeric(t),
                               ifelse(abs(pairs.df.pop$MedZ) <= thresh_upper_bound, 
                                      "outlier","outofbounds"),"inlier")) %>%
      group_by(Gene) %>% 
      mutate(
        OutlierGenePairs = sum(outlier == "outlier",na.rm = T),
        GeneKeep = OutlierGenePairs > 0) %>% # limit analysis to genes with at least 1 outlier individual
      ungroup() %>%
      filter(GeneKeep)
    # cbind(pairs.df.pop, all.pairs.thresh[[t]]) %>% 
    # filter(GeneKeep)
    if (save_thresh_files) {
      write.table(gktmp, file = paste0(output_dir,"/",output_prefix,"_",t,".tsv"), quote = F, sep = '\t', row.names = FALSE, col.names = TRUE)
    }
    
    gc()
    gktmp <- gktmp %>%
      filter(GeneKeep) 
    
    tmp1 <- gktmp %>%
      count(outlier) %>%
      pivot_wider(names_from=outlier,values_from=n) %>%
      rename_at(which(colnames(.) %in% c("outlier","inlier","outofbounds")), ~paste0(.,"_all_pairs"))
    
    tmp2 <- gktmp %>% 
      filter(has_rare == 1) %>%
      count(outlier) %>%
      pivot_wider(names_from=outlier,values_from=n) %>%
      rename_at(which(colnames(.) %in% c("outlier","inlier","outofbounds")), ~paste0(.,"_rare_pairs"))
    
    tmp <- cbind(tmp1,tmp2)
    rm(tmp1,tmp2)
    tmp$thresh <- t
    tmp$outlier_bin_width <- outlier_bin_width
    tmp$outlier_bin_width <- outlier_bin_percentile
    tmp$outlier_upper_thresh <- thresh_upper_bound
    gc()
    return(tmp)
  })
  # frequency table to count inliers and outliers with or without rare variants
  summary.df <- bind_rows(counts.list)
  
  return(summary.df)
}
# 
# computeEnrichment <- function(summary.df){
#   
#   # compute enrichment
#   enrichment.df  <- summary.df %>% mutate(numerator=outlier_rare_pairs/outlier_all_pairs, denominator=inlier_rare_pairs/inlier_all_pairs, enrichment=numerator/denominator) %>% 
#     # confidence intervals
#     mutate(log_se = sqrt(1/outlier_rare_pairs - 1/outlier_all_pairs + 1/inlier_rare_pairs - 1/inlier_all_pairs), lower_CI = enrichment * exp(-1.96*log_se), upper_CI = enrichment * exp(1.96*log_se)) %>%
#     mutate(sortcol = as.numeric(thresh)) %>%
#     arrange(sortcol)
#   
#   return(enrichment.df)
# }

computeEnrichment <- function(summary.df){
  
  # compute enrichment
  enrichment.df  <- summary.df %>% 
    mutate( numerator=outlier_rare_pairs/outlier_all_pairs,
            denominator=inlier_rare_pairs/inlier_all_pairs,
            enrichment=numerator/denominator) %>% 
    # confidence intervals
    mutate( log_se = sqrt(1/outlier_rare_pairs - 1/outlier_all_pairs + 1/inlier_rare_pairs - 1/inlier_all_pairs),
            lower_CI = enrichment * exp(-1.96*log_se),
            upper_CI = enrichment * exp(1.96*log_se)) %>%
    mutate(sortcol = as.numeric(thresh)) %>%
    arrange(sortcol)
  
  return(enrichment.df)
}

graphEnrichment <- function(enrichment.df, title.mod.pop){
  enrichment.plot = ggplot(enrichment.df, aes(x = sortcol, y = enrichment)) +
    geom_point(size = 2) + 
    geom_text(aes(label = outlier_all_pairs), hjust = 0.7, vjust = -4) + 
    geom_errorbar(aes(ymax = upper_CI, ymin = lower_CI)) + 
    labs(title=paste0("Enrichment of rare variants in outlier genes (",title.mod.pop,")"), x = "Z-score threshold")
  
  return(enrichment.plot)
}

enrichmentContingencies <- function(enrichment.df, thresholds){
  ctables <- lapply(thresholds, function(x){
    cont.tmp <- enrichment.df %>% filter(sortcol == x) %>%
      transmute(outlier_rare_pairs = outlier_rare_pairs, 
                inlier_rare_pairs = inlier_rare_pairs,
                outlier_not_rare = outlier_all_pairs - outlier_rare_pairs,
                inlier_not_rare = inlier_all_pairs - inlier_rare_pairs, threshold = sortcol )
    data.frame(Outlier = c(cont.tmp$outlier_rare_pairs,cont.tmp$outlier_not_rare),
               Inlier = c(cont.tmp$inlier_rare_pairs, cont.tmp$inlier_not_rare),
               rownames=c("Rare", "Not Rare"))
    # %>% pander(caption = paste0("Threshold = ",x))
  })
  return(ctables)
}



enrichmentTipTail <- 
  function(
    pop_subset_file, 
    z_scores_file, 
    rare_var_pairs_file, 
    thresholds=seq(1,5,0.5),
    outlier_bin_width = NULL, 
    outlier_bin_percentile=NULL,
    output_dir = "./", 
    output_prefix = NULL, 
    title.mod.pop="Population", 
    col.names = NULL,
    save_thresh_files = TRUE){
    
    counts.summary <- summarizeRareInlierOutlier(
      pop_subset_file = pop_subset_file, 
      rare_var_pairs_file =  rare_var_pairs_file,
      z_scores_file = z_scores_file,
      thresholds = thresholds,
      outlier_bin_width = outlier_bin_width,
      outlier_bin_percentile = outlier_bin_percentile,
      output_dir = output_dir,
      output_prefix = output_prefix,
      col.names = col.names,
      save_thresh_files = save_thresh_files)
    
    enrichment.df <- computeEnrichment(counts.summary)
    
    gg <- graphEnrichment(enrichment.df = enrichment.df,title.mod.pop = title.mod.pop)
    
    contingencies <- enrichmentContingencies(enrichment.df,thresholds)
    
    return(list(summary=counts.summary,enrichment=enrichment.df, plot=gg, contingencies=contingencies))
  }

