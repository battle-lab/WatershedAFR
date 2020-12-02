library(tidyr)
library(dplyr)
library(ggplot2)


summarizeRareInlierOutlier <- function(pop_subset_file, z_scores_file, rare_var_pairs_file, thresholds=seq(1,5,0.5), output_prefix = "default"){
  names(thresholds) <- thresholds
  
  pop_list <- as.character(read.table(pop_list_file, header = FALSE)$V1)
  
  rare_var_pairs <- read.table(rare_var_pairs_file, header=TRUE, check.names = FALSE) %>%
    mutate(easykey = paste0(Gene,Ind))
  # remove "duplicate" gene-individual pairs with rare variants NOTE: Only first occurance in table is kept
  rare_var_pairs <- rare_var_pairs[!duplicated(rare_var_pairs$easykey),] 
  
  # create a subtable that will be more useful ## DO NOT CHANGE ITS ORDER AFTER THIS STEP
  pairs.df.pop <- read.table(z_scores_file, header=TRUE, check.names = FALSE) %>%
    mutate(easykey = paste0(Gene,Ind)) %>% # create an id column for gene-indiv pairs
    select(easykey,Gene,Ind,MedZ) %>%
    filter(Ind %in% pop_list) %>% #limit pairs.df.pop to individuals within the population
    mutate(has_rare=ifelse(easykey %in% rare_var_pairs$easykey,1,0))  %>% # add column to pairs.df.pop indicating if the given gene-individual pair has at least 1 rare variant within 10kb of gene TSS
    ungroup()
  
  
  # determine inlier/outlier status based on Zscore thresholds
  # for (thresh in zscore_thresh_list) {
  #   pairs.df.pop[[paste0("outlier_status_at_thresh_",thresh)]] = ifelse(abs(pairs.df.pop$MedZ) > thresh, "outlier","inlier")
  # }
  # determine inlier/outlier status based on Zscore thresholds -- using apply for future parallelization cababilities jic
  # outlier_status <- as.data.frame(sapply(thresholds, function(x){
  #   tmp <- ifelse(abs(pairs.df.pop$MedZ) > as.numeric(x), "outlier","inlier")
  #   return(tmp)
  # }))
  # colnames(outlier_status) <- paste0("outlier_status_at_thresh_",thresholds)
  
  # all.pairs.thresh is a list titled by threshold that contains all inlier/outlier information for each pair at that threshold.
  all.pairs.thresh <- lapply(thresholds,function(t){
    outlier_stat <- ifelse(abs(pairs.df.pop$MedZ) > as.numeric(t), "outlier","inlier")
    #threshcol <- paste0("outlier_status_at_thresh_",t)
    tmp <- pairs.df.pop
    tmp$outlier <- outlier_stat
    tmp <- tmp %>%
      group_by(Gene) %>% 
      mutate(OutlierGenePairs = sum(outlier == "outlier",na.rm = T), GeneKeep = OutlierGenePairs > 0) %>%
      ungroup() 
    
    return(tmp[,c("outlier","GeneKeep")])
  })
  
  # # collect the indices of the rows pertinent to each threshold (meaning that all genes have at least one outlier at the given threshold)
  # outlier.pairs.thresh <- lapply(all.pairs.thresh, function(x){
  #   return(which(x$GeneKeep))
  #   
  # })
  
  
  # frequency table to count inliers and outliers with or without rare variants
  counts.list <- lapply(thresholds, function(t){
    gktmp <- 
      cbind(pairs.df.pop, all.pairs.thresh[[t]]) %>% 
      filter(GeneKeep)
    
    tmp1 <- gktmp %>%
      count(outlier)
    
    # if (all(c("outlier","inlier") %in% unique(tmp1$outlier))){
    tmp1 <- tmp1 %>%
      pivot_wider(names_from=outlier,values_from=n) %>%
      rename(outlier_all_pairs=outlier, inlier_all_pairs=inlier)
    # } else {
    # tmp1 <- tmp1 %>%
    # pivot_wider(names_from=outlier,values_from=n,names_prefix="ERROR")
    # }
    
    tmp2 <- gktmp %>% 
      filter(has_rare == 1) %>%
      count(outlier) 
    print("Here are the unique values in the outlier column after filtering for rare")
    print(unique(tmp2$outlier))
    print("Here's the head of the table:")
    print(head(tmp2))
    
    # if (all(c("outlier","inlier") %in% unique(tmp2$outlier))){
    tmp2 <- tmp2 %>%
      pivot_wider(names_from=outlier,values_from=n) %>%
      rename(outlier_rare_pairs=outlier, inlier_rare_pairs=inlier)
    # } else{
    # tmp2 <- tmp2 %>%
    # pivot_wider(names_from=outlier,values_from=n,names_prefix="ERROR")
    # }
    tmp <- cbind(tmp1,tmp2)
    tmp$thresh <- t
    
    return(tmp)
  })
  # frequency table to count inliers and outliers with or without rare variants
  summary.df <- bind_rows(counts.list)
  
  return(summary.df)
}


computeEnrichment <- function(summary.df){
  
  # compute enrichment
  enrichment.df  <- summary.df %>% mutate(numerator=outlier_rare_pairs/outlier_all_pairs, denominator=inlier_rare_pairs/inlier_all_pairs, enrichment=numerator/denominator) %>% 
    # confidence intervals
    mutate(log_se = sqrt(1/outlier_rare_pairs - 1/outlier_all_pairs + 1/inlier_rare_pairs - 1/inlier_all_pairs), lower_CI = enrichment * exp(-1.96*log_se), upper_CI = enrichment * exp(1.96*log_se)) %>%
    mutate(sortcol = as.numeric(thresh)) %>%
    arrange(sortcol)
  
  return(enrichment.df)
}

computeEnrichment <- function(summary.df){
  
  # compute enrichment
  enrichment.df  <- summary.df %>% mutate(numerator=outlier_rare_pairs/outlier_all_pairs, denominator=inlier_rare_pairs/inlier_all_pairs, enrichment=numerator/denominator) %>% 
    # confidence intervals
    mutate(log_se = sqrt(1/outlier_rare_pairs - 1/outlier_all_pairs + 1/inlier_rare_pairs - 1/inlier_all_pairs), lower_CI = enrichment * exp(-1.96*log_se), upper_CI = enrichment * exp(1.96*log_se)) %>%
    mutate(sortcol = as.numeric(thresh)) %>%
    arrange(sortcol)
  
  return(enrichment.df)
}

graphEnrichment <- function(enrichment.df, title.mod.pop){
  enrichment.plot = ggplot(enrichment.df, aes(x = sortcol, y = enrichment)) +
    geom_point(size = 2) + 
    geom_text(aes(label = outlier_all_pairs), hjust = 0.7, vjust = -4) + 
    geom_errorbar(aes(ymax = upper_CI, ymin = lower_CI)) + labs(title=paste0("Enrichment of rare variants in outlier genes (",title.mod.pop,")"), x = "Z-score threshold")
  
  return(enrichment.plot)
}


enrichmentTipTail <- 
  function(
    pop_subset_file, z_scores_file, 
    rare_var_pairs_file, thresholds=seq(1,5,0.5), 
    path.prefix = "default", title.mod.pop="Population"){
    
    counts.summary <- summarizeRareInlierOutlier(
      pop_subset_file = pop_subset_file, 
      rare_var_pairs_file =  rare_var_pairs_file,
      z_scores_file = z_scores_file, thresholds = thresholds)
    
    enrichment.df <- computeEnrichment(counts.summary)
    
    gg <- graphEnrichment(enrichment.df = enrichment.df,title.mod.pop = title.mod.pop)
    
    print(gg)
  }

