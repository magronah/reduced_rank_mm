load_models <- function(path, filenames) {
  file_paths <- paste0(path, filenames)
  mod_list  <- lapply(file_paths, readRDS)
}



my_aicc_fun  <- function(mod){
  
  loglik      =   as.numeric(logLik(mod))
  num_params  =   attr(logLik(mod), "df")  
  correction  =   2*num_params*(num_params+1)/(nobs(mod) - num_params -  1)
  
  if(is.na(loglik)){
    loglik      =    as.numeric(mod$obj$fn())
    aic         =    2*(loglik) +  2*(num_params)
    aicc        =    aic   +  correction
    return(aicc)
  }else{
    aic         =   -2*(loglik) +  2*(num_params)
    correction  =   2*num_params*(num_params+1)/(nobs(mod) - num_params -  1)
    aicc        =   aic   +  correction
    return(aicc)
  }
}
# filter_fun <- function(countdata, metadata,abund_thresh=10, sample_thresh=5){
#   # if( nrow(metadata) != ncol(countdata)){
#   #   countdata =  (countdata)
#   # }
#   keep <- rowSums(countdata >= abund_thresh) >= sample_thresh
#   countdata = countdata[keep,]
#   as.data.frame(t(countdata))
# }

#'
#' @param mod 
#' @param ntaxa 
#' @param conf_level 
#'
#' @return
#' @export
#' 
#' @examples
wald_confint = function(mod, conf_level = .95, in_sd = 1){
  pred    =   predict(mod, type = "latent", se.fit = TRUE) 
  est     =   pred$fit
  sd_err  =   pred$se.fit
  
  # Calculate z-score for the desired confidence level
  z_score    =    qnorm(conf_level + (1 - conf_level)/2)
  
  # Calculate confidence intervals
  dd_full    =    data.frame(est_param =   est,
                             lwr       =   est - z_score*in_sd*sd_err,
                             upr       =   est + z_score*in_sd*sd_err)
  
  # Calculate p-values
  z_stat = est / sd_err
  p_values = 2 * (1 - pnorm(abs(z_stat)))  # Two-tailed p-value
  
  dd_full$pvalue = p_values
  
  dd_full
}

##############################################
extract_grp_effect  =  function(dd_full,ntaxa){
  indx       =    seq(2, (2*ntaxa), 2)
  dd         =    dd_full[indx, ]
  dd$width    =    dd$upr  -  dd$lwr 
  # Add parameter names
  dd$param_name =   factor(paste0("taxon",1:ntaxa),
                           levels = paste0("taxon",1:ntaxa))
  dd
}

deseq_wald_confint = function(deseq_est, true_param_dd,conf_level = .95,model = "deseq_wald"){
  foldchange =  2^deseq_est$log2FoldChange
  dd         =  data.frame(foldchange = foldchange)
  z_score    =    qnorm(conf_level + (1 - conf_level)/2)
  sd_err     =   deseq_est$lfcSE
  dd$lwr     =   foldchange  -  z_score*sd_err
  dd$upr     =   foldchange  +  z_score*sd_err
  dd$model   =  rep(model,nrow(dd)) 
  dd$param_name  =  deseq_est$param_name
  dd$true_param  =   true_param_dd$true_param
  dd
}


custom_theme <- function(n) {
  theme_bw(base_size = n) +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = n, family = "Roboto"),
      axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      axis.text.y = element_text(family = "Roboto", size = n, color = "black")
    )
}


filter_fun <- function(countdata, metadata,abund_thresh=5, sample_thresh=3){
  
  ## sanity check
  if(all((metadata$subject)==colnames(countdata)) == FALSE){
    countdata = t(countdata)
  }
  
  ##########################################
  # filter
  dds <- DESeqDataSetFromMatrix(countdata,metadata, ~group)
  keep <- rowSums(counts(dds) >= abund_thresh) >= sample_thresh
  
  dds=dds[keep,]
  data.frame(counts(dds))
}

df_long = function(dd, otu_names = "sp", subject_name = "subject", ntaxa){
  if(ncol(dd) != ntaxa){
    dd  =  as.data.frame(t(dd))
  }else{
    dd  =  data.frame(dd)
  }
  df   =  dd %>%
    rownames_to_column(subject_name)
  ddd =   pivot_longer(df,
                       cols = starts_with(otu_names),
                       names_to  = otu_names,
                       values_to = "count")
  names(ddd)  = c(subject_name,"taxon","count")
  ddd
} 



# extract otu table and metadata
otu_meta_fun <- function(dd, ntaxa){
  wide_dd <- dd %>%
    dplyr::select(subject, taxon, count, group) %>%
    spread(key = taxon, value = count) %>%
    setNames(c("subject","group", paste0("taxon",1:ntaxa))) 
  
  otu_table  =  wide_dd %>%
    dplyr::select(paste0("taxon",1:ntaxa))
  
  metadata   =  wide_dd %>% 
    dplyr::select(c("subject","group"))
  
  otu_count  =   t(otu_table)
  dds        =   DESeqDataSetFromMatrix(otu_count,metadata, ~group)
  dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
  normalizer =   sizeFactors(dds) # one for each subject
  
  
  normalise_dd  =    data.frame(normalizer, subject = wide_dd$subject)
  dd            =    left_join(dd,normalise_dd, by ="subject")  
  dd
  
  #list(metadata = metadata, otu_table = otu_table,wide_dd=wide_dd)
}


deseqfun <- function(countdata,met_data,alpha_level=0.1,ref_name="NT",
                     minReplicatesForReplace = Inf, 
                     cooksCutoff = FALSE,
                     independentFiltering = FALSE, 
                     do_shrinkage =  "yes", 
                     shrinkage_method="normal",
                     design = ~group, 
                     ntaxa){
  

  #check otu table is in otu by samples format
  if(nrow(countdata) != ntaxa){
    countdata = t(countdata)
  }
  
  #remove samples with zeros for all taxa (if any such sample exist)
  keep <- (colSums(countdata) > 0)
  countdata = countdata[,keep]
  met_data  = met_data[keep, ]
  
  # call deseq
  dds <- DESeqDataSetFromMatrix(countdata,met_data, design = design)
  dds$group <- relevel(dds$group, ref = ref_name)
  
  dds <- DESeq(dds,sfType ="poscounts",
               minReplicatesForReplace = minReplicatesForReplace) 
  
  res <- results(dds, cooksCutoff=cooksCutoff, 
                 independentFiltering=independentFiltering,
                 alpha = alpha_level)
  
  if(do_shrinkage == "no"){
    reslt   <-   res
  }else{
    reslt <- lfcShrink(dds, res=res, coef=2, type=shrinkage_method)
  }
  
  deseq_est = data.frame(reslt)
  deseq_est$dispersion = dispersions(dds)
  deseq_est$intercept  = coef(dds)[, "Intercept"]
  
  deseq_dd   =  deseq_est  %>% 
               rownames_to_column(var = "param_name")
  
  list(result = deseq_dd, 
       data  =  list(countdata =  countdata,
                     meta_data = met_data),
       object = dds)
  
}
