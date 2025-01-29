# filter_fun <- function(countdata, metadata,abund_thresh=10, sample_thresh=5){
#   # if( nrow(metadata) != ncol(countdata)){
#   #   countdata =  (countdata)
#   # }
#   keep <- rowSums(countdata >= abund_thresh) >= sample_thresh
#   countdata = countdata[keep,]
#   as.data.frame(t(countdata))
# }

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
