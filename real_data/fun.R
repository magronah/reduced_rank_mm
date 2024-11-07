# filter_fun <- function(countdata, metadata,abund_thresh=10, sample_thresh=5){
#   # if( nrow(metadata) != ncol(countdata)){
#   #   countdata =  (countdata)
#   # }
#   keep <- rowSums(countdata >= abund_thresh) >= sample_thresh
#   countdata = countdata[keep,]
#   as.data.frame(t(countdata))
# }

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


df_long = function(dd, otu_names = "sp", subject_name = "subject"){
  if(nrow(dd) > ncol(dd)){
    dd  =  as.data.frame(t(dd))
  }
  df   =  dd %>%
    rownames_to_column(subject_name)
  ddd =   pivot_longer(df,
                      cols = starts_with(otu_names),
                      names_to  = otu_names,
                      values_to = "count")
  names(ddd)  = c("subject","taxon","count")
  ddd
}



# extract otu table and metadata
otu_meta_fun <- function(dd){
  
  wide_dd <- dd %>%
    dplyr::select(subject, taxon, count, group) %>%
    spread(key = taxon, value = count) %>%
    setnames(c("subject","group", paste0("taxon",1:ntaxa))) 
  
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
                     shrinkage_method="normal"){
  
  #check otu table is in otu by samples format
  if(all((met_data$subject)==colnames(countdata)) == FALSE){
    countdata = t(countdata)
  }
  
  #remove samples with zeros for all taxa (if any such sample exist)
  keep <- (colSums(countdata) > 0)
  countdata = countdata[,keep]
  met_data= met_data[keep, ]
  
  # call deseq
  dds <- DESeqDataSetFromMatrix(countdata,met_data, ~group)
  dds$group <- relevel(dds$group, ref = ref_name)
  
  dds <- DESeq(dds,sfType ="poscounts",
               minReplicatesForReplace = minReplicatesForReplace) 
  
  res <- results(dds, cooksCutoff=cooksCutoff, 
                 independentFiltering=independentFiltering,
                 alpha = alpha_level)
  
  reslt <- lfcShrink(dds, res=res, coef=2, type=shrinkage_method)
  
  deseq_est = data.frame(reslt)
  deseq_est$dispersion = dispersions(dds)
  deseq_dd   =  deseq_est  %>% 
               rownames_to_column(var = "param_name")
  
  list(result = deseq_dd, object = dds)
  
}
