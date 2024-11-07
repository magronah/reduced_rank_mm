## goal: pick 'theta' parameters for a reduced-rank model in a sensible way
##' @param d rank (dimension)
##' @param n full dimension (latent variables per group)
##' @param logsdvec vector of log-SDs of each factor
get_theta_corrRR <- function(d, n, logsdvec) {
  mat <- matrix(0, nrow=n, ncol=d)
  if (length(logsdvec) == 1) logsdvec <- rep(logsdvec, d)
  stopifnot(length(logsdvec) == d)
  for (i in 1:d) {
    r <- rnorm(n-i+1)
    rexp <- exp(r)
    mat[(i:n), i] <- sqrt(rexp/sum(rexp))
  }
  mat <- sweep(mat, 2, FUN = "*", exp(logsdvec))
  stopifnot(all.equal(sqrt(colSums(mat^2)), exp(logsdvec)))
  theta <- c(
    mat[row(mat)==col(mat)],  ## diagonal elements
    mat[row(mat)>col(mat)]    ## below-diagonal elements
  )
  return(theta)
}


# extract otu table and metadata
norm_fun <- function(dd){
  
  otu_tab  =   long_dd = meta_data = list()
  
  nsubj =  unique(dd$subject)
    t   =  unique(dd$time)
  
  for(i in 1:length(t)){
    wide_dd <- dd %>% 
      subset(time  == t[i]) %>% 
      dplyr::select(subject, taxon, count, group) %>%
      spread(key = taxon, value = count) %>%
      setNames(c("subject","group", paste0("taxon",1:ntaxa))) 
    
    otu_table  =  wide_dd %>%
      dplyr::select(paste0("taxon",1:ntaxa))
    
    metadata   =  wide_dd %>% 
      dplyr::select(c("subject","group"))
    
             otu_count  =   t(otu_table)
    colnames(otu_count) =   paste0("subject",nsubj)
    
    dds        =   DESeqDataSetFromMatrix(otu_count,metadata, ~group)
    dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
    normalizer =   sizeFactors(dds) # one for each subject
    
    
    normalise_dd  =    data.frame(normalizer, subject = wide_dd$subject)
    dd_res        =    left_join(dd,normalise_dd, by ="subject")  
    
    long_dd[[i]]    =   dd_res
    otu_tab[[i]]    =   otu_count
    meta_data[[i]]  =   metadata
    }

  names(otu_tab)  =  names(long_dd)  =  paste0("time", t)
  list(otu_tab    =  otu_tab, long_dd= long_dd, metadata = meta_data)
}



otu_meta_lst_fun = function(res_dd1){
  
  res_dd2   = list()
  for(i in 1:length(res_dd1)){
    dd  =  res_dd1[[i]]
    
    otu_table <- dd %>%
      dplyr::select(subject, taxon, count, group)
    
    otu_table <- spread(otu_table, key = taxon, value = count) 
    colnames(otu_table) <- c("subject","group", paste0("taxon",1:ntaxa))
    
    met_data    =   otu_table %>%
      dplyr::select(subject,group)
    
    countdata     =    otu_table %>%
      dplyr::select(-c(subject,group))
    
    rownames(countdata) = (met_data$subject)
    
    res_dd2[[i]]= lst(countdata,met_data)
  }
  
  names(res_dd2)  =   paste0("sim", 1:nsim)
  res_dd2
}


meta_data = function(ntaxa, nsubj){
  metadata <- data.frame(subject = factor(seq(nsubj)),
                         group = rep(c("control", "treat"), each = nsubj/2))
  df <- expand.grid(taxon = factor(seq(ntaxa)), subject = factor(seq(nsubj)))
  df <- merge(metadata, df, by = "subject")
  df
}


get_corr <- function(ntaxa, nsubject, sparse_thresh = 0.05, seed = NULL) {
  
  set.seed(seed)
  synthetic_data   <-   huge.generator(n = nsubject, d = ntaxa, graph = "scale-free")
  covariance_matrix  <- synthetic_data$sigma
  correlation_matrix <- cov2cor(covariance_matrix)
  
  correlation_matrix[abs(correlation_matrix) < sparse_thresh] <- 0
  return(correlation_matrix)
}


get_theta_corr <- function(ntaxa,nsubject, mat= NULL, seed = NULL) {
  if(!is.null(mat)){C  <- mat}
  else{set.seed(seed); C <- get_corr(ntaxa, nsubject,seed = seed)}
  C <- nearPD(C)$mat
  scale <- sqrt(fastmatrix::ldl(as.matrix(C))$d)
  cc2 <- chol(C) %*% diag(1/scale)
  cc2[upper.tri(cc2)]
}


metadata <- function(ntaxa, nIndiv, ntime){
  d <- expand.grid(taxon = factor(1:ntaxa),
                   subject = factor(rep(1:nIndiv,ntaxa)))[1:(ntaxa*nIndiv), ]
  expdes <- data.frame(subject = factor(1:nIndiv), 
                       group=rep(c("control","treatment"), 
                                 each = nIndiv/2))
  dat0 <- left_join(d, expdes, by = "subject")
  
  dd <- do.call("rbind", replicate(n=ntime, dat0, simplify = FALSE))
  dd$time = rep(1:ntime, each=ntaxa*nIndiv)
  dd$nugget <- factor(1:nrow(dd)) 
  dd
}


get_theta_logSD <- function(n,  meanlog = 0, sdlog = 1,seed = NULL, rank = NULL) {
  set.seed(seed)
  val <- rlnorm(n, meanlog, sdlog)  # Log-normal distribution
  logSD <- log(sqrt(val))
  
  if (!is.null(rank)) {
    logSD <- logSD[1:rank]
    return(logSD)
  } else {
    return(logSD)
  }
}



# extract otu table and metadata
otu_meta_fun <- function(dd){
  
  wide_dd <- dd %>%
    dplyr::select(subject, taxon, count, group) %>%
    spread(key = taxon, value = count) %>% 
    setNames(c("subject","group", paste0("taxon",1:ntaxa)))    
#  colnames(wide_dd) = c("subject","group", paste0("taxon",1:ntaxa))
  
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

 
