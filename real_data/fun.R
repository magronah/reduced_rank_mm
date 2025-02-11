## https://github.com/glmmTMB/glmmTMB/issues/1059#issuecomment-2187660510
library(nloptr)
library(glmmTMB)

nl_algs <- (nloptr.get.default.options()[1,"possible_values"]
            |> strsplit(",")
            |> unlist()
            |> trimws()
            |> grep(pattern="NLOPT_LD", value = TRUE)
            |> grep(pattern="AUGLAG", value = TRUE, invert = TRUE)
)

nloptwrap1 <- function(start, objective, gradient, control,
                       algorithm = "NLOPT_LD_LBFGS") {
  ## what are nlminb default tolerances?
  my_obj <- function(x) objective(x) ##
  my_grad <- function(x) gradient(x) ##
  fit <- nloptr(x0 = start, eval_f = my_obj, eval_grad_f = my_grad,
                opts = list(algorithm = algorithm, xtol_rel = 1e-8))
  return(fit)
}

optim_algs <- c("BFGS", "CG", "L-BFGS-B")
allFit.methods <- data.frame(
  optimizer = rep(c("nlminb", "optim", "nloptwrap1"),
                  times = c(1, length(optim_algs), length(nl_algs))),
  alg = c("default", optim_algs, nl_algs))

## FIXME: don't write over existing control ('smart_refit' branch has better
##        storage of control info)
## FIXME: don't refit with original alg?
## FIXME: allow for setting tolerances
## FIXME: factory etc. to store warnings
## FIXME: summary/tidy methods?
## FIXME: exclude terrible/sensitive optimizers
## FIXME: expand to use optimx methods?
## FIXME: we should only re-run the optimization, not update()! -> finalizeTMB()
##        if necessary do the modular steps again, ideally re-use info from fitted model
##        optional whether we should start from best fit or from default/specified starting values
##        optional whether to include original alg
allFit <- function(fit, methods = allFit.methods) {
  res <- vector("list", length = nrow(methods))
  names(res) <- paste(methods$optimizer,
                      methods$alg,
                      sep = ".")
  for (i in seq(nrow(methods))) {
    optimizer <- methods$optimizer[i]
    alg <- methods$alg[i]
    cat("** ", i, optimizer, alg, "\n")
    ctrl <- switch(optimizer,
                   nlminb = glmmTMBControl(optimizer=nlminb),
                   optim = glmmTMBControl(optimizer=optim,
                                          optArgs = list(method = alg)),
                   nloptwrap1 = glmmTMBControl(optimizer=nloptwrap1,
                                               optArgs = list(algorithm = alg)))
    time <- system.time(curfit <- try(update(fit, control = ctrl)))
    attr(curfit, "system.time") <- time
    res[[i]] <- curfit
  }
  res
}

# x1 <- update(mod, control = glmmTMBControl(optimizer = "nloptwrap1"))
# 
# modrr <- update(mod, control = glmmTMBControl(optimizer = nloptwrap1))
# x1 <- glmmTMB(formula(mod), data = getCall(mod)$data, 
#               family = getCall(mod)$family, 
#               control = glmmTMBControl(optimizer = "nloptwrap1"))

if (FALSE) {
  ## minimal example for testing allFit
  x0 <- glmmTMB(mpg ~ cyl, data = mtcars)
  x1 <- update(x0, control = glmmTMBControl(optimizer = nloptwrap1))
  x2 <- update(x0, control = glmmTMBControl(optimizer = nloptwrap2),
               optArgs = list(algorithm = nl_algs[1]))
  
  ## not identical (don't want them to be)
  fixef(x1)$cond-fixef(x2)$cond
}

load_models <- function(path, filenames) {
  file_paths <- paste0(path, filenames)
  mod_list  <- lapply(file_paths, readRDS)
}

filter_complete_separation = function(dd){
  df  = dd  %>% filter(abs(est_param) > 10)
  taxa_exclude =  df$param_name
  dd   %>%  filter(!param_name %in% taxa_exclude)
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
wald_confint = function(mod, conf_level = .95, 
                        in_sd = 1, ntaxa,
                        mean_count, mod_name,
                        path){
  
 grp_ind   =   seq(2, 2*(ntaxa), 2)
  #pred    =   predict(mod, type = "latent", se.fit = TRUE) 
  #est     =   pred$fit[grp_ind]
  #sd_err  =   pred$se.fit[grp_ind]
  ########################################################
  ss <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)
  saveRDS(ss, file = paste0(path,"sdreport_",mod_name,".rds"))
  #ss  <-  readRDS(paste0(path,"sdreport_",mod_name,".rds"))
  #########################################################
  ss$jointPrecision <-  as(ss$jointPrecision, "sparseMatrix")
  inverse_mat       <-  solve(ss$jointPrecision)
    
  if(any(diag(inverse_mat) < 0)){
    precision_nearPD <-  as.matrix(nearPD(ss$jointPrecision)$mat)
         precision   <-  as(precision_nearPD, "sparseMatrix")
    se_vec <- sqrt(diag(solve(precision)))
  }else{
    se_vec <- sqrt(diag(inverse_mat))
  }
  # Run nearPD
  #se_vec <- sqrt(diag(solve(nearPD(ss$jointPrecision))))
  # allFit is another option
  # https://github.com/glmmTMB/glmmTMB/blob/master/misc/allFit.R
  #se_vec <- sqrt(diag(solve(ss$jointPrecision)))
  #start_method in glmmTMBControl jitter.sd
  
  saveRDS(se_vec, file = paste0(path,"full_sdr_",mod_name,".rds"))
  # ########################################################
  full_est     =   mod$fit$parfull
     est       =   full_est[names(full_est) == "b"][grp_ind]
     sd_err    =   se_vec[grp_ind]
  
  # Calculate z-score for the desired confidence level
  z_score    =    qnorm(conf_level + (1 - conf_level)/2)
  
  # Calculate confidence intervals
     dd      =    data.frame(est_param =   est,
                             lwr       =   est - z_score*in_sd*sd_err,
                             upr       =   est + z_score*in_sd*sd_err)
  
  # Calculate p-values
  z_stat   =  est / sd_err
  p_values =  2 * (1 - pnorm(abs(z_stat)))  # Two-tailed p-value
  
  dd$pvalue  =  p_values
  dd$width   =   dd$upr  - dd$lwr
  dd$mean_count  =  as.numeric(mean_count)
  dd$param_name  =  names(mean_count)
  dd
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

deseq_wald_confint = function(deseq_est,mean_count, conf_level = .95){
  foldchange =   log(2)*deseq_est$log2FoldChange
  dd         =   data.frame(est_param = foldchange)
  z_score    =   qnorm(conf_level + (1 - conf_level)/2)
  sd_err     =   log(2)*deseq_est$lfcSE
  dd$lwr     =   foldchange  -  z_score*sd_err
  dd$upr     =   foldchange  +  z_score*sd_err
  dd$width   =   dd$upr  - dd$lwr
  dd$param_name  =  deseq_est$param_name
  dd$pvalue      =   deseq_est$padj  
  dd$mean_count  =  as.numeric(mean_count)
  dd$param_name  =  names(mean_count)
  dd
}

#mod =  autism_models$Zinbmm

zinbmm_confint = function(mod, mean_count,
                                group_label = "groupNT",
                                conf_level = .95){
  
  sd_err =  as.numeric(unlist(lapply(mod$fit,function(x) 
                       {summary(x)$tTable[group_label, "Std.Error"]})))
  
  est    =   as.numeric(unlist(lapply(mod$fit,function(x) 
                       {summary(x)$tTable[group_label, "Value"]})))
  
  pvalue =  as.numeric(unlist(lapply(mod$fit,function(x) 
                       {summary(x)$tTable[group_label, "p-value"]})))
  
  dd         =   data.frame(est_param = est)
  z_score    =   qnorm(conf_level + (1 - conf_level)/2)
  
  dd$lwr     =   est  -  z_score*sd_err
  dd$upr     =   est  +  z_score*sd_err
  dd$width   =   dd$upr  - dd$lwr
  
  dd$pvalue      =   pvalue  
  dd$mean_count  =  as.numeric(mean_count)
  dd$param_name  =  names(mean_count)
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
