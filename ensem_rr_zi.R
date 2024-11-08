library(RhpcBLASctl)
library(Matrix)
library(huge)
library(glmmTMB)
library(tidyverse)
library(dplyr)
library(DESeq2)
source("func2.R")
source("initial_param0.R")
path = paste0(getwd(),"/100_300/")
###########################################################
data	  =   readRDS(paste0(path,"sim_count_list_withzi_taxa.rds"))
form	  =   count ~ 1 + us(1 + group|taxon) +  rr(0 + taxon | subject,2)
################################################################
cc =   commandArgs(trailingOnly  = TRUE)
i	  =   as.integer(cc[1])

##Fit the model
par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)


 dd	=   data[[i]]
  ################################################################
  ##now add normalization constant 
  df  =  otu_meta_fun(dd)

  gprior  <- data.frame(prior = "gamma(2, 2.5)",
                        class = "theta_sd",
                        coef = "")
  
  form2  =   update(form,.~. + offset(normalizer))
  
  options(glmmTMB_openmp_debug = TRUE)
  blas_set_num_threads(1)
  
  system.time(
    fit <- glmmTMB(form2, data = df,
                   family = nbinom2,
                   ziformula = ~1 + (1|taxon),
                   control = par_ctrl  # Ensure par_ctrl is correctly configured for glmmTMB
    )
  )
  
  system.time(
    fit  <-  glmmTMB(form2, data = df,
                     family  = nbinom2,
                     ziformula  = ~1 + (1|taxon),
                     #prior   = gprior,
                     REML    = TRUE,
                     control = par_ctrl
    )
  )
  
#saveRDS(fit, file=paste0(path, "mod_fit/mod",i,".rds"))


saveRDS(fit, file=paste0("~/scratch/dataset/RR/100_300/rr/mod",i,".rds"))
# true =  readRDS(paste0(path,"true_param.rds"))

# res  =  ranef(fit, condVar=FALSE)
# est  = (res$cond$taxon$grouptreat)
# 
# plot(est,true$true_param)
# abline(0,1)
# nsim
                                                    
