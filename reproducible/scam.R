library(tibble)
library(dplyr)
library(scam)
library(tidyverse)
library(here)
source("reproducible/func.R")
############################################
nsubj = 150; ntaxa = 500 
path <- paste0(getwd(), "/", nsubj, "_", ntaxa, "/")
##################################################
## Read all the needed data
true_param =  readRDS(paste0(path,"true_param.rds"))

# Read mean count data
mean_count <- readRDS(paste0(path, "mean_count.rds"))

## Read estimates from fitting ensembles
dd   = list(
  rr     =  readRDS(paste0(path, "rr.rds")),
  rrzi   =  readRDS(paste0(path, "rrzi.rds")),
  us     =  readRDS(paste0(path, "us.rds")),
  uszi   =  readRDS(paste0(path, "uszi.rds"))
)
##########################################################
effect_size <-  true_param$true_param

# Compute p-values for all models
pvals <- lapply(dd, pvalue_cal)
names(pvals) <- c("rr", "rrzi", "us", "uszi")
##########################################################  
# Fit GAM models for each p-value set
mod_rr    <- gam_fit(pvals[["rr"]], effect_size, 
                     mean_count, alpha_level = 0.05)  # no error

mod_rrzi  <- gam_fit(pvals[["rrzi"]], effect_size, 
                     mean_count, alpha_level = 0.05)  # no error

mod_us    <- gam_fit(pvals[["us"]], effect_size, 
                     mean_count, alpha_level = 0.05)   # no error

mod_uszi  <- gam_fit(pvals[["uszi"]], effect_size, 
                     mean_count, alpha_level = 0.05)

##   Error in eigen(S.t, symmetric = TRUE) : infinite or missing values in 'x'
##   Called from: eigen(S.t, symmetric = TRUE)
#######################################################
any(is.na(pvals[["uszi"]]))
any(is.infinite(pvals[["uszi"]]))
#######################################################
#For nsubj = 150; ntaxa = 500, I get the same error for all the model 
## except for mod_us which works fine. 
