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
data      =   readRDS(paste0(path,"sim_count_list_withzi_taxa.rds"))
form      =   count ~ 1 + us(1 + group|taxon)  
################################################################
cc      =   commandArgs(trailingOnly  = TRUE)
i       =   as.integer(cc[1]) 
i =1 
dd      =   data[[i]]
################################################################
##now add normalization constant 
df  =  otu_meta_fun(dd)

##Fit the model
par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

form2  =   update(form,.~. + offset(normalizer))

system.time(
  fit  <-  glmmTMB(form2, data = df,
                   family  = nbinom2, 
                   ziformula  = ~ 1 + (1|taxon),
                   prior   = gprior,
                   REML    = TRUE,
                   control = par_ctrl
  )
)

saveRDS(fit, file=paste0("~/scratch/dataset/new_sim/100_300/us_mod/mod",i, ".rds"))

View(df)
