setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(glmmTMB)
library(ggplot2)
library(dplyr)
source("reproducible/leverage_funs.R")
library(Matrix)
##############################################################
path1   =   paste0("real_data/")
source(paste0(path1,"fun.R"))
#########################################################
# Define filenames to load
filenames <- c("rr_mod.rds")
######################################################################
# Load models for autism data
autism_path <- paste0("real_data/autism_data/results/")
soil_path   <- paste0("real_data/soil_data/results/")

# Load models for autism data
autism_models <- load_models(autism_path, filenames)
soil_models   <- load_models(soil_path, filenames)

# Assigning names 
mod_names <-  "RR"
names(autism_models)    =   names(soil_models)    =    mod_names
mod_list    =   lst(autism_models,soil_models)

rr_list <- lapply(mod_list, function(x){x$RR})
nam   =  c("a
#########################################################
par_ctrl <- glmmTMBControl(parallel = list(n = 10, autopar = TRUE),
                           optCtrl  = list(eval.max=1000, iter.max = 100))

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")
############################################################
cc     =   commandArgs(trailingOnly  = TRUE)
i      =   as.integer(cc[1])
############################################################
mmd  <-  mod
df   <-  model.frame(mod) 
dd   <-  df |> dplyr::rename(normalizer = "offset(normalizer)")

lev  <-  leverage_brute_modified(mmd, data = dd, inds = i, eps = 0.1,
                        fit_method = "update", scale = "response",
                        pred_method  = "predict",
                        progress = TRUE) 
############################################################
file_path  =  paste0(path,"leverage/rr/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(lev, file= paste0(file_path,"lev",i, ".rds"))
