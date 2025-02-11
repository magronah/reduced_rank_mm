library(huge)
library(glmmTMB)
library(Matrix)
library(here)
source("func2.R")
source("initial_param0.R")
###########################################################
path = paste0("~/scratch/data/",nsubj,"_",ntaxa,"/uszi")
path

files <-   list.files(path, full.names = TRUE)
############################################################
path_dir  =   paste0(nsubj,"_",ntaxa,"/")

cc        =   commandArgs(trailingOnly  = TRUE)
i         =   as.integer(cc[1])

for(i  in 1:10){
###########################################################
model     =   readRDS(files[i])
mu_count  =   readRDS(paste0(path_dir,"mean_count.rds"))
nam       =   paste0("uszi",i)
###########################################################
save_path  =  paste0(path_dir,"sdreport/")

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
  cat("Folder created at:", save_path, "\n")
} else {
  cat("Folder already exists at:", save_path, "\n")
}

confint   <-  wald_confint(mod   =  model,
                           ntaxa =  ntaxa,
                      mean_count =  mu_count,
                        mod_name =  nam,
                            path =  save_path)


confint_path  =  paste0(path_dir,"confint/")

if (!dir.exists(confint_path)) {
  dir.create(confint_path, recursive = TRUE)
  cat("Folder created at:", confint_path, "\n")
} else {
  cat("Folder already exists at:", confint_path, "\n")
}

saveRDS(confint, file = paste0(confint_path, "uszi",i,".rds"))
}
