library(tibble)
library(DESeq2)
library(Matrix)
library(huge)
library(foreach)
library(scam)
source("/project/6006158/agronahm/reduced_rank_mm/func2.R")
source("/project/6006158/agronahm/reduced_rank_mm/initial_param0.R")
############################################################
path  =  paste0("/project/6006158/agronahm/reduced_rank_mm/",nsubj,"_",ntaxa)

pvalue      =   readRDS(paste0(path,"/pvalues/deseq.rds"))
pvalue      =   rowMeans(pvalue)
####################################################################
data        =   readRDS(paste0(path,"/otu_meta_list_withzi_taxa.rds"))
effect      =   readRDS(paste0(path,"/true_param.rds"))
effect_size =   effect$true_param
mean_count  =   colMeans(do.call(rbind,lapply(data, function(x) (x$countdata))))
############################################################
mod    =  gam_fit(pvalue, effect_size, mean_count,grid_len = 500,alpha_level = 0.05)

file_path  =  paste0(path,"/GAM/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}



saveRDS(mod, file=paste0(file_path,"deseq.rds"))

#mm
