library(tibble)
library(DESeq2)
library(Matrix)
library(huge)
library(foreach)
source("/project/6006158/agronahm/reduced_rank_mm/func2.R")
source("/project/6006158/agronahm/reduced_rank_mm/initial_param0.R")
############################################################
path  =  paste0("/project/6006158/agronahm/reduced_rank_mm/",nsubj,"_",ntaxa)

pval  =  readRDS(paste0(path,"/pvalues/deseq.rds"))
####################################################################
data        =   readRDS(paste0(path,"/otu_meta_list_withzi_taxa.rds"))
effect      =   readRDS(paste0(path,"/true_param.rds"))
effect      =   effect$true_param
mean_count  =  colMeans(do.call(rbind,lapply(data, function(x) (x$countdata))))
############################################################

