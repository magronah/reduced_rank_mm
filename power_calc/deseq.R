library(tibble)
library(DESeq2)
library(Matrix)
library(huge)
library(foreach)

source("/project/6006158/agronahm/reduced_rank_mm/func2.R")
source("/project/6006158/agronahm/reduced_rank_mm/initial_param0.R")
############################################################
path = paste0("~/scratch/dataset/RR","/",nsubj,"_",ntaxa,"/deseq")
path

files   =   list.files(path, full.names = TRUE)
res = foreach(i = files, .combine = "cbind") %do% {
  res        =    readRDS(i)
  res$padj
}

dd            =    as.data.frame(res)
colnames(dd)  =    paste0("nsim", 1:ncol(dd))
rownames(dd)  =    paste0("taxon",1:ntaxa)

saveRDS(dd, file = paste0("/project/6006158/agronahm/reduced_rank_mm/",nsubj,"_",ntaxa,"/pvalues/deseq.rds"))
