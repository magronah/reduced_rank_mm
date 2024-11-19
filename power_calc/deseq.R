library(tibble)
library(DESeq2)
library(Matrix)
library(huge)
library(foreach)
source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/dataset/RR","/",nsubj,"_",ntaxa,"/deseq")
path

files   =   list.files(path, full.names = TRUE)
res = foreach(i = files, .combine = "cbind") %do% {
  res        =    readRDS(i)
  res$padj
}

res            =    as.data.frame(res)
colnames(res)  =    paste0("nsim", 1:ncol(res))
rownames(res)  =    paste0("taxon",1:ntaxa)

saveRDS(dd, file = paste0(getwd(),"/",nsubj,"_",ntaxa,"pvalues/deseq.rds"))
