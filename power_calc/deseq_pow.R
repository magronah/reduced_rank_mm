library(tibble)
library(DESeq2)
library(Matrix)
library(huge)
library(foreach)
source("func2.R")
source("initial_param0.R")
############################################################
path  =  "/project/6006158/agronahm/reduced_rank_mm/"
pval  =  readRDS(dd, file = paste0(path,nsubj,"_",ntaxa,"/pvalues/deseq.rds"))
####################################################################
dd        =   readRDS(paste0(path,nsubj,"_",ntaxa,"otu_meta_list_withzi_taxa.rds"))
effect    =   readRDS(paste0(path,nsubj,"_",ntaxa,"true_param.rds"))
mean_count  =  rowMeans(dd)
############################################################



files   =   list.files(path, full.names = TRUE)
res = foreach(i = files, .combine = "cbind") %do% {
  res        =    readRDS(i)
  res$padj
}

res            =    as.data.frame(res)
colnames(res)  =    paste0("nsim", 1:ncol(res))
rownames(res)  =    paste0("taxon",1:ntaxa)

saveRDS(dd, file = paste0(getwd(),"/",nsubj,"_",ntaxa,"pvalues/deseq.rds"))
