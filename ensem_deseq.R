setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/new_sim/func.R")
library(tibble)
library(DESeq2)
path = paste0(getwd(),"/reproducible/new_sim/10_50/")
###########################################################
data      =   readRDS(paste0(path,"otu_meta_list_nozi.rds"))
################################################################
cc      =   commandArgs(trailingOnly  = TRUE)
i       =   as.integer(cc[1]) 
dd      =   data[[i]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data

res        =   deseqfun(countdata,met_dd,ref_name="control")
saveRDS(res,file=paste0("~/scratch/dataset/new_sim/10_50/deseq_res/mod",i,".rds")) 

