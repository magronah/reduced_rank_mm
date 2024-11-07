setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)
source("load_glmmTMB.R")

#fl    =    paste0(ntaxa,"_",nsubj)
path  <-   paste0("~/scratch/dataset/new_sim/100_300/rr_mod/")
files <-   list.files(path, full.names = TRUE)

res = foreach(i = files, .combine ="cbind") %do% {
    mod  =   readRDS(i)
    dd   =  (ranef(mod)$cond$taxon$grouptreat)
    dd 
}

dd            =  as.data.frame(res)
rownames(dd)  =  paste0("taxon",1:nrow(dd))
colnames(dd)  =  paste0("sim",1:ncol(dd))

saveRDS(dd, file = paste0(getwd(),"/reproducible/new_sim/100_300/parametric_rr.rds"))


