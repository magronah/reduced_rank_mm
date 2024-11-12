library(foreach)
library(glmmTMB)
source("func2.R")
source("initial_param0.R")
path = paste0("~/scratch/dataset/RR","/",nsubj,"_",ntaxa,"/")
path
##########################################################
files <-   list.files(path, full.names = TRUE)

res = foreach(i = files, .combine ="cbind") %do% {
    mod  =   readRDS(i, )
    dd   =  (ranef(mod)$cond$taxon$grouptreat)
    dd
}

dd            =  as.data.frame(res)
rownames(dd)  =  paste0("taxon",1:nrow(dd))
colnames(dd)  =  paste0("sim",1:ncol(dd))

saveRDS(dd, file = paste0(getwd(),"/",nsubj,"_",ntaxa,"/rrzi.rds"))
# paste0(getwd(), "/reproducible/new_sim/100_300/rrzi_parametric.rds"))







