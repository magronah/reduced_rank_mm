setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")

library(foreach)
library(NBZIMM)


path  =   "~/scratch/dataset/new_sim/100_300/nbmm_mod/"
files =   list.files(path, full.names = TRUE)

res   = foreach(i = files,.packages = "NBZIMM") %do% {     
     mod   =   readRDS(i)
     est   =   unlist(lapply(mod$fit, function(x){coef(x)$fixed[["grouptreat"]]}))
     est
}


names_list    =  lapply(res, function(x){names(x)})
common_names  =  Reduce(intersect, names_list)
dd   =   as.data.frame(lapply(res, function(x){x[names(x) %in% common_names]}))


rownames(dd)  =    common_names
colnames(dd)  =    paste0("nsim",1:ncol(dd))
saveRDS(dd, file = paste0(getwd(), "/reproducible/new_sim/100_300/parametric_nbmm.rds"))
