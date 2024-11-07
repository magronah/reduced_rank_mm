setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(NBZIMM)
path = paste0(getwd(),"/reproducible/new_sim/100_300/")
####################################################################
data	=   readRDS(paste0(path,"otu_meta_list_withzi_taxa.rds"))
################################################################
cc	=   commandArgs(trailingOnly  = TRUE)
i	=   as.integer(cc[1])
dd	=   data[[i]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data
met_dd$dummy = factor(1)

mod        =   mms(y = countdata, fixed = ~group,
                      random = ~ 1|dummy,
                      zi_fixed = ~1,
                      data = met_dd, method = "zinb")

saveRDS(mod, file=paste0("~/scratch/dataset/new_sim/100_300/zinbmm/mod",i, ".rds"))





