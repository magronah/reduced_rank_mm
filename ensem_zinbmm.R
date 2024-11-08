library(NBZIMM)
library(DESeq2)
path = paste0(getwd(),"/100_300/")
####################################################################
data	=   readRDS(paste0(path,"otu_meta_list_withzi_taxa.rds"))
################################################################
cc	=   commandArgs(trailingOnly  = TRUE)
#i	=   as.integer(cc[1])
i=1
dd	=   data[[i]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data
met_dd$dummy = factor(1)

otu_count    =   t(countdata)
  dds        =   DESeqDataSetFromMatrix(otu_count,met_dd, ~group)
  dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
  normalizer =   sizeFactors(dds) 

mod        =   mms(y = countdata, fixed = ~group + offset(normalizer),
                      random =  ~ 1|dummy,
                      zi_fixed = ~1,
                      data = met_dd, method = "zinb")

saveRDS(mod, file=paste0("~/scratch/dataset/RR/100_300/zinbmm/mod",i,".rds"))





