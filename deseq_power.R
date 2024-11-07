setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)

directory_path <-   "~/scratch/dataset/power/deseq_res/" 
files          <-   list.files(directory_path, full.names = TRUE)

res = foreach(i = files, .combine = "cbind") %do% {
  res        =    readRDS(i)
  res$padj
}

res            =    as.data.frame(res)
colnames(res)  =    paste0("nsim", 1:ncol(res))
rownames(res)  =    paste0("taxon",1:nrow(res))

saveRDS(res, file = paste0(getwd(),"/reproducible/simp_single/data/deseq_power_padj.rds"))

