setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")


####################################################################
path   =   paste0(nsubj,"_",ntaxa,"/sim_data/")
path
####################################################################
deseq_dd   =   readRDS(paste0(path,"deseq_otu_meta_list_withzi_taxa.rds"))
deseq_dd   =   readRDS(paste0(path,"deseq_Shrink_otu_meta_list_withzi_taxa.rds"))
