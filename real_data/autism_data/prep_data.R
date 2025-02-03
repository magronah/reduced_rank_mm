library(tidyr)
library(tidyverse)
library(DESeq2)
# library(glmmTMB)
# library(NBZIMM)
##############################################################
path1   =   paste0(getwd(),"/real_data/")
source(paste0(path1,"/fun.R"))
##############################################################
data        =   readRDS(paste0(path1,"autism_data/aut_data.rds"))
metadata    =   readRDS(paste0(path1,"autism_data/aut_metadata.rds"))
taxa_ex     =   readRDS(paste0(path1,"autism_data/taxa_exclude.rds"))
taxa_ex     
##############################################################
nam        =  "PRJNA644763"
count_dd   =   data[[nam]]
meta_dd    =   metadata[[nam]] %>% 
  setNames(c("subject", "group"))

dd_fil      =   filter_fun(count_dd,meta_dd,abund_thresh=7, 
                            sample_thresh=4)
dd_fil     =  t(dd_fil) %>% 
               as.data.frame() %>% 
               dplyr::select(-all_of(taxa_ex))

dd_filt    =  t(dd_fil)  %>% 
              as.data.frame()
dd         =   t(dd_filt)
ntax       =   ncol(dd)
##############################################################
dd_long    =   df_long(dd, otu_names = "sp", 
                       subject_name = "subject", 
                       ntaxa=ntax)

dd_long    =   left_join(dd_long, meta_dd, by ="subject")  
###########################################################
pp   =  deseqfun(dd_filt,meta_dd,alpha_level=0.1,
                 ref_name="ASD", 
                 ntaxa = nrow(dd_filt))
###########################################################
countdata  =  pp$data$countdata
meta_dd    =  pp$data$meta_data
normalize_fac =    sizeFactors(pp$object)
meta_dd$normalizer  =   normalize_fac
###########################################################
normalizer =   data.frame(normalizer = normalize_fac) %>% 
  rownames_to_column("subject")

df    =   left_join(dd_long, normalizer, by ="subject")  
###########################################################
gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")
