library(phyloseq)
library(dplyr)
library(fido)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(microbiome)
library(vegan)
library(pheatmap)
##############################################################
path   =   paste0(getwd(),"/real_data/")
source(paste0(path,"/fun.R"))
##############################################################
data("atlas1006")

countdata  <-  data.frame(abundances(atlas1006)) 
colnames(countdata) <-  paste0("Sample-",1:ncol(countdata))
rownames(countdata) <-  paste0("taxon",1:nrow(countdata))

meta_data  <-  sample_data(atlas1006)
meta_dd    <-  meta_data %>%
               data.frame()  %>%
               dplyr::select(c("sample","bmi_group","age")) %>%
               dplyr::filter(meta_data[["time"]] == 0)  %>%
               na.omit() %>%
               setNames(c("subject","group", "age")) 

meta_dd$group <- factor(meta_dd$group, 
                        levels = c("lean", "underweight", "overweight", 
                                   "obese", "severeobese", "morbidobese"))


data  <-   countdata %>%  
              dplyr::select(rownames(meta_dd))
data  <-   data.frame(t(data))
###############################################################################
dd  <- data[, apply(data, 2, function(col) any(col == 0))]
dd  <- dd[, apply(dd, 2, function(col) sum(col) != 0)]

ntax    =  ncol(dd)
##############################################################
dd_long    =   df_long(dd, otu_names = "taxon", 
                       subject_name = "subject", 
                       ntaxa=ntax)

dd_long    =   left_join(dd_long, meta_dd, by ="subject")  
###########################################################
pp         =    deseqfun(dd, meta_dd, alpha_level=0.1,
                         design   = ~group + age,
                         ref_name = "lean",
                         ntaxa = ntax)
###########################################################
countdata     =    pp$data$countdata
meta_dd       =    pp$data$meta_data
normalize_fac =    sizeFactors(pp$object)
meta_dd$normalizer  =   normalize_fac
###########################################################
normalizer =   data.frame(normalizer = normalize_fac) %>% 
  rownames_to_column("subject")

df    =   left_join(dd_long, normalizer, by ="subject")  
################################################################
gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

