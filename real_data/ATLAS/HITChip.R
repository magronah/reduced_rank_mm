library(ANCOMBC)
library(tidyverse)
data(atlas1006, package = "microbiome")

x1         =   phyloseq(otu_table(atlas1006))
x2         =   phyloseq(tax_table(atlas1006))
meta_data  =   microbiome::meta(atlas1006)

meta_data$bmi = recode(meta_data$bmi_group,
                       obese = "obese",
                       severeobese = "obese",
                       morbidobese = "obese")

meta_data