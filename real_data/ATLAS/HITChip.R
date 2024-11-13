set.seed(123)
otu_table <- as.data.frame(matrix(sample(1:1000, 50, replace = TRUE), nrow = 10))
colnames(otu_table) <- paste0("OTU", 1:5)
rownames(otu_table) <- paste0("Subject", 1:10)

metadata <- data.frame(
  subject = rownames(otu_table),
  group = c(rep("A", 5), rep("B", 5)) # Two groups: A and B
)

rownames(metadata) <- rownames(otu_table)
library(ANCOMBC)
library(phyloseq)

OTU = otu_table(as.matrix(otu_table), taxa_are_rows = FALSE)
sample_data = sample_data(metadata)
physeq = phyloseq(OTU, sample_data)

# Run ANCOMBC
ancombc_results <- ancombc(phyloseq = physeq,
                           formula = "group", # Use group as the variable of interest
                           p_adj_method = "holm", # Adjust p-values using Holm's method
                           zero_cut = 0.90, # Threshold for filtering low-abundance taxa
                           lib_cut = 1000, # Minimum library size to retain samples
                           group = "group", # Specify group variable
                           struc_zero = TRUE, # Structural zeros for zero-inflation
                           global = TRUE)

# Extract results
res <- ancombc_results$res

# Estimated fold changes
fc_table <- res$beta
fc_table

 

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


meta_data$bmi = factor(meta_data$bmi, levels = c("obese", "overweight", "lean"))


