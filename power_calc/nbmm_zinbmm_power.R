library(tibble)
library(DESeq2)
library(Matrix)
library(huge)
library(foreach)
library(scam)
library(NBZIMM)
###############################################################
source("/project/6006158/agronahm/reduced_rank_mm/func2.R")
#source("/project/6006158/agronahm/reduced_rank_mm/initial_param0.R")

# Define sets of nsubj and ntaxa
parameter_sets <- list(
  list(nsubj = 50, ntaxa = 200),
  list(nsubj = 100, ntaxa = 300),
  list(nsubj = 150, ntaxa = 500)
)
