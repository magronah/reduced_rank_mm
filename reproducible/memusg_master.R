## remotes::install_github("mrc-ide/memprof")
source("reproducible/leverage_funs.R")
library(memprof)
library(glmmTMB, lib = "glmmTMB_lev")
library(peakRAM)
library(RTMB)
library(Matrix)

testframe <- expand.grid(nsubj = seq(5, 40, by = 5),
                         ntax = seq(25, 200, by = 25))
write.csv(testframe, "reproducible/testvals.csv")

res <- list()
for (i in seq(nrow(testframe))) {
    nsubj <- testframe$nsubj[i]
    ntax <- testframe$ntax[i]
    cat(i, nsubj, ntax, "\n")
    ## system(sprintf("./memusg Rscript <reproducible/memusg_worker.R >reproducible/mw%d.Rout 2>reproducible/mw%d_mem.Rout %d %d",
    ##                     i, i, nsubj, ntax))
    res[[i]] <- with_monitor(peakRAM_testfun(nsubj, ntax))
    saveRDS(res, "memusg.rds")
}


