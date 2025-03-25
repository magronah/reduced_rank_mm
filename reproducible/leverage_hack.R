library(glmmTMB)
mm <- readRDS("reproducible/rr_mod.rds")

## attempt leverage calculation

## reconstruct data from model frame
mf <- model.frame(mm)
newdata <- with(mf,
                data.frame(count, group, taxon, site,
                           normalizer = `offset(normalizer)`))
dim(newdata)
length(unique(newdata$taxon))

## get starting parameters 
p0 <- with(mm$obj$env, parList(last.par.best[-random]))
## exclude empty (unused) parameters
p0 <- p0[lengths(p0) > 0]
## exclude random effects
p0 <- p0[setdiff(names(p0), "b")]

## test (6 seconds)
system.time(
    lik1 <- with(mm$obj$env, mm$obj$fn(last.par.best[-random]))
)
lik2 <- with(mm$obj$env, mm$obj$fn(unlist(p0)))

## might not want threading if we're doing the leverage calculation in parallel ...
par_ctrl <- glmmTMBControl(
    parallel = list(n = 10, autopar = TRUE)
)

gprior <- data.frame(prior = "gamma(2, 2.5)",
                     class = "theta_sd",
                     coef = "")
## par_ctrl <- glmmTMBControl(optCtrl = list(trace=1))  ## default
par_ctrl <- glmmTMBControl()

options(glmmTMB_openmp_debug = FALSE)
calc_leverage <- function(i, eps = 0.001) {
    ## i <- 800; eps <- 0.001
    newdata$count[i]
    pert_data <- newdata
    pert_data$count[i] <- pert_data$count[i] + eps
    ## 2 minutes to optimize?
    st <- system.time(
        newfit <- update(mm, start = p0, data = pert_data, verbose = FALSE,
                         control = par_ctrl, priors = gprior)
    )
    ## newfit <- update(mm, start = p0, data = pert_data, verbose = FALSE,
    ##                     control = par_ctrl, priors = gprior)

    ## can we do this modularly and skip the slow stuff after the optimization?
    newfit0 <- update(mm, start = p0, data = pert_data, verbose = FALSE,
                      control = par_ctrl, priors = gprior, doFit = FALSE)
    system.time(newfit1 <- fitTMB(newfit0, doOptim = FALSE)) ## 1 second
    system.time(newfit2 <- with(newfit1, nlminb(par, fn, gr))) ## 80 seconds
    ## then we should be
    X <- getME(mm, "X")
    Z <- getME(mm, "Z")
    pp <- with(newfit1$env, parList(last.par.best[-random]))
    ## perturbed prediction
    pred <- drop(X[i,] %*% pp[["beta"]] + Z[i,] %*% pp[["b"]])
    ## original prediction
    pred0 <- fitted(mm)[i]
    (pred - pred0) / eps
}
