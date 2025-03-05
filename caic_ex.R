
## need to modify src/Makevars in glmmTMB directory to contain this
## (no fopenmp!)
c("PKG_CPPFLAGS = -DTMBAD_FRAMEWORK -DTMBAD_INDEX_TYPE=uint64_t",
  "## PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)",
  "## PKG_CXXFLAGS=$(SHLIB_OPENMP_CXXFLAGS)")

## need to install this way (with location adjusted to your liking)
## R CMD INSTALL -l ~/students/agronah/reduced_rank_mm/glmmTMB_lev glmmTMB

## from https://github.com/glmmTMB/glmmTMB/blob/leverage/misc/leverage.R
## @param diag Get diagonal only?
leverage <- function(fm, diag=TRUE) {
    has.random <- any(fm$obj$env$lrandom())
    obj <- fm$obj
    ## We mess with these... (cleanup on exit!)
    restore.on.exit <- c("ADreport",
                         "parameters",
                         "data")
    oldvars <- sapply(restore.on.exit, get, envir=obj$env, simplify=FALSE)
    restore.oldvars <- function(){
        for(var in names(oldvars)) assign(var, oldvars[[var]], envir=obj$env)
    }
    on.exit({restore.oldvars(); obj$retape()})
    ## #################################################################
    ## Get derivatives of prediction
    ##
    ##    mu_hat( b_hat( theta_hat(yobs), yobs) , theta_hat(yobs) )
    ##
    ## wrt. yobs.
    ## Note the three partial derivative 'paths' to consider.
    ## We can get this derivative by
    ##  1. yobs -> theta_hat(yobs)
    ##  2. theta -> mu_hat( b_hat( theta, yobs) , theta ) [fixed yobs]
    ##  3. yobs -> mu_hat( b_hat( theta, yobs) , theta ) [fixed theta]
    ## #################################################################
    ##parhat <- obj$env$last.par.best
    parhat <- fm$fit$parfull
    pl <- obj$env$parList(par=parhat)
    yobs <- obj$env$data$yobs
    obj$env$parameters <- pl
    theta <- parhat[!obj$env$lrandom()]  ## ALL top-level parameters
    b <- parhat[obj$env$lrandom()]
    if (!is.null(obj$env$spHess)) {
        Hbb <- obj$env$spHess(parhat, random=TRUE) ## Needed later for RE models
    }
    ## #################################################################
    ## 1. Get partial derivatives of theta_hat wrt to yobs
    ## Note: length(yobs) much greater that length(theta)
    ##       ==> Reverse mode AD is suitable !
    ## #################################################################
    ## Move 'yobs' from data -> parameters (preserve C++ template order!)
    obj$env$parameters <- c(list(yobs = yobs), obj$env$parameters)
    obj$env$data$yobs <- NULL
    obj$retape()
    ## New objective parameters: (yobs, b, theta)
    nobs <- length(obj$env$parameters$yobs)
    nb <- length(obj$env$random)
    ntheta <- length(obj$env$par) - nobs - nb
    TMB::config(tmbad.atomic_sparse_log_determinant=0, DLL="RTMB") ## TMB FIXME
    F <- GetTape(obj)
    r <- obj$env$random ## Which are random
    p <- tail(1:(nobs+ntheta), ntheta) ## Which are parameters *after* removing random
    ThetaHat <- F$laplace(r)$newton(p)
    J <- ThetaHat$jacobian(ThetaHat$par())
    ## Extra stuff we need in (3)
    F. <- F$jacfun() ## (yobs, [b], theta) -> (yobs, [b], theta)
    F. <- MakeTape(function(y) F.( c(y, parhat) ) [r] , yobs) ## yobs -> b
    Hby <- F.$jacfun(sparse=TRUE)(yobs)
    ## #################################################################
    ## 2. Get partial derivatives of mu_hat wrt to theta for *fixed* yobs
    ## Note: length(mu) much greater that length(theta)
    ##       ==> Forward mode AD is suitable !
    ## #################################################################
    obj$env$data$yobs <- yobs
    obj$env$parameters$yobs <- NULL
    obj$retape()
    r <- obj$env$random ## Which are now random
    F <- GetTape(obj)
    Bhat <- F$newton(r) ## theta -> bhat
    obj$env$data$doPredict <- as.double(1) ## Enable prediction of 'mu'
    obj$env$data$whichPredict <- as.double(1:nobs)
    obj$env$ADreport <- TRUE ## Modify return value from Cpp
    obj$retape()
    F <- GetTape(obj) ## (b, theta) -> mu
    ## This doesn't work:
    ## MuHat <- MakeTape(function(theta)F(c(Bhat(theta), theta)), theta)
    MuHat <- MakeTape(function(theta) {
        par <- advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
        r <- obj$env$lrandom()
        par[r] <- Bhat(theta)
        par[!r] <- theta
        F(par)
    } , theta)
    ## 'Adjoint trick'
    T2 <- MakeTape(function(weight) {
        WT <- MakeTape(function(th) sum(MuHat(th) * weight), theta)
        WT$jacfun()(advector(theta))
    }, rep(1,nobs))
    J2 <- T2$jacobian(yobs)
    if (diag) {
        term1 <- colSums(J*J2)
    } else {
        term1 <- t(J) %*% J2
    }
    ## #################################################################
    ## 3. Get partial derivatives of mu_hat wrt yobs for fixed theta
    ## Note: Tricky!
    ##       
    ## #################################################################
    term2 <- 0
    if (has.random) {
        F2 <- MakeTape(function(b) {
            par <- advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
            r <- obj$env$lrandom()
            par[r] <- b
            par[!r] <- theta
            F(par)
        }, b) ## (b) -> mu
        Hmb <- F2$jacfun(sparse=TRUE)(b) ## sparse deriv mu wrt b
        ## Implicit function theorem gives final partial deriv matrix:
        ##   - Hby %*% solve(Hbb) %*% Hbm
        ## of which we need the diagonal.
        ## Because mu and yobs link to the same random effects, all required b-cliques are part of Hbb !
        ## It follows that we can replace solve(Hbb) by its subset iH !
        if (diag) {
            if (length(b) == 0) {
                iH <- solve(Hbb)
            } else {
                iH <- TMB:::solveSubset(Hbb)
            }
            term2 <- -colSums( Hby * (  iH %*% t(Hmb) ) )
        } else {
            term2 <- -t(Hby) %*% solve(Hbb, t(Hmb))
        }
    }
    term1 + term2
}

library(glmmTMB, lib.loc = "glmmTMB_lev")
## both of these packages *must* be loaded
library(RTMB)
library(Matrix)
library(lme4)

data("sleepstudy", package = "lme4")
fm1 <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
fm0 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)

X <- getME(fm1, "X")
Z <- getME(fm1, "Z")
## compute _conditional_ log likelihood of glmmTMB object?
pp <- fm1$obj$env$last.par.best
pp0 <- with(fm1$obj$env, last.par.best[-random])
ppL <- split(pp, names(pp))
sigma(fm1)

## this successfully recovers the conditional log-likelihood as given by cAIC.
-sum(dnorm(sleepstudy$Reaction, drop(X %*% ppL[["beta"]] + Z%*% ppL[["b"]]), sigma(fm1), log = TRUE))
## so does this
-sum(dnorm(sleepstudy$Reaction, fitted(fm1), sigma(fm1), log = TRUE))


logLik(fm1)
plot(leverage(fm1), hatvalues(fm0))
abline(a=0, b=1, lty = 2)


## want to calculate the conditional log-likelihood (turn off Laplace approx)

fm1$obj$fn(pp0)
fm1$obj$env$random <- numeric(0)
m1$obj$retape()
## do something to zero out the log-likelihood of b with respect to Sigma?
pp2 <- ppL
pp2[["theta"]] <- c(100, 100, 0)
fm1$obj$fn(unlist(pp2))

## 6 top-level params
length(ppL[["b"]])  ## 36
## close, anyway ...
## NOTE leverage calc messes up the guts of the TMB object
sum(leverage(fm1))
sum(hatvalues(fm0))
## try setting up a sim with singular fit so equivalent to lm hatvalues??

library(cAIC4)
cAIC(fm0)

logLik(fm1)
logLik(fm0)


mm <- readRDS("reproducible/rr_mod.rds")
system.time(trace_hat <- sum(leverage(mm)))

## conditional AIC (no Z-I ...)

## can we get the conditional log-likelihood from the sum of squares of the deviance residuals???
## this *might* be it ... ??
-sum(dnbinom(model.frame(mm)$count, mu = fitted(mm), size = sigma(mm), log = TRUE))
