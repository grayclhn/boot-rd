# Verify that new versions of the functions in 'rdfunctions_ForMc.R' return
# the same results as functions in 'rdfunctions_old.R' in residual bootstrap.

# Verify that new versions of the functions in 'rdfunctions_ForMc.R' return
# the same results as functions in 'rdfunctions.R' in wild bootstrap.

library(rdrobust)
library(testthat)
library(rbenchmark)

source("rdfunctions.R")
source("rdfunctions_old.R")
source("rdfunctions_ForMc.R")

test_that("rdfunctions_ForMc.R gives the same results as rdfunctions_old.R in residual bootstrap.",{
  
  # generate a wrapper from rdfunctions_old.R, this allows
  # direct comparison with _ForMc results.
  old_estimator_wrapper <- function(y, x, a = 0.05, Nbc = 500, Nci = 999, p = 1, q = 2, 
                                    kernel = "uniform"){
    
    bw <- rdbwselect_2014(y, x, p=p, q=q, kernel=kernel)$bws
    h <- bw[1]
    b <- bw[2]
    
    estimate <- rd.estimate(data.frame(y=y, x=x), p, q, h, b, Nbc, kernel, T)
    ci <- rd.ci(data.frame(y=y, x=x), p, q, h, b, Nbc, Nci, 1-a, kernel, T)$ci
    
    return(c(estimate, ci))
  }
  
  dta <- generate.data(1)
  
  set.seed(798)
  old_uni <- old_estimator_wrapper(dta$y, dta$x, Nbc = 50, Nci = 50, kernel = "uniform")
  old_tri <- old_estimator_wrapper(dta$y, dta$x, Nbc = 50, Nci = 50, kernel = "triangular")
  old_epa <- old_estimator_wrapper(dta$y, dta$x, Nbc = 50, Nci = 50, kernel = "epanechnikov")
  
  set.seed(798)
  new_uni <- rdboot_ForMc(dta$y, dta$x, Nbc = 50, Nci = 50, bootstrap = "residual", 
                          kernel = "uniform", residual = "HC0")[2, ]
  new_tri <- rdboot_ForMc(dta$y, dta$x, Nbc = 50, Nci = 50, bootstrap = "residual", 
                          kernel = "triangular", residual = "HC0")[2, ]
  new_epa <- rdboot_ForMc(dta$y, dta$x, Nbc = 50, Nci = 50, bootstrap = "residual", 
                          kernel = "epanechnikov", residual = "HC0")[2, ]
  
  expect_equivalent(old_uni, new_uni)
  expect_equivalent(old_tri, new_tri)
  expect_equivalent(old_epa, new_epa)})


test_that("rdfunctions_ForMc.R gives the same results as rdfunctions.R in wild bootstrap.",{
  
  # rdfunctions.R take random draws for all observations in original data set,
  # rdrunctions_ForMc.R only take random draws for observations used in estimation,
  # this makes it impossible to compare in general. So we test the case h = b.
  
  # because both rdboot and rdboot_ForMc calculate bandwidth automatically and there
  # is no option to control that, we instead compare:
  # (1) point estimator:  boot_estimator VS boot_estimator_ForMc
  # (2) distribution of point estimator: boot_dist VS boot_dist_ForMc
  # boot_dist does not exist, we use boot_interval and return the object boots
  # from it.
  
  
  # first, we generate a dataset as usual, and drop all observations beyound b
  dta <- generate.data(1)
  y <- dta$y
  x <- dta$x
  p <- 1
  q <- 2
  kernel <- "uniform"
  
  bw <- rdbwselect_2014(y, x, p=p, q=q, kernel=kernel)$bws
  h <- bw[2]
  b <- bw[2]
  
  yql <- y[x > -max(bw) & x < 0]
  xql <- x[x > -max(bw) & x < 0]
  yqr <- y[x >= 0 & x < max(bw)]
  xqr <- x[x >= 0 & x < max(bw)]
  
  ihl <- xql > -h
  ihr <- xqr < h
  
  ypl <- y[x > -h & x < 0]
  xpl <- x[x > -h & x < 0]
  ypr <- y[x >= 0 & x < h]
  xpr <- x[x >= 0 & x < h]
  
  # a vector of weight for residual bootstrap
  wql <- kweight(xql, 0, b, kernel)
  wqr <- kweight(xqr, 0, b, kernel)
  wpl <- kweight(xpl, 0, b, kernel)
  wpr <- kweight(xpr, 0, b, kernel)
  
  # orthogonal polynomials
  xql.poly <- poly(xql, q)
  xqr.poly <- poly(xqr, q)
  xpl.poly <- poly(xpl, p)
  xpr.poly <- poly(xpr, p)
  
  # design matrix
  Xql <- cbind(1, poly(xql, q))
  Xqr <- cbind(1, poly(xqr, q))
  Xpl <- cbind(1, poly(xpl, p))
  Xpr <- cbind(1, poly(xpr, p))
  
  KXql <- kweight(xql, 0, b, kernel) * Xql
  KXqr <- kweight(xqr, 0, b, kernel) * Xqr
  KXpl <- kweight(xpl, 0, h, kernel) * Xpl
  KXpr <- kweight(xpr, 0, h, kernel) * Xpr
  
  # parameter maker
  WXql <- t(solve(crossprod(Xql, KXql), t(KXql)))
  WXqr <- t(solve(crossprod(Xqr, KXqr), t(KXqr)))
  WXpl <- t(solve(crossprod(Xpl, KXpl), t(KXpl)))
  WXpr <- t(solve(crossprod(Xpr, KXpr), t(KXpr)))
  
  # intercept maker
  coef.ql <- tcrossprod(c(1, predict(xql.poly, 0)), WXql)
  coef.qr <- tcrossprod(c(1, predict(xqr.poly, 0)), WXqr)
  coef.pl <- tcrossprod(c(1, predict(xpl.poly, 0)), WXpl)
  coef.pr <- tcrossprod(c(1, predict(xpr.poly, 0)), WXpr)
  
  # second, generate boot_dist from boot_interval
  boot_dist <- function(y, x, wp, wq, a, p = 1, q = p + 1,
                            nboot, bootfn, nboot2 = nboot, type = c("basic", "percentile", "both"),
                            m0 = lm(y ~ poly(x, q), subset = x <= 0, weights = wq),
                            m1 = lm(y ~ poly(x, q), subset = x > 0, weights = wq)) {
    
    type <- match.arg(type)
    
    yboot <- rep(NA, length(y))
    i0 <- x <= 0
    i1 <- !i0
    
    boots <- replicate(nboot2, {
      yboot[i1] <- bootfn(m1, y[i1], x[i1])
      yboot[i0] <- bootfn(m0, y[i0], x[i0])
      boot_estimator(yboot, x, wp, wq, p, q, nboot, bootfn)
    })
    
    return(boots)
  }
  
  # now we can compare
  set.seed(798)
  new_estimate <- boot_estimator_ForMc(ypl, ypr, yql, yqr,
                                       coef.ql, coef.qr, coef.pl, coef.pr,
                                       Xql, Xqr, Xpl, Xpr,
                                       WXql, WXqr, 1, 1, 
                                       wqr, wql, ihr, ihl, 50, bootstrap = "wild")
  
  new_dist <- boot_dist_ForMc(ypl, ypr, yql, yqr,
                              coef.ql, coef.qr, coef.pl, coef.pr,
                              Xql, Xqr, Xpl, Xpr,
                              WXql, WXqr, 1, 1,
                              wqr, wql, ihr, ihl, 10, 10, bootstrap = "wild")
  
  set.seed(798)
  old_estimate <- boot_estimator(c(yql, yqr), c(xql, xqr), c(wpl, wpr), c(wql, wqr),
                                 1, 2, 50, wild_bootstrap)
  
  old_dist <- boot_dist(c(yql, yqr), c(xql, xqr), c(wpl, wpr), c(wql, wqr),
                                     0.05, 1, 2, 10, wild_bootstrap, 10)
  
  expect_equivalent(new_estimate, old_estimate)
  expect_equivalent(new_dist, old_dist)})


test_that("rdfunctions_ForMc.R correctly calculates diagonal elements of projection matrix.",{
  
  # generate data
  dta <- generate.data(1)
  y <- dta$y
  x <- dta$x
  p <- 1
  q <- 2
  kernel <- "tri"
  
  # redo the steps in the function "rdboot_ForMc"
  bw <- rdbwselect_2014(y, x, p=p, q=q, kernel=kernel)$bws
  h <- bw[1]
  b <- bw[2]
  
  yql <- y[x > -max(bw) & x < 0]
  xql <- x[x > -max(bw) & x < 0]
  yqr <- y[x >= 0 & x < max(bw)]
  xqr <- x[x >= 0 & x < max(bw)]
  
  ihl <- xql > -h
  ihr <- xqr < h
  
  ypl <- y[x > -h & x < 0]
  xpl <- x[x > -h & x < 0]
  ypr <- y[x >= 0 & x < h]
  xpr <- x[x >= 0 & x < h]
  
  # a vector of weight for residual bootstrap
  wql <- kweight(xql, 0, b, kernel)
  wqr <- kweight(xqr, 0, b, kernel)
  
  # orthogonal polynomials
  xql.poly <- poly(xql, q)
  xqr.poly <- poly(xqr, q)
  xpl.poly <- poly(xpl, p)
  xpr.poly <- poly(xpr, p)
  
  # design matrix
  Xql <- cbind(1, poly(xql, q))
  Xqr <- cbind(1, poly(xqr, q))
  Xpl <- cbind(1, poly(xpl, p))
  Xpr <- cbind(1, poly(xpr, p))
  
  KXql <- kweight(xql, 0, b, kernel) * Xql
  KXqr <- kweight(xqr, 0, b, kernel) * Xqr
  KXpl <- kweight(xpl, 0, h, kernel) * Xpl
  KXpr <- kweight(xpr, 0, h, kernel) * Xpr
  
  # parameter maker
  WXql <- t(solve(crossprod(Xql, KXql), t(KXql)))
  WXqr <- t(solve(crossprod(Xqr, KXqr), t(KXqr)))
  WXpl <- t(solve(crossprod(Xpl, KXpl), t(KXpl)))
  WXpr <- t(solve(crossprod(Xpr, KXpr), t(KXpr)))
  
  # intercept maker
  coef.ql <- tcrossprod(c(1, predict(xql.poly, 0)), WXql)
  coef.qr <- tcrossprod(c(1, predict(xqr.poly, 0)), WXqr)
  coef.pl <- tcrossprod(c(1, predict(xpl.poly, 0)), WXpl)
  coef.pr <- tcrossprod(c(1, predict(xpr.poly, 0)), WXpr)
  
  # HC adjustment: diagnal vector of the projection matrix
  KXql.sqrt <- sqrt(kweight(xql, 0, b, kernel)) * Xql
  KXqr.sqrt <- sqrt(kweight(xqr, 0, b, kernel)) * Xqr
  
  # below are the fast version of diagonal elements
  hl <- diag(KXql.sqrt %*% solve(crossprod(Xql, KXql), t(KXql.sqrt)))
  hr <- diag(KXqr.sqrt %*% solve(crossprod(Xqr, KXqr), t(KXqr.sqrt)))
  
  ## below are the slower version
  # the weighted design matrix
  Xl <- sqrt(kweight(xql, 0, b, kernel)) * cbind(1, poly(xql, 2, raw = T))
  Xr <- sqrt(kweight(xqr, 0, b, kernel)) * cbind(1, poly(xqr, 2, raw = T))
  
  hl.slow <- diag(Xl %*% solve(t(Xl) %*% Xl) %*% t(Xl))
  hr.slow <- diag(Xr %*% solve(t(Xr) %*% Xr) %*% t(Xr))
  
  expect_equivalent(hl, hl.slow)
  expect_equivalent(hr, hr.slow)})

