## Copyright (c) 2016, Otávio Bartalotti, Gray Calhoun, and Yang He.
## Available under the MIT "Expat" License, see README.md

# As its name suggest, the functions in this file is for fast simulation
# Note: observations beyound max(h, b) will be dropped at the beginning.
# Notation: 
# p    indicates objectes associated with p'th order polynomial; q indicates q'th order.
# l    indicates objectes at left hand side; r indicates right hand side.
# fit  indicates matrix generating fitted values.
# coef indicates matrix generating coefficients.
# ih   indicates elements within bandwith h.
# w    indicates kernel weights.


## this function generates bias-corrected estimator
boot_estimator_ForMc <- function(ypl, ypr, yql, yqr, 
                                 coef.ql, coef.qr, coef.pl, coef.pr,
                                 Xql, Xqr, Xpl, Xpr,
                                 WXql, WXqr,
                                 wqr, wql, ihr, ihl, Nbc, bootstrap = "wild") {
  
  b.ql       <- crossprod(WXql, yql)
  fitted.l   <- Xql %*% b.ql
  residual.l <- yql - fitted.l
  
  b.qr       <- crossprod(WXqr, yqr)
  fitted.r   <- Xqr %*% b.qr
  residual.r <- yqr - fitted.r
  
  boot_parameter <- coef.qr %*% yqr - coef.ql %*% yql
  estimate       <- coef.pr %*% ypr - coef.pl %*% ypl

  if (bootstrap == "residual") {
    boots <- replicate(Nbc,
             coef.pr %*% (fitted.r[ihr] + sample(residual.r, length(ypr), T, wqr)) -
             coef.pl %*% (fitted.l[ihl] + sample(residual.l, length(ypl), T, wql)))
  }
  if (bootstrap == "wild") {
    boots <- replicate(Nbc,
             coef.pr %*% (fitted.r[ihr] + residual.r[ihr]*sample(wild_values, length(ypr), T, wild_weights)) -
             coef.pl %*% (fitted.l[ihl] + residual.l[ihl]*sample(wild_values, length(ypl), T, wild_weights)))
  }

  
  return(as.numeric(estimate - mean(boots) + boot_parameter))
}

## this function generates bootstrap distribution of bias-corrected estimator
boot_dist_ForMc <- function(ypl, ypr, yql, yqr, 
                            coef.ql, coef.qr, coef.pl, coef.pr,
                            Xql, Xqr, Xpl, Xpr,
                            WXql, WXqr,
                            wqr, wql, ihr, ihl, Nbc, Nci, bootstrap = "wild"){
  
  b.ql       <- crossprod(WXql, yql)
  fitted.l   <- Xql %*% b.ql
  residual.l <- yql - fitted.l
  
  b.qr       <- crossprod(WXqr, yqr)
  fitted.r   <- Xqr %*% b.qr
  residual.r <- yqr - fitted.r
  
  if (bootstrap == "residual"){
    boots <- replicate(Nci, {
      newyqr <- fitted.r + sample(residual.r, length(yqr), T, wqr)
      newyql <- fitted.l + sample(residual.l, length(yql), T, wql)
      boot_estimator_ForMc(newyql[ihl], newyqr[ihr], newyql, newyqr,
                           coef.ql, coef.qr, coef.pl, coef.pr,
                           Xql, Xqr, Xpl, Xpr,
                           WXql, WXqr,
                           wqr, wql, ihr, ihl, Nbc, bootstrap)})
  }
  if (bootstrap == "wild"){
    boots <- replicate(Nci, {
      newyqr <- fitted.r + residual.r*sample(wild_values, length(yqr), T, wild_weights)
      newyql <- fitted.l + residual.l*sample(wild_values, length(yql), T, wild_weights)
      boot_estimator_ForMc(newyql[ihl], newyqr[ihr], newyql, newyqr,
                           coef.ql, coef.qr, coef.pl, coef.pr,
                           Xql, Xqr, Xpl, Xpr,
                           WXql, WXqr,
                           wqr, wql, ihr, ihl, Nbc, bootstrap)})
  }
  
  return(boots)
}

## a wrapper for both point and interval estimator
rdboot_ForMc <- function(y, x, a = 0.05, Nbc = 500, Nci = 999, p = 1, q = 2, 
                         bootstrap = "residual", kernel = "uniform"){
  
  # to improve speed, (1) create all necessary objects only once and pass
  # them into functions for double bootstrap. (2) avoid multiplication of
  # matrix with large size.
  
  # drop all observations beyound max[h, b]
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
  
  estimate <- boot_estimator_ForMc(ypl, ypr, yql, yqr, 
                                   coef.ql, coef.qr, coef.pl, coef.pr,
                                   Xql, Xqr, Xpl, Xpr,
                                   WXql, WXqr,
                                   wqr, wql, ihr, ihl, Nbc, bootstrap)
  
  boots <- boot_dist_ForMc(ypl, ypr, yql, yqr, 
                           coef.ql, coef.qr, coef.pl, coef.pr,
                           Xql, Xqr, Xpl, Xpr,
                           WXql, WXqr,
                           wqr, wql, ihr, ihl, Nbc, Nci, bootstrap)
  
  boot_parameter <- coef.qr %*% yqr - coef.ql %*% yql
  
  ci_basic <- boot_ci_basic(estimate, boots, boot_parameter, a)
  ci_percentile <- quantile(boots, c(a/2, 1-a/2))
  
  result <- rbind(c(estimate, ci_basic), c(estimate, ci_percentile))
  row.names(result) <- c("Basic CI", "Percentile CI")
  
  return(result)
}
