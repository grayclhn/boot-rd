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
boot_estimator_ForMc <- function(ypl, ypr, yql, yqr, fit.ql, fit.qr, 
                           coef.qr, coef.ql, coef.pr, coef.pl,
                           wqr, wql, ihr, ihl, Nbc, bootstrap = "wild") {
  
  fitted.l <- fit.ql %*% yql
  fitted.r <- fit.qr %*% yqr
  residual.l <- yql - fitted.l
  residual.r <- yqr - fitted.r
  boot_parameter <- coef.qr %*% yqr - coef.ql %*% yql
  estimate <- coef.pr %*% ypr - coef.pl %*% ypl

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
boot_dist_ForMc <- function(ypl, ypr, yql, yqr, fit.ql, fit.qr, 
                      coef.qr, coef.ql, coef.pr, coef.pl,
                      wqr, wql, ihr, ihl, Nbc, Nci, bootstrap = "wild"){
  
  fitted.l <- fit.ql %*% yql
  fitted.r <- fit.qr %*% yqr
  residual.l <- yql - fitted.l
  residual.r <- yqr - fitted.r
  
  if (bootstrap == "residual"){
    boots <- replicate(Nci, {
      newyqr <- fitted.r + sample(residual.r, length(yqr), T, wqr)
      newyql <- fitted.l + sample(residual.l, length(yql), T, wql)
      boot_estimator_ForMc(newyql[ihl], newyqr[ihr], newyql, newyqr, fit.ql, fit.qr, 
                     coef.qr, coef.ql, coef.pr, coef.pl,
                     wqr, wql, ihr, ihl, Nbc, bootstrap)})
  }
  if (bootstrap == "wild"){
    boots <- replicate(Nci, {
      newyqr <- fitted.r + residual.r*sample(wild_values, length(yqr), T, wild_weights)
      newyql <- fitted.l + residual.l*sample(wild_values, length(yql), T, wild_weights)
      boot_estimator_ForMc(newyql[ihl], newyqr[ihr], newyql, newyqr, fit.ql, fit.qr, 
                     coef.qr, coef.ql, coef.pr, coef.pl,
                     wqr, wql, ihr, ihl, Nbc, bootstrap)})
  }
  
  return(boots)
}

## a wrapper for both point and interval estimator
rdboot_ForMc <- function(y, x, a = 0.05, Nbc = 500, Nci = 999, p = 1, q = 2, 
                   bootstrap = "residual", kernel = "uniform"){
  
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
  
  wpl <- kweight(xpl, 0, h, kernel)
  wpr <- kweight(xpr, 0, h, kernel)
  wql <- kweight(xql, 0, b, kernel)
  wqr <- kweight(xqr, 0, b, kernel)
  
  Xpl <- cbind(1, poly(xpl, p, raw = T))
  Xpr <- cbind(1, poly(xpr, p, raw = T))
  Xql <- cbind(1, poly(xql, q, raw = T))
  Xqr <- cbind(1, poly(xqr, q, raw = T))
  
  Wpl <- diag(wpl)
  Wpr <- diag(wpr)
  Wql <- diag(wql)
  Wqr <- diag(wqr)
  
  fit.ql <- Xql %*% solve(t(Xql) %*% Wql %*% Xql) %*% t(Xql) %*% Wql
  fit.qr <- Xqr %*% solve(t(Xqr) %*% Wqr %*% Xqr) %*% t(Xqr) %*% Wqr
  coef.ql <- t(c(1, rep(0, q))) %*% solve(t(Xql) %*% Wql %*% Xql) %*% t(Xql) %*% Wql
  coef.qr <- t(c(1, rep(0, q))) %*% solve(t(Xqr) %*% Wqr %*% Xqr) %*% t(Xqr) %*% Wqr
  coef.pl <- t(c(1, rep(0, p))) %*% solve(t(Xpl) %*% Wpl %*% Xpl) %*% t(Xpl) %*% Wpl
  coef.pr <- t(c(1, rep(0, p))) %*% solve(t(Xpr) %*% Wpr %*% Xpr) %*% t(Xpr) %*% Wpr
  
  estimate <- boot_estimator_ForMc(ypl, ypr, yql, yqr, fit.ql, fit.qr, 
                             coef.qr, coef.ql, coef.pr, coef.pl,
                             wqr, wql, ihr, ihl, Nbc, bootstrap)
  
  boots <- boot_dist_ForMc(ypl, ypr, yql, yqr, fit.ql, fit.qr, 
                     coef.qr, coef.ql, coef.pr, coef.pl,
                     wqr, wql, ihr, ihl, Nbc, Nci, bootstrap)
  
  boot_parameter <- coef.qr %*% yqr - coef.ql %*% yql
  
  ci_basic <- boot_ci_basic(estimate, boots, boot_parameter, a)
  ci_percentile <- quantile(boots, c(a/2, 1-a/2))
  
  result <- rbind(c(estimate, ci_basic), c(estimate, ci_percentile))
  row.names(result) <- c("Basic CI", "Percentile CI")
  
  return(result)
}
