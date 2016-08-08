
local_polynomial_regression <- function(dta, order, bw, kernel) {
  weight <- kweight(dta$x, 0, bw, kernel)
  fit    <- lm.wfit(cbind(1, poly(dta$x, order, raw = T)), dta$y, weight)
  return(list(cons = fit$coefficients[1],
              yhat = fit$fitted.values,
              res = fit$residuals,
              weight = weight))
}

rd.estimate <- function(dta, p, q, bw.p, bw.q, N.bc, kernel, bc = TRUE) {
  # dta:    a data set contains y and x.
  # p:      the order of polynomial for point estimation.
  # q:      the order of polynomial for bias correction.
  # bw.p:   bandwidth for point estimation.
  # bw.q:   bandwidth for bias correction.
  # N.bc:   number of bootstrap replications for bias correction.
  # kernel: epanechnikov, uniform or triangular.
  # w.q:    kernel weight for bias correction. Reduce repeated computation.
  # X:      the regressor matrix. Reduce repeated computation.
  # bc:     perform bias correction?

  if (length(bw.p) == 1) bw.p <- rep(bw.p, 2)
  if (length(bw.q) == 1) bw.q <- rep(bw.q, 2)

  # estimated effect from low order polynomial
  dta.p.l <- dta[dta$x > -bw.p[1] & dta$x < 0, ]
  dta.p.r <- dta[dta$x < bw.p[2] & dta$x >= 0, ]
  trd.hat <- local_polynomial_regression(dta.p.r, p, bw.p[2], kernel)$cons -
    local_polynomial_regression(dta.p.l, p, bw.p[1], kernel)$cons

  if (!bc) {
    return(trd.hat)
  } else {

    # the "true" effect from high order polynomial
    dta.q.l <- dta[dta$x > -bw.q[1] & dta$x < 0, ]
    dta.q.r <- dta[dta$x < bw.q[2] & dta$x >= 0, ]
    left  <- local_polynomial_regression(dta.q.l, q, bw.q[1], kernel)
    right <- local_polynomial_regression(dta.q.r, q, bw.q[2], kernel)
    trd   <- right$cons - left$cons

    # bootstrap
    trd.hat.star <- replicate(N.bc,
      local_polynomial_regression(data.frame(
        x = dta.p.r$x,
        y = right$yhat[dta.q.r$x < bw.p[2]] +
          sample(right$res, nrow(dta.p.r), T, right$weight)),
        p, bw.p[2], kernel)$cons -
      local_polynomial_regression(data.frame(
        x = dta.p.l$x,
        y = left$yhat[dta.q.l$x > -bw.p[1]] +
          sample(left$res, nrow(dta.p.l), T, left$weight)),
        p, bw.p[1], kernel)$cons
    )
    return(trd.hat - (mean(trd.hat.star) - trd))
  }
}

rd.ci <- function(dta, p, q, bw.p, bw.q, N.bc, N.ci, level, kernel, bc = T) {
  # dta:    a data set contains y and x.
  # p:      the order of polynomial for point estimation.
  # q:      the order of polynomial for bias correction.
  # bw.p:   bandwidth for point estimation.
  # bw.q:   bandwidth for bias correction.
  # N.bc:   number of bootstrap replications for bias correction.
  # N.ci:   number of bootstrap replications for confidence interval.
  # kernel: epanechnikov, uniform or triangular.
  # bc:     perform bias correction?

  if (length(bw.p) == 1) bw.p <- rep(bw.p, 2)
  if (length(bw.q) == 1) bw.q <- rep(bw.q, 2)

  # fit the high order polynomial on both sides
  dta.q.l <- dta[dta$x > -bw.q[1] & dta$x < 0, ]
  dta.q.r <- dta[dta$x < bw.q[2] & dta$x >= 0, ]
  left  <- local_polynomial_regression(dta.q.l, q, bw.q[1], kernel)
  right <- local_polynomial_regression(dta.q.r, q, bw.q[2], kernel)

  # bootstrap
  trd.star <- replicate(N.ci, rd.estimate(data.frame(
    x = c(dta.q.r$x, dta.q.l$x),
    y = c(right$yhat + sample(right$res, nrow(dta.q.r), T, right$weight),
          left$yhat + sample(left$res, nrow(dta.q.l), T, left$weight))),
    p, q, bw.p, bw.q, N.bc, kernel, bc))
  ci <- quantile(trd.star, c((1 - level)/2, 1 - (1 - level)/2), na.rm = T)
  sd <- sd(trd.star)
  return(list(ci = ci, sd = sd))
}

bw.select <- function(dta, r1, r2, N.bw, kernel, hstart = NULL) {
  # dta:    a data set contains y and x.
  # r1:     the order of polynomial on which MSE is minimized.
  # r2:     the order of polynomial which use initial bandwidth.
  # N.bw:   number of bootstrap to calculate MSE.
  # kernel: kernel.
  # hstart: initial bandwidth for polynomial of order r2.

  # initial bandwidth
  if (is.null(hstart)) {
    if (kernel == "epanechnikov" | kernel == "epa") {
      hstart <- rep(2.345*sd(dta$x)/nrow(dta)^(0.2), 2)
    } else if (kernel == "uniform" | kernel == "uni") {
      hstart <- rep(1.843*sd(dta$x)/nrow(dta)^(0.2), 2)
    } else {
      hstart <- rep(2.576*sd(dta$x)/nrow(dta)^(0.2), 2)
    }
  }

  # polynomial of roder r2
  dta.l       <- dta[dta$x > -hstart[1] & dta$x < 0, ]
  dta.r       <- dta[dta$x < hstart[2] & dta$x >= 0, ]
  left        <- local_polynomial_regression(dta.l, r2, hstart[1], kernel)
  right       <- local_polynomial_regression(dta.r, r2, hstart[2], kernel)

  # choose h to minimize MSE for polynomial of order r1
  boot.l <- replicate(N.bw, left$yhat + sample(left$res, nrow(dta.l), T, left$weight))
  boot.r <- replicate(N.bw, right$yhat + sample(right$res, nrow(dta.r), T, right$weight))

  boot.mse <- function(h, boot.y, x, true.value) {
    mean(apply(boot.y, 2, function(yb) {
      (lm.wfit(cbind(1, poly(x, r1, raw = T)), yb, kweight(x, 0, h, kernel))$coefficients[1] - true.value)^2
    }))
  }

  h.l <- optimize(boot.mse, c(0, hstart[1]), boot.y = boot.l, x = dta.l$x, true.value = left$cons)$minimum
  h.r <- optimize(boot.mse, c(0, hstart[2]), boot.y = boot.r, x = dta.r$x, true.value = right$cons)$minimum

  return(list(initial = hstart, update = c(h.l, h.r)))
}
