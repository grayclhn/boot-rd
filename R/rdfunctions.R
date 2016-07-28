## Copyright (c) 2016, Otávio Bartalotti, Gray Calhoun, and Yang He.
## Available under the MIT "Expat" License, see README.md

kweight <- rdrobust:::rdrobust_kweight

generate.data <- function(model.id) {
  x <- 2*rbeta(500, 2, 4) - 1
  e <- rnorm(500, 0, 0.1295)
  y <- e + switch(model.id,
    ifelse(x < 0, 0.48 +  1.27*x +  7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5,
                  0.52 +  0.84*x -  3.00*x^2 +  7.99*x^3 -  9.01*x^4 + 3.56*x^5),
    ifelse(x < 0, 3.71 +  2.30*x +  3.28*x^2 +  1.45*x^3 +  0.23*x^4 + 0.03*x^5,
                  0.26 + 18.49*x - 54.81*x^2 + 74.30*x^3 - 45.02*x^4 + 9.83*x^5),
    ifelse(x < 0, 0.48 +  1.27*x - 0.5*7.18*x^2 + 0.7*20.21*x^3 + 1.1*21.54*x^4 + 1.5*7.33*x^5,
                  0.52 +  0.84*x - 0.1*3.00*x^2 - 0.3*7.99*x^3 - 0.1*9.01*x^4 + 3.56*x^5))
  return(data.frame(y = y, x = x))}

# Functions that implement the residual and wild bootstraps.

residual_bootstrap <- function(m, y, x) {
  predict(m, data.frame(y = NA, x = x)) +
    sample(resid(m), size = length(y), replace = TRUE, prob = weights(m))
}

wild_values <- c(1 - sqrt(5), 1 + sqrt(5)) / 2
wild_weights <- c(sqrt(5) + 1, sqrt(5) - 1) / (2 * sqrt(5))

wild_bootstrap <- function(m, y, x) {
  yhat <- predict(m, data.frame(y = NA, x = x))
  ehat <- y - yhat
  yhat + (ehat * sample(wild_values, size = length(ehat),
    replace = TRUE, prob = wild_weights))
}

# Simple RD estimator of the treatment effect.

basic_estimator <- function(y, x, w, p = 1,
  m0 = lm(y ~ poly(x, p), subset = x <= 0, weights = w),
  m1 = lm(y ~ poly(x, p), subset = x > 0, weights = w)) {

  cutoff = data.frame(y = NA, x = 0)
  predict(m1, newdata = cutoff) - predict(m0, newdata = cutoff)
}

# The 'boot_estimator' function implements the RD bootstrap procedure
# to construct a de-biased point estimator
#
# Arguments:
# y: dependent variable as two vectors. 0 indicates that the running
#    variable is negative, 1 indicates that it is positive.
# x: running variable
# wp, wq: explicit weights imposed by the kernel near the cutoff
# a:      size of the implicit hypothesis test. The interval will have level 1-a
# nboot:  number of bootstrap replications used for debiasing
# bootfn: a function that generates new values of the depnedent model,
#         given an estimated local polynomial
# p, q: polynomial degree for the estimated model (p1) and the conditional
#         expectation used to implement the bootstrap (p2)
# m0, m1: fitted local regression models below and above the cutoff.

boot_estimator <- function(y, x, wp, wq, p = 1, q = p + 1, nboot, bootfn,
  m0 = lm(y ~ poly(x, q), subset = x <= 0, weights = wq),
  m1 = lm(y ~ poly(x, q), subset = x > 0, weights = wq)) {

  yboot <- rep(NA, length(y))
  i0 <- x <= 0
  i1 <- !i0

  # Generate replications of the naive estimator using the estimated
  # local regression as the bootstrap's true conditional expectation.

  boots <- replicate(nboot, {
    yboot[i1] <- bootfn(m1, y[i1], x[i1])
    yboot[i0] <- bootfn(m0, y[i0], x[i0])
    basic_estimator(yboot, x, wp, p)
  })

  # Subtract the bootstrap estimate of the bias (second and third terms)
  # from the basic estimator.
  basic_estimator(y, x, wp, p) -
    mean(boots) + basic_estimator(y, x, wq, q, m0, m1)
}

# The 'boot_interval' function implements the iterated bootstrap procedure
# to construct a confidence interval. We will use the 'basic bootstrap
# confidence interval' that inverts a bootstrap hypothesis test.
# I.e., if QL and QU are the lower and upper quantiles of the bootstrap
# distribution under the null of mean zero, the interval becomes
# [\tau-hat - QU, \tau-hat - QL]
#
# Arguments:
# y: dependent variable as two vectors. 0 indicates that the running
#    variable is negative, 1 indicates that it is positive.
# x: running variable
# wp, wq: explicit weights imposed by the kernel near the cutoff
# a: size of the implicit hypothesis test. The interval will have level 1-a
# nboot: number of bootstrap  replications for the "outer" bootstrap that
#         generates draws for the confidence interval construction
# nboot2: number of bootstrap replications used for debiasing
# bootfn: a function that generates new values of the depnedent model,
#         given an estimated local polynomial
# p1, p2: polynomial degree for the estimated model (p1) and the conditional
#         expectation used to implement the bootstrap (p2)
# m0, m1: fitted local regression models below and above the cutoff.

boot_interval <- function(y, x, wp, wq, a, p = 1, q = p + 1,
  nboot, bootfn, nboot2 = nboot, type = c("basic", "percentile", "both"),
  m0 = lm(y ~ poly(x, q), subset = x <= 0, weights = wq),
  m1 = lm(y ~ poly(x, q), subset = x > 0, weights = wq)) {

  type <- match.arg(type)

  yboot <- rep(NA, length(y))
  i0 <- x <= 0
  i1 <- !i0

  # Generate replications of the naive estimator using the estimated
  # local regression as the bootstrap's true conditional expectation.

  boots <- replicate(nboot2, {
    yboot[i1] <- bootfn(m1, y[i1], x[i1])
    yboot[i0] <- bootfn(m0, y[i0], x[i0])
    boot_estimator(yboot, x, wp, wq, p, q, nboot, bootfn)
  })

  switch(type,
    basic = boot_ci_basic(boot_estimator(y, x, wp, wq, p, q, nboot2, bootfn),
      quantile(boots, c(1 - a/2, a/2)), basic_estimator(y, x, wq, q)),
    percentile = quantile(boots, c(a/2, 1-a/2)),
    both = rbind(
      basic = boot_ci_basic(boot_estimator(y, x, wp, wq, p, q, nboot2, bootfn),
        quantile(boots, c(1 - a/2, a/2)), basic_estimator(y, x, wq, q)),
      percentile = quantile(boots, c(a/2, 1-a/2))),
    stop("The listed interval type is not supported."))
}

# Lower bound is the bias-corrected estimator (first term) minus the
# (1-a/2) quantile from the bootstrap distribtuion of the bias
# corrected estimator, centered to have mean zero (second and third terms).
# Upper bound is constructed similarly, but with the a/2 quantile.
#
# Arguments
# ---------
# estimate: value of the parameter estimator
# boots: vector of bootstrap draws of the estimator
# boot_parameter: true value of the estimand under the bootstrap distribution
# a: 1 - confidence level

boot_ci_basic <- function(estimate, boots, boot_parameter, a = 0.05) {
 ci <- estimate - quantile(boots, c(1 - a/2, a/2)) + boot_parameter
 names(ci) <- rev(names(ci))
 ci
}

## Simplified wrapper for the bootstrap confidence interval. This function
## is essentially the same as boot_interval except that it calculates the
## kernel weights inside the function using CCT's code.

rdboot <- function(y, x, a, nboot, bootfn, nboot2 = nboot,
  p = 1, q = p + 1, type = c("basic", "percentile", "both"),
  kernel = c("triangular", "uniform", "epanechnikov"),...) {

  type <- match.arg(type)
  kernel <- match.arg(kernel)
  bw <- rdbwselect_2014(y, x, p=p, q=q, kernel=kernel,...)$bws
  wp <- kweight(x, 0, bw[1], kernel)
  wq <- kweight(x, 0, bw[2], kernel)

  estimate <- basic_estimator(y, x, wp, p)
  names(estimate) <- "estimate"

  ci <- boot_interval(y, x, wp, wq, a, p, q, nboot, bootfn, nboot2, type)

  if (!is.matrix(ci)) {
    return(c(estimate, ci))
  } else {
    return(cbind(estimate, ci))
  }
}
