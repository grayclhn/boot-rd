# Verify that new versions of the functions (in 'rdfunctions.R') return
# the same values as our original versions (in 'rdfunctions_old.R')

library(rdrobust)
library(testthat)

source("rdfunctions.R")
source("rdfunctions_old.R")

test_that("Bias-corrected point estimates give the same result as our original version", {
  dta <- generate.data(1)
  common.seed <- 32

  ## Original approach
  ##
  ## Due to some differences in implemenation details, the two approaches
  ## only give identical results when the dataset and parameters are very
  ## constrained:
  ## - The same bandwidth must be used for the linear and quadratic models
  ## - All of the xs must fall in that bandwidth
  ## - All of the positive xs must be first, then the negative xs.
  ## Without these constraints, the two methods are essentially the same
  ## but draw pseudo random variables from the RNG in different orders so
  ## the results are not identical.
  cct <- rdrobust(dta$y, dta$x, kernel = "tri")
  bw <- rdbwselect_2014(dta$y, dta$x)$bws[2]

  smalldta <- dta[abs(dta$x) < bw,]
  smalldta <- smalldta[order(smalldta$x, decreasing = TRUE),]

  ## Original version of the point estimate
  set.seed(common.seed)
  old_estimate <- rd.estimate(smalldta, 1, 2, bw, bw, 17, "tri", T)
  old_interval <- rd.ci(smalldta, 1, 2, bw, bw, 17, 19, 0.9, "tri", T)

  ## New version
  set.seed(common.seed)
  new_estimate <- with(smalldta, {
    wq <- kweight(x, 0, bw, "tri")
    boot_estimator(y, x, wq, wq, p=1, q=2, nboot = 17,
      bootfn = residual_bootstrap)})
  new_interval <- with(smalldta, {
    wq <- kweight(x, 0, bw, "tri")
    boot_interval(y, x, wq, wq, a = 0.1, p = 1, q = 2,
      17, residual_bootstrap, 19, type = "percentile")})

  expect_equivalent(new_estimate, old_estimate)
  expect_equivalent(new_interval, old_interval$ci)})

test_that("Convenient version of interval constructions runs and matches explicit version.", {
  dta <- generate.data(1)
  common.seed <- 341
  set.seed(common.seed)
  basic_interval <- with(dta, {
    bw <- rdbwselect_2014(y, x)$bws
    wp <- kweight(x, 0, bw[1], "tri")
    wq <- kweight(x, 0, bw[2], "tri")
    estimate <- boot_estimator(y, x, wp, wq, 1, 2, 17, residual_bootstrap)
    boot_interval(y, x, wp, wq, a = 0.1, p = 1, q = 2,
      17, residual_bootstrap, 19, type = "percentile")})
  set.seed(common.seed)
  nice_interval <- rdboot(dta$y, dta$x, 0.1, 17, residual_bootstrap, 19,
    type = "percentile", kernel = "triangular")
  expect_equivalent(basic_interval, nice_interval[2:3])})

# Minimal testing for Wild bootstrap so far

"%isin%" <- function(x, y)
  all(sapply(x, function(xval)
    any(sapply(y, function(yval)
      isTRUE(all.equal(xval, yval, check.attributes = FALSE))))))

test_that("We didn't screw up the weights for the wild bootstrap", {
  expect_equal(sum(wild_weights), 1)
  expect_equal(sum(wild_weights * wild_values), 0)

  ## basic check that the wild bootstrap 'works'
  d <- data.frame(x = rep(c(0, 1), 4), y = rep(c(0, 1), each = 4))
  w <- wild_bootstrap(lm(y ~ x, data = d), d$y, d$x)
  expect_true(w %isin% c(.5 * (1 - wild_values), .5 * (1 + wild_values)))})
