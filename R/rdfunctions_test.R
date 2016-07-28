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
  cct <- rdrobust(dta$y, dta$x, kernel = "uni")
  bws <- cct$bws

  ## Due to some differences in implemenation details, the two approaches
  ## only give identical results when the dataset and parameters are very
  ## constrained:
  ## - The same bandwidth must be used for the linear and quadratic models
  ## - All of the xs must fall in that bandwidth
  ## - All of the positive xs must be first, then the negative xs.
  ## Without these constraints, the two methods are essentially the same
  ## but draw pseudo random variables from the RNG in different orders so
  ## the results are not identical.
  smalldta <- dta[abs(dta$x) < bws[2],]
  smalldta <- smalldta[order(smalldta$x, decreasing = TRUE),]

  ## Original version of the point estimate
  set.seed(common.seed)
  old_estimate <- rd.estimate(smalldta, 1, 2, bws[2], bws[2], 17, "tri", T)
  old_interval <- rd.ci(smalldta, 1, 2, bws[2], bws[2], 17, 19, 0.9, "tri", T)

  ## New version
  set.seed(common.seed)
  new_estimate <- with(smalldta, {
    wq <- kweight(x, 0, bws[2], "tri")
    boot_estimator(y, x, wq, wq, p=1, q=2, nboot = 17,
      bootfn = residual_bootstrap)})
  new_interval <- with(smalldta, {
    wq <- kweight(x, 0, bws[2], "tri")
    boot_interval(y, x, wq, wq, a = 0.1, p = 1, q = 2,
      17, residual_bootstrap, 19, type = "percentile")})

  expect_equivalent(new_estimate, old_estimate)
  expect_equivalent(new_interval, old_interval$ci)})
