rm(list = ls())
library(rdrobust)
library(xtable)
try(source("slackrISE.R"), silent = T)

source("functions.R")

setwd("C:/Users/Yang/OneDrive/QJE data and code PUBLIC/analysis data")

###################################################################################
## slightly modify the function rd.ci so that it exports standard error
###################################################################################

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
  left  <- lpreg(dta.q.l, q, bw.q[1], kernel)
  right <- lpreg(dta.q.r, q, bw.q[2], kernel)
  
  # bootstrap
  trd.star <- replicate(N.ci, rd.estimate(data.frame(x = c(dta.q.l$x, dta.q.r$x),
                                                     y = c(left$yhat + sample(left$res, nrow(dta.q.l), T, left$weight),
                                                           right$yhat + sample(right$res, nrow(dta.q.r), T, right$weight))),
                                          p, q, bw.p, bw.q, N.bc, kernel, bc))
  ci <- quantile(trd.star, c((1 - level)/2, 1 - (1 - level)/2), na.rm = T)
  sd <- sd(trd.star)
  return(list(ci = ci, sd = sd))
}

###################################################################################
## replicate part of their results ################################################
###################################################################################

cut <- 59.1984
census2 <- read.csv("census2.csv")
census3 <- read.csv("census3.csv")
mort <- read.csv("mort.csv")
mort4 <- read.csv("mort4.csv")
census_1990 <- read.csv("census_1990.csv")
census_2000 <- read.csv("census_2000.csv")

#### table 2 ####

## dependent variable: hsspend_per_kid_68

hsspend_per_kid_68 <- census2[ , c("hsspend_per_kid_68", "povrate60")]
hsspend_per_kid_68 <- hsspend_per_kid_68[complete.cases(hsspend_per_kid_68), ]
colnames(hsspend_per_kid_68) <- c("y", "x")
hsspend_per_kid_68$x <- hsspend_per_kid_68$x - cut

test_that("table 2 hsspend_per_kid_68",{
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_68, 1, 1, 9, 9, 0, "uni", F)), digits = 3), 137.251)    # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_68, 1, 1, 18, 18, 0, "uni", F)), digits = 3), 114.711)  # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_68, 1, 1, 36, 36, 0, "uni", F)), digits = 3), 134.491)  # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_68, 1, 1, 8, 8, 0, "uni", F)), digits = 3), 127.389)    # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_68, 2, 2, 16, 16, 0, "uni", F)), digits = 3), 112.706)  # Parametric quadratic h = 16
})


## dependent variable: hsspend_per_kid_72

hsspend_per_kid_72 <- census2[ , c("hsspend_per_kid_72", "povrate60")]
hsspend_per_kid_72 <- hsspend_per_kid_72[complete.cases(hsspend_per_kid_72), ]
colnames(hsspend_per_kid_72) <- c("y", "x")
hsspend_per_kid_72$x <- hsspend_per_kid_72$x - cut

test_that("table 2 hsspend_per_kid_72",{
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_72, 1, 1, 9, 9, 0, "uni", F)), digits = 3), 182.119)    # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_72, 1, 1, 18, 18, 0, "uni", F)), digits = 3), 88.959)   # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_72, 1, 1, 36, 36, 0, "uni", F)), digits = 3), 130.153)  # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_72, 1, 1, 8, 8, 0, "uni", F)), digits = 3), 175.773)    # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(hsspend_per_kid_72, 2, 2, 16, 16, 0, "uni", F)), digits = 3), 155.899)  # Parametric quadratic h = 16
})


## dependent variable: socspend_per_cap72

socspend_per_cap72 <- census2[ , c("socspend_per_cap72", "povrate60")]
socspend_per_cap72 <- socspend_per_cap72[complete.cases(socspend_per_cap72), ]
colnames(socspend_per_cap72) <- c("y", "x")
socspend_per_cap72$x <- socspend_per_cap72$x - cut

test_that("table 2 socspend_per_cap72",{
  expect_equal(round(as.numeric(rd.estimate(socspend_per_cap72, 1, 1, 9, 9, 0, "uni", F)), digits = 3), 14.474)     # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(socspend_per_cap72, 1, 1, 18, 18, 0, "uni", F)), digits = 3), 19.590)   # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(socspend_per_cap72, 1, 1, 36, 36, 0, "uni", F)), digits = 3), 14.506)   # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(socspend_per_cap72, 1, 1, 8, 8, 0, "uni", F)), digits = 3), -2.361)     # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(socspend_per_cap72, 2, 2, 16, 16, 0, "uni", F)), digits = 3), 5.659)    # Parametric quadratic h = 16
})

#### table 3 ####

## dependent variable: age5_9_sum2

age5_9_sum2 <- census3[ , c("age5_9_sum2", "povrate60")]
age5_9_sum2 <- age5_9_sum2[complete.cases(age5_9_sum2), ]
colnames(age5_9_sum2) <- c("y", "x")
age5_9_sum2$x <- age5_9_sum2$x - cut

test_that("table 2 age5_9_sum2",{
  expect_equal(round(as.numeric(rd.estimate(age5_9_sum2, 1, 1, 9, 9, 0, "uni", F)), digits = 3), -1.895)     # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(age5_9_sum2, 1, 1, 18, 18, 0, "uni", F)), digits = 3), -1.198)   # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(age5_9_sum2, 1, 1, 36, 36, 0, "uni", F)), digits = 3), -1.114)   # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(age5_9_sum2, 1, 1, 8, 8, 0, "uni", F)), digits = 3), -2.201)     # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(age5_9_sum2, 2, 2, 16, 16, 0, "uni", F)), digits = 3), -2.558)    # Parametric quadratic h = 16
})

## dependent variable: age5_9_injury_rate

age5_9_injury_rate <- census3[ , c("age5_9_injury_rate", "povrate60")]
age5_9_injury_rate <- age5_9_injury_rate[complete.cases(age5_9_injury_rate), ]
colnames(age5_9_injury_rate) <- c("y", "x")
age5_9_injury_rate$x <- age5_9_injury_rate$x - cut

test_that("table 2 age5_9_injury_rate",{
  expect_equal(round(as.numeric(rd.estimate(age5_9_injury_rate, 1, 1, 9, 9, 0, "uni", F)), digits = 3), 0.195)      # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(age5_9_injury_rate, 1, 1, 18, 18, 0, "uni", F)), digits = 3), 2.426)    # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(age5_9_injury_rate, 1, 1, 36, 36, 0, "uni", F)), digits = 3), 0.679)    # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(age5_9_injury_rate, 1, 1, 8, 8, 0, "uni", F)), digits = 3), -0.164)     # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(age5_9_injury_rate, 2, 2, 16, 16, 0, "uni", F)), digits = 3), 0.775)    # Parametric quadratic h = 16
})

## dependent variable: age5_9_rate

age5_9_rate <- census3[ , c("age5_9_rate", "povrate60")]
age5_9_rate <- age5_9_rate[complete.cases(age5_9_rate), ]
colnames(age5_9_rate) <- c("y", "x")
age5_9_rate$x <- age5_9_rate$x - cut

test_that("table 2 age5_9_rate",{
  expect_equal(round(as.numeric(rd.estimate(age5_9_rate, 1, 1, 9, 9, 0, "uni", F)), digits = 3), -3.416)      # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(age5_9_rate, 1, 1, 18, 18, 0, "uni", F)), digits = 3), 0.053)     # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(age5_9_rate, 1, 1, 36, 36, 0, "uni", F)), digits = 3), -1.537)    # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(age5_9_rate, 1, 1, 8, 8, 0, "uni", F)), digits = 3), -3.896)      # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(age5_9_rate, 2, 2, 16, 16, 0, "uni", F)), digits = 3), -2.927)    # Parametric quadratic h = 16
})

## dependent variable: age25plus_sum2

age25plus_sum2 <- census3[ , c("age25plus_sum2", "povrate60")]
age25plus_sum2 <- age25plus_sum2[complete.cases(age25plus_sum2), ]
colnames(age25plus_sum2) <- c("y", "x")
age25plus_sum2$x <- age25plus_sum2$x - cut

test_that("table 2 age25plus_sum2",{
  expect_equal(round(as.numeric(rd.estimate(age25plus_sum2, 1, 1, 9, 9, 0, "uni", F)), digits = 3), 2.204)      # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(age25plus_sum2, 1, 1, 18, 18, 0, "uni", F)), digits = 3), 6.016)    # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(age25plus_sum2, 1, 1, 36, 36, 0, "uni", F)), digits = 3), 5.872)    # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(age25plus_sum2, 1, 1, 8, 8, 0, "uni", F)), digits = 3), 2.091)      # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(age25plus_sum2, 2, 2, 16, 16, 0, "uni", F)), digits = 3), 2.574)    # Parametric quadratic h = 16
})

## dependent variable: age25plus_injury_rate

age25plus_injury_rate <- mort[ , c("age25plus_injury_rate", "povrate60")]
age25plus_injury_rate <- age25plus_injury_rate[complete.cases(age25plus_injury_rate), ]
colnames(age25plus_injury_rate) <- c("y", "x")
age25plus_injury_rate$x <- age25plus_injury_rate$x - cut

test_that("table 2 age25plus_injury_rate",{
  expect_equal(round(as.numeric(rd.estimate(age25plus_injury_rate, 1, 1, 9, 9, 0, "uni", F)), digits = 3), 5.697)       # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(age25plus_injury_rate, 1, 1, 18, 18, 0, "uni", F)), digits = 3), 7.276)     # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(age25plus_injury_rate, 1, 1, 36, 36, 0, "uni", F)), digits = 3), 4.398)     # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(age25plus_injury_rate, 1, 1, 8, 8, 0, "uni", F)), digits = 3), 2.65)        # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(age25plus_injury_rate, 2, 2, 16, 16, 0, "uni", F)), digits = 3), 4.276)     # Parametric quadratic h = 16
})

## dependent variable: rate_5964

rate_5964 <- mort4[ , c("rate_5964", "povrate60")]
rate_5964 <- rate_5964[complete.cases(rate_5964), ]
colnames(rate_5964) <- c("y", "x")
rate_5964$x <- rate_5964$x - cut

test_that("table 2 rate_5964",{
  expect_equal(round(as.numeric(rd.estimate(rate_5964, 1, 1, 9, 9, 0, "uni", F)), digits = 3), -3.327)      # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(rate_5964, 1, 1, 18, 18, 0, "uni", F)), digits = 3), -1.076)    # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(rate_5964, 1, 1, 36, 36, 0, "uni", F)), digits = 3), -0.066)    # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(rate_5964, 1, 1, 8, 8, 0, "uni", F)), digits = 3), -3.754)      # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(rate_5964, 2, 2, 16, 16, 0, "uni", F)), digits = 3), -4.869)    # Parametric quadratic h = 16
})

## dependent variable: white_age5_9_sum2

white_age5_9_sum2 <- census3[ , c("white_age5_9_sum2", "povrate60")]
white_age5_9_sum2 <- white_age5_9_sum2[complete.cases(white_age5_9_sum2), ]
colnames(white_age5_9_sum2) <- c("y", "x")
white_age5_9_sum2$x <- white_age5_9_sum2$x - cut

test_that("table 2 white_age5_9_sum2",{
  expect_equal(round(as.numeric(rd.estimate(white_age5_9_sum2, 1, 1, 9, 9, 0, "uni", F)), digits = 3), -1.105)      # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(white_age5_9_sum2, 1, 1, 18, 18, 0, "uni", F)), digits = 3), -0.865)    # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(white_age5_9_sum2, 1, 1, 36, 36, 0, "uni", F)), digits = 3), -0.749)    # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(white_age5_9_sum2, 1, 1, 8, 8, 0, "uni", F)), digits = 3), -1.334)      # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(white_age5_9_sum2, 2, 2, 16, 16, 0, "uni", F)), digits = 3), -1.746)    # Parametric quadratic h = 16
})

## dependent variable: black_age5_9_sum2

black_age5_9_sum2 <- census3[ , c("black_age5_9_sum2", "povrate60")]
black_age5_9_sum2 <- black_age5_9_sum2[complete.cases(black_age5_9_sum2), ]
colnames(black_age5_9_sum2) <- c("y", "x")
black_age5_9_sum2$x <- black_age5_9_sum2$x - cut

test_that("table 2 black_age5_9_sum2",{
  expect_equal(round(as.numeric(rd.estimate(black_age5_9_sum2, 1, 1, 9, 9, 0, "uni", F)), digits = 3), -2.275)      # Nonparametric h = 9
  expect_equal(round(as.numeric(rd.estimate(black_age5_9_sum2, 1, 1, 18, 18, 0, "uni", F)), digits = 3), -2.719)    # Nonparametric h = 18
  expect_equal(round(as.numeric(rd.estimate(black_age5_9_sum2, 1, 1, 36, 36, 0, "uni", F)), digits = 3), -1.589)    # Nonparametric h = 36
  expect_equal(round(as.numeric(rd.estimate(black_age5_9_sum2, 1, 1, 8, 8, 0, "uni", F)), digits = 3), -1.699)      # Parametric linear h = 8
  expect_equal(round(as.numeric(rd.estimate(black_age5_9_sum2, 2, 2, 16, 16, 0, "uni", F)), digits = 3), -1.93)     # Parametric quadratic h = 16
})

#### table 4 ####

## dependent variable: hsplus18_24

hsplus18_24 <- census_1990[ , c("hsplus18_24", "povrate60")]
hsplus18_24 <- hsplus18_24[complete.cases(hsplus18_24), ]
colnames(hsplus18_24) <- c("y", "x")
hsplus18_24$x <- hsplus18_24$x - cut

test_that("table 2 hsplus18_24",{
  expect_equal(round(as.numeric(rd.estimate(hsplus18_24, 1, 1, 7, 7, 0, "uni", F)), digits = 3), 0.030)    # Nonparametric h = 7
  expect_equal(round(as.numeric(rd.estimate(hsplus18_24, 1, 1, 4, 4, 0, "uni", F)), digits = 3), 0.043)    # Parametric linear h = 4
})

## dependent variable: some_clg18_24

some_clg18_24 <- census_1990[ , c("some_clg18_24", "povrate60")]
some_clg18_24 <- some_clg18_24[complete.cases(some_clg18_24), ]
colnames(some_clg18_24) <- c("y", "x")
some_clg18_24$x <- some_clg18_24$x - cut

test_that("table 2 some_clg18_24",{
  expect_equal(round(as.numeric(rd.estimate(some_clg18_24, 1, 1, 7, 7, 0, "uni", F)), digits = 3), 0.037)    # Nonparametric h = 7
  expect_equal(round(as.numeric(rd.estimate(some_clg18_24, 1, 1, 4, 4, 0, "uni", F)), digits = 3), 0.051)    # Parametric linear h = 4
})

## dependent variable: pcthsplus_18_24_00

pcthsplus_18_24_00 <- census_2000[ , c("pcthsplus_18_24_00", "povrate60")]
pcthsplus_18_24_00 <- pcthsplus_18_24_00[complete.cases(pcthsplus_18_24_00), ]
colnames(pcthsplus_18_24_00) <- c("y", "x")
pcthsplus_18_24_00$x <- pcthsplus_18_24_00$x - cut

test_that("table 2 pcthsplus_18_24_00",{
  expect_equal(round(as.numeric(rd.estimate(pcthsplus_18_24_00, 1, 1, 7, 7, 0, "uni", F)), digits = 3), 0.000)    # Nonparametric h = 7
  expect_equal(round(as.numeric(rd.estimate(pcthsplus_18_24_00, 1, 1, 4, 4, 0, "uni", F)), digits = 3), 0.024)    # Parametric linear h = 4
})

## dependent variable: pctsomecollege_18_24_00

pctsomecollege_18_24_00 <- census_2000[ , c("pctsomecollege_18_24_00", "povrate60")]
pctsomecollege_18_24_00 <- pctsomecollege_18_24_00[complete.cases(pctsomecollege_18_24_00), ]
colnames(pctsomecollege_18_24_00) <- c("y", "x")
pctsomecollege_18_24_00$x <- pctsomecollege_18_24_00$x - cut

test_that("table 2 pctsomecollege_18_24_00",{
  expect_equal(round(as.numeric(rd.estimate(pctsomecollege_18_24_00, 1, 1, 7, 7, 0, "uni", F)), digits = 3), 0.028)    # Nonparametric h = 7
  expect_equal(round(as.numeric(rd.estimate(pctsomecollege_18_24_00, 1, 1, 4, 4, 0, "uni", F)), digits = 3), 0.042)    # Parametric linear h = 4
})

## dependent variable: pcthsplus_25_34_00

pcthsplus_25_34_00 <- census_2000[ , c("pcthsplus_25_34_00", "povrate60")]
pcthsplus_25_34_00 <- pcthsplus_25_34_00[complete.cases(pcthsplus_25_34_00), ]
colnames(pcthsplus_25_34_00) <- c("y", "x")
pcthsplus_25_34_00$x <- pcthsplus_25_34_00$x - cut

test_that("table 2 pcthsplus_25_34_00",{
  expect_equal(round(as.numeric(rd.estimate(pcthsplus_25_34_00, 1, 1, 7, 7, 0, "uni", F)), digits = 3), 0.006)    # Nonparametric h = 7
  expect_equal(round(as.numeric(rd.estimate(pcthsplus_25_34_00, 1, 1, 4, 4, 0, "uni", F)), digits = 3), 0.015)    # Parametric linear h = 4
})

## dependent variable: pctsomecollege_25_34_00

pctsomecollege_25_34_00 <- census_2000[ , c("pctsomecollege_25_34_00", "povrate60")]
pctsomecollege_25_34_00 <- pctsomecollege_25_34_00[complete.cases(pctsomecollege_25_34_00), ]
colnames(pctsomecollege_25_34_00) <- c("y", "x")
pctsomecollege_25_34_00$x <- pctsomecollege_25_34_00$x - cut

test_that("table 2 pctsomecollege_25_34_00",{
  expect_equal(round(as.numeric(rd.estimate(pctsomecollege_25_34_00, 1, 1, 7, 7, 0, "uni", F)), digits = 3), 0.040)    # Nonparametric h = 7
  expect_equal(round(as.numeric(rd.estimate(pctsomecollege_25_34_00, 1, 1, 4, 4, 0, "uni", F)), digits = 3), 0.043)    # Parametric linear h = 4
})

###################################################################################
## compare with results from CCT and bootstrap ####################################
###################################################################################
N.bc <- 500
N.ci <- 999
level <- 0.95
set.seed(555)

robustestimate <- function(dta) {
  bws <- rdbwselect(dta$y, dta$x, kernel = "uni")$bws
  boots.estimate <- rd.estimate(dta, 1, 2, bws[1], bws[2], N.bc, "uni", T)
  boots.ci <- rd.ci(dta, 1, 2, bws[1], bws[2], N.bc, N.ci, level, "uni", T)
  cct <- rdrobust(dta$y, dta$x, kernel = "uni")
  return(list(bws = bws, 
              b.e = boots.estimate,
              b.c = boots.ci$ci,
              b.s = boots.ci$sd,
              c.e = cct$coef[3],
              c.c = cct$ci[3, ],
              c.s = cct$se[3]))
}

rd.hsspend_per_kid_68 <- robustestimate(hsspend_per_kid_68)
rd.hsspend_per_kid_72 <- robustestimate(hsspend_per_kid_72)
rd.socspend_per_cap72 <- robustestimate(socspend_per_cap72)
rd.age5_9_sum2 <- robustestimate(age5_9_sum2)
rd.hsplus18_24 <- robustestimate(hsplus18_24)
rd.some_clg18_24 <- robustestimate(some_clg18_24)
rd.pcthsplus_18_24_00 <- robustestimate(pcthsplus_18_24_00)
rd.pctsomecollege_18_24_00 <- robustestimate(pctsomecollege_18_24_00)
rd.pcthsplus_25_34_00 <- robustestimate(pcthsplus_25_34_00)
rd.pctsomecollege_25_34_00 <- robustestimate(pctsomecollege_25_34_00)







