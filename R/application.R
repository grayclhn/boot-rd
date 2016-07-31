## Copyright (c) 2016, Otávio Bartalotti, Gray Calhoun, and Yang He.
## Available under the MIT "Expat" License, see README.md

library(rdrobust)
library(xtable)
source("rdfunctions.R")

## Load Ludwig and Miller's (2007) dataset

cutoff <- 59.1984
census2 <- read.csv("./LM2007data/census2.csv")
census3 <- read.csv("./LM2007data/census3.csv")
mort <- read.csv("./LM2007data/mort.csv")
mort4 <- read.csv("./LM2007data/mort4.csv")
census_1990 <- read.csv("./LM2007data/census_1990.csv")
census_2000 <- read.csv("./LM2007data/census_2000.csv")

## Set parameterizations and functions for the analysis

N.bc <- 500
N.ci <- 999
level <- 0.95
set.seed(555)

original <- function(dta, bw) {
  cct <- rdrobust(dta$y, dta$x, h = bw, b = bw, kernel = "uni")
  t  <- cct$coef[1]
  ci <- cct$ci[1, ]
  return(c(t, ci, NA, bw, NA))
}

robust <- function(dta) {
  bws <- rdbwselect_2014(dta$y, dta$x, kernel = "uni")$bws
  cct <- rdrobust(dta$y, dta$x, kernel = "uni", h = bws[1], b = bws[2])
  t1  <- cct$coef[3]
  ci1 <- cct$ci[3, ]

  boot <- rdboot(dta$y, dta$x, 1 - level, N.bc, wild_bootstrap, N.ci,
    type = "basic", kernel = "uniform")
  t2 <- boot[1]
  ci2 <- boot[2:3]
  return(matrix(c(t1, ci1, NA, bws,
                  t2, ci2, NA, bws), nrow = 2, byrow = T))
}

## Table 1: spending per child

spending <- matrix(0, nrow = 10, ncol = 6)
colnames(spending) <- c("Effect", "Lower CI", "Upper CI", "P-value", "h", "b")

hsspend_per_kid_68 <- census2[ , c("hsspend_per_kid_68", "povrate60")]
hsspend_per_kid_68 <- hsspend_per_kid_68[complete.cases(hsspend_per_kid_68), ]
colnames(hsspend_per_kid_68) <- c("y", "x")
hsspend_per_kid_68$x <- hsspend_per_kid_68$x - cutoff

spending[1, ] <- original(hsspend_per_kid_68, 9)
spending[2, ] <- original(hsspend_per_kid_68, 18)
spending[3, ] <- original(hsspend_per_kid_68, 36)
spending[4:5, ] <- robust(hsspend_per_kid_68)

hsspend_per_kid_72 <- census2[ , c("hsspend_per_kid_72", "povrate60")]
hsspend_per_kid_72 <- hsspend_per_kid_72[complete.cases(hsspend_per_kid_72), ]
colnames(hsspend_per_kid_72) <- c("y", "x")
hsspend_per_kid_72$x <- hsspend_per_kid_72$x - cutoff

spending[6, ] <- original(hsspend_per_kid_72, 9)
spending[7, ] <- original(hsspend_per_kid_72, 18)
spending[8, ] <- original(hsspend_per_kid_72, 36)
spending[9:10, ] <- robust(hsspend_per_kid_72)

## Table 2: Mortality

mortality <- matrix(0, nrow = 5, ncol = 6)
colnames(mortality) <- c("Effect", "Lower CI", "Upper CI", "P-value", "h", "b")

age5_9_sum2 <- census3[ , c("age5_9_sum2", "povrate60")]
age5_9_sum2 <- age5_9_sum2[complete.cases(age5_9_sum2), ]
colnames(age5_9_sum2) <- c("y", "x")
age5_9_sum2$x <- age5_9_sum2$x - cutoff

mortality[1, ] <- original(age5_9_sum2, 9)
mortality[2, ] <- original(age5_9_sum2, 18)
mortality[3, ] <- original(age5_9_sum2, 36)
mortality[4:5, ] <- robust(age5_9_sum2)

## Table 3: Education

education <- matrix(0, nrow = 6, ncol = 6)
colnames(education) <- c("Effect", "Lower CI", "Upper CI", "P-value", "h", "b")

hsplus18_24 <- census_1990[ , c("hsplus18_24", "povrate60")]
hsplus18_24 <- hsplus18_24[complete.cases(hsplus18_24), ]
colnames(hsplus18_24) <- c("y", "x")
hsplus18_24$x <- hsplus18_24$x - cutoff

education[1, ] <- original(hsplus18_24, 7)
education[2:3, ] <- robust(hsplus18_24)

some_clg18_24 <- census_1990[ , c("some_clg18_24", "povrate60")]
some_clg18_24 <- some_clg18_24[complete.cases(some_clg18_24), ]
colnames(some_clg18_24) <- c("y", "x")
some_clg18_24$x <- some_clg18_24$x - cutoff

education[4, ] <- original(some_clg18_24, 7)
education[5:6, ] <- robust(some_clg18_24)

print(xtable(spending, digits = c(0,1,1,1,0,3,3)), include.rownames = F, type = "latex", file = "output/spending.tex")
print(xtable(mortality, digits = c(0,3,3,3,0,3,3)), include.rownames = F, type = "latex", file = "output/mortality.tex")
print(xtable(education, digits = c(0,3,3,3,0,3,3)), include.rownames = F, type = "latex", file = "output/education.tex")










