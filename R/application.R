## Copyright (c) 2016, Otávio Bartalotti, Gray Calhoun, and Yang He.
## Available under the MIT "Expat" License, see README.md

library(rdrobust)
library(xtable)
source("RDfunctions.R")

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
  cct <- rdrobust(dta$y, dta$x, kernel = "uni")
  bws <- cct$bws
  
  t1  <- cct$coef[3]
  ci1 <- cct$ci[3, ]
  
  t2  <- rd.estimate(dta, 1, 2, bws[1], bws[2], N.bc, "uni", T)
  ci2 <- rd.ci(dta, 1, 2, bws[1], bws[2], N.bc, N.ci, level, "uni", T)$ci
  
  return(matrix(c(t1, ci1, NA, bws,
                  t2, ci2, NA, bws), nrow = 2, byrow = T))
}

## Table 1: spending per child

spending <- matrix(0, nrow = 10, ncol = 6)
colnames(spending) <- c("Effect", "Lower CI", "Upper CI", "P-value", "h", "b")

spending[1, ] <- original(hsspend_per_kid_68, 9)
spending[2, ] <- original(hsspend_per_kid_68, 18)
spending[3, ] <- original(hsspend_per_kid_68, 36)
spending[4:5, ] <- robust(hsspend_per_kid_68)

spending[6, ] <- original(hsspend_per_kid_72, 9)
spending[7, ] <- original(hsspend_per_kid_72, 18)
spending[8, ] <- original(hsspend_per_kid_72, 36)
spending[9:10, ] <- robust(hsspend_per_kid_72)

## Table 2: Mortality

mortality <- matrix(0, nrow = 5, ncol = 6)
colnames(mortality) <- c("Effect", "Lower CI", "Upper CI", "P-value", "h", "b")

mortality[1, ] <- original(age5_9_sum2, 9)
mortality[2, ] <- original(age5_9_sum2, 18)
mortality[3, ] <- original(age5_9_sum2, 36)
mortality[4:5, ] <- robust(age5_9_sum2)

## Table 3: Education

education <- matrix(0, nrow = 6, ncol = 6)
colnames(education) <- c("Effect", "Lower CI", "Upper CI", "P-value", "h", "b")

education[1, ] <- original(hsplus18_24, 7)
education[2:3, ] <- robust(hsplus18_24)

education[4, ] <- original(some_clg18_24, 7)
education[5:6, ] <- robust(some_clg18_24)

print(xtable(spending, digits = c(0,1,1,1,0,3,3)), include.rownames = F, type = "latex", file = "output/spending.tex")
print(xtable(mortality, digits = c(0,3,3,3,0,3,3)), include.rownames = F, type = "latex", file = "output/mortality.tex")
print(xtable(education, digits = c(0,3,3,3,0,3,3)), include.rownames = F, type = "latex", file = "output/education.tex")










