## Copyright (c) 2016, Otávio Bartalotti, Gray Calhoun, and Yang He.
## Available under the MIT "Expat" License, see README.md

library(rdrobust)
library(foreach)
library(doParallel)
library(doRNG)
library(xtable)
try(source("slackrISE.R"), T)

source("rdfunctions_old.R")
source("rdfunctions.R")
source("rdfunctions_ForMc.R")

N.ci <- 999
a <- 0.05
N.simu <- 5000
N.core <- 35
N.bc <- 500

## naive bootstrap estimator
naiveboot <- function(y, x, kernel, Nci = 999) {
  
  ## point estimate
  h <- rdbwselect_2014(y, x, kernel=kernel, bwselect = "IK")$bws[1]
  tau <- rd.estimate(data.frame(y, x), 1, 2, h, h, 500, kernel, F)

  ## bootstrap
  Btau <- replicate(Nci, {
    i <- sample.int(length(y), length(y), replace = T)
    h <- rdbwselect_2014(y[i], x[i], kernel=kernel, bwselect = "IK")$bws[1]
    rd.estimate(data.frame(y[i], x[i]), 1, 2, h, h, 500, kernel, F)
  })
  
  return(c(tau, quantile(Btau, c(0.025, 0.975))))
}


## simulation 1: coverage of CI from bootstrap (basic CI)

simulation.1 <- function(model.id, kernel, heteroskedasticity = F) {
  
  if (heteroskedasticity == T) {
    avebw <- average_bw(model.id, kernel)
    gen_data <- function() generate.data.h(model.id, avebw)
  } else {
    gen_data <- function() generate.data(model.id)
  }
  
  simu <- function() {
    dta <- gen_data()
    rdboot_ForMc(dta$y, dta$x, a, N.bc, N.ci, 1, 2, "wild", kernel)[1, ]
  }
  
  cl <- makeCluster(N.core)
  registerDoParallel(cl)
  export.obj <- c("N.ci", "N.bc", "a", "generate.data", "kweight", "rdboot_ForMc",
                  "wild_values", "wild_weights", "boot_estimator_ForMc", "boot_dist_ForMc",
                  "boot_ci_basic", "generate.data.h")
  collect.simu <- foreach(i=1:N.simu, .combine="rbind", .packages="rdrobust", .export=export.obj, .inorder=F) %dorng%
    simu()
  stopCluster(cl)
  
  t.true    <- ifelse(model.id ==2, -3.45, 0.04)
  bias      <- t.true - mean(collect.simu[ , 1])
  SD        <- sd(collect.simu[ , 1])
  MSE       <- sqrt(mean((collect.simu[ , 1] - t.true)^2))
  coverage  <- mean(collect.simu[ , 2] <= t.true & collect.simu[ , 3] >= t.true)
  length    <- mean(collect.simu[ , 3] - collect.simu[ , 2])
  
  return(c(true = t.true, bias = bias, SD = SD, MSE = MSE, coverage = coverage, length = length))
}

## simulation 2: coverage of CI from CCT(2004)

simulation.2 <- function(model.id, kernel, heteroskedasticity = F) {
  
  if (heteroskedasticity == T) {
    avebw <- average_bw(model.id, kernel)
    gen_data <- function() generate.data.h(model.id, avebw)
  } else {
    gen_data <- function() generate.data(model.id)
  }
  
  simu <- function() {
    dta <- gen_data()
    bws <- rdbwselect_2014(dta$y, dta$x, kernel = kernel)$bws
    results <- rdrobust(dta$y, dta$x, kernel = kernel,
                        h = bws[1], b = bws[2], level = 100 * (1-a), vce = "hc3")
    t <- results$coef[3]
    ci <- results$ci[3, ]
    return(c(t, ci))
  }
  
  cl <- makeCluster(N.core)
  registerDoParallel(cl)
  export.obj <- c("generate.data", "generate.data.h", "a")
  collect.simu <- foreach(i=1:N.simu, .combine="rbind", .packages="rdrobust", .export=export.obj, .inorder=F) %dorng%
    simu()
  stopCluster(cl)
  
  t.true    <- ifelse(model.id ==2, -3.45, 0.04)
  bias      <- t.true - mean(collect.simu[ , 1])
  SD        <- sd(collect.simu[ , 1])
  MSE       <- sqrt(mean((collect.simu[ , 1] - t.true)^2))
  coverage  <- mean(collect.simu[ , 2] <= t.true & collect.simu[ , 3] >= t.true)
  length    <- mean(collect.simu[ , 3] - collect.simu[ , 2])
  
  return(c(true = t.true, bias = bias, SD = SD, MSE = MSE, coverage = coverage, length = length))
}

## simulation 3: coverage of CI from naive bootstrap

simulation.3 <- function(model.id, kernel, heteroskedasticity = F) {
  
  if (heteroskedasticity == T) {
    avebw <- average_bw(model.id, kernel)
    gen_data <- function() generate.data.h(model.id, avebw)
  } else {
    gen_data <- function() generate.data(model.id)
  }
  
  simu <- function() {
    dta <- gen_data()
    naiveboot(dta$y, dta$x, kernel, N.ci)
  }
  
  cl <- makeCluster(N.core)
  registerDoParallel(cl)
  export.obj <- c("N.ci", "a", "generate.data", "kweight", "generate.data.h",
                  "local_polynomial_regression", "naiveboot", "rd.estimate")
  collect.simu <- foreach(i=1:N.simu, .combine="rbind", .packages="rdrobust", .export=export.obj, .inorder=F) %dorng%
    simu()
  stopCluster(cl)
  
  t.true    <- ifelse(model.id ==2, -3.45, 0.04)
  bias      <- t.true - mean(collect.simu[ , 1])
  SD        <- sd(collect.simu[ , 1])
  MSE       <- sqrt(mean((collect.simu[ , 1] - t.true)^2))
  coverage  <- mean(collect.simu[ , 2] <= t.true & collect.simu[ , 3] >= t.true)
  length    <- mean(collect.simu[ , 3] - collect.simu[ , 2])
  
  return(c(true = t.true, bias = bias, SD = SD, MSE = MSE, coverage = coverage, length = length))
}

## simulation 4: coverage of traditional approach

simulation.4 <- function(model.id, kernel, heteroskedasticity = F) {
  
  if (heteroskedasticity == T) {
    avebw <- average_bw(model.id, kernel)
    gen_data <- function() generate.data.h(model.id, avebw)
  } else {
    gen_data <- function() generate.data(model.id)
  }
  
  simu <- function() {
    dta <- gen_data()
    bws <- rdbwselect_2014(dta$y, dta$x, kernel=kernel, bwselect = "IK")$bws[1]
    results <- rdrobust(dta$y, dta$x, kernel = kernel,
                        h = bws, b = bws, level = 100 * (1-a), vce = "hc3")
    t <- results$coef[1]
    ci <- results$ci[1, ]
    return(c(t, ci))
  }
  
  cl <- makeCluster(N.core)
  registerDoParallel(cl)
  export.obj <- c("generate.data", "generate.data.h", "a")
  collect.simu <- foreach(i=1:N.simu, .combine="rbind", .packages="rdrobust", .export=export.obj, .inorder=F) %dorng%
    simu()
  stopCluster(cl)
  
  t.true    <- ifelse(model.id ==2, -3.45, 0.04)
  bias      <- t.true - mean(collect.simu[ , 1])
  SD        <- sd(collect.simu[ , 1])
  MSE       <- sqrt(mean((collect.simu[ , 1] - t.true)^2))
  coverage  <- mean(collect.simu[ , 2] <= t.true & collect.simu[ , 3] >= t.true)
  length    <- mean(collect.simu[ , 3] - collect.simu[ , 2])
  
  return(c(true = t.true, bias = bias, SD = SD, MSE = MSE, coverage = coverage, length = length))
}

## generate tables

kernel <- "uni"
models <- c(1, 2, 3)
heteroskedasticity <- c(F, T)

rnames <- paste(c("m1", "h.m1", "m2", "h.m2","m3", "h.m3"), kernel, sep = ".")
boot.results <- matrix(0, nrow = 6, ncol = 6)
rownames(boot.results) <- rnames
cct.results <- matrix(0, nrow = 6, ncol = 6)
rownames(cct.results) <- rnames
naiveboot.results <- matrix(0, nrow = 6, ncol = 6)
rownames(naiveboot.results) <- rnames

i <- 1
for (m in models) {
  for (h in heteroskedasticity) {
    set.seed(798)
    boot.results[i, ] <- simulation.1(m, kernel, h)
    try(slackr(print(boot.results)), T)
    
    set.seed(798)
    cct.results[i, ] <- simulation.2(m, kernel, h)
    try(slackr(print(cct.results)), T)
    
    set.seed(798)
    naiveboot.results[i, ] <- simulation.3(m, kernel, h)
    try(slackr(print(naiveboot.results)), T)
    
    i <- i + 1
  }
}

print(xtable(boot.results, digits = 3), include.rownames = F, type = "latex", 
      file = paste("output/boot_results_", kernel, ".tex", sep = ""))
print(xtable(cct.results, digits = 3), include.rownames = F, type = "latex", 
      file = paste("output/cct_results_", kernel, ".tex", sep = ""))
print(xtable(naiveboot.results, digits = 3), include.rownames = F, type = "latex", 
      file = paste("output/naiveboot_results_", kernel, ".tex", sep = ""))
