## Copyright (c) 2016, Otávio Bartalotti, Gray Calhoun, and Yang He.
## Available under the MIT "Expat" License, see README.md

library(rdrobust)
library(RDbootstraps)
library(foreach)
library(doParallel)
library(doRNG)
library(xtable)

try(source("slackrISE.R"), T)

#### global parameters ####
Nci <- 999
Nbc <- 500
Nsimu <- 5000
Ncore <- 32

#### functions for simulation ####

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
  return(data.frame(y = y, x = x))
}

naiveboot <- function(y, x, kernel, Nci = 999) {
  ## point estimate
  tau <- rdrobust(y, x, kernel = kernel, bwselect = "IK")$coef[1]
  ## bootstrap
  tauboot <- function(y, x) {
    i <- sample.int(length(y), length(y), replace = T)
    rdrobust(y[i], x[i], kernel = kernel, bwselect = "IK")$coef[1]
  }
  Btau <- replicate(Nci, tauboot(y, x))

  return(c(tau,
           quantile(Btau, c(0.05, 0.95)),
           quantile(Btau, c(0.025, 0.975)),
           quantile(Btau, c(0.005, 0.995))))
}

rdsimu <- function(model.id, kernel) {

  ## a function for each random draw
  eachdraw <- function() {
    dta   <- generate.data(model.id)

    ## naive bootstraps + IK bandwidth
    boot.naive <- naiveboot(dta$y, dta$x, kernel = kernel, Nci = Nci)

    ## robust bootstraps + CCT bandwidth
    boot.robust <- rdboot(dta$y, dta$x, kernel = kernel, Nci = Nci, Nbc = Nbc)

    ## traditional + IK bandwidth
    result <- rdrobust(dta$y, dta$x, kernel = kernel, bwselect = "IK")
    tau <- result$coef[1]
    ci90 <- c(tau - qnorm(0.95)*result$se[1], tau + qnorm(0.95)*result$se[1])
    ci95 <- c(tau - qnorm(0.975)*result$se[1], tau + qnorm(0.975)*result$se[1])
    ci99 <- c(tau - qnorm(0.995)*result$se[1], tau + qnorm(0.995)*result$se[1])
    analytic.ik <- c(tau, ci90, ci95, ci99)

    ## robust + CCT bandwidth
    result <- rdrobust(dta$y, dta$x, kernel = kernel)
    tau <- result$coef[3]
    ci90 <- c(tau - qnorm(0.95)*result$se[3], tau + qnorm(0.95)*result$se[3])
    ci95 <- c(tau - qnorm(0.975)*result$se[3], tau + qnorm(0.975)*result$se[3])
    ci99 <- c(tau - qnorm(0.995)*result$se[3], tau + qnorm(0.995)*result$se[3])
    analytic.cct <- c(tau, ci90, ci95, ci99)

    return(c(boot.naive, boot.robust, analytic.ik, analytic.cct))
  }

  ## parallel computing
  cl <- makeCluster(Ncore)
  registerDoParallel(cl)
  export.obj <- c("Nci", "Nbc", "generate.data", "naiveboot")
  result <- foreach(i=1:Nsimu, .combine="rbind", .packages=c("rdrobust", "RDbootstraps"),
                    .export=export.obj, .inorder=F) %dorng% eachdraw()
  stopCluster(cl)

  ## clean the results
  tau    <- ifelse(model.id ==2, -3.45, 0.04)
  Nestimator <- 4
  table <- matrix(0, nrow = Nestimator, ncol = 9)
  colnames(table) <- c("bias", "SD", "RMSE",
                       "coverage90", "length90",
                       "coverage95", "length95",
                       "coverage99", "length99")
  rownames(table) <- c("boot.naive", "boot.robust",
                       "analytic.ik", "analytic.cct")
  for (i in 1:Nestimator) {
    table[i, 1] <- tau - mean(result[ , (7*(i - 1) + 1)])
    table[i, 2] <- sd(result[ , (7*(i - 1) + 1)])
    table[i, 3] <- sqrt(mean((result[ , (7*(i - 1) + 1)] - tau)^2))
    table[i, 4] <- mean(result[ , (7*(i - 1) + 2)] <= tau &
                        result[ , (7*(i - 1) + 3)] >= tau)
    table[i, 5] <- mean(result[ , (7*(i - 1) + 3)] - result[ , (7*(i - 1) + 2)])
    table[i, 6] <- mean(result[ , (7*(i - 1) + 4)] <= tau &
                          result[ , (7*(i - 1) + 5)] >= tau)
    table[i, 7] <- mean(result[ , (7*(i - 1) + 5)] - result[ , (7*(i - 1) + 4)])
    table[i, 8] <- mean(result[ , (7*(i - 1) + 6)] <= tau &
                          result[ , (7*(i - 1) + 7)] >= tau)
    table[i, 9] <- mean(result[ , (7*(i - 1) + 7)] - result[ , (7*(i - 1) + 6)])
  }
  return(round(table, digits = 3))
}

#### organize results ####

kernels <- c("uni", "tri", "epa")
models <- c(1, 2, 3)

for (k in kernels) {
  for (m in models) {
    set.seed(798)
    km <- paste(k, m, sep = ".")
    assign(km, rdsimu(m, k))
    try(slackr(print(km)), T)
    try(slackr(print(eval(parse(text = km)))), T)
  }
}

save.image("output/SRDsimuresult.RData")

