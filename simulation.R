rm(list = ls())
library(rdrobust)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)

# N_bc:   number of bootstraps for bias correction
# N_ci:   number of bootstraps for CI construction
# N_simu: number of simulations to calculate coverage
# N_core: number of cores used for parallel computing

# generate_data:  a function to generate data for three models in Calonic 2014
# correct_bias:   a function to perform residual sampling and correct bias
# rd_estimate:    a function to calculate RD estimate
# rd_ci:          a function to calculate CI for RD estimate
# simulate:       a function to perform one simulation
# coverage:       a function to collect simulation results from parallel computing

N_bc    <- 50
N_ci    <- 1000
N_simu  <- 100
N_core  <- 15


generate_data <- function(model_n) {  
# input is a scalar: which DGP we are using
  
  x <- 2*rbeta(500, 2, 4) - 1
  e <- rnorm(500, 0, 0.1295)
  
  if (model_n == 1) {
    y <- (0.48 + 1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5)*(x < 0) +
         (0.52 + 0.84*x - 3.00*x^2 + 7.99*x^3 - 9.01*x^4 + 3.56*x^5)*( x >= 0)
  }
  
  if (model_n == 2) {
    y <- (3.71 + 2.30*x + 3.28*x^2 + 1.45*x^3 + 0.23*x^4 + 0.03*x^5)*(x < 0) +
         (0.26 + 18.49*x - 54.81*x^2 + 74.30*x^3 - 45.02*x^4 + 9.83*x^5)*(x >= 0)
  }
  
  if (model_n == 3) {
    y <- (0.48 + 1.27*x - 0.5*7.18*x^2 + 0.7*20.21*x^3 + 1.1*21.54*x^4 + 1.5*7.33*x^5)*(x < 0) +
      (0.52 + 0.84*x - 0.1*3.00*x^2 - 0.3*7.99*x^3 - 0.1*9.01*x^4 + 3.56*x^5)*( x >= 0)
  }
  
  return(data.frame(y = y, x = x))
}



correct_bias <- function(dta, q) {  
# dta is a data set which contains x and y. q is the order of the
# lower order polynomial. The higher order is q+1. q=1 or 1=2.
  
  if (q == 1) {
    l_ols <- lm(y ~ x + r + x:r, data = dta)
    h_ols <- lm(y ~ x + x_sq + r + x:r + x_sq:r, data = dta)
  }
  
  if (q == 2) {
    l_ols <- lm(y ~ x + x_sq + r + x:r + x_sq:r, data = dta)
    h_ols <- lm(y ~ x + x_sq + x_cu + r + x:r + x_sq:r + x_cu:r, data = dta)
  }
  
  dta$y_hat <- predict(h_ols)
  dta$e_hat <- residuals(h_ols)
  trd_hat <- coefficients(l_ols)[["rTRUE"]]
  trd <- coefficients(h_ols)[["rTRUE"]]
  
  trd_hat_star <- vector(length = N_bc)
  for (i in 1:N_bc) {
    dta_star <- data.frame(y = dta$y_hat + dta$e_hat[sample(c(1:nrow(dta)),nrow(dta), replace = T)],
                           x = dta$x, x_sq = dta$x_sq,
                           r = dta$r)
    if (q == 1) {
      trd_hat_star[i] <- coefficients(lm(y ~ x + r + x:r, data = dta_star))[["rTRUE"]]
    }
    if (q == 2) {
      trd_hat_star[i] <- coefficients(lm(y ~ x + x_sq + r + x:r + x_sq:r, data = dta_star))[["rTRUE"]]
    }
  }
  
  # T_h is the "True" effect from higher order polynomial
  # T_l is the estimated effect from lower order polynomial
  # T_bc is the estimated effect after bias correction
  return(c(T_h = trd, T_l = trd_hat, T_bc = trd_hat - (mean(trd_hat_star) - trd)))
}



rd_estimate <- function(dta, q = 1) {
# dta is a data set contains x and y
# q is the order of lower order polynomial
  
  bw <- rdbwselect(dta$y, dta$x, bwselect = "CCT")$bws[1]/2
  dta <- dta[abs(dta$x) <= bw, ]
  
  # if the window contains too few observations, return NA
  if ((q == 1 & nrow(dta) < 6) | (q == 2 & nrow(dta) < 8)) {
    return(c(N = NA, T_h = NA, T_l = NA, T_bc = NA))
  }
  dta$x_sq <- dta$x^2
  dta$x_cu <- dta$x^3
  dta$r <- dta$x >= 0
  
  # also return how many observations are used
  # apply the function "correct_bias"
  return(c(N = nrow(dta), correct_bias(dta, q)))
}



rd_ci <- function(dta, q = 1, level = c(0.025, 0.975)) {
# dta is a data set contains x and y
# q is the order of lower order polynomial
# level defines the percentiles used for CI
  
  bw <- rdbwselect(dta$y, dta$x, bwselect = "CCT")$bws[1]/2
  dta <- dta[abs(dta$x) <= bw, ]
  
  # if the window contains too few observations, return NA
  if ((q == 1 & nrow(dta) < 6) | (q == 2 & nrow(dta) < 8)) {
    return(c(NA, NA))
  }
  dta$x_sq <- dta$x^2
  dta$x_cu <- dta$x^3
  dta$r <- dta$x >= 0
  
  if (q == 1) {
    l_ols <- lm(y ~ x + r + x:r, data = dta)
    h_ols <- lm(y ~ x + x_sq + r + x:r + x_sq:r, data = dta)
  }
  
  if (q == 2) {
    l_ols <- lm(y ~ x + x_sq + r + x:r + x_sq:r, data = dta)
    h_ols <- lm(y ~ x + x_sq + x_cu + r + x:r + x_sq:r + x_cu:r, data = dta)
  }
  
  dta$y_hat <- predict(h_ols)
  dta$e_hat <- residuals(h_ols)
  
  trd_bc_star <- vector(length = N_ci)
  for (j in 1:N_ci) {
    dta_star <- data.frame(y = dta$y_hat + dta$e_hat[sample(c(1:nrow(dta)),nrow(dta), replace = T)],
                           x = dta$x, x_sq = dta$x_sq, x_cu = dta$x_cu, r = dta$r)
    
    # apply the function "correct_bias"
    trd_bc_star[j] <- correct_bias(dta_star, q)["T_bc"]
  }
  ci <- quantile(trd_bc_star, level, na.rm = T)
  return(ci)
}



simulate <- function(model_n, q = 1) {
# model_n is a scalar: which DGP we are using
# q is the order of lower order polynomial
  
  dta <- generate_data(model_n)
  estimates <- rd_estimate(dta, q)
  
  # if there is not enough observatioins in the window, return NA
  if (is.na(estimates[1])) return(rep(NA, 7))
  ci <- rd_ci(dta, q)
  
  # determine the true effect for each DGP
  if (model_n == 1 | model_n == 3) t_true <- 0.04
  if (model_n == 2) t_true <- -3.45
  covered <- as.numeric(ci[1] <= t_true & ci[2] >= t_true)
  return(c(estimates, ci, covered = covered))
}



coverage <- function(model_n, q = 1) {
# model_n is a scalar: which DGP we are using
# q is the order of lower order polynomial
  
  cl <- makeCluster(N_core)
  registerDoParallel(cl)
  export_obj <- c("N_bc", "N_ci", "generate_data", "correct_bias", "rd_estimate", "rd_ci", "simulate")
  collect_simu <- foreach(i=1:N_simu, .combine="rbind", .packages="rdrobust", .export=export_obj, .inorder=F) %dopar% 
    simulate(model_n, q)
  stopCluster(cl)
  
  # return a N_simu X 7 matrix. Each row is result from a single simulation
  return(collect_simu)
}

set.seed(200)
results_1.1 <- coverage(1, 1)
results_2.1 <- coverage(2, 1)
results_3.1 <- coverage(3, 1)
results_1.2 <- coverage(1, 2)
results_2.2 <- coverage(2, 2)
results_3.2 <- coverage(3, 2)

summary(results_1.1)
summary(results_2.1)
summary(results_3.1)
summary(results_1.2)
summary(results_2.2)
summary(results_3.2)









