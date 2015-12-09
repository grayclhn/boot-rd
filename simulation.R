rm(list = ls())
library(rdrobust)
library(foreach)
library(doParallel)

# N_bc:   number of bootstraps for bias correction
# N_ci:   number of bootstraps for CI construction
# N_simu: number of simulations to calculate coverage
# N_core: number of cores used for parallel computing

# generate_data:  a function to generate data for three models in Calonic 2014
# rd_estimate:    a function to calculate RD estimate
# rd_ci:          a function to calculate CI for RD estimate
# simulate:       a function to perform one simulation
# coverage:       a function to collect results and calculate coverage

N_bc    <- 20
N_ci    <- 1000
N_simu  <- 100
N_core  <- 15

generate_data <- function(model_n) {
  
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

rd_estimate <- function(dta) {
  bw <- rdbwselect(dta$y, dta$x, bwselect = "CCT")$bws[1]
  dta <- dta[abs(dta$x) <= bw, ]
  N_obs <- nrow(dta)
  dta$x_sq <- dta$x^2
  dta$r <- dta$x >= 0
  
  linear_ols <- lm(y ~ x + r + x:r, data = dta)
  quadratic_ols <- lm(y ~ x + x_sq + r + x:r + x_sq:r, data = dta)
  
  dta$y_hat <- predict(quadratic_ols)
  dta$e_hat <- residuals(quadratic_ols)
  trd_hat <- coefficients(linear_ols)[["rTRUE"]]
  trd <- coefficients(quadratic_ols)[["rTRUE"]]
  
  ## correct bias
  trd_hat_star <- vector(length = N_bc)
  for (i in 1:N_bc) {
    dta_star <- data.frame(y = dta$y_hat + dta$e_hat[sample(c(1:N_obs),N_obs, replace = T)],
                           x = dta$x,
                           r = dta$r)
    trd_hat_star[i] <- coefficients(lm(y ~ x + r + x:r, data = dta_star))[["rTRUE"]]
  }
  
  trd_bc <- trd_hat - (mean(trd_hat_star) - trd)
  
  return(trd_bc)
}

rd_ci <- function(dta, level = c(0.025, 0.975)) {
  trd_bc_star <- vector(length = N_ci)
  for (j in 1:N_ci) {
    dta_star <- dta[sample(1:nrow(dta), nrow(dta), replace = T), ]
    trd_bc_star[j] <- rd_estimate(dta_star)
  }
  ci <- quantile(trd_bc_star, level, na.rm = T)
  return(ci)
}

simulate <- function(seed, model_n) {
  set.seed(seed)
  ci <- rd_ci(generate_data(model_n))
  if (model_n == 1 | model_n == 3) t_true <- 0.04
  if (model_n == 2) t_true <- -3.45
  covered <- as.numeric(ci[1] <= t_true & ci[2] >= t_true)
  return(covered)
}

coverage <- function(model_n) {
  cl <- makeCluster(N_core)
  registerDoParallel(cl)
  export_obj <- c("N_bc", "N_ci", "generate_data", "rd_estimate", "rd_ci", "simulate")
  collect_simu <- foreach(i=1:N_simu, .combine="rbind", .packages="rdrobust", .export=export_obj, .inorder=F) %dopar% 
    simulate(i, model_n)
  stopCluster(cl)
  mean(collect_simu)
}

coverage(1)
coverage(2)
coverage(3)












