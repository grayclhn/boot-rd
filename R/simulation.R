rm(list = ls())
library(rdrobust)
library(foreach)
library(doParallel)
library(xtable)
try(source("slackrISE.R"), silent = T)

source("functions.R")

##########################################################################

set.seed(555)
N.ci <- 999
level <- 0.95
N.simu <- 500
N.core <- 15
N.bc <- 250
N.bw <- 500


## simulation 1: coverage of CI from bootstraping corrected estimator

simulation.1 <- function(model.id, p, q, kernel) {
  
  simu <- function() {
    dta   <- generate.data(model.id)
    # CCT bandwidth from Calonico package
    bw    <- rdbwselect(dta$y, dta$x, p = p, q = q, kernel = kernel, bwselect = "CCT")
    bw.p  <- bw$bws[1,1]
    bw.q  <- bw$bws[1,2]
    t     <- rd.estimate(dta, p, q, bw.p, bw.q, N.bc, kernel, T)
    ci    <- rd.ci(dta, p, q, bw.p, bw.q, N.bc, N.ci, level, kernel, T)
    return(c(t, ci))
  }
  
  cl <- makeCluster(N.core)
  registerDoParallel(cl)
  export.obj <- c("N.ci", "N.bc", "level", "generate.data", "k.weight", "rd.estimate", "lpreg", "rd.ci")
  collect.simu <- foreach(i=1:N.simu, .combine="rbind", .packages="rdrobust", .export=export.obj, .inorder=F) %dopar% 
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

output.1 <- as.data.frame(matrix(nrow = 18, ncol = 10))
colnames(output.1) <- c("DGP", "q", "p", "Kernel", "True", "Bias", "SD", "MSE", "CI coverage", "CI length")

output.1[ , 1] <- rep(c(1,2,3), each = 6)
output.1[ , 2] <- rep(c(1,2), 9)
output.1[ , 3] <- rep(c(2,3), 9)
output.1[ , 4] <- rep(rep(c("uni", "epa", "tri"), each = 2), 3)

output.1[1, 5:10] <- simulation.1(1,1,2,"uni")
output.1[2, 5:10] <- simulation.1(1,2,3,"uni")
output.1[7, 5:10] <- simulation.1(2,1,2,"uni")
output.1[8, 5:10] <- simulation.1(2,2,3,"uni")
output.1[13, 5:10] <- simulation.1(3,1,2,"uni")
output.1[14, 5:10] <- simulation.1(3,2,3,"uni")

output.1[3, 5:10] <- simulation.1(1,1,2,"epa")
output.1[4, 5:10] <- simulation.1(1,2,3,"epa")
output.1[9, 5:10] <- simulation.1(2,1,2,"epa")
output.1[10, 5:10] <- simulation.1(2,2,3,"epa")
output.1[15, 5:10] <- simulation.1(3,1,2,"epa")
output.1[16, 5:10] <- simulation.1(3,2,3,"epa")

output.1[5, 5:10] <- simulation.1(1,1,2,"tri")
output.1[6, 5:10] <- simulation.1(1,2,3,"tri")
output.1[11, 5:10] <- simulation.1(2,1,2,"tri")
output.1[12, 5:10] <- simulation.1(2,2,3,"tri")
output.1[17, 5:10] <- simulation.1(3,1,2,"tri")
output.1[18, 5:10] <- simulation.1(3,2,3,"tri")

print(xtable(output.1, digits = c(0,0,0,0,0,2,3,3,3,3,3)), include.rownames = F)


## simulation 4: coverage of CI from CCT(2004)

simulation.4 <- function(model.id, p, q, kernel) {
  
  simu <- function() {
    dta   <- generate.data(model.id)
    results <- rdrobust(dta$y, dta$x, p = p, q = q, kernel = kernel, bwselect = "CCT", level = 100*level)
    t     <- results$coef[3]
    ci    <- results$ci[3, ]
    return(c(t, ci))
  }
  
  cl <- makeCluster(N.core)
  registerDoParallel(cl)
  export.obj <- c("generate.data", "level")
  collect.simu <- foreach(i=1:N.simu, .combine="rbind", .packages="rdrobust", .export=export.obj, .inorder=F) %dopar% 
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

output.4 <- as.data.frame(matrix(nrow = 18, ncol = 10))
colnames(output.4) <- c("DGP", "q", "p", "Kernel", "True", "Bias", "SD", "MSE", "CI coverage", "CI length")

output.4[ , 1] <- rep(c(1,2,3), each = 6)
output.4[ , 2] <- rep(c(1,2), 9)
output.4[ , 3] <- rep(c(2,3), 9)
output.4[ , 4] <- rep(rep(c("uni", "epa", "tri"), each = 2), 3)

set.seed(555)
output.4[1, 5:10] <- simulation.4(1,1,2,"uni")
output.4[2, 5:10] <- simulation.4(1,2,3,"uni")
output.4[7, 5:10] <- simulation.4(2,1,2,"uni")
output.4[8, 5:10] <- simulation.4(2,2,3,"uni")
output.4[13, 5:10] <- simulation.4(3,1,2,"uni")
output.4[14, 5:10] <- simulation.4(3,2,3,"uni")

output.4[3, 5:10] <- simulation.4(1,1,2,"epa")
output.4[4, 5:10] <- simulation.4(1,2,3,"epa")
output.4[9, 5:10] <- simulation.4(2,1,2,"epa")
output.4[10, 5:10] <- simulation.4(2,2,3,"epa")
output.4[15, 5:10] <- simulation.4(3,1,2,"epa")
output.4[16, 5:10] <- simulation.4(3,2,3,"epa")

output.4[5, 5:10] <- simulation.4(1,1,2,"tri")
output.4[6, 5:10] <- simulation.4(1,2,3,"tri")
output.4[11, 5:10] <- simulation.4(2,1,2,"tri")
output.4[12, 5:10] <- simulation.4(2,2,3,"tri")
output.4[17, 5:10] <- simulation.4(3,1,2,"tri")
output.4[18, 5:10] <- simulation.4(3,2,3,"tri")

print(xtable(output.4, digits = c(0,0,0,0,0,2,3,3,3,3,3)), include.rownames = F)


## simulation 2: coverage of CI from bootstraping corrected estimator

simulation.2 <- function(model.id, p, q, kernel) {
  
  simu <- function() {
    dta   <- generate.data(model.id)
    # CCT bandwidth from Calonico package
    bw    <- rdbwselect(dta$y, dta$x, p = p, q = q, kernel = kernel, bwselect = "CCT")
    bw.p  <- bw$bws[1,1]
    bw.q  <- bw$bws[1,2]
    t     <- rd.estimate(dta, p, q, bw.p, bw.q, N.bc, kernel)
    ci    <- rd.ci(dta, p, q, bw.p, bw.q, N.bc, N.ci, level, kernel)
    return(c(t, ci))
  }
  
  cl <- makeCluster(N.core)
  registerDoParallel(cl)
  export.obj <- c("N.ci", "N.bc", "level", "generate.data", "k.weight", "rd.estimate", "lpreg", "rd.ci")
  collect.simu <- foreach(i=1:N.simu, .combine="rbind", .packages="rdrobust", .export=export.obj, .inorder=F) %dopar% 
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

# > simulation.2(1,1,2,"uni")
# true         bias          SD         MSE         coverage      length 
# 0.04000000  -0.01141711  0.06604530  0.06695975  0.94000000   0.23919618 
# > simulation.2(2,1,2,"uni")
# true         bias          SD         MSE         coverage      length 
# -3.45000000 -0.01241777  0.08648155  0.08728288  0.94600000   0.32030309 
# > simulation.2(3,1,2,"uni")
# true         bias          SD         MSE         coverage      length 
# 0.040000000 -0.005055316 0.063537170 0.063674597  0.960000000 0.245990453

## simulation 3: check bandwidth from bootstrap method

simulation.3 <- function(model.id, p, q, kernel) {
  
  simu <- function() {
    dta   <- generate.data(model.id)
    # CCT bandwidth from Calonico package
    bw    <- rdbwselect(dta$y, dta$x, p = p, q = q, kernel = kernel, bwselect = "CCT")
    bw1.p  <- bw$bws[1,1]
    bw1.q  <- bw$bws[1,2]
    # bandwidth from bootstrap method
    bw2.q <- bw.select(dta, q, q + 1, N.bw, kernel)$update
    bw2.p <- bw.select(dta, p, q, N.bw, kernel, bw2.q)$update
    return(c(bw1.p, bw1.q, bw2.p, bw2.q))
  }
  
  cl <- makeCluster(N.core)
  registerDoParallel(cl)
  export.obj <- c("N.bw", "generate.data", "k.weight", "lpreg", "bw.select")
  collect.simu <- foreach(i=1:N.simu, .combine="rbind", .packages="rdrobust", .export=export.obj, .inorder=F) %dopar% 
    simu()
  stopCluster(cl)
  
  compare <- rbind(colMeans(collect.simu), apply(collect.simu, 2, sd))
  rownames(compare) <- c("Mean", "SD")
  colnames(compare) <- c("CCT(p)", "CCT(q)", "boot(pL)", "boot(pR)", "boot(qL)", "boot(qR)")
  
  return(compare)
}

simulation.3(1,1,2,"uni")
simulation.3(2,1,2,"uni")
simulation.3(3,1,2,"uni")

# > simulation.3(1,1,2,"uni")
#         CCT(p)     CCT(q)   boot(pL)   boot(pR)   boot(qL)   boot(qR)
# Mean 0.15940068 0.29566152 0.12012389 0.11395367 0.15601003 0.15212266
# SD   0.03398398 0.05644396 0.03118821 0.02903325 0.02212816 0.02049559
# > simulation.3(2,1,2,"uni")
#         CCT(p)     CCT(q)  boot(pL)   boot(pR)   boot(qL)   boot(qR)
# Mean 0.075783035 0.19557517 0.1197757 0.06765035 0.15678818 0.15141334
# SD   0.008665331 0.01996184 0.0284621 0.01307661 0.02056165 0.02139609
# > simulation.3(3,1,2,"uni")
#         CCT(p)     CCT(q)   boot(pL)   boot(pR)   boot(qL)   boot(qR)
# Mean 0.14329344 0.29167328 0.11560102 0.11269013 0.15724649 0.15104778
# SD   0.02398352 0.05380659 0.02925731 0.02866268 0.02178484 0.02081259








