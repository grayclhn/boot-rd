library(rdrobust)
library(foreach)
library(doParallel)
library(doRNG)
library(xtable)
try(source("slackrISE.R"), T)

source("RDfunctions.R")

#### for testing in my computer ####
# N.ci <- 20
# level <- 0.95
# N.simu <- 20
# N.core <- 2
# N.bc <- 20

#### to run in servor ####
N.ci <- 999
level <- 0.95
N.simu <- 1500
N.core <- 15
N.bc <- 500

## simulation 1: coverage of CI from bootstrap

simulation.1 <- function(model.id, p, q, kernel, bc) {
  
  simu <- function() {
    dta   <- generate.data(model.id)
    # CCT bandwidth from Calonico package
    bw    <- rdbwselect(dta$y, dta$x, p = p, q = q, kernel = kernel, bwselect = "CCT")
    bw.p  <- bw$bws[1,1]
    bw.q  <- bw$bws[1,2]
    t     <- rd.estimate(dta, p, q, bw.p, bw.q, N.bc, kernel, bc)
    ci    <- rd.ci(dta, p, q, bw.p, bw.q, N.bc, N.ci, level, kernel, bc)$ci
    return(c(t, ci))
  }
  
  cl <- makeCluster(N.core)
  registerDoParallel(cl)
  export.obj <- c("N.ci", "N.bc", "level", "generate.data", "kweight", "rd.estimate", "lpreg", "rd.ci")
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

simulation.2 <- function(model.id, p, q, kernel) {
  
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

gen.table <- function(model.id, p, q){
  table <- matrix(0, nrow = 4, ncol = 6)
  colnames(table) <- c("true", "bias", "SD", "MSE", "coverage", "length")
  rownames(table) <- c("cct.tri", "cct.uni", "boot.uni", "boot.uni(uncorrected)")

  try(environment(slackr) <- environment(), T)
  
  table[1, ] <- simulation.2(model.id, p, q, "tri")
  try(slackr(print(table[1,])), T)
  table[2, ] <- simulation.2(model.id, p, q, "uni")
  try(slackr(print(table[2,])), T)
  table[3, ] <- simulation.1(model.id, p, q, "uni", T)
  try(slackr(print(table[3,])), T)
  table[4, ] <- simulation.1(model.id, p, q, "uni", F)
  try(slackr(print(table[4,])), T)
  
  table
}

set.seed(798)
table.112 <- gen.table(1, 1, 2)
try(slackr("table.112 done 1/6"), T)

table.212 <- gen.table(2, 1, 2)
try(slackr("table.212 done 2/6"), T)

table.312 <- gen.table(3, 1, 2)
try(slackr("table.312 done 3/6"), T)

table.123 <- gen.table(1, 2, 3)
try(slackr("table.123 done 4/6"), T)

table.223 <- gen.table(2, 2, 3)
try(slackr("table.223 done 5/6"), T)

table.323 <- gen.table(3, 2, 3)
try(slackr("table.323 done 6/6"), T)

## table 1: p = 1, q = 2

table.1 <- rbind(table.112, table.212, table.312)
print(xtable(table.1, digits = 3), include.rownames = F, type = "latex", file = "output/simu_table_1.tex")

## table 2: p = 2, q = 3

table.2 <- rbind(table.123, table.223, table.323)
print(xtable(table.2, digits = 3), include.rownames = F, type = "latex", file = "output/simu_table_2.tex")




