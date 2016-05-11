#' The main function to perform bootstraps in sharp regression discontinuity
#' design.
#'
#' @description This package was created to speed up the simulations using
#' multiple cores. Functions compilded from C++ code cannot be loaded to each
#' worker unless they are bundled into a package. This package is primarily only
#' for sharp RD simulations in this paper and is not designed for general use.
#' @param y is the dependent variable.
#' @param x is the running variable (score or forcing variable).
#' @param c specifies the RD cutoff in x; default is c = 0.
#' @param p specifies the order of the local-polynomial used to construct the
#'   point-estimator; default is p = 1 (local linear regression).
#' @param q specifies the order of the local-polynomial used to construct the
#'   bias-correction; default is q = 2 (local quadratic regression).
#' @param h directly sets the main bandwidth. If not specified, it is computed
#'   by the command rdrobust::rdbwselect.
#' @param b directly sets the pilot bandwidth. If not specified, it is computed
#'   by the command rdrobust::rdbwselect.
#' @param Nbc is the number of bootstraps for bias correction.
#' @param Nci is the number of bootstraps to construct percentile confidence
#'   interval.
#' @param kernel is the kernel function used to construct the local-polynomial
#'   estimator(s). Options are triangular (default option), epanechnikov and
#'   uniform.
#' @param level sets the confidence levels, a vector of length 3.
#' @return A vector containing point and interval estimates.
#' @export

rdboot <- function(y, x, c = 0,  p = 1, q = 2,
                   h = NULL, b = NULL, Nbc = 500, Nci = 999,
                   kernel = "uni", level = c(90, 95, 99)) {

  ## prepare the data
  data <- data.frame(y = y, x = x - c)

  if (is.null(h) | is.null(b)) {
    bws <- rdrobust::rdbwselect(data$y, data$x, p = p, q = q, kernel = kernel)$bws
    h  <- bws[1,1]
    b  <- bws[1,2]
  }

  ## keep data within distance b around 0
  data <- data[abs(data$x) < b, ]
  data$wh <- rdrobust::kweight(data$x, 0, h, kernel)
  data$wb <- rdrobust::kweight(data$x, 0, b, kernel)

  left <- as.matrix(data[data$x < 0, ])
  right <- as.matrix(data[data$x >= 0, ])

  ## estimation
  tau <- srdbc(left, right, Nbc, p, q, 1)[3]      # point estimate
  Btau <- srdbcboot(left, right, Nbc, Nci, p, q)  # bootstrap point estimate

  ci1 <- quantile(Btau, c((1 - level[1]/100)/2, 1 - (1 - level[1]/100)/2))
  ci2 <- quantile(Btau, c((1 - level[2]/100)/2, 1 - (1 - level[2]/100)/2))
  ci3 <- quantile(Btau, c((1 - level[3]/100)/2, 1 - (1 - level[3]/100)/2))

  return(c(tau, ci1, ci2, ci3))
}
