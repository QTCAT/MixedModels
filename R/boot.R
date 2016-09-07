#' Parametric bootstrapping
#'
#' Parametric bootstrapping for mixed model with relationship matrix
#'
#' @param model a \code{\link[pedigreemm]{pedigreemm}} model object.
#' @param func a function which extracts the parameter of interest from the 
#' model. This function should have one argument - the model object.
#' @param R the number of bootstrapping runs.
#' @param cores the number of cores for parallel computing at Unix-like 
#' machines.
#'
#' @importFrom boot boot
#' @importFrom stats update simulate
#' @export
relMMboot <- function(model, func, R = 1000, cores = 1) {
  out <- boot(list(data = model@frame, model = model, func = func),
              statistic = function(x) {
                fit <- update(x$model, data = x$data)
                x$func(fit)
              }, 
              ran.gen = function(x, ...) {
                x$data[, 1] <- unlist(simulate(x$model, use.u = TRUE))
                x
              },
              sim = "parametric", R = R,
              parallel = "multicore", ncpus = cores)
  out
}


#' Beta Approximation Confidence Intervals
#'
#' Using the normal approximation to a statistic, calculate beta distrinution
#' two-sided confidence intervals.
#'
#' @param object a bootstrap output object returned from a call to boot.
#' @param conf	a scalar or vector containing the confidence level of the 
#' required interval.
#' @param index the index of the statistic of interest within the output of a 
#' call to boot.out$statistic.
#'
#' @importFrom stats qbeta var
#' @export
beta.ci <- function(object, conf = 0.95, index = 1) {
  stopifnot(class(object) == "boot")
  p <- (1 + conf) / 2
  x <- object$t[, index[1]]
  mu0 <- object$t0[index[1]]
  # var (including mean bias) 
  var0 <- var(x) + (mean(x) - mu0)^2 * length(x) / (length(x) - 1)
  # from mean and var to alpha and beta
  a <- ((1 - mu0) / var0 - 1 / mu0) * mu0^2
  b <- a * (1 / mu0 - 1)
  c(conf, qbeta(1 - p, a, b), qbeta(p, a, b))
}

