#' Genetic relationship matrix
#'
#' Computers a genetic relationship matrix (GRM) from SNP data.
#'
#' @param x a matrix (n x p) of n individuals (rows) and p SNPs (columns). The
#' alleles have to be coded by AA = 0, AB = 0.5, and BB = 1. Missing data are
#' not yet allowed.
#' @param checkPD logickal, if true (default) the nearest positive definite 
#' matrix to the GRM is compute if the GRM does not fulfill these 
#' requirements anyways.
#' @param ... additional arguments to \code{\link[Matrix]{nearPD}}.
#'
#' @importFrom Matrix nearPD
#' @export
grm <- function(x, checkPD = TRUE, ...) {
  stopifnot(is.matrix(x))
  stopifnot(!is.na(x))
  if (!all(na.omit(unique(as.vector(x))) %in% c(0, .5, 1)))
    stop("alleles have to be coded by AA = 0, AB = 0.5, and BB = 1")
  f <- colMeans(x)
  v <- mean(f * (1 - f))
  w <- x - matrix(rep(f, each = nrow(x)), nrow(x))
  A <- tcrossprod(w) / (ncol(x) / 2 * v)
  if (checkPD)
    A <- tryCatch(nearPD(A, ...)$mat,
                  error = function(x) {
                    message("it was not possible to find a positive definite matrix")
                    message(x)
                    return(NULL)
                  })
  A
}


#' Rescale a relationship matrix
#'
#' Rescale a relationship matrix for a realistic variance component estimation.
#'
#' @param x a symmetric matrix.
#' @param tolerance numeric value giving the tolerance for the rescaling.
#'
#' @export
rescale <- function(x, tolerance = 1e-3) {
  n <- ncol(x)
  scpar <- ((sum(diag(x)) - sum(x) / n) / n)
  if (abs(1 - scpar) > tolerance)
    x <- x / scpar
  x
}
