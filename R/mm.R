#' Genetic relationship matrix
#'
#' Computers a genetic relationship matrix (GRM) from SNP data.
#'
#' @param x a matrix (n x p) of n individuals (rows) and p SNPs (columns). The
#' alleles have to be coded by AA = 0, AB = 0.5, and BB = 1. Missing data are
#' not yet allowed.
#'
#' @importFrom Matrix nearPD
#' @export
grm <- function (x) {
  stopifnot(is.matrix(x))
  stopifnot(!is.na(x))
  if (!all(na.omit(unique(as.vector(x))) %in% c(0, .5, 1)))
    stop("alleles have to be coded by AA = 0, AB = 0.5, and BB = 1")
  f <- colMeans(x)
  v <- mean(f * (1 - f))
  w <- x - matrix(rep(f, each = nrow(x)), nrow(x))
  A <- tcrossprod(w) / (ncol(x) / 2 * v)
  tryCatch(expr = nearPD(A)$mat,
           error = message("it was not possible to find a positive definite matrix"))
}
