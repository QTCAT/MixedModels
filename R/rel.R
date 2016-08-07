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
#' @examples
#' # simulate data
#' set.seed(123)
#' x <- matrix(rbinom(1000, 2, .7)/2, 10, 100)
#' rownames(x) <- paste0("indiv", 1:nrow(x))
#' # estimated relationship matrix
#' genorel <- grm(x)
#'
#' @importFrom Matrix nearPD
#' @export
grm <- function(x, checkPD = TRUE, ...) {
  stopifnot(is.matrix(x))
  if (any(is.na(x)))
    stop("in 'x' missing data are not allowed")
  if (is.null(rownames(x)))
    stop("row names of 'x' are missing")
  if (!all(na.omit(unique(as.vector(x))) %in% c(0, .5, 1)))
    stop("alleles have to be coded by AA = 0, AB = 0.5, and BB = 1")
  n <- nrow(x)
  p <- ncol(x)
  f <- colMeans(x)
  w <- x - matrix(rep(f, each = n), n, p)
  A <- tcrossprod(w) / (p / 2 * mean(f * (1 - f)))
  if (checkPD)
    A <- tryCatch(nearPD(A, ...)$mat,
                  error = function(x) {
                    message("it was not possible to find a positive definite matrix")
                    message(x)
                    return(NULL)
                  })
  A
}
