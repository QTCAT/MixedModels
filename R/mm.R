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


#' @title Mixed Models incopareting relationship matrices
#' 
#' @description Fit (Generalized) Linear Mixed Models incopareting user 
#' defined 'G-site' covariance relationship matrices. Derifed from 
#' \code{\link[pedigreemm]{pedigreemm}}.
#' 
#' @inheritParams lme4::lmer 
#' @param family a GLM family, see \code{\link[stats]{glm}} and 
#' \code{\link[stats]{family}}.
#' @param covarrel a named list of relationship matrices. The names must 
#' correspond to the names of grouping factors for random-effects terms in the 
#' \code{formula} argument.
#' 
#' 
#' @details All arguments to this function are the same as those to the 
#' function \code{\link[lme4]{lmer}} (or in case of family 
#' \code{\link[lme4]{glmer}}) except \code{covarrel} which must be a 
#' named list of relationship matrices.  Each name (frequently 
#' there is only one) must correspond to the name of a grouping factor in a 
#' random-effects term in the \code{formula}.  The observed levels of that 
#' factor must be contained in the column and row names of the relationship 
#' matrix.  For each relationship matrix the (left) Cholesky factor of the 
#' observed levels is calculated and applied to the model matrix for that term.
#' 
#' @return a \code{\linkS4class{pedigreemm}} object.
#' 
#' @references Vazquez, A.I., D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. 
#' Weigel. (2010). Technical Note: An R package for fitting generalized linear 
#' mixed models in animal breeding. Journal of Animal Science, 88:497-504.
#' 
#' @keywords models
#' @export
relMM <-  function(formula, data, family = NULL, REML = TRUE, 
                   covarrel = list(), control = list(), start = NULL, 
                   verbose = FALSE, subset, weights, na.action, 
                   offset, contrasts = NULL, devFunOnly = FALSE, ...) {
  gaus <- FALSE
  if (is.null(family)) {
    gaus <- TRUE
  } else {
    ## copied from glm()
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (!inherits(family, "family")) stop("unknown family type")
    gaus <- family$family == "gaussian" && family$link == "identity"
  }
  mc <- match.call()
  # create a call to lmer
  lmerc <- mc
  lmerc[[1]] <- if (gaus) as.name("lmer") else as.name("glmer")
  lmerc$covarrel <- NULL
  if (!gaus) lmerc$REML <- NULL
  # call [g]lmer instead
  if (!length(covarrel))
    return(eval.parent(lmerc))
  # check the covarrel argument
  stopifnot(is.list(covarrel),
            length(names(covarrel)) == length(covarrel))
  if (!all(sapply(covarrel, function(x) {
    nrow(x) == ncol(x) && 
      all(rownames(x) == colnames(x)) && 
      (inherits(x, "Matrix") || inherits(x, "matrix"))
  })))
    stop("every relatinships matrix must be a square matrix which identical column and row names")
  lmf <- eval(lmerc, parent.frame())
  # copy the pedigree list for relfactor
  relfac <- covarrel
  pnms <- names(covarrel)
  stopifnot(all(pnms %in% names(lmf@flist)))
  asgn <- attr(lmf@flist, "assign")
  Zt <- lmf@pp$Zt
  for (i in seq_along(covarrel)) {
    tn <- which(match(pnms[i], names(lmf@flist)) == asgn)
    if (length(tn) > 1)
      stop("each relationship matrix must be associated with only one random-effect term")
    ind <- (lmf@Gp)[tn:(tn + 1L)]
    Zti <- (ind[1] + 1L):ind[2]
    if (!all(rownames(Zt)[Zti] %in% rownames(covarrel[[i]])))
      stop("levels of a random-effect term are not part of the according relationship matrix")
    relfaci <- match(rownames(covarrel[[i]]), rownames(Zt)[Zti])
    relfac[[i]] <- tryCatch(chol(covarrel[[i]][relfaci, relfaci]),
                            error = function(x) {
                              message("each relatinships matrices must be positive semidefinite")
                              message(x)
                              return(NULL)
                            })
    Zt[Zti, ] <- relfac[[i]] %*% Zt[Zti, ]
  }
  reTrms <- list(Zt = Zt, theta = lmf@theta, Lambdat = lmf@pp$Lambdat,
                 Lind = lmf@pp$Lind, lower = lmf@lower, flist = lmf@flist,
                 cnms = lmf@cnms, Gp = lmf@Gp)
  dfl <- list(fr = lmf@frame, X = lmf@pp$X, reTrms = reTrms, start = lmf@theta)
  if (gaus) {
    dfl$REML <- lmf@resp$REML > 0L
    devfun <- do.call(mkLmerDevfun, dfl)
    opt <- optimizeLmer(devfun, optimizer = "Nelder_Mead", ...)
  } else {
    dfl$family <- family
    devfun <- do.call(mkGlmerDevfun, dfl)
    opt <- optimizeGlmer(devfun, optimizer = "Nelder_Mead", ...)
  }
  mm <- mkMerMod(environment(devfun), opt, reTrms, lmf@frame, mc)
  cls <- if (gaus) "lmerpedigreemm" else "glmerpedigreemm"
  ans <- do.call(new, list(Class = cls, relfac = relfac, frame = mm@frame, 
                           flist = mm@flist, cnms = mm@cnms, Gp = mm@Gp,
                           theta = mm@theta, beta = mm@beta, u = mm@u, 
                           lower = mm@lower, devcomp = mm@devcomp, pp = mm@pp, 
                           resp = mm@resp, optinfo = mm@optinfo))
  ans@call <- evalq(mc)
  ans
}


# setMethod("ranef", signature(object = "pedigreemm"),
#           function(object, postVar = FALSE, drop = FALSE, 
#                    whichel = names(ans), pedigree = TRUE, ...) {
#             if ((postVar <- as.logical(postVar)) && (pedigree <- as.logical(pedigree)))
#               stop("code for applying pedigree and posterior variances not yet written")
#             ans <- ranef(as(object, "merMod"), postVar, drop = FALSE)
#             ans <- ans[whichel]
#             if (pedigree) {
#               if (postVar)
#                 stop("postVar and pedigree cannot both be true")
#               rf <- object@relfac
#               for (nm in names(rf)) {
#                 dm <- data.matrix(ans[[nm]])
#                 cn <- colnames(dm)
#                 rn <- rownames(dm)
#                 dm <- as.matrix(rf[[nm]] %*% dm)
#                 colnames(dm) <- cn
#                 rownames(dm) <- rn
#                 ans[[nm]] <- data.frame(dm, check.names = FALSE)
#               }
#             }
#             if (drop)
#               ans <- lapply(ans, function(el)
#               {
#                 if (ncol(el) > 1) return(el)
#                 pv <- drop(attr(el, "postVar"))
#                 el <- drop(as.matrix(el))
#                 if (!is.null(pv))
#                   attr(el, "postVar") <- pv
#                 el
#               })
#             ans
#           })
