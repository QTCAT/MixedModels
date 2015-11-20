
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
  if (!all(sapply(covarrel, function(x, tol = 1e-3) {
    nrow(x) == ncol(x) &&
      identical(rownames(x), colnames(x)) &&
      inherits(x, "Matrix") || inherits(x, "matrix") &&
      isTRUE(all.equal(x[upper.tri(x)], rev(x[lower.tri(x)]), tolerance = tol))
  })))
    stop("every relatinships matrix must be a symmetric matrix which identical column and row names")
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
    relfaci <- match(rownames(Zt)[Zti], rownames(covarrel[[i]]))
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
