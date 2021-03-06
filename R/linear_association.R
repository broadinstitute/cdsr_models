require(magrittr)
require(ashr)
require(limma)
require(corpcor)
require(WGCNA)
require(plyr)
require(tidyverse)

#' Calculates the linear association between a matrix of features A and a vector y.
#'
#' @param X n x m numerical matrix of features (missing values are allowed).
#' @param Y n x p numerical matrix of responses (missing values are allowed)
#' @param W n x p numerical matrix of confounders (missing values are not allowed).
#' @param n.min Numeric value controlling the minimum degrees of freedom required for p-values.
#' @param shrinkage Boolean to control whether adaptive shrinkage should be applied.
#' @param alpha Numeric value controling the alpha value for adaptive shrinkage.
#' @param MHC_direction String ("x" or "y") indicating which variable is independent, defaults to matrix with more columns.
#'
#' @return A list with components:
#' \describe{
#'   \item{N}{Vector of degrees of freedom for each feature.}
#'   \item{rho}{Vector of uncorrected beta estimates controlled for W.}
#'   \item{beta}{Vector of corrected beta estimates controlled for W.}
#'   \item{beta.se}{Vector of standard errors of the beta estimates.}
#'   \item{p.val}{Vector of p-values (without adaptive shrinkage).}
#'   \item{q.val}{Vector of q-values (without adaptive shrinkage).}
#'   \item{res.table}{A table with the following columns (values calculated with adaptive shrinkage): \cr
#'                      betahat: estimated linear coefficient controlled for W. \cr
#'                      sebetahat: standard error for the estimate of beta. \cr
#'                      NegativeProb: posterior probabilities that beta is negative. \cr
#'                      PositiveProb: posterior probabilities that beta is positive \cr
#'                      lfsr: local and global FSR values. \cr
#'                      svalue: s-values.\cr
#'                      lfdr: local and global FDR values. \cr
#'                      qvalue: q-values for the effect size estimates. \cr
#'                      PosteriorMean: moderated effect size estimates. \cr
#'                      PosteriorSD: standard deviations of moderated effect size estimates. \cr
#'                      dep.var: dependent variables (column names of Y). \cr
#'                      ind.var: independent variables (column names of X).
#'                   }
#'}
#'
#' @export
#'
lin_associations = function (X,
                             Y,
                             W = NULL,
                             n.min = 4,
                             shrinkage = T,
                             alpha = 0,
                             MHC_direction = NULL) {
  if (is.null(MHC_direction)) {
    MHC_direction = ifelse(length(Y) >= length(X), "x", "y")
  }
  X.NA <- !is.finite(X); X[X.NA] = NA
  Y.NA <- !is.finite(Y); Y[Y.NA] = NA
  N <- (t(!X.NA) %*% (!Y.NA)) - 2
  if (!is.null(W)) {
    P <- W %*% corpcor::pseudoinverse(t(W) %*% W) %*% t(W)
    X[X.NA] <- 0
    X <- X - (P %*% X)
    X[X.NA] <- NA
    Y[Y.NA] <- 0
    Y <- Y - (P %*% Y)
    Y[Y.NA] <- NA
    N = N - dim(W)[2]
  }

  f = function(A) {
    if (is.matrix(A)) {
      res = apply(A, 2, function(x)
        sd(x, na.rm = T))
    } else{
      res = sd(A, na.rm = T)
    }
    return(res)
  }
  sx <- f(X)
  sy <- f(Y)

  rho <- WGCNA::cor(X, Y, use = "pairwise.complete.obs")
  beta <- t(t(rho / sx) * sy)
  beta.se <- t(t(sqrt(1 - rho ^ 2) / sx) * sy) / sqrt(N)
  p.val <- 2 * pt(-abs(beta / beta.se), N)
  p.val[N < n.min] <- NA
  p.val[(sx == 0) | !is.finite(sx), ] <- NA
  p.val[, (sy == 0) | !is.finite(sy)] <- NA

  if (MHC_direction == "y") {
    q.val = apply(p.val, 2, function(x)
      p.adjust(x, method = "BH"))
  } else{
    q.val = t(apply(p.val, 1, function(x)
      p.adjust(x, method = "BH")))
  }
  if (shrinkage) {
    res.table <- list()
    if (MHC_direction == "y") {
      for (ix in 1:dim(p.val)[2]) {
        fin <- is.finite(p.val[, ix])
        res <- ashr::ash(beta[fin, ix],
                         pmax(beta.se[fin, ix], 1e-10),
                         mixcompdist = "halfuniform",
                         alpha = alpha)$result
        if (!is.null(res)) {
          res$dep.var <- colnames(Y)[ix]
          res$ind.var <- rownames(res)
          res$p.val <- p.val[fin, ix]
          res.table[[ix]] <- res
        }
      }
    }
    else {
      for (ix in 1:dim(p.val)[1]) {
        fin <- is.finite(p.val[ix, ])
        res <- tryCatch(
          ashr::ash(
            beta[ix, fin],
            pmax(beta.se[ix, fin], 1e-10),
            mixcompdist = "halfuniform",
            alpha = alpha
          )$result,
          error = function(e)
            NULL
        )
        if (!is.null(res)) {
          res$dep.var <- rownames(res)
          res$ind.var <- colnames(Y)[ix]
          res$p.val <- p.val[ix, fin]
          res.table[[ix]] <- res
        }
      }
    }
    res.table <- dplyr::bind_rows(res.table)
  }
  else {
    res.table <- NULL
  }
  return(
    list(
      N = N,
      rho = rho,
      beta = beta,
      beta.se = beta.se,
      p.val = p.val,
      q.val = q.val,
      res.table = res.table
    )
  )
}




#' Estimate linear-model stats for a matrix of data using limma with empirical Bayes moderated t-stats for p-values
#'
#' @param mat: Nxp data matrix with N cell lines and p genes
#' @param vec: N vector of independent variables. Can be two-group labels as factors, bools, or can be numeric
#' @param covars: Optional Nxk matrix of covariates
#' @param weights: Optional N vector of precision weights for each data point
#' @param target_type: Name of the column variable in the data (default 'Gene')
#' @param limma_trend: Whether to fit an intensity trend with the empirical Bayes variance model
#'
#' @return: data frame of stats
#' @export
#'
#' @examples
#' CRISPR = load.from.taiga(data.name='avana-2-0-1-d98f',
#' data.version=1,
#' data.file='ceres_gene_effects',
#' transpose = T)
#' is_panc <- load.from.taiga(data.name = 'ccle-lines-lineages') %>% .[, 'pancreas']
#' ulines <- intersect(rownames(CRISPR), names(is_panc))
#' lim_res <- run_lm_stats_limma(CRISPR[ulines,], is_panc[ulines])
#' @export run_lm_stats_limma
run_lm_stats_limma <- function(mat, vec, covars = NULL, weights = NULL,
                               target_type = 'Gene', limma_trend = FALSE) {

  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  } else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata,]))
  }
  #there must be more than one unique value of the independent variable
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  #if using covariates add them as additional predictors to the model
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata,, drop = FALSE]
    combined[['pred']] <- pred
    form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  } else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata,])
    } else{
      weights <- weights[udata]
    }
  }

  fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- grep('pred', colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf)

  if (colnames(results)[1] == 'ID') {
    colnames(results)[1] <- target_type
  } else {
    results %<>% rownames_to_column(var = target_type)
  }

  results$min_samples <- min_samples[results[[target_type]]]

  two_to_one_sided <- function(two_sided_p, stat, test_dir) {
    #helper function for converting two-sided p-values to one-sided p-values
    one_sided_p <- two_sided_p / 2
    if (test_dir == 'right') {
      one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
    } else {
      one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
    }
    return(one_sided_p)
  }

  results %<>%
    set_colnames(plyr::revalue(colnames(.), c('logFC' = 'EffectSize', 'AveExpr' = 'Avg',
                                        't' = 't_stat', 'B' = 'log_odds',
                                        'P.Value' = 'p.value', 'adj.P.Val' = 'q.value',
                                        'min_samples' = 'min_samples'))) %>%
    na.omit()
  results %<>%
    dplyr::mutate(p.left = two_to_one_sided(p.value, EffectSize, 'left'),
                  p.right = two_to_one_sided(p.value, EffectSize, 'right'),
                  q.left = p.adjust(p.left, method = 'BH'),
                  q.right = p.adjust(p.right, method = 'BH'))
  return(results)
}

