#' @title Efficiently compute scalar balance statistics
#' @name bal.compute
#' 
#' @description These are functions primarily designed for programmers who want to be able to quickly compute one of several scalar (single number) sample balance statistics, e.g., for use in selecting a tuning parameter when estimating balancing weights. `bal.compute()` computes a scalar balance statistics from the supplied inputs. `bal.init()` initializes the input so that when `bal.compute()` is used on the output along with a set of weights, the computation of the balance statistic is fast. `vignette("optimizing-balance")` provides an overview and more examples of how to use these functions. `available.stats()` returns the balance statistics available for the given treatment type.
#' 
#' @param x for `bal.compute()`, a `bal.init` object created by `bal.init()` or a matrix or data frame containing the covariates. For `bal.init()`, a matrix or data frame containing the covariates.
#' @param weights a vector of balancing weights to compute the weighted statistics.
#' @param treat a vector containing the treatment variable. Can be `NULL` to compute target balance statistics.
#' @param stat string; the name of the statistic to compute. See Details.
#' @param s.weights optional; a vector of sampling weights.
#' @param ... other arguments used to specify options for the balance statistic. See Details for which arguments are allowed with each balance statistic. Ignored for the `bal.init` method of `bal.compute()`.
#' @param treat.type string; the treatment type, either `"binary"`, `"multinomial"`, `"continuous"`, or `"target"`. Abbreviations allowed.
#' 
#' @returns For `bal.compute()`, a single numeric value. For `bal.init()`, a `bal.init` object containing the components created in the initialization and the function used to compute the balance statistic. For `available.stats()`, a character vector of available statistics.
#' 
#' @details 
#' The following list contains the allowable balance statistics that can be supplied to `bal.init()` or the default method of `bal.compute()`, the additional arguments that can be used with each one, and the treatment types allowed with each one. For all balance statistics, lower values indicate better balance. Target balance refers to the similarity between a sample and the same sample but weighted.
#' \describe{
#'     \item{`smd.mean`, `smd.max`, `smd.rms`}{
#'         The mean, maximum, or root-mean-squared absolute standardized mean difference, computed using [col_w_smd()]. The other allowable arguments include `estimand` (`"ATE"`, `"ATT"`, or `"ATC"`) to select the estimand (default is `"ATE"`), `focal` to identify the focal treatment group when the ATT is the estimand and the treatment has more than two categories, and `pairwise` to select whether mean differences should be computed between each pair of treatment groups or between each treatment group and the target group identified by `estimand` (default `TRUE`). Can be used with binary and multi-category treatments and for target balance.
#'     }
#'     \item{`ks.mean`, `ks.max`, `ks.rms`}{
#'         The mean, maximum, or root-mean-squared Kolmogorov-Smirnov statistic, computed using [col_w_ks()]. The other allowable arguments include `estimand` (`"ATE"`, `"ATT"`, or `"ATC"`) to select the estimand (default is `"ATE"`), `focal` to identify the focal treatment group when the ATT is the estimand and the treatment has more than two categories, and `pairwise` to select whether statistics should be computed between each pair of treatment groups or between each treatment group and the target group identified by `estimand` (default `TRUE`). Can be used with binary and multi-category treatments and for target balance.
#'     }
#'     \item{`ovl.mean`, `ovl.max`, `ovl.rms`}{
#'         The mean, maximum, or root-mean-squared overlapping coefficient complement, computed using [col_w_ovl()]. The other allowable arguments include `estimand` (`"ATE"`, `"ATT"`, or `"ATC"`) to select the estimand (default is `"ATE"`), `integrate` to select whether integration is done using using [integrate()] (`TRUE`) or a Riemann sum (`FALSE`, the default), `focal` to identify the focal treatment group when the ATT is the estimand and the treatment has more than two categories, `pairwise` to select whether statistics should be computed between each pair of treatment groups or between each treatment group and the target group identified by `estimand` (default `TRUE`). Can be used with binary and multi-category treatments and for target balance.
#'     }
#'     \item{`mahalanobis`}{
#'         The Mahalanobis distance between the treatment group means. This is similar to `smd.rms` but the covariates are standardized to remove correlations between them and de-emphasize redundant covariates. The other allowable arguments include `estimand` (`"ATE"`, `"ATT"`, or `"ATC"`) to select the estimand (default is `"ATE"`) and `focal` to identify the focal treatment group when the ATT is the estimand. Can only be used with binary treatments and for target balance.
#'     }
#'     \item{`energy.dist`}{
#'         The total energy distance between each treatment group and the target sample, which is a scalar measure of the similarity between two multivariate distributions. The other allowable arguments include `estimand` (`"ATE"`, `"ATT"`, `"ATC"`, or `NULL`) to select the estimand (default is `NULL`), `focal` to identify the focal treatment group when the ATT is the estimand and the treatment has more than two categories, and `improved` to select whether the "improved" energy distance should be used when `estimand = "ATE"`, which emphasizes difference between treatment groups in addition to difference between each treatment group and the target sample (default `TRUE`). When `estimand = NULL`, only the energy distance between the treatment groups will be computed (i.e., as opposed to the energy distance between each treatment groups and the target sample). Can be used with binary and multi-category treatments and for target balance.
#'     }
#'     \item{`kernel.dist`}{
#'         The kernel distance between the treatment groups, which is a scalar measure of the similarity between two multivariate distributions. Can only be used with binary treatments.
#'     }
#'     \item{`l1.med`}{
#'         The median L1 statistic computed across a random selection of possible coarsening of the data. The other allowable arguments include `estimand` (`"ATE"`, `"ATT"`, or `"ATC"`) to select the estimand (default is `"ATE"`), `focal` to identify the focal treatment group when the ATT is the estimand and the treatment has more than two categories, `l1.min.bin` (default 2) and `l1.max.bin` default (12) to select the minimum and maximum number of bins with which to bin continuous variables and `l1.n` (default 101) to select the number of binnings used to select the binning at the median. `covs` should be supplied without splitting factors into dummies to ensure the binning works correctly; for simplicity, the `.covs` argument can be supplied, which will override `covs` but isn't used by other statistics. Can be used with binary and multi-category treatments.
#'     }
#'     \item{`r2`, `r2.2`, `r2.3`}{
#'         The post-weighting \eqn{R^2} of a model for the treatment. The other allowable arguments include `poly` to add polynomial terms of the supplied order to the model and `int` (default `FALSE`) to add two-way interaction between covariates into the model. Using `r2.2` is a shortcut to requesting squares, and using `r2.3` is a shortcut to requesting cubes. Can be used with binary and continuous treatments. For binary treatments, the McKelvey and Zavoina \eqn{R^2} from a logistic regression is used; for continuous treatments, the \eqn{R^2} from a linear regression is used.
#'     }
#'     \item{`p.mean`, `p.max`, `p.rms`}{
#'         The mean, maximum, or root-mean-squared absolute Pearson correlation between the treatment and covariates, computed using [col_w_corr()]. Can only be used with continuous treatments.
#'     }
#'     \item{`s.mean`, `s.max`, `s.rms`}{
#'         The mean, maximum, or root-mean-squared absolute Spearman correlation between the treatment and covariates, computed using [col_w_corr()]. Can only be used with continuous treatments.
#'     }
#'     \item{`distance.cov`, `distance.cor`}{
#'         The distance covariance or distance correlation, respectively, between the scaled covariates and treatment, which is a scalar measure of the independence of two possibly multivariate distributions. The distance correlation is scale-free and ranges from 0 (completely independent) to 1 (perfectly associated). Can only be used with continuous treatments.
#'     }
#' }
#' 
#' Although statistics can be computed directly using `bal.compute()` alone, the intended workflow is to use `bal.init()` to initialize a `bal.init` object, which can then be passed to `bal.compute()` many times with different sets of weights, thereby minimizing the processing that `bal.init()` does because it is only done once. In contrast, using `bal.compute()` on covariates directly (i.e., using the default method) calls `bal.init()` internally each time, which can slow down evaluation. When speed isn't of interest or to calculate a balance statistic outside the context of balance optimization, the default method of `bal.compute()` can be a quick shortcut to avoid having to create a `bal.init` object first.
#' 
#' @seealso [`balance-summary`], [bal.tab()]
#' 
#' See `vignette("optimizing-balance")` for references and definitions of some of the above quantities.
#' 
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE)
#' # Select the optimal number of subclasses for
#' # subclassification:
#' data("lalonde")
#' covs <- c("age", "educ", "race", "married",
#'           "nodegree", "re74", "re75")
#' 
#' # Estimate propensity score
#' p <- glm(reformulate(covs, "treat"),
#'          data = lalonde, 
#'          family = "binomial")$fitted.values
#' 
#' # Function to compute subclassification weights
#' subclass_ATE <- function(treat, p, nsub) {
#'     m <- MatchIt::matchit(treat ~ 1,
#'                           data = lalonde,
#'                           distance = p,
#'                           method = "subclass",
#'                           estimand = "ATE",
#'                           subclass = nsub)
#'     m$weights
#' }
#' 
#' # Initialize balance statistic; largest KS statistic
#' init <- bal.init(lalonde[covs],
#'                  treat = lalonde$treat, 
#'                  stat = "ks.max",
#'                  estimand = "ATE")
#' 
#' # Statistic prior to subclassification:
#' bal.compute(init)
#' 
#' # Testing 4 to 50 subclasses
#' nsubs <- 4:50
#' stats <- vapply(nsubs, function(n) {
#'     w <- subclass_ATE(lalonde$treat, p, n)
#'     bal.compute(init, w)
#' }, numeric(1L))
#' 
#' plot(stats ~ nsubs)
#' 
#' # 6 subclass gives lowest ks.max value (.238)
#' nsubs[which.min(stats)]
#' stats[which.min(stats)]
#' 
#' # See which statistics are available
#' available.stats("binary")
#' available.stats("multinomial")
#' available.stats("continuous")
#' available.stats("target")


#' @export
bal.compute <- function(x, ...) {
  UseMethod("bal.compute")
}

#' @rdname bal.compute
#' @exportS3Method bal.compute bal.init
bal.compute.bal.init <- function(x,
                                 weights = NULL,
                                 ...) {
  fun <- attr(x, "fun", TRUE)
  
  fun(init = x, weights = weights)
}

#' @rdname bal.compute
#' @exportS3Method bal.compute default
bal.compute.default <- function(x,
                                treat = NULL,
                                stat,
                                s.weights = NULL,
                                weights = NULL,
                                ...) {
  .chk_not_missing(x, "`x`")
  .chk_not_missing(stat, "`stat`")
  
  bal.init(x = x, treat = treat, stat = stat, s.weights = s.weights, ...) |>
    bal.compute.bal.init(weights = weights)
}

#' @rdname bal.compute
#' @export
bal.init <- function(x,
                     treat = NULL,
                     stat,
                     s.weights = NULL,
                     ...) {
  .chk_not_missing(x, "`x`")
  .chk_not_missing(stat, "`stat`")
  
  .chk_string(stat)
  
  if (is_null(treat)) {
    treat.type <- "target"
  }
  else {
    .chk_vector(treat)
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
  }
  
  stat <- match_arg(stat, available.stats(treat.type))
  
  .chk_null_or(s.weights, chk = .chk_numeric)
  
  init <- bal_criterion(treat.type, stat)
  
  out <- init$init(x = x, treat = treat, s.weights = s.weights, ...)
  
  attr(out, "fun") <- init$fun
  attr(out, "treat.type") <- treat.type
  attr(out, "stat") <- stat
  
  set_class(out, "bal.init", .replace = FALSE, .last = FALSE)
}

#' @rdname bal.compute
#' @export
available.stats <- function(treat.type = "binary") {
  .chk_string(treat.type)
  if (identical(treat.type, "multi-category")) {
    treat.type <- "multinomial"
  }
  treat.type <- match_arg(treat.type, c("binary", "multinomial", "continuous", "target"))
  
  switch(
    treat.type,
    binary = c(
      "smd.mean",
      "smd.max",
      "smd.rms",
      "ks.mean",
      "ks.max",
      "ks.rms",
      "ovl.mean",
      "ovl.max",
      "ovl.rms",
      "mahalanobis",
      "energy.dist",
      "kernel.dist",
      "l1.med",
      "r2",
      "r2.2",
      "r2.3"
    ),
    multinomial = c(
      "smd.mean",
      "smd.max",
      "smd.rms",
      "ks.mean",
      "ks.max",
      "ks.rms",
      "ovl.mean",
      "ovl.max",
      "ovl.rms",
      "energy.dist",
      "l1.med"
    ),
    continuous = c(
      "p.mean",
      "p.max",
      "p.rms",
      "s.mean",
      "s.max",
      "s.rms",
      "r2",
      "r2.2",
      "r2.3",
      "distance.cov",
      "distance.cor"
    ),
    target = c(
      "smd.mean",
      "smd.max",
      "smd.rms",
      "ks.mean",
      "ks.max",
      "ks.rms",
      "ovl.mean",
      "ovl.max",
      "ovl.rms",
      "mahalanobis",
      "energy.dist"
    )
  )
}

#' @exportS3Method print bal.init
print.bal.init <- function(x, ...) {
  cat("A `bal.init` object\n")
  cat(sprintf("  treatment type: %s\n", attr(x, "treat.type")))
  cat(sprintf("  statistic: %s (%s)\n", attr(x, "stat"), bal_stat.to.phrase(attr(x, "stat"))))
  cat("Use `bal.compute()` to compute the balance statistic.\n")
  invisible(x)
}

bal_stat.to.phrase <- function(stat) {
  phrase <- switch(stat,
                   "smd.mean" = "average absolute standardized mean difference",
                   "smd.max" = "maximum absolute standardized mean difference",
                   "smd.rms" = "root-mean-square absolute standardized mean difference",
                   "ks.mean" = "average Kolmogorov-Smirnov statistic",
                   "ks.max" = "maximum Kolmogorov-Smirnov statistic",
                   "ks.rms" = "root-mean-square Kolmogorov-Smirnov statistic",
                   "ovl.mean" = "average overlapping coefficient complement",
                   "ovl.max" = "maximum overlapping coefficient complement",
                   "ovl.rms" = "root-mean-square overlapping coefficient complement",
                   "mahalanobis" = "sample Mahalanobis distance",
                   "energy.dist" = "energy distance",
                   "kernel.dist" = "kernel distance",
                   "l1.med" = "L1 median",
                   "r2" = "post-weighting treatment R-squared",
                   "p.mean" = "average Pearson correlation",
                   "p.max" = "maximum Pearson correlation",
                   "p.rms" = "root-mean-square Pearson correlation",
                   "s.mean" = "average Spearman correlation",
                   "s.max" = "maximum Spearman correlation",
                   "s.rms" = "root-mean-square Spearman correlation",
                   "distance.cov" = "distance covariance",
                   "distance.cor" = "distance correlation",
                   NA_character_
  )
  
  if (anyNA(phrase)) {
    .err(sprintf("%s is not an allowed statistic",
                 add_quotes(stat, 2L)))
  }
  
  phrase
}

process_init_covs <- function(covs) {
  nm <- deparse1(substitute(covs))
  
  if ((!is.numeric(covs) && !is.data.frame(covs)) || length(dim(covs)) > 2L) {
    .err(sprintf("`%s` must be a data.frame or numeric matrix", nm))
  }
  
  needs.splitting <- FALSE
  
  if (is.data.frame(covs)) {
    to.split <- vapply(covs, chk::vld_character_or_factor, logical(1L))
    
    if (any(to.split)) {
      needs.splitting <- TRUE
    }
    else {
      covs <- as.matrix(covs)
    }
  }
  else if (!is.matrix(covs)) {
    covs <- matrix(covs, ncol = 1L)
  }
  
  bin.vars <- if_null_then(attr(covs, "bin"), .process_bin_vars(mat = covs))
  
  if (needs.splitting) {
    bin.vars[to.split] <- TRUE
    covs <- splitfactor(covs, drop.first = "if2", split.with = bin.vars)
    bin.vars <- attr(covs, "split.with")[[1L]]
  }
  
  covs <- as.matrix(covs)
  
  attr(covs, "bin") <- bin.vars
  
  covs
}

#init functions
init_smd <- function(x, treat = NULL, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, ...) {
  
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  if (is_null(treat)) {
    n <- length(s.weights)
    
    id2 <- rep.int(seq_len(n), 2L)
    
    treat <- rep(0:1, each = n)
    s.weights <- s.weights[id2]
    x <- x[id2, , drop = FALSE]
    
    denoms <- .compute_s.d.denom(x, treat = treat,
                                 s.d.denom = "0",
                                 s.weights = s.weights,
                                 bin.vars = bin.vars)
    
    focal <- pairwise <- treatment.pairs <- NULL
  }
  else {
    .chk_atomic(treat)
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
      .err("`treat` must be a binary or multi-category variable")
    }
    
    .chk_null_or(estimand, .chk_string)
    if (is_null(estimand)) estimand <- "ATE"
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    s.d.denom <- .get_s.d.denom(estimand = estimand, treat = treat,
                                focal = focal, quietly = TRUE)
    
    denoms <- .compute_s.d.denom(x, treat = treat,
                                 s.d.denom = s.d.denom,
                                 s.weights = s.weights,
                                 bin.vars = bin.vars)
    
    if (treat.type == "binary" && is_not_null(focal)) {
      pairwise <- TRUE
    }
    
    .chk_flag(pairwise)
    
    unique.treats <- {
      if (is.factor(treat)) levels(droplevels(treat))
      else unique(treat)
    }
    
    if (pairwise) {
      treatment.pairs <- {
        if (is_null(focal))
          utils::combn(unique.treats, 2L, simplify = FALSE)
        else 
          lapply(setdiff(unique.treats, focal), c, focal)
      }
    }
    else {
      treat.all <- last(make.unique(c(unique.treats, "All")))
      treat <- factor(c(as.character(treat), rep_with(treat.all, treat)),
                      levels = c(unique.treats, treat.all))
      x <- x[rep(seq_row(x), 2L), , drop = FALSE]
      s.weights <- rep.int(s.weights, 2L)
      
      treatment.pairs <- lapply(unique.treats, c, treat.all)
    }
  }
  
  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              s.d.denom = denoms,
              focal = focal,
              pairwise = pairwise,
              treatment.pairs = treatment.pairs)
  
  set_class(out, "init_smd")
}
init_ks <- function(x, treat = NULL, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, ...) {
  
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  if (is_null(treat)) {
    n <- length(s.weights)
    
    id2 <- rep.int(seq_len(n), 2L)
    
    treat <- rep(0:1, each = n)
    s.weights <- s.weights[id2]
    x <- x[id2, , drop = FALSE]
    
    focal <- pairwise <- treatment.pairs <- NULL
  }
  else {
    .chk_atomic(treat)
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
      .err("`treat` must be a binary or multi-category variable")
    }
    
    .chk_null_or(estimand, .chk_string)
    if (is_null(estimand)) {
      estimand <- "ATE"
    }
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    if (treat.type == "binary" && is_not_null(focal)) {
      pairwise <- TRUE
    }
    
    .chk_flag(pairwise)
    
    unique.treats <- {
      if (is.factor(treat)) levels(droplevels(treat))
      else unique(treat)
    }
    
    if (!pairwise) {
      treat.all <- last(make.unique(c(unique.treats, "All")))
      treat <- factor(c(as.character(treat), rep_with(treat.all, treat)),
                      levels = c(unique.treats, treat.all))
      x <- x[rep(seq_row(x), 2L), , drop = FALSE]
      s.weights <- rep.int(s.weights, 2L)
      focal <- treat.all
      
      treatment.pairs <- lapply(unique.treats, c, treat.all)
    }
    else if (treat.type == "multinomial") {
      treatment.pairs <- {
        if (is_null(focal))
          utils::combn(unique.treats, 2L, simplify = FALSE)
        else 
          lapply(setdiff(unique.treats, focal), c, focal)
      }
    }
    else {
      treatment.pairs <- list(unique.treats)
    }
  }
  
  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              focal = focal,
              pairwise = pairwise,
              treatment.pairs = treatment.pairs)
  
  set_class(out, "init_ks")
}
init_ovl <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE,
                     integrate = TRUE, ...) {
  
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  .chk_flag(integrate)
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  if (is_null(treat)) {
    n <- length(s.weights)
    
    id2 <- rep.int(seq_len(n), 2L)
    
    treat <- rep(0:1, each = n)
    s.weights <- s.weights[id2]
    x <- x[id2, , drop = FALSE]
    
    focal <- pairwise <- treatment.pairs <- NULL
  }
  else {
    .chk_atomic(treat)
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
      .err("`treat` must be a binary or multi-category variable")
    }
    
    .chk_null_or(estimand, .chk_string)
    if (is_null(estimand)) estimand <- "ATE"
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    if (treat.type == "binary" && is_not_null(focal)) {
      pairwise <- TRUE
    }
    
    .chk_flag(pairwise)
    
    unique.treats <- {
      if (is.factor(treat)) levels(droplevels(treat))
      else unique(treat)
    }
    
    if (!pairwise) {
      treat.all <- last(make.unique(c(unique.treats, "All")))
      treat <- factor(c(as.character(treat), rep_with(treat.all, treat)),
                      levels = c(unique.treats, treat.all))
      x <- x[rep(seq_row(x), 2L), , drop = FALSE]
      s.weights <- rep.int(s.weights, 2L)
      focal <- treat.all
      
      treatment.pairs <- lapply(unique.treats, c, treat.all)
    }
    else if (treat.type == "multinomial") {
      treatment.pairs <- {
        if (is_null(focal))
          utils::combn(unique.treats, 2L, simplify = FALSE)
        else 
          lapply(setdiff(unique.treats, focal), c, focal)
      }
    }
    else {
      treatment.pairs <- list(unique.treats)
    }
  }
  
  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              focal = focal,
              pairwise = pairwise,
              treatment.pairs = treatment.pairs,
              integrate = integrate)
  
  set_class(out, "init_ovl")
}
init_mahalanobis <- function(x, treat = NULL, s.weights = NULL, estimand = NULL, focal = NULL, ...) {
  
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  if (anyNA(x)) {
    .err('"mahalanobis" cannot be used when there are missing values in the covariates')
  }
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  if (is_null(treat)) {
    sigma <- cov.wt(x, s.weights)$cov
    
    n <- length(s.weights)
    
    id2 <- rep.int(seq_len(n), 2L)
    
    treat <- rep(0:1, each = n)
    s.weights <- s.weights[id2]
    x <- x[id2, , drop = FALSE]
    
    s.d.denom <- NULL
  }
  else {
    .chk_atomic(treat)
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
      .err("`treat` must be a binary or multi-category variable")
    }
    
    .chk_null_or(estimand, .chk_string)
    if (is_null(estimand)) estimand <- "ATE"
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    s.d.denom <- .get_s.d.denom(estimand = estimand, treat = treat, focal = focal, quietly = TRUE)
    
    if (!all(bin.vars)) {
      x[, !bin.vars] <- scale(x[, !bin.vars])
    }
    
    if (s.d.denom %in% as.character(treat)) {
      sigma <- cov.wt(x[treat == s.d.denom, , drop = FALSE], s.weights[treat == s.d.denom])$cov
      zeros <- diag(sigma) == 0
      if (any(zeros)) {
        sigma_all <- cov.wt(x, s.weights)$cov
        sigma[zeros, ] <- sigma_all[zeros, ]
        sigma[, zeros] <- sigma_all[, zeros]
        
      }
    }
    else if (s.d.denom == "pooled") {
      sigma <- .5 * (cov.wt(x[treat == treat[1L], , drop = FALSE], s.weights[treat == treat[1L]])$cov +
                       cov.wt(x[treat != treat[1L], , drop = FALSE], s.weights[treat != treat[1L]])$cov)
    }
    else {
      sigma <- cov.wt(x, s.weights)$cov
    }
  }
  
  zeros <- diag(sigma) == 0
  if (any(zeros)) {
    diag(sigma)[zeros] <- 1
  }
  
  #MASS::ginv
  sigmasvd <- svd(sigma)
  pos <- sigmasvd$d > max(1e-8 * sigmasvd$d[1L], 0)
  sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% ((1 / sigmasvd$d[pos]) * t(sigmasvd$u[, pos, drop = FALSE]))
  
  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              s.d.denom = s.d.denom,
              sigma_inv = sigma_inv)
  
  set_class(out, "init_mahalanobis")
}
init_energy.dist <- function(x, treat = NULL, s.weights = NULL, estimand = NULL, focal = NULL, improved = TRUE, ...) {
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  if (anyNA(x)) {
    .err('"energy.dist" cannot be used when there are missing values in the covariates')
  }
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  s.weights_n <- s.weights / sum(s.weights)
  
  if (is_null(treat)) {
    dist.covs <- scale(x, scale = sqrt(col.w.v(x, s.weights, bin.vars)))
    
    d <- unname(as.matrix(dist(dist.covs)))
    
    n <- length(s.weights)
    
    P <- -d / (n^2)
    q <- drop(s.weights_n %*% d) * 2 / n
    k <- -sum(q * s.weights_n) * n / 2
    
    unique.treats <- focal <- NULL
  }
  else {
    .chk_atomic(treat)
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
      .err("`treat` must be a binary or multi-category variable")
    }
    
    treat <- factor(treat)
    
    if (is_not_null(estimand)) {
      .chk_string(estimand)
      f.e <- process_focal_and_estimand(focal, estimand, treat)
      focal <- f.e[["focal"]]
      estimand <- f.e[["estimand"]]
    }
    
    dist.covs <- scale(x, scale = sqrt(col.w.v(x, s.weights, bin.vars)))
    
    d <- unname(as.matrix(dist(dist.covs)))
    
    n <- length(treat)
    unique.treats <- levels(treat)
    
    treat_t <- vapply(unique.treats, function(t) treat == t, logical(n))
    n_t <- colSums(treat_t)
    
    treat_n_t <- vapply(unique.treats, function(t) treat_t[, t] / n_t[t],
                        numeric(n))
    
    if (is_null(estimand)) {
      .col_diff <- function(x) x[, 1L] - x[, 2L]
      all_pairs <- utils::combn(unique.treats, 2L, simplify = FALSE)
      nn <- tcrossprod(vapply(all_pairs, function(p) .col_diff(treat_n_t[, p, drop = FALSE]),
                              numeric(n)))
      
      P <- -d * nn
      q <- rep.int(0, n)
      k <- 0
    }
    else if (is_null(focal)) {
      .chk_flag(improved)
      
      nn <- tcrossprod(treat_n_t)
      
      if (improved) {
        .col_diff <- function(x) x[, 1L] - x[, 2L]
        all_pairs <- utils::combn(unique.treats, 2L, simplify = FALSE)
        nn <- nn + tcrossprod(vapply(all_pairs, function(p) .col_diff(treat_n_t[, p, drop = FALSE]),
                                     numeric(n)))
      }
      
      P <- -d * nn
      q <- drop(s.weights_n %*% d) * 2 * rowSums(treat_n_t)
      k <- - length(n_t) * drop(s.weights_n %*% d %*% s.weights_n)
    }
    else {
      nn <- tcrossprod(treat_n_t)
      in_focal <- treat == focal
      
      s.weights_n_f <- s.weights[in_focal]
      s.weights_n_f <- s.weights_n_f / sum(s.weights_n_f)
      
      P <- -d * nn
      q <- drop(s.weights_n_f %*% d[in_focal, , drop = FALSE]) * 2 * rowSums(treat_n_t)
      k <- -length(n_t) * drop(s.weights_n_f %*% d[in_focal, in_focal, drop = FALSE] %*% s.weights_n_f)
    }
  }
  
  out <- list(P = P,
              q = q,
              k = k,
              s.weights = s.weights,
              treat = treat,
              unique.treats = unique.treats)
  
  set_class(out, "init_energy.dist")
}
init_kernel.dist <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL, ...) {
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  if (anyNA(x)) {
    .err('"kernel.dist" cannot be used when there are missing values in the covariates')
  }
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)
  
  if (treat.type %nin% c("binary")) {
    .err("`treat` must be a binary variable")
  }
  
  treat <- as.numeric(treat == treat[1L])
  
  for (t in 0:1) {
    s.weights[treat == t] <- s.weights[treat == t] / mean_fast(s.weights[treat == t])
  }
  
  dist.covs <- scale(x, scale = sqrt(col.w.v(x, s.weights, bin.vars)))
  
  d <- unname(as.matrix(dist(dist.covs)))
  
  K <- exp(-(d^2) / median(d))
  
  T_star <- numeric(length(treat))
  T_star[treat == 1] <- 1 / sum(treat == 1)
  T_star[treat == 0] <- -1 / sum(treat == 0)
  
  out <- list(K = K,
              T_star = T_star,
              s.weights = s.weights,
              treat = treat)
  
  set_class(out, "init_kernel.dist")
}
init_p <- function(x, treat, s.weights = NULL, ...) {
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)
  
  if (treat.type %nin% c("continuous")) {
    .err("`treat` must be a continuous (numeric) variable")
  }
  
  s.d.denom <- .get_s.d.denom.cont(quietly = TRUE)
  
  denoms <- .compute_s.d.denom(x, treat = treat,
                               s.d.denom = s.d.denom, s.weights = s.weights,
                               bin.vars = bin.vars)
  
  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              s.d.denom = denoms)
  
  set_class(out, "init_p")
}
init_s <- function(x, treat, s.weights = NULL, ...) {
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)
  
  if (treat.type %nin% c("continuous")) {
    .err("`treat` must be a continuous (numeric) variable")
  }
  
  for (i in seq_col(x)[!bin.vars[i]]) {
    x[, i] <- rank(x[, i], na.last = "keep")
  }
  treat <- rank(treat, na.last = "keep")
  
  s.d.denom <- .get_s.d.denom.cont(quietly = TRUE)
  
  denoms <- .compute_s.d.denom(x, treat = treat,
                               s.d.denom = s.d.denom, s.weights = s.weights,
                               bin.vars = bin.vars)
  
  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              s.d.denom = denoms)
  
  set_class(out, "init_s")
}
init_r2 <- function(x, treat, s.weights = NULL, poly = 1, int = FALSE, ...) {
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  
  if (anyNA(x)) {
    .err('"r2" cannot be used when there are missing values in the covariates')
  }
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)
  
  if (treat.type %nin% c("binary", "continuous")) {
    .err("`treat` must be a binary or continuous (numeric) variable")
  }
  
  if (treat.type == "binary") {
    treat <- as.numeric(treat == treat[1L])
  }
  
  x <- cbind(`(Intercept)` = 1, x, .int_poly_f2(x, poly = poly, int = int))
  
  out <- list(treat = treat,
              x = x,
              s.weights = s.weights)
  
  set_class(out, "init_r2")
}
init_distance.cov <- function(x, treat, s.weights = NULL, std = FALSE, ...) {
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")
  
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  
  if (anyNA(x)) {
    .err(sprintf("%s cannot be used when there are missing values in the covariates",
                 add_quotes(if (std) "distance.cor" else "distance.cov", 2L)))
  }
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  s.weights <- s.weights / sum(s.weights)
  
  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)
  
  if (treat.type %nin% c("continuous")) {
    .err("`treat` must be a continuous (numeric) variable")
  }
  
  dist.covs <- scale(x, scale = sqrt(col.w.v(x, s.weights, bin.vars)))
  
  Xdist <- dist(dist.covs) |>
    as.matrix() |>
    unname()
  
  n <- length(treat)
  
  Adist <- dist(treat / sqrt(col.w.v(treat, s.weights))) |>
    as.matrix() |>
    unname()
  
  Xmeans <- colMeans(Xdist)
  Xgrand_mean <- mean(Xmeans)
  XX <- Xdist + Xgrand_mean - outer(Xmeans, Xmeans, "+")
  
  Ameans <- colMeans(Adist)
  Agrand_mean <- mean(Ameans)
  AA <- Adist + Agrand_mean - outer(Ameans, Ameans, "+")
  
  if (std) {
    Adenom <- sqrt(drop(t(s.weights) %*% (AA^2) %*% s.weights))
    Xdenom <- sqrt(drop(t(s.weights) %*% (XX^2) %*% s.weights))
    denom <- Adenom * Xdenom
  }
  else {
    denom <- 1
  }
  
  out <- list(P = XX * AA,
              denom = denom,
              s.weights = s.weights,
              treat = treat)
  
  set_class(out, "init_distance.cov")
}
init_l1.med <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL,
                        .covs = NULL, l1.min.bin = 2, l1.max.bin = 12, l1.n = 101, ...) {
  
  if (is_not_null(.covs)) {
    x <- .covs
  }
  
  if (!is.data.frame(x)) {
    if (is.atomic(x) && is_null(dim(x))) {
      x <- data.frame(x)
    }
    else if (!is.matrix(x)) {
      .err("`x` must be a data.frame or matrix")
    }
  }
  x <- as.data.frame(x)
  
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  
  if (anyNA(x)) {
    .err('"l1.med" cannot be used when there are missing values in the covariates')
  }
  
  check_arg_lengths(x, treat, s.weights)
  
  if (is_null(s.weights)) {
    s.weights <- rep.int(1, NROW(x))
  }
  
  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)
  
  if (treat.type %nin% c("binary", "multinomial")) {
    .err("`treat` must be a binary or multi-category variable")
  }
  
  coarsen <- function(covs, cutpoints = NULL, grouping = NULL) {
    is.numeric.cov <- setNames(vapply(covs, is.numeric, logical(1L)), names(covs))
    for (i in names(cutpoints)) {
      if (cutpoints[[i]] == 0) {
        is.numeric.cov[i] <- FALSE #Will not be binned
      }
    }
    
    #Process grouping
    if (is_not_null(grouping) && is_not_null(names(grouping))) {
      covs[names(grouping)] <- lapply(names(grouping), function(g) {
        x <- covs[[g]]
        groups <- grouping[[g]]
        
        for (i in seq_along(groups)) {
          x[x %in% groups[[i]]] <- groups[[i]][1L]
        }
        x
      })
      cutpoints[names(grouping)] <- NULL
    }
    
    #Create bins for numeric variables
    for (i in names(covs)[is.numeric.cov]) {
      bins <- cutpoints[[i]]
      
      #cutpoints is number of bins, unlike in cem
      breaks <- seq(min(covs[[i]]), max(covs[[i]]), length.out = bins + 1L)
      breaks[c(1L, bins + 1L)] <- c(-Inf, Inf)
      
      covs[[i]] <- findInterval(covs[[i]], breaks)
    }
    
    #Reduce to strata
    factor(do.call("paste", c(covs, sep = " | ")))
  }
  
  .chk_null_or(estimand, .chk_string)
  if (is_null(estimand)) {
    estimand <- "ATE"
  }
  
  f.e <- process_focal_and_estimand(focal, estimand, treat)
  focal <- f.e[["focal"]]
  estimand <- f.e[["estimand"]]
  
  unique.treats <- unique(treat)
  for (t in unique.treats) {
    s.weights[treat == t] <- s.weights[treat == t] / sum(s.weights[treat == t])
  }
  
  is.numeric.cov <- setNames(vapply(x, is.numeric, logical(1L)), names(x))
  nunique.covs <- vapply(x, nunique, integer(1L))
  
  coarsenings <- lapply(seq_len(l1.n), function(i) {
    cutpoints <- setNames(lapply(nunique.covs[is.numeric.cov], function(nu) {
      sample(seq(min(l1.min.bin, nu), min(l1.max.bin, nu)), 1L)
    }), names(x)[is.numeric.cov])
    
    grouping <- setNames(lapply(seq_along(x)[!is.numeric.cov], function(i) {
      nu <- nunique.covs[i]
      u <- unique(x[[i]], nmax = nu)
      
      #Randomly select number of bins
      nbins <- sample(seq(min(l1.min.bin, nu), min(l1.max.bin, nu)), 1L)
      
      #Randomly assign bin numbers to levels of covariate
      bin.assignments <- sample.int(nbins, nu, replace = TRUE)
      
      #Group levels with same bin number
      lapply(unique(bin.assignments, nmax = nbins),
             function(b) u[bin.assignments == b])
      
    }), names(x)[!is.numeric.cov])
    
    list(cutpoints = cutpoints, grouping = grouping, treat_cutpoints = NULL)
  })
  
  l1s <- unlist(lapply(coarsenings, function(co) {
    x <- coarsen(x, cutpoints = co[["cutpoints"]], grouping = co[["grouping"]])
    
    if (treat.type == "binary" || is_null(focal)) {
      sum(vapply(levels(x), function(l) {
        in_l <- which(x == l)
        abs(diff(range(vapply(unique.treats, function(t) sum(s.weights[in_l][treat[in_l] == t]),
                              numeric(1L)))))
      }, numeric(1L))) / length(unique.treats)
    }
    else {
      sum(vapply(levels(x), function(l) {
        in_l <- which(x == l)
        sum.s.weights.focal <- sum(s.weights[in_l][treat[in_l] == focal])
        max(abs(vapply(unique.treats[unique.treats != focal],
                       function(t) sum(s.weights[in_l][treat[in_l] == t]) - sum.s.weights.focal,
                       numeric(1L))))
      }, numeric(1L))) / length(unique.treats)
    }
  }))
  
  l1.med <- sort(l1s, partial = ceiling(l1.n / 2))[ceiling(l1.n / 2)]
  
  l1.med.coarsening <- coarsenings[[which(l1s == l1.med)[1L]]]
  
  out <- list(coarsened.covs = coarsen(x,
                                       cutpoints = l1.med.coarsening[["cutpoints"]],
                                       grouping = l1.med.coarsening[["grouping"]]),
              s.weights = s.weights,
              treat = treat,
              unique.treats = unique.treats,
              focal = focal)
  
  set_class(out, "init_l1.med")
}

# Statistics
## binary
smd.binary <- function(init, weights = NULL) {
  .check_init(init, "init_smd")
  
  if (!init$pairwise) {
    return(smd.multinomial(init, weights))
  }
  
  col_w_smd(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
            bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE)
}
ks.binary <- function(init, weights = NULL) {
  .check_init(init, "init_ks")
  
  if (!init$pairwise) {
    return(ks.multinomial(init, weights))
  }
  
  col_w_ks(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
           bin.vars = init$bin.vars)
}
ovl.binary <- function(init, weights = NULL) {
  .check_init(init, "init_ovl")
  
  if (!init$pairwise) {
    return(ovl.multinomial(init, weights))
  }
  
  col_w_ovl(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
            bin.vars = init$bin.vars, integrate = init$integrate)
}
mahalanobis.binary <- function(init, weights = NULL) {
  .check_init(init, "init_mahalanobis")
  mean.diffs <- col_w_smd(init$covs, init$treat,
                          weights = weights,
                          s.weights = init$s.weights,
                          bin.vars = init$bin.vars,
                          std = FALSE, abs = TRUE)
  drop(sqrt(t(mean.diffs) %*% init$sigma_inv %*% mean.diffs))
}
energy.dist.binary <- function(init, weights = NULL) {
  .check_init(init, "init_energy.dist")
  
  weights <- {
    if (is_null(weights)) init[["s.weights"]]
    else weights * init[["s.weights"]]
  }
  
  for (t in init[["unique.treats"]]) {
    weights[init[["treat"]] == t] <- weights[init[["treat"]] == t] / mean_fast(weights[init[["treat"]] == t])
  }
  
  drop(weights %*% init[["P"]] %*% weights) + sum(init[["q"]] * weights) + init[["k"]]
}
kernel.dist.binary <- function(init, weights = NULL) {
  .check_init(init, "init_kernel.dist")
  
  weights <- {
    if (is_null(weights)) init[["s.weights"]]
    else weights * init[["s.weights"]]
  }
  
  for (t in 0:1) {
    weights[init[["treat"]] == t] <- weights[init[["treat"]] == t] / mean_fast(weights[init[["treat"]] == t])
  }
  
  T_weights <- init[["T_star"]] * weights
  
  drop(sqrt(T_weights %*% init[["K"]] %*% T_weights))
}
r2.binary <- function(init, weights = NULL) {
  .check_init(init, "init_r2")
  
  weights <- {
    if (is_null(weights)) init[["s.weights"]]
    else weights * init[["s.weights"]]
  }
  
  fit <- glm.fit(x = init$x, y = init$treat, weights = weights,
                 family = quasibinomial())
  
  wmtreat <- w.m(fit$linear.predictors, weights)
  
  SSmodel <- sum(weights * (fit$linear.predictors - wmtreat)^2)
  
  SSmodel / (sum(weights) * pi^2 / 3 + SSmodel)
}
l1.med.binary <- function(init, weights = NULL) {
  .check_init(init, "init_l1.med")
  
  weights <- {
    if (is_null(weights)) init[["s.weights"]]
    else weights * init[["s.weights"]]
  }
  
  for (t in init[["unique.treats"]]) {
    weights[init[["treat"]] == t] <- weights[init[["treat"]] == t] / sum(weights[init[["treat"]] == t])
  }
  
  x <- init[["coarsened.covs"]]
  
  sum(vapply(levels(x), function(l) {
    in_l <- which(x == l)
    abs(diff(vapply(init[["unique.treats"]],
                    function(t) sum(weights[in_l][init[["treat"]][in_l] == t]),
                    numeric(1L))))
  }, numeric(1L))) / 2
  
}
## multinomial
smd.multinomial <- function(init, weights = NULL) {
  .check_init(init, "init_smd")
  
  if (!init$pairwise && is_not_null(weights)) {
    weights <- c(weights, rep_with(1, weights))
  }
  
  unlist(lapply(init$treatment.pairs, function(x) {
    t_in_x <- which(init$treat %in% x)
    
    col_w_smd(init$covs[t_in_x, , drop = FALSE],
              treat = as.integer(init$treat[t_in_x] == x[1L]),
              weights = weights[t_in_x],
              s.weights = init$s.weights[t_in_x],
              bin.vars = init$bin.vars,
              s.d.denom = init$s.d.denom,
              abs = TRUE)
  }))
}
ks.multinomial <- function(init, weights = NULL) {
  .check_init(init, "init_ks")
  
  if (!init$pairwise && is_not_null(weights)) {
    weights <- c(weights, rep_with(1, weights))
  }
  
  unlist(lapply(init$treatment.pairs, function(x) {
    t_in_x <- which(init$treat %in% x)
    
    col_w_ks(init$covs[t_in_x, , drop = FALSE],
             treat = as.integer(init$treat[t_in_x] == x[1L]),
             weights = weights[t_in_x],
             s.weights = init$s.weights[t_in_x],
             bin.vars = init$bin.vars)
  }))
}
ovl.multinomial <- function(init, weights = NULL) {
  .check_init(init, "init_ovl")
  
  if (!init$pairwise && is_not_null(weights)) {
    weights <- c(weights, rep_with(1, weights))
  }
  
  unlist(lapply(init$treatment.pairs, function(x) {
    t_in_x <- which(init$treat %in% x)
    
    col_w_ovl(init$covs[t_in_x, , drop = FALSE],
              treat = as.integer(init$treat[t_in_x] == x[1L]),
              weights = weights[t_in_x],
              s.weights = init$s.weights[t_in_x],
              bin.vars = init$bin.vars,
              integrate = init$integrate)
  }))
}
energy.dist.multinomial <- function(init, weights = NULL) {
  energy.dist.binary(init, weights)
}
l1.med.multinomial <- function(init, weights = NULL) {
  .check_init(init, "init_l1.med")
  
  weights <- {
    if (is_null(weights)) init[["s.weights"]]
    else weights * init[["s.weights"]]
  }
  
  for (t in init[["unique.treats"]]) {
    weights[init[["treat"]] == t] <- weights[init[["treat"]] == t] / sum(weights[init[["treat"]] == t])
  }
  
  x <- init[["coarsened.covs"]]
  
  if (is_null(init[["focal"]])) {
    sum(vapply(levels(x), function(l) {
      in_l <- which(x == l)
      abs(diff(range(vapply(init[["unique.treats"]], function(t) {
        sum(weights[in_l][init[["treat"]][in_l] == t])
      }, numeric(1L)))))
    }, numeric(1L))) / length(init[["unique.treats"]])
  }
  else {
    sum(vapply(levels(x), function(l) {
      in_l <- which(x == l)
      sum.weights.focal <- sum(weights[in_l][init[["treat"]][in_l] == init[["focal"]]])
      max(abs(vapply(init[["unique.treats"]][init[["unique.treats"]] != init[["focal"]]],
                     function(t) {
                       sum(weights[in_l][init[["treat"]][in_l] == t]) - sum.weights.focal
                     }, numeric(1L))))
    }, numeric(1L))) / length(init[["unique.treats"]])
  }
  
}
## continuous
pearson.corr.continuous <- function(init, weights = NULL) {
  .check_init(init, "init_p")
  col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
            bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
            std = TRUE)
}
spearman.corr.continuous <- function(init, weights = NULL) {
  .check_init(init, "init_s")
  col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
            bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
            std = TRUE)
}
r2.continuous <- function(init, weights = NULL) {
  .check_init(init, "init_r2")
  
  weights <- {
    if (is_null(weights)) init[["s.weights"]]
    else weights * init[["s.weights"]]
  }
  
  fit <- lm.wfit(x = init$x, y = init$treat, w = weights)
  
  SSresid <- sum(weights * (fit$residuals - w.m(fit$residuals, weights))^2)
  SStreat <- sum(weights * (init$treat - w.m(init$treat, weights))^2)
  
  1 - SSresid / SStreat
}
distance.cov.continuous <- function(init, weights = NULL) {
  .check_init(init, "init_distance.cov")
  
  weights <- {
    if (is_null(weights)) init[["s.weights"]]
    else weights * init[["s.weights"]]
  }
  
  weights <- weights / sum(weights)
  
  dcov <- drop(t(weights) %*% init[["P"]] %*% weights)
  
  if (dcov <= 0) {
    return(0)
  }
  
  sqrt(dcov / init[["denom"]])
}
## target
smd.target <- function(init, weights = NULL) {
  .check_init(init, "init_smd")
  
  col_w_smd(init$covs, treat = init$treat,
            weights = c(rep_with(1, weights), weights),
            s.weights = init$s.weights,
            bin.vars = init$bin.vars,
            s.d.denom = init$s.d.denom,
            abs = TRUE)
}
ks.target <- function(init, weights = NULL) {
  .check_init(init, "init_ks")
  col_w_ks(init$covs, treat = init$treat,
           weights = c(rep_with(1, weights), weights),
           s.weights = init$s.weights,
           bin.vars = init$bin.vars)
}
ovl.target <- function(init, weights = NULL) {
  .check_init(init, "init_ovl")
  col_w_ovl(init$covs, treat = init$treat,
            weights = c(rep_with(1, weights), weights),
            s.weights = init$s.weights,
            bin.vars = init$bin.vars, integrate = init$integrate)
}
mahalanobis.target <- function(init, weights = NULL) {
  .check_init(init, "init_mahalanobis")
  
  mean.diffs <- col_w_smd(init$covs, treat = init$treat,
                          weights = c(rep_with(1, weights), weights),
                          s.weights = init$s.weights,
                          bin.vars = init$bin.vars,
                          std = FALSE, abs = TRUE)
  
  drop(sqrt(t(mean.diffs) %*% init$sigma_inv %*% mean.diffs))
}
energy.dist.target <- function(init, weights = NULL) {
  .check_init(init, "init_energy.dist")
  
  weights <- {
    if (is_null(weights)) init[["s.weights"]]
    else weights * init[["s.weights"]]
  }
  
  weights <- weights / mean_fast(weights)
  
  drop(weights %*% init[["P"]] %*% weights) + sum(init[["q"]] * weights) + init[["k"]]
}

bal_criterion <- function(treat.type, criterion) {
  .chk_not_missing(treat.type, "`treat.type`")
  .chk_string(treat.type)
  .chk_not_missing(criterion, "`criterion`")
  .chk_string(criterion)
  
  switch(
    treat.type,
    binary = {
      switch(criterion,
             smd.mean = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_smd(covs, treat, s.weights, estimand)
                 }
                 mean_fast(smd.binary(init, weights))
               },
               init = init_smd
             ),
             smd.max = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_smd(covs, treat, s.weights, estimand)
                 }
                 max(smd.binary(init, weights))
               },
               init = init_smd
             ),
             smd.rms = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_smd(covs, treat, s.weights, estimand)
                 }
                 rms(smd.binary(init, weights))
               },
               init = init_smd
             ),
             ks.mean = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ks(covs, treat, s.weights)
                 }
                 mean_fast(ks.binary(init, weights))
               },
               init = init_ks
             ),
             ks.max = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ks(covs, treat, s.weights)
                 }
                 max(ks.binary(init, weights))
               },
               init = init_ks
             ),
             ks.rms = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ks(covs, treat, s.weights)
                 }
                 rms(ks.binary(init, weights))
               },
               init = init_ks
             ),
             ovl.mean = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, integrate = TRUE,
                              init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ovl(covs, treat, s.weights, integrate = integrate)
                 }
                 mean_fast(ovl.binary(init, weights))
               },
               init = init_ovl
             ),
             ovl.max = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, integrate = TRUE,
                              init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ovl(covs, treat, s.weights, integrate = integrate)
                 }
                 max(ovl.binary(init, weights))
               },
               init = init_ovl
             ),
             ovl.rms = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, integrate = TRUE,
                              init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ovl(covs, treat, s.weights, integrate = integrate)
                 }
                 rms(ovl.binary(init, weights))
               },
               init = init_ovl
             ),
             mahalanobis = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_mahalanobis(covs, treat, s.weights, estimand)
                 }
                 mahalanobis.binary(init, weights)
               },
               init = init_mahalanobis
             ),
             energy.dist = list(
               fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, improved = TRUE, init = NULL, ...) {
                 
                 if (is_null(init)) {
                   init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                 }
                 energy.dist.binary(init, weights)
               },
               init = init_energy.dist
             ),
             kernel.dist = list(
               fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, init = NULL, ...) {
                 
                 if (is_null(init)) {
                   init <- init_kernel.dist(covs, treat, s.weights, estimand, focal)
                 }
                 kernel.dist.binary(init, weights)
               },
               init = init_kernel.dist
             ),
             l1.med = list(
               fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_l1.med(covs, treat, s.weights, estimand, focal, ...)
                 }
                 l1.med.binary(init, weights)
               },
               init = init_l1.med
             ),
             r2 = list(
               fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_r2(covs, treat, s.weights)
                 }
                 r2.binary(init, weights)
               },
               init = init_r2
             ),
             r2.2 = list(
               fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_r2(covs, treat, s.weights, poly = 2)
                 }
                 r2.binary(init, weights)
               },
               init = function(...) init_r2(..., poly = 2)
             ),
             r2.3 = list(
               fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_r2(covs, treat, s.weights, poly = 3)
                 }
                 r2.binary(init, weights)
               },
               init = function(...) init_r2(..., poly = 3)
             )
      )},
    multinomial = {
      switch(criterion,
             smd.mean = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_smd(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                 }
                 mean_fast(smd.multinomial(init, weights))
               },
               init = init_smd
             ),
             smd.max = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_smd(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                 }
                 max(smd.multinomial(init, weights))
               },
               init = init_smd
             ),
             smd.rms = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_smd(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                 }
                 rms(smd.multinomial(init, weights))
               },
               init = init_smd
             ),
             ks.mean = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                              pairwise = TRUE, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ks(covs, treat, s.weights, focal = focal, pairwise = pairwise)
                 }
                 mean_fast(ks.multinomial(init, weights))
               },
               init = init_ks
             ),
             ks.max = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                              pairwise = TRUE, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ks(covs, treat, s.weights, focal = focal, pairwise = pairwise)
                 }
                 max(ks.multinomial(init, weights))
               },
               init = init_ks
             ),
             ks.rms = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                              pairwise = TRUE, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ks(covs, treat, s.weights, focal = focal, pairwise = pairwise)
                 }
                 rms(ks.multinomial(init, weights))
               },
               init = init_ks
             ),
             ovl.mean = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                              pairwise = TRUE, integrate = TRUE, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ovl(covs, treat, s.weights, focal = focal, pairwise = pairwise,
                                    integrate = integrate)
                 }
                 mean_fast(ovl.multinomial(init, weights))
               },
               init = init_ovl
             ),
             ovl.max = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                              pairwise = TRUE, integrate = TRUE, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ovl(covs, treat, s.weights, focal = focal, pairwise = pairwise,
                                    integrate = integrate)
                 }
                 max(ovl.multinomial(init, weights))
               },
               init = init_ovl
             ),
             ovl.rms = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                              pairwise = TRUE, integrate = TRUE, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ovl(covs, treat, s.weights, focal = focal, pairwise = pairwise,
                                    integrate = integrate)
                 }
                 rms(ovl.multinomial(init, weights))
               },
               init = init_ovl
             ),
             energy.dist = list(
               fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, improved = TRUE, init = NULL, ...) {
                 
                 if (is_null(init)) {
                   init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                 }
                 energy.dist.multinomial(init, weights)
               },
               init = init_energy.dist
             ),
             l1.med = list(
               fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_l1.med(covs, treat, s.weights, estimand, focal, ...)
                 }
                 l1.med.multinomial(init, weights)
               },
               init = init_l1.med
             )
      )},
    continuous = {
      switch(criterion,
             p.mean = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_p(covs, treat, s.weights)
                 }
                 mean_fast(pearson.corr.continuous(init, weights))
               },
               init = init_p
             ),
             p.max = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_p(covs, treat, s.weights)
                 }
                 max(pearson.corr.continuous(init, weights))
               },
               init = init_p
             ),
             p.rms = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_p(covs, treat, s.weights)
                 }
                 rms(pearson.corr.continuous(init, weights))
               },
               init = init_p
             ),
             s.mean = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_s(covs, treat, s.weights)
                 }
                 mean_fast(spearman.corr.continuous(init, weights))
               },
               init = init_s
             ),
             s.max = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_s(covs, treat, s.weights)
                 }
                 max(spearman.corr.continuous(init, weights))
               },
               init = init_s
             ),
             s.rms = list(
               fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_s(covs, treat, s.weights)
                 }
                 rms(spearman.corr.continuous(init, weights))
               },
               init = init_s
             ),
             r2 = list(
               fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_r2(covs, treat, s.weights)
                 }
                 r2.continuous(init, weights)
               },
               init = init_r2
             ),
             r2.2 = list(
               fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_r2(covs, treat, s.weights, poly = 2)
                 }
                 r2.continuous(init, weights)
               },
               init = function(...) init_r2(..., poly = 2)
             ),
             r2.3 = list(
               fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_r2(covs, treat, s.weights, poly = 3)
                 }
                 r2.continuous(init, weights)
               },
               init = function(...) init_r2(..., poly = 3)
             ),
             distance.cov = list(
               fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_distance.cov(covs, treat, s.weights)
                 }
                 distance.cov.continuous(init, weights)
               },
               init = init_distance.cov
             ),
             distance.cor = list(
               fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_distance.cov(covs, treat, s.weights, std = TRUE)
                 }
                 distance.cov.continuous(init, weights)
               },
               init = function(...) init_distance.cov(..., std = TRUE)
             )
      )},
    target = {
      switch(criterion,
             smd.mean = list(
               fun = function(covs, weights = NULL, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_smd(covs, s.weights)
                 }
                 mean_fast(smd.target(init, weights))
               },
               init = init_smd
             ),
             smd.max = list(
               fun = function(covs, weights = NULL, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_smd(covs, s.weights)
                 }
                 max(smd.target(init, weights))
               },
               init = init_smd
             ),
             smd.rms = list(
               fun = function(covs, weights = NULL, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_smd(covs, s.weights)
                 }
                 rms(smd.target(init, weights))
               },
               init = init_smd
             ),
             ks.mean = list(
               fun = function(covs, weights = NULL, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ks(covs, s.weights)
                 }
                 mean_fast(ks.target(init, weights))
               },
               init = init_ks
             ),
             ks.max = list(
               fun = function(covs, weights = NULL, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ks(covs, s.weights)
                 }
                 max(ks.target(init, weights))
               },
               init = init_ks
             ),
             ks.rms = list(
               fun = function(covs, weights = NULL, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ks(covs, s.weights)
                 }
                 rms(ks.target(init, weights))
               },
               init = init_ks
             ),
             ovl.mean = list(
               fun = function(covs, weights = NULL, s.weights = NULL, integrate = TRUE,
                              init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ovl(covs, s.weights, integrate = integrate)
                 }
                 mean_fast(ovl.target(init, weights))
               },
               init = init_ovl
             ),
             ovl.max = list(
               fun = function(covs, weights = NULL, s.weights = NULL, integrate = TRUE,
                              init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ovl(covs, s.weights, integrate = integrate)
                 }
                 max(ovl.target(init, weights))
               },
               init = init_ovl
             ),
             ovl.rms = list(
               fun = function(covs, weights = NULL, s.weights = NULL, integrate = TRUE,
                              init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_ovl(covs, s.weights, integrate = integrate)
                 }
                 rms(ovl.target(init, weights))
               },
               init = init_ovl
             ),
             mahalanobis = list(
               fun = function(covs, weights = NULL, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_mahalanobis(covs, s.weights)
                 }
                 mahalanobis.target(init, weights)
               },
               init = init_mahalanobis
             ),
             energy.dist = list(
               fun = function(covs, weights = NULL, s.weights = NULL, init = NULL, ...) {
                 if (is_null(init)) {
                   init <- init_energy.dist(covs, s.weights)
                 }
                 energy.dist.target(init, weights)
               },
               init = init_energy.dist
             )
      )}
  )
}

.check_init <- function(init, init_class) {
  .chk_not_missing(init, "`init`")
  .chk_not_missing(init_class, "`init_class`")
  .chk_is(init, init_class)
}
