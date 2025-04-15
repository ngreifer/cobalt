#' @title Compute Balance and Summary Statistics for Covariates
#' 
#' @description These functions quickly compute balance statistics for the given covariates. These functions are used in [bal.tab()], but they are available for use in programming without having to call `bal.tab()` to get them.
#' \itemize{
#'     \item{`col_w_mean()` computes the (weighted) means for a set of covariates and weights and is essentially a weighted version of [colMeans()].}
#'     \item{`col_w_sd()` computes the (weighted) standard deviations for a set of covariates and weights.}
#'     \item{`col_w_smd()` computes the (weighted) (absolute) (standardized) difference in means for a set of covariates, a binary treatment, and weights.}
#'     \item{`col_w_vr()` computes the (weighted) variance ratio for a set of covariates, a binary treatment, and weights.}
#'     \item{`col_w_ks()` computes the (weighted) Kolmogorov-Smirnov (KS) statistic for a set of covariates, a binary treatment, and weights.}
#'     \item{`col_w_ovl()` computes the complement of the (weighted) overlapping coefficient compliment for a set of covariates, a binary treatment, and weights (based on Franklin et al, 2014).}
#'     \item{`col_w_cov()` and `col_w_corr()` compute the (weighted) (absolute) treatment-covariate covariance or correlation for a set of covariates, a continuous treatment, and weights.}
#'     \item{`col_w_dcov()` and `col_w_dcorr()` compute the (weighted) treatment-covariate distance covariance or distance correlation for a set of covariates, a continuous treatment, and weights.}
#' }
#' 
#' @param mat a numeric matrix or a data frame containing the covariates for which the statistic is to be computed. If a data frame, [splitfactor()] with `drop.first = "if2"` will be called if any character or factor variables are present. This can slow down the function, so it's generally best to supply a numeric matrix. If a numeric vector is supplied, it will be converted to a 1-column matrix first.
#' @param weights `numeric`; an optional set of weights used to compute the weighted statistics. If sampling weights are supplied through `s.weights`, the `weights` should not incorporate these weights, as `weights` and `s.weights` will be multiplied together prior to computing the weighted statistics.
#' @param s.weights `numeric`; an optional set of sampling weights used to compute the weighted statistics. If weights are supplied through `weights`, `weights` and `s.weights` will be multiplied together prior to computing the weighted statistics. Some functions use `s.weights` in a particular way; for others, supplying `weights` and `s.weights` is equivalent to supplying their product to either `weights` or `s.weights`. See Details.
#' @param subset a `logical` vector with length equal to the number of rows of `mat` used to subset the data. See Details for notes on its use with `col_w_smd()`, `col_w_cov()`, and `col_w_corr()`.
#' @param na.rm `logical`; whether `NA`s should be ignored or not. If `FALSE`, any variable with any `NA`s will have its corresponding statistic returned as `NA`. If `TRUE`, any variable with any `NA`s will have its corresponding statistic computed as if the missing value were not there.
#' @param treat a vector of treatment status for each individual. For `col_w_smd()`, `col_w_vr()`, `col_w_ks()`, and `col_w_ovl()`, `treat` should have exactly two unique values. For `col_w_cov()`, `col_w_corr()`, `col_w_dcov()`, and `col_w_dcorr()`, `treat` should be a many-valued numeric vector.
#' @param std `logical`; for `col_w_smd()`, whether the computed mean differences for each variable should be standardized; for `col_w_cov()`, whether treatment-covariate correlations should be computed (`TRUE`) rather than covariances (`FALSE`); for `col_w_dcov()`, whether treatment-covariate distance correlations should be computed (`TRUE`) rather than distance covariances (`FALSE`). Can have either length 1, whereby all variables will be standardized or not, or length equal to the number of columns of `mat`, whereby only variables with a value of `TRUE` will be standardized. See Details.
#' @param s.d.denom for `col_w_smd()`, `col_w_cov()`, and `col_w_dcov()` when `std` is `TRUE` for some variables, and for `col_w_corr()` and `col_w_dcorr()`, how the standardization factor should be computed. For `col_w_smd()` (i.e., when computing standardized mean differences), allowable options include 
#' \itemize{
#'     \item{`"treated"` - uses the standard deviation of the variable in the treated group}
#'     \item{`"control"` - uses the standard deviation of the variable in the control group}
#'     \item{`"pooled"` - uses the square root of the average of the variances of the variable in the treated and control groups}
#'     \item{`"all"` - uses the standard deviation of the variable in the full sample}
#'     \item{`"weighted"` - uses the standard deviation of the variable in the full sample weighted by `weighted.weights`}
#'     \item{`"hedges"` - uses the small-sample corrected version of Hedge's G described in the WWC Procedures Handbook (see References)}
#' \item{the name of one of the treatment values - uses the standard deviation of the variable in that treatment group.}
#' }
#' For `col_w_cov()`, `col_w_corr()`, `col_w_dcov()`, and `col_w_dcorr()`, only `"all"` and `"weighted"` are allowed. Abbreviations allowed. This can also be supplied as a numeric vector of standard deviations with length equal to the number of columns of `mat`; the values will be used as the standardization factors.
#' @param abs `logical`; for `col_w_smd()`, `col_w_cov()`, and `col_w_corr()`, whether the returned statistics should be in absolute value (`TRUE`) or not. For `col_w_vr()`, whether the ratio should always include the larger variance in the numerator, so that the ratio is always greater than or equal to 1. Default is `FALSE`.
#' @param bin.vars a vector used to denote whether each variable is binary or not. Can be a `logical` vector with length equal to the number of columns of `mat` or a vector of numeric indices or character names of the binary variables. If missing (the default), the function will figure out which covariates are binary or not, which can increase computation time. If `NULL`, it will be assumed no variables are binary. All functions other than `col_w_mean()` treat binary variables different from continuous variables. If a factor or character variable is in `mat`, all the dummies created will automatically be marked as binary, but it should still receive an entry when `bin.vars` is supplied as `logical`.
#' @param weighted.weights for `col_w_smd()`, `col_w_cov()`, `col_w_corr()`, `col_w_dcov()`, and `col_w_dcorr()`, when `std = TRUE` and `s.d.denom = "weighted"`, a vector of weights to be applied to the computation of the denominator standard deviation. If not specified, will use the argument to `weights`. When `s.d.denom` is not `"weighted"`, this is ignored. The main purpose of this is to allow `weights` to be `NULL` while weighting the denominator standard deviations for assessing balance in the unweighted sample but using the standard deviations of the weighted sample.
#' @param type for `col_w_cov()` and `col_w_corr()`, the type of covariance/correlation to be computed. Allowable options include `"pearson"` and `"spearman"`. When `"spearman"` is requested, the covariates and treatment are first turned into ranks using [rank()] with `na.last = "keep"`.
#' @param integrate `logical`; for `col_w_ovl()`, whether to use [integrate()] to calculate the area of overlap for continuous variables. If `FALSE`, a midpoint Riemann sum will be used instead. The Riemann sum is a little slower and very slightly imprecise (unnoticibly in most contexts). When `TRUE`, `integrate()` will be tried, and if it fails, the Riemann sum will be used as a fallback. The default (`TRUE`) is to use `integrate()` when possible.
#' @param steps for `col_w_ovl()` when `integrate = FALSE`, the number of points to use to compute the Riemann sum to approximate the integral. Default is 1001 for 1000 partitions.
#' @param ... for all functions, additional arguments supplied to [splitfactor()] when `mat` is a data.frame. `data`, `var.name`, `drop.first`, and `drop.level` are ignored; `drop.first` is automatically set to `"if2"`. For `col_w_ovl()`, other arguments passed to [density()] besides `x` and `weights`. Note that the default value for `bw` when unspecified is `"nrd"` rather than the default in `density()`, which is `"nrd0"`.
#' @returns A vector of balance statistics, one for each variable in `mat`. If `mat` has column names, the output will be named as well.
#' 
#' @details
#' `col_w_mean()` computes column weighted means for a matrix of variables. It is similar to [colMeans()] but (optionally) incorporates weights. `weights` and `s.weights` are multiplied together prior to being used, and there is no distinction between them. This could be used to compute the weighted means of each covariate in the general population to examine the degree to which a weighting method has left the weighted samples resembling the original population.
#' 
#' `col_w_sd()` computes column weighted standard deviations for a matrix of variables. `weights` and `s.weights` are multiplied together prior to being used, and there is no distinction between them. The variance of binary variables is computed as \eqn{p(1-p)}, where \eqn{p} is the (weighted) proportion of 1s, while the variance of continuous variables is computed using the standard formula; the standard deviation is the square root of this variance.
#' 
#' `col_w_smd()` computes the mean difference for each covariate between treatment groups defined by `treat`. These mean differences can optionally be weighted, standardized, and/or in absolute value. The standardization factor is computed using the unweighted standard deviation or variance when `s.weights` are absent, and is computed using the `s.weights`-weighted standard deviation or variance when `s.weights` are present, except when `s.d.denom = "weighted"`, in which case the product of `weighted.weights` and `s.weights` (if present) are used to weight the standardization factor. The standardization factor is computed using the whole sample even when `subset` is used. Note that unlike `bal.tab()`, `col_w_smd()` requires the user to specify whether each individual variable should be standardized using `std` rather than relying on `continuous` or `binary`. The weighted mean difference is computed using the product of `weights` and `s.weights`, if specified. The variance of binary variables is computed as \eqn{p(1-p)}, where \eqn{p} is the (weighted) proportion of 1s, while the variance of continuous variables is computed using the standard formula. 
#' 
#' `col_w_vr()` computes the variance ratio for each covariate between treatment groups defined by `treat`. When `abs = TRUE`, `pmax(out, 1/out)` is applied to the output so that the ratio is always greater than or equal to 1. For binary variables, the variance is computed as \eqn{p(1-p)}, where \eqn{p} is the (weighted) proportion of 1s, while the variance of continuous variables is computed using the standard formula. Note that in `bal.tab()`, variance ratios are not computed for binary variables, while here, they are (but likely should not be interpreted). `weights` and `s.weights` are multiplied together prior to being used, and there is no distinction between them. Because of how the weighted variance is computed, exactly balanced groups may have variance ratios that differ slightly from 1.
#' 
#' `col_w_ks()` computes the KS statistic for each covariate using the method implemented in \pkg{twang}. The KS statistics can optionally be weighted. For binary variables, the KS statistic is just the difference in proportions. `weights` and `s.weights` are multiplied together prior to being used, and there is no distinction between them.
#' 
#' `col_w_ovl()` computes the complement of the overlapping coefficient as described by Franklin et al. (2014). It does so by computing the density of the covariate in the treated and control groups, then finding the area where those density overlap, and subtracting that number from 1, yielding a value between 0 and 1 where 1 indicates complete imbalance, and 0 indicates perfect balance. [density()] is used to model the density in each group. The bandwidth of the covariate in the smaller treatment group is used for both groups. The area of overlap can be computed using `integrate`, which quickly and accurately computes the integral, or using a midpoint Riemann sum with 1000 partitions, which approximates the area more slowly. A reason to prefer the Riemann sum is that `integrate` can fail for unknown reasons, though Riemann sums will fail with some extreme distributions. When either method fails, the resulting value will be `NA`. For binary variables, the complement of the overlapping coefficient is just the difference in proportions. `weights` and `s.weights` are multiplied together prior to being used, and there is no distinction between them. The weights are used to compute the weighted density by supplying them to the `weights` argument of `density`.
#' 
#' `col_w_cov()` computes the covariances between a continuous treatment and the covariates to assess balance for a continuous treatment as recommended in Austin (2019). These covariances can optionally be weighted or in absolute value or can be requested as correlations (i.e., standardized covariances). Each correlations is computed as the covariance between the treatment and covariate divided by a standardization factor, which is equal to the square root of the product of the variance of treatment and the variance of the covariate. The standardization factor is computed using the unweighted variances when `s.weights` are absent, and is computed using the sampling weighted variances when `s.weights` are present, except when `s.d.denom = "weighted"`, in which case the product of `weighted.weights` and `s.weights` (if present) are used to weight the standardization factor. For this reason, the computed correlation can be greater than 1 or less than -1. The standardization factor is always computed using the whole sample even when `subset` is used. The covariance is computed using the product of `weights` and `s.weights`, if specified. The variance of binary variables is computed as \eqn{p(1-p)}, where \eqn{p} is the (weighted) proportion of 1s, while the variance of continuous variables is computed using the standard formula. 
#' 
#' `col_w_corr()` is a wrapper for `col_w_cov` with `std` set to `TRUE`. 
#' 
#' `col_w_dcov()` computes the distance covariances between a continuous treatment and the covariates to assess balance for a continuous treatment. A multivariate version is described by Huling et al. (2023) for computing a scalar value that represents the balance for all covariates simultaneously; the statistic computed here is for one covariate at a time. The distance covariances can optionally be weighted or can be requested as distance correlations (i.e., standardized distance covariances). The distance correlations are computed as the distance covariance between the treatment and covariate divided by a standardization factor, which is equal to the square root of the product of the distance variance of treatment and the distance variance of the covariate, where the distance variance is the distance covariance of a variable with itself. The standardization factor is computed using the unweighted distance variances when `s.weights` are absent, and is computed using the sampling weighted distance variances when `s.weights` are present, except when `s.d.denom = "weighted"`, in which case the product of `weighted.weights` and `s.weights` (if present) are used to weight the standardization factor. For this reason, the computed distance correlation can be greater than 1. The standardization factor is always computed using the whole sample even when `subset` is used. The distance covariance is computed using the product of `weights` and `s.weights`, if specified.
#' 
#' `col_w_dcorr()` is a wrapper for `col_w_dcov` with `std` set to `TRUE`. 
#' 
#' @references 
#' Austin, P. C. (2019). Assessing covariate balance when using the generalized propensity score with quantitative or continuous exposures. *Statistical Methods in Medical Research*, 28(5), 1365–1377. \doi{10.1177/0962280218756159}
#' 
#' Franklin, J. M., Rassen, J. A., Ackermann, D., Bartels, D. B., & Schneeweiss, S. (2014). Metrics for covariate balance in cohort studies of causal effects. *Statistics in Medicine*, 33(10), 1685–1699. \doi{10.1002/sim.6058}
#' 
#' Huling, J. D., Greifer, N., & Chen, G. (2023). Independence Weights for Causal Inference with Continuous Treatments. *Journal of the American Statistical Association*, 0(0), 1–14. \doi{10.1080/01621459.2023.2213485}
#' 
#' What Works Clearinghouse. (2020). WWC Procedures Handbook (Version 4.1). Retrieved from
#' <https://ies.ed.gov/ncee/wwc/Handbooks>
#' 
#' @seealso 
#' * [bal.tab()]
#' * [bal.compute()]
#' * [balance-statistics]
#' @examplesIf requireNamespace("WeightIt", quietly = TRUE)
#' data("lalonde", package = "cobalt")
#' 
#' treat <- lalonde$treat
#' covs <- subset(lalonde, select = -c(treat, re78))
#' covs0 <- splitfactor(covs, drop.first = "if2")
#' bin.vars <- c(FALSE, FALSE, TRUE, TRUE, TRUE,
#'               TRUE, TRUE, FALSE, FALSE)
#' W <- WeightIt::weightit(treat ~ covs, method = "glm", 
#'                         estimand = "ATE")
#' weights <- W$weights
#' 
#' round(data.frame(
#'     m0 = col_w_mean(covs0, weights = weights, subset = treat == 0),
#'     sd0 = col_w_sd(covs0, weights = weights,
#'                    bin.vars = bin.vars, subset = treat == 0),
#'     m1 = col_w_mean(covs0, weights = weights, subset = treat == 1),
#'     sd1 = col_w_sd(covs0, weights = weights,
#'                    bin.vars = bin.vars, subset = treat == 1),
#'     smd = col_w_smd(covs0, treat = treat, weights = weights,
#'                     std = TRUE, bin.vars = bin.vars),
#'     vr = col_w_vr(covs0, treat = treat, weights = weights,
#'                   bin.vars = bin.vars),
#'     ks = col_w_ks(covs0, treat = treat, weights = weights,
#'                   bin.vars = bin.vars),
#'     ovl = col_w_ovl(covs0, treat = treat, weights = weights,
#'                   bin.vars = bin.vars),
#'     row.names = colnames(covs0)
#' ), 4)
#' 
#' # Compare to bal.tab():
#' bal.tab(covs, treat = treat, weights = weights,
#'         disp = c("m", "sd"),
#'         stats = c("m", "v", "ks", "ovl"),
#'         estimand = "ATE", method = "weighting",
#'         binary = "std")
#' 
#' 

#' @name balance-summary
#' @export 
col_w_mean <- function(mat, weights = NULL, s.weights = NULL, subset = NULL, na.rm = TRUE, ...) {
  
  mat <- process_mat1(mat, ...)
  
  check_arg_lengths(mat, weights, s.weights, subset)
  
  if (is_null(weights)) weights <- rep.int(1, NROW(mat))
  if (is_null(s.weights)) s.weights <- rep.int(1, NROW(mat))
  
  .chk_null_or(subset, .chk_logical)
  if (is_null(subset)) subset <- rep.int(TRUE, NROW(mat))
  
  weights <- weights * s.weights
  
  if (all(weights == 0)) {
    .err("at least 1 unit must have a nonzero weight to compute weighted means")
  }
  
  col.w.m(mat[subset, , drop = FALSE], w = weights[subset], na.rm = na.rm)
}

#' @rdname balance-summary
#' @export 
col_w_sd <- function(mat, weights = NULL, s.weights = NULL, bin.vars, subset = NULL, na.rm = TRUE, ...) {
  
  mat <- process_mat2(mat, .bin.vars = bin.vars, ...)
  bin.vars <- attr(mat, "bin")
  
  check_arg_lengths(mat, weights, s.weights, subset)
  
  if (is_null(weights)) weights <- rep.int(1, NROW(mat))
  if (is_null(s.weights)) s.weights <- rep.int(1, NROW(mat))
  
  .chk_null_or(subset, .chk_logical)
  if (is_null(subset)) subset <- rep.int(TRUE, NROW(mat))
  
  weights <- weights * s.weights
  
  if (sum(weights != 0) < 2L) {
    .err("at least 2 units must have nonzero weights to compute weighted standard deviations")
  }
  
  sqrt(col.w.v(mat[subset, , drop = FALSE], w = weights[subset], 
               bin.vars = bin.vars, na.rm = na.rm))
}

#' @rdname balance-summary
#' @export 
col_w_smd <- function(mat, treat, weights = NULL, std = TRUE, s.d.denom = "pooled", abs = FALSE,
                      s.weights = NULL, bin.vars, subset = NULL, weighted.weights = weights, na.rm = TRUE, ...) {
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  .chk_not_any_na(treat)
  
  mat <- process_mat2(mat, ..., .bin.vars = bin.vars)
  bin.vars <- attr(mat, "bin")
  
  .chk_logical(std)
  .chk_not_any_na(std)
  if (length(std) %nin% c(1L, NCOL(mat))) {
    .err("`std` must have length equal to 1 or the number of columns of `mat`")
  }
  
  .chk_flag(abs)
  
  check_arg_lengths(mat, treat, weights, s.weights, subset)
  
  if (is_null(weights)) weights <- rep.int(1, NROW(mat))
  if (is_null(s.weights)) s.weights <- rep.int(1, NROW(mat))
  
  .chk_null_or(subset, .chk_logical)
  if (is_null(subset)) subset <- rep.int(TRUE, NROW(mat))
  
  if (!is_binary(treat[subset])) {
    .err("`treat` must be a binary variable")
  }
  
  weights <- weights * s.weights
  
  if (length(std) == 1L) {
    std <- rep.int(std, NCOL(mat))
  }
  
  tval1_0 <- treat[1L]
  
  if (all(weights[treat == tval1_0] == 0) || 
      all(weights[treat != tval1_0] == 0)) {
    .err("at least 1 unit in each level of `treat` must have a nonzero weight to compute weighted SMDs")
  }
  
  m1 <- col.w.m(mat[treat == tval1_0 & subset, , drop = FALSE], weights[treat == tval1_0 & subset], na.rm = na.rm)
  m0 <- col.w.m(mat[treat != tval1_0 & subset, , drop = FALSE], weights[treat != tval1_0 & subset], na.rm = na.rm)
  diffs <- m1 - m0
  zeros <- check_if_zero(diffs)
  
  to.sd <- std & !is.na(zeros) & !zeros
  if (any(to.sd)) {
    if (!is.numeric(s.d.denom) || length(s.d.denom) != ncol(mat)) {
      s.d.denom <- .get_s.d.denom(s.d.denom, weights = list(weights), treat = treat,
                                  quietly = TRUE)
      
      denoms <- .compute_s.d.denom(mat, treat = treat, 
                                   s.d.denom = s.d.denom, s.weights = s.weights, 
                                   bin.vars = bin.vars, subset = subset, to.sd = to.sd,
                                   weighted.weights = weighted.weights, na.rm = na.rm)
    }
    else {
      denoms <- s.d.denom
    }
    
    diffs[to.sd] <- diffs[to.sd] / denoms[to.sd]
  }
  
  if (abs) {
    diffs <- abs(diffs)
  }
  else {
    tval1 <- {
      if (is_0_1(treat))  1
      else get_treated_level(treat[subset])
    }
    
    if (tval1 != tval1_0) {
      diffs <- -diffs
    }
  }
  
  setNames(diffs, colnames(mat))
}

#' @rdname balance-summary
#' @export 
col_w_vr <- function(mat, treat, weights = NULL, abs = FALSE, s.weights = NULL, bin.vars,
                     subset = NULL, na.rm = TRUE, ...) {
  
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  .chk_not_any_na(treat)
  
  mat <- process_mat2(mat, ..., .bin.vars = bin.vars)
  bin.vars <- attr(mat, "bin")
  
  .chk_flag(abs)
  
  check_arg_lengths(mat, treat, weights, s.weights, subset)
  
  if (is_null(weights)) weights <- rep.int(1, NROW(mat))
  if (is_null(s.weights)) s.weights <- rep.int(1, NROW(mat))
  
  .chk_null_or(subset, .chk_logical)
  if (is_null(subset)) subset <- rep.int(TRUE, NROW(mat))
  
  if (!is_binary(treat[subset])) {
    .err("`treat` must be a binary variable")
  }
  
  weights <- weights * s.weights
  
  weights <- weights[subset]
  treat <- treat[subset]
  mat <- mat[subset, , drop = FALSE]
  
  tval1_0 <- treat[1L]
  
  if (sum(weights[treat == tval1_0] != 0) < 2L || 
      sum(weights[treat != tval1_0] != 0) < 2L) {
    .err("at least 2 units in each level of `treat` must have nonzero weights to compute weighted variance ratios")
  }
  
  v1 <- col.w.v(mat[treat == tval1_0, , drop = FALSE], weights[treat == tval1_0],
                bin.vars = bin.vars, na.rm = na.rm)
  v0 <- col.w.v(mat[treat != tval1_0, , drop = FALSE], weights[treat != tval1_0],
                bin.vars = bin.vars, na.rm = na.rm)
  
  v.ratios <- v1 / v0
  
  if (abs) {
    v.ratios <- abs_(v.ratios, ratio = TRUE)
  }
  else {
    tval1 <- {
      if (is_0_1(treat))  1
      else get_treated_level(treat[subset])
    }
    
    if (tval1 != tval1_0) {
      v.ratios <- 1 / v.ratios
    }
  }
  
  setNames(v.ratios, colnames(mat))
}

#' @rdname balance-summary
#' @export 
col_w_ks <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars, subset = NULL,
                     na.rm = TRUE, ...) {
  
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  .chk_not_any_na(treat)
  
  mat <- process_mat2(mat, ..., .bin.vars = bin.vars)
  bin.vars <- attr(mat, "bin")
  
  check_arg_lengths(mat, treat, weights, s.weights, subset)
  
  if (is_null(weights)) weights <- rep.int(1, NROW(mat))
  if (is_null(s.weights)) s.weights <- rep.int(1, NROW(mat))
  
  .chk_null_or(subset, .chk_logical)
  if (is_null(subset)) subset <- rep.int(TRUE, NROW(mat))
  
  if (!is_binary(treat[subset])) {
    .err("`treat` must be a binary variable")
  }
  
  weights <- weights * s.weights
  
  weights <- weights[subset]
  treat <- treat[subset]
  mat <- mat[subset, , drop = FALSE]
  
  tval1 <- treat[1L]
  ks <- rep.int(NA_real_, NCOL(mat))
  
  if (all(weights[treat == tval1] == 0) || 
      all(weights[treat != tval1] == 0)) {
    .err("at least 1 unit in each level of `treat` must have a nonzero weight to compute weighted KS statistics")
  }
  
  if (!all(bin.vars)) {
    weights_ <- weights
    weights_[treat == tval1] <-  weights[treat == tval1] / sum(weights[treat == tval1])
    weights_[treat != tval1] <- -weights[treat != tval1] / sum(weights[treat != tval1])
    
    ks[!bin.vars] <- apply(mat[, !bin.vars, drop = FALSE], 2L, function(x) {
      if (anyNA(x)) {
        if (na.rm) x <- na.rem(x)
        else return(NA_real_)
      }
      
      ordered.index <- order(x)
      cumv <- abs(cumsum(weights_[ordered.index]))[c(diff(x[ordered.index]) != 0, TRUE)]
      
      if (is_null(cumv)) 0 else max(cumv)
    })
  }
  
  if (any(bin.vars)) {
    ks[bin.vars] <- abs(col.w.m(mat[treat == tval1, bin.vars, drop = FALSE],
                                weights[treat == tval1], na.rm = na.rm) - 
                          col.w.m(mat[treat != tval1, bin.vars, drop = FALSE],
                                  weights[treat != tval1], na.rm = na.rm))
  }
  
  setNames(ks, colnames(mat))
}

#' @rdname balance-summary
#' @export 
col_w_ovl <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars, 
                      subset = NULL, na.rm = TRUE, integrate = TRUE, steps = 1001L, ...) {
  
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  .chk_not_any_na(treat)
  
  mat <- process_mat2(mat, ..., .bin.vars = bin.vars)
  bin.vars <- attr(mat, "bin")
  
  check_arg_lengths(mat, treat, weights, s.weights, subset)
  
  if (is_null(weights)) weights <- rep.int(1, NROW(mat))
  if (is_null(s.weights)) s.weights <- rep.int(1, NROW(mat))
  
  .chk_null_or(subset, .chk_logical)
  if (is_null(subset)) subset <- rep.int(TRUE, NROW(mat))
  
  if (!is_binary(treat[subset])) {
    .err("`treat` must be a binary variable")
  }
  
  weights <- weights * s.weights
  
  weights <- weights[subset]
  treat <- treat[subset]
  mat <- mat[subset, , drop = FALSE]
  
  tval1 <- treat[1L]
  
  .chk_gte(weights, 0)
  
  if (all(weights[treat == tval1] == 0) || 
      all(weights[treat != tval1] == 0)) {
    .err("at least 1 unit in each level of `treat` must have a nonzero weight to compute weighted OVL statistics")
  }
  
  if (check_if_zero(sum(weights[treat == tval1])) || 
      check_if_zero(sum(weights[treat != tval1]))) {
    .err("the sum of weights in each treatment group must be nonzero to compute weighted OVL statistics")
  }
  
  all_pos_w <- all(weights >= 0)
  steps <- 1001L
  
  unique.treat <- unique(treat, nmax = 2L)
  t.sizes <- setNames(vapply(unique.treat, function(x) sum(treat == x), numeric(1L)),
                      unique.treat)
  smallest.t <- names(t.sizes)[which.min(t.sizes)]
  ovl <- setNames(numeric(ncol(mat)), colnames(mat))
  
  if (!all(bin.vars)) {
    .chk_flag(integrate)
    
    if (!integrate) {
      .chk_count(steps)
      .chk_gte(steps, 5)
    }
    
    bw <- ...get("bw", "nrd")
    
    A <- ...mget(setdiff(names(formals(density_neg_w_safe)),
                         c("x", "weights", "bw")))
    
    bw_fun <- get0(paste.("bw", bw))
    if (!is.function(bw_fun)) {
      .err(sprintf("%s is not an acceptable entry to `bw`. See `?stats::density` for allowable options",
                   add_quotes(bw)))
    }
    
    .w_ovl <- function(x) {
      if (anyNA(x)) {
        if (!na.rm) {
          return(NA_real_)
        }
        t <- treat[!is.na(x)]
        w <- weights[!is.na(x)]
        x <- x[!is.na(x)]
      }
      else {
        t <- treat
        w <- weights
      }
      
      #If perfect separation, return 1
      if (all(all_pos_w) &&
          (min(x[t == tval1]) > max(x[t != tval1]) ||
           min(x[t != tval1]) > max(x[t == tval1]))) {
        return(1)
      }
      
      x <- center(x) / sd(x)
      
      A[["bw"]] <- bw_fun(x[t == smallest.t])
      
      min.c <- min(x) - 4 * A[["bw"]]
      max.c <- max(x) + 4 * A[["bw"]]
      
      A[["x"]] <- x[t == tval1]
      A[["weights"]] <- w[t == tval1] / sum(w[t == tval1])
      
      f1_ <- approxfun(do.call(density_neg_w_safe, A))
      
      f1 <- function(z) {
        y <- f1_(z)
        if (anyNA(y)) {
          y[is.na(y)] <- 0
        }
        y
      }
      
      A[["x"]] <- x[t != tval1]
      A[["weights"]] <- w[t != tval1] / sum(w[t != tval1])
      
      f0_ <- approxfun(do.call(density_neg_w_safe, A))
      
      f0 <- function(z) {
        y <- f0_(z)
        if (anyNA(y)) {
          y[is.na(y)] <- 0
        }
        y
      }
      
      fn <- function(z) {
        pmin(f1(z), f0(z))
      }
      
      s <- NULL
      if (integrate) {
        s <- try(integrate(fn, lower = min.c, upper = max.c)$value,
                 silent = TRUE)
      }
      
      if (null_or_error(s)) {
        s <- intapprox(fn, min.c, max.c, steps = steps, method = "midpoint")
      }
      
      min(max(1 - s, 0), 1) #Reverse: measure imbalance
    }
    
    ovl[!bin.vars] <- apply(mat[, !bin.vars, drop = FALSE], 2L, .w_ovl)
  }
  
  if (any(bin.vars)) {
    ovl[bin.vars] <- abs(col.w.m(mat[treat == tval1, bin.vars, drop = FALSE],
                                 weights[treat == tval1], na.rm = na.rm) - 
                           col.w.m(mat[treat != tval1, bin.vars, drop = FALSE],
                                   weights[treat != tval1], na.rm = na.rm))
  }
  
  ovl
}

#' @rdname balance-summary
#' @export 
col_w_cov <- function(mat, treat, weights = NULL, type = "pearson", std = FALSE,
                      s.d.denom = "all", abs = FALSE, s.weights = NULL, bin.vars,
                      subset = NULL, weighted.weights = weights, na.rm = TRUE, ...) {
  
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  .chk_not_any_na(treat)
  
  mat <- process_mat2(mat, ..., .bin.vars = bin.vars)
  bin.vars <- attr(mat, "bin")
  
  .chk_logical(std)
  .chk_not_any_na(std)
  if (length(std) %nin% c(1L, NCOL(mat))) {
    .err("`std` must have length equal to 1 or the number of columns of `mat`")
  }
  
  .chk_flag(abs)
  
  check_arg_lengths(mat, treat, weights, s.weights, subset)
  
  if (is_null(weights)) weights <- rep.int(1, NROW(mat))
  if (is_null(s.weights)) s.weights <- rep.int(1, NROW(mat))
  
  .chk_null_or(subset, .chk_logical)
  if (is_null(subset)) subset <- rep.int(TRUE, NROW(mat))
  
  if (length(std) == 1L) {
    std <- rep.int(std, NCOL(mat))
  }
  
  .chk_string(type)
  type <- tolower(type)
  type <- match_arg(type, c("pearson", "spearman"))
  if (type == "spearman") {
    for (i in which(!bin.vars)) {
      mat[,i] <- rank(mat[,i], na.last = "keep")
    }
    treat <- rank(treat, na.last = "keep")
  }
  
  weights <- weights * s.weights
  
  if (sum(weights != 0) <= 1) {
    .err("at least 2 units must have nonzero weights to compute weighted covariances")
  }
  
  covars <- col.w.cov(mat[subset, , drop = FALSE], y = treat[subset],
                      w = weights[subset], na.rm = na.rm)
  
  zeros <- check_if_zero(covars)
  to.sd <- std & !is.na(zeros) & !zeros
  
  if (any(to.sd)) {
    if (is.numeric(s.d.denom) && length(s.d.denom) == ncol(mat)) {
      denoms <- s.d.denom
    }
    else {
      s.d.denom <- .get_s.d.denom.cont(s.d.denom, weights = list(weights),
                                       quietly = TRUE)
      
      denoms <- .compute_s.d.denom(mat, treat = treat, 
                                   s.d.denom = s.d.denom, s.weights = s.weights, 
                                   bin.vars = bin.vars, subset = subset, to.sd = to.sd,
                                   weighted.weights = weighted.weights, na.rm = na.rm)
    }
    
    covars[to.sd] <- covars[to.sd] / denoms[to.sd]
  }
  
  if (abs) covars <- abs(covars)
  
  setNames(covars, colnames(mat))
}

#' @rdname balance-summary
#' @export 
col_w_corr <- function(mat, treat, weights = NULL, type = "pearson", s.d.denom = "all",
                       abs = FALSE, s.weights = NULL, bin.vars, subset = NULL,
                       weighted.weights = weights, na.rm = TRUE, ...) {
  .call <- match.call(expand.dots = TRUE)
  .call[[1L]] <- quote(cobalt::col_w_cov)
  .call[["std"]] <- quote(TRUE)
  eval.parent(.call)
}

#' @rdname balance-summary
#' @export 
col_w_dcov <- function(mat, treat, weights = NULL, std = FALSE, s.d.denom = "all",
                       s.weights = NULL, subset = NULL, weighted.weights = weights,
                       na.rm = TRUE, ...) {
  .chk_not_missing(treat, "`treat`")
  .chk_atomic(treat)
  .chk_not_any_na(treat)
  
  .chk_logical(std)
  .chk_not_any_na(std)
  if (length(std) %nin% c(1L, NCOL(mat))) {
    .err("`std` must have length equal to 1 or the number of columns of `mat`")
  }
  
  mat <- process_mat2(mat, ...)
  
  check_arg_lengths(mat, treat, weights, s.weights, subset)
  
  if (is_null(weights)) weights <- rep.int(1, NROW(mat))
  if (is_null(s.weights)) s.weights <- rep.int(1, NROW(mat))
  
  .chk_numeric(weights)
  .chk_numeric(s.weights)
  
  s.weights <- s.weights / sum(s.weights)
  
  .chk_null_or(subset, .chk_logical)
  if (is_null(subset)) subset <- rep.int(TRUE, NROW(mat))
  
  weights <- weights * s.weights
  
  if (!all(subset)) {
    weights[!subset] <- 0
  }
  
  weights <- weights / sum(weights)
  
  Adist <- abs(outer(treat, treat, "-"))
  Ameans <- colMeans(Adist)
  AA <- Adist + mean(Ameans) - outer(Ameans, Ameans, "+")
  
  if (any(std)) {
    s.d.denom <- .get_s.d.denom.cont(s.d.denom, weights = list(weights),
                                     quietly = TRUE)
    
    weighted.weights <- weighted.weights * s.weights
    weighted.weights <- weighted.weights / sum(weighted.weights)
    
    dvarA <- switch(s.d.denom,
                    "all" = drop(t(s.weights) %*% (AA^2) %*% s.weights),
                    "weighted" = drop(t(weighted.weights) %*% (AA^2) %*% weighted.weights))
  }
  
  out <- vapply(seq_col(mat), function(i) {
    Xdist <- abs(outer(mat[,i], mat[,i], "-"))
    Xmeans <- colMeans(Xdist)
    Xgrand_mean <- mean(Xmeans)
    XX <- Xdist + Xgrand_mean - outer(Xmeans, Xmeans, "+")
    
    if (!std[i]) {
      return(sqrt(drop(t(weights) %*% (XX * AA) %*% weights)))
    }
    
    dvarX <- switch(s.d.denom,
                    "all" = drop(t(s.weights) %*% (XX^2) %*% s.weights),
                    "weighted" = drop(t(weighted.weights) %*% (XX^2) %*% weighted.weights))
    sqrt(drop(t(weights) %*% (XX * AA) %*% weights) / sqrt(dvarX * dvarA))
  }, numeric(1L))
  
  setNames(out, colnames(mat))
}

#' @rdname balance-summary
#' @export 
col_w_dcorr <- function(mat, treat, weights = NULL, s.d.denom = "all",
                        s.weights = NULL, subset = NULL, weighted.weights = weights,
                        na.rm = TRUE, ...) {
  .call <- match.call(expand.dots = TRUE)
  .call[[1L]] <- quote(cobalt::col_w_dcov)
  .call[["std"]] <- quote(TRUE)
  eval.parent(.call)
}

process_mat1 <- function(mat, ...) {
  needs.splitting <- FALSE
  
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (!any_apply(mat, chk::vld_character_or_factor)) {
        return(as.matrix(mat))
      }
      
      needs.splitting <- TRUE
    }
    else if (is.numeric(mat)) {
      return(matrix(mat, ncol = 1L))
    }
    else {
      .err("`mat` must be a data.frame or numeric matrix")
    }
  }
  else if (!is.numeric(mat)) {
    .err("`mat` must be a data.frame or numeric matrix")
  }
  
  if (!needs.splitting) {
    return(mat)
  }
  
  A <- ...mget(setdiff(names(formals(splitfactor)),
                       c("data", "var.name", "drop.first", "drop.level", "split.with")))
  A[["data"]] <- mat
  A[["drop.first"]] <- "if2"
  
  as.matrix(do.call("splitfactor", A))
}
process_mat2 <- function(mat, ..., .bin.vars) {
  if ((!is.numeric(mat) && !is.data.frame(mat)) || length(dim(mat)) > 2L) {
    .err("`mat` must be a data.frame or numeric matrix")
  }
  
  needs.splitting <- FALSE
  
  if (is.data.frame(mat)) {
    to.split <- vapply(mat, chk::vld_character_or_factor, logical(1L))
    
    if (any(to.split)) {
      needs.splitting <- TRUE
    }
    else {
      mat <- as.matrix(mat)
    }
  }
  else if (!is.matrix(mat)) {
    mat <- matrix(mat, ncol = 1L)
  }
  
  bin.vars <- .process_bin_vars(.bin.vars, mat)
  
  if (needs.splitting) {
    bin.vars[to.split] <- TRUE
    
    A <- ...mget(setdiff(names(formals(splitfactor)),
                         c("data", "var.name", "drop.first",
                           "drop.level", "split.with")))
    A[["data"]] <- mat
    A[["drop.first"]] <- "if2"
    A[["split.with"]] <- bin.vars
    
    mat <- do.call("splitfactor", A)
    bin.vars <- attr(mat, "split.with")[[1L]]
    mat <- as.matrix(mat)
  }
  
  attr(mat, "bin") <- bin.vars
  
  mat
}

.process_bin_vars <- function(bin.vars, mat) {
  if (missing(bin.vars)) {
    return(apply(mat, 2L, is_0_1))
  }
  
  if (is_null(bin.vars)) {
    return(rep.int(FALSE, ncol(mat)))
  }
  
  if (is.logical(bin.vars)) {
    
    if (length(bin.vars) != ncol(mat)) {
      .err("if `bin.vars` is logical, it must have length equal to the number of columns of `mat`")
    }
    
    bin.vars[is.na(bin.vars)] <- FALSE
  }
  else if (is.numeric(bin.vars)) {
    bin.vars <- bin.vars[!is.na(bin.vars) & bin.vars != 0]
    
    if (any(bin.vars < 0) && any(bin.vars > 0)) {
      .err("positive and negative indices cannot be mixed with `bin.vars`")
    }
    
    if (any(abs(bin.vars) > ncol(mat))) {
      .err("if `bin.vars` is numeric, none of its values can exceed the number of columns of `mat`")
    }
    
    logical.bin.vars <- rep.int(any(bin.vars < 0), ncol(mat))
    logical.bin.vars[abs(bin.vars)] <- !logical.bin.vars[abs(bin.vars)]
    bin.vars <- logical.bin.vars
  }
  else if (is.character(bin.vars)) {
    bin.vars <- bin.vars[!is.na(bin.vars) & nzchar(bin.vars)]
    
    if (is_null(colnames(mat))) {
      .err("if `bin.vars` is character, `mat` must have column names")
    }
    
    if (!all(bin.vars %in% colnames(mat))) {
      .err("if `bin.vars` is character, all its values must be column names of `mat`")
    }
    
    bin.vars <- colnames(mat) %in% bin.vars
  }
  else {
    .err("`bin.vars` must be a logical, numeric, or character vector")
  }
  
  bin.vars
}
