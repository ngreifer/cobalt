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
#' }
#' 
#' @param mat a numeric matrix or a data frame containing the covariates for which the statistic is to be computed. If a data frame, [splitfactor()] with `drop.first = "if2"` will be called if any character or factor variables are present. This can slow down the function, so it's generally best to supply a numeric matrix. If a numeric vector is supplied, it will be converted to a 1-column matrix first.
#' @param weights `numeric`; an optional set of weights used to compute the weighted statistics. If sampling weights are supplied through `s.weights`, the `weights` should not incorporate these weights, as `weights` and `s.weights` will be multiplied together prior to computing the weighted statistics.
#' @param s.weights `numeric`; an optional set of sampling weights used to compute the weighted statistics. If weights are supplied through `weights`, `weights` and `s.weights` will be multiplied together prior to computing the weighted statistics. Some functions use `s.weights` in a particular way; for others, supplying `weights` and `s.weights` is equivalent to supplying their product to either `weights` or `s.weights`. See Details.
#' @param subset a `logical` vector with length equal to the number of rows of `mat` used to subset the data. See Details for notes on its use with `col_w_smd()`, `col_w_cov()`, and `col_w_corr()`.
#' @param na.rm `logical`; whether `NA`s should be ignored or not. If `FALSE`, any variable with any `NA`s will have its corresponding statistic returned as `NA`. If `TRUE`, any variable with any `NA`s will have its corresponding statistic computed as if the missing value were not there.
#' @param treat a vector of treatment status for each individual. For `col_w_smd()`, `col_w_vr()`, `col_w_ks()`, and `col_w_ovl()`, `treat` should have exactly two unique values. For `col_w_cov()` and `col_w_corr()`, `treat` should be a many-valued numeric vector.
#' @param std `logical`; for `col_w_smd()`, whether the computed mean differences for each variable should be standardized; for `col_w_cov()`, whether treatment-covariate correlations should be computed (`TRUE`) rather than covariances (`FALSE`). Can be either length 1, whereby all variables will be standardized or not, or length equal to the number of columns of `mat`, whereby only variables with a value of `TRUE` will be standardized. See Details.
#' @param s.d.denom for `col_w_smd()` and `col_w_cov()` when `std` is `TRUE` for some variables, and for `col_w_corr()`, how the standardization factor should be computed. For `col_w_smd()` (i.e., when computing standardized mean differences), allowable options include 
#' \itemize{
#'     \item{`"treated"` - uses the standard deviation of the variable in the treated group}
#'     \item{`"control"` - uses the standard deviation of the variable in the control group}
#'     \item{`"pooled"` - uses the square root of the average of the variances of the variable in the treated and control groups}
#'     \item{`"all"` - uses the standard deviation of the variable in the full sample}
#'     \item{`"weighted"` - uses the standard deviation of the variable in the full sample weighted by `weighted.weights`}
#'     \item{`"hedges"` - uses the small-sample corrected version of Hedge's G described in the WWC Procedures Handbook (see References)}
#' \item{the name of one of the treatment values - uses the standard deviation of the variable in that treatment group.}
#' }
#' For `col_w_cov()` and `col_w_corr()`, only `"all"` and `"weighted"` are allowed. Abbreviations allowed. This can also be supplied as a numeric vector of standard deviations with length equal to the number of columns of `mat`; the values will be used as the standardization factors.
#' @param abs `logical`; for `col_w_smd()`, `col_w_cov()`, and `col_w_corr()`, whether the returned statistics should be in absolute value (`TRUE`) or not. For `col_w_vr()`, whether the ratio should always include the larger variance in the numerator, so that the ratio is always greater than or equal to 1. Default is `FALSE`.
#' @param bin.vars a vector used to denote whether each variable is binary or not. Can be a `logical` vector with length equal to the number of columns of `mat` or a vector of numeric indices or character names of the binary variables. If missing (the default), the function will figure out which covariates are binary or not, which can increase computation time. If `NULL`, it will be assumed no variables are binary. All functions other than `col_w_mean()` treat binary variables different from continuous variables. If a factor or character variable is in `mat`, all the dummies created will automatically be marked as binary, but it should still receive an entry when `bin.vars` is supplied as `logical`.
#' @param weighted.weights for `col_w_smd()`, `col_w_cov()`, and `col_w_corr()`, when `std = TRUE` and `s.d.denom = "weighted"`, a vector of weights to be applied to the computation of the denominator standard deviation. If not specified, will use the argument to `weights`. When `s.d.denom` is not `"weighted"`, this is ignored. The main purpose of this is to allow `weights` to be `NULL` while weighting the denominator standard deviations for assessing balance in the unweighted sample but using the standard deviations of the weighted sample.
#' @param type for `col_w_cov()` and `col_w_corr()`, the type of covariance/correlation to be computed. Allowable options include `"pearson"` and `"spearman"`. When `"spearman"` is requested, the covariates and treatment are first turned into ranks using [rank()] with `na.last = "keep"`.
#' @param integrate `logical`; for `col_w_ovl()`, whether to use [integrate()] to calculate the area of overlap or the distance between the densities, respectively. If `FALSE`, a midpoint Riemann sum with 1000 partitions will be used instead. The Riemann sum is a little slower and very slightly imprecise (unnoticibly in most contexts), but the integral can fail sometimes and thus is less stable. The default is to use the Riemann sum.
#' @param ... for all functions, additional arguments supplied to [splitfactor()] when `mat` is a data.frame. `data`, `var.name`, `drop.first`, and `drop.level` are ignored; `drop.first` is automatically set to `"if2"`. For `col_w_ovl()` and `col_w_ent()`, other arguments passed to [density()] besides `x` and `weights`. Note that the default value for `bw` when unspecified is `"nrd"` rather than the default in `density()`, which is `"nrd0"`.
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
#' `col_w_cov()` computes the covariance between a continuous treatment and the covariates to assess balance for continuous treatments as recommended in Austin (2019). These covariance can optionally be weighted or in absolute value or can be requested as correlations (i.e., standardized covariances). The correlations are computed as the covariance between the treatment and covariate divided by a standardization factor, which is equal to the square root of the product of the variance of treatment and the variance of the covariate. The standardization factor is computed using the unweighted variances when `s.weights` are absent, and is computed using the sampling weighted variances when `s.weights` are present, except when `s.d.denom = "weighted"`, in which case the product of `weighted.weights` and `s.weights` (if present) are used to weight the standardization factor. For this reason, the computed correlation can be greater than 1 or less than -1. The standardization factor is always computed using the whole sample even when `subset` is used. The covariance is computed using the product of `weights` and `s.weights`, if specified. The variance of binary variables is computed as \eqn{p(1-p)}, where \eqn{p} is the (weighted) proportion of 1s, while the variance of continuous variables is computed using the standard formula. 
#' 
#' `col_w_corr()` is a wrapper for `col_w_cov` with `std` set to `TRUE`. 
#' 
#' @references 
#' Austin, P. C. (2019). Assessing covariate balance when using the generalized propensity score with quantitative or continuous exposures. *Statistical Methods in Medical Research*, 28(5), 1365–1377. \doi{10.1177/0962280218756159}
#' 
#' Franklin, J. M., Rassen, J. A., Ackermann, D., Bartels, D. B., & Schneeweiss, S. (2014). Metrics for covariate balance in cohort studies of causal effects. *Statistics in Medicine*, 33(10), 1685–1699. \doi{10.1002/sim.6058}
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
#' W <- WeightIt::weightit(treat ~ covs, method = "ps", 
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
#'     row.names = colnames(covs0)
#' ), 4)
#' 
#' # Compare to bal.tab():
#' bal.tab(covs, treat = treat, weights = weights,
#'         disp = c("m", "sd"), stats = c("m", "v", "ks"),
#'         estimand = "ATE", method = "weighting",
#'         binary = "std")
#' 
#' 

#' @name balance-summary
#' @export 
col_w_mean <- function(mat, weights = NULL, s.weights = NULL, subset = NULL, na.rm = TRUE, ...) {
    
    mat <- process_mat1(mat, ...)
    
    check_arg_lengths(mat, weights, s.weights, subset)
    
    if (is_null(weights)) weights <- rep(1, NROW(mat))
    if (is_null(s.weights)) s.weights <- rep(1, NROW(mat))
    
    .chk_null_or(subset, .chk_logical)
    if (is_null(subset)) subset <- rep(TRUE, NROW(mat))
    
    weights <- weights * s.weights
    
    if (sum(weights > 0) < 1) {
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
    
    if (is_null(weights)) weights <- rep(1, NROW(mat))
    if (is_null(s.weights)) s.weights <- rep(1, NROW(mat))
    
    .chk_null_or(subset, .chk_logical)
    if (is_null(subset)) subset <- rep(TRUE, NROW(mat))
    
    weights <- weights * s.weights
    
    if (sum(weights > 0) < 2) {
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
    
    if (is_null(weights)) weights <- rep(1, NROW(mat))
    if (is_null(s.weights)) s.weights <- rep(1, NROW(mat))
    
    .chk_null_or(subset, .chk_logical)
    if (is_null(subset)) subset <- rep(TRUE, NROW(mat))
    
    if (!is_binary(treat[subset])) .err("`treat` must be a binary variable")
    
    weights <- weights * s.weights
    
    if (length(std) == 1L) std <- rep(std, NCOL(mat))
    
    tval1_0 <- treat[1]
    
    if (sum(weights[treat==tval1_0] > 0) < 1 || 
        sum(weights[treat!=tval1_0] > 0) < 1) {
        .err("at least 1 unit in each level of `treat` must have a nonzero weight to compute weighted SMDs")
    }
    
    m1 <- col.w.m(mat[treat == tval1_0 & subset, , drop = FALSE], weights[treat == tval1_0 & subset], na.rm = na.rm)
    m0 <- col.w.m(mat[treat != tval1_0 & subset, , drop = FALSE], weights[treat != tval1_0 & subset], na.rm = na.rm)
    diffs <- m1 - m0
    zeros <- check_if_zero(diffs)
    
    if (any(to.sd <- std & !is.na(zeros) & !zeros)) {
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
    
    if (abs) diffs <- abs(diffs)
    else {
        tval1 <- {
            if (is_0_1(treat))  1
            else get_treated_level(treat[subset])
        }
        if (tval1 != tval1_0) diffs <- -1 * diffs
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
    
    if (is_null(weights)) weights <- rep(1, NROW(mat))
    if (is_null(s.weights)) s.weights <- rep(1, NROW(mat))
    
    .chk_null_or(subset, .chk_logical)
    if (is_null(subset)) subset <- rep(TRUE, NROW(mat))
    
    if (!is_binary(treat[subset])) .err("`treat` must be a binary variable")
    
    weights <- weights * s.weights
    
    weights <- weights[subset]
    treat <- treat[subset]
    mat <- mat[subset, , drop = FALSE]
    
    tval1_0 <- treat[1]
    
    if (sum(weights[treat == tval1_0] > 0) < 2 || 
        sum(weights[treat != tval1_0] > 0) < 2) {
        .err("at least 2 units in each level of `treat` must have nonzero weights to compute weighted variance ratios.")
    }
    
    v1 <- col.w.v(mat[treat==tval1_0, , drop = FALSE], weights[treat==tval1_0], bin.vars = bin.vars, na.rm = na.rm)
    v0 <- col.w.v(mat[treat!=tval1_0, , drop = FALSE], weights[treat!=tval1_0], bin.vars = bin.vars, na.rm = na.rm)
    
    v.ratios <- v1/v0
    
    if (abs) v.ratios <- abs_(v.ratios, ratio = TRUE)
    else {
        tval1 <- {
            if (is_0_1(treat))  1
            else get_treated_level(treat[subset])
        }
        if (tval1 != tval1_0) v.ratios <- 1/v.ratios
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
    
    if (is_null(weights)) weights <- rep(1, NROW(mat))
    if (is_null(s.weights)) s.weights <- rep(1, NROW(mat))
    
    .chk_null_or(subset, .chk_logical)
    if (is_null(subset)) subset <- rep(TRUE, NROW(mat))
    
    if (!is_binary(treat[subset])) .err("`treat` must be a binary variable")
    
    weights <- weights * s.weights
    
    weights <- weights[subset]
    treat <- treat[subset]
    mat <- mat[subset, , drop = FALSE]
    
    tval1 <- treat[1]
    ks <- rep(NA_real_, NCOL(mat))
    
    if (sum(weights[treat == tval1] > 0) < 1 || 
        sum(weights[treat != tval1] > 0) < 1) {
        .err("at least 1 unit in each level of `treat` must have a nonzero weight to compute weighted KS statistics")
    }
    
    if (any(!bin.vars)) {
        weights_ <- weights
        weights_[treat == tval1] <-  weights[treat == tval1]/sum(weights[treat == tval1])
        weights_[treat != tval1] <- -weights[treat != tval1]/sum(weights[treat != tval1])
        ks[!bin.vars] <- apply(mat[, !bin.vars, drop = FALSE], 2, function(x) {
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
        ks[bin.vars] <- abs(col.w.m(mat[treat == tval1, bin.vars, drop = FALSE], weights[treat == tval1], na.rm = na.rm) - 
                                col.w.m(mat[treat != tval1, bin.vars, drop = FALSE], weights[treat != tval1], na.rm = na.rm))
    }
    
    setNames(ks, colnames(mat))
}

#' @rdname balance-summary
#' @export 
col_w_ovl <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars, integrate = FALSE,
                      subset = NULL, na.rm = TRUE, ...) {
    
    A <- list(...)
    
    .chk_not_missing(treat, "`treat`")
    .chk_atomic(treat)
    .chk_not_any_na(treat)
    
    mat <- process_mat2(mat, ..., .bin.vars = bin.vars)
    bin.vars <- attr(mat, "bin")
    
    check_arg_lengths(mat, treat, weights, s.weights, subset)
    
    if (is_null(weights)) weights <- rep(1, NROW(mat))
    if (is_null(s.weights)) s.weights <- rep(1, NROW(mat))
    
    .chk_null_or(subset, .chk_logical)
    if (is_null(subset)) subset <- rep(TRUE, NROW(mat))
    
    if (!is_binary(treat[subset])) .err("`treat` must be a binary variable")
    
    weights <- weights * s.weights
    
    weights <- weights[subset]
    treat <- treat[subset]
    mat <- mat[subset, , drop = FALSE]
    
    tval1 <- treat[1]
    
    .chk_gte(weights, 0)
    
    if (sum(weights[treat == tval1] > 0) < 1 || 
        sum(weights[treat != tval1] > 0) < 1) {
        .err("at least 1 unit in each level of `treat` must have a nonzero weight to compute weighted OVL statistics")
    }
    
    unique.treat <- unique(treat, nmax = 2)
    t.sizes <- setNames(vapply(unique.treat, function(x) sum(treat == x), numeric(1L)),
                        unique.treat)
    smallest.t <- names(t.sizes)[which.min(t.sizes)]
    ovl <- setNames(numeric(ncol(mat)), colnames(mat))
    if (any(!bin.vars)) {
        if (is_null(A[["bw"]])) A[["bw"]] <- "nrd"
        if (!is.function(get0(paste0("bw.", A[["bw"]])))) {
            .err(sprintf("%s is not an acceptable entry to `bw`. See `?stats::density` for allowable options",
                         add_quotes(A[["bw"]])))
        }
        
        A[names(A) %nin% names(formals(density.default))] <- NULL
        
        ovl[!bin.vars] <- apply(mat[, !bin.vars, drop = FALSE], 2, function(cov) {
            if (na.rm) {
                t <- treat[!is.na(cov)]
                w <- weights[!is.na(cov)]
                cov <- cov[!is.na(cov)]
            }
            else if (anyNA(cov)) return(NA_real_)
            else {
                t <- treat
                w <- weights
            }
            
            if (min(cov[t == tval1]) > max(cov[t != tval1]) ||
                min(cov[t != tval1]) > max(cov[t == tval1])) return(1)
            
            cov <- center(cov)/sd(cov)
            
            bw <- get0(paste0("bw.", A[["bw"]]))(cov[t == smallest.t])
            
            f1_ <- approxfun(do.call(density.default,
                                     c(list(cov[t == tval1], 
                                            weights = w[t == tval1] / sum(w[t == tval1]),
                                            bw = bw),
                                       A[setdiff(names(A), "bw")])))
            f1 <- function(x) {
                y <- f1_(x)
                y[is.na(y)] <- 0
                y
            }
            f0_ <- approxfun(do.call(density.default,
                                     c(list(cov[t != tval1], 
                                            weights = w[t != tval1] / sum(w[t != tval1]),
                                            bw = bw),
                                       A[setdiff(names(A), "bw")])))
            f0 <- function(x) {
                y <- f0_(x)
                y[is.na(y)] <- 0
                y
            }
            fn <- function(x) {
                pmin(f1(x), f0(x))
            }
            min.c <- min(cov) - 4 * bw
            max.c <- max(cov) + 4 * bw
            # range <- max.c - min.c
            # min.c.ext <- min.c - .01 * range
            # max.c.ext <- max.c + .01 * range
            if (isTRUE(integrate)) {
                s <- try(integrate(fn, lower = min.c,
                                   upper = max.c)$value,
                         silent = TRUE)
            }
            else {
                s <- intapprox(fn, min.c, max.c, 1001, method = "midpoint")
            }
            
            if (inherits(s, "try-error"))  return(NA_real_)
            
            min(max(1 - s, 0), 1) #Reverse: measure imbalance
            
        })
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
    
    if (is_null(weights)) weights <- rep(1, NROW(mat))
    if (is_null(s.weights)) s.weights <- rep(1, NROW(mat))
    
    .chk_null_or(subset, .chk_logical)
    if (is_null(subset)) subset <- rep(TRUE, NROW(mat))
    
    if (length(std) == 1L) std <- rep(std, NCOL(mat))
    
    .chk_string(type)
    type <- tolower(type)
    type <- match_arg(type, c("pearson", "spearman"))
    if (type == "spearman") {
        for (i in seq_len(ncol(mat))) if (!bin.vars[i]) mat[,i] <- rank(mat[,i], na.last = "keep")
        treat <- rank(treat, na.last = "keep")
    }
    
    weights <- weights * s.weights
    
    if (sum(weights > 0) <= 1) {
        .err("at least 2 units must have nonzero weights to compute weighted covariances")
    }
    
    covars <- col.w.cov(mat[subset, , drop = FALSE], y = treat[subset], w = weights[subset], na.rm = na.rm)
    
    zeros <- check_if_zero(covars)
    
    if (any(to.sd <- std & !is.na(zeros) & !zeros)) {
        if (!is.numeric(s.d.denom) || length(s.d.denom) != ncol(mat)) {
            s.d.denom <- .get_s.d.denom.cont(s.d.denom, weights = list(weights),
                                             quietly = TRUE)
            
            denoms <- .compute_s.d.denom(mat, treat = treat, 
                                         s.d.denom = s.d.denom, s.weights = s.weights, 
                                         bin.vars = bin.vars, subset = subset, to.sd = to.sd,
                                         weighted.weights = weighted.weights, na.rm = na.rm)
        }
        else {
            denoms <- s.d.denom
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
    .call[[1]] <- quote(col_w_cov)
    .call[["std"]] <- quote(TRUE)
    eval.parent(.call)
}


process_mat1 <- function(mat, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else .err("`mat` must be a data.frame or numeric matrix")
    }
    else if (!is.numeric(mat)) {
        .err("`mat` must be a data.frame or numeric matrix")
    }
    
    if (needs.splitting) {
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call("splitfactor", c(list(mat, drop.first ="if2"),
                                        A))
        mat <- as.matrix(mat)
    }
    
    mat
}
process_mat2 <- function(mat, ..., .bin.vars) {
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(to.split <- vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) {
            mat <- matrix(mat, ncol = 1)
        }
        else {
            .err("`mat` must be a data.frame or numeric matrix")
        }
    }
    else if (!is.numeric(mat)) .err("`mat` must be a data.frame or numeric matrix")
    
    bin.vars <- process.bin.vars(.bin.vars, mat)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call("splitfactor", c(list(mat, drop.first ="if2",
                                             split.with = bin.vars),
                                        A))
        bin.vars <- attr(mat, "split.with")[[1]]
        mat <- as.matrix(mat)
    }
    attr(mat, "bin") <- bin.vars
    mat
}

process.bin.vars <- function(bin.vars, mat) {
    if (missing(bin.vars)) return(is_binary_col(mat))
    if (is_null(bin.vars)) return(rep(FALSE, ncol(mat)))
    
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
        logical.bin.vars <- rep(any(bin.vars < 0), ncol(mat))
        logical.bin.vars[abs(bin.vars)] <- !logical.bin.vars[abs(bin.vars)]
        bin.vars <- logical.bin.vars
    }
    else if (is.character(bin.vars)) {
        bin.vars <- bin.vars[!is.na(bin.vars) & bin.vars != ""]
        if (is_null(colnames(mat))) {
            .err("if `bin.vars` is character, `mat` must have column names")
        }
        if (any(bin.vars %nin% colnames(mat))) {
            .err("if `bin.vars` is character, all its values must be column names of `mat`")
        }
        bin.vars <- colnames(mat) %in% bin.vars
    }
    else {
        .err("`bin.vars` must be a logical, numeric, or character vector")
    }
    
    bin.vars
}
