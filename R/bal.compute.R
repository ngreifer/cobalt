#' @title Initialize and efficiently compute scalar balance statistics
#' @name bal.compute
#' 
#' @description These are functions primarily designed for programmers who want to be able to quickly compute one of several scalar (single number) sample balance statistics, e.g., for use in selecting a tuning parameter when estimating balancing weights. `bal.init()` initializes the input so that when `bal.compute()` is used on the output along with a set of weights, the computation of the balance statistic is fast. `vignette("optimizing-balance")` provides an overview and more examples of how to use these functions.
#' 
#' @param stat the name of the statistic to compute. See Details.
#' @param treat a vector containing the treatment variable.
#' @param covs a matrix or data frame containing the covariates.
#' @param s.weights optional; a vector of sampling weights.
#' @param ... other arguments used to specify options for the balance statistic. See Details for which arguments are allowed with each balance statistic.
#' @param init a `bal.init` object created by `bal.init()`.
#' @param weights a vector of balancing weights to compute the weighted statistics
#' 
#' @returns For `bal.init()`, a `bal.init` object containing the components created in the initialization and the function used to compute the balance statistic. For `bal.compute()`, a single numeric value.
#' 
#' @details 
#' The following list contains the allowable balance statistics that can be supplied to `bal.init()`, the additional arguments that can be used with each one, and the treatment types allowed with each one. For all balance statistics, lower values indicate better balance.
#' \describe{
#'     \item{`smd.mean`, `smd.max`, `smd.rms`}{
#'         The mean, maximum, or root-mean-squared absolute standardized mean difference, computed using [col_w_smd()]. The other allowable arguments include `estimand` (ATE, ATC, or ATT) to select the estimand, `focal` to identify the focal treatment group when the ATT is the estimand and the treatment has more than two categories, and `pairwise` to select whether mean differences should be computed between each pair of treatment groups or between each treatment group and the target group identified by `estimand` (default `TRUE`). Can be used with binary and multi-category treatments.
#'     }
#'     \item{`ks.mean`, `ks.max`, `ks.rms`}{
#'         The mean, maximum, or root-mean-squared Kolmogorov-Smirnov statistic, computed using [col_w_ks()]. The other allowable arguments include `estimand` (ATE, ATC, or ATT) to select the estimand, `focal` to identify the focal treatment group when the ATT is the estimand and the treatment has more than two categories, and `pairwise` to select whether statistics should be computed between each pair of treatment groups or between each treatment group and the target group identified by `estimand` (default `TRUE`). Can be used with binary and multi-category treatments.
#'     }
#'     \item{`ovl.mean`, `ovl.max`, `ovl.rms`}{
#'         The mean, maximum, or root-mean-squared overlapping coefficient complement, computed using [col_w_ovl()]. The other allowable arguments include `estimand` (ATE, ATC, or ATT) to select the estimand, `integrate` to select whether integration is done using using [integrate()] (`TRUE`) or a Riemann sum (`FALSE`, the default), `focal` to identify the focal treatment group when the ATT is the estimand and the treatment has more than two categories, `pairwise` to select whether statistics should be computed between each pair of treatment groups or between each treatment group and the target group identified by `estimand` (default `TRUE`). Can be used with binary and multi-category treatments.
#'     }
#'     \item{`mahalanobis`}{
#'         The Mahalanobis distance between the treatment group means. This is similar to `smd.rms` but the covariates are standardized to remove correlations between them and de-emphasize redundant covariates. The other allowable arguments include `estimand` (ATE, ATC, or ATT) to select the estimand and `focal` to identify the focal treatment group when the ATT is the estimand. Can only be used with binary treatments.
#'     }
#'     \item{`energy.dist`}{
#'         The total energy distance between each treatment group and the target sample, which is a scalar measure of the similarity between two multivariate distributions. The other allowable arguments include `estimand` (ATE, ATC, or ATT) to select the estimand, `focal` to identify the focal treatment group when the ATT is the estimand and the treatment has more than two categories, and `improved` to select whether the "improved" energy distance should be used, which emphasizes difference between treatment groups in addition to difference between each treatment group and the target sample (default `TRUE`). Can be used with binary and multi-category treatments.
#'     }
#'     \item{`l1.med`}{
#'         The median L1 statistic computed across a random selection of possible coarsening of the data. The other allowable arguments include `l1.min.bin` (default 2) and `l1.max.bin` default (12) to select the minimum and maximum number of bins with which to bin continuous variables and `l1.n` (default 101) to select the number of binnings used to select the binning at the median. `covs` should be supplied without splitting factors into dummies to ensure the binning works correctly. Can be used with binary and multi-category treatments.
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
#'     \item{`distance.cov`}{
#'         The distance covariance between the scaled covariates and treatment, which is a scalar measure of the independence of two possibly multivariate distributions. Can only be used with continuous treatments.
#'     }
#' }
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
#'     return(m$weights)
#' }
#' 
#' # Initialize balance statistic; largest KS statistic
#' init <- bal.init("ks.max", treat = lalonde$treat,
#'                  covs = lalonde[covs],
#'                  estimand = "ATE")
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

#' @export
bal.init <- function(stat, treat, covs, s.weights = NULL, ...) {
    chk::chk_not_missing(stat)
    chk::chk_not_missing(treat)
    chk::chk_not_missing(covs)

    chk::chk_string(stat)
    chk::chk_vector(treat)
    
    chk::chk_null_or(s.weights, vld = chk::chk_numeric)

    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    init <- bal_criterion(treat.type, stat)
    
    out <- init$init(treat = treat, covs = covs, s.weights = s.weights, ...)
    
    attr(out, "fun") <- init$fun
    attr(out, "treat.type") <- attr(init, "treat.type")
    attr(out, "stat") <- attr(init, "stat")
    
    class(out) <- c("bal.init", class(out))
    
    out
}

#' @rdname bal.compute
#' @export
bal.compute <- function(init, weights = NULL) {
    chk::chk_is(init, "bal.init")

    fun <- attr(init, "fun")
    
    fun(init = init, weights = weights)
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
                     # "kernel.dist" = "kernel distance",
                     "l1.med" = "L1 median",
                     "r2" = "post-weighting treatment R-squared",
                     "p.mean" = "average Pearson correlation",
                     "p.max" = "maximum Pearson correlation",
                     "p.rms" = "root-mean-square Pearson correlation",
                     "s.mean" = "average Spearman correlation",
                     "s.max" = "maximum Spearman correlation",
                     "s.rms" = "root-mean-square Spearman correlation",
                     "distance.cov" = "distance covariance",
                     NA_character_
    )
    
    if (anyNA(phrase)) {
        .err(sprintf("%s is not an allowed statistic", add_quotes(stat, 2)))
    }
    
    return(phrase)
}

process_init_covs <- function(covs) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else .err("`covs` must be a data.frame or numeric matrix")
    }
    else if (!is.numeric(covs)) .err("`covs` must be a data.frame or numeric matrix")
    
    bin.vars <- process.bin.vars(mat = covs)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call(splitfactor, list(covs, drop.first ="if2",
                                          split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }
    
    attr(covs, "bin") <- bin.vars
    covs
}

#init functions
init_smd <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, ...) {
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    chk::chk_flag(pairwise)
    
    covs <- process_init_covs(covs)
    bin.vars <- attr(covs, "bin")

    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
        .err("`treat` must be a binary or multi-category variable")
    }
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    unique.treats <- unique(treat)
    
    if (treat.type == "multinomial") {
        if (is_null(focal) && !pairwise) {
            treat.all <- last(make.unique(unique.treats, "All"))
            treat <- factor(c(as.character(treat), rep(treat.all, length(treat))),
                            levels = c(unique.treats, treat.all))
            covs <- rbind(covs, covs)
            s.weights <- rep(s.weights, 2)
            focal <- treat.all
        }
        
        if (is_null(focal) || pairwise) {
            treatment.pairs <- combn(unique.treats, 2, simplify = FALSE)
        }
        else {
            treatment.pairs <- lapply(setdiff(unique.treats, focal), function(x) c(x, focal))
        }
    }
    else {
        treatment.pairs <- list(unique.treats)
        pairwise <- TRUE
    }
    
    s.d.denom <- get.s.d.denom(estimand = estimand, treat = treat, focal = focal, quietly = TRUE)
    
    denoms <- compute_s.d.denom(covs, treat = treat,
                                s.d.denom = s.d.denom, s.weights = s.weights,
                                bin.vars = bin.vars)
    
    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                s.d.denom = denoms,
                focal = focal,
                pairwise = pairwise,
                treatment.pairs = treatment.pairs)
    class(out) <- "init_smd"
    out
}
init_ks <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, ...) {
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    chk::chk_flag(pairwise)
    
    covs <- process_init_covs(covs)
    bin.vars <- attr(covs, "bin")
    
    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
        .err("`treat` must be a binary or multi-category variable")
    }
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    unique.treats <- unique(treat)
    
    if (treat.type == "multinomial") {
        if (is_null(focal) && !pairwise) {
            treat.all <- last(make.unique(unique.treats, "All"))
            treat <- factor(c(as.character(treat), rep(treat.all, length(treat))),
                            levels = c(unique.treats, treat.all))
            covs <- rbind(covs, covs)
            s.weights <- rep(s.weights, 2)
            focal <- treat.all
        }
        
        if (is_null(focal) || pairwise) {
            treatment.pairs <- combn(unique.treats, 2, simplify = FALSE)
        }
        else {
            treatment.pairs <- lapply(setdiff(unique.treats, focal), function(x) c(x, focal))
        }
    }
    else {
        treatment.pairs <- list(unique.treats)
        pairwise <- TRUE
    }
    
    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                focal = focal,
                pairwise = pairwise,
                treatment.pairs = treatment.pairs)
    class(out) <- "init_ks"
    out
}
init_ovl <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE,
                     integrate = FALSE, ...) {
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    chk::chk_flag(pairwise)
    
    covs <- process_init_covs(covs)
    bin.vars <- attr(covs, "bin")

    chk::chk_flag(integrate)

    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
        .err("`treat` must be a binary or multi-category variable")
    }
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    unique.treats <- unique(treat)
    
    if (treat.type == "multinomial") {
        if (is_null(focal) && !pairwise) {
            treat.all <- last(make.unique(unique.treats, "All"))
            treat <- factor(c(as.character(treat), rep(treat.all, length(treat))),
                            levels = c(unique.treats, treat.all))
            covs <- rbind(covs, covs)
            s.weights <- rep(s.weights, 2)
            focal <- treat.all
        }
        
        if (is_null(focal) || pairwise) {
            treatment.pairs <- combn(unique.treats, 2, simplify = FALSE)
        }
        else {
            treatment.pairs <- lapply(setdiff(unique.treats, focal), function(x) c(x, focal))
        }
    }
    else {
        treatment.pairs <- list(unique.treats)
        pairwise <- TRUE
    }
    
    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                focal = focal,
                pairwise = pairwise,
                treatment.pairs = treatment.pairs,
                integrate = integrate)
    class(out) <- "init_ovl"
    out
}
init_mahalanobis <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, ...) {
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    covs <- process_init_covs(covs)
    bin.vars <- attr(covs, "bin")
    
    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
        .err("`treat` must be a binary or multi-category variable")
    }
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    s.d.denom <- get.s.d.denom(estimand = estimand, treat = treat, focal = focal, quietly = TRUE)
    
    if (any(!bin.vars)) covs[,!bin.vars] <- scale(covs[,!bin.vars])
    
    if (s.d.denom %in% as.character(treat)) {
        sigma <- cov.wt(covs[treat == s.d.denom,,drop = FALSE], s.weights[treat == s.d.denom])$cov
        if (any(zeros <- diag(sigma) == 0)) {
            sigma_all <- cov.wt(covs, s.weights)$cov
            sigma[zeros,] <- sigma_all[zeros,]
            sigma[,zeros] <- sigma_all[,zeros]
            if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros]  <- 1
        }
    }
    else if (s.d.denom == "pooled") {
        sigma <- .5 * (cov.wt(covs[treat == treat[1],,drop = FALSE], s.weights[treat == treat[1]])$cov +
                           cov.wt(covs[treat != treat[1],,drop = FALSE], s.weights[treat != treat[1]])$cov)
        if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros] <- 1
    }
    else {
        sigma <- cov.wt(covs, s.weights)$cov
        if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros] <- 1
    }
    
    #MASS::ginv
    sigmasvd <- svd(sigma)
    pos <- sigmasvd$d > max(1e-8 * sigmasvd$d[1L], 0)
    sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% ((1/sigmasvd$d[pos]) *
                                                          t(sigmasvd$u[, pos, drop = FALSE]))
    
    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                s.d.denom = s.d.denom,
                sigma_inv = sigma_inv)
    class(out) <- "init_mahalanobis"
    out
}
init_energy.dist <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE, ...) {
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    covs <- process_init_covs(covs)
    bin.vars <- attr(covs, "bin")
    
    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
        .err("`treat` must be a binary or multi-category variable")
    }
    
    treat <- factor(treat)
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    dist.covs <- scale(covs, scale = sqrt(col.w.v(covs, s.weights, bin.vars)))
    
    d <- unname(as.matrix(dist(dist.covs)))
    
    n <- length(treat)
    unique.treats <- levels(treat)
    treat.seq <- seq_along(unique.treats)
    diagn <- diag(n)
    
    for (t in unique.treats) {
        s.weights[treat == t] <- s.weights[treat == t]/mean_fast(s.weights[treat == t])
    }

    treat_t <- vapply(unique.treats, function(t) treat == t, logical(n))
    n_t <- colSums(treat_t)
    
    s.weights_n_t <- setNames(lapply(unique.treats, function(t) treat_t[,t] * s.weights / n_t[t]),
                              unique.treats)
   
    if (is_null(focal)) {
        
        P <- -d * Reduce("+", lapply(s.weights_n_t, tcrossprod))
        
        q <- ((s.weights * 2 / n) %*% d) * Reduce("+", s.weights_n_t)
        
        if (improved) {
            all_pairs <- combn(unique.treats, 2, simplify = FALSE)
            P <- P - d * Reduce("+", lapply(all_pairs, function(p) {
                tcrossprod(s.weights_n_t[[p[1]]] - s.weights_n_t[[p[2]]])
            }))
        }
    }
    else {
        non_focal <- setdiff(unique.treats, focal)
        in_focal <- treat == focal
        
        P <- -d[!in_focal, !in_focal] *
            Reduce("+", lapply(s.weights_n_t[non_focal], function(x) tcrossprod(x[!in_focal])))
        
        q <- 2 * (s.weights_n_t[[focal]][in_focal] %*% d[in_focal, !in_focal]) *
            Reduce("+", lapply(s.weights_n_t[non_focal], function(x) x[!in_focal]))
    }
    
    out <- list(q = q,
                P = P,
                s.weights = s.weights,
                treat = treat,
                unique.treats = unique.treats,
                focal = focal)
    class(out) <- "init_energy.dist"
    out
}
init_p <- function(covs, treat, s.weights = NULL, ...) {
    covs <- process_init_covs(covs)
    bin.vars <- attr(covs, "bin")
    
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("continuous")) {
        .err("`treat` must be a continuous (numeric) variable")
    }
    
    s.d.denom <- get.s.d.denom.cont(quietly = TRUE)
    
    denoms <- compute_s.d.denom(covs, treat = treat,
                                s.d.denom = s.d.denom, s.weights = s.weights,
                                bin.vars = bin.vars)
    
    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                s.d.denom = denoms)
    class(out) <- "init_p"
    out
}
init_s <- function(covs, treat, s.weights = NULL, ...) {
    covs <- process_init_covs(covs)
    bin.vars <- attr(covs, "bin")
    
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("continuous")) {
        .err("`treat` must be a continuous (numeric) variable")
    }
    
    for (i in seq_len(ncol(covs))[!bin.vars[i]]) {
        covs[,i] <- rank(covs[,i], na.last = "keep")
    }
    treat <- rank(treat, na.last = "keep")
    
    s.d.denom <- get.s.d.denom.cont(quietly = TRUE)
    
    denoms <- compute_s.d.denom(covs, treat = treat,
                                s.d.denom = s.d.denom, s.weights = s.weights,
                                bin.vars = bin.vars)
    
    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                s.d.denom = denoms)
    class(out) <- "init_s"
    out
}
init_r2 <- function(covs, treat, s.weights = NULL, poly = 1, int = FALSE, ...) {
    covs <- process_init_covs(covs)
    bin.vars <- attr(covs, "bin")
    
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "continuous")) {
        .err("`treat` must be a binary or continuous (numeric) variable")
    }
    
    if (treat.type == "binary") treat <- as.numeric(treat == treat[1])
    
    x <- cbind(`(Intercept)` = 1, covs, int.poly.f2(covs, poly = poly, int = int))
    
    out <- list(treat = treat,
                x = x,
                s.weights = s.weights)
    class(out) <- "init_r2"
    out
}
init_distance.cov <- function(covs, treat, s.weights = NULL, ...) {
    covs <- process_init_covs(covs)
    bin.vars <- attr(covs, "bin")
    
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("continuous")) {
        .err("`treat` must be a continuous (numeric) variable")
    }
    
    dist.covs <- scale(covs, scale = sqrt(col.w.v(covs, s.weights, bin.vars)))
    
    Xdist <- unname(as.matrix(dist(dist.covs)))
    
    n <- length(treat)
    
    Adist <- unname(as.matrix(dist(treat/sqrt(col.w.v(treat, s.weights)))))
    
    Xmeans <- colMeans(Xdist)
    Xgrand_mean <- mean(Xmeans)
    XA <- Xdist + Xgrand_mean - outer(Xmeans, Xmeans, "+")
    
    Ameans <- colMeans(Adist)
    Agrand_mean <- mean(Ameans)
    AA <- Adist + Agrand_mean - outer(Ameans, Ameans, "+")
    
    P <- XA * AA/n^2
    
    out <- list(P = P,
                s.weights = s.weights,
                treat = treat)
    class(out) <- "init_distance.cov"
    out
}
init_l1.med <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL,
                        .covs = NULL, l1.min.bin = 2, l1.max.bin = 12, l1.n = 101, ...) {
    
    if (is_not_null(.covs)) covs <- .covs
    if (!is.data.frame(covs)) {
        if (is.atomic(covs) && is_null(dim(covs))) covs <- data.frame(covs)
        else if (!is.matrix(covs)) .err("`covs` must be a data.frame or matrix.")
    }
    covs <- as.data.frame(covs)
    
    chk::chk_not_missing(treat)
    chk::chk_atomic(treat)
    
    check_arg_lengths(covs, treat, s.weights)
    
    if (is_null(s.weights)) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %nin% c("binary", "multinomial")) {
        .err("`treat` must be a binary or multi-category variable")
    }
    
    coarsen <- function(covs, cutpoints = NULL, grouping = NULL) {
        is.numeric.cov <- setNames(vapply(covs, is.numeric, logical(1L)), names(covs))
        for (i in names(cutpoints)) {
            if (cutpoints[[i]] == 0) is.numeric.cov[i] <- FALSE #Will not be binned
        }
        
        #Process grouping
        if (!is.null(grouping) && !is.null(names(grouping))) {
            covs[names(grouping)] <- lapply(names(grouping), function(g) {
                x <- covs[[g]]
                groups <- grouping[[g]]
                
                for (i in seq_along(groups)) {
                    x[x %in% groups[[i]]] <- groups[[i]][1]
                }
                x
            })
            cutpoints[names(cutpoints) %in% names(grouping)] <- NULL
        }
        
        #Create bins for numeric variables
        for (i in names(covs)[is.numeric.cov]) {
            bins <- cutpoints[[i]]
            
            #cutpoints is number of bins, unlike in cem
            breaks <- seq(min(covs[[i]]), max(covs[[i]]), length = bins + 1)
            breaks[c(1, bins + 1)] <- c(-Inf, Inf)
            
            covs[[i]] <- findInterval(covs[[i]], breaks)
        }
        
        #Reduce to strata
        factor(do.call("paste", c(covs, sep = " | ")))
    }
    
    f.e <- process_focal_and_estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
    
    unique.treats <- unique(treat)
    for (t in unique.treats) s.weights[treat == t] <- s.weights[treat == t]/sum(s.weights[treat == t])
    
    is.numeric.cov <- setNames(vapply(covs, is.numeric, logical(1L)), names(covs))
    nunique.covs <- vapply(covs, nunique, integer(1L))
    
    coarsenings <- lapply(1:l1.n, function(i) {
        cutpoints <- setNames(lapply(nunique.covs[is.numeric.cov], function(nu) {
            sample(seq(min(l1.min.bin, nu), min(l1.max.bin, nu)), 1)
        }), names(covs)[is.numeric.cov])
        grouping <- setNames(lapply(seq_along(covs)[!is.numeric.cov], function(x) {
            nu <- nunique.covs[x]
            u <- unique(covs[[x]], nmax = nu)
            
            #Randomly select number of bins
            nbins <- sample(seq(min(l1.min.bin, nu), min(l1.max.bin, nu)), 1)
            
            #Randomly assign bin numbers to levels of covariate
            bin.assignments <- sample(seq_len(nbins), nu, replace = TRUE)
            
            #Group levels with same bin number
            lapply(unique(bin.assignments, nmax = nbins),
                   function(b) u[bin.assignments == b])
            
        }), names(covs)[!is.numeric.cov])
        list(cutpoints = cutpoints, grouping = grouping, treat_cutpoints = NULL)
    })
    
    l1s <- unlist(lapply(coarsenings, function(co) {
        x <- coarsen(covs, cutpoints = co[["cutpoints"]], grouping = co[["grouping"]])
        
        if (treat.type == "binary" || is_null(focal)) {
            sum(vapply(levels(x), function(l) {
                in_l <- which(x == l)
                abs(diff(range(vapply(unique.treats, function(t) sum(s.weights[in_l][treat[in_l] == t]), numeric(1L)))))
            }, numeric(1L))) / length(unique.treats)
        }
        else {
            sum(vapply(levels(x), function(l) {
                in_l <- which(x == l)
                sum.s.weights.focal <- sum(s.weights[in_l][treat[in_l] == focal])
                max(abs(vapply(unique.treats[unique.treats != focal],
                               function(t) sum(s.weights[in_l][treat[in_l] == t]) - sum.s.weights.focal, numeric(1L))))
            }, numeric(1L))) / length(unique.treats)
        }
    }))
    
    l1.med <- sort(l1s, partial = ceiling(l1.n/2))[ceiling(l1.n/2)]
    
    l1.med.coarsening <- coarsenings[[which(l1s == l1.med)[1]]]
    
    out <- list(coarsened.covs = coarsen(covs,
                                         cutpoints = l1.med.coarsening[["cutpoints"]],
                                         grouping = l1.med.coarsening[["grouping"]]),
                s.weights = s.weights,
                treat = treat,
                unique.treats = unique.treats,
                focal = focal)
    class(out) <- "init_l1.med"
    out
    
}

#Statistics
smd.binary <- function(init, weights = NULL) {
    check_init(init, "init_smd")
    col_w_smd(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
              bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE)
}
ks.binary <- function(init, weights = NULL) {
    check_init(init, "init_ks")
    col_w_ks(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
             bin.vars = init$bin.vars)
}
ovl.binary <- function(init, weights = NULL) {
    check_init(init, "init_ovl")
    col_w_ovl(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
              bin.vars = init$bin.vars, integrate = init$integrate)
}
mahalanobis.binary <- function(init, weights = NULL) {
    check_init(init, "init_mahalanobis")
    mean.diffs <- col_w_smd(init$covs, init$treat, weights, s.weights = init$s.weights,
                            bin.vars = init$bin.vars, std = FALSE)
    drop(sqrt(t(mean.diffs) %*% init$sigma_inv %*% mean.diffs))
}
energy.dist.binary <- function(init, weights = NULL) {
    check_init(init, "init_energy.dist")
    
    if (is_null(weights)) weights <- init[["s.weights"]]
    else weights <- weights * init[["s.weights"]]
    
    for (t in init[["unique.treats"]]) {
        weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/mean_fast(weights[init[["treat"]] == t])
    }
    
    if (is_not_null(init[["focal"]])) {
        weights <- weights[init[["treat"]] != init[["focal"]]]
    }
    
    return(drop(weights %*% init[["P"]] %*% weights + init[["q"]] %*% weights))
}
r2.binary <- function(init, weights = NULL) {
    check_init(init, "init_r2")
    
    if (is_null(weights)) weights <- init[["s.weights"]]
    else weights <- weights * init[["s.weights"]]
    
    fit <- glm.fit(init$x, init$treat, weights, family = quasibinomial())
    
    wmtreat <- sum(weights*fit$linear.predictors)/sum(weights)
    
    SSmodel <- sum(weights * (fit$linear.predictors - wmtreat)^2)
    
    SSmodel / (sum(weights) * pi^2/3 + SSmodel)
}
l1.med.binary <- function(init, weights = NULL) {
    check_init(init, "init_l1.med")
    
    if (is_null(weights)) weights <- init[["s.weights"]]
    else weights <- weights * init[["s.weights"]]
    
    for (t in init[["unique.treats"]]) {
        weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/sum(weights[init[["treat"]] == t])
    }
    
    x <- init[["coarsened.covs"]]
    
    sum(vapply(levels(x), function(l) {
        in_l <- which(x == l)
        abs(diff(vapply(init[["unique.treats"]], function(t) sum(weights[in_l][init[["treat"]][in_l] == t]), numeric(1L))))
    }, numeric(1L))) / 2
    
}
smd.multinomial <- function(init, weights = NULL) {
    check_init(init, "init_smd")
    
    if (!init$pairwise) {
        weights <- c(weights, rep(1, length(weights)))
    }
    
    unlist(lapply(init$treatment.pairs, function(x) {
        col_w_smd(init$covs[init$treat %in% x,,drop = FALSE],
                  treat = init$treat[init$treat %in% x],
                  weights = weights[init$treat %in% x],
                  s.weights = init$s.weights[init$treat %in% x],
                  bin.vars = init$bin.vars,
                  s.d.denom = init$s.d.denom, abs = TRUE)
    }))
}
ks.multinomial <- function(init, weights = NULL) {
    check_init(init, "init_ks")
    
    if (!init$pairwise) {
        weights <- c(weights, rep(1, length(weights)))
    }
    
    unlist(lapply(init$treatment.pairs, function(x) {
        col_w_ks(init$covs[init$treat %in% x,,drop = FALSE],
                 treat = init$treat[init$treat %in% x],
                 weights = weights[init$treat %in% x],
                 s.weights = init$s.weights[init$treat %in% x],
                 bin.vars = init$bin.vars)
    }))
}
ovl.multinomial <- function(init, weights = NULL) {
    check_init(init, "init_ks")
    
    if (!init$pairwise) {
        weights <- c(weights, rep(1, length(weights)))
    }
    
    unlist(lapply(init$treatment.pairs, function(x) {
        col_w_ovl(init$covs[init$treat %in% x,,drop = FALSE],
                 treat = init$treat[init$treat %in% x],
                 weights = weights[init$treat %in% x],
                 s.weights = init$s.weights[init$treat %in% x],
                 bin.vars = init$bin.vars,
                 integrate = init$integrate)
    }))
}
energy.dist.multinomial <- energy.dist.binary
l1.med.multinomial <- function(init, weights = NULL) {
    check_init(init, "init_l1.med")
    
    if (is_null(weights)) weights <- init[["s.weights"]]
    else weights <- weights * init[["s.weights"]]
    
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
pearson.corr.continuous <- function(init, weights = NULL) {
    check_init(init, "init_p")
    col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
              bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
              std = TRUE)
}
spearman.corr.continuous <- function(init, weights = NULL) {
    check_init(init, "init_s")
    col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
              bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
              std = TRUE)
}
r2.continuous <- function(init, weights = NULL) {
    check_init(init, "init_r2")
    
    if (is_null(weights)) weights <- init[["s.weights"]]
    else weights <- weights * init[["s.weights"]]
    
    fit <- lm.wfit(init$x, init$treat, weights)
    
    SSresid <- sum(weights * (fit$residuals - w.m(fit$residuals, weights))^2)
    SStreat <- sum(weights * (init$treat - w.m(init$treat, weights))^2)
    
    1 - SSresid/SStreat
}
distance.cov.continuous <- function(init, weights = NULL) {
    check_init(init, "init_distance.cov")
    
    if (is_null(weights)) weights <- init[["s.weights"]]
    else weights <- weights * init[["s.weights"]]
    
    weights <- weights/mean_fast(weights)
    
    return(drop(t(weights) %*% init[["P"]] %*% weights))
}

bal_criterion <- function(treat.type, criterion) {
    chk::chk_not_missing(criterion)
    chk::chk_not_missing(treat.type)

    treat.type <- match_arg(treat.type, c("binary", "multinomial", "continuous"))
    criteria <- switch(
        treat.type,
        binary = list(
            smd.mean = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_smd(covs, treat, s.weights, estimand)
                    }
                    mean_fast(smd.binary(init, weights))
                },
                init = init_smd
            ),
            smd.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_smd(covs, treat, s.weights, estimand)
                    }
                    max(smd.binary(init, weights))
                },
                init = init_smd
            ),
            smd.rms = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
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
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, integrate = FALSE,
                               init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ovl(covs, treat, s.weights, integrate = integrate)
                    }
                    mean_fast(ovl.binary(init, weights))
                },
                init = init_ovl
            ),
            ovl.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, integrate = FALSE,
                               init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ovl(covs, treat, s.weights, integrate = integrate)
                    }
                    max(ovl.binary(init, weights))
                },
                init = init_ovl
            ),
            ovl.rms = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, integrate = FALSE,
                               init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ovl(covs, treat, s.weights, integrate = integrate)
                    }
                    rms(ovl.binary(init, weights))
                },
                init = init_ovl
            ),
            mahalanobis = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_mahalanobis(covs, treat, s.weights, estimand)
                    }
                    mahalanobis.binary(init, weights)
                },
                init = init_mahalanobis
            ),
            energy.dist = list(
                fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE, init = NULL, ...) {
                    
                    if (is_null(init)) {
                        init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                    }
                    energy.dist.binary(init, weights)
                },
                init = init_energy.dist
            ),
            # kernel.dist = list(
            #     fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, init = NULL, ...) {
            #
            #         if (is_null(init)) {
            #             init <- init_kernel.dist(covs, treat, s.weights, estimand, focal)
            #         }
            #         kernel.dist.binary(init, weights)
            #     },
            #     init = init_kernel.dist
            # ),
            l1.med = list(
                fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, init = NULL, ...) {
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
        ),
        multinomial = list(
            smd.mean = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_smd(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                    }
                    mean_fast(smd.multinomial(init, weights))
                },
                init = init_smd
            ),
            smd.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_smd(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                    }
                    max(smd.multinomial(init, weights))
                },
                init = init_smd
            ),
            smd.rms = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_smd(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                    }
                    rms(smd.multinomial(init, weights))
                },
                init = init_smd
            ),
            ovl.mean = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                               pairwise = TRUE, integrate = FALSE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ovl(covs, treat, s.weights, focal = focal, pairwise = pairwise,
                                         integrate = integrate)
                    }
                    mean_fast(ovl.multinomial(init, weights))
                },
                init = init_ks
            ),
            ovl.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                               pairwise = TRUE, integrate = FALSE, init = NULL, ...) {
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
                               pairwise = TRUE, integrate = FALSE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ovl(covs, treat, s.weights, focal = focal, pairwise = pairwise,
                                         integrate = integrate)
                    }
                    rms(ovl.multinomial(init, weights))
                },
                init = init_ovl
            ),
            energy.dist = list(
                fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE, init = NULL, ...) {
                    
                    if (is_null(init)) {
                        init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                    }
                    energy.dist.multinomial(init, weights)
                },
                init = init_energy.dist
            ),
            l1.med = list(
                fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_l1.med(covs, treat, s.weights, estimand, focal, ...)
                    }
                    l1.med.multinomial(init, weights)
                },
                init = init_l1.med
            )
        ),
        
        continuous = list(
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
            )
        )
    )
    criterion <- match_arg(criterion, names(criteria))
    
    criteria[[criterion]]
}

check_init <- function(init, init_class) {
    chk::chk_not_missing(init)
    chk::chk_not_missing(init_class)
    chk::chk_is(init, init_class)
}

