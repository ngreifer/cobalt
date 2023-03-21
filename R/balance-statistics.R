#' Balance Statistics in `bal.tab` and `love.plot`
#' @name balance-statistics
#' 
#' @description [bal.tab()] and [love.plot()] display balance statistics for the included covariates. The `stats` argument in each of these functions controls which balance statistics are to be displayed. The argument to `stats` should be a character vector with the names of the desired balance statistics.
#' 
#' This page describes all of the available balance statistics and how to request them. Abbreviations are allowed, so you can use the first few letters of each balance statistics to request it instead of typing out its whole name. That convention is used throughout the documentation. For example, to request mean differences and variance ratios in `bal.tab()` or `love.plot()`, you could include `stats = c("m", "v")`. In addition, the `thresholds` argument uses the same naming conventions and can be used to request balance thresholds on each statistic. For example, to request a balance threshold of .1 for mean differences, you could include `thresholds = c(m = .1)`.
#'     
#' Below, each allowable entry to `stats` and `thresholds` are described, along with other details or option that accompany them.
#'     
#' ## Binary/Multi-Category Treatments
#' \describe{
#'    \item{`"mean.diffs"`}{Mean differences as computed by [col_w_smd()]. Can be abbreviated as `"m"`. Setting the arguments `continuous` and `binary` to either `"std"` or `"raw"` will determine whether standardized mean differences or raw mean differences are calculated for continuous and categorical variables, respectively. When standardized mean differences are requested, the `s.d.denom` argument controls how the standardization occurs. When `abs = TRUE`, negative values become positive. Mean differences are requested by default when no entry to `stats` is provided.}
#'             
#'   \item{`"variance.ratios"`}{Variance ratios as computed by [col_w_vr()]. Can be abbreviated as `"v"`. Will not be computed for binary variables. When `abs = TRUE`, values less than 1 will have their inverse taken. When used with `love.plot`, the x-axis scaled will be logged so that, e.g., .5 is as far away from 1 as 2 is.}
#'             
#'   \item{`"ks.statistics"`}{Kolmogorov-Smirnov (KS) statistics as computed by [col_w_ks()].}
#'             
#'   \item{`"ovl.coefficients"`}{Overlapping (OVL) statistics as computed by [col_w_ovl()]. Can be abbreviated as `"ovl"`. Additional arguments passed to `col_w_ovl()`, such as `integrate` or `bw`, can be supplied to `bal.tab()` or `love.plot()`.}
#' }
#'     
#' ## Continuous Treatments
#' \describe{
#'             \item{`"correlations"`}{Pearson correlations as computed by [col_w_cov()]. Can be abbreviated as `"cor"`. Setting the arguments `continuous` and `binary` to either `"std"` or `"raw"` will determine whether correlations or covariances are calculated for continuous and categorical variables, respectively (they are both `"std"` by default). When correlations are requested, the `s.d.denom` argument controls how the standardization occurs. When `abs = TRUE`, negative values become positive. Pearson correlations are requested by default when no entry to `stats` is provided.}
#'             
#' \item{`"spearman.correlations"`}{Spearman correlations as computed by [col_w_cov()]. Can be abbreviated as `"sp"`. All arguments are the same as those for `"correlations"`. When `abs = TRUE`, negative values become positive.}
#'             
#' \item{`"mean.diffs.target"`}{Mean differences computed between the weighted and unweighted sample to ensure the weighted sample is representative of the original population. Can be abbreviated as `"m"`. Setting the arguments `continuous` and `binary` to either `"std"` or `"raw"` will determine whether standardized mean differences or raw mean differences are calculated for continuous and categorical variables, respectively. The standardization factor will be computed in the unweighted sample. When `abs = TRUE`, negative values become positive. This statistic is only computed for the adjusted samples.}
#'             
#' \item{`"ks.statistics.target"`}{KS-statistics computed between the weighted and unweighted sample to ensure the weighted sample is representative of the original population. Can be abbreviated as `"ks"`. This statistic is only computed for the adjusted samples.}
#' }
#'     
#' If a statistic is requested in `thresholds`, it will automatically be placed in `stats`. For example, `bal.tab(..., stats = "m", thresholds = c(v = 2))` will display both mean differences and variance ratios, and the variance ratios will have a balance threshold set to 2.
#' 
#' @examples
#' data(lalonde)
#' 
#' #Binary treatments
#' bal.tab(treat ~ age + educ + married + re74, data = lalonde,
#'         stats = c("m", "v", "ks"))
#' love.plot(treat ~ age + educ + married + re74, data = lalonde,
#'           stats = c("m", "v", "ks"), binary = "std",
#'           thresholds = c(m = .1, v = 2))
#' 
#' #Continuous treatments
#' bal.tab(re75 ~ age + educ + married + re74, data = lalonde,
#'         stats = c("cor", "sp"))
#' love.plot(re75 ~ age + educ + married + re74, data = lalonde,
#'           thresholds = c(cor = .1, sp = .1))
#' 
NULL

get_from_STATS <- function(what) {
    setNames(sapply(STATS, function(s) s[[what]], USE.NAMES = FALSE),
             names(STATS))
}

STATS <- list()

STATS[["mean.diffs"]] <- {list(
    type = "bin",
    threshold = "m.threshold",
    Threshold = "M.Threshold",
    disp_stat = "disp.diff",
    adj_only = FALSE,
    abs = function(x) abs_(x),
    bal.tab_column_prefix = "Diff", #Also which.stat in love.plot
    threshold_range = c(0, Inf),
    balance_tally_for = "mean differences",
    variable_with_the_greatest = "mean difference", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        A <- list(...)
        binary <- A$binary #attr(x, "print.options")$binary
        continuous <- A$continuous #attr(x, "print.options")$continuous
        abs <- A$abs
        var_type <- A$var_type #B[["type"]]
        stars <- A$stars
        
        #All std, no std, some std
        if ((binary == "std" || !any(var_type == "Binary")) && 
            (continuous == "std" || !any(var_type == "Contin."))) {
            xlab.diff <- "Standardized Mean Differences"
        } 
        else if ((binary == "raw" || !any(var_type == "Binary")) && 
                 (continuous == "raw" || !any(var_type == "Contin."))) {
            xlab.diff <- "Mean Differences"
        }
        else {
            stars <- match_arg(stars, c("none", "std", "raw"))
            if (stars == "none") {
                xlab.diff <- "Mean Differences"
            }
            else if (stars == "std") {
                xlab.diff <- "Mean Differences"
            }
            else if (stars == "raw") {
                xlab.diff <- "Standardized Mean Differences"
            }
        }
        
        xlab <- if (abs) paste("Absolute", xlab.diff) else xlab.diff
        return(xlab)
    },
    love.plot_add_stars = function(SS.var, variable.names, ...) {
        A <- list(...)
        binary <- A$binary #attr(x, "print.options")$binary
        continuous <- A$continuous #attr(x, "print.options")$continuous
        var_type <- A$var_type #B[["Type"]]
        stars <- A$stars
        star_char = A$star_char #args$star_char
        
        #All std, no std, some std
        if (!((binary == "std" || sum(var_type == "Binary") == 0) && 
              (continuous == "std" || sum(var_type != "Binary") == 0)) 
            &&
            !((binary == "raw" || sum(var_type == "Binary") == 0) && 
              (continuous == "raw" || sum(var_type != "Binary") == 0))) {
            
            stars <- match_arg(stars, c("none", "std", "raw"))
            if (stars == "none") {
                warning("Standardized mean differences and raw mean differences are present in the same plot. \nUse the 'stars' argument to distinguish between them and appropriately label the x-axis.", call. = FALSE)
            }
            else {
                if (!rlang::is_string(star_char)) star_char <- "*"
                
                vars_to_star <- setNames(rep(FALSE, length(variable.names)), variable.names)
                if (stars == "std") {
                    if (binary == "std") vars_to_star[variable.names[var_type == "Binary"]] <- TRUE
                    if (continuous == "std") vars_to_star[variable.names[var_type != "Binary"]] <- TRUE
                }
                else if (stars == "raw") {
                    if (binary == "raw") vars_to_star[variable.names[var_type == "Binary"]] <- TRUE
                    if (continuous == "raw") vars_to_star[variable.names[var_type != "Binary"]] <- TRUE
                }
                new.variable.names <- setNames(variable.names, variable.names)
                names(new.variable.names)[vars_to_star[variable.names]] <- paste0(variable.names[vars_to_star[variable.names]], star_char)
                SS.var <- do.call(f.recode, c(list(SS.var), new.variable.names))
            }
        }
        
        return(SS.var)
    },
    baseline.xintercept = 0,
    threshold.xintercepts = function(threshold, abs) {
        if (abs) c(lower = base::abs(threshold))
        else c(lower = -base::abs(threshold), upper = base::abs(threshold))
    },
    love.plot_axis_scale = ggplot2::scale_x_continuous,
    fun = function(C, treat, weights, std, s.d.denom, abs, s.weights, bin.vars, weighted.weights = weights, subset = NULL, ...) {
        col_w_smd(C, treat = treat, weights = weights, 
                  std = std, s.d.denom = s.d.denom,
                  abs = abs, s.weights = s.weights, bin.vars = bin.vars,
                  weighted.weights = weighted.weights,
                  subset = subset)
    }
)}

STATS[["variance.ratios"]] <- {list(
    type = "bin",
    threshold = "v.threshold",
    Threshold = "V.Threshold",
    disp_stat = "disp.v.ratio",
    adj_only = FALSE,
    abs = function(x) abs_(x, ratio = TRUE),
    bal.tab_column_prefix = "V.Ratio", #Also which.stat in love.plot
    threshold_range = c(1, Inf),
    balance_tally_for = "variance ratios",
    variable_with_the_greatest = "variance ratio", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        "Variance Ratios"
    },
    love.plot_add_stars = function(SS.var, variable.names, ...) {
        return(SS.var)
    },
    baseline.xintercept = 1,
    threshold.xintercepts = function(threshold, abs) {
        if (abs) c(lower = abs_(threshold, ratio = TRUE))
        else c(lower = abs_(threshold, ratio = TRUE)^-1, upper = abs_(threshold, ratio = TRUE))
    },
    love.plot_axis_scale = ggplot2::scale_x_log10,
    fun = function(C, treat, weights, abs, s.weights, bin.vars, subset = NULL, ...) {
        vrs <- rep(NA_real_, ncol(C))
        if (any(!bin.vars)) {
            vrs[!bin.vars] <- col_w_vr(C[, !bin.vars, drop = FALSE], treat = treat, 
                                       weights = weights, abs = abs, 
                                       s.weights = s.weights, bin.vars = bin.vars[!bin.vars],
                                       subset = subset)
        }
        vrs
    }
)}

STATS[["ks.statistics"]] <- {list(
    type = "bin",
    threshold = "ks.threshold",
    Threshold = "KS.Threshold",
    disp_stat = "disp.ks",
    adj_only = FALSE,
    abs = function(x) abs_(x),
    bal.tab_column_prefix = "KS", #Also which.stat in love.plot
    threshold_range = c(0, 1),
    balance_tally_for = "KS statistics",
    variable_with_the_greatest = "KS statistic", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        "Kolmogorov-Smirnov Statistics"
    },
    love.plot_add_stars = function(SS.var, variable.names, ...) {
        return(SS.var)
    },
    baseline.xintercept = 0,
    threshold.xintercepts = function(threshold, abs) {
        c(lower = base::abs(threshold))
    },
    love.plot_axis_scale = ggplot2::scale_x_continuous,
    fun = function(C, treat, weights, s.weights, bin.vars, subset = NULL, ...) {
        A <- list(...)
        do.call("col_w_ks", c(list(C, treat = treat, weights = weights, s.weights = s.weights, bin.vars = bin.vars,
                                   subset = subset), A))
    }
)}

STATS[["ovl.coefficients"]] <- {list(
    type = "bin",
    threshold = "ovl.threshold",
    Threshold = "OVL.Threshold",
    disp_stat = "disp.ovl",
    adj_only = FALSE,
    abs = function(x) abs_(x),
    bal.tab_column_prefix = "OVL", #Also which.stat in love.plot
    threshold_range = c(0, 1),
    balance_tally_for = "overlapping coefficients",
    variable_with_the_greatest = "overlapping coefficient", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        "Overlapping Coefficients"
    },
    love.plot_add_stars = function(SS.var, variable.names, ...) {
        return(SS.var)
    },
    baseline.xintercept = 0,
    threshold.xintercepts = function(threshold, abs) {
        c(lower = base::abs(threshold))
    },
    love.plot_axis_scale = ggplot2::scale_x_continuous,
    fun = function(C, treat, weights, s.weights, bin.vars, integrate = FALSE, subset = NULL, ...) {
        A <- list(...)
        do.call("col_w_ovl", c(list(C, treat = treat, weights = weights, s.weights = s.weights, bin.vars = bin.vars,
                                    subset = subset, integrate = integrate), A))
    }
)}

STATS[["correlations"]] <- {list(
    type = "cont",
    threshold = "r.threshold",
    Threshold = "R.Threshold",
    disp_stat = "disp.corr",
    adj_only = FALSE,
    abs = function(x) abs_(x),
    bal.tab_column_prefix = "Corr", #Also which.stat in love.plot
    threshold_range = c(0, 1),
    balance_tally_for = "treatment correlations",
    variable_with_the_greatest = "treatment correlation", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        A <- list(...)
        abs <- A$abs
        if (abs) return("Absolute Treatment-Covariate Correlations")
        else return("Treatment-Covariate Correlations")
    },
    love.plot_add_stars = function(SS.var, variable.names, ...) {
        return(SS.var)
    },
    baseline.xintercept = 0,
    threshold.xintercepts = function(threshold, abs) {
        if (abs) c(lower = base::abs(threshold))
        else c(lower = -base::abs(threshold), upper = base::abs(threshold))
    },
    love.plot_axis_scale = ggplot2::scale_x_continuous,
    fun = function(C, treat, weights, abs, s.weights, std, s.d.denom, bin.vars, weighted.weights = weights, subset = NULL, ...) {
        col_w_cov(C, treat = treat, weights = weights, abs = abs, s.weights = s.weights, 
                  std = std, type = "pearson",
                  s.d.denom = s.d.denom,
                  bin.vars = bin.vars, weighted.weights = weighted.weights, na.rm = TRUE,
                  subset = subset)
    }
)}

STATS[["spearman.correlations"]] <- {list(
    type = "cont",
    threshold = "s.threshold",
    Threshold = "S.Threshold",
    disp_stat = "disp.spear",
    adj_only = FALSE,
    abs = function(x) abs_(x),
    bal.tab_column_prefix = "S.Corr", #Also which.stat in love.plot
    threshold_range = c(0, 1),
    balance_tally_for = "treatment Spearman correlations",
    variable_with_the_greatest = "treatment Spearman correlation", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        A <- list(...)
        abs <- A$abs
        if (abs) return("Absolute Treatment-Covariate Spearman Correlations")
        else return("Treatment-Covariate Spearman Correlations")
    },
    love.plot_add_stars = function(SS.var, variable.names, ...) {
        return(SS.var)
    },
    baseline.xintercept = 0,
    threshold.xintercepts = function(threshold, abs) {
        if (abs) c(lower = base::abs(threshold))
        else c(lower = -base::abs(threshold), upper = base::abs(threshold))
    },
    love.plot_axis_scale = ggplot2::scale_x_continuous,
    fun = function(C, treat, weights, abs, s.weights, std, s.d.denom, bin.vars, weighted.weights = weights, subset = NULL, ...) {
        col_w_cov(C, treat = treat, weights = weights, abs = abs, s.weights = s.weights, 
                  std = std, type = "spearman",
                  s.d.denom = s.d.denom,
                  bin.vars = bin.vars, weighted.weights = weighted.weights, na.rm = TRUE,
                  subset = subset)
    }
)}

STATS[["mean.diffs.target"]] <- {list(
    type = "cont",
    threshold = "m.threshold",
    Threshold = "M.Threshold",
    disp_stat = "disp.diff",
    adj_only = TRUE,
    abs = function(x) abs_(x),
    bal.tab_column_prefix = "Diff", #Also which.stat in love.plot
    threshold_range = c(0, Inf),
    balance_tally_for = "target mean differences",
    variable_with_the_greatest = "target mean difference", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        A <- list(...)
        binary <- A$binary #attr(x, "print.options")$binary
        continuous <- A$continuous #attr(x, "print.options")$continuous
        abs <- A$abs
        var_type <- A$var_type #B[["type"]]
        stars <- A$stars
        
        #All std, no std, some std
        if ((binary == "std" || !any(var_type == "Binary")) && 
            (continuous == "std" || !any(var_type == "Contin."))) {
            xlab.diff <- "Standardized Target Mean Differences"
        } 
        else if ((binary == "raw" || !any(var_type == "Binary")) && 
                 (continuous == "raw" || !any(var_type == "Contin."))) {
            xlab.diff <- "Target Mean Differences"
        }
        else {
            stars <- match_arg(stars, c("none", "std", "raw"))
            if (stars == "none") {
                xlab.diff <- "Target Mean Differences"
            }
            else if (stars == "std") {
                xlab.diff <- "Target Mean Differences"
            }
            else if (stars == "raw") {
                xlab.diff <- "Standardized Target Mean Differences"
            }
        }
        
        xlab <- if (abs) paste("Absolute", xlab.diff) else xlab.diff
        return(xlab)
    },
    love.plot_add_stars = function(SS.var, variable.names, ...) {
        A <- list(...)
        binary <- A$binary #attr(x, "print.options")$binary
        continuous <- A$continuous #attr(x, "print.options")$continuous
        var_type <- A$var_type #B[["Type"]]
        stars <- A$stars
        star_char = A$star_char #args$star_char
        
        #All std, no std, some std
        if (!((binary == "std" || sum(var_type == "Binary") == 0) && 
              (continuous == "std" || sum(var_type != "Binary") == 0)) 
            &&
            !((binary == "raw" || sum(var_type == "Binary") == 0) && 
              (continuous == "raw" || sum(var_type != "Binary") == 0))) {
            
            stars <- match_arg(stars, c("none", "std", "raw"))
            if (stars == "none") {
                warning("Standardized mean differences and raw mean differences are present in the same plot. \nUse the 'stars' argument to distinguish between them and appropriately label the x-axis.", call. = FALSE)
            }
            else {
                if (!rlang::is_string(star_char)) star_char <- "*"
                
                vars_to_star <- setNames(rep(FALSE, length(variable.names)), variable.names)
                if (stars == "std") {
                    if (binary == "std") vars_to_star[variable.names[var_type == "Binary"]] <- TRUE
                    if (continuous == "std") vars_to_star[variable.names[var_type != "Binary"]] <- TRUE
                }
                else if (stars == "raw") {
                    if (binary == "raw") vars_to_star[variable.names[var_type == "Binary"]] <- TRUE
                    if (continuous == "raw") vars_to_star[variable.names[var_type != "Binary"]] <- TRUE
                }
                new.variable.names <- setNames(variable.names, variable.names)
                names(new.variable.names)[vars_to_star[variable.names]] <- paste0(variable.names[vars_to_star[variable.names]], star_char)
                SS.var <- do.call(f.recode, c(list(SS.var), new.variable.names))
            }
        }
        
        return(SS.var)
    },
    baseline.xintercept = 0,
    threshold.xintercepts = function(threshold, abs) {
        if (abs) c(lower = base::abs(threshold))
        else c(lower = -base::abs(threshold), upper = base::abs(threshold))
    },
    love.plot_axis_scale = ggplot2::scale_x_continuous,
    fun = function(C, treat, weights, std, s.d.denom, abs, s.weights, bin.vars, weighted.weights = weights, subset = NULL, ...) {
        n <- nrow(C)
        C <- rbind(C, C)
        treat <- rep(c(0,1), each = n)
        if (is_not_null(weights)) weights <- c(weights, rep(1, n))
        if (is_not_null(s.weights)) s.weights <- c(s.weights, s.weights)
        if (is_not_null(subset)) subset <- c(subset, subset)
        s.d.denom <- "1"
        
        if (is_not_null(weights)) {
        col_w_smd(C, treat = treat, weights = weights, 
                  std = std, s.d.denom = s.d.denom,
                  abs = abs, s.weights = s.weights, bin.vars = bin.vars,
                  weighted.weights = weighted.weights,
                  subset = subset)
        }
        else rep(NA_real_, ncol(C))
    }
)}

STATS[["ks.statistics.target"]] <- {list(
    type = "cont",
    threshold = "ks.threshold",
    Threshold = "KS.Threshold",
    disp_stat = "disp.ks",
    adj_only = TRUE,
    abs = function(x) abs_(x),
    bal.tab_column_prefix = "KS", #Also which.stat in love.plot
    threshold_range = c(0, 1),
    balance_tally_for = "target KS statistics",
    variable_with_the_greatest = "target KS statistic", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        "Target Kolmogorov-Smirnov Statistics"
    },
    love.plot_add_stars = function(SS.var, variable.names, ...) {
        return(SS.var)
    },
    baseline.xintercept = 0,
    threshold.xintercepts = function(threshold, abs) {
        c(lower = base::abs(threshold))
    },
    love.plot_axis_scale = ggplot2::scale_x_continuous,
    fun = function(C, treat, weights, s.weights, bin.vars, subset = NULL, ...) {
        A <- list(...)
        
        n <- nrow(C)
        C <- rbind(C, C)
        treat <- rep(c(0,1), each = n)
        if (is_not_null(weights)) weights <- c(weights, rep(1, n))
        if (is_not_null(s.weights)) s.weights <- c(s.weights, s.weights)
        if (is_not_null(subset)) subset <- c(subset, subset)
        
        do.call("col_w_ks", c(list(C, treat = treat, weights = weights, s.weights = s.weights, bin.vars = bin.vars,
                                   subset = subset), A))
    }
)}

all_STATS <- function(type) {
    if (missing(type)) out <- names(STATS)
    else {
        type <- match_arg(type, c("bin", "cont"))
        out <- names(STATS)[get_from_STATS("type") == type]
    }
    unique(out)
}
