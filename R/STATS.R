#Trying to make a general stats object

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
    bin_only = FALSE,
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
                  subset = NULL)
    }
)}

STATS[["variance.ratios"]] <- {list(
    type = "bin",
    threshold = "v.threshold",
    Threshold = "V.Threshold",
    disp_stat = "disp.v.ratio",
    bin_only = TRUE,
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
                                       subset = NULL)
        }
        vrs
    }
)}

STATS[["ks.statistics"]] <- {list(
    type = "bin",
    threshold = "ks.threshold",
    Threshold = "KS.Threshold",
    disp_stat = "disp.ks",
    bin_only = FALSE,
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
                                   subset = NULL), A))
    }
)}

STATS[["ovl.coefficients"]] <- {list(
    type = "bin",
    threshold = "ovl.threshold",
    Threshold = "OVL.Threshold",
    disp_stat = "disp.ovl",
    bin_only = FALSE,
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
                                    subset = NULL, integrate = integrate), A))
    }
)}

STATS[["correlations"]] <- {list(
    type = "cont",
    threshold = "r.threshold",
    Threshold = "R.Threshold",
    disp_stat = "disp.corr",
    bin_only = FALSE,
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
                  subset = NULL)
    }
)}

STATS[["spearman.correlations"]] <- {list(
    type = "cont",
    threshold = "s.threshold",
    Threshold = "S.Threshold",
    disp_stat = "disp.spear",
    bin_only = FALSE,
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
                  subset = NULL)
    }
)}

all_STATS <- function(type) {
    if (missing(type)) names(STATS)
    else {
        type <- match_arg(type, c("bin", "cont"))
        names(STATS)[get_from_STATS("type") == type]
    }
}