#STATS that rely on other packages

STATS[["distance.correlations"]] <- {list(
    type = "cont",
    threshold = "d.threshold",
    Threshold = "D.Threshold",
    disp_stat = "disp.dcorr",
    bin_only = FALSE,
    abs = function(x) abs_(x),
    bal.tab_column_prefix = "D.Corr", #Also which.stat in love.plot
    threshold_range = c(0, 1),
    balance_tally_for = "treatment distance correlations",
    variable_with_the_greatest = "treatment distance correlation", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        A <- list(...)
        abs <- A$abs
        if (abs) return("Treatment-Covariate Distance Correlations")
        else return("Treatment-Covariate Distance Correlations")
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
        if (length(treat) > 20000) stop("distance correlation cannot be used with sample sizes greater than 20000.", call. = FALSE)
        check.package("extracat")
        if (is_null(subset)) subset <- rep(TRUE, length(treat))
        if (is_null(weights)) weights <- rep(1, length(treat))
        apply(C[subset, , drop = FALSE], 2, extracat::wdcor, y = treat[subset], w = weights[subset]*s.weights[subset])
    }
)}