#' @title Display Balance Statistics in a Love Plot
#' 
#' @description Generates a "Love" plot graphically displaying covariate balance before and after adjusting. Options are available for producing publication-ready plots. Detailed examples are available in `vignette("love.plot")`.
#' 
#' @param x the valid input to a call to [bal.tab()] (e.g., the output of a preprocessing function). Other arguments that would be supplied to `bal.tab()` can be entered with `...`. Can also be a `bal.tab` object, i.e., the output of a call to `bal.tab()`. See Examples. If `x` is not a `bal.tab` object, `love.plot()` calls `bal.tab()` with the arguments supplied.
#' @param stats `character`; which statistic(s) should be reported. See [`stats`][balance-statistics] for allowable options. For binary and multi-category treatments, "mean.diffs" (i.e., mean differences) is the default. For continuous treatments, "correlations" (i.e., treatment-covariate Pearson correlations) is the default. Multiple options are allowed.
#' @param abs `logical`; whether to present the statistic in absolute value or not. For variance ratios, this will force all ratios to be greater than or equal to 1. If `x` is a `bal.tab` object, `love.plot()` might ignore `abs` depending on the original `bal.tab()` call. If unspecified, uses whatever was used in the call to `bal.tab()`.
#' @param agg.fun if balance is to be displayed across clusters or imputations rather than within a single cluster or imputation, which summarizing function ("mean", "max", or "range") of the balance statistics should be used. If "range" is entered, `love.plot()` will display a line from the min to the max with a point at the mean for each covariate. Abbreviations allowed; "range" is default. Remember to set `which.<ARG> = .none` (where `<ARG>` is the grouping argument, such as `cluster` or `imp`) to use `agg.fun`. See Details.
#' @param var.order a `character` or `love.plot` object; how to order the variables in the plot. See Details. 
#' @param drop.missing `logical`; whether to drop rows for variables for which the statistic has a value of `NA`, for example, variance ratios for binary variables. If `FALSE`, there will be rows for these variables but no points representing their value. Default is `TRUE`, so that variables with missing balance statistics are absent. When multiple `stats` are requested, only variables with `NA`s for all `stats` will be dropped if `drop.missing = TRUE`. This argument used to be called `no.missing`, and that name still works (but has been deprecated).
#' @param drop.distance `logical`; whether to ignore the distance measure (if there are any) in plotting.
#' @param thresholds `numeric`; an optional value to be used as a threshold marker in the plot. Should be a named vector where each name corresponds to the statistic for which the threshold is to be applied. See example at [`stats`][balance-statistics]. If `x` is a `bal.tab` object and a threshold was set in it (e.g., with `thresholds`), its threshold will be used unless overridden using the `threshold` argument in `love.plot()`.
#' @param line `logical`; whether to display a line connecting the points for each sample.
#' @param stars when mean differences are to be displayed, which variable names should have a star (i.e., an asterisk) next to them. Allowable values are "none", "std" (for variables with mean differences that have been standardized), or "raw" (for variables with mean differences that have not been standardized). If "raw", the x-axis title will be "Standardized Mean Differences". Otherwise, it will be "Mean Differences". Ignored when mean difference are not displayed. See Details for an explanation of the purpose of this option.
#' @param grid `logical`; whether gridlines should be shown on the plot. Default is `FALSE`.
#' @param limits `numeric`; the bounds for the x-axis of the plot. Must a (named) list of vectors of length 2 in ascending order, one for each value of `stats` that is to have limits; e.g., `list(m = c(-.2, .2))`. If values exceed the limits, they will be plotted at the edge.
#' @param colors the colors of the points on the plot. See 'Color Specification' at [graphics::par()] or the `ggplot2` [aesthetic specifications](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#colour-and-fill) page. The first value corresponds to the color for the unadjusted sample, and the second color to the adjusted sample. If only one is specified, it will apply to both. Defaults to the default \pkg{ggplot2} colors.
#' @param shapes the shapes of the points on the plot. Must be one or two numbers between 1 and 25 or the name of a valid shape. See the `ggplot2` [aesthetic specifications](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#point) page for valid options. Values 15 to 25 are recommended. The first value corresponds to the shape for the unadjusted sample, and the second color to the adjusted sample. If only one is specified, it will apply to both. Defaults to 19 (`"circle filled"`).
#' @param alpha `numeric`; the transparency of the points. See [ggplot2::scale_alpha()].
#' @param size `numeric`; the size of the points on the plot. Defaults to 3. In previous versions, the size was scaled by a factor of 3. Now `size` corresponds directly to the `size` aesthetic in [ggplot2::geom_point()].
#' @param wrap `numeric`; the number of characters at which to wrap axis labels to the next line. Defaults to 30. Decrease this if the axis labels are excessively long.
#' @param var.names an optional object providing alternate names for the variables in the plot, which will otherwise be the variable names as they are stored. This may be useful when variables have ugly names. See Details on how to specify `var.names`. [var.names()] can be a useful tool for extracting and editing the names from the `bal.tab` object.
#' @param title `character`; the title of the plot.
#' @param sample.names `character`; new names to be given to the samples (i.e., in place of "Unadjusted" and "Adjusted"). For example, when matching it used, it may be useful to enter `c("Unmatched", "Matched")`.
#' @param labels `logical` or `character`; labels to give the plots when multiple `stats` are requested. If `TRUE`, the labels will be capital letters. Otherwise, must be a string with the same length as `stats`. This can be useful when the plots are to be used in an article.
#' @param position the position of the legend. When `stats` has length 1, this can be any value that would be appropriate as an argument to `legend.position` in [ggplot2::theme()]. When `stat` has length greater than 1, can be one of "none", "left", "right", "bottom", or "top".
#' @param themes an optional list of `theme` objects to append to each individual plot. Each entry should be the output of a call to [ggplot2::theme()] in \pkg{ggplot2}. This is a way to customize the individual plots when multiple `stats` are requested since the final output is not a manipulable `ggplot` object. It can be used with length-1 `stats`, but it probably makes more sense to just add the `theme()` call after `love.plot()`.
#' @param ... additional arguments passed to `bal.tab()` or options for display of the plot. The following related arguments are currently accepted:
#' \describe{
#'     \item{`use.grid`}{whether to use [gridExtra::arrangeGrob()] in `gridExtra` to make the plot when `stats` has length 1. See section Value.}
#'     \item{`disp.subclass`}{whether to display individual subclasses if subclassification is used. Overrides the `disp.subclass` option in the original `bal.tab()` call if `x` is a `bal.tab` object.}
#'     \item{`star_char`}{`character`; when `stars` are used, the character that should be the "star" next to the starred variables. The default is `"*"`. `"†"` or `"\u2020"` (i.e., dagger) might be appealing as well.}
#' }
#' Additionally, any of the `which.` arguments used with clustered or multiply imputed data or longitudinal or multi-category treatments can be specified to display balance on selected groupings. Set to `.none` to aggregate across groups (in which `agg.fun` comes into effect) and set to `.all` to view all groups. See [display-options] for options, and see `vignette("segmented-data")` for details and examples.
#' 
#' @returns When only one type of balance statistic is requested, the returned object is a standard `ggplot` object that can be manipulated using \pkg{ggplot2} syntax. This facilitates changing fonts, background colors, and features of the legend outside of what `love.plot()` provides automatically. 
#' 
#' When more than one type of balance statistic is requested, the plot is constructed using [gridExtra::arrangeGrob()] in `gridExtra`, which arranges multiple plots and their shared legend into one plot. Because the output of `arrangeGrob` is a `gtable` object, its features cannot be manipulated in the standard way. Use the `themes` argument to change theme elements of the component plots. The original plots are stored in the `"plots"` attribute of the output object.
#' 
#' @details `love.plot` can be used with clusters, imputations, and multi-category and longitudinal treatments in addition to the standard case. Setting the corresponding `which.` argument to `.none` will aggregate across that dimension. When aggregating, an argument should be specified to `agg.fun` referring to whether the mean, minimum ("min"), or maximum ("max") balance statistic or range ("range", the default) of balance statistics for each covariate should be presented in the plot. See `vignette("segmented-data")` for examples.
#' 
#' With subclasses, balance will be displayed for the unadjusted sample and the aggregated subclassified sample. If `disp.subclass` is `TRUE`, each subclass will be displayed additionally as a number on the plot. 
#' 
#' ### Variable order using `var.order`
#' 
#' The order that the variables are presented in depends on the argument to `var.order`. If `NULL`, the default, they will be displayed in the same order as in the call to `bal.tab()`, which is the order of the underlying data set. If "alphabetical", they will be displayed in alphabetical order. If "unadjusted", they will be ordered by the balance statistic of the unadjusted sample. To order by the values of the adjusted sample, "adjusted" can be supplied if only one set of weights (or subclasses) are specified; otherwise, the name of the set of weights should be specified.
#' 
#' If multiple `stats` are requested, the order will be determined by the first entry to `stats` (e.g., if both "mean.diffs" and "ks.statistics" are requested, and `var.order = "unadjusted"`, the variables will be displayed in order of the unadjusted mean differences for both plots). If multiple plots are produced simultaneously (i.e., for individual clusters or imputations), `var.order` can only be `NULL` or "alphabetical".
#' 
#' If a `love.plot` object is supplied, the plot being drawn will use the variable order in the supplied `love.plot` object. This can be useful when making more than one plot and the variable order should be the same across plots.
#' 
#' ### Variable names using `var.names`
#' 
#' The default in `love.plot()` is to present variables as they are named in the output of the call to `bal.tab()`, so it is important to know this output before specifying alternate variable names when using `var.names`, as the displayed variable names may differ from those in the original data.
#' 
#' There are several ways to specify alternate names for presentation in the displayed plot using the `var.names` argument by specifying a list of old and new variable names, pairing the old name with the new name. You can do this in three ways: 1) use a vector or list of new variable names, with the `names` of the values the old variable names; 2) use a data frame with exactly one column containing the new variable names and the row names containing the old variable names; or 3) use a data frame with two columns, the first (or the one named "old") containing the old variable names and the second (or the one named "new") containing the new variable names. If a variable in the output from `bal.tab()` is not provided in the list of old variable names, `love.plot()` will use the original old variable name.
#' 
#' `love.plot()` can replace old variables names with new ones based on exact matching for the name strings or matching using the variable name components. For example, if a factor variable `"X"` with levels `"a"`, `"b"`, and `"c"` is displayed with `love.plot()`, the variables `"X_a"`, `"X_b"`, and `"X_c"` will be displayed. You can enter replacement names for all three variables individually with `var.names`, or you can simply specify a replacement name for `"X"`, and `"X"` will be replaced by the given name in all instances it appears, including not just factor expansions, but also polynomials and interactions in `int = TRUE` in the original `bal.tab()` call. In an interaction with another variable, say `"Y"`, there are several ways to replace the name of the interaction term `"X_a * Y"`. If the entire string (`"X_a * Y"`) is included in `var.names`, the entire string will be replaced. If `"X_a"` is included in `var.names`, only it will be replaced (and it will be replaced everywhere else it appears). If `"X"` is included in `var.names`, only it will be replaced (and it will be replaced everywhere else it appears). See example at [var.names()].
#' 
#' ### Stars and the x-axis label with mean differences
#' 
#' When mean differences are to be displayed, `love.plot()` attempts to figure out the appropriate label for the x-axis. If all mean differences are standardized, the x-axis label will be "Standardized Mean Differences". If all mean differences are raw (i.e., unstandardized), the x-axis label will be "Mean Differences". Otherwise, `love.plot()` turns to the `stars` argument. If "raw", the x-axis label will be "Standardized Mean Differences" (i.e., because un-starred variables have standardized mean differences displayed). If "std", the x-axis label will be "Mean Differences" (i.e., because un-starred variables have raw mean differences displayed). If "none", the x-axis label will be "Mean Differences" and a warning will be issued recommending the use of `stars`. 
#' 
#' The default is to display standardized mean differences for continuous variables, raw mean differences for binary variables, and no stars, so this warning will be issued in most default uses of `love.plot()`. The purpose of this is to correct behavior of previous versions of \pkg{cobalt} in which the default x-axis label was "Mean Differences", even when standardized mean differences were displayed, yielding a potentially misleading plot. This warning requires the user to think about what values are being displayed. The idea of using `stars` is that the user can, in a caption for the plot, explain that variables with an asterisk have standardized (or raw) mean differences display, in contrast to un-starred variables.
#' 
#' @note `love.plot` can also be called by using `plot()` or `autoplot()` on a `bal.tab` object. If used in this way, some messages may appear twice. It is recommended that you just use `love.plot()` instead.
#' 
#' @seealso 
#' [bal.tab()], `vignette("love.plot")`
#' 
#' @examplesIf requireNamespace("WeightIt", quietly = TRUE)
#' data("lalonde", package = "cobalt")
#' 
#' ## Propensity score weighting
#' w.out1 <- WeightIt::weightit(treat ~ age + educ + race + married +
#'                                  nodegree + re74 + re75, 
#'                              data = lalonde)
#' 
#' love.plot(w.out1, thresholds = c(m = .1), var.order = "unadjusted")
#' 
#' ## Using alternate variable names
#' v <- data.frame(old = c("age", "educ", "race_black", "race_hispan", 
#'                         "race_white", "married", "nodegree", "re74", 
#'                         "re75", "distance"),
#'                 new = c("Age", "Years of Education", "Black", 
#'                         "Hispanic", "White", "Married", "No Degree", 
#'                         "Earnings 1974", "Earnings 1975", 
#'                         "Propensity Score"))
#' 
#' love.plot(w.out1, stats = "m", threshold = .1, 
#'           var.order = "unadjusted", var.names = v)
#' 
#' #Using multiple stats
#' love.plot(w.out1, stats = c("m", "ks"), 
#'           thresholds = c(m = .1, ks = .05), 
#'           var.order = "unadjusted", var.names = v, stars = "raw",
#'           position = "bottom", wrap = 20)
#' 
#' #Changing visual elements
#' love.plot(w.out1, thresholds = c(m = .1), 
#'           var.order = "unadjusted", var.names = v, abs = TRUE,
#'           shapes = c("triangle filled", "circle"), 
#'           colors = c("red", "blue"), line = TRUE,
#'           grid = FALSE, sample.names = c("Original", "Weighted"),
#'           stars = "raw", position = "top")
#' 

#' @rdname love.plot
#' @export 
love.plot <- function(x, stats, abs, agg.fun = NULL, 
                      var.order = NULL, drop.missing = TRUE, drop.distance = FALSE, 
                      thresholds = NULL, line = FALSE, stars = "none", grid = FALSE, 
                      limits = NULL, colors = NULL, shapes = NULL, alpha = 1, size = 3, 
                      wrap = 30, var.names = NULL, title, sample.names, labels = FALSE,
                      position = "right", themes = NULL, ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call()
    .alls <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.all)), logical(1L))
    .nones <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.none)), logical(1L))
    if (any(c(.alls, .nones))) {
        .call[.alls] <- expression(NULL)
        .call[.nones] <- expression(NA)
        return(eval.parent(.call))
    }
    
    if (missing(stats)) stats <- NULL
    
    #Re-call bal.tab with disp.v.ratio or disp.ks if stats = "v" or "k".
    if (typeof(.call[["x"]]) == "language") { #if x is not an object (i.e., is a function call)
        
        replace.args <- function(m) {
            #m is bal.tab call or list (for do.call)
            m[["un"]] <- TRUE
            m[["subclass.summary"]] <- TRUE
            
            if (is_not_null(stats)) m[["stats"]] <- stats
            
            if (any(names(m) == "agg.fun")) m[["agg.fun"]] <- NULL
            
            if (any(names(m) %pin% "abs")) m[["abs"]] <- abs
            
            if (any(names(m) %pin% "thresholds")) m["thresholds"] <- list(NULL)
            
            m
        }
        
        if (deparse1(.call[["x"]][[1]]) %in% c("bal.tab", "cobalt::bal.tab", methods("bal.tab"))) { #if x i bal.tab call
            .call[["x"]] <- replace.args(.call[["x"]])
            x <- eval.parent(.call[["x"]])
            
        }
        else if (deparse1(.call[["x"]][[1]]) == "do.call") { #if x is do.call
            d <- match.call(eval(.call[["x"]][[1]]), .call[["x"]])
            if (deparse1(d[["what"]]) %in% c("bal.tab", "cobalt::bal.tab", methods("bal.tab"))) {
                d[["args"]] <- replace.args(d[["args"]])
                x <- eval.parent(d)
            }
        }
    }
    
    tryCatch(force(x), error = function(e) .err(conditionMessage(e)))
    
    if (!inherits(x, "bal.tab")) {
        #Use bal.tab on inputs first, then love.plot on that
        .call2 <- .call
        .call2[[1]] <- quote(cobalt::bal.tab)
        .call2[["x"]] <- x
        
        .call2["thresholds"] <- list(NULL)
        
        .call[["x"]] <- .call2
        
        return(eval.parent(.call))
    }
    
    args <- list(...)
    
    #shape (deprecated)
    #un.color (deprecated)
    #adj.color (deprecated)
    #cluster.fun (deprecated)
    #star_char
    
    p.ops <- c("which.cluster", "which.imp", "which.treat", "which.time", "disp.subclass")
    for (i in p.ops) {
        if (rlang::has_name(args, i)) attr(x, "print.options")[[i]] <- args[[i]]
    }
    
    #Using old argument names
    if (is_not_null(args$cluster.fun) && is_null(agg.fun)) agg.fun <- args$cluster.fun
    if (is_not_null(args$no.missing)) drop.missing <- args$no.missing
    
    Agg.Fun <- NULL
    subtitle <- NULL
    
    #Process abs
    if (missing(abs)) {
        abs <- if_null_then(attr(x, "print.options")[["abs"]], TRUE)
    }
    
    #Process stats
    if (is_null(stats)) stats <- attr(x, "print.options")$stats
    stats <- match_arg(stats, all_STATS(attr(x, "print.options")$type), several.ok = TRUE)
    
    #Get B and config
    if (inherits(x, "bal.tab.subclass")) {
        if (is_null(x[["Balance.Across.Subclass"]])) {
            .err("`subclass.summary` must be set to `TRUE` in the original call to `bal.tab()`")
        }
        B <- cbind(x[["Balance.Across.Subclass"]], variable.names = row.names(x[["Balance.Across.Subclass"]]))
        
        disp.subclass <- isTRUE(attr(x, "print.options")$disp.subclass)
        if (disp.subclass) {
            subclass.names <- names(x[["Subclass.Balance"]])
            sub.B <- do.call("cbind", c(
                lapply(subclass.names, function(s) {
                    sub <- x[["Subclass.Balance"]][[s]]
                    setNames(sub[endsWith(names(sub), ".Adj")],
                                       gsub(".Adj", paste0(".", s), names(sub)[endsWith(names(sub), ".Adj")]))
                }),
                list(variable.names = row.names(x[["Balance.Across.Subclass"]]))))
        }
        else {
            subclass.names <- sub.B <- NULL
        }
        
        attr(x, "print.options")$weight.names <- "Adj"
        subtitle <- "Across Subclasses"
        config <- "agg.none"
        facet <- NULL
    }
    else {
        B_list <- unpack_bal.tab(x)
        namesep <- attr(B_list, "namesep")
        class_sequence <- attr(B_list, "class_sequence")
        
        if (is_not_null(class_sequence)) {
            #Multiple layers present
            facet_mat <- as.matrix(do.call(rbind, strsplit(names(B_list), namesep, fixed = TRUE)))
            facet <- unname(vapply(class_sequence, switch, character(1L),
                                   bal.tab.cluster = "cluster",
                                   bal.tab.msm = "time",
                                   bal.tab.multi = "treat",
                                   bal.tab.imp = "imp", NULL))
            dimnames(facet_mat) <- list(names(B_list), facet)
            
            for (b in seq_along(B_list)) {
                B_list[[b]][["variable.names"]] <- factor(rownames(B_list[[b]]), levels = rownames(B_list[[b]]))
                for (i in facet) {
                    B_list[[b]][[i]] <- {
                        if (i == "imp") factor(paste("Imputation:", facet_mat[b, i]),
                                               levels = paste("Imputation:", sort(unique(as.numeric(facet_mat[b, i])))))
                        else facet_mat[b, i]
                    }
                }
            }
            
            #Process which. so that B_list can be shortened
            agg.over <- character(0)
            for (i in facet) {
                which. <- attr(x, "print.options")[[paste0("which.", i)]]
                if (is_null(which.)) {
                    #All levels; facet_mat stays the same.
                }
                else if (anyNA(which.)) {
                    agg.over <- c(agg.over, i)
                }
                else {
                    if (i == "treat") {
                        treat_levels <- attr(x, "print.options")$treat_vals_multi
                        if (is.numeric(which.)) which. <- treat_levels[which.]
                        if (!all(which. %in% treat_levels)) {
                            .err("all values in `which.treat` must be names or indices of treatment levels")
                        }
                        if (attr(x, "print.options")$pairwise) {
                            vs.combs <- cbind(vs.tmp <- as.matrix(expand.grid(treat_levels, treat_levels, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)), 
                                              apply(vs.tmp, 1, paste, collapse = " vs. "))
                            vs.combs <- vs.combs[vs.combs[,3] %in% facet_mat[, i],]
                            if (length(which.) == 1) facet_mat <- facet_mat[facet_mat[,i] %in% vs.combs[,3][vs.combs[,1] == which. | vs.combs[,2] == which.], , drop = FALSE]
                            else facet_mat <- facet_mat[facet_mat[,i] %in% vs.combs[,3][vs.combs[,1] %in% which. & vs.combs[,2] %in% which.], , drop = FALSE]
                        }
                        else {
                            vs.combs <- cbind(vs.tmp <- as.matrix(data.frame("Others", treat_levels, stringsAsFactors = FALSE)), 
                                              apply(vs.tmp, 1, paste, collapse = " vs. "))
                            vs.combs <- vs.combs[vs.combs[,3] %in% facet_mat[, i],]
                            facet_mat <- facet_mat[facet_mat[,i] %in% vs.combs[,3][vs.combs[,2] %in% which.], , drop = FALSE]
                        }
                    }
                    else {
                        if (is.numeric(which.) && max(which.) <= nunique(facet_mat[,i])) {
                            if (i == "imp") facet_mat <- facet_mat[facet_mat[,i] %in% as.character(which.), ,drop = FALSE]
                            facet_mat <- facet_mat[facet_mat[,i] %in% sort(unique(facet_mat[,i]))[which.], ,drop = FALSE]
                        }
                        else if (is.character(which.) && all(which. %in% unique(facet_mat[,i]))) {
                            facet_mat <- facet_mat[facet_mat[,i] %in% which., ,drop = FALSE]
                        }
                        else .err(sprintf("The argument to `which.%s` must be `.none`, `.all`, or the desired levels or indices of %s", i, switch(i, time = "time points", i)))
                    }
                }
                
            }
            B_list <- B_list[rownames(facet_mat)]
            B_names <- names(B_list[[1]])
            
            stat.cols <- expand.grid_string(vapply(stats, function(s) STATS[[s]]$bal.tab_column_prefix, character(1L)),
                                            c("Un", attr(x, "print.options")[["weight.names"]]),
                                            collapse = ".")
            stat.cols <- stat.cols[stat.cols %in% B_names]
            cols.to.keep <- c("variable.names", "Type", facet, stat.cols)
            
            for (b in seq_along(B_list)) {
                B_list[[b]] <- B_list[[b]][cols.to.keep]
            }
            
            B_stack <- do.call("rbind", c(B_list, list(make.row.names = FALSE)))
            
            if (is_not_null(agg.over)) {
                if (is_null(agg.fun)) {
                    agg.fun <- {
                        if (any(c("treat", "time") %in% agg.over)) "max"
                        else "range"
                    }
                }
                agg.fun <- tolower(agg.fun)
                Agg.Fun <- firstup(agg.fun <- match_arg(agg.fun, c("range", "max", "mean")))
                if (agg.fun == "max") abs <- TRUE
                
                if (abs) {
                    B_stack[stat.cols] <- lapply(stat.cols, function(sc) {
                        abs_(B_stack[[sc]], ratio = startsWith(sc, "V.Ratio"))
                    })
                }
                
                facet <- setdiff(facet, agg.over)
                
                aggregate_B <- function(FUN, B) {
                    B_agged <- aggregate(B[stat.cols], 
                                         by = B[c("variable.names", "Type", facet)], 
                                         FUN = FUN)
                    names(B_agged)[names(B_agged) %in% stat.cols] <- paste.(firstup(FUN), names(B_agged)[names(B_agged) %in% stat.cols])
                    B_agged
                }
                
                if (agg.fun == "range") {
                    B <- Reduce(function(x, y) merge(x, y, by = c("variable.names", "Type", facet), 
                                                     sort = FALSE),
                                lapply(c("min", "mean", "max"), aggregate_B, B_stack))
                }
                else {
                    B <- aggregate_B(agg.fun, B_stack)
                }
                
                B <- B[order(B[["variable.names"]]),]
                
                subtitle1 <- paste0(Agg.Fun, " across ", word_list(vapply(agg.over, switch, character(1L),
                                                                          "cluster" = "clusters",
                                                                          "time" = "time points",
                                                                          "treat" = "treatment pairs",
                                                                          "imp" = "imputations")))
                config <- paste.("agg", agg.over)
            }
            else {
                B <- B_stack
                subtitle1 <- NULL
                config <- "agg.none"
            }
            
            one.level.facet <- facet[vapply(B[facet], all_the_same, logical(1L))]
            subtitle2 <- {
                if (is_null(one.level.facet)) NULL
                else paste(vapply(one.level.facet, function(olf) {
                    paste(firstup(olf), B[1, olf], sep = ": ")
                }, character(1L)), collapse = ", ")
            }
            
            B[names(B) %in% one.level.facet] <- NULL
            
            if (sum(facet %nin% one.level.facet) > 1) {
                .err(sprintf("At least one of %s must be `.none` or of length 1",
                             word_list(paste.("which", facet), "or", quotes = "`")))
            }
            
            facet <- setdiff(facet, one.level.facet)
            
            subtitle <- paste(c(subtitle1, subtitle2), collapse = "\n")
            
            #one.level.facet - go in subtitle
            #facet - go in facet
            #agg.over - aggregated (e.g., averaged) over
            
        }
        else {
            #Single-layer bal.tab
            B <- cbind(B_list, 
                       variable.names = factor(rownames(B_list), levels = rownames(B_list)))
            
            facet <- one.level.facet <- agg.over <- NULL
            
            B_names <- names(B)
            
            stat.cols <- expand.grid_string(vapply(stats, function(s) STATS[[s]]$bal.tab_column_prefix, character(1L)),
                                            c("Un", attr(x, "print.options")[["weight.names"]]),
                                            collapse = ".")
            stat.cols <- stat.cols[stat.cols %in% B_names]
            
            cols.to.keep <- c("variable.names", "Type", stat.cols)
            B <- B[cols.to.keep]
            
            config <- "agg.none"
            
            subtitle <- NULL
        }
        sub.B <- NULL
        disp.subclass <- NULL
    }
    
    if (is_not_null(facet) && length(stats) > 1) {
        .err("`stats` can only have a length of 1 when faceting by other dimension (e.g., cluster, treatment)")
    }
    
    if (is_not_null(agg.fun) && config == "agg.none") {
        .wrn("no aggregation will take place, so `agg.fun` will be ignored. Remember to set `which.<ARG> = .none` to aggregate across <ARG>")
    }
    
    #Process variable names
    if (is_not_null(var.names)) {
        if (is.data.frame(var.names)) {
            if (ncol(var.names)==1) {
                if (is_not_null(row.names(var.names))) {
                    new.labels <- setNames(unlist(as.character(var.names[,1])), rownames(var.names))
                }
                else .wrn("`var.names` is a data frame, but its rows are unnamed")
            }
            else {
                if (all(c("old", "new") %in% names(var.names))) {
                    new.labels <- setNames(unlist(as.character(var.names[,"new"])), var.names[,"old"])
                }
                else {
                    if (ncol(var.names)>2) .wrn("only using first 2 columns of `var.names`")
                    new.labels <- setNames(unlist(as.character(var.names[,2])), var.names[,1])
                }
            } 
        }
        else if (is.atomic(var.names)) {
            if (is_not_null(names(var.names))) {
                new.labels <- setNames(as.character(var.names), names(var.names))
            }
            else .wrn("`var.names` is a vector, but its values are unnamed")
        }
        else if (is.list(var.names)) {
            if (!all(vapply(var.names, chk::vld_character_or_factor, logical(1L)))) {
                .wrn("`var.names` is a list, but its values are not the new names of the variables")
            }
            else if (is_null(names(var.names))) {
                .wrn("`var.names` is a list, but its values are unnamed")
            }
            else{
                new.labels <- unlist(var.names) #already a list
            }
        }
        else {
            .wrn("the argument to `var.names` is not one of the accepted structures and will be ignored.\n  See `?love.plot` for details")
        }
        
        co.names <- attr(x, "print.options")[["co.names"]]
        seps <- attr(co.names, "seps")
        for (i in names(co.names)) {
            comp <- co.names[[i]][["component"]]
            type <- co.names[[i]][["type"]]
            
            if (i %in% names(new.labels) && !is.na(new.labels[i])) {
                co.names[[i]][["component"]] <- new.labels[i]
                co.names[[i]][["type"]] <- "base"
            }
            else {
                if ("isep" %in% type) {
                    named.vars <- character(sum(type == "isep") + 1)
                    sep.inds <- c(which(type == "isep"), length(comp) + 1)
                    named.vars <- lapply(seq_along(sep.inds), function(k) {
                        inds <- (if (k == 1) seq(1, sep.inds[k] - 1) 
                                 else seq(sep.inds[k-1] + 1, sep.inds[k] - 1))
                        var <- comp[inds]
                        var.is.base <- type[inds] == "base"
                        pasted.var <- paste(var, collapse = "")
                        if (pasted.var %in% names(new.labels)) return(new.labels[pasted.var])
                       
                        paste(ifelse(var.is.base & var %in% names(new.labels) & !is.na(new.labels[var]),
                                     new.labels[var], var), collapse = "")
                    })
                    co.names[[i]][["component"]] <- do.call("paste", c(unname(named.vars), list(sep = seps["int"])))
                }
                else co.names[[i]][["component"]] <- ifelse(type == "base" & comp %in% names(new.labels) & !is.na(new.labels[comp]), new.labels[comp], comp)
            }
        }
        
        recode.labels <- setNames(names(co.names), 
                                  vapply(co.names, function(x) paste0(x[["component"]], collapse = ""), character(1L)))
        
        B[["variable.names"]] <- do.call(f.recode, c(list(B[["variable.names"]]), recode.labels))
    }
    
    distance.names <- as.character(unique(B[["variable.names"]][B[["Type"]] == "Distance"], nmax = sum(B[["Type"]] == "Distance")))
    
    if (drop.distance) {
        B <- B[B[["variable.names"]] %nin% distance.names, , drop = FALSE]
    }
    
    #Process variable order
    if (is_not_null(var.order) && !inherits(var.order, "love.plot")) {
        if (!inherits(x, "bal.tab.subclass") && 
            (is_null(attr(x, "print.options")$nweights) ||
             attr(x, "print.options")$nweights == 0)) {
            ua <- c("Unadjusted", "Alphabetical")
            names(ua) <- c("unadjusted", "alphabetical")
        }
        else if (inherits(x, "bal.tab.subclass") ||
                 attr(x, "print.options")$nweights == 1) {
            ua <- c("Adjusted", "Unadjusted", "Alphabetical")
            names(ua) <- c("adjusted", "unadjusted", "alphabetical")
        }
        else {
            ua <- c("Unadjusted", attr(x, "print.options")$weight.names, "Alphabetical")
            names(ua) <- c("unadjusted", attr(x, "print.options")$weight.names, "alphabetical")
        }
        if (get_from_STATS("adj_only")[stats[1]]) ua <- ua[names(ua) != "unadjusted"]
        var.order <- ua[match_arg(var.order, tolower(ua))]
    }
    
    #Process sample names
    
    ntypes <- length(attr(x, "print.options")$weight.names) + 1
    
    original.sample.names <- c("Unadjusted", attr(x, "print.options")$weight.names)
    if (length(original.sample.names) == 2) original.sample.names[2] <- "Adjusted"
    
    if (!missing(sample.names)) {
        if (!is.character(sample.names)) {
            .wrn("the argument to `sample.names` must be a character vector. Ignoring `sample.names`")
            sample.names <- NULL
        }
        else if (length(sample.names) %nin% c(ntypes, ntypes - 1)) {
            .wrn("the argument to `sample.names` must contain as many names as there are sample types, or one fewer. Ignoring `sample.names`")
            sample.names <- NULL
        }
    }
    else sample.names <- NULL
    
    if (is_null(sample.names)) {
        sample.names <- original.sample.names
    }
    else if (length(sample.names) == ntypes - 1) {
        sample.names <- c("Unadjusted", sample.names)
    }
    names(sample.names) <- original.sample.names
    
    #Process limits
    if (is_not_null(limits)) {
        if (!is.list(limits)) {
            limits <- list(limits)
        }
        if (any(vapply(limits, 
                       function(l) !is.numeric(l) || length(l) %nin% c(0L, 2L), 
                       logical(1L)))) {
            .wrn("`limits` must be a list of numeric vectors of legnth 2. Ignoring `limits`")
            limits <- NULL
        }
        
        if (is_not_null(names(limits))) {
            names(limits) <- stats[pmatch(names(limits), stats, duplicates.ok = TRUE)]
            limits <- limits[!is.na(names(limits))]
        }
        else {
            names(limits) <- stats[seq_along(limits)]
        }
    }
    
    #Setting up appearance
    
    #Alpha (transparency)
    if (is.numeric(alpha[1]) && !anyNA(alpha[1]) && 
        between(alpha[1], c(0,1))) {
        alpha <- alpha[1]
    }
    else {
        .wrn("the argument to `alpha` must be a number between 0 and 1. Using 1 instead")
        alpha <- 1
    }
    
    #Color
    if (is_not_null(args[["colours"]])) colors <- args[["colours"]]
    
    if (is_null(colors)) {
        colors <- {
            if (shapes.ok(shapes, ntypes) && length(shapes) > 1 && length(shapes) == ntypes) {
                rep("black", ntypes)
            }
            else gg_color_hue(ntypes)
        }
    }
    else {
        if (length(colors) == 1) {
            colors <- rep(colors, ntypes)
        }
        else if (length(colors) > ntypes) {
            colors <- colors[seq_len(ntypes)]
            .wrn(sprintf("only using first %s value%%s in `colors`", ntypes),
                 n = ntypes)
        }
        else if (length(colors) < ntypes) {
            .wrn("not enough colors were specified. Using default colors instead")
            colors <- gg_color_hue(ntypes)
        }
        
        if (!all(vapply(colors, isColor, logical(1L)))) {
            .wrn("the argument to `colors` contains at least one value that is not a recognized color. Using default colors instead")
            colors <- gg_color_hue(ntypes)
        }
        
    }
    # colors[] <- vapply(colors, col_plus_alpha, character(1L), alpha = alpha)
    names(colors) <- sample.names
    fill <- colors
    
    #Shapes
    if (is_null(shapes)) {
        shapes <- assign.shapes(colors)
    }
    else if (!shapes.ok(shapes, ntypes)) {
        .wrn(sprintf("the argument to `shape` must be %s valid shape%%s. See `?love.plot` for more information.\nUsing default shapes instead", ntypes), n = ntypes)
        shapes <- assign.shapes(colors)
    }
    else if (length(shapes) == 1) {
        shapes <- rep(shapes, ntypes)
    }
    names(shapes) <- sample.names
    
    #Size
    if (is.numeric(size)) size <- size[1]
    else {
        .wrn("the argument to `size` must be a number. Using 3 instead")
        size <- 3
    }
    
    stroke <- rep(0, ntypes)
    size0 <- size <- rep(size, ntypes)
    
    shapes.with.fill <- grepl("filled", shapes, fixed = TRUE)
    stroke[shapes.with.fill] <- size[shapes.with.fill]/3
    size[shapes.with.fill] <- size[shapes.with.fill]* .58
    
    # stroke <- .8*size
    
    if (is_not_null(facet) && is_not_null(var.order) &&
        !inherits(var.order, "love.plot") && tolower(var.order) != "alphabetical") {
            .wrn("`var.order` cannot be set with faceted plots (unless \"alphabetical\"). Ignoring `var.order`")
            var.order <- NULL
    }
    
    agg.range <- isTRUE(Agg.Fun == "Range")
    
    #Process thresholds
    thresholds <- if_null_then(attr(x, "print.options")$thresholds[stats], 
                               process_thresholds(thresholds, stats))
    
    #Title
    if (missing(title)) title <- "Covariate Balance"
    else title <- as.character(title)
    # if (missing(subtitle)) subtitle <- as.character(subtitle)
    
    #Process themes
    if (is_not_null(themes)) {
        if (!is.vector(themes, "list")) {
            themes <- list(themes)
        }
        
        if (any(vapply(themes, 
                       function(t) !inherits(t, "theme") || !inherits(t, "gg"), 
                       logical(1L)))) {
            .wrn("`themes` must be a list of `theme` objects. Ignoring `themes`")
            themes <- NULL
        }
        
        if (is_not_null(names(themes))) {
            names(themes) <- stats[pmatch(names(themes), stats, duplicates.ok = TRUE)]
            themes <- themes[!is.na(names(themes))]
        }
        else {
            names(themes) <- stats[1:length(themes)]
        }
    }
    
    variable.names <- as.character(B[["variable.names"]])
    
    plot.list <- make_list(stats)
    for (s in stats) {
        adj_only <- get_from_STATS("adj_only")[s]
        col.sample.names <- c("Un"[!adj_only], attr(x, "print.options")$weight.names)
        
        #Get SS
        if (agg.range) {
            SS <- do.call("rbind", 
                          lapply(col.sample.names,
                                 function(w) data.frame(var = variable.names,
                                                        type = B[["Type"]],
                                                        min.stat = B[[paste.("Min", STATS[[s]]$bal.tab_column_prefix, w)]],
                                                        max.stat = B[[paste.("Max", STATS[[s]]$bal.tab_column_prefix, w)]],
                                                        mean.stat = B[[paste.("Mean", STATS[[s]]$bal.tab_column_prefix, w)]],
                                                        Sample = switch(w, "Un"= "Unadjusted", 
                                                                        "Adj" = "Adjusted", w),
                                                        B[facet],
                                                        row.names = NULL,
                                                        stringsAsFactors = TRUE)))
            
            sample.vals <- sample.names[levels(SS[["Sample"]])]
            SS[["Sample"]] <- factor(SS[["Sample"]], levels = original.sample.names, labels = sample.names)
            
            if (all(sapply(SS[c("min.stat", "max.stat", "mean.stat")], is.na))) 
                .err(sprintf("no balance statistics to display. This can occur when `%s = FALSE` and `quick = TRUE` in the original call to `bal.tab()`",
                             STATS[[s]]$disp_stat))
            
            missing.stat <- all(is.na(SS[["mean.stat"]]))
            if (missing.stat) {
                .err(sprintf("%s cannot be displayed. This can occur when %s `FALSE` and `quick = TRUE` in the original call to `bal.tab()`",
                             word_list(firstup(STATS[[s]]$balance_tally_for)),
                             word_list(STATS[[s]]$disp_stat, and.or = "and", is.are = TRUE,
                                       quotes = "`")))
            }
            
            gone <- character(0)
            for (i in sample.vals) {
                if (all(sapply(SS[SS[["Sample"]] == i, c("min.stat", "max.stat", "mean.stat")], is.na))) {
                    gone <- c(gone, i)
                    if (i == sample.names["Unadjusted"] && !adj_only) {
                        .wrn("unadjusted values are missing. This can occur when `un = FALSE` and `quick = TRUE` in the original call to `bal.tab()`")
                    }
                    SS <- SS[SS[["Sample"]] != i,]
                }
            }
            
            dec <- FALSE
            
            if (is_not_null(plot.list[[1]])) var.order <- plot.list[[1]]
            
            if (is_not_null(var.order)) {
                if (inherits(var.order, "love.plot")) {
                    old.vars <- levels(var.order$data$var)
                    old.vars[endsWith(old.vars, "*")] <- substr(old.vars[endsWith(old.vars, "*")], 1, nchar(old.vars[endsWith(old.vars, "*")])-1)
                    if (!all(SS[["var"]] %in% old.vars)) {
                        .wrn("the `love.plot` object in `var.order` doesn't have the same variables as the current input. Ignoring `var.order`")
                        var.order <- NULL
                    }
                    else {
                        SS[["var"]] <- factor(SS[["var"]], levels = old.vars[old.vars %in% SS[["var"]]])
                    }
                }
                else if (tolower(var.order) == "alphabetical") {
                    if ("time" %in% facet) {
                        covnames0 <- make_list(length(unique(SS[["time"]])))
                        for (i in seq_along(covnames0)) {
                            covnames0[[i]] <- {
                                if (i == 1) sort(levels(SS[["var"]][SS[["time"]] == i]))
                                else sort(setdiff(levels(SS[["var"]][SS[["time"]] == i]),
                                                  unlist(covnames0[seq_along(covnames0) < i])))
                            }
                        }
                        covnames <- unlist(covnames0)
                    }
                    else {
                        covnames <- sort(levels(SS[["var"]]))
                    }
                    
                    SS[["var"]] <- factor(SS[["var"]], levels = c(rev(setdiff(covnames, distance.names)),
                                                                  sort(distance.names, decreasing = TRUE)))
                }
                else if (var.order %in% ua) {
                    if (var.order %in% gone) {
                        .wrn(sprintf("`var.order` was set to %s but no %s %s were calculated. Ignoring `var.order`",
                                     add_quotes(tolower(var.order)),
                                     tolower(var.order),
                                     STATS[[s]]$balance_tally_for))
                        var.order <- NULL
                    }
                    else {
                        v <- as.character(SS[["var"]][order(SS[["mean.stat"]][SS[["Sample"]] == sample.names[var.order]], 
                                                            decreasing = dec, na.last = FALSE)])
                        
                        SS[["var"]] <- factor(SS[["var"]], 
                                              levels = c(setdiff(v, distance.names), 
                                                         sort(distance.names, decreasing = TRUE)))
                    }
                }
                
            }
            if (is_null(var.order)) {
                covnames <- as.character(unique(SS[["var"]]))
                SS[["var"]] <- factor(SS[["var"]], levels = c(rev(setdiff(covnames, distance.names)),
                                                              sort(distance.names, decreasing = TRUE)))
            }
            
            if (s == "mean.diffs" && any(base::abs(SS[["max.stat"]]) > 5, na.rm = TRUE)) {
                .wrn("large mean differences detected; you may not be using standardized mean differences for continuous variables")
            }
            if (length(stats) == 1 && drop.missing) SS <- SS[!is.na(SS[["min.stat"]]),]
            SS[["stat"]] <- SS[["mean.stat"]]
        }
        else {
            SS <- do.call("rbind", 
                          lapply(col.sample.names,
                                 function(w) data.frame(var = variable.names,
                                                        type = B[["Type"]],
                                                        stat = B[[ifelse(is_null(Agg.Fun), paste.(STATS[[s]]$bal.tab_column_prefix, w),
                                                                         paste.(Agg.Fun, STATS[[s]]$bal.tab_column_prefix, w))]],
                                                        Sample = switch(w, "Un"= "Unadjusted",
                                                                        "Adj" = "Adjusted", w),
                                                        B[facet],
                                                        row.names = NULL,
                                                        stringsAsFactors = TRUE)
                          ))
            
            
            
            sample.vals <- sample.names[levels(SS[["Sample"]])]
            SS[["Sample"]] <- factor(SS[["Sample"]], levels = original.sample.names, labels = sample.names)
            
            missing.stat <- all(is.na(SS[["stat"]]))
            if (missing.stat) {
                .err(sprintf("%s cannot be displayed. This can occur when %s `FALSE` and `quick = TRUE` in the original call to `bal.tab()`",
                             word_list(firstup(STATS[[s]]$balance_tally_for)), 
                             word_list(STATS[[s]]$disp_stat, and.or = "and", is.are = TRUE,
                                       quotes = "`")))
            }
            
            gone <- character(0)
            for (i in sample.vals) {
                if (all(is.na(SS[["stat"]][SS[["Sample"]] == i]))) {
                    gone <- c(gone, i)
                    if (!adj_only && i == sample.names["Unadjusted"]) {
                        .wrn("unadjusted values are missing. This can occur when `un = FALSE` and `quick = TRUE` in the original call to `bal.tab()`")
                    }
                    SS <- SS[SS[["Sample"]] !=i ,]
                }
            }
            
            if (abs) {
                SS[["stat"]] <- abs_(SS[["stat"]], ratio = s == "variance.ratios")
            }
            dec <- FALSE
            
            if (is_not_null(plot.list[[1]])) var.order <- plot.list[[1]]
            
            #Apply var.order
            if (is_not_null(var.order)) {
                if (inherits(var.order, "love.plot")) {
                    old.vars <- levels(var.order$data$var)
                    old.vars[endsWith(old.vars, "*")] <- substr(old.vars[endsWith(old.vars, "*")], 1, nchar(old.vars[endsWith(old.vars, "*")])-1)
                    if (!all(SS[["var"]] %in% old.vars)) {
                        .wrn("the `love.plot` object in `var.order` doesn't have the same variables as the current input. Ignoring `var.order`")
                        var.order <- NULL
                    }
                    else {
                        SS.var.levels <- old.vars[old.vars %in% SS[["var"]]]
                    }
                }
                else if (tolower(var.order) == "alphabetical") {
                    if ("time" %in% facet) {
                        covnames0 <- make_list(length(unique(SS[["time"]])))
                        for (i in seq_along(covnames0)) {
                            covnames0[[i]] <- {
                                if (i == 1) sort(levels(SS[["var"]][SS[["time"]] == i]))
                                else sort(setdiff(levels(SS[["var"]][SS[["time"]] == i]),
                                                  unlist(covnames0[seq_along(covnames0) < i])))
                            }
                        }
                        covnames <- unlist(covnames0)
                    }
                    else {
                        covnames <- sort(levels(SS[["var"]]))
                    }
                    
                    SS.var.levels <- c(rev(setdiff(covnames, distance.names)), sort(distance.names, decreasing = TRUE))
                    
                }
                else if (var.order %in% ua) {
                    if (var.order %in% gone) {
                        .wrn(sprintf("`var.order` was set to %s, but no %s %s were calculated. Ignoring `var.order`",
                                     add_quotes(tolower(var.order)),
                                     tolower(var.order),
                                     STATS[[s]]$balance_tally_for))
                        var.order <- NULL
                    }
                    else {
                        v <- as.character(SS[["var"]][order(SS[["stat"]][SS[["Sample"]] == sample.names[var.order]], 
                                                            decreasing = dec, na.last = FALSE)])
                        SS.var.levels <- c(setdiff(v,  distance.names), sort(distance.names, decreasing = TRUE))
                    }
                }
                
            }
            if (is_null(var.order)) {
                covnames <- as.character(unique(SS[["var"]])) #Don't use levels here to preserve original order
                SS.var.levels <- c(rev(setdiff(covnames, distance.names)), sort(distance.names, decreasing = TRUE))
            }
            SS[["var"]] <- factor(SS[["var"]], levels = SS.var.levels)
            
            SS[["Sample"]] <- SS[["Sample"]][, drop = TRUE]
            
            if (s == "mean.diffs" && any(base::abs(SS[["stat"]]) > 5, na.rm = TRUE)) {
                .wrn("large mean differences detected; you may not be using standardized mean differences for continuous variables")
            }
            if (length(stats) == 1 && drop.missing) SS <- SS[!is.na(SS[["stat"]]),]
            
            if (is_not_null(sub.B)) {
                #Add subclass statistics when disp.subclass = TRUE
                SS.sub <- do.call("rbind", 
                                  lapply(subclass.names,
                                         function(w) data.frame(var = variable.names,
                                                                type = B[["Type"]],
                                                                stat = sub.B[[paste.(STATS[[s]]$bal.tab_column_prefix, w)]],
                                                                Sample = w,
                                                                row.names = NULL,
                                                                stringsAsFactors = TRUE)
                                  ))
                SS.sub[["Sample"]] <- factor(SS.sub[["Sample"]], levels = subclass.names, labels = subclass.names)
                if (abs) {
                    SS.sub[["stat"]] <- abs_(SS.sub[["stat"]], ratio = s == "variance.ratios")
                }
                SS <- rbind(SS, SS.sub)
            }
        }
        
        SS <- SS[order(SS[["var"]], na.last = FALSE),]
        SS[["var"]] <- SS[["var"]][, drop = TRUE]
        
        #Make the plot
        baseline.xintercept <- STATS[[s]]$baseline.xintercept
        threshold.xintercepts <- {
            if (is_null(thresholds[[s]])) NULL
            else STATS[[s]]$threshold.xintercepts(thresholds[[s]], abs)
        }
        xlab <- STATS[[s]]$love.plot_xlab(abs = abs, binary = attr(x, "print.options")$binary,
                                          continuous = attr(x, "print.options")$continuous,
                                          var_type = B[["Type"]],
                                          stars = stars)
        SS[["var"]] <- STATS[[s]]$love.plot_add_stars(SS[["var"]], 
                                                      variable.names = variable.names,
                                                      binary = attr(x, "print.options")$binary,
                                                      continuous = attr(x, "print.options")$continuous,
                                                      var_type = B[["Type"]],
                                                      stars = stars,
                                                      star_char = args$star_char)
        scale_Statistics <- STATS[[s]]$love.plot_axis_scale
        
        apply.limits <- FALSE
        SS[["on.border"]] <- FALSE
        if (is_not_null(limits[[s]])) {
            if (limits[[s]][2] < limits[[s]][1]) {
                limits[[s]] <- c(limits[[s]][2], limits[[s]][1])
            }
            
            if (limits[[s]][1] >= baseline.xintercept) limits[[s]][1] <- baseline.xintercept - .05*limits[[s]][2]
            if (limits[[s]][2] <= baseline.xintercept) limits[[s]][2] <- baseline.xintercept - .05*limits[[s]][1]
            
            if (identical(scale_Statistics, ggplot2::scale_x_log10)) limits[[s]][limits[[s]] <= 1e-2] <- 1e-2
            
            if (agg.range) {
                if (any(SS[["mean.stat"]] < limits[[s]][1], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["mean.stat"]] < limits[[s]][1]] <- TRUE
                    SS[["mean.stat"]][SS[["mean.stat"]] < limits[[s]][1]] <- limits[[s]][1]
                    SS[["max.stat"]][SS[["max.stat"]] < limits[[s]][1]] <- limits[[s]][1]
                    SS[["min.stat"]][SS[["min.stat"]] < limits[[s]][1]] <- limits[[s]][1]
                }
                if (any(SS[["mean.stat"]] > limits[[s]][2], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["mean.stat"]] > limits[[s]][2]] <- TRUE
                    SS[["mean.stat"]][SS[["mean.stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    SS[["max.stat"]][SS[["max.stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    SS[["min.stat"]][SS[["min.stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
                # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
            }
            else {
                if (any(SS[["stat"]] < limits[[s]][1], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["stat"]] < limits[[s]][1]] <- TRUE
                    SS[["stat"]][SS[["stat"]] < limits[[s]][1]] <- limits[[s]][1]
                }
                if (any(SS[["stat"]] > limits[[s]][2], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["stat"]] > limits[[s]][2]] <- TRUE
                    SS[["stat"]][SS[["stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
            }
            
            apply.limits <- TRUE
        }
        
        lp <- ggplot2::ggplot(aes(y = .data$var, x = .data$stat, group = .data$Sample), data = SS) +
            ggplot2::theme(panel.background = element_rect(fill = "white"),
                           axis.text.x = element_text(color = "black"),
                           axis.text.y = element_text(color = "black"),
                           panel.border = element_rect(fill = NA, color = "black"),
                           plot.background = element_blank(),
                           legend.background = element_blank(),
                           legend.key = element_blank()
            ) +
            ggplot2::scale_shape_manual(values = shapes) +
            ggplot2::scale_size_manual(values = size) +
            ggplot2::scale_discrete_manual(aesthetics = "stroke", values = stroke) +
            ggplot2::scale_fill_manual(values = fill) +
            ggplot2::scale_color_manual(values = colors) +
            ggplot2::labs(y = NULL, x = wrap(xlab, wrap))
        
        lp <- lp + ggplot2::geom_vline(xintercept = baseline.xintercept,
                                       linetype = 1, color = "gray5")
        
        if (is_not_null(threshold.xintercepts)) {
            lp <- lp + ggplot2::geom_vline(xintercept = threshold.xintercepts,
                                           linetype = 2, color = "gray8")
        }
        
        if (agg.range) {
            position.dodge <- ggplot2::position_dodge(.5*(size0[1]/3))
            if (line) { #Add line except to distance
                f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
                lp <- lp + ggplot2::layer(geom = "path", data = f, 
                                          position = position.dodge, 
                                          stat = "identity",
                                          mapping = aes(x = .data$mean.stat, color = .data$Sample), 
                                          params = list(linewidth = size0[1]*.8/3, na.rm = TRUE,
                                                        alpha = alpha))
            }
            
            lp <- lp +
                ggplot2::geom_linerange(aes(y = .data$var, xmin = .data$min.stat, xmax = .data$max.stat,
                                            color = .data$Sample), position = position.dodge,
                                        linewidth = size0[1]*.8/3,
                                        alpha = alpha, 
                                        orientation = "y",
                                        show.legend = FALSE,
                                        na.rm = TRUE) +
                ggplot2::geom_point(aes(y = .data$var, 
                                        x = .data$mean.stat, 
                                        shape = .data$Sample,
                                        size = .data$Sample,
                                        stroke = .data$Sample,
                                        color = .data$Sample),
                                    fill = "white", na.rm = TRUE,
                                    alpha = alpha,
                                    position = position.dodge)
            
        }
        else {
            if (is_not_null(sub.B)) {
                SS.sub <- SS[SS[["Sample"]] %in% subclass.names,]
                SS.sub[["Sample"]] <- SS.sub[["Sample"]][, drop = TRUE]
                
                SS <- SS[SS[["Sample"]] %nin% subclass.names,]
                SS[["Sample"]] <- SS[["Sample"]][, drop = TRUE]
            }
            if (isTRUE(line)) { #Add line except to distance
                f <- function(q) {is.na(q[["stat"]])[q$type == "Distance"] <- TRUE; q}
                lp <- lp + ggplot2::layer(geom = "path", data = f(SS),
                                          position = "identity", stat = "identity",
                                          mapping = aes(color = .data$Sample),
                                          params = list(linewidth = size0[1]*.8/3,
                                                        na.rm = TRUE,
                                                        alpha = alpha))
            }
            lp <- lp + ggplot2::geom_point(data = SS, aes(shape = .data$Sample,
                                                          size = .data$Sample,
                                                          stroke = .data$Sample,
                                                          color = .data$Sample),
                                           fill = "white", 
                                           na.rm = TRUE,
                                           alpha = alpha)
            if (is_not_null(sub.B)) {
                #Add subclass label text
                lp <- lp + ggplot2::geom_text(data = SS.sub,
                                              mapping = aes(label = .data$Sample),
                                              size = 2.5*size0[1]/3, na.rm = TRUE)
            }
        }
        
        if (!drop.distance && is_not_null(distance.names)) {
            lp <- lp + ggplot2::geom_hline(linetype = 1, color = "black",
                                           yintercept = nunique(SS[["var"]]) - length(distance.names) + .5)
        }
        
        if (apply.limits) {
            lp <- lp + scale_Statistics(limits = limits[[s]], expand = c(0, 0))
        }
        else {
            lp <- lp + scale_Statistics()
        }
        
        if (isFALSE(grid)) {
            lp <- lp + ggplot2::theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank())
        }
        else {
            lp <- lp + ggplot2::theme(panel.grid.major = element_line(color = "gray87"),
                                      panel.grid.minor = element_line(color = "gray90"))
        }
        
        if (is_not_null(facet)) {
            lp <- lp + ggplot2::facet_grid(reformulate(facet, "."), drop = FALSE) + ggplot2::labs(x = xlab)
        }
        
        class(lp) <- c(class(lp), "love.plot")
        plot.list[[s]] <- lp
    }
    
    if (length(stats) > 1 || isTRUE(args$use.grid)) {
        
        position <- {
            if (!chk::vld_string(position)) NA_character_
            else match_arg(position, 
                           c("right", "left", "top", "bottom", "none"))
        }
        
        #Process labels
        if (isTRUE(labels)) labels <- LETTERS[seq_along(plot.list)]
        else if (is_null(labels) || isFALSE(labels)) labels <- NULL
        else if (!is.atomic(labels) || length(labels) != length(plot.list)) {
            .wrn("`labels` must be `TRUE` or a string with the same length as `stats`. Ignoring `labels`")
            labels <- NULL
        }
        else labels <- as.character(labels)
        
        plots.to.combine <- plot.list
        for (i in seq_along(plots.to.combine)) {
            plots.to.combine[[i]] <- {
                if (i > 1) {
                    plots.to.combine[[i]] + 
                        ggplot2::theme(axis.text.y=element_blank(),
                                       axis.ticks.y=element_blank(),
                                       legend.position = "none")
                }
                else {
                    plots.to.combine[[i]] + ggplot2::theme(legend.position = "none")
                }
            }
            
            if (is_not_null(labels)) {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + ggplot2::labs(title = labels[i])
            }
            
            if (is_not_null(themes[[stats[i]]])) {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + themes[[stats[i]]]
            }
        }
        
        g <- ggarrange_simple(plots = plots.to.combine, nrow = 1)
        title.grob <- grid::textGrob(title, gp = grid::gpar(fontsize=13.2))
        subtitle.grob <- grid::textGrob(subtitle, gp = grid::gpar(fontsize=13.2))
        
        if (position == "none") {
            p <- gridExtra::arrangeGrob(grobs = list(g), nrow = 1)
        }
        else {
            if (any(!get_from_STATS("adj_only")[stats])) {
                legend.to.get <- which(!get_from_STATS("adj_only")[stats])[1]
            }
            else legend.to.get <- 1
            
            legg <- ggplot2::ggplotGrob(plots.to.combine[[legend.to.get]] + ggplot2::theme(legend.position = position))
            leg <- legg$grobs[[which(legg$layout$name == "guide-box")]]
            
            if (position == "left") {
                p <- gridExtra::arrangeGrob(grobs = list(leg, g), nrow = 1, 
                                            widths = grid::unit.c(sum(leg$widths), grid::unit(1, "npc") - sum(leg$widths)))
            }
            else if (position == "right") {
                p <- gridExtra::arrangeGrob(grobs = list(g, leg), nrow = 1, 
                                            widths = grid::unit.c(grid::unit(1, "npc") - sum(leg$widths), sum(leg$widths)))
            }
            else if (position == "top") {
                p <- gridExtra::arrangeGrob(grobs = list(leg, g), nrow = 2,
                                            heights = grid::unit.c(sum(leg$heights), grid::unit(1, "npc") - sum(leg$heights)))
            }
            else if (position == "bottom") {
                p <- gridExtra::arrangeGrob(grobs = list(g, leg), nrow = 2,
                                            heights = grid::unit.c(grid::unit(1, "npc") - sum(leg$heights), sum(leg$heights)))
            }
        }
        
        if (is_not_null(subtitle)) {
            p <- gridExtra::arrangeGrob(p, top = subtitle.grob)
        }
        p <- gridExtra::arrangeGrob(p, top = title.grob)
        
        grid::grid.newpage()
        grid::grid.draw(p)
        
        attr(p, "plots") <- plot.list
        class(p) <- c(class(p), "love.plot")
        
        return(invisible(p))
    }
    
    p <- plot.list[[1]] + 
        ggplot2::labs(title = title, subtitle = subtitle) +
        ggplot2::theme(plot.title = element_text(hjust = 0.5),
                       plot.subtitle = element_text(hjust = 0.5),
                       legend.position = position)
    
    if (is_not_null(themes[[1]])) {
        p <- p + themes[[1]]
    }
    
    p
}

#' @exportS3Method autoplot bal.tab
autoplot.bal.tab <- function(object, ...) {
    love.plot(object, ...)
}

#' @exportS3Method plot bal.tab
plot.bal.tab <- function(x, ...) {
    love.plot(x, ...)
}

# Helper functions
isColor <- function(x) {
    tryCatch(is.matrix(grDevices::col2rgb(x)), 
             error = function(e) FALSE)
}

f.recode <- function(f, ...) {
    #Simplified version of forcats::fct_recode
    f <- factor(f)
    new_levels <- unlist(list(...), use.names = TRUE)
    old_levels <- levels(f)
    idx <- match(new_levels, old_levels)
    
    old_levels[idx] <- names(new_levels)
    
    levels(f) <- old_levels
    
    f
}

seq_int_cycle <- function(begin, end, max) {
    seq(begin, end, by = 1) - max * (seq(begin - 1, end - 1, by = 1) %/% max)
}

assign.shapes <- function(colors, default.shape = "circle") {
    if (nunique(colors) < length(colors)) {
        seq_int_cycle(19, 19 + length(colors) - 1, max = 25)
    }
    else {
        rep.int(default.shape, length(colors))
    }
}

shapes.ok <- function(shapes, nshapes) {
    shape_names <- c(
        "circle", paste("circle", c("open", "filled", "cross", "plus", "small")), "bullet",
        "square", paste("square", c("open", "filled", "cross", "plus", "triangle")),
        "diamond", paste("diamond", c("open", "filled", "plus")),
        "triangle", paste("triangle", c("open", "filled", "square")),
        paste("triangle down", c("open", "filled")),
        "plus", "cross", "asterisk"
    )
    shape_nums <- 1:25
    (length(shapes) == 1 || length(shapes) == nshapes) &&
        ((is.numeric(shapes) && all(shapes %in% shape_nums)) ||
             (is.character(shapes) && all(shapes %in% shape_names)))
}

gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}

ggarrange_simple <- function(plots, nrow = NULL, ncol = NULL) {
    #A thin version of egg:ggarrange
    
    gtable_frame <- function (g, width = grid::unit(1, "null"), height = grid::unit(1, "null")) {
        panels <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]),]
        pargins <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]),]
        ll <- unique(panels$l)
        margins <- if (length(ll) == 1) grid::unit(0, "pt") else g$widths[ll[-length(ll)] + 2]
        tt <- unique(panels$t)
        fixed_ar <- g$respect
        if (fixed_ar) {
            ar <- as.numeric(g$heights[tt[1]])/as.numeric(g$widths[ll[1]])
            height <- width * (ar/length(ll))
            g$respect <- FALSE
        }
        core <- g[seq(min(tt), max(tt)), seq(min(ll), max(ll))]
        top <- g[seq(1, min(tt) - 1), seq(min(ll), max(ll))]
        bottom <- g[seq(max(tt) + 1, nrow(g)), seq(min(ll), max(ll))]
        left <- g[seq(min(tt), max(tt)), seq(1, min(ll) - 1)]
        right <- g[seq(min(tt), max(tt)), seq(max(ll) + 1, ncol(g))]
        fg <- grid::nullGrob()
        if (is_not_null(left)) {
            lg <- gtable::gtable_add_cols(left, grid::unit(1, "null"), 0)
            lg <- gtable::gtable_add_grob(lg, fg, 1, l = 1)
        }
        else {
            lg <- fg
        }
        
        if (is_not_null(right)) {
            rg <- gtable::gtable_add_cols(right, grid::unit(1, "null"))
            rg <- gtable::gtable_add_grob(rg, fg, 1, l = ncol(rg))
        }
        else {
            rg <- fg
        }
        
        if (is_not_null(top)) {
            tg <- gtable::gtable_add_rows(top, grid::unit(1, "null"), 0)
            tg <- gtable::gtable_add_grob(tg, fg, t = 1, l = 1)
        }
        else {
            tg <- fg
        }
        
        if (is_not_null(bottom)) {
            bg <- gtable::gtable_add_rows(bottom, grid::unit(1, "null"), 
                                          -1)
            bg <- gtable::gtable_add_grob(bg, fg, t = nrow(bg), l = 1)
        }
        else {
            bg <- fg
        }
        
        grobs <- list(fg, tg, fg, lg, core, rg, fg, bg, fg)
        widths <- grid::unit.c(sum(left$widths), width, sum(right$widths))
        heights <- grid::unit.c(sum(top$heights), height, sum(bottom$heights))
        all <- gtable::gtable_matrix("all", grobs = matrix(grobs, ncol = 3, nrow = 3, byrow = TRUE), 
                                     widths = widths, heights = heights)
        
        all[["layout"]][5, "name"] <- "panel"
        if (fixed_ar) 
            all$respect <- TRUE
        all
    }
    
    n <- length(plots)
    
    grobs <- lapply(plots, ggplot2::ggplotGrob)
    
    if (is_null(nrow) && is_null(ncol)) {
        nm <- grDevices::n2mfrow(n)
        nrow <- nm[1]
        ncol <- nm[2]
    }
    
    hw <- lapply(rep(1, n), grid::unit, "null")
    
    fg <- lapply(seq_along(plots), function(i) gtable_frame(g = grobs[[i]], 
                                                            width = hw[[i]], height = hw[[i]]))
    
    spl <- split(fg, rep(1, n))
    
    rows <- lapply(spl, function(r) do.call(gridExtra::gtable_cbind, r))
    
    gt <- do.call(gridExtra::gtable_rbind, rows)
    
    invisible(gt)
}

bal.tab_class_sequence <- function(b) {
    if (inherits(b, "bal.tab.bin") || inherits(b, "bal.tab.cont")) return(NULL)
    
    b_ <- b[[which(endsWith(names(b), ".Balance"))]][[1]]
    c(class(b)[1], bal.tab_class_sequence(b_))
}

unpack_bal.tab <- function(b) {
    unpack_bal.tab_internal <- function(b) {
        if (inherits(b, "bal.tab.bin") || inherits(b, "bal.tab.cont")) return(b[["Balance"]])
        
        b_ <- b[[which(endsWith(names(b), ".Balance"))]]
        
        b_list <- lapply(b_, function(i) {
            if (inherits(b, "bal.tab.bin") || inherits(b, "bal.tab.cont")) return(i[["Balance"]])
            
            unpack_bal.tab_internal(i)
        })
        
        b_list
    }
    LinearizeNestedList <- function(NList, NameSep) {
        # LinearizeNestedList:
        #
        # https://sites.google.com/site/akhilsbehl/geekspace/
        #         articles/r/linearize_nested_lists_in_r
        #
        # Akhil S Bhel
        # 
        
        if (is.data.frame(NList)) return(NList)
        
        A <- 1
        B <- length(NList)
        
        while (A <= B) {
            Element <- NList[[A]]
            EName <- names(NList)[A]
            if (is.list(Element)) {
                
                Before <- {
                    if (A == 1) NULL
                    else NList[1:(A - 1)]
                }
                
                After <- {
                    if (A == B) NULL
                    else NList[(A + 1):B]
                }
                
                if (is.data.frame(Element)) {
                    Jump <- 1
                }
                else {
                    NList[[A]] <- NULL
                    
                    Element <- LinearizeNestedList(Element, NameSep)
                    names(Element) <- paste(EName, names(Element), sep = NameSep)
                    Jump <- length(Element)
                    NList <- c(Before, Element, After)
                }
            }
            else {
                Jump <- 1
            }
            
            A <- A + Jump
            B <- length(NList)
        }
       
        NList
    }
    
    namesep <- paste(c("|", unlist(lapply(1:20, function(i) sample(LETTERS, 1))), "|"), collapse = "")
    
    out_ <- unpack_bal.tab_internal(b)
    out <- LinearizeNestedList(out_, NameSep = namesep)
    
    attr(out, "namesep") <- namesep
    attr(out, "class_sequence") <- bal.tab_class_sequence(b)
    
    out
}
