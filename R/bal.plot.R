#' @title Visualize Distributional Balance
#' 
#' @description Generates density plots, bar graphs, or scatterplots displaying distributional balance between treatment and covariates using \pkg{ggplot2}.
#' 
#' @param x the object for which balance is to be assessed; can be any object for which there is support in [bal.tab()].
#' @param var.name `character`; the name of the variable whose values are to be plotted. To view distributions of the distance measure (e.g., propensity score), if any, use `"distance"` as the argument unless the distance variable has been named. If there are duplicate variable names across inputs, `bal.plot()` will first look in the covariate `data.frame` from `x`, followed by `addl`, and then `distance`, if any. If not specified, will use the first covariate available with a warning.
#' @param ... other arguments to define the variable, treatment, and weights. Some inputs are required depending on the method. See Additional Arguments. Can also be used to supply the `bw`, `adjust`, `kernel`, and `n` arguments for [ggplot2::geom_density()] and the `bins` argument for [ggplot2::geom_histogram()].
#' @param which whether to display distributional balance for the adjusted (`"adjusted"`) or unadjusted sample (`"unadjusted"`) or both at the same time (`"both"`). When multiple weights are present, the names of the weights can be supplied, too. The default is to display balance for the adjusted sample only unless no weights, subclasses, or matching strata are specified. Multiple values and abbreviations allowed.
#' @param which.sub `numeric`; if subclassification is used, a vector corresponding to the subclass(es) for which the distributions are to be displayed. If `.all` (the default), distributions from all subclasses are displayed in a grid.
#' @param cluster optional; a vector of cluster membership, or the name of a variable in an available data set passed to `bal.plot()` that contains cluster membership.
#' @param which.cluster if clusters are used, which cluster(s) to display. Can be cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names. If `.all` (the default), all clusters are displayed. If `.none`, cluster information is ignored and the marginal distribution of the covariates is displayed.
#' @param imp optional; a vector of imputation indices, or the name of a variable in an available data set passed to `bal.plot()` that contains imputation indices.
#' @param which.imp if imputations are used, which imputations(s) to display. Must be numerical indices for which to display balance. If `.all` (the default), all imputations are displayed. If `.none`, data from all imputations are combined into one distribution.
#' @param which.treat which treatment groups to display. If `NULL` (the default) or `NA`, all treatment groups are displayed.
#' @param which.time for longitudinal treatments, which time points to display. Can be treatment names or time period indices. If `NULL` (the default) or `NA`, all time points are displayed.
#' @param mirror `logical`; if the treatment is binary, the covariate is continuous, and densities or histograms are requested, whether to display mirrored densities/histograms or overlapping densities/histograms. Ignored otherwise.
#' @param type `character`; for binary and multi-category treatments with a continuous covariate, whether to display densities (`"density"`), histograms  (`"histogram"`), or empirical cumulative density function plots (`"ecdf"`). The default is to display densities. Abbreviations are allowed.
#' @param colors a vector of colors for the plotted densities/histograms. See 'Color Specification' at [graphics::par()]. Defaults to the default \pkg{ggplot2} colors.
#' @param grid `logical`; whether gridlines should be shown on the plot. Default is `TRUE`.
#' @param sample.names `character`; new names to be given to the samples (i.e., in place of "Unadjusted Sample" and "Adjusted Sample"). For example, when matching it used, it may be useful to enter `c("Unmatched", "Matched")`.
#' @param position the position of the legend. This can be any value that would be appropriate as an argument to `legend.position` in [ggplot2::theme()].
#' @param facet.formula a `formula` designating which facets should be on the rows and columns. This should be of the "historical" formula interface to [ggplot2::facet_grid()]. If of the form `a ~ b`, `a` will be faceted on the rows and `b` on the columns. To only facet on the rows, provide a one-sided formula with an empty left-hand side. To only facet on the columns, the formula should be of the form `a ~ .` (i.e., with only `.` on the right-hand side). The allowable facets depend on which arguments have been supplied to `bal.plot()`; possible values include `which`, `cluster`, `imp`, and (for longitudinal treatments) `time`. If `NULL`, `bal.plot()` will decide what looks best; this argument exists in case you disagree with its choice.
#' @param disp.means `logical`; for a categorical treatment with a continuous covariate, whether a line should be drawn for each treatment level denoting the (weighted) mean of the covariate. Ignored if `type` is not "density" or "histogram". Default is `FALSE`.
#' @param alpha.weight `logical`; if both the treatment and the covariate are continuous, whether points should be shaded according to their weight. Fainter points are those that have smaller weights. Default is `TRUE`.
#' 
#' @section Additional Arguments:
#' `bal.plot()` works like [bal.tab()] in that it can take a variety of types of inputs and yield the same output for each. Depending on what kind of input is given, different additional parameters are required in `...`. For details on what is required and allowed for each additional input and their defaults, see the help file for the [bal.tab()] method associated with the input. The following are the required additional arguments based on each input type:
#'     
#' * For `matchit` objects: None
#' * For `weightit` objects: None
#' * For `ps`, `ps.cont`, `mnps`, and `iptw` objects: (`stop.method`; see [defaults][bal.tab.ps]).
#' * For `Match` objects: `formula` and `data` or `covs` and `treat`.
#' * For `optmatch` objects: `formula` and `data` or `covs` (`treat` is not required).
#' * For `CBPS` objects: None
#' * For `ebalance` objects: `formula` and `data` or `covs` and `treat`.
#' * For `formula`s: `data`
#' * For `data.frame`s: `treat`
#' * For `designmatch` objects: `formula` and `data` or `covs` and `treat`.
#' * For `sbw` objects: None
#' * For `mimids` and `wimids` objects: None, but an argument to `which.imp` should be specified.
#' * For other objects processed through `bal.tab()`'s default method, whichever arguments are required to identify treatment, variables, and a conditioning method (if any).
#'
#' @returns A `"ggplot"` object, returned invisibly.
#' 
#' @details 
#' `bal.plot()` uses [ggplot2::ggplot()] from the \pkg{ggplot2} package, and (invisibly) returns a `"ggplot"` object. For categorical treatments with continuous covariates or continuous treatments with categorical covariates, density plots are created using [ggplot2::geom_density()], histograms are created using [ggplot2::geom_histogram()], and empirical CDF plots are created using [ggplot2::geom_step()]; for categorical treatments with categorical covariates, bar graphs are created using [ggplot2::geom_bar()]; for continuous treatments with continuous covariates, scatterplots are created using [ggplot2::geom_point()].
#' 
#' For continuous treatments with continuous covariates, four additional lines are presented for aid in balance assessment. The red line is the linear fit line. The blue line is a smoothing curve generated with \pkg{ggplot2}'s [ggplot2::geom_smooth()] with `method = "auto"`. The horizontal black line is a horizontal reference line intercepting the (unweighted) treatment mean. The vertical black line is a reference line intercepting the (unweighted) treatment mean. Balance is indicated by the flatness of both fit lines and whether they pass through the intersection of the two black reference lines.
#' 
#' When multiple plots are to be displayed (i.e., when requesting subclass balance, cluster balance, or imputation balance, or when multiple sets of weights are provided or `which = "both"`, or when treatment is longitudinal), the plots will be displayed in a grid using \pkg{ggplot2}'s [ggplot2::facet_grid()]. Subclassification cannot be used with clusters or multiply imputed data.
#' 
#' To change the plot and axis titles, use [ggplot2::labs()]. Because the output is a `ggplot` object, other elements can be changed using \pkg{ggplot2} functions; see [here](https://stackoverflow.com/questions/61255335/change-legend-generated-by-bal-plot) for an example.
#' 
#' @seealso [bal.tab()], [love.plot()]
#' 
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE)
#' data("lalonde", package = "cobalt")
#' 
#' #Nearest Neighbor Matching
#' m.out <- MatchIt::matchit(treat ~ age + educ + race +  married +
#'                               nodegree + re74 + re75, 
#'                           data = lalonde)
#' 
#' bal.plot(m.out, "age", which = "both")
#' bal.plot(m.out, "re74", which = "both", type = "ecdf")
#' bal.plot(m.out, "race", which = "both")
#' bal.plot(m.out, "distance", which = "both", mirror = TRUE,
#'          type = "histogram", colors = c("white", "black"))
#' 
#' @examplesIf requireNamespace("WeightIt", quietly = TRUE)
#' #PS weighting with a continuous treatment
#' w.out <- WeightIt::weightit(re75 ~ age + I(age^2) + educ + 
#'                                 race + married + nodegree,
#'                             data = lalonde)
#' 
#' bal.plot(w.out, "age", which = "both")
#' bal.plot(w.out, "married", which = "both")
 
#' @rdname bal.plot
#' @export 
bal.plot <- function(x, var.name, ..., which, which.sub = NULL, cluster = NULL, which.cluster = NULL, 
                     imp = NULL, which.imp = NULL, which.treat = NULL, which.time = NULL, 
                     mirror = FALSE, type = "density", colors = NULL, grid = FALSE, sample.names,
                     position = "right", facet.formula = NULL, disp.means = getOption("cobalt_disp.means", FALSE), 
                     alpha.weight = TRUE) {
    
    tryCatch(identity(x), error = function(e) .err(conditionMessage(e)))
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    .alls <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.all)), logical(1L))
    .nones <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.none)), logical(1L))
    if (any(c(.alls, .nones))) {
        .call[.alls] <- expression(NULL)
        .call[.nones] <- expression(NA)
        return(eval.parent(.call))
    }
    
    tryCatch(force(x), error = function(e) .err(conditionMessage(e)))
    
    args <- list(...)
    
    x <- process_obj(x)
    
    X <- x2base(x, ..., cluster = cluster, imp = imp)
    
    if (is_null(X$covs.list)) {
        #Point treatment
        X$covs <- .get_C2(X$covs, addl = X$addl, distance = X$distance, cluster = X$cluster, treat = X$treat,
                         drop = FALSE)
        co.names <- attr(X$covs, "co.names")
        if (missing(var.name)) {
            var.name <- NULL; k = 1
            while (is_null(var.name)) {
                x <- co.names[[k]]
                if ("isep" %nin% x[["type"]]) var.name <- x[["component"]][x[["type"]] == "base"][1]
                else if (k < length(co.names)) k <- k + 1
                else {
                    .err("please specify an argument to `var.name`")
                }
            }
            .msg(sprintf("no `var.name` was provided. Displaying balance for %s", var.name))
        }
        
        var.name_in_name <- vapply(co.names, function(x) {
            var.name %in% x[["component"]][x[["type"]] == "base"] &&
                "isep" %nin% x[["type"]]
        }, logical(1L))
        
        var.name_in_name_and_factor <- vapply(seq_along(co.names), function(x) {
            var.name_in_name[x] && "fsep" %in% co.names[[x]][["type"]]
        }, logical(1L))
        
        if (any(var.name_in_name_and_factor)) {
            X$var <- unsplitfactor(as.data.frame(X$covs[,var.name_in_name_and_factor, drop = FALSE]), 
                                   var.name, sep = attr(co.names, "seps")["factor"])[[1]]
        }
        else if (any(var.name_in_name)) {
            X$var <- X$covs[,var.name]
        }
        else {
            .err(sprintf("%s is not the name of an available covariate", add_quotes(var.name, 2)))
        }
        
        if (get.treat.type(X$treat) != "continuous") {
            X$treat <- treat_vals(X$treat)[X$treat]
        }
    }
    else {
        #Longitudinal
        X$covs.list <- lapply(seq_along(X$covs.list), function(i) {
            .get_C2(X$covs.list[[i]], addl = X$addl.list[[i]], distance = X$distance.list[[i]], cluster = X$cluster,
                   treat = X$treat.list[[i]], drop = FALSE)
        })
        co.names.list <- lapply(X$covs.list, attr, "co.names")
        ntimes <- length(X$covs.list)
        
        if (missing(var.name)) {
            var.name <- NULL; k <- 1; time <- 1
            while (is_null(var.name)) {
                x <- co.names.list[[time]][[k]]
                if ("isep" %nin% x[["type"]]) {
                    var.name <- x[["component"]][x[["type"]] == "base"][1]
                }
                else if (time < ntimes) {
                    if (k <  length(co.names.list[[time]])) {
                        k <- k + 1
                    }
                    else {
                        k <- 1
                        time <- time + 1
                    }
                }
                else if (k < length(co.names.list[[time]]))  {
                    k <- k + 1
                }
                else {
                    .err("please specify an argument to `var.name`")
                }
            }
            .msg(sprintf("no `var.name` was provided. Displaying balance for %s", var.name))
        }
        
        var.list <- make_list(length(X$covs.list))
        appears.in.time <- rep.int(TRUE, length(X$covs.list))
        for (i in seq_along(X$covs.list)) {
            var.name_in_name <- vapply(co.names.list[[i]], function(x) {
                var.name %in% x[["component"]][x[["type"]] == "base"] &&
                    "isep" %nin% x[["type"]]
            }, logical(1L))
            
            var.name_in_name_and_factor <- var.name_in_name & vapply(co.names.list[[i]], function(x) {
                "fsep" %in% x[["type"]]
            }, logical(1L))
            
            if (any(var.name_in_name_and_factor)) {
                var.list[[i]] <- unsplitfactor(as.data.frame(X$covs.list[[i]][,var.name_in_name_and_factor, drop = FALSE]), 
                                               var.name, sep = attr(co.names.list[[i]], "seps")["factor"])[[1]]
            }
            else if (any(var.name_in_name)) {
                var.list[[i]] <- X$covs.list[[i]][,var.name]
            }
            else {
                appears.in.time[i] <- FALSE
            }
        }
        
        if (all(vapply(var.list, is_null, logical(1L)))) {
            .err(sprintf("%s is not the name of an available covariate.", add_quotes(var.name)))
        }
        
        X$var <- unlist(var.list[appears.in.time])
        
        X$time <- rep(which(appears.in.time), times = lengths(var.list[appears.in.time]))
        
        X$treat.list[appears.in.time] <- lapply(X$treat.list[appears.in.time], function(t) {
            if (get.treat.type(t) != "continuous") treat_vals(t)[t] else t
        })
        X$treat <- unlist(X$treat.list[appears.in.time])
        
        treat.names <- {
            if (is_null(names(X$treat.list)[appears.in.time])) which(appears.in.time)
            else names(X$treat.list)[appears.in.time]
        }
        
        if (is_not_null(X$weights)) {
            X$weights <- do.call("rbind", lapply(seq_len(sum(appears.in.time)), function(x) X$weights))
        }
        
        if (is_not_null(X$s.weights)) {
            X$s.weights <- do.call("c", lapply(seq_len(sum(appears.in.time)), function(x) X$s.weights))
        }
        
        if (is_not_null(X$cluster)) {
            X$cluster <- rep(X$cluster, sum(appears.in.time))
        }
        
        if (is_not_null(X$imp)) {
            X$imp <- rep(X$imp, sum(appears.in.time))
        }
    }
    
    weight.names <- {
        if (is_null(X$subclass) && NCOL(X$weights) != 1L) names(X$weights)
        else "adjusted"
    }

    if (is_null(X$s.weights)) {
        X$s.weights <- rep(1, length(X$treat))
    }
    
    if (missing(which)) {
        if (is_not_null(args$un)) {
            .msg("note: `un` is deprecated; please use `which` for the same and added functionality")
            which <- {
                if (args$un) "unadjusted"
                else weight.names
            }
        }
        else {
            which <- {
                if (is_null(X$weights) && is_null(X$subclass)) "unadjusted"
                else weight.names
            }
        }
    }
    else {
        if (is_null(X$weights) && is_null(X$subclass)) {
            which <- "unadjusted"
        }
        else {
            which <- tolower(which)
            which <- match_arg(which, unique(c("adjusted", "unadjusted", "both", weight.names)),
                               several.ok = TRUE)
        }
    }
    
    if (is_not_null(args$size.weight)) {
        .msg("note: `size.weight` is no longer allowed; please use `alpha.weight` for similar functionality")
    }
    
    title <- sprintf("Distributional Balance for %s", add_quotes(var.name))
    subtitle <- NULL
    
    facet <- NULL
    is.categorical.var <- is.factor(X$var) || is.character(X$var) || is_binary(X$var) 
    
    if (is_null(X$subclass) || (length(which) == 1 && which == "unadjusted")) {
        if (is_not_null(which.sub) && !all(is.na(which.sub))) {
            if (is_null(X$subclass)) {
                .wrn("`which.sub` was specified but no subclasses were supplied. Ignoring `which.sub`")
            }
            else {
                .wrn("`which.sub` was specified but only the unadjusted sample was requested. Ignoring `which.sub`")
            }
        }
        
        facet <- "which"
        
        if ("both" %in% which) {
            which <- c("Unadjusted Sample", weight.names)
        }
        else {
            if ("adjusted" %in% which) which <- c(which[which != "adjusted"], weight.names)
            if ("unadjusted" %in% which) which <- c("Unadjusted Sample", which[which != "unadjusted"])
        }
        which <- unique(which)
        
        if (is_null(X$weights)) {
            X$weights <- setNames(data.frame(rep.int(1, length(X$treat))),
                                  "Unadjusted Sample")
        }
        else {
            if (ncol(X$weights) == 1) {
                which[which != "Unadjusted Sample"] <- "Adjusted Sample"
                names(X$weights) <- "Adjusted Sample"
            }
            
            if ("Unadjusted Sample" %in% which) {
                X$weights <- setNames(data.frame(rep.int(max(X$weights), length(X$treat)),
                                                 X$weights),
                                      c("Unadjusted Sample", names(X$weights)))
            }
            
            
        }
        
        X$weights <- X$weights[which]
        
        #Process sample names
        ntypes <- length(which)
        nadj <- sum(which != "Unadjusted Sample")
        if (!missing(sample.names)) {
            if (!is.vector(sample.names, "character")) {
                .wrn("the argument to `sample.names` must be a character vector. Ignoring `sample.names`")
                sample.names <- NULL
            }
            else if (length(sample.names) %nin% c(ntypes, nadj)) {
                .wrn("the argument to `sample.names` must contain as many names as there are sample types, or one fewer. Ignoring `sample.names`")
                sample.names <- NULL
            }
        }
        else sample.names <- NULL
        
        #NULL: all; NA: none
        in.imp <- rep.int(TRUE, length(X$var))
        if (is_not_null(X$imp)) {
            if (is_null(which.imp)) {
                in.imp <- !is.na(X$imp)
                facet <- c("imp", facet)
            }
            else if (!all(is.na(which.imp))) {
                if (!is.numeric(which.imp)) {
                    .err("the argument to `which.imp` must be the indices corresponding to the imputations for which distributions are to be displayed")
                }
                
                if (!all(which.imp %in% seq_len(nlevels(X$imp)))) {
                    .err(paste0("The following inputs to `which.imp` do not correspond to given imputations:\n\t", word_list(which.imp[!which.imp %in% seq_len(nlevels(X$imp))])), tidy = FALSE)
                }
                
                in.imp <- !is.na(X$imp) & X$imp %in% levels(X$imp)[which.imp]
                facet <- c("imp", facet)
                
            }
            
        }
        else if (is_not_null(which.imp) && !all(is.na(which.imp))) {
            .wrn("`which.imp` was specified but no `imp` values were supplied. Ignoring `which.imp`")
        }
        
        in.cluster <- rep.int(TRUE, length(X$var))
        if (is_not_null(X$cluster)) {
            if (is_null(which.cluster)) {
                in.cluster <- !is.na(X$cluster)
                facet <- c("cluster", facet)
            }
            else if (!all(is.na(which.cluster))) {
                if (is.numeric(which.cluster)) {
                    if (!all(which.cluster %in% seq_len(nlevels(X$cluster)))) {
                        .err(paste0("The following inputs to `which.cluster` do not correspond to given clusters:\n\t", word_list(which.cluster[!which.cluster %in% seq_len(nlevels(X$cluster))])), tidy = FALSE)
                    }
                    
                    in.cluster <- !is.na(X$cluster) & X$cluster %in% levels(X$cluster)[which.cluster]
                }
                else if (is.character(which.cluster)) {
                    if (!all(which.cluster %in% levels(X$cluster))) {
                        .err(paste0("The following inputs to `which.cluster` do not correspond to given clusters:\n\t", word_list(which.cluster[which.cluster %nin% levels(X$cluster)])), tidy = FALSE)
                    }
                    
                    in.cluster <- !is.na(X$cluster) & X$cluster %in% which.cluster
                }
                else {
                    .err("the argument to `which.cluster` must be the names or indices corresponding to the clusters for which distributions are to be displayed")
                }
                
                facet <- c("cluster", facet)
            }
        }
        else if (is_not_null(which.cluster)) {
            .wrn("`which.cluster` was specified but no `cluster` values were supplied. Ignoring `which.cluster`")
        }
        
        in.time <- rep.int(TRUE, length(X$var))
        if (is_not_null(X$time)) {
            if (is_null(which.time) || all(is.na(which.time))) {
                in.time <- !is.na(X$time)
            }
            else {
                if (is.numeric(which.time)) {
                    if (!all(which.time %in% seq_along(X$covs.list))) {
                        .err(paste0("The following inputs to `which.time` do not correspond to given time periods:\n\t", word_list(which.time[!which.time %in% seq_along(X$covs.list)])), tidy = FALSE)
                    }
                    
                    if (all(which.time %in% seq_along(X$covs.list)[appears.in.time])) {
                        #nothing; which.time is good
                    }
                    else if (any(which.time %in% seq_along(X$covs.list)[appears.in.time])) {
                        .wrn(sprintf("%s does not appear in time period %s",
                                     add_quotes(var.name), word_list(which.time[!which.time %in% seq_along(X$covs.list)[appears.in.time]], "or")))
                        which.time <- which.time[which.time %in% seq_along(X$covs.list)[appears.in.time]]
                    }
                    else {
                        .err(sprintf("%s does not appear in time period %s",
                                     add_quotes(var.name), word_list(which.time, "or")))
                    }
                    in.time <- !is.na(X$time) & X$time %in% which.time
                }
                else if (is.character(which.time)) {
                    if (!all(which.time %in% treat.names)) {
                        .err(paste0("The following inputs to `which.time` do not correspond to given time periods:\n\t", word_list(which.time[!which.time %in% treat.names])), tidy = FALSE)
                    }
                    
                    if (all(which.time %in% treat.names[appears.in.time])) {
                        #nothing; which.time is good
                    }
                    else if (any(which.time %in% treat.names[appears.in.time])) {
                        time.periods <- word_list(which.time[!which.time %in% treat.names[appears.in.time]], "and")
                        .wrn(sprintf("%s does not appear in the time period%%s corresponding to treatment%%s %s",
                                     add_quotes(var.name), time.periods), n = sum(!which.time %in% treat.names[appears.in.time]))
                        which.time <- which.time[which.time %in% treat.names[appears.in.time]]
                    }
                    else {
                        time.periods <- word_list(which.time, "and")
                        .err(sprintf("%s does not appear in the time period%%s corresponding to treatment%%s %s",
                                     add_quotes(var.name), time.periods), n = length(which.time))
                    }
                    in.time <- !is.na(X$time) & treat.names[X$time] %in% which.time
                }
                else {
                    .err("the argument to `which.time` must be the names or indices corresponding to the time periods for which distributions are to be displayed")
                }
            }
            facet <- c("time", facet)
        }
        else if (is_not_null(which.time)) {
            .wrn("`which.time` was specified but a point treatment was supplied. Ignoring `which.time`")
        }
        
        nobs <- sum(in.imp & in.cluster & in.time)
        if (nobs == 0) .err("no observations to display")
        
        D <- make_list(which)
        for (i in which) {
            D[[i]] <- make_df(c("treat", "var", "weights", "s.weights", "which"), nobs)
            D[[i]]$treat <- X$treat[in.imp & in.cluster & in.time]
            D[[i]]$var <- X$var[in.imp & in.cluster & in.time]
            D[[i]]$weights <-  X$weights[in.imp & in.cluster & in.time, i]
            D[[i]]$s.weights <- X$s.weights[in.imp & in.cluster & in.time]
            D[[i]]$which <- i
            
            #Add columns for additional facets
            if ("imp" %in% facet) D[[i]]$imp <- factor(paste("Imputation", X$imp[in.imp & in.cluster & in.time]))
            if ("cluster" %in% facet) D[[i]]$cluster <- factor(X$cluster[in.imp & in.cluster & in.time])
            if ("time" %in% facet) D[[i]]$time <- factor(paste("Time", X$time[in.imp & in.cluster & in.time]))
        }
        D <- do.call(rbind, D)
        
        D$which <- factor(D$which, levels = which)
        
        if (is_not_null(sample.names)) {
            if (length(sample.names) == nadj) {
                levels(D$which)[levels(D$which) != "Unadjusted Sample"] <- sample.names
            }
            else if (length(sample.names) == ntypes) {
                levels(D$which) <- sample.names
            }
        }
        
    }
    else {
        if (is_not_null(X$cluster)) .err("subclasses are not supported with clusters")
        if (is_not_null(X$imp)) .err("subclasses are not supported with multiple imputations")
        
        #Process sample names
        # ntypes <- as.numeric("both" %in% which)
        if (missing(sample.names)) {
            sample.names <- NULL
        }
        else if (which %nin% c("both", "unadjusted")) {
            .wrn("`sample.names` can only be used with `which = \"both\"` or `\"unadjusted\"` to rename the unadjusted sample when called with subclasses. Ignoring `sample.names`")
            sample.names <- NULL
        }
        else if (!is.vector(sample.names, "character")) {
            .wrn("the argument to `sample.names` must be a character vector. Ignoring `sample.names`")
            sample.names <- NULL
        }
        else if (length(sample.names) != 1) {
            .wrn("the argument to `sample.names` must be of length 1 when called with subclasses. Ignoring `sample.names`")
            sample.names <- NULL
        }
        
        #Get sub.names.good
        sub.names <- levels(X$subclass)
        
        if (is_null(which.sub)) {
            which.sub <- sub.names
        }
        else {
            which.sub.original <- which.sub
            if (anyNA(which.sub)) which.sub <- which.sub[!is.na(which.sub)]
            
            if (is_null(which.sub)) {
                .err(sprintf("the argument to `which.sub` cannot be `.none` or `NA` when `which = %s`",
                             add_quotes(which)))
            }
            
            if (is.character(which.sub)) {
                which.sub <- which.sub[which.sub %in% sub.names]
            }
            else if (is.numeric(which.sub)) {
                which.sub <- sub.names[which.sub[which.sub %in% seq_along(sub.names)]]
            }
            
            if (is_null(which.sub)) {
                .err("the argument to `which.sub` must be `.none`, `.all`, or the valid names or indices of subclasses")
            }
            
            if (!all(which.sub.original %in% which.sub)) {
                w.l <- word_list(setdiff(which.sub.original,  which.sub))
                .wrn(sprintf("%s %s not correspond to any subclass in the object and will be ignored",
                             w.l, if (attr(w.l, "plural")) "do" else "does"))
            }
        }
        
        in.sub <- !is.na(X$subclass) & X$subclass %in% which.sub
        D <- make_df(c("weights", "s.weights", "treat", "var", "subclass"), sum(in.sub))
        D$weights <- 1
        D$s.weights <- X$s.weights[in.sub]
        D$treat <- X$treat[in.sub]
        D$var <- X$var[in.sub]
        D$subclass <- paste("Subclass", X$subclass[in.sub])
        
        if (which == "both") {
            #Make unadjusted sample
            D2 <- make_df(c("weights", "s.weights", "treat", "var", "subclass"), length(X$treat))
            D2$weights <- 1
            D$s.weights <- X$s.weights
            D2$treat <- X$treat
            D2$var <- X$var
            D2$subclass <- rep("Unadjusted Sample", length(X$treat))
            D <- rbind(D2, D, stringsAsFactors = TRUE)
            D$subclass <- relevel(factor(D$subclass), "Unadjusted Sample")
        }
        
        facet <- "subclass"
        
        if (is_not_null(sample.names)) {
            levels(D$subclass)[levels(D$subclass) == "Unadjusted Sample"] <- sample.names
        }
        
    }
    
    treat.type <- get.treat.type(assign.treat.type(D$treat))
    
    D <- na.omit(D[order(D$var),])
    # D <- D[D$weights > 0,]
    
    if (treat.type == "continuous") { #Continuous treatments
        
        D$treat <- as.numeric(D$treat)
        
        if (is.categorical.var) { #Categorical vars
            #Density arguments supplied through ...
            bw <- if_null_then(args$bw, "nrd0")
            adjust <- if_null_then(args$adjust, 1)
            kernel <- if_null_then(args$kernel, "gaussian")
            n <- if_null_then(args$n, 512)
            
            D$var <- factor(D$var)
            cat.sizes <- tapply(rep(1, NROW(D)), D$var, sum)
            smallest.cat <- names(cat.sizes)[which.min(cat.sizes)]
            
            if (is.character(bw)) {
                if (!is.function(get0(paste0("bw.", bw)))) {
                    .err(sprintf("%s is not an acceptable entry to `bw`. See `?stats::density` for allowable options",
                                 add_quotes(bw, "`")))
                }
                bw <- get0(paste0("bw.", bw))(D$treat[D$var == smallest.cat])
            }
            
            #Color
            ntypes <- length(cat.sizes)
            if (is_not_null(args$colours)) colors <- args$colours
            
            if (is_null(colors)) {
                colors <- gg_color_hue(ntypes)
            }
            else {
                if (length(colors) > ntypes) {
                    colors <- colors[seq_len(ntypes)]
                    .wrn(sprintf("only using first %s values in `colors`", ntypes))
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
            
            D$weights <- ave(D[["weights"]] * D[["s.weights"]], 
                             D[c("var", facet)], 
                             FUN = function(x) x/sum(x))
            
            bp <- ggplot2::ggplot(D, mapping = aes(x = .data$treat, fill = .data$var, weight = .data$weights)) + 
                ggplot2::geom_density(alpha = .4, bw = bw, adjust = adjust, kernel = kernel, 
                                      n = n, trim = TRUE, outline.type = "full", stat = StatDensity2) + 
                ggplot2::labs(fill = var.name, y = "Density", x = "Treat", title = title, subtitle = subtitle) +
                ggplot2::scale_fill_manual(values = colors) + 
                ggplot2::geom_hline(yintercept = 0) +
                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .05)))
        }
        else { #Continuous vars
            D$var.mean <- ave_w.m(D[["var"]], D[facet], w = D[["s.weights"]])
            D$treat.mean <- ave_w.m(D[["treat"]], D[facet], w = D[["s.weights"]])
            
            bp <- ggplot2::ggplot(D, mapping = aes(x = .data$var, y = .data$treat, weight = .data$weights * .data$s.weights))
            
            if (identical(which, "Unadjusted Sample") || isFALSE(alpha.weight)) bp <- bp + ggplot2::geom_point(alpha = .9)
            else bp <- bp + ggplot2::geom_point(aes(alpha = .data$weights), show.legend = FALSE) + 
                ggplot2::scale_alpha(range = c(.04, 1))
            
            bp <- bp + 
                ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "firebrick2",
                                     alpha = .4, linewidth = 1.5) + 
                {
                    if (nrow(D) <= 1000)
                        ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "royalblue1",
                                             alpha = .1, linewidth = 1.5) 
                    else
                        ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                                             se = FALSE, color = "royalblue1", alpha = .1, linewidth = 1.5)
                } +
                ggplot2::geom_hline(aes(yintercept = .data$treat.mean), linetype = 1, alpha = .9) + 
                ggplot2::geom_vline(aes(xintercept = .data$var.mean), linetype = 1, alpha = .8) +
                ggplot2::labs(y = "Treat", x = var.name, title = title, subtitle = subtitle)
        }
    }
    else { #Categorical treatments (multinomial supported)
        D$treat <- factor(D$treat)
        
        if (is_null(which.treat)) {
            which.treat <- character(0)
        }
        else if (is.numeric(which.treat)) {
            which.treat <- levels(D$treat)[seq_along(levels(D$treat)) %in% which.treat]
            if (is_null(which.treat)) {
                .wrn("no numbers in `which.treat` correspond to treatment values. All treatment groups will be displayed")
                which.treat <- character(0)
            }
        }
        else if (is.character(which.treat)) {
            which.treat <- levels(D$treat)[levels(D$treat) %in% which.treat]
            if (is_null(which.treat)) {
                .wrn("no names in `which.treat` correspond to treatment values. All treatment groups will be displayed")
                which.treat <- character(0)
            }
        }
        else if (anyNA(which.treat)) {
            which.treat <- character(0)
        }
        else {
            .wrn("the argument to `which.treat` must be `NA`, `NULL`, or a vector of treatment names or indices. All treatment groups will be displayed")
            which.treat <- character(0)
        }
        
        if (is_not_null(which.treat) && !anyNA(which.treat)) D <- D[D$treat %in% which.treat,]
        
        for (i in which(vapply(D, is.factor, logical(1L)))) D[[i]] <- factor(D[[i]])
        
        D$weights <- ave(D[["weights"]] * D[["s.weights"]], 
                         D[c("treat", facet)], 
                         FUN = function(x) x/sum(x))
        
        #Color
        ntypes <- nlevels.treat <- nlevels(D$treat)
        if (is_not_null(args$colours)) colors <- args$colours
        
        if (is_null(colors)) {
            colors <- gg_color_hue(ntypes)
        }
        else {
            if (length(colors) > ntypes) {
                colors <- colors[seq_len(ntypes)]
                .wrn(sprintf("only using first %s values in `colors`", ntypes))
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
        names(colors) <- levels(D$treat)
        
        if (is.categorical.var) { #Categorical vars
            D$var <- factor(D$var)
            bp <- ggplot2::ggplot(D, mapping = aes(x = .data$var, fill = .data$treat, weight = .data$weights)) + 
                ggplot2::geom_bar(position = "dodge", alpha = .8, color = "black") + 
                ggplot2::labs(x = var.name, y = "Proportion", fill = "Treatment", title = title, subtitle = subtitle) + 
                ggplot2::scale_x_discrete(drop=FALSE) + 
                ggplot2::scale_fill_manual(drop=FALSE, values = colors) +
                ggplot2::geom_hline(yintercept = 0) + 
                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .05)))
        }
        else { #Continuous vars
            
            type <- match_arg(type, c("histogram", "density", "ecdf"))
            
            if (type %in% c("ecdf")) {
                mirror <- FALSE
                alpha <- 1
                legend.alpha <- alpha
                expandScale <- ggplot2::expansion(mult = .005)
            }
            else if (nlevels.treat <= 2 && isTRUE(mirror)) {
                posneg <- c(1, -1)
                alpha <- .8
                legend.alpha <- alpha
                expandScale <- ggplot2::expansion(mult = .05)
            }
            else {
                mirror <- FALSE
                posneg <- rep(1, nlevels.treat)
                alpha <- .4
                legend.alpha <- alpha/nlevels.treat
                expandScale <- ggplot2::expansion(mult = c(0, .05))
            }
            
            if (type == "histogram") {
                if (isTRUE(disp.means)) {
                    D$var.mean <- ave_w.m(D[["var"]], D[c("treat", facet)], w = D[["weights"]])
                }
                
                if (!is.numeric(args$bins)) args$bins <- 12
                geom_fun <- function(t) {
                    out <- list(ggplot2::geom_histogram(data = D[D$treat == levels(D$treat)[t],],
                                                        mapping = aes(x = .data$var, y = posneg[t]*ggplot2::after_stat(`count`), 
                                                                      weight = .data$weights,
                                                                      fill = names(colors)[t]),
                                                        alpha = alpha, bins = args$bins, color = "black"),
                                NULL)
                    if (isTRUE(disp.means)) out[[2]] <-
                            ggplot2::geom_segment(data = unique(D[D$treat == levels(D$treat)[t], c("var.mean", facet), drop = FALSE]),
                                                  mapping = aes(x = .data$var.mean, xend = .data$var.mean, y = 0, yend = posneg[t]*Inf), 
                                                  color = if (isTRUE(mirror)) "black" else colors[[t]])
                    clear_null(out)
                }
                ylab <- "Proportion"
                
            }
            else if (type %in% c("ecdf")) {
                
                D$cum.pt <- ave(D[["weights"]], 
                                D[c("treat", facet)], 
                                FUN = function(x) cumsum(x)/sum(x))
                
                #Pad 0 and 1
                extra <- setNames(do.call(expand.grid, c(list(c("top", "bottom")),
                                                         lapply(c("treat", facet), function(i) levels(D[[i]])))),
                                  c("pos_", "treat", facet))
                
                merge.extra <- data.frame(pos_ = c("top", "bottom"),
                                          var = c(-Inf, Inf),
                                          cum.pt = c(0, 1))
                extra <- merge(extra, merge.extra)
                extra[["pos_"]] <- NULL
                
                D[names(D) %nin% names(extra)] <- NULL
                
                D <- rbind(extra[names(D)], D)
                
                geom_fun <- function(t) {
                    
                    ggplot2::geom_step(data = D[D$treat == levels(D$treat)[t],],
                                       mapping = aes(x = .data$var, y = .data$cum.pt, color = names(colors)[t]))
                    
                }
                ylab <- "Cumulative Proportion"
            }
            else {
                #Density arguments supplied through ...
                bw <- if_null_then(args$bw, "nrd0")
                adjust <- if_null_then(args$adjust, 1)
                kernel <- if_null_then(args$kernel, "gaussian")
                n <- if_null_then(args$n, 512)
                
                if (is.character(bw)) {
                    t.sizes <- tapply(rep(1, NROW(D)), D$treat, sum)
                    smallest.t <- names(t.sizes)[which.min(t.sizes)]
                    if (is.function(get0(paste0("bw.", bw)))) {
                        bw <- get0(paste0("bw.", bw))(D$var[D$treat == smallest.t])
                    }
                    else {
                        .err(sprintf("%s is not an acceptable entry to `bw`. See `?stats::density` for allowable options",
                                     add_quotes(bw, "`")))
                    }
                }
                
                if (isTRUE(disp.means)) {
                    D$var.mean <- ave_w.m(D[["var"]], D[c("treat", facet)], w = D[["weights"]])
                }
                
                geom_fun <- function(t) {
                    out <- list(
                        ggplot2::geom_density(data = D[D$treat == levels(D$treat)[t],],
                                              mapping = aes(x = .data$var, y = posneg[t]*ggplot2::after_stat(`density`),
                                                            weight = .data$weights, fill = names(colors)[t]),
                                              alpha = alpha, bw = bw, adjust = adjust,
                                              kernel = kernel, n = n, trim = TRUE,
                                              outline.type = "full", stat = StatDensity2),
                        NULL)
                    if (isTRUE(disp.means)) out[[2]] <-
                            ggplot2::geom_segment(data = unique(D[D$treat == levels(D$treat)[t], c("var.mean", facet), drop = FALSE]),
                                                  mapping = aes(x = .data$var.mean, xend = .data$var.mean, y = 0, yend = posneg[t]*Inf), 
                                                  color = if (isTRUE(mirror)) "black" else colors[[t]])
                    clear_null(out)
                }
                ylab <- "Density"
                
            }
            
            bp <- Reduce("+", c(list(ggplot2::ggplot()),
                                lapply(seq_len(nlevels.treat), geom_fun))) +
                ggplot2::scale_fill_manual(values = colors, guide = ggplot2::guide_legend(override.aes = list(alpha = legend.alpha))) +
                ggplot2::scale_color_manual(values = colors) +
                ggplot2::labs(x = var.name, y = ylab, title = title, subtitle = subtitle,
                              fill = "Treatment", color = "Treatment") + 
                ggplot2::scale_y_continuous(expand = expandScale)
            
            if (isTRUE(mirror)) bp <- bp + ggplot2::geom_hline(yintercept = 0)
        }
    }
    
    #Themes
    bp <- bp + ggplot2::theme(panel.background = element_rect(fill = "white", color = "black"),
                              panel.border = element_rect(fill = NA, color = "black"),
                              plot.background = element_blank(),
                              legend.position = position,
                              legend.key = element_rect(fill = "white", color = "black", linewidth = .25))
    if (isTRUE(grid)) {
        bp <- bp + ggplot2::theme(panel.grid.major = element_line(color = "gray87"),
                                  panel.grid.minor = element_line(color = "gray90"))
    }
    else {
        bp <- bp + ggplot2::theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())
    }
    
    if (is_not_null(facet)) {
        if (is_not_null(facet.formula)) {
            if (!rlang::is_formula(facet.formula)) {
                .err("`facet.formula` must be a formula.")
            }
            test.facet <- invisible(ggplot2::facet_grid(facet.formula))
            if (any(c(names(test.facet$params$rows), names(test.facet$params$cols)) %nin% facet)) {
                .err(sprintf("only %s %%r allowed in `facet.formula`",
                             word_list(facet, quotes = 2)), n = length(facet))
            }
            if ("which" %nin% c(names(test.facet$params$rows), names(test.facet$params$cols))) {
                if (length(which) > 1) {
                    .err("\"which\" must be in the facet formula when the `which` argument refers to more than one sample")
                }
                
                .msg(sprintf("Displaying balance for the %s",
                             if (which %in% c("Adjusted Sample", "Unadjusted Sample")) tolower(which)
                             else paste(which, "sample")))
            }
        }
        else if (length(facet) >= 2) {
            if ("which" %in% facet) {
                facet.formula <- reformulate("which", facet[facet %nin% "which"])
            }
            else if ("imp" %in% facet) {
                facet.formula <- reformulate(facet[facet %nin% "imp"], "imp")
            }
            else {
                facets <- data.frame(facet = facet, length = vapply(facet, function(x) nlevels(D[[x]]), numeric(1L)),
                                     stringsAsFactors = FALSE)
                facets <- facets[with(facets, order(length, facet, decreasing = c(FALSE, TRUE))), ]
                facet.formula <- formula(facets)
            }
        }
        else {
            facet.formula <- reformulate(facet, ".")
        }
        
        bp <- bp + ggplot2::facet_grid(facet.formula, drop = FALSE,
                                       scales = if ("subclass" %in% facet) "free_x" else "fixed")
    }
    
   bp
}

# Helper functions

#Similar to StatDensity but allows for negative weights in density
StatDensity2 <- ggplot2::ggproto("StatDensity2", ggplot2::StatDensity,
                                 compute_group = function(data, scales, bw = "nrd0", adjust = 1, kernel = "gaussian",
                                                          n = 512, trim = FALSE, na.rm = FALSE, flipped_aes = FALSE) {
                                     data <- ggplot2::flip_data(data, flipped_aes)
                                     if (trim) {
                                         range <- range(data$x, na.rm = TRUE)
                                     }
                                     else {
                                         range <- scales[[flipped_names(flipped_aes)$x]]$dimension()
                                     }
                                     
                                     density <- compute_density2(data$x, w = data$weight, from = range[1],
                                                                 to = range[2], bw = bw, adjust = adjust, kernel = kernel, n = n)
                                     density$flipped_aes <- flipped_aes
                                     ggplot2::flip_data(density, flipped_aes)
                                 }
                                 
)

compute_density2 <- function(x, w, from, to, bw = "nrd0", adjust = 1,
                             kernel = "gaussian", n = 512) {
    nx <- length(x)
    if (is.null(w)) {
        w <- rep(1 / nx, nx)
    }
    else {
        w <- w / sum(w)
    }
    
    # if less than 2 points return data frame of NAs and a warning
    if (nx < 2) {
        rlang::warn("Groups with fewer than two data points have been dropped.")
        return(data.frame(
            x = NA_real_,
            density = NA_real_,
            scaled = NA_real_,
            ndensity = NA_real_,
            count = NA_real_,
            n = NA_integer_
        ))
    }
    
    dens <- density_neg_w_safe(x, weights = w, bw = bw, adjust = adjust,
                               kernel = kernel, n = n, from = from, to = to)
    
    data.frame(
        x = dens$x,
        density = dens$y,
        scaled =  dens$y / max(dens$y, na.rm = TRUE),
        ndensity = dens$y / max(dens$y, na.rm = TRUE),
        count =   dens$y * nx,
        n = nx
    )
}
density_neg_w_safe <- function(x, weights, bw = "nrd0", adjust = 1,
                               kernel = "gaussian", n = 512, from, to, ...) {
    
    if (any(weights < 0)) wpos <- pmax(0, weights)
    else wpos <- weights
    
    outpos <- density(x = x, bw = bw, adjust = adjust, kernel = kernel,
                      n = n, from = from,
                      to = to, weights = wpos/sum(wpos), ...)
    
    out <- data.frame(x = outpos$x, y = sum(wpos)*outpos$y)
    
    if (any(weights < 0)) {
        wneg <- -pmin(0, weights)
        
        outneg <- density(x = x, bw = bw, adjust = adjust, kernel = kernel,
                          n = n, from = from,
                          to = to, weights = wneg/sum(wneg), ...)
        out$y <- out$y - sum(wneg)*outneg$y
    }
    
    out
}
