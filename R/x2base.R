#Functions to convert object to base.bal.tab input

x2base <- function(...) UseMethod("x2base")

x2base.matchit <- function(m, ...) {
    A <- list(...)
    
    #Process matchit
    
    #Process data and get imp
    m.data <- if (NROW(m[["model"]][["data"]]) != length(m[["treat"]])) NULL else m[["model"]][["data"]]
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data)) {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data, m.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(m[["treat"]], datalist = list(data, m.data))
    
    #Process covs
    if (is.data.frame(m[["X"]])) {
        covs <- get_covs_from_formula(data = m[["X"]])
    }
    else if (is_not_null(m[["model"]][["model"]])) {
        if (nrow(m[["model"]][["model"]]) == length(treat)) {
            covs <- get_covs_from_formula(m[["formula"]], data = m[["model"]][["model"]])
        }
        else {
            #Recreating covs from model object and m[["X"]]. Have to do this because when 
            #discard != NULL and reestimate = TRUE, cases are lost. This recovers them.
            
            # if (is_not_null(data)) {
            #     covs <- get_covs_from_formula(m[["formula"]], data = m[["model"]][["model"]])
            # }
            # else {
            order <- setNames(attr(m[["model"]][["terms"]], "order"),
                              attr(m[["model"]][["terms"]], "term.labels"))
            assign <- setNames(attr(m[["X"]], "assign"), colnames(m[["X"]]))
            assign1 <- assign[assign %in% which(order == 1)] #Just main effects
            
            dataClasses <- attr(m[["model"]][["terms"]], "dataClasses")
            factors.to.unsplit <- names(dataClasses)[dataClasses %in% c("factor", "character", "logical")]
            f0 <- setNames(lapply(factors.to.unsplit, 
                                  function(x) {
                                      if (dataClasses[x] == "factor")
                                          list(levels = levels(m[["model"]][["model"]][[x]]),
                                               faclev = paste0(x, levels(m[["model"]][["model"]][[x]])))
                                      else 
                                          list(levels = unique(m[["model"]][["model"]][[x]]),
                                               faclev = paste0(x, unique(m[["model"]][["model"]][[x]])))
                                  }),
                           factors.to.unsplit)
            covs <- as.data.frame(m[["X"]][, names(assign1)])
            for (i in factors.to.unsplit) {
                covs <- unsplitfactor(covs, i, sep = "",
                                      dropped.level = f0[[i]][["levels"]][f0[[i]][["faclev"]] %nin% colnames(m[["X"]])])
                if (dataClasses[i] == "logical") covs[[i]] <- as.logical(covs[[i]])
            }
            covs <- get_covs_from_formula(m[["formula"]], data = covs)
            # }
        }
        
    }
    else if (is_not_null(m[["formula"]]) && is_not_null(data)) {
        covs <- get_covs_from_formula(m[["formula"]], data = data)
    }
    else {
        covs <- get_covs_from_formula(data = m[["X"]])
    }
    
    #Get estimand
    if (is_null(estimand <- A[["estimand"]])) {
        estimand <- m[["estimand"]]
    }
    
    #Get method
    if (inherits(m, "matchit.subclass")) {
        if (is_not_null(method <- A[["method"]]) && rlang::is_string(method)) {
            method <- match_arg(method, c("weighting", "subclassification"))
        }
        else method <- "subclassification"
    }
    else {
        method <- "matching"
    }
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, m.data, covs))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, m.data, covs),
                                 obj.distance = m[["distance"]], 
                                 obj.distance.name = "distance")
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with matchit objects.", call. = FALSE)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (method == "subclassification") {
        subclass <- as.factor(m[["subclass"]])
    }
    else subclass <- NULL
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with matchit objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(m[["weights"]]) && !all_the_same(m[["weights"]])) {
        weights <- process_weights(m, A, treat, covs, method, addl.data = list(data, m.data))
        method <- attr(weights, "method")
    }
    else weights <- NULL
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A[["s.weights"]], m[["s.weights"]]))) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, m.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, m.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    discarded <- m[["discarded"]]
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "matchit()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- m[["call"]]
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.ps <- function(ps, ...) {
    A <- list(...)
    
    #Process ps
    if (is_not_null(A) && names(A)[1]=="" && is_null(A[["stop.method"]])) A[["stop.method"]] <- A[[1]]
    if (is_null(A[["stop.method"]]) && is_not_null(A[["full.stop.method"]])) A[["stop.method"]] <- A[["full.stop.method"]]
    
    if (is_not_null(A[["stop.method"]])) {
        if (is.character(A[["stop.method"]])) {
            rule1 <- names(ps[["w"]])[sapply(t(sapply(tolower(A[["stop.method"]]), function(x) startsWith(tolower(names(ps[["w"]])), x))), any)]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(names(ps[["w"]]), and.or = "or", quotes = 2), ".\nUsing all available stop methods instead."))
                rule1 <- names(ps[["w"]])
            }
        }
        else if (is.numeric(A[["stop.method"]]) && any(A[["stop.method"]] %in% seq_along(names(ps[["w"]])))) {
            if (any(!A[["stop.method"]] %in% seq_along(names(ps[["w"]])))) {
                message(paste0("Warning: There are ", length(names(ps[["w"]])), " stop methods available, but you requested ", 
                               word_list(A[["stop.method"]][!A[["stop.method"]] %in% seq_along(names(ps[["w"]]))], and.or = "and"),"."))
            }
            rule1 <- names(ps[["w"]])[A[["stop.method"]] %in% seq_along(names(ps[["w"]]))]
        }
        else {
            warning("stop.method should be ", word_list(names(ps[["w"]]), and.or = "or", quotes = 2), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- names(ps[["w"]])
        }
    }
    else {
        rule1 <- names(ps[["w"]])
    }
    
    s <- names(ps[["w"]])[match(tolower(rule1), tolower(names(ps[["w"]])))]
    
    #Process data and get imp
    ps.data <- ps[["data"]]
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data, ps.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(ps[["treat"]], datalist = list(data, ps.data))
    
    #Process covs
    f <- reformulate(ps[["gbm.obj"]][["var.names"]])
    covs <- get_covs_from_formula(f, data = ps.data)
    
    #Get estimand
    estimand <- ps[["estimand"]]
    
    #Get method
    method <- rep("weighting", length(s))
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, ps.data))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, ps.data, covs),
                                 obj.distance = ps[["ps"]][s], 
                                 obj.distance.name = if (length(s) > 1) paste.("prop.score", substr(s, 1, nchar(s) - 4)) else "prop.score")
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with ps objects.", call. = FALSE)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("'subclasses' are not allowed with ps objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with ps objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(ps, A, treat, covs, method, addl.data = list(data, ps.data), 
                               stop.method = s, estimand = estimand)
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A[["s.weights"]], ps[["sampw"]]))) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, ps.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, ps.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "ps()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    # call <- ps[["parameters"]]
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.mnps <- function(mnps, ...) {
    A <- list(...)
    
    #Process mnps
    if (is_not_null(A)&& names(A)[1]=="" && is_null(A[["stop.method"]])) A[["stop.method"]] <- A[[1]]
    if (is_null(A[["stop.method"]]) == 0 && is_not_null(A[["full.stop.method"]])) A[["stop.method"]] <- A[["full.stop.method"]]
    
    if (is_not_null(A[["stop.method"]])) {
        if (any(is.character(A[["stop.method"]]))) {
            rule1 <- mnps[["stopMethods"]][sapply(t(sapply(tolower(A[["stop.method"]]), function(x) startsWith(tolower(mnps[["stopMethods"]]), x))), any)]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(mnps[["stopMethods"]], and.or = "or", quotes = 2), ".\nUsing all available stop methods instead."))
                rule1 <- mnps[["stopMethods"]]
            }
        }
        else if (is.numeric(A[["stop.method"]]) && any(A[["stop.method"]] %in% seq_along(mnps[["stopMethods"]]))) {
            if (any(!A[["stop.method"]] %in% seq_along(mnps[["stopMethods"]]))) {
                message(paste0("Warning: There are ", length(mnps[["stopMethods"]]), " stop methods available, but you requested ", 
                               word_list(A[["stop.method"]][!A[["stop.method"]] %in% seq_along(mnps[["stopMethods"]])], and.or = "and"),"."))
            }
            rule1 <- mnps[["stopMethods"]][A[["stop.method"]] %in% seq_along(mnps[["stopMethods"]])]
        }
        else {
            warning("stop.method should be ", word_list(mnps[["stopMethods"]], and.or = "or", quotes = 2), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- mnps[["stopMethods"]]
        }
    }
    else {
        rule1 <- mnps[["stopMethods"]]
    }
    
    s <- mnps[["stopMethods"]][match(tolower(rule1), tolower(mnps[["stopMethods"]]))]
    
    #Process data and get imp
    mnps.data <- mnps[["data"]]
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data, mnps.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(mnps[["treatVar"]], datalist = list(data, mnps.data))
    
    #Process covs
    f <- reformulate(mnps[["psList"]][[1]][["gbm.obj"]][["var.names"]])
    covs <- get_covs_from_formula(f, mnps.data)
    
    #Get estimand
    estimand <- mnps[["estimand"]]
    
    #Get method
    method <- rep("weighting", length(s))
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, mnps.data))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, mnps.data))
    
    #Process focal
    focal <- mnps[["treatATT"]]
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with mnps objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with mnps objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(mnps, A, treat, covs, method, addl.data = list(data, mnps.data), 
                               stop.method = s)
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A[["s.weights"]], mnps[["sampw"]]))) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, mnps.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, mnps.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "mnps()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- NULL
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "multi"
    
    return(X)
}
x2base.ps.cont <- function(ps.cont, ...) {
    A <- list(...)
    
    #Process data and get imp
    ps.data <- ps.cont[["data"]]
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data, ps.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(ps.cont[["treat"]], datalist = list(data, ps.data))
    
    #Process covs
    f <- reformulate(ps.cont[["gbm.obj"]][["var.names"]])
    covs <- get_covs_from_formula(f, ps.data)
    
    #Get estimand
    
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, ps.data))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, ps.data))
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with ps.cont objects.", call. = FALSE)
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with ps.cont objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with ps.cont objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(ps.cont, A, treat, covs, method, addl.data = list(data, ps.data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A[["s.weights"]], ps.cont[["sampw"]]))) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, ps.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, ps.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "ps.cont()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("correlations" %in% stats) {
            s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    # call <- ps.cont[["parameters"]]
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "cont"
    
    return(X)
}
x2base.Match <- function(Match, ...) {
    A <- list(...)
    
    #Process Match
    if (is_not_null(Match) && !is.list(Match)) stop("'Match' object contains no valid matches")
    
    #Process data and get imp
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    t.c <- use.tc.fd(A[["formula"]], data, A[["treat"]], A[["covs"]])
    treat <- process_treat(t.c[["treat"]], datalist = list(data))
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- Match[["estimand"]]
    
    #Get method
    method <- "matching"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, covs))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, covs))
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with Match objects.", call. = FALSE)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with Match objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with Match objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(Match, A, treat, covs, method, addl.data = list(data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["s.weights"]])) {
        stop("Sampling weights are not allowed with Match objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    discarded <- rep(FALSE, length(treat))
    if (is_not_null(Match[["index.dropped"]])) discarded[Match[["index.dropped"]]] <- TRUE
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "Match()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- NULL
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.formula <- function(formula, ...) {
    A <- list(...)

    #Pass to x2base.data.frame, which processes covs as a formula
    A[["covs"]] <- NULL

    X <- do.call(x2base.data.frame, c(list(covs = formula), A))
    return(X)
}
x2base.data.frame <- function(covs, ...) {
    A <- list(...)
    
    #Process data.frame
    
    #Process data and get imp
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    if (is.formula(covs)) A[["treat"]] <- get_treat_from_formula(covs, data, treat = A[["treat"]])
    treat <- process_treat(A[["treat"]], datalist = list(data))
    
    #Process covs
    if (is_null(covs)) {
        stop("'covs' data.frame must be specified.", call. = FALSE)
    }
    if (is.formula(covs)) {
        covs <- get_covs_from_formula(covs, data = data)
        if (is_null(covs)) {
            stop("The right hand side of the formula must contain covariates for which balance is to be assessed.", call. = FALSE)
        }
    }
    if (is_null(attr(covs, "co.names"))) {
        if (is.matrix(covs)) covs <- as.data.frame.matrix(covs)
        covs <- get_covs_from_formula(data = covs)
    }
    # is_(covs, "data.frame", stop = TRUE)
    
    #Get estimand
    estimand <- A[["estimand"]]
    
    #Get method
    specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
    if (is_not_null(A[["weights"]])) {
        specified["weights"] <- TRUE
    }
    if (is_not_null(A[["subclass"]])){
        specified["subclass"] <- TRUE
    }
    if (is_not_null(A[["match.strata"]])) {
        specified["match.strata"] <- TRUE
    }
    
    if (is_null(A[["method"]])) {
        if (specified["match.strata"]) {
            if (sum(specified) > 1) {
                message(word_list(names(specified)[specified]), " are specified. Assuming \"matching\" and using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                A[["weights"]] <- A[["subclass"]] <- NULL
            }
            method <- "matching"
        }
        else if (specified["subclass"]) {
            if (sum(specified) > 1) {
                message(word_list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                A[["weights"]] <- A[["match.strata"]] <- NULL
            }
            method <- "subclassification"
            #weights <- rep(1, nrow(covs))
        }
        else if (specified["weights"]) {
            if (sum(specified) > 1) {
                message(word_list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                A[["match.strata"]] <- A[["subclass"]] <- NULL
            }
            method <- "weighting"
        }
        else {
            method <- "matching"
        }
    }
    else if (length(A[["method"]]) == 1) {
        specified.method <- match_arg(A[["method"]], c("weighting", "matching", "subclassification"))
        if (specified.method == "weighting") {
            if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
                    A[["match.strata"]] <- A[["subclass"]] <- NULL
                }
                method <- "weighting"
            }
            else if (specified["match.strata"]) {
                message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
                A[["subclass"]] <- NULL
                method <- "matching"
            }
            else if (specified["subclass"]) {
                message("method = \"weighting\" is specified, but no weights are present. Assuming \"subclassification\" and using subclass instead.")
                method <- "subclassification"
            }
            else {
                method <- "matching"
            }
        }
        else if (specified.method == "matching") {
            if (specified["match.strata"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                    A[["weights"]] <- A[["subclass"]] <- NULL
                }
                method <- "matching"
            }
            else if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
                    A[["match.strata"]] <- A[["subclass"]] <- NULL
                }
                method <- "matching"
            }
            else if (specified["subclass"]) {
                message("method = \"matching\" is specified, but no weights or match.strata are present. Assuming \"subclassification\" and using subclass instead.")
                method <- "subclassification"
            }
            else {
                method <- "matching"
            }
        }
        else if (specified.method == "subclassification") {
            if (specified["subclass"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    A[["weights"]] <- A[["match.strata"]] <- NULL
                }
                method <- "subclassification"
            }
            else if (specified["match.strata"]) {
                message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
                A[["weights"]] <- NULL
                method <- "matching"
            }
            else if (specified["weights"]) {
                message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"weighting\" and using weights instead.")
                method <- "weighting"
            }
        }
    }
    else {
        specified.method <- match_arg(A[["method"]], c("weighting", "matching", "subclassification"), several.ok = TRUE)
        if (any(specified.method == "subclassification") || specified["subclass"]) {
            stop("Subclassification cannot be specified along with other methods.", call. = FALSE)
        }
        else if (specified["match.strata"]) {
            stop("Only weights can be specified with multiple methods.", call. = FALSE)
        }
        else if (!specified["weights"]) {
            warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
            method <- "matching"
        }
        else {
            #Matching and/or weighting with various weights
            method <- specified.method
            A[["match.strata"]] <- A[["subclass"]] <- NULL
        }
    }
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, covs))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, covs))
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]]) && get.treat.type(treat) != "continuous") {
        focal <- process_focal(focal, treat)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        subclass <- vector.process(subclass, 
                                   datalist = list(data),
                                   name = "subclass", 
                                   which = "subclass membership",
                                   missing.okay = TRUE)
        subclass <- factor(subclass)
        weights <- NULL
    }
    
    #Process match.strata
    else if (is_not_null(match.strata <- A[["match.strata"]])) {
        match.strata <- vector.process(match.strata, 
                                       datalist = list(data),
                                       name = "match.strata", 
                                       which = "matching strata membership",
                                       missing.okay = TRUE)
        weights <- data.frame(weights = strata2weights(match.strata,
                                                       treat = treat,
                                                       estimand = estimand,
                                                       focal = focal))
    }
    #Process weights
    else if (is_not_null(A[["weights"]])) {
        weights <- process_weights(NULL, A, treat, covs, method, addl.data = list(data))
        method <- attr(weights, "method")
    }
    else {
        weights <- NULL
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["s.weights"]])) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    discarded <- A[["discarded"]]
    
    #Process length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp)
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, 
                                       weights = weights, subclass = subclass, 
                                       treat = treat, focal = focal)
        }
        else if ("correlations" %in% stats) {
            s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    return(X)
}
x2base.CBPS <- function(cbps.fit, ...) {
    A <- list(...)
    
    #Process CBPS
    
    #Process data and get imp
    c.data <- cbps.fit[["data"]]
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data, c.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- get_treat_from_formula(cbps.fit[["formula"]], c.data)
    treat <- process_treat(treat, datalist = list(data, c.data))
    
    #Process covs
    covs <- get_covs_from_formula(cbps.fit[["formula"]], c.data)
    
    #Get estimand
    if (is_not_null(estimand <- A[["estimand"]])) {
        stop("'estimand' is not allowed with CBPS objects.", call. = FALSE)
    }
    
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, c.data))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, c.data),
                                 obj.distance = if (get.treat.type(treat) == "binary") cbps.fit[["fitted.values"]], 
                                 obj.distance.name = "prop.score")
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with CBPS objects.", call. = FALSE)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with CBPS objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with CBPS objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(cbps.fit, A, treat, covs, method, addl.data = list(data, c.data), 
                               use.weights = A[["use.weights"]])
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["sampw"]])) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, c.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
        weights <- weights/s.weights #Because CBPS weights contain s.weights in them
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, c.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "CBPS()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal,
                                       quietly = TRUE)
        }
        else if ("correlations" %in% stats) {
            s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- cbps.fit[["call"]]
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    return(X)
}
x2base.ebalance <- function(ebalance, ...) {
    A <- list(...)
    
    #Process ebalance
    
    #Process data and get imp
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    t.c <- use.tc.fd(A[["formula"]], data, A[["treat"]], A[["covs"]])
    treat <- process_treat(t.c[["treat"]], datalist = list(data))
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- "ATT"
    
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, covs))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, covs))
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with ebalance objects.", call. = FALSE)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with ebalance objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with ebalance objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(ebalance, A, treat, covs, method, addl.data = list(data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["sampw"]])) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "ebalance()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.optmatch <- function(optmatch, ...) {
    A <- list(...)
    
    #Process optmatch
    if (all(is.na(optmatch))) stop("The 'optmatch' object contains no valid matches.", call. = FALSE)
    
    #Process data and get imp
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    t.c <- use.tc.fd(A[["formula"]], data = data, covs = A[["covs"]],
                     treat = if_null_then(A[["treat"]], as.numeric(attr(optmatch, "contrast.group"))))
    treat <- process_treat(t.c[["treat"]], datalist = list(data))
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- A[["estimand"]]
    
    #Get method
    method <- "matching"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, covs))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, covs))
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with optmatch objects.", call. = FALSE)
    }
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with optmatch objects.", call. = FALSE)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with optmatch objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(optmatch, A, treat, covs, method, addl.data = list(data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["s.weights"]])) {
        stop("Sampling weights are not allowed with optmatch objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = paste0(deparse1(attr(optmatch, "call")[[1]]), "()"))
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- attr(optmatch, "call")
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.cem.match <- function(cem.match, ...) {
    A <- list(...)
    
    #Process cem.match
    if (is_(cem.match, "cem.match.list")) {
        cem.match[["vars"]] <- cem.match[[1]][["vars"]]
        cem.match[["baseline.group"]] <- cem.match[[1]][["baseline.group"]]
        cem.match[["groups"]] <- unlist(lapply(cem.match[vapply(cem.match, is_, logical(1L), "cem.match")], `[[`, "groups"))
        cem.match[["w"]] <- get.w.cem.match(cem.match)
    }
    if (all(check_if_zero(cem.match[["w"]]))) stop("The 'cem.match' object contains no valid matches.", call. = FALSE)
    
    #Process data and get imp
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    if (is_null(data)) {
        stop("An argument to 'data' must be specified with cem.match objects.", call. = FALSE)
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    else if (is_(cem.match, "cem.match.list") && sum(vapply(cem.match, is_, logical(1L), "cem.match")) != 1) {
        stop("An argument to 'imp' must be specified or the argument to data must be a mids object.", call. = FALSE)
    }
    
    #Process treat
    t.c <- use.tc.fd(data = data, treat = cem.match[["groups"]], 
                     covs = cem.match[["vars"]])
    treat <- process_treat(t.c[["treat"]], datalist = list(data))
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- A[["estimand"]]
    
    #Get method
    method <- "matching"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, covs))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, covs))
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with cem.match objects.", call. = FALSE)
    }
    
    #Process focal
    focal <- cem.match[["baseline.group"]]
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with cem.match objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(cem.match, A, treat, covs, method, addl.data = list(data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["s.weights"]])) {
        stop("Sampling weights are not allowed with cem.match objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "cem()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.weightit <- function(weightit, ...) {
    A <- list(...)
    
    #Process CBPS
    
    #Process data and get imp
    d.e.in.w <- vapply(c("covs", "exact", "by", "moderator"), function(x) is_not_null(weightit[[x]]), logical(1L))
    if (any(d.e.in.w)) weightit.data <- do.call("data.frame", unname(weightit[c("covs", "exact", "by", "moderator")[d.e.in.w]]))
    else weightit.data <- NULL
    
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data, weightit.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(weightit[["treat"]], datalist = list(data, weightit.data))
    
    #Process covs
    if (is_null(covs <- weightit[["covs"]])) stop("No covariates were specified in the weightit object.", call. = FALSE)
    covs <- get_covs_from_formula(data = covs)
    
    #Get estimand
    estimand <- weightit[["estimand"]]
    
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, weightit.data))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, weightit.data),
                                 obj.distance = if (get.treat.type(treat) == "binary") weightit[["ps"]], 
                                 obj.distance.name = "prop.score")
    
    #Process focal
    focal <- weightit[["focal"]]
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with weightit objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with weightit objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(weightit, A, treat, covs, method, addl.data = list(data, weightit.data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A[["s.weights"]], weightit[["s.weights"]]))) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, weightit.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, weightit.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    discarded <- weightit[["discarded"]]
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "weightit()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
        else if ("correlations" %in% stats) {
            s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- weightit[["call"]]
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    return(X)
}
x2base.designmatch <- function(dm, ...) {
    A <- list(...)
    
    #Process designmatch
    if (all(c("id_1", "id_2") %in% names(dm))) {
        stop("Balance cannot currently be checked on a nonbipartite match.", call. = FALSE)
    }
    
    #Process data and get imp
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    t.c <- use.tc.fd(A[["formula"]], data, A[["treat"]], A[["covs"]])
    treat <- process_treat(t.c[["treat"]], datalist = list(data))
    if (is.unsorted(rev(treat))) warning("designmatch requires the input data to be sorted by treatment; the data supplied to bal.tab() was not, indicating a possible coding error.", call. = FALSE)
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- A[["estimand"]]
    
    #Get method
    method <- "matching"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, covs))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, covs))
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with designmatch objects.", call. = FALSE)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with designmatch objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with designmatch objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(dm, A, treat, covs, method, addl.data = list(data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["s.weights"]])) {
        stop("Sampling weights are not allowed with designmatch objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "the matching function in designmatch")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- NULL
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.mimids <- function(mimids, ...) {
    A <- list(...)
    
    #Process mimids
    old_version <- !all(c("object", "models", "approach") %in% names(mimids))
    models <- if (old_version) mimids[["models"]][-1] else mimids[["models"]]
    
    #Process data and get imp
    if (old_version) {
        if (is_(mimids[["original.datasets"]], "mids")) m.data <- imp.complete(mimids[["original.datasets"]])
        else m.data <- imp.complete(mimids[["others"]][["source"]])
    }
    else {
        m.data <- imp.complete(mimids[["object"]])
    }
    
    imp <- m.data[[".imp"]]
    
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data, m.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(unlist(lapply(models, function(m) m[["treat"]])))
    
    #Process covs
    covs <- do.call("rbind", lapply(models, function(m) m[["X"]]))
    covs <- get_covs_from_formula(data = covs)
    
    #Get estimand
    estimand <- models[[1]][["estimand"]]
    
    #Get method
    method <- "matching"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, m.data))
    
    #Process distance
    m.distance <- unlist(lapply(models, function(m) m[["distance"]]))
    
    if (all(is.na(m.distance))) m.distance <- NULL
    
    distance <- process_distance(A[["distance"]], datalist = list(data, m.data),
                                 obj.distance = m.distance, 
                                 obj.distance.name = "distance")
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with mimids objects.", call. = FALSE)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with mimids objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with mimids objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(mimids, A, treat, covs, method, addl.data = list(data, m.data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["s.weights"]])) {
        stop("Sampling weights are not allowed with mimids objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, m.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, min(table(imp)))
    }
    
    #Process discarded
    discarded <- unlist(lapply(models, function(m) m[["discarded"]]))
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "matchthem()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- NULL
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "imp"
    
    return(X)
}
x2base.wimids <- function(wimids, ...) {
    A <- list(...)
    
    #Process wimids
    old_version <- !all(c("object", "models", "approach") %in% names(wimids))
    models <- if (old_version) wimids[["models"]][-1] else wimids[["models"]]
    
    #Process data and get imp
    if (old_version) {
        if (is_(wimids[["original.datasets"]], "mids")) w.data <- imp.complete(wimids[["original.datasets"]])
        else w.data <- imp.complete(wimids[["others"]][["source"]])
    }
    else {
        w.data <- imp.complete(wimids[["object"]])
    }
    
    imp <- w.data[[".imp"]]
    
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data, w.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(unlist(lapply(models, function(w) w[["treat"]])))
    
    #Process covs
    covs <- do.call("rbind", lapply(models, function(w) w[["covs"]]))
    covs <- get_covs_from_formula(data = covs)
    
    #Get estimand
    estimand <- unique(unlist(lapply(models, function(w) w[["estimand"]])))
 
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, w.data))
    
    #Process distance
    w.distance <- unlist(lapply(models, function(m) m[["ps"]]))
    if (all(is.na(w.distance))) w.distance <- NULL
    
    distance <- process_distance(A[["distance"]], datalist = list(data, w.data),
                                 obj.distance = if (get.treat.type(treat) == "binary") w.distance, 
                                 obj.distance.name = "prop.score")
    
    #Process focal
    focal <- unique(unlist(lapply(models, function(w) w[["focal"]])))
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with wimids objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with wimids objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(wimids, A, treat, covs, method, addl.data = list(data, w.data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A[["s.weights"]], unlist(lapply(models, function(w) w[["s.weights"]]))))) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, w.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, w.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, min(table(imp)))
    }
    
    #Process discarded
    discarded <- unlist(lapply(models, function(w) w[["discarded"]]))
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "weightthem()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
        else if ("correlations" %in% stats) {
            s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- NULL
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "imp"
    
    return(X)
}
x2base.sbwcau <- function(sbwcau, ...) {
    A <- list(...)
    
    #Process matchit
    
    #Process data and get imp
    sbw.data <- sbwcau[["dat_weights"]][names(sbwcau[["dat_weights"]]) != "weights"]
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              datalist = list(data, sbw.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(sbwcau[["ind"]], datalist = list(data, sbw.data))
    
    #Process covs
    f <- reformulate(sbwcau[["bal"]][["bal_cov"]])
    covs <- get_covs_from_formula(f, data = sbw.data)
    
    #Get estimand
    estimand <- sbwcau[["par"]][["par_est"]]
    
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- process_addl(A[["addl"]], datalist = list(data, sbw.data))
    
    #Process distance
    distance <- process_distance(A[["distance"]], datalist = list(data, sbw.data))
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with sbwcau objects.", call. = FALSE)
    } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
        focal <- switch(toupper(estimand), 
                        "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                        "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                        NULL)
    }
    
    #Process pairwise
    if (get.treat.type(treat) == "binary" && is_null(focal)) {
        if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
        if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with sbwcau objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with sbwcau objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(sbwcau, A, treat, covs, method, addl.data = list(data, sbw.data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["sampw"]])) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, sbw.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, sbw.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, length(treat))
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "sbw()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat, focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs) || anyNA(addl)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    
    #Process output
    X <- initialize_X()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}

#MSMs wth multiple time points
x2base.iptw <- function(iptw, ...) {
    A <- list(...)
    
    #Process iptw
    if (is_not_null(A) && names(A)[1]=="" && is_null(A[["stop.method"]])) A[["stop.method"]] <- A[[1]] #for bal.plot
    if (is_null(A[["stop.method"]]) && is_not_null(A[["full.stop.method"]])) A[["stop.method"]] <- A[["full.stop.method"]]
    available.stop.methods <- names(iptw[["psList"]][[1]][["ps"]])
    if (is_not_null(A[["stop.method"]])) {
        if (any(is.character(A[["stop.method"]]))) {
            rule1 <- available.stop.methods[vapply(available.stop.methods, function(x) any(startsWith(tolower(x), tolower(A[["stop.method"]]))), logical(1L))]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(available.stop.methods, and.or = "or", quotes = 2), ".\nUsing all available stop methods instead."))
                rule1 <- available.stop.methods
            }
        }
        else if (is.numeric(A[["stop.method"]]) && any(A[["stop.method"]] %in% seq_along(available.stop.methods))) {
            if (any(!A[["stop.method"]] %in% seq_along(available.stop.methods))) {
                message(paste0("Warning: There are ", length(available.stop.methods), " stop methods available, but you requested ", 
                               word_list(A[["stop.method"]][!A[["stop.method"]] %in% seq_along(available.stop.methods)], and.or = "and"),"."))
            }
            rule1 <- available.stop.methods[A[["stop.method"]] %in% seq_along(available.stop.methods)]
        }
        else {
            warning("stop.method should be ", word_list(available.stop.methods, and.or = "or", quotes = 2), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- available.stop.methods
        }
    }
    else {
        rule1 <- available.stop.methods
    }
    
    s <- available.stop.methods[match(tolower(rule1), tolower(available.stop.methods))]
    
    #Process data and get imp
    ps.data <- iptw[["psList"]][[1]][["data"]]
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", datalist = list(data), missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat.list
    treat.list <- process_treat.list(lapply(iptw[["psList"]], function(x) x[["treat"]]), datalist = list(data, ps.data))
    
    #Process covs.list
    covs.list <- lapply(iptw[["psList"]], function(x) get_covs_from_formula(reformulate(x[["gbm.obj"]][["var.names"]]), data = x[["data"]]))
    
    #Get estimand
    estimand <- substr(toupper(s), nchar(s)-2, nchar(s))
    
    #Get method
    method <- rep("weighting", length(s))
    
    #Process addl.list 
    addl.list <- process_addl.list(if_null_then(A[["addl.list"]], A[["addl"]]),
                                   datalist = list(data, ps.data),
                                   covs.list = covs.list)
    
    #Process distance
    # ntimes <- iptw[["nFits"]]
    # distance.list <- list.process("distance.list", A[["distance.list"]], ntimes, 
    #                               "the original call to iptw()",
    #                               treat.list,
    #                               covs.list,
    #                               list(data, ps.data))
    # if (is_not_null(distance.list)) {
    #     for (ti in seq_along(distance.list)) {
    #         if (length(s) == 1) {
    #             distance.list[[ti]] <- data.frame(distance[[ti]], prop.score = iptw[["psList"]][[ti]][["ps"]][[s]])
    #         }
    #         else {
    #             distance.list[[ti]] <- data.frame(distance[[ti]], prop.score = iptw[["psList"]][[ti]][["ps"]][s])
    #         }
    #     }
    #     
    # }
    # else {
    #     distance.list <- make_list(ntimes)
    #     for (ti in seq_along(distance.list)) {
    #         if (length(s) == 1) {
    #             distance.list[[ti]] <- data.frame(prop.score = iptw[["psList"]][[ti]][["ps"]][[s]])
    #         }
    #         else {
    #             distance.list[[ti]] <- data.frame(prop.score = iptw[["psList"]][[ti]][["ps"]][s])
    #         }
    #     }
    # }
    # if (is_not_null(distance.list)) distance.list <- lapply(distance.list, function(x) get_covs_from_formula(~x))
    # 
    distance.list <- process_distance.list(if_null_then(A[["distance.list"]], A[["distance"]]),
                                           datalist = list(data, ps.data),
                                           covs.list = covs.list, obj.distance = lapply(iptw[["psList"]], function(x) x[["ps"]][,s,drop = FALSE]),
                                           obj.distance.name = if (length(s) > 1) paste.("prop.score", substr(s, 1, nchar(s) - 4)) else "prop.score")
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with iptw objects.", call. = FALSE)
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with iptw objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with iptw objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(iptw, A, treat.list[[1]], covs.list[[1]], method, addl.data = list(data, ps.data), 
                               stop.method = s)
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A[["s.weights"]], iptw[["psList"]][[1]][["sampw"]]))) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, ps.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, ps.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat.list)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, min(lengths(treat.list)))
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("weights"),
                       lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                       imp = imp,
                       original.call.to = "iptw()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat.list)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat.list)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
        }
    }
    
    #Missing values warning
    if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- NULL
    
    #Process output
    X <- initialize_X_msm()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    return(X)
}
x2base.data.frame.list <- function(covs.list, ...) {
    A <- list(...)
    
    #Process iptw
    
    #Process data and get imp
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", datalist = list(data), missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat.list
    treat.list <- process_treat.list(A[["treat.list"]], datalist = list(data))
    
    #Process covs.list
    if (is_null(covs.list)) {
        stop("'covs.list' must be specified.", call. = FALSE)
    }
    if (!is_(covs.list, "list")) {
        stop("'covs.list' must be a list of covariates for which balance is to be assessed at each time point.", call. = FALSE)
    }
    if (any(!vapply(covs.list, is_mat_like, logical(1L)))) {
        stop("Each item in 'covs.list' must be a data frame.", call. = FALSE)
    }
    if (any(vapply(covs.list, function(x) is_null(attr(x, "co.names")), logical(1L)))) {
        covs.list <- lapply(covs.list, function(x) get_covs_from_formula(data = x))
    }
    
    if (length(treat.list) != length(covs.list)) {
        stop("'treat.list' must be a list of treatment statuses at each time point.", call. = FALSE)
    }
    
    #Get estimand
    estimand <- NULL
    
    #Get method
    specified <- setNames(rep(FALSE, 1), "weights")
    if (is_not_null(A[["weights"]])) {
        if (!is_(A[["weights"]], c("character", "numeric", "data.frame", "list"))) {
            stop("The argument to 'weights' must be a vector, list, or data frame of weights or the (quoted) names of variables in 'data' that contain weights.", call. = FALSE)
        }
        specified["weights"] <- TRUE
    }
    
    if (is_null(method <- A[["method"]])) {
        if (specified["weights"]) {
            method <- "weighting"
        }
        else {
            method <- "matching"
        }
    }
    else if (length(method) == 1) {
        specified.method <- match_arg(method, c("weighting", "matching", "subclassification"))
        if (specified.method == "weighting") {
            if (specified["weights"]) {
                method <- "weighting"
            }
            else {
                method <- "matching"
            }
        }
        else {
            if (specified["weights"]) {
                warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
                method <- "matching"
            }
            else {
                method <- "matching"
            }
        }
    }
    else {
        specified.method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
        if (any(specified.method == "subclassification") || specified["subclass"]) {
            warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
            method <- "matching"
        }
        else if (specified["match.strata"]) {
            warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
            method <- "matching"
        }
        else if (!specified["weights"]) {
            warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
            method <- "matching"
        }
        else {
            #Matching and/or weighting with various weights
            method <- specified.method
        }
    }
    
    #Process addl.list 
    addl.list <- process_addl.list(if_null_then(A[["addl.list"]], A[["addl"]]),
                                   datalist = list(data),
                                   covs.list = covs.list)
    
    #Process distance
    distance.list <- process_distance.list(if_null_then(A[["distance.list"]], A[["distance"]]),
                                           datalist = list(data),
                                           covs.list = covs.list)
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with longitudinal treatments.", call. = FALSE)
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with longitudinal treatments.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with longitudinal treatments.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- A[["weights"]])) {
        weights <- process_weights(NULL, A, treat.list[[1]], covs.list[[1]], method, addl.data = list(data))
        method <- attr(weights, "method")
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["s.weights"]])) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat.list)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, min(lengths(treat.list)))
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("weights"),
                       lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                       imp = imp)
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat.list)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat.list)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom("pooled", estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
        }
        else if ("correlations" %in% stats){
            s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
        }
    }
    
    #Missing values warning
    if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- NULL
    
    #Process output
    X <- initialize_X_msm()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "msm"
    
    return(X)
}
x2base.formula.list <- function(formula.list, ...) {
    A <- list(...)
    A[["covs.list"]] <- NULL
    A[["treat.list"]] <- NULL
    
    treat.list <- covs.list <- make_list(length(formula.list))
    for (i in seq_along(formula.list)) {
        treat.list[[i]] <- get_treat_from_formula(formula.list[[i]], A[["data"]])
        covs.list[[i]] <- get_covs_from_formula(formula.list[[i]], A[["data"]])
        names(treat.list)[i] <- attr(treat.list[[i]], "treat.name")
    }
    
    X <- do.call("x2base.data.frame.list", c(list(covs.list, treat.list = treat.list), A))
    return(X)
}
x2base.CBMSM <- function(cbmsm, ...) {
    A <- list(...)
    
    #Process CBMSM
    ID <- sort(unique(cbmsm[["id"]]))
    times <- sort(unique(cbmsm[["time"]]))
    cbmsm[["data"]] <- cbmsm[["data"]][order(cbmsm[["id"]], cbmsm[["time"]]),,drop = FALSE]
    
    #Process data and get imp
    cbmsm.data <- cbmsm[["data"]][cbmsm[["time"]] == 1, , drop = FALSE]
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", datalist = list(data), missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat.list
    treat.list <- process_treat.list(lapply(times, function(x) cbmsm[["treat.hist"]][ID, x]), 
                                     datalist = list(data, cbmsm.data))
    
    #Process covs.list
    covs.list <- make_list(times)
    for (i in seq_along(times)) {
        ti <- times[i]
        cov_i <- get_covs_from_formula(cbmsm[["formula"]], data = cbmsm[["data"]][cbmsm[["time"]] == ti, , drop = FALSE])
        for (co in seq_along(attr(cov_i, "co.names"))) {
            attr(cov_i, "co.names")[[co]][["component"]][attr(cov_i, "co.names")[[co]][["type"]] == "base"] <-
                paste0(attr(cov_i, "co.names")[[co]][["component"]][attr(cov_i, "co.names")[[co]][["type"]] == "base"], "_T", ti)
        }
        names(attr(cov_i, "co.names")) <- vapply(attr(cov_i, "co.names"), function(x) paste0(x[["component"]], collapse = ""), character(1L))
        colnames(cov_i) <- names(attr(cov_i, "co.names"))
        if (i == 1) {
            covs.list[[i]] <- cov_i
        }
        else {
            covs.list[[i]] <- co.cbind(covs.list[[i-1]], cov_i)
        }
    }
    
    #Get estimand
    estimand <- NULL
    
    #Get method
    method <- "weighting"
    
    #Process addl.list 
    addl.list <- process_addl.list(if_null_then(A[["addl.list"]], A[["addl"]]),
                                   datalist = list(data, cbmsm.data),
                                   covs.list = covs.list)
    
    #Process distance
    distance.list <- process_distance.list(if_null_then(A[["distance.list"]], A[["distance"]]),
                                           datalist = list(data, cbmsm.data),
                                           covs.list = covs.list, obj.distance = cbmsm[["fitted.values"]],
                                           obj.distance.name = "prop.score")
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with CBMSM objects.", call. = FALSE)
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with CBMSM objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with CBMSM objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(cbmsm, A, treat.list[[1]], covs.list[[1]], method, addl.data = list(data, cbmsm.data))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- A[["s.weights"]])) {
        stop("Sampling weights are not allowed with CBMSM objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, cbmsm.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat.list)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, min(lengths(treat.list)))
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("weights"),
                       lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                       imp = imp,
                       original.call.to = "CBMSM()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat.list)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat.list)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom("pooled", estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
        }
        else if ("correlations" %in% stats) {
            s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
        }
    }
    
    #Missing values warning
    if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- cbmsm[["call"]]
    
    #Process output
    X <- initialize_X_msm()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    return(X)
}
x2base.weightitMSM <- function(weightitMSM, ...) {
    A <- list(...)
    
    #Process weightitMSM
    
    #Process data and get imp
    weightitMSM.data <- weightitMSM[["data"]]
    d.e.in.w <- vapply(c("exact", "by", "moderator"), function(x) is_not_null(weightitMSM[[x]]), logical(1L))
    if (any(d.e.in.w)) weightitMSM.data2 <- do.call("data.frame", unname(weightitMSM[c("exact", "by", "moderator")[d.e.in.w]]))
    else weightitMSM.data2 <- NULL
    
    imp <- A[["imp"]]
    if (is_not_null(data <- A[["data"]])) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is.data.frame(data))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", datalist = list(data), missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat.list
    treat.list <- process_treat.list(weightitMSM[["treat.list"]],
                                     datalist = list(data, weightitMSM.data, weightitMSM.data2))    
    #Process covs.list
    covs.list <- lapply(weightitMSM[["covs.list"]], function(x) get_covs_from_formula(data = x))
    
    #Get estimand
    estimand <- weightitMSM[["estimand"]]
    
    #Get method
    method <- "weighting"
    
    #Process addl.list 
    addl.list <- process_addl.list(if_null_then(A[["addl.list"]], A[["addl"]]), 
                                   datalist = list(data, weightitMSM.data,
                                                   weightitMSM.data2),
                                   covs.list = covs.list)
    
    #Process distance
    # ntimes <- length(covs.list)
    # distance.list <- list.process("distance.list", A[["distance.list"]], ntimes, 
    #                               "the original call to weightitMSM()",
    #                               treat.list,
    #                               covs.list,
    #                               list(data, weightitMSM.data,
    #                                    weightitMSM.data2))
    # if (is_not_null(distance.list)) distance.list <- lapply(seq_along(distance.list), function(x) data.frame(distance.list[[x]], prop.score = weightitMSM[["ps.list"]][[x]]))
    # else if (is_not_null(weightitMSM[["ps.list"]])) distance.list <- lapply(seq_along(weightitMSM[["ps.list"]]), function(x) data.frame(prop.score = weightitMSM[["ps.list"]][[x]]))
    # else distance.list <- NULL
    # if (is_not_null(distance.list)) distance.list <- lapply(distance.list, function(x) get_covs_from_formula(~x))
    distance.list <- process_distance.list(if_null_then(A[["distance.list"]], A[["distance"]]),
                                           datalist = list(data, weightitMSM.data, weightitMSM.data2),
                                           covs.list = covs.list, obj.distance = weightitMSM[["ps.list"]],
                                           obj.distance.name = "prop.score")
    
    #Process focal
    if (is_not_null(focal <- A[["focal"]])) {
        stop("'focal' is not allowed with weightitMSM objects.", call. = FALSE)
    }
    
    #Process subclass
    if (is_not_null(subclass <- A[["subclass"]])) {
        stop("subclasses are not allowed with weightitMSM objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A[["match.strata"]])) {
        stop("Matching strata are not allowed with weightitMSM objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- process_weights(weightitMSM, A, treat.list[[1]], covs.list[[1]], method, 
                               addl.data = list(data, weightitMSM.data, weightitMSM.data2))
    method <- attr(weights, "method")
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A[["s.weights"]], weightitMSM[["s.weights"]]))) {
        s.weights <- vector.process(s.weights, 
                                    datalist = list(data, weightitMSM.data, weightitMSM.data2),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A[["cluster"]])) {
        cluster <- vector.process(cluster, 
                                  datalist = list(data, weightitMSM.data, weightitMSM.data2),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
        cluster.check(cluster, treat.list)
    }
    
    #Process subset
    if (is_not_null(subset <- A[["subset"]])) {
        subset <- process_subset(subset, min(lengths(treat.list)))
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("weights"),
                       lists = c("treat.list", "covs.list", "addl.list", "distance.list"),
                       imp = imp,
                       original.call.to = "weightitMSM()")
    
    #Process stats and thresholds
    if (!check_if_call_from_fun(bal.plot)) {
        stats <- process_stats(A[["stats"]], treat = treat.list)
        type <- attr(stats, "type")
        
        if (is_not_null(thresholds <- A[["thresholds"]])) {
            thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
            if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
        }
        else thresholds <- list()
        
        for (s in all_STATS(type)) {
            #If disp.stat is TRUE, add stat to stats
            if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- unique(c(stats, s))
            }
            else if (isFALSE(A[[STATS[[s]][["disp_stat"]]]])) {
                stats <- setdiff(stats, s)
            }
            
            #Process and check thresholds
            if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
            }
            if (is_not_null(thresholds[[s]])) {
                thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                    thresholds[[s]] <- NULL
                    warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                   "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                }
                else stats <- unique(c(stats, s))
            }
        }
        
        stats <- process_stats(stats, treat = treat.list)
        
        #Get s.d.denom
        if ("mean.diffs" %in% stats) {
            s.d.denom <- get.s.d.denom("pooled", estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
        }
        else if ("correlations" %in% stats) {
            s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
        }
    }
    
    #Missing values warning
    if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- weightitMSM[["call"]]
    
    #Process output
    X <- initialize_X_msm()
    X.names <- names(X)
    
    for (i in X.names) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    return(X)
}

x2base.default <- function(obj, ...) {
    
    A <- list(...)
    
    if (is_not_null(A) && (is_null(names(A)) || "" %in% names(A))) {
        stop("All arguments to '...' must be named.", call. = FALSE)
    }
    
    if (!is.list(obj)) stop("The input object must be an appropriate list, data.frame, formula, or the output of one of the supported packages.", call. = FALSE)
    
    Q <- list(treat = list(name = c("treat", "tr"), 
                           type = c("numeric", "character", "factor", "logical")),
              treat.list = list(name = c("treat.list", "treat", "tr"),
                                type = c("list", "data.frame")),
              covs = list(name = c("covs", "covariates"), 
                          type = c("data.frame")),
              covs.list = list(name = c("covs.list", "covs", "covariates"),
                               type = c("list")),
              formula = list(name = c("formula", "form"), 
                             type = c("formula")),
              formula.list = list(name = c("formula.list", "formula", "form"),
                                  type = c("list")),
              data = list(name = c("data"),
                          type = c("data.frame", "mids")),
              weights = list(name = c("weights", "w", "wts"),
                             type = c("data.frame", "matrix", "numeric")),
              distance = list(name = c("distance", "distance.list", "ps", "pscore", "p.score", "propensity.score"),
                              type = c("data.frame", "matrix", "numeric", "list")),
              subclass = list(name = c("subclass", "strata"),
                              type = c("factor", "character", "numeric")),
              match.strata = list(name = c("match.strata"),
                                  type = c("factor", "character", "numeric")),
              estimand = list(name = c("estimand", "target", "att", "ate"),
                              type = c("character", "numeric")),
              s.weights = list(name = c("s.weights", "sw", "sweights", "sampw"),
                               type = c("numeric")),
              focal = list(name = c("focal", "treatATT"), 
                           type = c("character", "numeric")),
              call = list(name = c("call"),
                          type = c("call")))
    
    P <- make_list(names(Q))
    names(obj) <- tolower(names(obj))
    
    for (i in names(Q)) {
        if (i %nin% names(A)) {
            for (j in Q[[i]][["name"]]) {
                if (is_null(P[[i]])) {
                    if (is_not_null(obj[[j]])) {
                        if (any(which.type <- vapply(Q[[i]][["type"]], function(x) is_(obj[[j]], x), logical(1L)))) {
                            P[[i]] <- obj[[j]]
                            attr(P[[i]], "name") <- j
                            attr(P[[i]], "type") <- Q[[i]][["type"]][which.type]
                        }
                    }
                }
            }
            if (is_not_null(P[[i]])) {
                assign(i, P[[i]])
                A[[i]] <- P[[i]]
            }
        }
    }
    
    msm <- FALSE
    
    #treat OK
    
    #treat.list
    if (is_not_null(A[["treat.list"]])) {
        if (!all(sapply(A[["treat.list"]], function(x) any(vapply(Q[["treat"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
            A[["treat.list"]] <- NULL
        }
        else msm <- TRUE
    }
    
    #covs 
    if (is_not_null(A[["covs"]])) A[["covs"]] <- as.data.frame(A[["covs"]])
    
    #covs.list
    if (is_not_null(A[["covs.list"]])) {
        if (!all(sapply(A[["covs.list"]], function(x) any(vapply(Q[["covs"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
            A[["covs.list"]] <- NULL
        }
        else msm <- TRUE
    }
    
    #formula OK
    
    #formula.list
    if (is_not_null(A[["formula.list"]])) {
        if (!all(sapply(A[["formula.list"]], function(x) any(vapply(Q[["formula"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
            A[["formula.list"]] <- NULL
        }
        else msm <- TRUE
    }
    
    #data
    if (is_not_null(A[["data"]])) {
        if (is_(A[["data"]], "mids")) {
            A[["data"]] <- imp.complete(A[["data"]])
            if ("imp" %nin% names(A)) A[["imp"]] <- A[["data"]][[".imp"]]
        }
        A[["data"]] <- as.data.frame(A[["data"]])
    }
    
    #weights
    if (is_not_null(A[["weights"]])) {
        # if (is.vector(A[["weights"]], "numeric")) A[["weights"]] <- data.frame(weights = A[["weights"]])
        # else A[["weights"]] <- as.data.frame(A[["weights"]])
    }
    
    #distance
    if (is_not_null(A[["distance"]])) {
        if (is.list(A[["distance"]]) && !is.data.frame(A[["distance"]])) {
            if (!all(sapply(A[["distance"]], function(x) any(vapply(Q[["distance"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
                A[["distance"]] <- NULL
            }
        }
        else if (is.numeric(A[["distance"]])) {
            if (is_not_null(attr(A[["distance"]], "name"))) A[["distance"]] <- setNames(data.frame(A[["distance"]]),
                                                                                        attr(A[["distance"]], "name"))
            else A[["distance"]] <- data.frame(distance = A[["distance"]])
        }
        else A[["distance"]] <- as.data.frame(A[["distance"]])
    }
    
    #distance.list
    if (is_not_null(A[["distance.list"]])) {
        if (!all(sapply(A[["distance.list"]], function(x) any(vapply(Q[["distance"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
            A[["distance.list"]] <- NULL
        }
        #msm <- TRUE
    }
    
    #subclass
    if (is_not_null(A[["subclass"]])) A[["subclass"]] <- factor(A[["subclass"]])
    
    #match.strata
    if (is_not_null(A[["match.strata"]])) A[["match.strata"]] <- factor(A[["match.strata"]])
    
    #estimand
    if (is_not_null(A[["estimand"]])) {
        estimand.name <- attr(A[["estimand"]], "name")
        if (is_not_null(estimand.name) && toupper(estimand.name) == "ATT") {
            if (A[["estimand"]] == 0) A[["estimand"]] <- "ATE"
            else A[["estimand"]] <- "ATT"
        }
        else if (is_not_null(estimand.name) && toupper(estimand.name) == "ATE") {
            if (A[["estimand"]] == 0) A[["estimand"]] <- "ATT"
            else A[["estimand"]] <- "ATE"
        }
        else {
            if (tolower(A[["estimand"]]) %in% c("att", "treat", "treated", "tr", "t", "atet")) A[["estimand"]] <- "ATT"
            else if (tolower(A[["estimand"]]) %in% c("ate", "all")) A[["estimand"]] <- "ATE"
            else if (tolower(A[["estimand"]]) %in% c("atc", "control", "untreated", "u", "c", "ctrl", "atu", "atec", "ateu")) A[["estimand"]] <- "ATC"
            else A[["estimand"]] <- NULL
        }
    }
    
    #s.weights OK
    
    #focal OK
    
    #call OK
    
    #model (only to extract data)
    if (is_not_null(obj[["model"]])) {
        if (is_null(A[["data"]]) && "data" %in% names(obj[["model"]])) {
            A[["data"]] <- obj[["model"]][["data"]]
        }
    }
    
    if (!msm) {
        
        #Process data and get imp
        imp <- A[["imp"]]
        if (is_not_null(data <- A[["data"]])) {
            if (is_(data, "mids")) {
                data <- imp.complete(data)
                if (is_null(imp)) imp <- data[[".imp"]]
            }
            else if (!is.data.frame(data))
            {
                # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
                data <- NULL
            }
        }
        
        #Process imp
        if (is_not_null(imp)) {
            imp <- vector.process(imp, 
                                  name = "imp", 
                                  which = "imputation identifiers", 
                                  datalist = list(data), 
                                  missing.okay = FALSE)
            imp <- factor(imp)
        }
        
        #Process treat
        t.c <- use.tc.fd(A[["formula"]], data, A[["treat"]], A[["covs"]])
        treat <- process_treat(t.c[["treat"]], datalist = list(data))
        
        #Process covs
        covs <- t.c[["covs"]]
        if (is_null(covs)) {
            stop("Covariates must be specified using 'covs' or 'formula'.", call. = FALSE)
        }
        
        #Get estimand
        estimand <- A[["estimand"]]
        
        #Get method
        specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
        if (is_not_null(A[["weights"]])) {
            specified["weights"] <- TRUE
        }
        if (is_not_null(A[["subclass"]])){
            specified["subclass"] <- TRUE
        }
        if (is_not_null(A[["match.strata"]])) {
            specified["match.strata"] <- TRUE
        }
        
        if (is_null(method <- A[["method"]])) {
            if (specified["match.strata"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Assuming \"matching\" and using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                    A[["weights"]] <- A[["subclass"]] <- NULL
                }
                method <- "matching"
            }
            else if (specified["subclass"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    A[["weights"]] <- A[["match.strata"]] <- NULL
                }
                method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    A[["match.strata"]] <- A[["subclass"]] <- NULL
                }
                method <- "weighting"
            }
            else {
                method <- "matching"
            }
        }
        else if (length(method) == 1) {
            specified.method <- match_arg(method, c("weighting", "matching", "subclassification"))
            if (specified.method == "weighting") {
                if (specified["weights"]) {
                    if (sum(specified) > 1) {
                        message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
                        A[["match.strata"]] <- A[["subclass"]] <- NULL
                    }
                    method <- "weighting"
                }
                else if (specified["match.strata"]) {
                    message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
                    A[["subclass"]] <- NULL
                    method <- "matching"
                }
                else if (specified["subclass"]) {
                    message("method = \"weighting\" is specified, but no weights are present. Assuming \"subclassification\" and using subclass instead.")
                    method <- "subclassification"
                }
                else {
                    method <- "matching"
                }
            }
            else if (specified.method == "matching") {
                if (specified["match.strata"]) {
                    if (sum(specified) > 1) {
                        message(word_list(names(specified)[specified]), " are specified. Using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                        A[["weights"]] <- A[["subclass"]] <- NULL
                    }
                    method <- "matching"
                }
                else if (specified["weights"]) {
                    if (sum(specified) > 1) {
                        message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
                        A[["match.strata"]] <- A[["subclass"]] <- NULL
                    }
                    method <- "matching"
                }
                else if (specified["subclass"]) {
                    message("method = \"matching\" is specified, but no weights or match.strata are present. Assuming \"subclassification\" and using subclass instead.")
                    method <- "subclassification"
                }
                else {
                    method <- "matching"
                }
            }
            else if (specified.method == "subclassification") {
                if (specified["subclass"]) {
                    if (sum(specified) > 1) {
                        message(word_list(names(specified)[specified]), " are specified. Using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                        A[["weights"]] <- A[["match.strata"]] <- NULL
                    }
                    method <- "subclassification"
                }
                else if (specified["match.strata"]) {
                    message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
                    A[["weights"]] <- NULL
                    method <- "matching"
                }
                else if (specified["weights"]) {
                    message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"weighting\" and using weights instead.")
                    method <- "weighting"
                }
            }
        }
        else {
            specified.method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
            if (any(specified.method == "subclassification") || specified["subclass"]) {
                stop("Subclassification cannot be specified along with other methods.", call. = FALSE)
            }
            else if (specified["match.strata"]) {
                stop("Only weights can be specified with multiple methods.", call. = FALSE)
            }
            else if (!specified["weights"]) {
                warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
                method <- "matching"
            }
            else {
                #Matching and/or weighting with various weights
                method <- specified.method
                A[["match.strata"]] <- A[["subclass"]] <- NULL
            }
        }
        
        #Process addl 
        addl <- process_addl(A[["addl"]], datalist = list(data, covs))
        
        #Process distance
        distance <- process_distance(A[["distance"]], datalist = list(data, covs))
        
        #Process subclass
        if (is_not_null(subclass <- A[["subclass"]])) {
            subclass <- vector.process(subclass, 
                                       datalist = list(data),
                                       name = "subclass", 
                                       which = "subclass membership",
                                       missing.okay = TRUE)
            subclass <- factor(subclass)
        }
        
        #Process match.strata
        else if (is_not_null(match.strata <- A[["match.strata"]])) {
            match.strata <- vector.process(match.strata, 
                                           datalist = list(data),
                                           name = "match.strata", 
                                           which = "matching strata membership",
                                           missing.okay = TRUE)
            weights <- data.frame(weights = strata2weights(match.strata,
                                                           treat = treat,
                                                           estimand = estimand,
                                                           focal = focal))
        }
        
        #Process weights
        else if (is_not_null(A[["weights"]])) {
            weights <- process_weights(NULL, A, treat, covs, method, addl.data = list(data))
            method <- attr(weights, "method")
        }
        else {
            weights <- NULL
        }
        
        #Process s.weights
        if (is_not_null(s.weights <- A[["s.weights"]])) {
            s.weights <- vector.process(s.weights, 
                                        datalist = list(data),
                                        name = "s.weights", 
                                        which = "sampling weights",
                                        missing.okay = FALSE)
            weight.check(s.weights)
        }
        
        #Process cluster
        if (is_not_null(cluster <- A[["cluster"]])) {
            cluster <- vector.process(cluster, 
                                      datalist = list(data),
                                      name = "cluster", 
                                      which = "cluster membership",
                                      missing.okay = FALSE)
            cluster <- factor(cluster)
            cluster.check(cluster, treat)
        }
        
        #Process subset
        if (is_not_null(subset <- A[["subset"]])) {
            subset <- process_subset(subset, length(treat))
        }
        
        #Process discarded
        discarded <- A[["discarded"]]
        
        #Process length
        length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                           data.frames = c("covs", "weights", "distance", "addl"),
                           imp = imp)
        
        #Process focal
        if (is_not_null(focal <- A[["focal"]]) && get.treat.type(treat) != "continuous") {
            focal <- process_focal(focal, treat)
        } else if (get.treat.type(treat) == "binary" && is_not_null(estimand)) {
            focal <- switch(toupper(estimand), 
                            "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                            "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                            NULL)
        }
        
        #Process pairwise
        if (get.treat.type(treat) == "binary" && is_null(focal)) {
            if (is_null(A[["pairwise"]])) A[["pairwise"]] <- TRUE
            if (isFALSE(A[["pairwise"]])) attr(treat, "treat.type") <- "multinomial"
        }
        
        #Process stats and thresholds
        if (!check_if_call_from_fun(bal.plot)) {
            stats <- process_stats(A[["stats"]], treat = treat)
            type <- attr(stats, "type")
            
            if (is_not_null(thresholds <- A[["thresholds"]])) {
                thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
                if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
            }
            else thresholds <- list()
            
            for (s in all_STATS(type)) {
                #If disp.stat is TRUE, add stat to stats
                if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                    stats <- unique(c(stats, s))
                }
                
                #Process and check thresholds
                if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                    thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
                }
                if (is_not_null(thresholds[[s]])) {
                    thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                    if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                        thresholds[[s]] <- NULL
                        warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                       "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                    }
                    else stats <- unique(c(stats, s))
                }
            }
            
            stats <- process_stats(stats, treat = treat)
            
            #Get s.d.denom
            if ("mean.diffs" %in% stats) {
                s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = estimand, 
                                           weights = weights, subclass = subclass, 
                                           treat = treat, focal = focal)
            }
            else if ("correlations" %in% stats) {
                s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
            }
        }
        
        #Missing values warning
        if (anyNA(covs) || anyNA(addl)) {
            warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
        }
        
        #Get call
        call <- A[["call"]]
        
        #Process output
        X <- initialize_X()
        X.names <- names(X)
        
        for (i in X.names) {
            X[[i]] <- get0(i, inherits = FALSE)
        }
        
        X <- subset_X(X, subset)
        X <- setNames(X[X.names], X.names)
        
    }
    else {
        
        #Process input
        initial.list.lengths <- c(length(A[["formula.list"]]), length(A[["covs.list"]]), length(A[["treat.list"]]))
        if (!all_the_same(initial.list.lengths[initial.list.lengths != 0])) stop("The lists in the object were not the same length.", call. = FALSE)
        ntimes.guess <- max(initial.list.lengths)
        if (is_null(A[["treat.list"]])) A[["treat.list"]] <- make_list(length(ntimes.guess)) 
        if (is_null(A[["covs.list"]])) A[["covs.list"]] <- make_list(length(ntimes.guess)) 
        
        #Process data and get imp
        imp <- A[["imp"]]
        if (is_not_null(data <- A[["data"]])) {
            if (is_(data, "mids")) {
                data <- imp.complete(data)
                if (is_null(imp)) imp <- data[[".imp"]]
            }
            else if (!is.data.frame(data))
            {
                # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
                data <- NULL
            }
        }
        
        #Process imp
        if (is_not_null(imp)) {
            imp <- vector.process(imp, "imp", "imputation identifiers", datalist = list(data), missing.okay = FALSE)
            imp <- factor(imp)
        }
        
        #Process treat.list
        for (i in seq_len(ntimes.guess)) {
            t.c <- use.tc.fd(A[["formula.list"]][[i]], data, A[["treat.list"]][[i]], A[["covs.list"]][[i]])
            
            A[["treat.list"]][[i]] <- t.c[["treat"]]
            A[["covs.list"]][[i]]  <- t.c[["covs"]]
            if (is_not_null(t.c[["treat.name"]])) names(A[["treat.list"]])[i] <- t.c[["treat.name"]]
        }
        treat.list <- process_treat.list(A[["treat.list"]], datalist = list(data))
        
        #Process covs.list
        if (is_null(covs.list <- A[["covs.list"]])) {
            stop("'covs.list' must be specified.", call. = FALSE)
        }
        if (!is_(covs.list, "list")) {
            stop("'covs.list' must be a list of covariates for which balance is to be assessed at each time point.", call. = FALSE)
        }
        if (any(!vapply(covs.list, is.data.frame, logical(1L)))) {
            stop("Each item in 'covs.list' must be a data frame.", call. = FALSE)
        }
        if (length(treat.list) != length(covs.list)) {
            stop("'treat.list' must be a list of treatment statuses at each time point.", call. = FALSE)
        }
        
        #Get estimand
        estimand <- NULL
        
        #Get method
        specified <- setNames(rep(FALSE, 1), "weights")
        if (is_not_null(A[["weights"]])) {
            if (!is_(A[["weights"]], c("character", "numeric", "data.frame", "list"))) {
                stop("The argument to 'weights' must be a vector, list, or data frame of weights or the (quoted) names of variables in 'data' that contain weights.", call. = FALSE)
            }
            specified["weights"] <- TRUE
        }
        
        if (is_null(method <- A[["method"]])) {
            if (specified["weights"]) {
                method <- "weighting"
            }
            else {
                method <- "matching"
            }
        }
        else if (length(method) == 1) {
            specified.method <- match_arg(method, c("weighting", "matching", "subclassification"))
            if (specified.method == "weighting") {
                if (specified["weights"]) {
                    method <- "weighting"
                }
                else {
                    method <- "matching"
                }
            }
            else {
                if (specified["weights"]) {
                    warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
                    method <- "matching"
                }
                else {
                    method <- "matching"
                }
            }
        }
        else {
            specified.method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
            if (any(specified.method == "subclassification") || specified["subclass"]) {
                warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
                method <- "matching"
            }
            else if (specified["match.strata"]) {
                warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
                method <- "matching"
            }
            else if (!specified["weights"]) {
                warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
                method <- "matching"
            }
            else {
                #Matching and/or weighting with various weights
                method <- specified.method
            }
        }
        
        #Process addl.list 
        addl.list <- process_addl.list(if_null_then(A[["addl.list"]], A[["addl"]]), 
                                       datalist = list(data),
                                       covs.list = covs.list)
        
        #Process distance
        # ntimes <- length(covs.list)
        # distance.list <- list.process("distance.list", A[["distance.list"]], ntimes, 
        #                               "covs.list",
        #                               treat.list,
        #                               covs.list,
        #                               list(data))
        # if (is_not_null(distance.list)) distance.list <- lapply(distance.list, function(x) get_covs_from_formula(~x))
        distance.list <- process_distance.list(if_null_then(A[["distance.list"]], A[["distance"]]),
                                               datalist = list(data),
                                               covs.list = covs.list)
        #Process focal
        if (is_not_null(focal <- A[["focal"]])) {
            stop("'focal' is not allowed with longitudinal treatments.", call. = FALSE)
        }
        
        #Process subclass
        if (is_not_null(subclass <- A[["subclass"]])) {
            stop("subclasses are not allowed with longitudinal treatments.", call. = FALSE)
        }
        
        #Process match.strata
        if (is_not_null(match.strata <- A[["match.strata"]])) {
            stop("Matching strata are not allowed with longitudinal treatments.", call. = FALSE)
        }
        
        #Process weights
        weights <- process_weights(NULL, A, treat.list[[1]], covs.list[[1]], method, addl.data = list(data))
        method <- attr(weights, "method")
        
        #Process s.weights
        if (is_not_null(s.weights <- A[["s.weights"]])) {
            s.weights <- vector.process(s.weights, 
                                        datalist = list(data),
                                        name = "s.weights", 
                                        which = "sampling weights",
                                        missing.okay = FALSE)
        }
        
        #Process cluster
        if (is_not_null(cluster <- A[["cluster"]])) {
            cluster <- vector.process(cluster, 
                                      datalist = list(data),
                                      name = "cluster", 
                                      which = "cluster membership",
                                      missing.okay = FALSE)
            cluster <- factor(cluster)
            cluster.check(cluster, treat.list)
        }
        
        #Process subset
        if (is_not_null(subset <- A[["subset"]])) {
            subset <- process_subset(subset, min(lengths(treat.list)))
        }
        
        #Process discarded
        
        #Process length
        length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                           data.frames = c("weights"),
                           lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                           imp = imp)
        
        #Process stats and thresholds
        if (!check_if_call_from_fun(bal.plot)) {
            stats <- process_stats(A[["stats"]], treat = treat.list)
            type <- attr(stats, "type")
            
            if (is_not_null(thresholds <- A[["thresholds"]])) {
                thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
                if (any(names(thresholds) %nin% stats)) stats <- unique(c(stats, names(thresholds)))
            }
            else thresholds <- list()
            
            for (s in all_STATS(type)) {
                #If disp.stat is TRUE, add stat to stats
                if (isTRUE(A[[STATS[[s]][["disp_stat"]]]])) {
                    stats <- unique(c(stats, s))
                }
                
                #Process and check thresholds
                if (is_not_null(A[[STATS[[s]][["threshold"]]]])) {
                    thresholds[[s]] <- A[[STATS[[s]][["threshold"]]]]
                }
                if (is_not_null(thresholds[[s]])) {
                    thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])
                    if (!between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
                        thresholds[[s]] <- NULL
                        warning(paste0(STATS[[s]][["threshold"]], " must be between ", word_list(STATS[[s]][["threshold_range"]]),
                                       "; ignoring ", STATS[[s]][["threshold"]], "."), call. = FALSE)
                    }
                    else stats <- unique(c(stats, s))
                }
            }
            
            stats <- process_stats(stats, treat = treat.list)
            
            #Get s.d.denom
            if ("mean.diffs" %in% stats) {
                s.d.denom <- get.s.d.denom("pooled", estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
            }
            else if ("correlations" %in% stats) {
                s.d.denom <- get.s.d.denom.cont(A[["s.d.denom"]], weights = weights, subclass = subclass)
            }
        }
        
        #Missing values warning
        if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
            warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
        }
        
        #Get call
        call <- NULL
        
        #Process output
        X <- initialize_X_msm()
        X.names <- names(X)
        
        for (i in X.names) {
            X[[i]] <- get0(i, inherits = FALSE)
        }
        
        X <- subset_X(X, subset)
        X <- setNames(X[X.names], X.names)
        
        class(X) <- "msm"
    }
    
    return(X)
}
