#Functions to convert object to base.bal.tab input

x2base <- function(...) UseMethod("x2base")

x2base.matchit <- function(m, ...) {
    A <- list(...)
    
    #Process matchit
    
    #Process data and get imp
    m.data <- m$model$data
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data, m.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(m$treat, data = list(data, m.data))
    
    #Process covs
    if (is_not_null(m$model$model)) {
        if (nrow(m$model$model) == length(treat)) {
            covs <- data.frame(m$model$model[, names(m$model$model) %in% attributes(terms(m$model))$term.labels])
        }
        else {
            #Recreating covs from model object and m$X. Have to do this because when 
            #drop != NULL and reestimate = TRUE, cases are lost. This recovers them.
            
            order <- setNames(attr(m$model$terms, "order"),
                              attr(m$model$terms, "term.labels"))
            assign <- setNames(attr(m$X, "assign"), colnames(m$X))
            assign1 <- assign[assign %in% which(order == 1)] #Just main effects
            
            dataClasses <- attr(m$model$terms, "dataClasses")
            factors.to.unsplit <- names(dataClasses)[dataClasses %in% c("factor", "character", "logical")]
            f0 <- setNames(lapply(factors.to.unsplit, 
                                  function(x) {
                                      if (dataClasses[x] == "factor")
                                          list(levels = levels(m$model$model[[x]]),
                                               faclev = paste0(x, levels(m$model$model[[x]])))
                                      else 
                                          list(levels = unique(m$model$model[[x]]),
                                               faclev = paste0(x, unique(m$model$model[[x]])))
                                  }),
                           factors.to.unsplit)
            covs <- as.data.frame(m$X[, names(assign1)])
            for (i in factors.to.unsplit) {
                covs <- unsplitfactor(covs, i, sep = "",
                                      dropped.level = f0[[i]]$levels[f0[[i]]$faclev %nin% colnames(m$X)])
                if (dataClasses[i] == "logical") covs[[i]] <- as.logical(covs[[i]])
            }
        }
        
    }
    else if ("matchit.mahalanobis" %nin% class(m) && is_not_null(data)) {
        t.c <- get.covs.and.treat.from.formula(m$formula, data = data)
        covs <- t.c[["reported.covs"]]
    }
    else {
        covs <- data.frame(m$X)
    }
    
    #Get estimand
    estimand <- A$estimand
    
    #Get method
    if (any(class(m) == "matchit.subclass")) {
        method <- "subclassification"
    }
    else if (any(class(m) == "matchit.full")) {
        method <- "weighting"
    }
    else {
        method <- "matching"
    }
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data, m.data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data, m.data)
    if (any(is.finite(m$distance))) {
        if (is_not_null(distance)) distance <- cbind(distance, distance = m$distance)
        else distance <- data.frame(distance = m$distance)
    }
    
    #Process subclass
    if ("matchit.subclass" %in% class(m)) {
        subclass <- factor(m$subclass)
    }
    else subclass <- NULL
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with matchit objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- data.frame(weights = m$weights)
    weight.check(weights)
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        stop("s.weights are not allowed with matchit objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, m.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    discarded <- m$discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "matchit()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with matchit objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- m$call
    
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
    if (is_not_null(A) && names(A)[1]=="" && is_null(A$stop.method)) A$stop.method <- A[[1]]
    if (is_null(A$stop.method) && is_not_null(A$full.stop.method)) A$stop.method <- A$full.stop.method
    
    if (is_not_null(A$stop.method)) {
        if (is.character(A$stop.method)) {
            rule1 <- names(ps$w)[sapply(t(sapply(tolower(A$stop.method), function(x) startsWith(tolower(names(ps$w)), x))), any)]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- names(ps$w)
            }
        }
        else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(names(ps$w)))) {
            if (any(!A$stop.method %in% seq_along(names(ps$w)))) {
                message(paste0("Warning: There are ", length(names(ps$w)), " stop methods available, but you requested ", 
                               word_list(A$stop.method[!A$stop.method %in% seq_along(names(ps$w))], and.or = "and"),"."))
            }
            rule1 <- names(ps$w)[A$stop.method %in% seq_along(names(ps$w))]
        }
        else {
            warning("stop.method should be ", word_list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- names(ps$w)
        }
    }
    else {
        rule1 <- names(ps$w)
    }
    
    s <- names(ps$w)[match(tolower(rule1), tolower(names(ps$w)))]

    #Process data and get imp
    ps.data <- ps$data
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data, ps.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(ps$treat, data = list(data, ps.data))
    
    #Process covs
    covs <- ps.data[, ps$gbm.obj$var.names, drop = FALSE]
    
    #Get estimand
    estimand <- ps$estimand
    
    #Get method
    method <- rep("weighting", length(s))
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data, ps.data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data, ps.data)
    if (is_not_null(distance)) {
        if (length(s) == 1) {
            distance <- cbind(distance, prop.score = ps$ps[[s]])
        }
        else {
            distance <- cbind(distance, prop.score = ps$ps[s])
        }
    }
    else {
        if (length(s) == 1) {
            distance <- data.frame(prop.score = ps$ps[[s]])
        }
        else {
            distance <- data.frame(prop.score = ps$ps[s])
        }
    }
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with ps objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with ps objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(get.w(ps, s, estimand)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A$s.weights, ps$sampw))) {
        s.weights <- vector.process(s.weights, 
                                    data = list(data, ps.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }

    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, ps.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "ps()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with ps objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)

    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- ps$parameters
    
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
    if (is_not_null(A)&& names(A)[1]=="" && is_null(A$stop.method)) A$stop.method <- A[[1]]
    if (is_null(A$stop.method) == 0 && is_not_null(A$full.stop.method)) A$stop.method <- A$full.stop.method
    
    if (is_not_null(A$stop.method)) {
        if (any(is.character(A$stop.method))) {
            rule1 <- mnps$stopMethods[sapply(t(sapply(tolower(A$stop.method), function(x) startsWith(tolower(mnps$stopMethods), x))), any)]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- mnps$stopMethods
            }
        }
        else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(mnps$stopMethods))) {
            if (any(!A$stop.method %in% seq_along(mnps$stopMethods))) {
                message(paste0("Warning: There are ", length(mnps$stopMethods), " stop methods available, but you requested ", 
                               word_list(A$stop.method[!A$stop.method %in% seq_along(mnps$stopMethods)], and.or = "and"),"."))
            }
            rule1 <- mnps$stopMethods[A$stop.method %in% seq_along(mnps$stopMethods)]
        }
        else {
            warning("stop.method should be ", word_list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- mnps$stopMethods
        }
    }
    else {
        rule1 <- mnps$stopMethods
    }
    
    s <- mnps$stopMethods[match(tolower(rule1), tolower(mnps$stopMethods))]
    
    #Process data and get imp
    mnps.data <- mnps$data
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data, mnps.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(mnps$treatVar, data = list(data, mnps.data))
    
    #Process covs
    covs <- mnps$data[mnps$psList[[1]]$gbm.obj$var.names]
    
    #Get estimand
    estimand <- mnps$estimand
    
    #Get method
    method <- rep("weighting", length(s))
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data, mnps.data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data, mnps.data)
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with mnps objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with mnps objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(get.w(mnps, s)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A$s.weights, mnps$sampw))) {
        s.weights <- vector.process(s.weights, 
                                    data = list(data, mnps.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, mnps.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "mnps()")
    
    #Process focal
    focal <- mnps$treatATT
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
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
    
    #Process ps
    if (is_not_null(A) && names(A)[1]=="" && is_null(A$stop.method)) A$stop.method <- A[[1]]
    if (is_null(A$stop.method) && is_not_null(A$full.stop.method)) A$stop.method <- A$full.stop.method
    
    if (is_not_null(A$stop.method)) {
        if (is.character(A$stop.method)) {
            rule1 <- names(ps.cont$w)[sapply(t(sapply(tolower(A$stop.method), function(x) startsWith(tolower(names(ps.cont$w)), x))), any)]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(names(ps.cont$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- names(ps.cont$w)
            }
        }
        else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(names(ps.cont$w)))) {
            if (any(!A$stop.method %in% seq_along(names(ps.cont$w)))) {
                message(paste0("Warning: There are ", length(names(ps.cont$w)), " stop methods available, but you requested ", 
                               word_list(A$stop.method[!A$stop.method %in% seq_along(names(ps.cont$w))], and.or = "and"),"."))
            }
            rule1 <- names(ps.cont$w)[A$stop.method %in% seq_along(names(ps.cont$w))]
        }
        else {
            warning("stop.method should be ", word_list(names(ps.cont$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- names(ps.cont$w)
        }
    }
    else {
        rule1 <- names(ps.cont$w)
    }
    
    s <- names(ps.cont$w)[match(tolower(rule1), tolower(names(ps.cont$w)))]
    
    #Process data and get imp
    ps.data <- ps.cont$data
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data, ps.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(ps.cont$treat, data = list(data, ps.data))
    
    #Process covs
    covs <- ps.cont$data[, ps.cont$gbm.obj$var.names, drop = FALSE]
    
    #Get estimand

    #Get method
    method <- rep("weighting", length(s))
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data, ps.data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data, ps.data)
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with ps.cont objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with ps.cont objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(get.w(ps.cont, s)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A$s.weights, ps.cont$sampw))) {
        s.weights <- vector.process(s.weights, 
                                    data = list(data, ps.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, ps.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "ps.cont()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with ps.cont objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- ps.cont$parameters
    
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
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    treat <- process_treat(t.c[["treat"]], data = data)
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- Match$estimand
    
    #Get method
    method <- "matching"
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data)
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with Match objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with Match objects.", call. = FALSE)
    }
    
    #Process weights
    weights <- data.frame(weights = get.w.Match(Match))
    weight.check(weights)
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        stop("s.weights are not allowed with Match objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = data,
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    discarded <- rep(FALSE, length(treat))
    if (is_not_null(Match$index.dropped)) discarded[Match$index.dropped] <- TRUE
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "Match()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with Match objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
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
    
    if ("data" %in% names(A) && is_(A[["data"]], "mids")) {
        A[["data"]] <- imp.complete(A[["data"]])
        if (is_null(A[["imp"]])) A[["imp"]] <- A[["data"]][[".imp"]]
    }
    
    t.c <- get.covs.and.treat.from.formula(formula, A[["data"]], treat = A[["treat"]])
    covs <- t.c[["reported.covs"]]
    treat <- t.c[["treat"]]
    
    if (is_null(covs)) stop("The right hand side of the formula must contain covariates for which balance is to be assessed.", call. = FALSE)
    
    A[["covs"]] <- NULL
    A[["treat"]] <- NULL
    
    X <- do.call(x2base.data.frame, c(list(covs = covs, treat = treat), A))
    return(X)
}
x2base.data.frame <- function(covs, ...) {
    A <- list(...)
    
    #Process data.frame
    
    #Process data and get imp
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(A$treat, data = data)
    
    #Process covs
    if (is_null(covs)) {
        stop("covs data.frame must be specified.", call. = FALSE)
    }
    is_(covs, "data.frame", stop = TRUE)
    
    #Get estimand
    estimand <- A$estimand
    
    #Get method
    specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
    if (is_not_null(A$weights)) {
        specified["weights"] <- TRUE
    }
    if (is_not_null(A$subclass)){
        specified["subclass"] <- TRUE
    }
    if (is_not_null(A$match.strata)) {
        specified["match.strata"] <- TRUE
    }
    
    if (is_null(A$method)) {
        if (specified["match.strata"]) {
            if (sum(specified) > 1) {
                message(word_list(names(specified)[specified]), " are specified. Assuming \"matching\" and using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                A$weights <- A$subclass <- NULL
            }
            method <- "matching"
        }
        else if (specified["subclass"]) {
            if (sum(specified) > 1) {
                message(word_list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                A$weights <- A$match.strata <- NULL
            }
            method <- "subclassification"
            #weights <- rep(1, nrow(covs))
        }
        else if (specified["weights"]) {
            if (sum(specified) > 1) {
                message(word_list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                A$match.strata <- A$subclass <- NULL
            }
            method <- "weighting"
        }
        else {
            method <- "matching"
        }
    }
    else if (length(A$method) == 1) {
        specified.method <- match_arg(A$method, c("weighting", "matching", "subclassification"))
        if (specified.method == "weighting") {
            if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
                    A$match.strata <- A$subclass <- NULL
                }
                method <- "weighting"
            }
            else if (specified["match.strata"]) {
                message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
                A$subclass <- NULL
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
                    A$weights <- A$subclass <- NULL
                }
                method <- "matching"
            }
            else if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
                    A$match.strata <- A$subclass <- NULL
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
                    A$weights <- A$match.strata <- NULL
                }
                method <- "subclassification"
            }
            else if (specified["match.strata"]) {
                message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
                A$weights <- NULL
                method <- "matching"
            }
            else if (specified["weights"]) {
                message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"weighting\" and using weights instead.")
                method <- "weighting"
            }
        }
    }
    else {
        specified.method <- match_arg(A$method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
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
            A$match.strata <- A$subclass <- NULL
        }
    }
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data)

    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        subclass <- vector.process(subclass, 
                                   data = data,
                                   name = "subclass", 
                                   which = "subclass membership",
                                   missing.okay = TRUE)
        subclass <- factor(subclass)
        weights <- NULL
    }
    
    #Process match.strata
    else if (is_not_null(match.strata <- A$match.strata)) {
        match.strata <- vector.process(match.strata, 
                                       data = data,
                                       name = "match.strata", 
                                       which = "matching strata membership",
                                       missing.okay = TRUE)
        weights <- data.frame(weights = strata2weights(match.strata,
                                                       treat = treat))
    }
    #Process weights
    else if (is_not_null(weights <- data.frame.process("weights", A[["weights"]], treat, covs, data))) {
        
        weight.check(weights)
        
        if (length(method) == 1) {
            method <- rep.int(method, ncol(weights))
        }
        else if (length(method) != ncol(weights)) {
            stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
        }
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        s.weights <- vector.process(s.weights, 
                                    data = data,
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = data,
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    discarded <- A$discarded
    
    #Process length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp)
    
    #Process focal
    if (is_not_null(focal <- A$focal) && get.treat.type(treat) != "continuous") {
        if (is.numeric(focal)) {
            if (focal <= nunique(treat)) focal <- levels(treat)[focal]
            else 
                stop(paste0("focal was specified as ", focal, 
                            ", but there are only ", nunique(treat), " treatment groups."), call. = FALSE)
        }
        else {
            if (focal %nin% levels(treat)) 
                stop(paste0("The name specified to focal is not the name of any treatment group."), call. = FALSE)
        }
    }
    
    #Get s.d.denom
    if (get.treat.type(treat) != "continuous") {
        s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, 
                                   weights = weights, subclass = subclass, 
                                   treat = treat, focal = focal)
    }
    else {
        s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
    }
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
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
    c.data <- cbps.fit$data
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data, c.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    t.c <- use.tc.fd(cbps.fit$formula, c.data)
    treat <- process_treat(t.c[["treat"]], data = list(data, c.data))
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- A$estimand
    
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data, c.data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data, c.data)
    if (get.treat.type(treat) == "binary") {
        if (all_the_same(cbps.fit$fitted.values)) {
            if (is_null(distance)) distance <- NULL
        }
        else {
            if (is_null(distance)) distance <- data.frame(prop.score = cbps.fit$fitted.values)
            else distance <- cbind(distance, prop.score = cbps.fit$fitted.values)
        }
    }
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with CBPS objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with CBPS objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w(cbps.fit, use.weights = A$use.weights)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$sampw)) {
        s.weights <- vector.process(s.weights, 
                                  data = list(data, c.data),
                                  name = "s.weights", 
                                  which = "sampling weights",
                                  missing.okay = FALSE)
        weight.check(s.weights)
        weights <- weights/s.weights #Because CBPS weights contain s.weights in them
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, c.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "CBPS()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with CBPS objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    if (get.treat.type(treat) != "continuous") {
        s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    }
    else {
        s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
    }
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- cbps.fit$call
    
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
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    treat <- process_treat(t.c[["treat"]], data = data)
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- "ATT"
    
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data)
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with ebalance objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with ebalance objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w(ebalance, treat)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$sampw)) {
        s.weights <- vector.process(s.weights, 
                                    data = data,
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = data,
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "ebalance()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with ebalance objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
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
x2base.ebalance.trim <- x2base.ebalance
x2base.optmatch <- function(optmatch, ...) {
    A <- list(...)
    
    #Process optmatch
    if (all(is.na(optmatch))) stop("'optmatch' object contains no valid matches")

    #Process data and get imp
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    treat <- process_treat(t.c[["treat"]], data = data)
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- A$estimand
    
    #Get method
    method <- "matching"
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data)
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with optmatch objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with optmatch objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w.optmatch(optmatch)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        stop("s.weights are not allowed with optmatch objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = data,
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = paste0(deparse(attr(optmatch, "call")[[1]]), "()"))
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with optmatch objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
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
x2base.weightit <- function(weightit, ...) {
    A <- list(...)
    
    #Process CBPS
    
    #Process data and get imp
    d.e.in.w <- vapply(c("covs", "exact", "by", "moderator"), function(x) is_not_null(weightit[[x]]), logical(1L))
    if (any(d.e.in.w)) weightit.data <- do.call("cbind", unname(weightit[c("covs", "exact", "by", "moderator")[d.e.in.w]]))
    else weightit.data <- NULL
    
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
        {
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, 
                              name = "imp", 
                              which = "imputation identifiers", 
                              data = list(data, weightit.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(weightit$treat, data = list(data, weightit.data))
    
    #Process covs
    if (is_null(covs <- weightit$covs)) stop("No covariates were specified in the weightit object.", call. = FALSE)
    
    #Get estimand
    estimand <- weightit$estimand
    
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data, weightit.data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data, weightit.data)
    if (is_not_null(distance)) distance <- cbind(distance, prop.score = weightit$ps)
    else if (is_not_null(weightit$ps)) distance <- data.frame(prop.score = weightit$ps)
    else distance <- NULL
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with weightit objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with weightit objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w(weightit)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A$s.weights, weightit$s.weights))) {
        s.weights <- vector.process(s.weights, 
                                    data = list(data, weightit.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, weightit.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    discarded <- weightit$discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "weightit()")
    
    #Process focal
    focal <- weightit$focal
    
    #Get s.d.denom
    if (get.treat.type(treat) != "continuous") {
        s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    }
    else {
        s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
    }
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- weightit$call
    
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
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    treat <- process_treat(t.c[["treat"]], data = data)
    
    #Process covs
    covs <- t.c[["covs"]]
    
    #Get estimand
    estimand <- A$estimand
    
    #Get method
    method <- "matching"
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data)
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with designmatch objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with designmatch objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w.designmatch(dm, treat)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        stop("s.weights are not allowed with designmatch objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = data,
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "the matching function in designmatch")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with designmatch objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
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
    
    
    #Process data and get imp
    if (is_(mimids[["original.datasets"]], "mids")) m.data <- imp.complete(mimids[["original.datasets"]])
    else m.data <- imp.complete(mimids$others$source)
    
    imp <- m.data[[".imp"]]
    
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data, m.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(unlist(lapply(mimids[["models"]][-1], function(m) m[["treat"]])))
    
    #Process covs
    covs <- do.call("rbind", lapply(levels(imp), function(i) {
        m <- mimids[["models"]][-1][[as.numeric(i)]]
        if (is_not_null(m$model$model)) {
            if (nrow(m$model$model) == length(m$treat)) {
                covs <- data.frame(m$model$model[, names(m$model$model) %in% attributes(terms(m$model))$term.labels])
            }
            else {
                #Recreating covs from model object and m$X. Have to do this because when 
                #drop != NULL and reestimate = TRUE, cases are lost. This recovers them.
                
                order <- setNames(attr(m$model$terms, "order"),
                                  attr(m$model$terms, "term.labels"))
                assign <- setNames(attr(m$X, "assign"), colnames(m$X))
                assign1 <- assign[assign %in% which(order == 1)] #Just main effects
                
                dataClasses <- attr(m$model$terms, "dataClasses")
                factors.to.unsplit <- names(dataClasses)[dataClasses %in% c("factor", "character", "logical")]
                f0 <- setNames(lapply(factors.to.unsplit, 
                                      function(x) {
                                          if (dataClasses[x] == "factor")
                                              list(levels = levels(m$model$model[[x]]),
                                                   faclev = paste0(x, levels(m$model$model[[x]])))
                                          else 
                                              list(levels = unique(m$model$model[[x]]),
                                                   faclev = paste0(x, unique(m$model$model[[x]])))
                                      }),
                               factors.to.unsplit)
                covs <- as.data.frame(m$X[, names(assign1)])
                for (j in factors.to.unsplit) {
                    covs <- unsplitfactor(covs, j, sep = "",
                                          dropped.level = f0[[j]]$levels[f0[[j]]$faclev %nin% colnames(m$X)])
                    if (dataClasses[j] == "logical") covs[[j]] <- as.logical(covs[[j]])
                }
            }
        }
        else  {
            t.c <- get.covs.and.treat.from.formula(m$formula, data = m.data[imp == i, , drop = FALSE])
            covs <- t.c[["reported.covs"]]
        }
        return(covs)
    }))
    
    #Get estimand
    estimand <- A$estimand
    
    #Get method
    method <- "matching"
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data, m.data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data, m.data)
    m.distance <- unlist(lapply(mimids[["models"]][-1], function(m) m[["distance"]]))
    if (all(is.na(m.distance))) m.distance <- NULL
    else m.distance <- data.frame(distance = m.distance)
    if (is_not_null(distance)) distance <- cbind(distance, m.distance)
    else distance <- m.distance
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with mimids objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with mimids objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w.mimids(mimids)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        stop("s.weights are not allowed with mimids objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, m.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    discarded <- unlist(lapply(mimids[["models"]][-1], function(m) m[["discarded"]]))
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "matchthem()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with mimids objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
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
     
    #Process data and get imp
    if (is_(wimids[["original.datasets"]], "mids")) w.data <- imp.complete(wimids[["original.datasets"]])
    else w.data <- imp.complete(wimids$others$source)
    
    imp <- w.data[[".imp"]]
    
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
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
                              data = list(data, w.data), 
                              missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(unlist(lapply(wimids[["models"]][-1], function(w) w[["treat"]])))
    
    #Process covs
    covs <- do.call("rbind", lapply(wimids[["models"]][-1], function(w) w[["covs"]]))
    
    #Get estimand
    estimand <- unique(unlist(lapply(wimids[["models"]][-1], function(w) w[["estimand"]])))
    
    #Get method
    method <- "weighting"
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data, w.data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data, w.data)
    w.distance <- unlist(lapply(wimids[["models"]][-1], function(m) m[["ps"]]))
    if (all(is.na(w.distance))) w.distance <- NULL
    else w.distance <- data.frame(distance = w.distance)
    if (is_not_null(distance)) distance <- cbind(distance, w.distance)
    else distance <- w.distance
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with wimids objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with wimids objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w.wimids(wimids)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A$s.weights, unlist(lapply(wimids[["models"]][-1], function(w) w[["s.weights"]]))))) {
        s.weights <- vector.process(s.weights, 
                                    data = list(data, w.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, w.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    discarded <- unlist(lapply(wimids[["models"]][-1], function(w) w[["discarded"]]))
    
    #Process imp and length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "weightthem()")
    
    #Process focal
    focal <- unique(unlist(lapply(wimids[["models"]][-1], function(w) w[["focal"]])))
    
    #Get s.d.denom
    if (get.treat.type(treat) != "continuous") {
        s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    }
    else {
        s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
    }
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
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

#MSMs wth multiple time points
x2base.iptw <- function(iptw, ...) {
    A <- list(...)
    
    #Process iptw
    if (is_not_null(A) && names(A)[1]=="" && is_null(A$stop.method)) A$stop.method <- A[[1]] #for bal.plot
    if (is_null(A$stop.method) && is_not_null(A$full.stop.method)) A$stop.method <- A$full.stop.method
    available.stop.methods <- names(iptw$psList[[1]]$ps)
    if (is_not_null(A$stop.method)) {
        if (any(is.character(A$stop.method))) {
            rule1 <- available.stop.methods[vapply(available.stop.methods, function(x) any(startsWith(tolower(x), tolower(A$stop.method))), logical(1L))]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(available.stop.methods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- available.stop.methods
            }
        }
        else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(available.stop.methods))) {
            if (any(!A$stop.method %in% seq_along(available.stop.methods))) {
                message(paste0("Warning: There are ", length(available.stop.methods), " stop methods available, but you requested ", 
                               word_list(A$stop.method[!A$stop.method %in% seq_along(available.stop.methods)], and.or = "and"),"."))
            }
            rule1 <- available.stop.methods[A$stop.method %in% seq_along(available.stop.methods)]
        }
        else {
            warning("stop.method should be ", word_list(available.stop.methods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- available.stop.methods
        }
    }
    else {
        rule1 <- available.stop.methods
    }
    
    s <- available.stop.methods[match(tolower(rule1), tolower(available.stop.methods))]
    
    #Process data and get imp
    ps.data <- iptw$psList[[1]]$data
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", data = data, missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat.list
    treat.list <- process_treat.list(lapply(iptw$psList, function(x) x$treat), data = list(data, ps.data))
    
    #Process covs.list
    covs.list <- lapply(iptw$psList, function(x) x$data[x$gbm.obj$var.names])
    all.covs <- unique(unlist(lapply(covs.list, names)))
    covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
    
    #Get estimand
    estimand <- substr(toupper(s), nchar(s)-2, nchar(s))
    
    #Get method
    method <- rep("weighting", length(s))
    
    #Process addl.list 
    ntimes <- iptw$nFits
    addl.list <- list.process("addl.list", A[["addl.list"]], ntimes, 
                              "the original call to iptw()",
                              treat.list,
                              covs.list,
                              data,
                              ps.data)
    
    #Process distance
    distance.list <- list.process("distance.list", A[["distance.list"]], ntimes, 
                                  "the original call to iptw()",
                                  treat.list,
                                  covs.list,
                                  data,
                                  ps.data)
    if (is_not_null(distance.list)) {
        for (ti in seq_along(distance.list)) {
            if (length(s) == 1) {
                distance.list[[ti]] <- data.frame(distance[[ti]], prop.score = iptw$psList[[ti]]$ps[[s]])
            }
            else {
                distance.list[[ti]] <- data.frame(distance[[ti]], prop.score = iptw$psList[[ti]]$ps[s])
            }
        }
        
    }
    else {
        distance.list <- vector("list", ntimes)
        for (ti in seq_along(distance.list)) {
            if (length(s) == 1) {
                distance.list[[ti]] <- data.frame(prop.score = iptw$psList[[ti]]$ps[[s]])
            }
            else {
                distance.list[[ti]] <- data.frame(prop.score = iptw$psList[[ti]]$ps[s])
            }
        }
    }
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with iptw objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with iptw objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w(iptw, s)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A$s.weights, iptw$psList[[1]]$sampw))) {
        s.weights <- vector.process(s.weights, 
                                    data = list(data, ps.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, ps.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("weights"),
                       lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                       imp = imp,
                       original.call.to = "iptw()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with iptw objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
    
    #Missing values warning
    if (any(c(any(vapply(covs.list, anyNA, logical(1L))), any(vapply(addl.list, anyNA, logical(1L)))))) {
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
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", data = data, missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat.list
    treat.list <- process_treat.list(A$treat.list, data)
    
    #Process covs.list
    if (is_null(covs.list)) {
        stop("covs.list must be specified.", call. = FALSE)
    }
    if (!is.vector(covs.list, "list")) {
        stop("covs.list must be a list of covariates for which balance is to be assessed at each time point.", call. = FALSE)
    }
    if (any(!vapply(covs.list, is.data.frame, logical(1L)))) {
        stop("Each item in covs.list must be a data frame.", call. = FALSE)
    }
    all.covs <- unique(unlist(lapply(covs.list, names)))
    covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
    
    if (length(treat.list) != length(covs.list)) {
        stop("treat.list must be a list of treatment statuses at each time point.", call. = FALSE)
    }
    
    #Get estimand
    estimand <- NULL
    
    #Get method
    specified <- setNames(rep(FALSE, 1), "weights")
    if (is_not_null(A$weights)) {
        if (!is_(A$weights, c("character", "numeric", "data.frame", "list"))) {
            stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
        }
        specified["weights"] <- TRUE
    }
    
    if (is_null(method <- A$method)) {
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
    ntimes <- length(covs.list)
    addl.list <- list.process("addl.list", A[["addl.list"]], ntimes, 
                              "covs.list",
                              treat.list,
                              covs.list,
                              data)
    
    #Process distance
    distance.list <- list.process("distance.list", A[["distance.list"]], ntimes, 
                                  "covs.list",
                                  treat.list,
                                  covs.list,
                                  data)
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with longitudinal treatments.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with longitudinal treatments.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame.process("weights", 
                                  A[["weights"]], 
                                  do.call("cbind", treat.list), 
                                  do.call("cbind", covs.list), 
                                  data))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        s.weights <- vector.process(s.weights, 
                                  data = data,
                                  name = "s.weights", 
                                  which = "sampling weights",
                                  missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = data,
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("weights"),
                       lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                       imp = imp)
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with longitudinal treatments.", call. = FALSE)
    }
    
    #Get s.d.denom
    if (all(vapply(treat.list, get.treat.type, character(1L)) != "continuous")) {
        s.d.denom <- get.s.d.denom("pooled", estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
    }
    else {
        s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
    }
    
    #Missing values warning
    if (any(c(any(vapply(covs.list, anyNA, logical(1L))), any(vapply(addl.list, anyNA, logical(1L)))))) {
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
    
    treat.list <- covs.list <- vector("list", length(formula.list))
    for (i in seq_along(formula.list)) {
        t.c <- get.covs.and.treat.from.formula(formula.list[[i]], A[["data"]])
        covs.list[[i]] <- t.c[["reported.covs"]]
        treat.list[[i]] <- t.c[["treat"]]
        names(treat.list)[i] <- t.c[["treat.name"]]
    }
    
    X <- do.call("x2base.data.frame.list", c(list(covs.list, treat.list = treat.list), A))
    return(X)
}
x2base.CBMSM <- function(cbmsm, ...) {
    A <- list(...)
    
    #Process CBMSM
    ID <- sort(unique(cbmsm$id))
    times <- sort(unique(cbmsm$time))
    
    #Process data and get imp
    cbmsm.data <- cbmsm$data[cbmsm$time == 1, , drop = FALSE][ID, , drop = FALSE]
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", data = data, missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat.list
    treat.list <- process_treat.list(lapply(times, function(x) cbmsm$treat.hist[ID, x]), 
                                     data = list(data, cbmsm.data))
    
    #Process covs.list
    covs <- cbmsm$data[names(cbmsm$data) %in% attributes(terms(cbmsm$model))$term.labels]
    covs.list <- lapply(times, function(ti) {
        if (ti == 1) {
            out <- setNames(data.frame(covs[cbmsm$time == ti, , drop = FALSE][ID, , drop = FALSE]),
                            paste0(names(covs), "0"))
        }
        else {
            out <- cbind(covs.list[[ti - 1]], 
                         setNames(data.frame(cbmsm$y[cbmsm$time == ti - 1][ID], 
                                             covs[cbmsm$time == ti, , drop = FALSE][ID, , drop = FALSE]),
                                  c(paste0("treat", ti - 1), paste0(names(covs), ti))))
        }
    })
    all.covs <- unique(unlist(lapply(covs.list, names)))
    covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
    
    #Get estimand
    estimand <- NULL
    
    #Get method
    method <- "weighting"
    
    #Process addl.list 
    ntimes <- length(times)
    addl.list <- list.process("addl.list", A[["addl.list"]], ntimes, 
                              "the original call to CBMSM()",
                              treat.list,
                              covs.list,
                              data,
                              cbmsm.data)
    
    #Process distance
    distance.list <- list.process("distance.list", A[["distance.list"]], ntimes, 
                                  "the original call to CBMSM()",
                                  treat.list,
                                  covs.list,
                                  data,
                                  cbmsm.data)
    if (is_not_null(distance.list)) distance.list <- lapply(times, function(x) data.frame(distance.list[[x]], prop.score = cbmsm$fitted.values))
    else if (is_not_null(cbmsm$fitted.values)) distance.list <- lapply(times, function(x) data.frame(prop.score = cbmsm$fitted.values))
    else distance.list <- NULL
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with CBMSM objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with CBMSM objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w(cbmsm)[ID]))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        stop("sampling weights are not allowed with CBMSM objects.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, cbmsm.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("weights"),
                       lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                       imp = imp,
                       original.call.to = "CBMSM()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with CBMSM objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    if (all(vapply(treat.list, get.treat.type, character(1L)) != "continuous")) {
        s.d.denom <- get.s.d.denom("pooled", estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
    }
    else {
        s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
    }
    
    #Missing values warning
    if (any(c(any(vapply(covs.list, anyNA, logical(1L))), any(vapply(addl.list, anyNA, logical(1L)))))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- cbmsm$call
    
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
    weightitMSM.data <- weightitMSM$data
    d.e.in.w <- vapply(c("covs.list", "exact", "by", "moderator"), function(x) is_not_null(weightitMSM[[x]]), logical(1L))
    if (any(d.e.in.w)) weightitMSM.data2 <- do.call("cbind", c(do.call(cbind, weightitMSM$covs.list), weightitMSM[c("exact", "by", "moderator")])[d.e.in.w])
    else weightitMSM.data2 <- NULL
    
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", data = data, missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat.list
    treat.list <- process_treat.list(weightitMSM$treat.list,
                                     data = list(data, weightitMSM.data, weightitMSM.data2))    
    #Process covs.list
    covs.list <- weightitMSM$covs.list
    all.covs <- unique(unlist(lapply(covs.list, names)))
    covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
    
    #Get estimand
    estimand <- weightitMSM$estimand
    
    #Get method
    method <- "weighting"
    
    #Process addl.list 
    ntimes <- length(covs.list)
    addl.list <- list.process("addl.list", A[["addl.list"]], ntimes, 
                              "the original call to weightitMSM()",
                              treat.list,
                              covs.list,
                              data,
                              weightitMSM.data,
                              weightitMSM.data2)
    
    #Process distance
    distance.list <- list.process("distance.list", A[["distance.list"]], ntimes, 
                                  "the original call to weightitMSM()",
                                  treat.list,
                                  covs.list,
                                  data,
                                  weightitMSM.data,
                                  weightitMSM.data2)
    
    if (is_not_null(distance.list)) distance.list <- lapply(seq_along(distance.list), function(x) data.frame(distance.list[[x]], prop.score = weightitMSM$ps.list[[x]]))
    else if (is_not_null(weightitMSM$ps.list)) distance.list <- lapply(seq_along(weightitMSM$ps.list), function(x) data.frame(prop.score = weightitMSM$ps.list[[x]]))
    else distance.list <- NULL
    
    #Process subclass
    if (is_not_null(subclass <- A$subclass)) {
        stop("subclasses are not allowed with weightitMSM objects.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- A$match.strata)) {
        stop("match.strata are not allowed with weightitMSM objects.", call. = FALSE)
    }
    
    #Process weights
    if (is_not_null(weights <- data.frame(weights = get.w(weightitMSM)))) {
        weight.check(weights)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- if_null_then(A$s.weights, weightitMSM$s.weights))) {
        s.weights <- vector.process(s.weights, 
                                    data = list(data, weightitMSM.data, weightitMSM.data2),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
        weight.check(s.weights)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, weightitMSM.data, weightitMSM.data2),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("weights"),
                       lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                       imp = imp,
                       original.call.to = "weightitMSM()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        stop("focal is not allowed with weightitMSM objects.", call. = FALSE)
    }
    
    #Get s.d.denom
    if (all(vapply(treat.list, get.treat.type, character(1L)) != "continuous")) {
        s.d.denom <- get.s.d.denom("pooled", estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
    }
    else {
        s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
    }
    
    #Missing values warning
    if (any(c(any(vapply(covs.list, anyNA, logical(1L))), any(vapply(addl.list, anyNA, logical(1L)))))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- weightitMSM$call
    
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
    X.names <- c("covs",
                 "treat",
                 "weights",
                 "distance",
                 "addl",
                 "s.d.denom",
                 "call",
                 "cluster",
                 "imp",
                 "s.weights",
                 "focal",
                 "discarded",
                 "method",
                 "subclass",
                 "covs.list",
                 "treat.list",
                 "distance.list",
                 "addl.list")
    
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
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
              distance = list(name = c("distance", "ps", "pscore","p.score", "propensity.score"),
                              type = c("data.frame", "matrix", "numeric")),
              distance.list = list(name = c("distance.list", "ps.list", "distance", "ps"),
                                   type = c("list")),
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
    
    P <- setNames(vector("list", length(Q)), names(Q))
    names(obj) <- tolower(names(obj))
    
    for (i in names(Q)) {
        if (is_null(A[[i]])) {
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
        assign(i, A[[i]])
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
    if (is_not_null(covs)) A[["covs"]] <- as.data.frame(A[["covs"]])
    
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
        if (is.vector(A[["weights"]], "numeric")) A[["weights"]] <- data.frame(weights = A[["weights"]])
        else A[["weights"]] <- as.data.frame(A[["weights"]])
    }

    #distance
    if (is_not_null(A[["distance"]])) {
        if (is.numeric(A[["distance"]])) {
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
        imp <- A$imp
        if (is_not_null(data <- A$data)) {
            if (is_(data, "mids")) {
                data <- imp.complete(data)
                if (is_null(imp)) imp <- data[[".imp"]]
            }
            else if (!is_(data, "data.frame"))
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
                                  data = list(data), 
                                  missing.okay = FALSE)
            imp <- factor(imp)
        }
        
        #Process treat
        t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
        treat <- process_treat(t.c[["treat"]], 
                               data = data)
        
        #Process covs
        covs <- t.c[["covs"]]
        if (is_null(covs)) {
            stop("covs data.frame must be specified.", call. = FALSE)
        }
        is_(covs, "data.frame", stop = TRUE)
        
        #Get estimand
        estimand <- A$estimand
        
        #Get method
        specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
        if (is_not_null(A$weights)) {
            specified["weights"] <- TRUE
        }
        if (is_not_null(A$subclass)){
            specified["subclass"] <- TRUE
        }
        if (is_not_null(A$match.strata)) {
            specified["match.strata"] <- TRUE
        }
        
        if (is_null(method <- A$method)) {
            if (specified["match.strata"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Assuming \"matching\" and using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                    A$weights <- A$subclass <- NULL
                }
                method <- "matching"
            }
            else if (specified["subclass"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    A$weights <- A$match.strata <- NULL
                }
                method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    A$match.strata <- A$subclass <- NULL
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
                        A$match.strata <- A$subclass <- NULL
                    }
                    method <- "weighting"
                }
                else if (specified["match.strata"]) {
                    message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
                    A$subclass <- NULL
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
                        A$weights <- A$subclass <- NULL
                    }
                    method <- "matching"
                }
                else if (specified["weights"]) {
                    if (sum(specified) > 1) {
                        message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
                        A$match.strata <- A$subclass <- NULL
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
                        A$weights <- A$match.strata <- NULL
                    }
                    method <- "subclassification"
                }
                else if (specified["match.strata"]) {
                    message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
                    A$weights <- NULL
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
                A$match.strata <- A$subclass <- NULL
            }
        }
        
        #Process addl 
        addl <- data.frame.process("addl", A[["addl"]], treat, covs, data)
        
        #Process distance
        distance <- data.frame.process("distance", A[["distance"]], treat, covs, data)
        
        #Process subclass
        if (is_not_null(subclass <- A$subclass)) {
            subclass <- vector.process(subclass, 
                                       data = data,
                                       name = "subclass", 
                                       which = "subclass membership",
                                       missing.okay = TRUE)
            subclass <- factor(subclass)
        }
        
        #Process match.strata
        if (is_not_null(match.strata <- A$match.strata)) {
            match.strata <- vector.process(match.strata, 
                                           data = data,
                                           name = "match.strata", 
                                           which = "matching strata membership",
                                           missing.okay = TRUE)
            weights <- data.frame(weights = strata2weights(match.strata,
                                                           treat = treat))
        }
        
        #Process weights
        if (is_not_null(weights <- data.frame.process("weights", A[["weights"]], treat, covs, data))) {
            
            weight.check(weights)
            
            if (length(method) == 1) {
                method <- rep.int(method, ncol(weights))
            }
            else if (length(method) != ncol(weights)) {
                stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
            }
        }
        
        #Process s.weights
        if (is_not_null(s.weights <- A$s.weights)) {
            s.weights <- vector.process(s.weights, 
                                        data = data,
                                        name = "s.weights", 
                                        which = "sampling weights",
                                        missing.okay = FALSE)
            weight.check(s.weights)
        }
        
        #Process cluster
        if (is_not_null(cluster <- A$cluster)) {
            cluster <- vector.process(cluster, 
                                      data = data,
                                      name = "cluster", 
                                      which = "cluster membership",
                                      missing.okay = FALSE)
            cluster <- factor(cluster)
        }
        
        #Process subset
        if (is_not_null(subset <- A$subset)) {
            if (!is.logical(subset)) {
                stop("The argument to subset must be a logical vector.", call. = FALSE)
            }
            if (anyNA(subset)) {
                warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
                subset[is.na(subset)] <- FALSE
            }
        }
        
        #Process discarded
        discarded <- A$discarded
        
        #Process length
        length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                           data.frames = c("covs", "weights", "distance", "addl"),
                           imp = imp)
        
        #Process focal
        if (is_not_null(focal <- A$focal) && get.treat.type(treat) != "continuous") {
            if (is.numeric(focal)) {
                if (focal <= nunique(treat)) focal <- levels(treat)[focal]
                else 
                    stop(paste0("focal was specified as ", focal, 
                                ", but there are only ", nunique(treat), " treatment groups."), call. = FALSE)
            }
            else {
                if (focal %nin% levels(treat)) 
                    stop(paste0("The name specified to focal is not the name of any treatment group."), call. = FALSE)
            }
        }
        
        #Get s.d.denom
        if (get.treat.type(treat) != "continuous") {
            s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, 
                                         weights = weights, subclass = subclass, 
                                         treat = treat, focal = focal)
        }
        else {
            s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
        }
        
        #Missing values warning
        if (any(c(anyNA(covs), anyNA(addl)))) {
            warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
        }
        
        #Get call
        call <- A$call
        
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
        initial.list.lengths <- c(length(A$formula.list), length(A$covs.list), length(A$treat.list))
        if (!all_the_same(initial.list.lengths[initial.list.lengths != 0])) stop("The lists in the object were not the same length.", call. = FALSE)
        ntimes.guess <- max(initial.list.lengths)
        if (is_null(A$treat.list)) A$treat.list <- vector("list", length(ntimes.guess)) 
        if (is_null(A$covs.list)) A$covs.list <- vector("list", length(ntimes.guess)) 

        #Process data and get imp
        imp <- A$imp
        if (is_not_null(data <- A$data)) {
            if (is_(data, "mids")) {
                data <- imp.complete(data)
                if (is_null(imp)) imp <- data[[".imp"]]
            }
            else if (!is_(data, "data.frame"))
            {
                # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
                data <- NULL
            }
        }
        
        #Process imp
        if (is_not_null(imp)) {
            imp <- vector.process(imp, "imp", "imputation identifiers", data = data, missing.okay = FALSE)
            imp <- factor(imp)
        }
        
        #Process treat.list
        for (i in seq_len(ntimes.guess)) {
            t.c <- use.tc.fd(A$formula.list[[i]], data, A$treat.list[[i]], A$covs.list[[i]])
            
            A$treat.list[[i]] <- t.c[["treat"]]
            A$covs.list[[i]]  <- t.c[["covs"]]
            if (is_not_null(t.c[["treat.name"]])) names(A$treat.list)[i] <- t.c[["treat.name"]]
        }
        treat.list <- process_treat.list(A$treat.list, data)
        
        #Process covs.list
        if (is_null(covs.list <- A$covs.list)) {
            stop("covs.list must be specified.", call. = FALSE)
        }
        if (!is.vector(covs.list, "list")) {
            stop("covs.list must be a list of covariates for which balance is to be assessed at each time point.", call. = FALSE)
        }
        if (any(!vapply(covs.list, is.data.frame, logical(1L)))) {
            stop("Each item in covs.list must be a data frame.", call. = FALSE)
        }
        all.covs <- unique(unlist(lapply(covs.list, names)))
        covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
        
        if (length(treat.list) != length(covs.list)) {
            stop("treat.list must be a list of treatment statuses at each time point.", call. = FALSE)
        }
        
        #Get estimand
        estimand <- NULL
        
        #Get method
        specified <- setNames(rep(FALSE, 1), "weights")
        if (is_not_null(A$weights)) {
            if (!is_(A$weights, c("character", "numeric", "data.frame", "list"))) {
                stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
            }
            specified["weights"] <- TRUE
        }
        
        if (is_null(method <- A$method)) {
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
        ntimes <- length(covs.list)
        addl.list <- list.process("addl.list", A[["addl.list"]], ntimes, 
                                  "covs.list",
                                  treat.list,
                                  covs.list,
                                  data)
        
        #Process distance
        distance.list <- list.process("distance.list", A[["distance.list"]], ntimes, 
                                      "covs.list",
                                      treat.list,
                                      covs.list,
                                      data)
        
        #Process subclass
        if (is_not_null(subclass <- A$subclass)) {
            stop("subclasses are not allowed with longitudinal treatments.", call. = FALSE)
        }
        
        #Process match.strata
        if (is_not_null(match.strata <- A$match.strata)) {
            stop("match.strata are not allowed with longitudinal treatments.", call. = FALSE)
        }
        
        #Process weights
        weights <- data.frame.process("weights", 
                                      A[["weights"]], 
                                      do.call("cbind", treat.list), 
                                      do.call("cbind", covs.list), 
                                      data)
        weight.check(weights)
        
        #Process s.weights
        if (is_not_null(s.weights <- A$s.weights)) {
            s.weights <- vector.process(s.weights, 
                                        data = data,
                                        name = "s.weights", 
                                        which = "sampling weights",
                                        missing.okay = FALSE)
        }
        
        #Process cluster
        if (is_not_null(cluster <- A$cluster)) {
            cluster <- vector.process(cluster, 
                                      data = data,
                                      name = "cluster", 
                                      which = "cluster membership",
                                      missing.okay = FALSE)
            cluster <- factor(cluster)
        }
        
        #Process subset
        if (is_not_null(subset <- A$subset)) {
            if (!is.logical(subset)) {
                stop("The argument to subset must be a logical vector.", call. = FALSE)
            }
            if (anyNA(subset)) {
                warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
                subset[is.na(subset)] <- FALSE
            }
        }
        
        #Process discarded
        
        #Process length
        length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                           data.frames = c("weights"),
                           lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                           imp = imp)
        
        #Process focal
        if (is_not_null(focal <- A$focal)) {
            stop("focal is not allowed with longitudinal treatments.", call. = FALSE)
        }
        
        #Get s.d.denom
        if (all(vapply(treat.list, get.treat.type, character(1L)) != "continuous")) {
            s.d.denom <- get.s.d.denom("pooled", estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
        }
        else {
            s.d.denom <- get.s.d.denom.cont(A$s.d.denom, weights = weights, subclass = subclass)
        }
        
        #Missing values warning
        if (any(c(any(vapply(covs.list, anyNA, logical(1L))), any(vapply(addl.list, anyNA, logical(1L)))))) {
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
