#Functions to convert object to base.bal.tab input

x2base <- function(...) UseMethod("x2base")

x2base.matchit <- function(m, ...) {
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Initializing variables
    
    if (any(class(m) == "matchit.subclass")) {
        subclass <- factor(m$subclass)
        method <- "subclassification"
    }
    else if (any(class(m) == "matchit.full")) {
        subclass <- NULL
        method <- "weighting"
    }
    else {
        subclass <- NULL
        method <- "matching"
    }
    
    weights <- data.frame(weights = m$weights)
    treat <- m$treat
    data <- A$data
    subset <- A$subset
    cluster <- A$cluster
    s.d.denom <- A$s.d.denom
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    
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
    m.data <- m$model$data
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(m.data) == cluster)) {
                cluster <- m.data[[cluster]]
            }
        }
        else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is_(subset, "logical")) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data, m.data))
    }
    
    if (any(is.finite(m$distance))) {
        if (is_not_null(distance)) distance <- cbind(distance, distance = m$distance)
        else distance <- data.frame(distance = m$distance)
    }
    
    #Get s.d.denom
    estimand <- "ATT"
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method = method)
    
    ensure.equal.lengths <- TRUE
    vectors <- c("cluster", "treat", "subset", "subclass")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors, ifnotfound = list(NULL))), 
                          vapply(data.frames, 
                                 function(x) NROW(get0(x)),
                                 numeric(1L))), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to matchit()."), call. = FALSE)
    }
    
    if (any(c(is.na(covs), is.na(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$method <- method
    X$treat <- treat
    X$weights <- weights
    X$discarded <- m$discarded
    X$covs <- covs
    X$distance <- distance
    X$addl <- addl
    X$cluster <- factor(cluster)
    X$call <- m$call
    X$subclass <- factor(subclass)
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.ps <- function(ps, ...) {
    #stop.method
    #s.d.denom
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Initializing variables
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
    estimand <- ps$estimand
    
    weights <- data.frame(get.w(ps, s, estimand))
    treat <- ps$treat
    covs <- ps$data[, ps$gbm.obj$var.names, drop = FALSE]
    data <- A$data
    ps.data <- ps$data
    cluster <- A$cluster
    subset <- A$subset
    s.weights <- ps$sampw
    s.d.denom <- A$s.d.denom
    method <- rep("weighting", ncol(weights))
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(ps.data) == cluster)) {
                cluster <- ps.data[[cluster]]
            }
            else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
        }
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data, ps.data))
    }
    
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
    
    ensure.equal.lengths <- TRUE
    vectors <- c("s.weights", "cluster", "subset")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to ps()."), call. = FALSE)
    }
    
    #Get s.d.denom
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)

    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$addl <- addl
    X$covs <- covs
    X$call <- ps$parameters
    X$cluster <- factor(cluster)
    X$method <- method
    X$s.weights <- s.weights
    
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.mnps <- function(mnps, ...) {
    #stop.method
    #s.d.denom
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Initializing variables
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
    
    weights <- data.frame(get.w(mnps, s))
    treat <- mnps$treatVar
    covs <- mnps$data[mnps$psList[[1]]$gbm.obj$var.names]
    data <- A$data
    cluster <- A$cluster
    subset <- A$subset
    s.weights <- mnps$sampw
    s.d.denom <- A$s.d.denom
    focal <- mnps$treatATT
    method <- rep("weighting", ncol(weights))
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    mnps.data <- mnps$data
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(mnps.data) == cluster)) {
                cluster <- mnps.data[[cluster]]
            }
        }
        else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data, mnps.data))
    }    
    # model.ps <- setNames(as.data.frame(lapply(mnps$psList, function(x) x$ps)), names(mnps$psList))
    # model.ps.combined <- numeric(length(treat))
    # for (i in levels(treat)) {
    #     model.ps.combined[treat == i] <- model.ps[treat == i, i]
    # }
    # if (length(distance) > 0) distance <- cbind(distance, prop.score = model.ps.combined)
    # else distance <- data.frame(prop.score = model.ps.combined)
    
    ensure.equal.lengths <- TRUE
    vectors <- c("s.weights", "cluster", "subset")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to ps()."), call. = FALSE)
    }
    
    #Get s.d.denom
    estimand <- mnps$estimand
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal, method)
    
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$addl <- addl
    X$covs <- covs
    X$call <- NULL
    X$cluster <- factor(cluster)
    X$s.weights <- mnps$sampw
    X$focal <- focal
    X$method <- method
    
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- get.X.class(X)
    
    return(X)
}
x2base.ps.cont <- function(ps.cont, ...) {
    #stop.method
    #s.d.denom
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Initializing variables
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
    
    weights <- data.frame(get.w(ps.cont, s))
    treat <- ps.cont$treat
    covs <- ps.cont$data[, ps.cont$gbm.obj$var.names, drop = FALSE]
    data <- A$data
    ps.data <- ps.cont$data
    cluster <- A$cluster
    subset <- A$subset
    s.weights <- ps.cont$sampw
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(ps.data) == cluster)) {
                cluster <- ps.data[[cluster]]
            }
            else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
        }
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data, ps.data))
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("s.weights", "cluster", "subset")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to ps.cont()."), call. = FALSE)
    }
    
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$addl <- addl
    X$covs <- covs
    X$call <- ps.cont$parameters
    X$cluster <- factor(cluster)
    X$method <- rep("weighting", ncol(weights))
    X$s.weights <- s.weights
    
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "cont"
    
    return(X)
}
x2base.Match <- function(Match, ...) {
    
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Checks
    if (is_not_null(Match) && !is.list(Match)) {
        stop("'Match' object contains no valid matches")}
    
    #Get treat and covs
    data <- A$data
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    
    #Initializing variables
    m <- Match
    estimand <- m$estimand
    method <- "matching"
    s.d.denom <- A$s.d.denom
    
    treat <- t.c[["treat"]]
    covs  <- t.c[["covs"]]
    
    weights <- data.frame(weights = get.w.Match(m))
    
    cluster <- A$cluster
    subset <- A$subset
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
            cluster <- data[[cluster]]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
        
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    dropped <- rep(FALSE, length(treat))
    if (is_not_null(m$index.dropped)) dropped[m$index.dropped] <- TRUE
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data))
    }
    
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("treat", "cluster", "subset")
    data.frames <- c(covs.data, "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in names(lengths)[names(lengths) != "weights"]) {
            if (lengths[i] > 0 && lengths[i] != lengths["weights"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original call to Match()."), call. = FALSE)
    }
    
    #Get s.d.denom
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
    
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$treat <- treat
    X$weights <- weights
    X$discarded <- dropped
    X$distance <- NULL #NAs in distance bcause of incomplete list in Match object
    X$addl <- addl
    X$covs <- covs
    X$call <- NULL
    X$method <- method
    X$cluster <- factor(cluster)
    
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
    
    X <- do.call("x2base.data.frame", c(list(covs = covs, treat = treat), A))
    return(X)
}
x2base.data.frame <- function(covs, ...) {
    #treat
    #data
    #weights
    #distance
    #subclass
    #match.strata
    #addl
    #s.d.denom
    #method
    #cluster
    #estimand
    #imp
    
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    treat <- A$treat
    data <- A$data
    weights <- A$weights
    distance <- A$distance
    subclass <- A$subclass
    match.strata <- A$match.strata
    cluster <- A$cluster
    addl <- A$addl
    s.d.denom <- A$s.d.denom
    method <- A$method
    estimand <- A$estimand
    imp <- A$imp
    s.weights <- A$s.weights
    subset <- A$subset
    focal <- A$focal
    
    #Checks
    if (is_null(covs)) {
        stop("covs data.frame must be specified.", call. = FALSE)
    }
    is_(covs, "data.frame", stop = TRUE)
    
    #Process data
    if (is_not_null(data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if ("imp" %nin% names(A)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
        {
            warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }

    specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
    if (is_not_null(weights)) {
        if (!is_(weights, c("character", "numeric", "data.frame", "list"))) {
            stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
        }
        specified["weights"] <- TRUE
    }
    if (is_not_null(subclass)){
        if (!is_(subclass, c("character", "numeric", "factor"))) {
            stop("The argument to subclass must be a vector of subclass membership or the (quoted) name of a variable in data that contains subclass membership.", call. = FALSE)
        }
        specified["subclass"] <- TRUE
    }
    if (is_not_null(match.strata)) {
        if (!is_(match.strata, c("character", "numeric", "factor"))) {
            stop("The argument to match.strata must be a vector of match stratum membership or the (quoted) name of a variable in data that contains match stratum membership.", call. = FALSE)
        }
        specified["match.strata"] <- TRUE
    }
    
    #Getting method
    if (is_null(method)) {
        if (specified["match.strata"]) {
            if (sum(specified) > 1) {
                message(word_list(names(specified)[specified]), " are specified. Assuming \"matching\" and using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                weights <- subclass <- NULL
            }
            method <- "matching"
        }
        else if (specified["subclass"]) {
            if (sum(specified) > 1) {
                message(word_list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                weights <- match.strata <- NULL
            }
            method <- "subclassification"
            #weights <- rep(1, nrow(covs))
        }
        else if (specified["weights"]) {
            if (sum(specified) > 1) {
                message(word_list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                match.strata <- subclass <- NULL
            }
            else {
                message("Assuming \"weighting\". If not, specify with an argument to method.")
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
                    match.strata <- subclass <- NULL
                }
                method <- "weighting"
            }
            else if (specified["match.strata"]) {
                message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
                subclass <- NULL
                method <- "matching"
            }
            else if (specified["subclass"]) {
                message("method = \"weighting\" is specified, but no weights are present. Assuming \"subclassification\" and using subclass instead.")
                method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else {
                method <- "matching"
            }
        }
        else if (specified.method == "matching") {
            if (specified["match.strata"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                    weights <- subclass <- NULL
                }
                method <- "matching"
            }
            else if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
                    match.strata <- subclass <- NULL
                }
                method <- "matching"
            }
            else if (specified["subclass"]) {
                message("method = \"matching\" is specified, but no weights or match.strata are present. Assuming \"subclassification\" and using subclass instead.")
                method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else {
                method <- "matching"
            }
        }
        else if (specified.method == "subclassification") {
            if (specified["subclass"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    weights <- match.strata <- NULL
                }
                method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else if (specified["match.strata"]) {
                message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
                weights <- NULL
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
            stop("Only weights can be specified with mutiple methods.", call. = FALSE)
        }
        else if (!specified["weights"]) {
            warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
            method <- "matching"
        }
        else {
            #Matching and/or weighting with various weights
            method <- specified.method
            match.strata <- subclass <- NULL
        }
    }
    
    if (is_not_null(cluster) && !is_(cluster, c("character", "numeric", "factor"))) {
        stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (is_not_null(imp) && !is_(imp, c("character", "numeric", "factor"))) {
        stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
    }
    
    #Process treat
    if (is_null(treat)) stop("treat must be specified.", call. = FALSE)
    
    
    if (is.character(treat) && length(treat)==1 && treat %in% names(data)) {
        treat <- data[[treat]]
    }
    else if (is_(treat, c("numeric", "logical", "factor")) || (is.character(treat) && length(treat) > 1)) {
        treat <- treat
    }
    else stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
    
    if (sum(is.na(treat)) > 0)
        stop("Missing values exist in treat.", call. = FALSE)
    
    if (is_binary(treat)) {
        treat <- binarize(treat)
    }
    else if (is.character(treat)) {
        treat <- factor(treat)
    }
    
    #Process weights, addl, and distance
    for (i in c("weights", "addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data))
    }
    
    #Process subclass
    if (is_not_null(subclass)) {
        if (is_(subclass, c("numeric", "factor")) || (is.character(subclass) && length(subclass) > 1)) {
            subclass <- factor(subclass)
        }
        else if (is.character(subclass) && length(subclass)==1 && any(names(data) == subclass)) {
            subclass <- factor(data[[subclass]])
        }
        else stop("The name supplied to subclass is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata)) {
        if (is.character(match.strata) && length(match.strata)==1) {
            if (any(names(data) == match.strata)) {
                match.strata <- data[[match.strata]]
            }
            else stop("The name supplied to match.strata is not the name of a variable in data.", call. = FALSE)
        }
    }
    
    #Process sampling weights
    if (is_not_null(s.weights)) {
        if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
            stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
        }
        if (is.character(s.weights) && length(s.weights)==1) {
            if (any(names(data) == s.weights)) {
                s.weights <- data[[s.weights]]
            }
            else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
        }
        if (anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is_(cluster, c("numeric", "factor")) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
            cluster <- data[[cluster]]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    #Process imp
    if (is_not_null(imp)) {
        if (is_(imp, c("numeric", "factor")) || (is.character(imp) && length(imp) > 1)) {
            imp <- imp
        }
        else if (is.character(imp) && length(imp)==1 && any(names(data) == imp)) {
            imp <- data[[imp]]
        }
        else stop("The name supplied to imp is not the name of a variable in data.", call. = FALSE)
        
        imp.lengths <- vapply(unique(imp), function(i) sum(imp == i), numeric(1L))
        
        if (all_the_same(imp.lengths)) { #all the same
            for (i in vectors) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) { 
                    if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        
                        temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, temp.merge[[4]][order(temp.merge[[3]])])
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
            for (i in data.frames) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) {
                    if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, setNames(temp.merge[order(temp.merge[[3]]), -c(1:3), drop = FALSE], names(get(i))))
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
        }
        else {
            problematic <- lengths > 0 & lengths != length(imp)
        }
        if (any(problematic)) {
            stop(paste0(word_list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
        }
        else ensure.equal.lengths <- FALSE
    }
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
    }
    
    #Turn match.strata into weights
    if (is_not_null(match.strata)) {
        weights <- data.frame(weights = match.strata2weights(match.strata = match.strata,
                                                             treat = treat,
                                                             covs = covs
        ))
    }
   
    if (is_not_null(weights)) {
        if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
        if (any(vapply(weights, function(x) !is.numeric(x), logical(1L)))) stop("All weights must be numeric.", call. = FALSE)
        if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
        
        if (length(method) == 1) {
            method <- rep(method, ncol(weights))
        }
        else if (length(method) != ncol(weights)) {
            stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
        }
        
    }
    
    #Check focal
    if (is_not_null(focal) && is.factor(treat)) {
        if (is.numeric(focal)) {
            if (focal <= nunique(treat)) focal <- levels(treat)[focal]
            else 
                stop(paste0("focal was specified as ", focal, 
                            ", but there are only ", levels(treat), " treatment groups."), call. = FALSE)
        }
        else {
            if (!any(levels(treat) == focal)) 
                stop(paste0("The name specified to focal is not the name of any treatment group."), call. = FALSE)
        }
    }
    
    #Get s.d.denom
    if (is_binary(treat) || !is.numeric(treat)) { #non-continuous
        X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal, method = method)
    }
    
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$method <- method
    X$covs <- covs
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$subclass <- subclass
    X$cluster <- factor(cluster)
    X$call <- NULL
    X$addl <- addl
    X$imp <- factor(imp)
    X$s.weights <- s.weights
    X$focal <- focal
    
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- get.X.class(X)
    
    return(X)
}
x2base.CBPS <- function(cbps.fit, ...) {
    #s.d.denom
    #cluster
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    #Checks
    
    treat <- model.response(model.frame(cbps.fit$terms, cbps.fit$data))
    covs <- cbps.fit$data[names(cbps.fit$data) %in% attributes(terms(cbps.fit))$term.labels]
    data <- A$data
    s.weights <- A$s.weights
    subset <- A$subset
    weights <- data.frame(weights = get.w(cbps.fit, use.weights = A$use.weights))
    cluster <- A$cluster
    s.d.denom <- A$s.d.denom
    estimand <- A$estimand
    method <- "weighting"
    c.data <- cbps.fit$data
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    #Process sampling weights
    if (is_not_null(s.weights)) {
        if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
            stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
        }
        if (is.character(s.weights) && length(s.weights)==1) {
            if (any(names(data) == s.weights)) {
                s.weights <- data[[s.weights]]
            }
            else if (any(names(c.data) == s.weights)) {
                s.weights <- data[[s.weights]]
            }
            else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
        }
    }
    else s.weights <- rep(1, length(treat))
    
    weights <- weights/s.weights #Because CBPS weights contain s.weights in them
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(c.data) == cluster)) {
                cluster <- c.data[[cluster]]
            }
        }
        else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data, c.data))
    }
    
    if (!any(class(cbps.fit) == "CBPSContinuous") && !is_binary(treat)) {
        #Multinomial
        # model.ps <- setNames(data.frame(cbps.fit$fitted.values), levels(treat))
        # model.ps.combined <- numeric(length(treat))
        # for (i in levels(treat)) {
        #     model.ps.combined[treat == i] <- model.ps[treat == i, i]
        # }
        # if (is_not_null(distance)) distance <- cbind(distance, prop.score = model.ps.combined)
        # else distance <- data.frame(prop.score = model.ps.combined)
    }
    else {
        if (all_the_same(cbps.fit$fitted.values)) {
            if (is_null(distance)) distance <- NULL
        }
        else {
            if (is_null(distance)) distance <- data.frame(prop.score = cbps.fit$fitted.values)
            else distance <- cbind(distance, prop.score = cbps.fit$fitted.values)
        }
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("cluster", "s.weights", "subset")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in names(lengths)[names(lengths) != "covs"]) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to CBPS()."), call. = FALSE)
    }
    
    #Get s.d.denom
    if (!any(class(cbps.fit) == "CBPSContinuous")) {
        if (is_binary(treat)) {
            X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
        }
        else {
            X$s.d.denom <- "pooled"
        }
    }
    
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$distance <- distance
    X$addl <- addl
    X$weights <- weights
    X$treat <- treat
    X$covs <- covs
    X$cluster <- factor(cluster)
    X$call <- cbps.fit$call
    X$s.weights <- s.weights
    X$method <- method
    
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- get.X.class(X)
    
    return(X)
}
x2base.ebalance <- function(ebalance, ...) {
    #formula
    #data
    #treat
    #covs
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Get treat and covs
    data <- A$data
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    
    #Initializing variables
    
    treat <- t.c[["treat"]]
    covs  <- t.c[["covs"]]
    cluster <- A$cluster
    subset <- A$subset
    weights <- data.frame(weights = get.w(ebalance, treat))
    method <- "weighting"
    s.d.denom <- A$s.d.denom
    estimand <- "ATT"
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
            cluster <- data[[cluster]]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
        
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data))
    }
    
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("treat", "cluster", "subset")
    data.frames <- c(covs.data, "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in names(lengths)[names(lengths) != "weights"]) {
            if (lengths[i] > 0 && lengths[i] != lengths["weights"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original call to ebalance()."), call. = FALSE)
    }
    
    #Get s.d.denom
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
    
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$treat <- treat
    X$weights <- weights
    X$covs <- covs
    X$distance <- distance
    X$addl <- addl
    X$call <- NULL
    X$method <- method
    X$cluster <- factor(cluster)
    
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
}
x2base.ebalance.trim <- x2base.ebalance
x2base.optmatch <- function(optmatch, ...) {
    #formula
    #data
    #treat
    #covs
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Get treat and covs
    data <- A$data
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    
    #Initializing variables
    treat <- binarize(t.c$treat)
    covs  <- t.c$covs
    distance <- A$distance
    subset <- A$subset
    cluster <- A$cluster
    s.d.denom <- A$s.d.denom
    estimand <- "ATT"
    method <- "matching"
    
    #Process match.strata (optmatch)
    if (length(optmatch) != length(treat) || length(optmatch) != nrow(covs)) {
        stop(paste0("The optmatch object must have the same length as ", ifelse(attr(t.c, "which")=="fd", "data", "covs"), "."), call. = FALSE)
    }
    a <- attributes(optmatch)
    
    d.reordered <- setNames(seq_len(nrow(covs)), rownames(covs))[a$names]
    
    weights <- data.frame(weights = get.w(optmatch))
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
            cluster <- data[[cluster]]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data))
    }    
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("treat", "cluster", "subset")
    data.frames <- c(covs.data, "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in names(lengths)[names(lengths) != "weights"]) {
            if (lengths[i] > 0 && lengths[i] != lengths["weights"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original call to optmatch()."), call. = FALSE)
    }
    
    #Get s.d.denom
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
    
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$treat <- treat[d.reordered]
    X$distance <- distance[d.reordered, , drop = FALSE]
    X$covs <- covs[d.reordered, , drop = FALSE]
    X$weights <- weights
    X$addl <- addl[d.reordered, , drop = FALSE]
    X$call <- attr(optmatch, "call")
    X$method <- "matching"
    X$cluster <- factor(cluster[d.reordered])
    
    
    subset <- subset[d.reordered]
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
    
}
x2base.weightit <- function(weightit, ...) {
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Initializing variables
    estimand <- weightit$estimand
    weights <- data.frame(weights = get.w(weightit))
    treat <- weightit$treat
    covs <- weightit$covs
    s.weights <- weightit$s.weights
    data <- A$data
    cluster <- A$cluster
    imp <- A$imp
    subset <- A$subset
    s.d.denom <- A$s.d.denom
    focal <- weightit$focal
    method <- rep("weighting", ncol(weights))
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    d.e.in.w <- vapply(c("covs", "exact", "by"), function(x) is_not_null(weightit[[x]]), logical(1L))
    if (any(d.e.in.w)) weightit.data <- do.call("cbind", unname(weightit[c("covs", "exact", "by")[d.e.in.w]]))
    else weightit.data <- NULL
    
    if (has.treat.type(treat)) {
        treat.type <- get.treat.type(treat)
    }
    else if (is_not_null(weightit$treat.type)) {
        treat.type <- weightit$treat.type
    }
    else {
        treat.type <- get.treat.type(assign.treat.type(treat))
    }
    
    if (is_null(covs)) stop("No covariates were specified in the weightit object.", call. = FALSE)
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (!is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
            stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
        }
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(weightit.data) == cluster)) {
                cluster <- weightit.data[[cluster]]
            }
            else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
        }
    }
    
    #Process imp
    if (is_not_null(imp) && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
        stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data, weightit.data))
    }    
    
    if (is_not_null(distance)) distance <- cbind(distance, prop.score = weightit$ps)
    else if (is_not_null(weightit$ps)) distance <- data.frame(prop.score = weightit$ps)
    else distance <- NULL
    
    ensure.equal.lengths <- TRUE
    vectors <- c("s.weights", "cluster", "subset")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    #Process imp further
    if (is_not_null(imp)) {
        if (is.numeric(imp) || is.factor(imp) || (is.character(imp) && length(imp)>1)) {
            imp <- imp
        }
        else if (is.character(imp) && length(imp)==1 && any(names(data) == imp)) {
            imp <- data[[imp]]
        }
        else if (is.character(imp) && length(imp)==1 && any(names(weightit.data) == imp)) {
            imp <- weightit.data[[imp]]
        }
        else stop("The name supplied to imp is not the name of a variable in data.", call. = FALSE)
        
        imp.lengths <- vapply(unique(imp), function(i) sum(imp == i), numeric(1L))
        
        if (all_the_same(imp.lengths)) { #all the same
            for (i in vectors) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) { 
                    if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        
                        temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, temp.merge[[4]][order(temp.merge[[3]])])
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
            for (i in data.frames) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) {
                    if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        temp.var <- data.frame(sort(imp),rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, setNames(temp.merge[order(temp.merge[[3]]), -c(1:3), drop = FALSE], names(get(i))))
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
        }
        else {
            problematic <- lengths > 0 & lengths != length(imp)
        }
        if (any(problematic)) {
            stop(paste0(word_list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
        }
        else ensure.equal.lengths <- FALSE
    }
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
    }
    
    #Get s.d.denom
    if (treat.type != "continuous") {
        X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal, method)
    }
    
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$addl <- addl
    X$covs <- covs
    X$cluster <- factor(cluster)
    X$method <- method
    X$imp <- factor(imp)
    X$s.weights <- weightit$s.weights
    X$discarded <- weightit$discarded
    X$focal <- focal
    X$call <- weightit$call
    
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- get.X.class(X)
    
    return(X)
}
x2base.designmatch <- function(dm, ...) {
    #formula
    #data
    #treat
    #covs
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Get treat and covs
    data <- A$data
    
    if (all(c("id_1", "id_2") %in% names(dm))) {
        t.c <- use.tc.fd(A$formula, data, A$treat, A$covs, needs.treat = FALSE)
        covs  <- t.c$covs
        
        if (anyDuplicated(c(dm[["id_1"]], dm[["id_2"]])) != 0) {
            stop("Some units are used more than once. Balance cannot be checked.", call. = FALSE)
        }
        treat <- rep(NA_real_, nrow(covs))
        treat[dm[["id_1"]]] <- 1
        treat[dm[["id_2"]]] <- 0
        
        in.matched <- !is.na(treat)
        weights <- NULL
    }
    else {
        t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
        covs  <- t.c$covs
        treat <- binarize(t.c$treat)
        in.matched <- rep(TRUE, length(treat))
        
        weights <- data.frame(weights = get.w.designmatch(dm, treat))
        
    }
    
    #Initializing variables
    distance <- A$distance
    subset <- A$subset
    cluster <- A$cluster
    estimand <- A$estimand
    s.d.denom <- A$s.d.denom
    method <- "matching"
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster) > 1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
            cluster <- data[[cluster]]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data))
    }    
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("treat", "cluster", "subset")
    data.frames <- c(covs.data, "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is_null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in names(lengths)[names(lengths) != "treat"]) {
            if (lengths[i] > 0 && lengths[i] != lengths["treat"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original call to designmatch()."), call. = FALSE)
    }
    
    #Get s.d.denom
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
    
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    
    X$treat <- treat
    X$distance <- distance
    X$covs <- covs
    X$weights <- weights
    X$addl <- addl
    X$call <- NULL
    X$method <- method
    X$cluster <- factor(cluster)
    
    X <- subset_X(X, in.matched)
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "binary"
    
    return(X)
    
}
.x2base.mimids <- function(mimids, ...) {
    
    nimp <- length(mimids[[2]]) - 1
    if (nimp == 1) {
        X <- x2base.matchit(mimids[[2]][-1][[1]])
    }
    else {
        Xs <- vector("list", nimp)
        for (i in 1:nimp) {
            Xs[[i]] <- x2base.matchit(mimids[[2]][-1][[i]])
            Xs[[i]][["imp"]] <- factor(rep(i, length(Xs[[i]][["treat"]])),
                                       levels = as.character(seq_along(Xs)))
        }
        n <- length(Xs[[1]][["treat"]])
        
        X <- setNames(vector("list", length(Xs[[1]])),
                      names(Xs[[1]]))
        for (x in names(X)) {
            if (is.data.frame(Xs[[1]][[x]]) || is.matrix(Xs[[1]][[x]])) {
                X[[x]] <- do.call("rbind", lapply(Xs, function(x_) x_[[x]]))
            }
            else if (length(Xs[[1]][[x]]) == n) {
                X[[x]] <- unlist(lapply(Xs, function(x_) x_[[x]]))
            }
            else {
                X[[x]] <- Xs[[1]][[x]]
            }
        }
        
        class(X) <- "imp"
    }
    
    if ("wimids" %in% class(mimids)) {
        X$weights[] <- unlist(lapply(mimids[[4]][-1], function(i) i[["inverse.weights"]]))
        X$s.d.denom <- "pooled"
        X$method <- "weighting"
    }
    
    return(X)
}
x2base.mimids <- function(mimids, ...) {
    
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
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Initializing variables
    nimp <- length(mimids[[4]]) - 1
    
    if (any(class(mimids) == "wimids")) {
        subclass <- NULL
        method <- "weighting"
        estimand <- "ATE"
        weights <- do.call("rbind", lapply(mimids[[4]][-1], function(x) x["inverse.weights"]))
        m.distance <- data.frame(distance = unlist(lapply(mimids[[2]][-1], function(m) m[["distance"]])))
        if ("average.distance" %in% names(mimids[[4]][-1][[1]])) {
            m.distance$average.distance <- rep(rowMeans(do.call("cbind", lapply(mimids[[2]][-1], function(m) m[["distance"]]))), nimp)
        }
        discarded <- NULL
    } else {
        subclass <- NULL
        method <- "matching"
        estimand <- "ATT"
        weights <- do.call("rbind", lapply(mimids[[2]][-1], function(x) data.frame(weights = x[["weights"]])))
        m.distance <- data.frame(distance = unlist(lapply(mimids[[2]][-1], function(m) m[["distance"]])))
        if ("average.distance" %in% names(mimids[[4]][-1][[1]])) names(m.distance) <- "ave.distance"
        discarded <- if ("discarded" %in% names(mimids[[2]][-1][[1]])) unlist(lapply(mimids[[2]][-1], function(x) x[["discarded"]])) else NULL
    }
    
    #Process data
    data <- A$data
    
    if (is_not_null(data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if ("imp" %nin% names(A)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
        {
            warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    m.data1 <- do.call("rbind", lapply(mimids[[2]][-1], function(m) m$model$data))
    m.data2 <- do.call("rbind", mimids[[4]][-1])
    
    subset <- A$subset
    cluster <- A$cluster
    s.d.denom <- A$s.d.denom
    imp <- m.data2[[".imp"]]
    weights[is.na(weights),] <- 0
    
    treat <- unlist(lapply(mimids[[2]][-1], function(m) m[["treat"]]))
    covs <- do.call("rbind", lapply(1:nimp, function(i) {
        m <- mimids[[2]][-1][[i]]
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
                for (i in factors.to.unsplit) {
                    covs <- unsplitfactor(covs, i, sep = "",
                                          dropped.level = f0[[i]]$levels[f0[[i]]$faclev %nin% colnames(m$X)])
                    if (dataClasses[i] == "logical") covs[[i]] <- as.logical(covs[[i]])
                }
            }
            
        }
        else if (!anyNA(mimids[[4]][-1][[i]])) {
            t.c <- get.covs.and.treat.from.formula(m$formula, data = mimids[[4]][-1][[i]])
            covs <- t.c[["reported.covs"]]
        }
        else if (is_not_null(data)) {
            if (nrow(data) == length(mimids[[2]][-1][[i]][["treat"]])) {
                t.c <- get.covs.and.treat.from.formula(m$formula, data = data)
                covs <- t.c[["reported.covs"]]
            }
            else if (nrow(data) == length(treat)) {
                t.c <- get.covs.and.treat.from.formula(m$formula, 
                                                       data = data[((i-1)*length(mimids[[2]][-1][[i]][["treat"]]) + 1):(i*length(mimids[[2]][-1][[i]][["treat"]])), , drop = FALSE])
                covs <- t.c[["reported.covs"]]
            }
            else {
                covs <- data.frame(m$X)
            }
        }
        else {
            covs <- data.frame(m$X)
        }
        return(covs)
    }))
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(m.data1) == cluster)) {
                cluster <- m.data1[[cluster]]
            }
        }
        else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is_(subset, "logical")) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data, m.data1, m.data2))
    }

    if (is_not_null(distance)) distance <- cbind(distance, m.distance)
    else distance <- m.distance
    
    #Get s.d.denom
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method = method)
    
    ensure.equal.lengths <- TRUE
    vectors <- c("treat", "cluster", "subset")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is_null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L))), c(vectors, data.frames))
    #Process imp further
    if (is_not_null(imp)) {
        
        imp.lengths <- vapply(unique(imp, nmax = nimp), function(i) sum(imp == i), numeric(1L))
        
        if (all_the_same(imp.lengths)) { #all the same
            for (i in vectors) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) { 
                    if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        
                        temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, temp.merge[[4]][order(temp.merge[[3]])])
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
            for (i in data.frames) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) {
                    if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        temp.var <- data.frame(sort(imp),rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, setNames(temp.merge[order(temp.merge[[3]]), -c(1:3), drop = FALSE], names(get(i))))
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
        }
        else {
            problematic <- lengths > 0 & lengths != length(imp)
        }
        if (any(problematic)) {
            stop(paste0(word_list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
        }
        else ensure.equal.lengths <- FALSE
    }
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
    }
    
    if (any(c(is.na(covs), is.na(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$method <- method
    X$treat <- treat
    X$weights <- weights
    X$discarded <- discarded
    X$covs <- covs
    X$distance <- distance
    X$addl <- addl
    X$cluster <- factor(cluster)
    X$imp <- factor(imp)
    X$call <- NULL
    X$subclass <- factor(subclass)
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- "imp"
    
    return(X)
}
x2base.wimids <- x2base.mimids

#MSMs wth multiple time points
x2base.iptw <- function(iptw, ...) {
    A <- list(...)
    
    X.names <- c("covs.list",
                 "treat.list",
                 "weights",
                 "distance.list",
                 "addl.list",
                 "s.d.denom",
                 "call",
                 "cluster",
                 "imp",
                 "s.weights",
                 "focal",
                 "discarded",
                 "method",
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
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
    estimand <- substr(toupper(s), nchar(s)-2, nchar(s))
    
    weights <- data.frame(get.w(iptw, s))
    treat.list <- lapply(iptw$psList, function(x) x$treat)
    covs.list <- lapply(iptw$psList, function(x) x$data[x$gbm.obj$var.names])
    subset <- A$subset
    data <- A$data
    cluster <- A$cluster
    ps.data <- iptw$psList[[1]]$data
    s.weights <- iptw$psList[[1]]$sampw
    ntimes <- iptw$nFits
    s.d.denom <- A$s.d.denom
    method <- rep("weighting", ncol(weights))
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    #Order covs.list
    all.covs <- unique(unlist(lapply(covs.list, names)))
    covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(ps.data) == cluster)) {
                cluster <- ps.data[[cluster]]
            }
            else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
        }
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl.list", "distance.list")) {
        assign(i, list.process(i, A[[i]], ntimes, 
                               "the original call to iptw()",
                               treat.list,
                               covs.list,
                               data,
                               ps.data))
    }
    
    
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
    
    ensure.equal.lengths <- TRUE
    vectors <- c("s.weights", "cluster", "subset")
    data.frames <- c("weights")
    lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L)),
                          vapply(lists, function(x) {
                              if (is.null(get0(x))) 0 
                              else if (is.vector(get(x))) {
                                  if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                                  else max(lengths(get(x)))
                              }
                              else max(vapply(get(x), function(y) if (is_not_null(y)) nrow(y) else 0, numeric(1L)))
                          }, numeric(1L))), c(vectors, data.frames, lists))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames, lists)) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to iptw()."), call. = FALSE)
    }
    
    #Get s.d.denom
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat.list[[1]], focal = NULL, method)
    
    if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$weights <- weights
    X$treat.list <- lapply(treat.list, function(x) x)
    X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
    X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
    X$covs.list <- lapply(covs.list, function(x) x)
    X$call <- NULL
    X$cluster <- factor(cluster)
    X$method <- rep("weighting", ncol(weights))
    X$s.weights <- s.weights
    
    if (is_null(subset)) subset <- rep(TRUE, lengths["treat.list"])
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- get.X.class(X)
    
    return(X)
}
x2base.data.frame.list <- function(covs.list, ...) {
    A <- list(...)
    X.names <- c("covs.list",
                 "treat.list",
                 "weights",
                 "distance.list",
                 "addl.list",
                 "s.d.denom",
                 "call",
                 "cluster",
                 "imp",
                 "s.weights",
                 "focal",
                 "discarded",
                 "method",
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    covs.list <- covs.list
    treat.list <- A$treat.list
    data <- A$data
    weights <- A$weights
    distance.list <- A$distance.list
    cluster <- A$cluster
    addl.list <- A$addl.list
    s.d.denom <- A$s.d.denom
    method <- A$method
    imp <- A$imp
    s.weights <- A$s.weights
    subset <- A$subset
    focal <- A$focal
    ntimes <- length(covs.list)
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (is_not_null(s.weights) && any(vapply(s.weights, anyNA, logical(1L)))) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    #Checks
    if (is_null(covs.list)) {
        stop("covs.list must be specified.", call. = FALSE)
    }
    if (!is.list(covs.list)) {
        stop("covs.list must be a list of covariates for which balanced is to be assessed at each time point.", call. = FALSE)
    }
    if (any(!vapply(covs.list, is.data.frame, logical(1L)))) {
        stop("Each item in covs.list must be a data frame.", call. = FALSE)
    }
    
    if (is_not_null(data) && !is.data.frame(data)) {
        warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execuction will halt.")
        data <- NULL
    }
    
    specified <- setNames(rep(FALSE, 1), "weights")
    if (is_not_null(weights)) {
        if (!is.character(weights) && !is.numeric(weights) && !is.data.frame(weights) && !is.list(weights)) {
            stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
        }
        specified["weights"] <- TRUE
    }
    
    #Getting method
    if (is_null(method)) {
        if (specified["weights"]) {
            
            #message("Assuming \"weighting\". If not, specify with an argument to method.")
            
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
    
    if (is_not_null(cluster) && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
        stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (is_not_null(imp) && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
        stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
    }
    
    #Order covs.list
    all.covs <- unique(unlist(lapply(covs.list, names)))
    covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
    
    #Process treat
    if (is_null(treat.list)) stop("treat.list must be specified.", call. = FALSE)
    if (!is.vector(treat.list)) {
        treat.list <- as.list(treat.list)
    }
    if (length(treat.list) != length(covs.list)) {
        stop("treat.list must be a list of treatment statuses at each time point.", call. = FALSE)
    }
    
    for (ti in seq_along(treat.list)) {
        if (is.numeric(treat.list[[ti]]) || is.factor(treat.list[[ti]]) || (is.character(treat.list[[ti]]) && length(treat.list[[ti]]) > 1)) {
            #treat.list[[ti]] <- treat.list[[ti]]
        }
        else if (is.character(treat.list[[ti]]) && length(treat.list[[ti]])==1 && any(names(data) == treat.list[[ti]])) {
            names(treat.list)[ti] <- treat.list[[ti]]
            treat.list[[ti]] <- data[[treat.list[[ti]]]]
        }
        else stop("Each item in treat.list must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
        
        if (sum(is.na(treat.list[[ti]])) > 0)
            stop("Missing values exist in treat.list", call. = FALSE)
        
        if (is_binary(treat.list[[ti]])) {
            treat.list[[ti]] <- binarize(treat.list[[ti]])
        }
        else if (is.character(treat.list[[ti]])) {
            treat.list[[ti]] <- factor(treat.list[[ti]])
        }
    }
    
    #Process weights
    for (i in c("weights")) {
        assign(i, data.frame.process(i, A[[i]], do.call("cbind", treat.list), do.call("cbind", covs.list), data))
    }
    
    #Process addl and distance
    for (i in c("addl.list", "distance.list")) {
        assign(i, list.process(i, A[[i]], ntimes, 
                               "covs.list",
                               treat.list,
                               covs.list,
                               data
        ))
    }
    
    #Process sampling weights
    if (is_not_null(s.weights)) {
        if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
            stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
        }
        if (is.character(s.weights) && length(s.weights)==1) {
            if (any(names(data) == s.weights)) {
                s.weights <- data[[s.weights]]
            }
            else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
        }
        if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
        
    }
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
            cluster <- data[[cluster]]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("s.weights", "cluster", "subset")
    data.frames <- c("weights")
    lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is_null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L)),
                          vapply(lists, function(x) {
                              if (is_null(get0(x))) 0 
                              else if (is.vector(get(x))) {
                                  if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                                  else max(lengths(get(x)))
                              }
                              else max(vapply(get(x), function(y) if (is_not_null(y)) {if (is_not_null(nrow(y))) nrow(y) else length(y)} else 0, numeric(1L)))
                          }, numeric(1L))), c(vectors, data.frames, lists))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames, lists)) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs.list."), call. = FALSE)
    }
    
    if (is_not_null(weights)) {
        if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
        if (any(vapply(weights, function(x) any(!is.finite(x)), logical(1L)))) stop("All weights must be numeric.", call. = FALSE)
        if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
        
        if (length(method) == 1) {
            method <- rep(method, ncol(weights))
        }
        else if (length(method) != ncol(weights)) {
            stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
        }
        
    }
    
    #Get s.d.denom
    X$s.d.denom <- rep("pooled", max(1, ncol(weights)))
    
    if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    if (is_null(s.weights)) s.weights <- rep(1, length(treat.list[[1]]))
    
    X$method <- method
    X$covs.list <- lapply(covs.list, function(x) x)
    X$weights <- weights
    X$treat.list <- lapply(treat.list, function(x) x)
    X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
    X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
    X$cluster <- factor(cluster)
    X$call <- NULL
    X$imp <- factor(imp)
    X$s.weights <- s.weights
    
    if (is_null(subset)) subset <- rep(TRUE, length(treat.list[[1]]))
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- get.X.class(X)
    
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
    
    X.names <- c("covs.list",
                 "treat.list",
                 "weights",
                 "distance.list",
                 "addl.list",
                 "s.d.denom",
                 "call",
                 "cluster",
                 "imp",
                 "s.weights",
                 "focal",
                 "discarded",
                 "method",
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    ID <- sort(unique(cbmsm$id))
    times <- sort(unique(cbmsm$time))
    treat.list <- lapply(times, function(x) cbmsm$treat.hist[ID, x]) 
    covs <- cbmsm$data[names(cbmsm$data) %in% attributes(terms(cbmsm$model))$term.labels]
    weights <- data.frame(weights = get.w(cbmsm)[ID])
    ntimes <- length(times)
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    
    covs.list <- vector("list", ntimes)
    for (ti in times) {
        if (ti == 1) {
            covs.list[[ti]] <- setNames(data.frame(covs[cbmsm$time == ti, , drop = FALSE][ID, , drop = FALSE]),
                                        paste0(names(covs), "0"))
        }
        else {
            covs.list[[ti]] <- cbind(covs.list[[ti - 1]], 
                                     setNames(data.frame(cbmsm$y[cbmsm$time == ti - 1][ID], 
                                                         covs[cbmsm$time == ti, , drop = FALSE][ID, , drop = FALSE]),
                                              c(paste0("treat", ti - 1), paste0(names(covs), ti))))
        }
    }
    
    cluster <- A$cluster
    subset <- A$subset
    data <- A$data
    
    cbmsm.data <- cbmsm$data[cbmsm$time == 1, , drop = FALSE][ID, , drop = FALSE]
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (!is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
            stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
        }
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(cbmsm.data) == cluster)) {
                cluster <- cbmsm.data[[cluster]]
            }
            else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
        }
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl.list", "distance.list")) {
        assign(i, list.process(i, A[[i]], ntimes, 
                               "the original call to CBMSM()",
                               treat.list,
                               covs.list,
                               data,
                               cbmsm.data))
    }
    
    if (is_not_null(distance.list)) distance.list <- lapply(times, function(x) data.frame(distance.list[[x]], prop.score = cbmsm$fitted.values))
    else if (is_not_null(cbmsm$fitted.values)) distance.list <- lapply(times, function(x) data.frame(prop.score = cbmsm$fitted.values))
    else distance.list <- NULL
    
    ensure.equal.lengths <- TRUE
    vectors <- c("cluster", "subset")
    data.frames <- c("weights")
    lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L)),
                          vapply(lists, function(x) {
                              if (is.null(get0(x))) 0 
                              else if (is.vector(get(x))) {
                                  if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                                  else max(lengths(get(x)))
                              }
                              else max(vapply(get(x), function(y) if (is_not_null(y)) nrow(y) else 0, numeric(1L)))
                          }, numeric(1L))), c(vectors, data.frames, lists))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames, lists)) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to CBMSM()."), call. = FALSE)
    }
    
    if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$weights <- weights
    X$treat.list <- lapply(treat.list, function(x) x)
    X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
    X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
    X$covs.list <- lapply(covs.list, function(x) x)
    X$call <- cbmsm$call
    X$cluster <- factor(cluster)
    X$method <- rep("weighting", ncol(weights))
    X$s.weights <- NULL
    X$s.d.denom <- "pooled"
    
    if (is_null(subset)) subset <- rep(TRUE, lengths["treat.list"])
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- get.X.class(X)
    
    return(X)
}
x2base.weightitMSM <- function(weightitMSM, ...) {
    A <- list(...)
    
    X.names <- c("covs.list",
                 "treat.list",
                 "weights",
                 "distance.list",
                 "addl.list",
                 "s.d.denom",
                 "call",
                 "cluster",
                 "imp",
                 "s.weights",
                 "focal",
                 "discarded",
                 "method",
                 "subclass")
    X <- setNames(vector("list", length(X.names)),
                  X.names)
    
    #Initializing variables
    estimand <- weightitMSM$estimand
    weights <- data.frame(weights = get.w(weightitMSM))
    treat.list <- weightitMSM$treat.list
    covs.list <- weightitMSM$covs.list
    if (is_null(covs.list)) stop("No covariates were specified in the weightit object.", call. = FALSE)
    covs <- do.call("cbind", covs.list)
    s.weights <- weightitMSM$s.weights
    data <- A$data
    cluster <- A$cluster
    imp <- A$imp
    subset <- A$subset
    ntimes <- length(treat.list)
    
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    weightitMSM.data <- weightitMSM$data
    d.e.in.w <- vapply(c("covs.list", "exact", "by"), function(x) is_not_null(weightitMSM[[x]]), logical(1L))
    if (any(d.e.in.w)) weightitMSM.data <- do.call("cbind", c(list(covs), weightitMSM[c("exact", "by")])[d.e.in.w])
    else weightitMSM.data <- NULL
    
    if (all(vapply(treat.list, has.treat.type, logical(1L)))) {
        treat.type <- vapply(treat.list, get.treat.type, character(1L))
    }
    else if (length(weightitMSM$treat.type) == length(treat.list)) {
        treat.type <- weightitMSM$treat.type
    }
    else {
        treat.type <- vapply(treat.list, function(treat) {
            get.treat.type(assign.treat.type(treat))
        }, character(1L))
    } 
    
    #Order covs.list
    all.covs <- unique(unlist(lapply(covs.list, names)))
    covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (!is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
            stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
        }
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else if (any(names(weightitMSM.data) == cluster)) {
                cluster <- weightitMSM.data[[cluster]]
            }
            else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
        }
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process addl and distance
    for (i in c("addl.list", "distance.list")) {
        assign(i, list.process(i, A[[i]], ntimes, 
                               "the original call to weightitMSM()",
                               treat.list,
                               covs.list,
                               data,
                               weightitMSM.data))
    }
    
    if (is_not_null(distance.list)) distance.list <- lapply(seq_along(distance.list), function(x) data.frame(distance.list[[x]], prop.score = weightitMSM$ps.list[[x]]))
    else if (is_not_null(weightitMSM$ps.list)) distance.list <- lapply(seq_along(weightitMSM$ps.list), function(x) data.frame(prop.score = weightitMSM$ps.list[[x]]))
    else distance.list <- NULL
    
    ensure.equal.lengths <- TRUE
    vectors <- c("s.weights", "cluster", "subset")
    data.frames <- c("weights")
    lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L)),
                          vapply(lists, function(x) {
                              if (is.null(get0(x))) 0 
                              else if (is.vector(get(x))) {
                                  if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                                  else max(lengths(get(x)))
                              }
                              else max(vapply(get(x), function(y) if (is_not_null(y)) nrow(y) else 0, numeric(1L)))
                          }, numeric(1L))), c(vectors, data.frames, lists))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames, lists)) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to weightitMSM()."), call. = FALSE)
    }
    
    #Get s.d.denom
    X$s.d.denom <- rep("pooled", ncol(weights))
    
    if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    X$weights <- weights
    X$treat.list <- lapply(treat.list, function(x) x)
    X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
    X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
    X$covs.list <- lapply(covs.list, function(x) x)
    X$call <- NULL
    X$cluster <- factor(cluster)
    X$method <- rep("weighting", ncol(weights))
    X$s.weights <- s.weights
    
    if (is_null(subset)) subset <- rep(TRUE, lengths["treat.list"])
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- get.X.class(X)
    
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
    if (is_not_null(treat.list)) {
        if (!all(sapply(treat.list, function(x) any(vapply(Q[["treat"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
            treat.list <- A[["treat.list"]]
        }
        msm <- TRUE
    }
    
    #covs 
    if (is_not_null(covs)) covs <- as.data.frame(covs)
    
    #covs.list
    if (is_not_null(covs.list)) {
        if (!all(sapply(covs.list, function(x) any(vapply(Q[["covs"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
            covs.list <- A[["covs.list"]]
        }
        msm <- TRUE
    }
    
    #formula OK
    
    #formula.list
    if (is_not_null(formula.list)) {
        if (!all(sapply(formula.list, function(x) any(vapply(Q[["formula"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
            formula.list <- A[["formula.list"]]
        }
        msm <- TRUE
    }
    
    #data
    if (is_not_null(data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if ("imp" %nin% names(A)) A[["imp"]] <- data[[".imp"]]
        }
        data <- as.data.frame(data)
    }
    
    #weights
    if (is_not_null(weights)) {
        if (is.vector(weights, "numeric")) weights <- data.frame(weights = weights)
        else weights <- as.data.frame(weights)
    }

    #distance
    if (is_not_null(distance)) {
        if (is.numeric(distance)) {
            if (is_not_null(attr(distance, "name"))) distance <- setNames(data.frame(distance),
                                                                          attr(distance, "name"))
            else distance <- data.frame(distance = distance)
        }
        else distance <- as.data.frame(distance)
    }
    
    #distance.list
    if (is_not_null(distance.list)) {
        if (!all(sapply(distance.list, function(x) any(vapply(Q[["distance"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
            distance.list <- A[["distance.list"]]
        }
        #msm <- TRUE
    }
    
    #subclass
    if (is_not_null(subclass)) subclass <- factor(subclass)
    
    #match.strata
    if (is_not_null(match.strata)) match.strata <- factor(match.strata)
    
    #estimand
    if (is_not_null(estimand)) {
        estimand.name <- attr(estimand, "name")
        if (is_not_null(estimand.name) && toupper(estimand.name) == "ATT") {
            if (estimand == 0) estimand <- "ATE"
            else estimand <- "ATT"
        }
        else if (is_not_null(estimand.name) && toupper(estimand.name) == "ATE") {
            if (estimand == 0) estimand <- "ATT"
            else estimand <- "ATE"
        }
        else {
            if (tolower(estimand) %in% c("att", "treat", "treated", "tr", "t", "atet")) estimand <- "ATT"
            else if (tolower(estimand) %in% c("ate", "all")) estimand <- "ATE"
            else if (tolower(estimand) %in% c("atc", "control", "untreated", "u", "c", "ctrl", "atu", "atec", "ateu")) estimand <- "ATC"
            else estimand <- NULL
        }
    }
    
    #s.weights OK
    
    #focal OK
    
    #call OK
    
    #model (only to extract data)
    if (is_not_null(obj[["model"]])) {
        if (is_null(data) && "data" %in% names(obj[["model"]])) {
            data <- obj[["model"]][["data"]]
        }
    }
    
    if (!msm) {
        
        cluster <- A$cluster
        addl <- A[["addl"]]
        s.d.denom <- A$s.d.denom
        method <- A$method
        imp <- A$imp
        subset <- A$subset
        
        # if (length(distance) > 0 && !is.character(distance) && !is.numeric(distance) && !is.data.frame(distance)) {
        #     stop("The argument to distance must be a vector of distance scores or the (quoted) name of a variable in data that contains distance scores.", call. = FALSE)
        # }
        
        if (is_not_null(cluster) && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
            stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
        }
        if (is_not_null(imp) && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
            stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
        }
        
        t.c <- use.tc.fd(formula, data, treat, covs)
        
        treat <- t.c[["treat"]]
        covs  <- t.c[["covs"]]
        
        #Checks
        if (is_null(covs)) {
            stop("covs data.frame must be specified.", call. = FALSE)
        }
        if (!is.data.frame(covs)) {
            stop("covs must be a data.frame.", call. = FALSE)
        }
        if (is_not_null(data) && !is.data.frame(data)) {
            warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
        
        #Process treat
        if (is_null(treat)) stop("treat must be specified.", call. = FALSE)
        
        if (is.numeric(treat) || is.factor(treat) || (is.character(treat) && length(treat) > 1)) {
            treat <- treat
        }
        else if (is.character(treat) && length(treat)==1 && any(names(data) == treat)) {
            treat <- data[[treat]]
        }
        else stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
        
        if (sum(is.na(treat)) > 0)
            stop("Missing values exist in treat.", call. = FALSE)
        
        if (nunique(treat) == 2) {
            treat <- binarize(treat)
        }
        else if (is.character(treat)) {
            treat <- factor(treat)
        }
        
        #Process weights, addl, and distance
        for (i in c("weights", "addl", "distance")) {
            assign(i, data.frame.process(i, A[[i]], treat, covs, data))
        }
        
        #Process subclass
        if (is_not_null(subclass)) {
            if (is.numeric(subclass) || is.factor(subclass) || (is.character(subclass) && length(subclass) > 1)) {
                subclass <- factor(subclass)
            }
            else if (is.character(subclass) && length(subclass)==1 && any(names(data) == subclass)) {
                subclass <- factor(data[[subclass]])
            }
            else stop("The name supplied to subclass is not the name of a variable in data.", call. = FALSE)
            
            subclass.weights <- match.strata2weights(match.strata = subclass,
                                                     treat = treat,
                                                     covs = covs)
        }
        
        #Process match.strata
        if (is_not_null(match.strata)) {
            if (is.numeric(match.strata) || is.factor(match.strata) || (is.character(match.strata) && length(match.strata) > 1)) {
                match.strata <- factor(match.strata)
            }
            else if (is.character(match.strata) && length(match.strata)==1 && any(names(data) == match.strata)) {
                match.strata <- factor(data[[match.strata]])
            }
            else stop("The name supplied to match.strata is not the name of a variable in data.", call. = FALSE)
            
            match.strata.weights <- match.strata2weights(match.strata = match.strata,
                                                         treat = treat,
                                                         covs = covs)
        }
        
        specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
        if (is_not_null(weights)) {
            if (!is.character(weights) && !is.numeric(weights) && !is.data.frame(weights) && !is.list(weights)) {
                stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
            }
            specified["weights"] <- TRUE
        }
        if (is_not_null(subclass)){
            if (!is.character(subclass) && !is.factor(subclass) && !is.numeric(subclass)) {
                stop("The argument to subclass must be a vector of subclass membership or the (quoted) name of a variable in data that contains subclass membership.", call. = FALSE)
            }
            specified["subclass"] <- TRUE
        }
        if (is_not_null(match.strata)) {
            if (!is.character(match.strata) && !is.numeric(match.strata) && !is.factor(match.strata)) {
                stop("The argument to match.strata must be a vector of match stratum membership or the (quoted) name of a variable in data that contains match stratum membership.", call. = FALSE)
            }
            specified["match.strata"] <- TRUE
        }
        
        #Getting method
        if (is_null(method)) {
            if (specified["match.strata"]) {
                if (specified["subclass"]) {
                    if (isTRUE(all.equal(match.strata.weights, subclass.weights, 
                                         check.attributes = FALSE))) {
                        subclass.weights <- subclass <- NULL
                        specified["subclass"] <- FALSE
                    }
                }
                if (specified["weights"]) {
                    if (ncol(weights) == 1 && isTRUE(all.equal(match.strata.weights, weights, 
                                                               check.attributes = FALSE))) {
                        weights <- NULL
                        specified["weights"] <- FALSE
                    }
                }
                method <- "matching"
            }
            else if (specified["subclass"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    weights <- match.strata <- NULL
                }
                method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word_list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    match.strata <- subclass <- NULL
                }
                else {
                    if (!any(c("optweight", "weightit") %in% class(obj))) {
                        message("Assuming \"weighting\". If not, specify with an argument to method.")}
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
                        match.strata <- subclass <- NULL
                    }
                    method <- "weighting"
                }
                else if (specified["match.strata"]) {
                    message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
                    subclass <- NULL
                    method <- "matching"
                }
                else if (specified["subclass"]) {
                    message("method = \"weighting\" is specified, but no weights are present. Assuming \"subclassification\" and using subclass instead.")
                    method <- "subclassification"
                    #weights <- rep(1, nrow(covs))
                }
                else {
                    method <- "matching"
                }
            }
            else if (specified.method == "matching") {
                if (specified["match.strata"]) {
                    if (sum(specified) > 1) {
                        message(word_list(names(specified)[specified]), " are specified. Using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                        weights <- subclass <- NULL
                    }
                    method <- "matching"
                }
                else if (specified["weights"]) {
                    if (sum(specified) > 1) {
                        message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
                        match.strata <- subclass <- NULL
                    }
                    method <- "matching"
                }
                else if (specified["subclass"]) {
                    message("method = \"matching\" is specified, but no weights or match.strata are present. Assuming \"subclassification\" and using subclass instead.")
                    method <- "subclassification"
                    #weights <- rep(1, nrow(covs))
                }
                else {
                    method <- "matching"
                }
            }
            else if (specified.method == "subclassification") {
                if (specified["subclass"]) {
                    if (sum(specified) > 1) {
                        message(word_list(names(specified)[specified]), " are specified. Using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
                        weights <- match.strata <- NULL
                    }
                    method <- "subclassification"
                    #weights <- rep(1, nrow(covs))
                }
                else if (specified["match.strata"]) {
                    message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
                    weights <- NULL
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
                stop("Only weights can be specified with mutiple methods.", call. = FALSE)
            }
            else if (!specified["weights"]) {
                warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
                method <- "matching"
            }
            else {
                #Matching and/or weighting with various weights
                method <- specified.method
                match.strata <- subclass <- NULL
            }
        }
        
        #Process sampling weights
        if (is_not_null(s.weights)) {
            if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
                stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
            }
            if (is.character(s.weights) && length(s.weights)==1) {
                if (any(names(data) == s.weights)) {
                    s.weights <- data[[s.weights]]
                }
                else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
            }
            if (anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
        }
        
        #Process cluster
        if (is_not_null(cluster)) {
            if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
                cluster <- cluster
            }
            else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
        }
        
        #Process subset
        if (is_not_null(subset)) {
            if (!is.logical(subset)) {
                stop("The argument to subset must be a logical vector.", call. = FALSE)
            }
            if (anyNA(subset)) {
                warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
                subset[is.na(subset)] <- FALSE
            }
        }
        
        ensure.equal.lengths <- TRUE
        vectors <- c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset")
        data.frames <- c("covs", "weights", "distance", "addl")
        problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
        lengths <- setNames(c(lengths(mget(vectors, ifnotfound = list(NULL))), 
                              vapply(data.frames, 
                                     function(x) NROW(get0(x)), 
                                     numeric(1L))), c(vectors, data.frames))
        #Process imp
        if (is_not_null(imp)) {
            if (is_(imp, c("numeric", "factor")) || (is.character(imp) && length(imp) > 1)) {
                imp <- imp
            }
            else if (is.character(imp) && length(imp)==1 && any(names(data) == imp)) {
                imp <- data[[imp]]
            }
            else stop("The name supplied to imp is not the name of a variable in data.", call. = FALSE)
            
            imp.lengths <- vapply(unique(imp), function(i) sum(imp == i), numeric(1L))
            
            if (all_the_same(imp.lengths)) { #all the same
                for (i in vectors) {
                    if (lengths[i] > 0 && lengths[i] != length(imp)) { 
                        if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                        if (lengths[i] == imp.lengths[1]) {
                            temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                                   order2 = seq_along(imp))
                            
                            temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                                   get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
                            )
                            temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                                by.y = 1:2, sort = FALSE)
                            assign(i, temp.merge[[4]][order(temp.merge[[3]])])
                        }
                        else {
                            problematic[i] <- TRUE
                        }
                    }
                }
                for (i in data.frames) {
                    if (lengths[i] > 0 && lengths[i] != length(imp)) {
                        if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                        if (lengths[i] == imp.lengths[1]) {
                            temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                                   order2 = seq_along(imp))
                            temp.var <- data.frame(sort(imp),rep(seq_len(lengths[i]), length(imp.lengths)),
                                                   get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
                            )
                            temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                                by.y = 1:2, sort = FALSE)
                            assign(i, setNames(temp.merge[order(temp.merge[[3]]), -c(1:3), drop = FALSE], names(get(i))))
                        }
                        else {
                            problematic[i] <- TRUE
                        }
                    }
                }
            }
            else {
                problematic <- lengths > 0 & lengths != length(imp)
            }
            if (any(problematic)) {
                stop(paste0(word_list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
            }
            else ensure.equal.lengths <- FALSE
        }
        
        #Ensure all input lengths are the same.
        if (ensure.equal.lengths) {
            for (i in c(vectors, data.frames[data.frames!="covs"])) {
                if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                    problematic[i] <- TRUE
                }
            }
        }
        if (any(problematic)) {
            stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
        }
        
        #Turn match.strata into weights
        if (is_not_null(get0("match.strata.weights"))) {
            weights <- data.frame(weights = match.strata.weights)
        }
        
        if (is_not_null(weights)) {
            if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
            if (any(vapply(weights, function(x) !is.numeric(x), logical(1L)))) {
                stop("All weights must be numeric.", call. = FALSE)
            }
            if (length(method) == 1) {
                method <- rep(method, ncol(weights))
            }
            else if (length(method) != ncol(weights)) {
                stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
            }
            
        }
        
        #Check focal
        if (is_not_null(focal) && is.factor(treat)) {
            if (is.numeric(focal)) {
                if (focal <= nunique(treat)) focal <- levels(treat)[focal]
                else 
                    stop(paste0("focal was specified as ", focal, 
                                ", but there are only ", levels(treat), " treatment groups."), call. = FALSE)
            }
            else {
                if (!any(levels(treat) == focal)) 
                    stop(paste0("The name specified to focal is not the name of any treatment group."), call. = FALSE)
            }
        }
        
        #Get s.d.denom
        if (is_binary(treat) || !is.numeric(treat)) { #non-continuous
            X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal, method)
        }
        
        if (any(c(is.na(covs), is.na(addl)))) {
            warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
        }
        
        if (is_null(subset)) subset <- rep(TRUE, length(treat))
        
        X$method <- method
        X$covs <- covs
        X$weights <- weights
        X$treat <- treat
        X$distance <- distance
        X$subclass <- subclass
        X$cluster <- factor(cluster)
        X$call <- call
        X$addl <- addl
        X$imp <- factor(imp)
        X$s.weights <- s.weights
        X$focal <- focal
        
    }
    else {
        cluster <- A$cluster
        addl.list <- A$addl.list
        s.d.denom <- A$s.d.denom
        method <- A$method
        imp <- A$imp
        subset <- A$subset
        
        if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
        if (is_not_null(s.weights) && any(vapply(s.weights, anyNA, logical(1L)))) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
        
        initial.list.lengths <- c(length(formula.list), length(covs.list), length(treat.list))
        if (!all_the_same(initial.list.lengths[initial.list.lengths != 0])) stop("The lists in the object were not the same length.", call. = FALSE)
        ntimes.guess <- max(initial.list.lengths)
        
        if (is_null(treat.list)) treat.list <- vector("list", length(ntimes.guess)) 
        if (is_null(covs.list)) covs.list <- vector("list", length(ntimes.guess)) 
        for (i in seq_len(ntimes.guess)) {
            t.c <- use.tc.fd(formula.list[[i]], data, treat.list[[i]], covs.list[[i]])
            
            treat.list[[i]] <- t.c[["treat"]]
            covs.list[[i]]  <- t.c[["covs"]]
            if (is_not_null(t.c[["treat.name"]])) names(treat.list)[i] <- t.c[["treat.name"]]
        }
        
        ntimes <- length(covs.list)
        
        #Checks
        if (is_null(covs.list)) {
            stop("A covariate list must be specified.", call. = FALSE)
        }
        if (!is.list(covs.list)) {
            stop("covs.list must be a list of covariates for which balanced is to be assessed at each time point.", call. = FALSE)
        }
        if (any(!vapply(covs.list, function(x) is.data.frame(x), logical(1L)))) {
            stop("Each item in covs.list must be a data frame.", call. = FALSE)
        }
        
        if (is_not_null(data) && !is.data.frame(data)) {
            warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execuction will halt.")
            data <- NULL
        }
        
        specified <- setNames(rep(FALSE, 1), "weights")
        if (is_not_null(weights)) {
            if (!is.character(weights) && !is.numeric(weights) && !is.data.frame(weights) && !is.list(weights)) {
                stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
            }
            specified["weights"] <- TRUE
        }
        
        #Getting method
        if (is_null(method)) {
            if (specified["weights"]) {
                
                #message("Assuming \"weighting\". If not, specify with an argument to method.")
                
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
        
        if (is_not_null(cluster) && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
            stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
        }
        if (is_not_null(imp) && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
            stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
        }
        
        #Order covs.list
        all.covs <- unique(unlist(lapply(covs.list, names)))
        covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
        
        #Process treat
        if (is_null(treat.list)) stop("A treatment list must be specified.", call. = FALSE)
        if (!is.vector(treat.list)) {
            treat.list <- as.list(treat.list)
        }
        if (length(treat.list) != length(covs.list)) {
            stop("treat.list must be a list of treatment statuses at each time point.", call. = FALSE)
        }
        
        for (ti in seq_along(treat.list)) {
            if (is.numeric(treat.list[[ti]]) || is.factor(treat.list[[ti]]) || (is.character(treat.list[[ti]]) && length(treat.list[[ti]]) > 1)) {
                #treat.list[[ti]] <- treat.list[[ti]]
            }
            else if (is.character(treat.list[[ti]]) && length(treat.list[[ti]])==1 && any(names(data) == treat.list[[ti]])) {
                names(treat.list)[ti] <- treat.list[[ti]]
                treat.list[[ti]] <- data[[treat.list[[ti]]]]
            }
            else stop("Each item in treat.list must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
            
            if (sum(is.na(treat.list[[ti]])) > 0)
                stop("Missing values exist in treat.list", call. = FALSE)
            
            if (is_binary(treat.list[[ti]])) {
                treat.list[[ti]] <- binarize(treat.list[[ti]])
            }
            else if (is.character(treat.list[[ti]])) {
                treat.list[[ti]] <- factor(treat.list[[ti]])
            }
        }
        
        #Process weights
        for (i in c("weights")) {
            assign(i, data.frame.process(i, A[[i]], do.call("cbind", treat.list), do.call("cbind", covs.list), data))
        }
        
        #Process addl and distance
        for (i in c("addl.list", "distance.list")) {
            assign(i, list.process(i, A[[i]], ntimes, 
                                   "covs.list",
                                   treat.list,
                                   covs.list,
                                   data
            ))
        }
        
        #Process sampling weights
        if (is_not_null(s.weights)) {
            if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
                stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
            }
            if (is.character(s.weights) && length(s.weights)==1) {
                if (any(names(data) == s.weights)) {
                    s.weights <- data[[s.weights]]
                }
                else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
            }
        }
        
        #Process cluster
        if (is_not_null(cluster)) {
            if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
                cluster <- cluster
            }
            else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
                cluster <- data[[cluster]]
            }
            else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
        }
        
        #Process subset
        if (is_not_null(subset)) {
            if (!is.logical(subset)) {
                stop("The argument to subset must be a logical vector.", call. = FALSE)
            }
            if (anyNA(subset)) {
                warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
                subset[is.na(subset)] <- FALSE
            }
        }
        
        ensure.equal.lengths <- TRUE
        vectors <- c("s.weights", "cluster", "subset")
        data.frames <- c("weights")
        lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
        problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
        lengths <- setNames(c(lengths(mget(vectors)), 
                              vapply(data.frames, 
                                     function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                     }, numeric(1L)),
                              vapply(lists, function(x) {
                                  if (is.null(get0(x))) 0 
                                  else if (is.vector(get(x))) {
                                      if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                                      else max(lengths(get(x)))
                                  }
                                  else max(vapply(get(x), function(y) if (is_not_null(y)) {if (is_not_null(nrow(y))) nrow(y) else length(y)} else 0, numeric(1L)))
                              }, numeric(1L))), c(vectors, data.frames, lists))
        
        #Ensure all input lengths are the same.
        if (ensure.equal.lengths) {
            for (i in c(vectors, data.frames, lists)) {
                if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
                    problematic[i] <- TRUE
                }
            }
        }
        if (any(problematic)) {
            stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs.list."), call. = FALSE)
        }
        
        if (is_not_null(weights)) {
            if (any(vapply(weights, function(x) any(!is.finite(x)), logical(1L)))) {
                stop("All weights must be numeric.", call. = FALSE)
            }
            if (length(method) == 1) {
                method <- rep(method, ncol(weights))
            }
            else if (length(method) != ncol(weights)) {
                stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
            }
            
        }
        
        #Get s.d.denom
        X$s.d.denom <- rep("pooled", max(1, ncol(weights)))
        
        if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
            warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
        }
        
        if (is_null(s.weights)) s.weights <- rep(1, length(treat.list[[1]]))
        if (is_null(subset)) subset <- rep(TRUE, length(treat.list[[1]]))
        
        X$method <- method
        X$covs.list <- lapply(covs.list, function(x) x)
        X$weights <- weights
        X$treat.list <- lapply(treat.list, function(x) x)
        X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
        X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
        X$cluster <- factor(cluster)
        X$call <- call
        X$imp <- factor(imp)
        X$s.weights <- s.weights
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[X.names], X.names)
    
    class(X) <- get.X.class(X)
    
    return(X)
}
