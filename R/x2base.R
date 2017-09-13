#Functions to convert object to base.bal.tab input

x2base <- function(obj, ...) UseMethod("x2base")

x2base.matchit <- function(m, ...) {
    A <- list(...)
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              subclass=NA,
              method=NA,
              addl=NA,
              distance=NA,
              call=NA,
              cluster=NA,
              discarded=NA)
    
    #Initializing variables
    
    if (any(class(m) == "matchit.subclass")) {
        X$subclass <- factor(m$subclass)
        X$method <- "subclassification"
    }
    else if (any(class(m) == "matchit.full")) {
        X$subclass <- NULL
        X$method <- "weighting"
    }
    else {
        X$subclass <- NULL
        X$method <- "matching"
    }
    weights <- data.frame(weights = m$weights)
    treat <- m$treat
    
    if (length(m$model$model) > 0) {
        o.data <- m$model$model #data used in the PS formula, including treatment and covs
        covs <- data.frame(o.data[, !is.na(match(names(o.data), attributes(terms(m$model))$term.labels))])
        if (identical(o.data, data)) o.data <- NULL
    }
    else {
        o.data <- NULL
        covs <- data.frame(m$X)
    }
    m.data <- m$model$data
    
    #Process cluster
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (cluster %in% names(data)) {
                cluster <- data[, cluster]
            }
            else if (cluster %in% names(m.data)) {
                cluster <- m.data[, cluster]
            }
        }
        else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.vector(val, mode = "list")) {
                val.list <- lapply(val, function(x) process.val(x, i, treat, covs, data))
                val.list <- lapply(seq_along(val.list), function(x) {
                    if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                    val.list[[x]]})
                if (length(unique(sapply(val.list, nrow))) > 1) {
                    stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                }
                
                val.df <- setNames(do.call("cbind", val.list),
                                   c(sapply(val.list, names)))
            }
            else {
                val.df <- process.val(val, i, treat, covs, data)
            }
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    if (any(is.finite(m$distance))) {
        if (length(distance) > 0) distance <- cbind(distance, distance = m$distance)
        else distance <- data.frame(distance = m$distance)
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("cluster", "treat")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word.list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to matchit()."), call. = FALSE)
    }
    
    X$treat <- treat
    X$weights <- weights
    X$discarded <- m$discarded
    X$covs <- covs
    X$distance <- distance
    X$addl <- addl
    X$cluster <- factor(cluster)
    X$call <- m$call
    return(X)
}
x2base.ps <- function(ps, ...) {
    #stop.method
    #s.d.denom
    A <- list(...)
    
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              distance=NA,
              addl=NA,
              s.d.denom=NA,
              call=NA,
              cluster = NA,
              s.weights = NA)
    
    #Initializing variables
    if (length(A) > 0 && names(A)[1]=="" && length(A$stop.method)==0) A$stop.method <- A[[1]]
    if (length(A$stop.method) == 0 && length(A$full.stop.method) > 0) A$stop.method <- A$full.stop.method
    
    if (length(A$stop.method) > 0) {
        if (is.character(A$stop.method)) {
            rule1 <- names(ps$w)[apply(sapply(tolower(A$stop.method), function(x) startsWith(tolower(names(ps$w)), x)), 1, any)]
            if (length(rule1) == 0) {
                message(paste0("Warning: stop.method should be ", word.list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- names(ps$w)
            }
        }
        else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(names(ps$w)))) {
            if (any(!A$stop.method %in% seq_along(names(ps$w)))) {
                message(paste0("Warning: There are ", length(names(ps$w)), " stop methods available, but you requested ", 
                               word.list(A$stop.method[!A$stop.method %in% seq_along(names(ps$w))], and.or = "and"),"."))
            }
            rule1 <- names(ps$w)[A$stop.method %in% seq_along(names(ps$w))]
        }
        else {
            warning("stop.method should be ", word.list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- names(ps$w)
        }
    }
    else {
        rule1 <- names(ps$w)
    }

    s <- names(ps$w)[match(tolower(rule1), tolower(names(ps$w)))]
    estimand <- substr(tolower(s), nchar(s)-2, nchar(s))
    
    if (length(A$s.d.denom>0) && is.character(A$s.d.denom)) {
        X$s.d.denom <- tryCatch(match.arg(A$s.d.denom, c("treated", "control", "pooled")),
                                error = function(cond) {
                                    new.s.d.denom <- switch(substr(tolower(s), nchar(s)-2, nchar(s)), att = "treated", ate = "pooled")
                                    message(paste0("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\nUsing ", deparse(new.s.d.denom), " instead."))
                                    return(new.s.d.denom)})
    }
    else X$s.d.denom <- sapply(estimand, switch, att = "treated", ate = "pooled")
    
    weights <- get.w(ps, s, estimand)
    treat <- ps$treat
    covs <- ps$data[, ps$gbm.obj$var.names, drop = FALSE]
    
    ps.data <- ps$data
    
    #Process cluster
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (cluster %in% names(data)) {
                cluster <- data[, cluster]
            }
            else if (cluster %in% names(ps.data)) {
                cluster <- ps.data[, cluster]
            }
        }
        else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.vector(val, mode = "list")) {
                val.list <- lapply(val, function(x) process.val(x, i, treat, covs, data))
                val.list <- lapply(seq_along(val.list), function(x) {
                    if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                    val.list[[x]]})
                if (length(unique(sapply(val.list, nrow))) > 1) {
                    stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                }
                
                val.df <- setNames(do.call("cbind", val.list),
                                   c(sapply(val.list, names)))
            }
            else {
                val.df <- process.val(val, i, treat, covs, data)
            }
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    if (length(distance) > 0) distance <- cbind(distance, prop.score = ps$ps[s][, ])
    else distance <- data.frame(prop.score = ps$ps[s][, ])
    
    ensure.equal.lengths <- TRUE
    vectors <- c("cluster")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word.list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to ps()."), call. = FALSE)
    }
    
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$addl <- addl
    X$covs <- covs
    X$call <- ps$parameters
    X$cluster <- factor(cluster)
    X$method <- rep("weighting", ncol(weights))
    X$s.weights <- ps$sampw

    return(X)
}
x2base.mnps <- function(mnps, ...) {
    #stop.method
    #s.d.denom
    A <- list(...)
    
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              distance=NA,
              addl=NA,
              s.d.denom=NA,
              call=NA,
              cluster = NA,
              s.weights = NA)
    
    #Initializing variables
    if (length(A) > 0 && names(A)[1]=="" && length(A$stop.method)==0) A$stop.method <- A[[1]]
    if (length(A$stop.method) == 0 && length(A$full.stop.method) > 0) A$stop.method <- A$full.stop.method
    
    if (length(A$stop.method) > 0) {
        if (any(is.character(A$stop.method))) {
            rule1 <- mnps$stopMethods[sapply(t(sapply(tolower(A$stop.method), function(x) startsWith(tolower(mnps$stopMethods), x))), any)]
            if (length(rule1) == 0) {
                message(paste0("Warning: stop.method should be ", word.list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- mnps$stopMethods
            }
            # rule1 <- tryCatch(match.arg(tolower(stop.method), tolower(names(ps$w)), several.ok = TRUE),
            #                   error = function(cond) {message(paste0("Warning: stop.method should be ", word.list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."));
            #                       return(names(ps$w))})
        }
        else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(mnps$stopMethods))) {
            if (any(!A$stop.method %in% seq_along(mnps$stopMethods))) {
                message(paste0("Warning: There are ", length(mnps$stopMethods), " stop methods available, but you requested ", 
                               word.list(A$stop.method[!A$stop.method %in% seq_along(mnps$stopMethods)], and.or = "and"),"."))
            }
            rule1 <- mnps$stopMethods[A$stop.method %in% seq_along(mnps$stopMethods)]
        }
        else {
            warning("stop.method should be ", word.list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- mnps$stopMethods
        }
    }
    else {
        rule1 <- mnps$stopMethods
    }
    
    s <- mnps$stopMethods[match(tolower(rule1), tolower(mnps$stopMethods))]
    
    estimand <- setNames(mnps$estimand, s)
    
    if (length(A$s.d.denom>0) && is.character(A$s.d.denom)) {
        X$s.d.denom <- tryCatch(match.arg(A$s.d.denom, c("treated", "control", "pooled")),
                                error = function(cond) {
                                    new.s.d.denom <- switch(substr(tolower(s), nchar(s)-2, nchar(s)), att = "treated", ate = "pooled")
                                    message(paste0("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\nUsing ", deparse(new.s.d.denom), " instead."))
                                    return(new.s.d.denom)})
    }
    else X$s.d.denom <- sapply(estimand, switch, att = "treated", ate = "pooled")
    
    weights <- get.w(mnps, s)
    treat <- mnps$treatVar
    covs <- mnps$data[, mnps$psList[[1]]$gbm.obj$var.names, drop = FALSE]
    
    mnps.data <- mnps$data
    
    #Process cluster
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (cluster %in% names(data)) {
                cluster <- data[, cluster]
            }
            else if (cluster %in% names(mnps.data)) {
                cluster <- mnps.data[, cluster]
            }
        }
        else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.vector(val, mode = "list")) {
                val.list <- lapply(val, function(x) process.val(x, i, treat, covs, data))
                val.list <- lapply(seq_along(val.list), function(x) {
                    if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                    val.list[[x]]})
                if (length(unique(sapply(val.list, nrow))) > 1) {
                    stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                }
                
                val.df <- setNames(do.call("cbind", val.list),
                                   c(sapply(val.list, names)))
            }
            else {
                val.df <- process.val(val, i, treat, covs, data)
            }
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    # if (length(distance) > 0) distance <- cbind(distance, prop.score = ps$ps[s][, ])
    # else distance <- data.frame(prop.score = ps$ps[s][, ])
    
    ensure.equal.lengths <- TRUE
    vectors <- c("cluster")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word.list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to ps()."), call. = FALSE)
    }
    
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$addl <- addl
    X$covs <- covs
    X$call <- NULL
    X$cluster <- factor(cluster)
    X$method <- rep("weighting", ncol(weights))
    X$s.weights <- mnps$sampw
    
    return(X)
}
x2base.Match <- function(Match, ...) {
    #formula
    #data
    #treat
    #covs
    #addl
    #distance
    #s.d.denom
    A <- list(...)
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              method=NA,
              distance=NA,
              addl=NA,
              call=NA,
              s.d.denom=NA,
              cluster = NA)
    #Checks
    if (!is.list(Match) & length(Match) > 0) {
        stop("'Match' object contains no valid matches")}
    
    #Get treat and covs
    data <- A$data
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    
    if (sum(is.na(t.c$covs))>0)
        stop("Missing values exist in the covariates.", call. = FALSE)
    
    #Initializing variables
    m <- Match
    s <- m$estimand
    if (length(A$s.d.denom>0) && is.character(A$s.d.denom)) {
        X$s.d.denom <- tryCatch(match.arg(A$s.d.denom, c("treated", "control", "pooled")),
                                error = function(cond) {
                                    new.s.d.denom <- switch(toupper(s), ATT = "treated", ATE = "treated", ATC = "control")
                                    message(paste0("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\nUsing ", deparse(new.s.d.denom), " instead."))
                                    return(new.s.d.denom)})
    }
    else X$s.d.denom <- switch(toupper(s), ATT = "treated", ATE = "treated", ATC = "control")
    treat0 <- t.c$treat
    covs0  <- t.c$covs
    nobs <- nrow(covs0)
    
    #distance <- NULL
    
    data.list <- covs.list <- treat.list <- weights.list <- distance.list <- list(control=NA, treated=NA, unmatched=NA, dropped=NA)
    
    covs.list$control <- cbind(covs0[m$index.control, ], index=m$index.control)
    covs.list$treated <- cbind(covs0[m$index.treat, ], index=m$index.treat)
    covs.list$unmatched <- cbind(covs0[!seq_len(nobs) %in% c(m$index.treated, m$index.control, m$index.dropped), ], index=as.numeric(row.names(covs0)[!(1:nobs) %in% c(m$index.treated, m$index.control, m$index.dropped)]))
    covs.list$dropped <- cbind(covs0[m$index.dropped, ], index=m$index.dropped)
    
    treat.list$control <- treat0[m$index.control]
    treat.list$treated <- treat0[m$index.treat]
    treat.list$unmatched <- treat0[!seq_len(nobs) %in% c(m$index.treated, m$index.control, m$index.dropped)]
    treat.list$dropped <- treat0[m$index.dropped]
    
    weights.list$control <- weights.list$treated <- m$weights
    weights.list$unmatched <- rep(0, length(treat0[!seq_len(nobs) %in% c(m$index.treated, m$index.control, m$index.dropped)]))
    weights.list$dropped <- rep(0, length(m$index.dropped))
    
    data.list <- lapply(1:4, function(x) cbind(data.frame(treat=treat.list[[x]]), data.frame(weights=weights.list[[x]]), covs.list[[x]]))
    o.data <- do.call(rbind, data.list)
    o.data2 <- merge(unique(o.data[, is.na(match(names(o.data), "weights"))]), aggregate(weights~index, data=o.data, FUN=sum), by="index")
    
    #Process cluster
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && cluster %in% names(data)) {
            cluster <- data[, cluster]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
        
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.vector(val, mode = "list")) {
                val.list <- lapply(val, function(x) process.val(x, i, treat, covs, data))
                val.list <- lapply(seq_along(val.list), function(x) {
                    if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                    val.list[[x]]})
                if (length(unique(sapply(val.list, nrow))) > 1) {
                    stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                }
                
                val.df <- setNames(do.call("cbind", val.list),
                                   c(sapply(val.list, names)))
            }
            else {
                val.df <- process.val(val, i, treat, covs, data)
            }
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    treat <- o.data2$treat
    weights <- data.frame(weights = o.data2$weights)
    covs <- o.data2[, is.na(match(names(o.data2), c("treat", "weights", "index")))]
    dropped <- rep(0, length(treat))
    if (length(m$index.dropped) > 0) dropped[m$index.dropped] <- 1
    
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("treat", "cluster")
    data.frames <- c(covs.data, "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in names(lengths)[names(lengths) != "weights"]) {
            if (lengths[i] > 0 && lengths[i] != lengths["weights"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word.list(names(problematic[problematic])), " must have the same number of observations as the original call to Match()."), call. = FALSE)
    }
    
    X$treat <- treat
    X$weights <- weights
    X$discarded <- dropped
    X$distance <- NULL #NAs in distance bcause of incomplete list in Match object
    X$addl <- addl
    X$covs <- covs
    X$call <- NULL
    X$method <- "matching"
    X$cluster <- factor(cluster)

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
    X <- list(covs = NA,
              weights = NA,
              treat = NA,
              distance = NA,
              subclass = NA,
              match.strata = NA,
              addl = NA,
              method = NA,
              call = NA,
              cluster = NA,
              imp = NA,
              s.weights = NA)
    
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
    
    
    #Checks
    if (length(covs) == 0) {
        stop("covs dataframe must be specified.", call. = FALSE)
    }
    if (!is.data.frame(covs)) {
        stop("covs must be a dataframe.", call. = FALSE)
    }
    if (sum(is.na(covs)) > 0) {
        stop("Missing values exist in the covariates.", call. = FALSE)
    }
    if (length(data) > 0 && !is.data.frame(data)) {
        warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execuction will halt.")
        data <- NULL
    }
    # if (length(distance) > 0 && !is.character(distance) && !is.numeric(distance) && !is.data.frame(distance)) {
    #     stop("The argument to distance must be a vector of distance scores or the (quoted) name of a variable in data that contains distance scores.", call. = FALSE)
    # }
    
    specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
    if (length(weights) > 0) {
        if (!is.character(weights) && !is.numeric(weights) && !is.data.frame(weights) && !is.list(weights)) {
            stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
        }
        specified["weights"] <- TRUE
    }
    if (length(subclass) > 0){
        if (!is.character(subclass) && !is.numeric(subclass)) {
            stop("The argument to subclass must be a vector of subclass membership or the (quoted) name of a variable in data that contains subclass membership.", call. = FALSE)
        }
        specified["subclass"] <- TRUE
    }
    if (length(match.strata) > 0) {
        if (!is.character(match.strata) && !is.numeric(match.strata) && !is.factor(match.strata)) {
            stop("The argument to match.strata must be a vector of match stratum membership or the (quoted) name of a variable in data that contains match stratum membership.", call. = FALSE)
        }
        specified["match.strata"] <- TRUE
    }
    
    #Getting method
    if (length(method) == 0) {
        if (specified["match.strata"]) {
            if (sum(specified) > 1) {
                message(word.list(names(specified)[specified]), " are specified. Assuming \"matching\" and using match.strata and ignoring ", word.list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                weights <- subclass <- NULL
            }
            X$method <- "matching"
        }
        else if (specified["subclass"]) {
            if (sum(specified) > 1) {
                message(word.list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word.list(names(specified)[specified & names(specified)!="subclass"]), ".")
                weights <- match.strata <- NULL
            }
            X$method <- "subclassification"
            #weights <- rep(1, nrow(covs))
        }
        else if (specified["weights"]) {
            if (sum(specified) > 1) {
                message(word.list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word.list(names(specified)[specified & names(specified)!="subclass"]), ".")
                match.strata <- subclass <- NULL
            }
            else {
                message("Assuming \"weighting\". If not, specify with an argument to method.")
            }
            X$method <- "weighting"
        }
        else {
            X$method <- "matching"
        }
    }
    else if (length(method) == 1) {
        specified.method <- match.arg(method, c("weighting", "matching", "subclassification"))
        if (specified.method == "weighting") {
            if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word.list(names(specified)[specified]), " are specified. Using weights and ignoring ", word.list(names(specified)[specified & names(specified)!="weights"]), ".")
                    match.strata <- subclass <- NULL
                }
                X$method <- "weighting"
            }
            else if (specified["match.strata"]) {
                message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
                subclass <- NULL
                X$method <- "matching"
            }
            else if (specified["subclass"]) {
                message("method = \"weighting\" is specified, but no weights are present. Assuming \"subclassification\" and using subclass instead.")
                X$method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else {
                X$method <- "matching"
            }
        }
        else if (specified.method == "matching") {
            if (specified["match.strata"]) {
                if (sum(specified) > 1) {
                    message(word.list(names(specified)[specified]), " are specified. Using match.strata and ignoring ", word.list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                    weights <- subclass <- NULL
                }
                X$method <- "matching"
            }
            else if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word.list(names(specified)[specified]), " are specified. Using weights and ignoring ", word.list(names(specified)[specified & names(specified)!="weights"]), ".")
                    match.strata <- subclass <- NULL
                }
                X$method <- "matching"
            }
            else if (specified["subclass"]) {
                message("method = \"matching\" is specified, but no weights or match.strata are present. Assuming \"subclassification\" and using subclass instead.")
                X$method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else {
                X$method <- "matching"
            }
        }
        else if (specified.method == "subclassification") {
            if (specified["subclass"]) {
                if (sum(specified) > 1) {
                    message(word.list(names(specified)[specified]), " are specified. Using subclass and ignoring ", word.list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    weights <- match.strata <- NULL
                }
                X$method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else if (specified["match.strata"]) {
                message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
                weights <- NULL
                X$method <- "matching"
            }
            else if (specified["weights"]) {
                message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"weighting\" and using weights instead.")
                X$method <- "weighting"
            }
        }
    }
    else {
        specified.method <- match.arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
        if (any(specified.method == "subclassification") || specified["subclass"]) {
            stop("Subclassification cannot be specified along with other methods.", call. = FALSE)
        }
        else if (specified["match.strata"]) {
            stop("Only weights can be specified with mutiple methods.", call. = FALSE)
        }
        else if (!specified["weights"]) {
            warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
            X$method <- "matching"
        }
        else if (length(specified.method) != ncol(weights)) {
            stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
        }
        else {
            #Matching and/or weighting with various weights
            X$method <- specified.method
            match.strata <- subclass <- NULL
        }
    }
    
    if (length(cluster) > 0 && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
        stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (length(imp) > 0 && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
        stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
    }
    
    #Process treat
    if (length(treat) == 0) stop("treat must be specified.", call. = FALSE)
    
    if (is.numeric(treat) || is.factor(treat) || (is.character(treat) && length(treat) > 1)) {
        treat <- treat
    }
    else if (is.character(treat) && length(treat)==1 && treat %in% names(data)) {
        treat <- data[, treat]
    }
    else stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
    
    if (sum(is.na(treat)) > 0)
        stop("Missing values exist in treat.", call. = FALSE)
    
    if (length(unique(treat)) == 2) {
        treat <- binarize(treat)
    }
    
    #Process weights, addl, and distance
    for (i in c("weights", "addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.vector(val, mode = "list")) {
                val.list <- lapply(val, function(x) process.val(x, i, treat, covs, data))
                val.list <- lapply(seq_along(val.list), function(x) {
                    if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                    val.list[[x]]})
                if (length(unique(sapply(val.list, nrow))) > 1) {
                    stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                }
                
                val.df <- setNames(do.call("cbind", val.list),
                                   c(sapply(val.list, names)))
            }
            else {
                val.df <- process.val(val, i, treat, covs, data)
            }
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    #Process subclass
    if (length(subclass) > 0) {
        if (is.numeric(subclass) || is.factor(subclass) || (is.character(subclass) && length(subclass) > 1)) {
            subclass <- subclass
        }
        else if (is.character(subclass) && length(subclass)==1 && subclass %in% names(data)) {
            subclass <- data[, subclass]
        }
        else stop("The name supplied to subclass is not the name of a variable in data.", call. = FALSE)
        #weights <- data.frame(weights = rep(1, length(treat)))
    }
    
    #Process match.strata
    if (length(match.strata) > 0) {
        if (is.character(match.strata) && length(match.strata)==1) {
            if (match.strata %in% names(data)) {
                match.strata <- data[, match.strata]
            }
            else stop("The name supplied to match.strata is not the name of a variable in data.", call. = FALSE)
        }
    }
    #Process sampling weights
    if (length(s.weights) > 0) {
        if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
            stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
        }
        if (is.character(s.weights) && length(s.weights)==1) {
            if (s.weights %in% names(data)) {
                s.weights <- data[, s.weights]
            }
            else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
        }
    }
    else s.weights <- rep(1, length(treat))

    #Process cluster
    if (length(cluster) > 0) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && cluster %in% names(data)) {
            cluster <- data[, cluster]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("treat", "subclass", "match.strata", "cluster")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    #Process imp
    if (length(imp) > 0) {
        if (is.numeric(imp) || is.factor(imp) || (is.character(imp) && length(imp)>1)) {
            imp <- imp
        }
        else if (is.character(imp) && length(imp)==1 && imp %in% names(data)) {
            imp <- data[, imp]
        }
        else stop("The name supplied to imp is not the name of a variable in data.", call. = FALSE)
        
        imp.lengths <- sapply(unique(imp), function(i) sum(imp == i))
        
        if (abs(max(imp.lengths) - min(imp.lengths)) < sqrt(.Machine$double.eps)) { #all the same
            for (i in vectors) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) { 
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        
                        temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, temp.merge[order(temp.merge[,3]), -c(1:3), drop = TRUE])
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
            for (i in data.frames) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) {
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        temp.var <- data.frame(sort(imp),rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, setNames(temp.merge[order(temp.merge[,3]), -c(1:3), drop = FALSE], names(get(i))))
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
            stop(paste0(word.list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
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
        stop(paste0(word.list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
    }
    
    #Turn match.strata into weights
    if (length(match.strata) > 0) {
        weights <- data.frame(weights = match.strata2weights(match.strata = match.strata,
                                                             treat = treat,
                                                             covs = covs
                                                             ))
    }
    
    if (length(weights) > 0) {
        if (any(sapply(weights, function(x) !is.numeric(x)))) {
            stop("All weights must be numeric.", call. = FALSE)
        }
        if (length(X$method) == 1) {
            X$method <- rep(X$method, ncol(weights))
        }
        
    }
    

    
    #Get s.d.denom
    check.estimand <- check.weights <- bad.s.d.denom <- bad.estimand <- FALSE
    if (length(unique(treat)) <= 2 || !is.numeric(treat)) { #non-continuous
        if (length(s.d.denom) > 0) {
            try.s.d.denom <- tryCatch(match.arg(s.d.denom, c("treated", "control", "pooled"), several.ok = TRUE),
                                      error = function(cond) FALSE)
            if (any(try.s.d.denom == FALSE)) {
                check.estimand <- TRUE
                bad.s.d.denom <- TRUE
            }
            else {
                if (length(try.s.d.denom) > 1 && length(try.s.d.denom) != ncol(weights)) {
                    stop("s.d.denom must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
                }
                else X$s.d.denom <- try.s.d.denom
            }
        }
        else {
            check.estimand <- TRUE
        }
        
        if (check.estimand == TRUE) {
            if (length(estimand) > 0) {
                try.estimand <- tryCatch(match.arg(tolower(estimand), c("att", "atc", "ate"), several.ok = TRUE),
                                         error = function(cond) FALSE)
                if (any(try.estimand == FALSE)) {
                    check.weights <- TRUE
                    bad.estimand <- TRUE
                }
                else {
                    if (length(try.estimand) > 1 && length(try.estimand) != ncol(weights)) {
                        stop("estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
                    }
                    else X$s.d.denom <- sapply(try.estimand, function(x) switch(x, att = "treated", atc = "control", ate = "pooled"))
                }
            }
            else {
                check.weights <- TRUE
            }
        }
        
        if (check.weights == TRUE) {
            if (length(weights) == 0) {
                X$s.d.denom <- "treated"
                estimand <- "ATT"
            }
            else {
                X$s.d.denom <- estimand <- character(ncol(weights))
                for (i in seq_len(ncol(weights))) {
                    if (X$method[i] == "weighting") {
                        if (max(weights[treat==1 & weights[, i] > sqrt(.Machine$double.eps), i]) - min(weights[treat==1 & weights[, i] > sqrt(.Machine$double.eps), i]) < sqrt(.Machine$double.eps) &&
                            max(weights[treat==0 & weights[, i] > sqrt(.Machine$double.eps), i]) - min(weights[treat==0 & weights[, i] > sqrt(.Machine$double.eps), i]) >= sqrt(.Machine$double.eps)
                        ) { #if treated weights are only all either 0 the same; ATT
                            estimand[i] <- "att"
                            X$s.d.denom[i] <- "treated"
                        }
                        else if (max(weights[treat==0 & weights[, i] > sqrt(.Machine$double.eps), i]) - min(weights[treat==0 & weights[, i] > sqrt(.Machine$double.eps), i]) < sqrt(.Machine$double.eps) &&
                                 max(weights[treat==1 & weights[, i] > sqrt(.Machine$double.eps), i]) - min(weights[treat==1 & weights[, i] > sqrt(.Machine$double.eps), i]) >= sqrt(.Machine$double.eps)
                        ) { #if control weights are only all either 0 the same; ATC
                            estimand[i] <- "atc"
                            X$s.d.denom[i] <- "control"
                        }
                        else {
                            estimand[i] <- "ate"
                            X$s.d.denom[i] <- "pooled"
                        }
                    }
                    
                    else {
                        estimand[i] <- "att"
                        X$s.d.denom[i] <- "treated"
                    }
                }
            }
        }
        if (length(weights) > 0 && length(X$s.d.denom) == 1) X$s.d.denom <- rep(X$s.d.denom, ncol(weights))
        
        if (bad.s.d.denom && bad.estimand) {
            message("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\n         Using \"", word.list(X$s.d.denom), "\" instead.")
        }
        else if (bad.estimand) {
            message("Warning: estimand should be one of \"ATT\", \"ATC\", or \"ATE\". Using \"", ifelse(length(unique(estimand)) == 1, toupper(estimand), word.list(toupper(estimand))), "\" instead.")
        }
        else if (check.weights) {
            message("Note: estimand and s.d.denom not specified; assuming ", word.list(toupper(estimand)), " and ", word.list(X$s.d.denom), ".")
        }
        
        if (all(X$method %in% c("weighting", "matching"))) {
            if (length(weights) > 0 && length(X$s.d.denom) != ncol(weights)) {
                stop("Valid inputs to s.d.denom or estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
            }
        }
    }
    
    X$covs <- covs
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$subclass <- factor(subclass)
    X$cluster <- factor(cluster)
    X$call <- NULL
    X$addl <- addl
    X$imp <- factor(imp)
    X$s.weights <- s.weights

    return(X)
}
x2base.formula <- function(formula, ...) {
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
    
    A <- list(...)
    
    #Checks
    if (length(A$data) == 0) {
        stop("Dataframe must be specified.", call. = FALSE)}
    if (!is.data.frame(A$data)) {
        stop("Data must be a dataframe.", call. = FALSE)}
    
    #Initializing variables
    tt <- terms(formula)
    attr(tt, "intercept") <- 0
    if (is.na(match(rownames(attr(tt, "factors"))[1], names(A$data)))) {
        stop(paste0("The given response variable, \"", rownames(attr(tt, "factors"))[1], "\", is not a variable in data."))
    }
    m.try <- try({mf <- model.frame(tt, A$data)}, TRUE)
    if (class(m.try) == "try-error") {
        stop(paste0(c("All variables of formula must be variables in data.\nVariables not in data: ",
                      paste(attr(tt, "term.labels")[is.na(match(attr(tt, "term.labels"), names(A$data)))], collapse=", "))), call. = FALSE)}
    treat <- model.response(mf)
    covs <- A$data[, !is.na(match(names(A$data), attr(tt, "term.labels"))), drop = FALSE]
    X <- x2base.data.frame(covs, treat = treat, ...)
    return(X)
}
x2base.CBPS <- function(cbps.fit, ...) {
    #s.d.denom
    #cluster
    A <- list(...)
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              distance=NA,
              addl=NA,
              s.d.denom=NA,
              call=NA,
              cluster=NA)
    #Checks
    
    treat <- cbps.fit$y
    covs <- cbps.fit$data[, !is.na(match(names(cbps.fit$data), attributes(terms(cbps.fit))$term.labels))]
    weights <- data.frame(weights = get.w(cbps.fit, 
                                                use.weights = ifelse(length(A$use.weights) == 0, TRUE,
                                                                               A$use.weights)))
    
    if (!(any(class(cbps.fit) == "CBPSContinuous") || nlevels(as.factor(treat)) > 2)) {
        if (length(A$s.d.denom > 0) && is.character(A$s.d.denom)) {
            X$s.d.denom <- tryCatch(match.arg(A$s.d.denom, c("treated", "control", "pooled")),
                                    error = function(cond) {
                                        new.s.d.denom <- switch(tolower(A$estimand), att = "treated", ate = "pooled")
                                        message(paste0("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\nUsing ", deparse(new.s.d.denom), " instead."))
                                        return(new.s.d.denom)})
        }
        else {
            if (abs(max(weights[treat == 1,], na.rm = TRUE) - 
                    min(weights[treat == 1,], na.rm = TRUE)) < 
                sqrt(.Machine$double.eps)) {
                X$s.d.denom <- "treated"
            }
            else {
                X$s.d.denom <- "pooled"
            }
        }
    }
    
    c.data <- cbps.fit$data
    
    #Process cluster
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1) {
            if (cluster %in% names(data)) {
                cluster <- data[, cluster]
            }
            else if (cluster %in% names(c.data)) {
                cluster <- c.data[, cluster]
            }
        }
        else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.vector(val, mode = "list")) {
                val.list <- lapply(val, function(x) process.val(x, i, treat, covs, data))
                val.list <- lapply(seq_along(val.list), function(x) {
                    if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                    val.list[[x]]})
                if (length(unique(sapply(val.list, nrow))) > 1) {
                    stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                }
                
                val.df <- setNames(do.call("cbind", val.list),
                                   c(sapply(val.list, names)))
            }
            else {
                val.df <- process.val(val, i, treat, covs, data)
            }
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    if (abs(max(cbps.fit$fitted.values) - min(cbps.fit$fitted.values)) < sqrt(.Machine$double.eps)) {
        if (length(distance) == 0) distance <- NULL
    }
    else {
        if (length(distance) == 0) distance <- data.frame(prop.score = cbps.fit$fitted.values)
        else distance <- cbind(distance, prop.score = cbps.fit$fitted.values)
    }
    
    
    ensure.equal.lengths <- TRUE
    vectors <- c("cluster")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in names(lengths)[names(lengths) != "covs"]) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word.list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to CBPS()."), call. = FALSE)
    }
    
    X$distance <- distance
    X$addl <- addl
    X$weights <- weights
    X$treat <- treat
    X$covs <- covs
    X$cluster <- factor(cluster)
    X$call <- cbps.fit$call

    return(X)
}
x2base.ebalance <- function(ebalance, ...) {
    #formula
    #data
    #treat
    #covs
    A <- list(...)
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              method=NA,
              distance=NA,
              call=NA,
              cluster = NA,
              addl=NA)
    
    #Get treat and covs
    data <- A$data
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    
    if (sum(is.na(t.c$covs))>0)
        stop("Missing values exist in the covariates.", call. = FALSE)
    
    #Initializing variables
    
    treat <- t.c$treat
    covs  <- t.c$covs
    # weights <- rep(1, length(treat))
    # 
    # if (length(ebalance$w) != sum(treat == 0)) {
    #     stop("There are more control units in treat than weights in the ebalance object.", call. = FALSE)
    # }
    # weights[treat == 0] <- ebalance$w
    weights <- data.frame(weights = get.w(ebalance, treat))
    
    #Process cluster
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && cluster %in% names(data)) {
            cluster <- data[, cluster]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
        
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.vector(val, mode = "list")) {
                val.list <- lapply(val, function(x) process.val(x, i, treat, covs, data))
                val.list <- lapply(seq_along(val.list), function(x) {
                    if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                    val.list[[x]]})
                if (length(unique(sapply(val.list, nrow))) > 1) {
                    stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                }
                
                val.df <- setNames(do.call("cbind", val.list),
                                   c(sapply(val.list, names)))
            }
            else {
                val.df <- process.val(val, i, treat, covs, data)
            }
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("treat", "cluster")
    data.frames <- c(covs.data, "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in names(lengths)[names(lengths) != "weights"]) {
            if (lengths[i] > 0 && lengths[i] != lengths["weights"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word.list(names(problematic[problematic])), " must have the same number of observations as the original call to ebalance()."), call. = FALSE)
    }
    
    X$treat <- treat
    X$weights <- weights
    X$covs <- covs
    X$distance <- distance
    X$addl <- addl
    X$call <- NULL
    X$method <- "weighting"
    X$cluster <- factor(cluster)

    return(X)
}
x2base.optmatch <- function(optmatch, ...) {
    #formula
    #data
    #treat
    #covs
    A <- list(...)
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              method=NA,
              distance=NA,
              call=NA,
              cluster = NA)
    
    #Get treat and covs
    data <- A$data
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    
    if (sum(is.na(t.c$covs))>0)
        stop("Missing values exist in the covariates.", call. = FALSE)
    
    #Initializing variables
    treat <- t.c$treat
    covs  <- t.c$covs
    distance <- A$distance
    
    #Process match.strata (optmatch)
    if (length(optmatch) != length(treat) || length(optmatch) != nrow(covs)) {
        stop(paste0("The optmatch object must have the same length as ", ifelse(attr(t.c, "which")=="fd", "data", "covs"), "."), call. = FALSE)
    }
    weights <- data.frame(weights = get.w(optmatch, treat = treat))
    
    #Process cluster
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && cluster %in% names(data)) {
            cluster <- data[, cluster]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.vector(val, mode = "list")) {
                val.list <- lapply(val, function(x) process.val(x, i, treat, covs, data))
                val.list <- lapply(seq_along(val.list), function(x) {
                    if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                    val.list[[x]]})
                if (length(unique(sapply(val.list, nrow))) > 1) {
                    stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                }
                
                val.df <- setNames(do.call("cbind", val.list),
                                   c(sapply(val.list, names)))
            }
            else {
                val.df <- process.val(val, i, treat, covs, data)
            }
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("treat", "cluster")
    data.frames <- c(covs.data, "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in names(lengths)[names(lengths) != "weights"]) {
            if (lengths[i] > 0 && lengths[i] != lengths["weights"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word.list(names(problematic[problematic])), " must have the same number of observations as the original call to optmatch()."), call. = FALSE)
    }
    
    X$treat <- treat
    X$distance <- distance
    X$covs <- covs
    X$weights <- weights
    X$addl <- addl
    X$call <- NULL
    X$method <- "matching"
    X$cluster <- factor(cluster)

    return(X)
    
}
