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
              obj=NA,
              call=NA,
              cluster=NA)
    
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
    X$weights <- m$weights
    X$treat <- m$treat
    
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
            if (is.numeric(val)) {
                val.df <- setNames(data.frame(val), i)
            }
            else if (is.character(val)) {
                if (length(data) > 0 && any(val %in% names(data))) {
                    val.df <- data[, val[val %in% names(data)], drop = FALSE]
                    val <- val[!val %in% names(data)]
                    
                }
                if (length(m.data) > 0 && any(val %in% names(m.data))) {
                    if (length(val.df) > 0) val.df <- cbind(val.df, m.data[, val[val %in% names(m.data)], drop = FALSE])
                    else val.df <- m.data[, val[val %in% names(m.data)], drop = FALSE]
                    val <- val[!val %in% names(m.data)]
                    
                }
                if (length(val) > 0) {
                    
                    warning(paste("The following variable(s) named in", i, "are not in any given data set and will be ignored: ",
                                  paste(val)))
                }
            }
            else if (is.data.frame(val)) {
                val.df <- val
            }
            else stop(paste("The names supplied to", i, "are not the name of a variable in any given data set."), call. = FALSE)
            
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
    vectors <- c("cluster")
    data.frames <- c("covs", "distance", "addl")
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
    
    X$covs <- covs
    X$distance <- distance
    X$addl <- addl
    X$cluster <- factor(cluster)
    X$obj <- m
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
    if (length(X$subclass) > 0) X$obj$subclass <- X$subclass
    X$call <- m$call
    return(X)
}
x2base.ps <- function(ps, ...) {
    #full.stop.method
    #s.d.denom
    A <- list(...)
    
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              distance=NA,
              addl=NA,
              s.d.denom=NA,
              call=NA,
              obj=NA,
              cluster = NA)
    
    #Initializing variables
    
    if (length(A$full.stop.method) > 1) {
        warning("The argument to full.stop.method must have length 1; using the first available stop method instead.", call. = FALSE)
        rule1 <- names(ps$w)[1]
    }
    else if (length(A$full.stop.method)==1) {
        if (is.character(A$full.stop.method)) {
            rule1 <- tryCatch(match.arg(tolower(A$full.stop.method), tolower(names(ps$w))),
                              error = function(cond) {message(paste0("Warning: full.stop.method should be one of ", paste(deparse(names(ps$w)), collapse = ", "), ".\nUsing ", deparse(names(ps$w)[1]), " instead."));
                                  return(names(ps$w)[1])})
        }
        else if (is.numeric(A$full.stop.method) && A$full.stop.method %in% seq_along(names(ps$w))) {
            rule1 <- names(ps$w)[A$full.stop.method]
        }
        else {
            warning("full.stop.method should be ", word.list(sapply(names(ps$w), deparse), "or"), ".\nUsing ", deparse(names(ps$w)[1]), " instead.", call. = FALSE)
            rule1 <- names(ps$w)[1]
        }
    }
    else {
        rule1 <- names(ps$w)[1]
    }
    
    s <- names(ps$w)[match(tolower(rule1), tolower(names(ps$w)))]
    
    attr(ps, "which") <- s
    
    if (length(A$s.d.denom>0) && is.character(A$s.d.denom)) {
        X$s.d.denom <- tryCatch(match.arg(A$s.d.denom, c("treated", "control", "pooled")),
                                error = function(cond) {
                                    new.s.d.denom <- switch(substr(tolower(s), nchar(s)-2, nchar(s)), att = "treated", ate = "pooled")
                                    message(paste0("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\nUsing ", deparse(new.s.d.denom), " instead."))
                                    return(new.s.d.denom)})
    }
    else X$s.d.denom <- switch(substr(tolower(s), nchar(s)-2, nchar(s)), att = "treated", ate = "pooled")
    
    weights <- ps$w[, s]
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
            if (is.numeric(val)) {
                val.df <- setNames(data.frame(val), i)
            }
            else if (is.character(val)) {
                if (length(data) > 0 && any(val %in% names(data))) {
                    val.df <- data[, val[val %in% names(data)], drop = FALSE]
                    val <- val[!val %in% names(data)]
                    
                }
                if (length(ps.data) > 0 && any(val %in% names(ps.data))) {
                    if (length(val.df) > 0) val.df <- cbind(val.df, ps.data[, val[val %in% names(ps.data)], drop = FALSE])
                    else val.df <- ps.data[, val[val %in% names(ps.data)], drop = FALSE]
                    val <- val[!val %in% names(ps.data)]
                    
                }
                if (length(val) > 0) {
                    
                    warning(paste("The following variable(s) named in", i, "are not in any given data set and will be ignored: ",
                                  paste(val)))
                }
            }
            else if (is.data.frame(val)) {
                val.df <- val
            }
            else stop(paste("The names supplied to", i, "are not the name of a variable in any given data set."), call. = FALSE)
            
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
    data.frames <- c("covs", "distance", "addl")
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
    X$obj <- list(treat=treat, weights=weights)
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
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
              obj=NA,
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
            if (is.numeric(val)) {
                val.df <- setNames(data.frame(val), i)
            }
            else if (is.character(val)) {
                if (length(A$data) > 0 && any(val %in% names(data))) {
                    val.df <- data[, val[val %in% names(data)], drop = FALSE]
                    val <- val[!val %in% names(data)]
                }
                if (length(val) > 0) {
                    warning(paste("The following variable(s) named in", i, "are not in data and will be ignored: ",
                                  paste(val)))
                }
            }
            else if (is.data.frame(val)) {
                val.df <- val
            }
            else stop(paste("The names supplied to", i, "are not the name of a variable in data."), call. = FALSE)
            
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    treat <- o.data2$treat
    weights <- o.data2$weights
    covs <- o.data2[, is.na(match(names(o.data2), c("treat", "weights", "index")))]
    
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("weights", "treat", "cluster")
    data.frames <- c(covs.data, "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors[vectors!="weights"], data.frames)) {
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
    X$distance <- NULL #NAs in distance bcause of incomplete list in Match object
    X$addl <- addl
    X$covs <- covs
    X$call <- NULL
    X$method <- "matching"
    X$cluster <- factor(cluster)
    X$obj <- list(treat=treat, weights=weights)
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
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
              obj = NA,
              cluster = NA,
              imp = NA)
    
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
    if (length(distance) > 0 && !is.character(distance) && !is.numeric(distance) && !is.data.frame(distance)) {
        stop("The argument to distance must be a vector of distance scores or the (quoted) name of a variable in data that contains distance scores.", call. = FALSE)
    }
    
    specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
    if (length(weights) > 0) {
        if (!is.character(weights) && !is.numeric(weights)) {
            stop("The argument to weights must be a vector of weights or the (quoted) name of a variable in data that contains weights.", call. = FALSE)
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
            weights <- rep(1, nrow(covs))
        }
        else if (specified["weights"]) {
            if (sum(specified) > 1) {
                message(word.list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word.list(names(specified)[specified & names(specified)!="subclass"]), ".")
                match.strata <- match.strata <- NULL
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
    else {
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
                weights <- rep(1, nrow(covs))
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
                weights <- rep(1, nrow(covs))
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
                weights <- rep(1, nrow(covs))
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
    
    if (length(cluster) > 0 && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
        stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (length(imp) > 0 && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
        stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
    }
    if (length(treat) == 0) stop("treat must be specified.", call. = FALSE)
    
    #Process treat
    if (is.numeric(treat) || is.factor(treat) || (is.character(treat) && length(treat) > 1)) {
        treat <- treat
    }
    else if (is.character(treat) && length(treat)==1 && treat %in% names(data)) {
        treat <- data[, treat]
    }
    else stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
    
    if (sum(is.na(treat)) > 0)
        stop("Missing values exist in treat.", call. = FALSE)
    
    #Process weights
    if (length(weights) > 0) {
        if (is.numeric(weights)) {
            weights <- weights
        }
        else if (is.character(weights) && length(weights)==1 && weights %in% names(data)) {
            weights <- data[, weights]
        }
        else stop("The name supplied to weights is not the name of a variable in data.", call. = FALSE)
        
        if (sum(is.na(weights)) > 0)
            stop("Missing values exist in weights.", call. = FALSE)
    }
    
    #Process subclass
    if (length(subclass) > 0) {
        if (is.numeric(subclass) || is.factor(subclass) || (is.character(subclass) && length(subclass)>1)) {
            subclass <- subclass
        }
        else if (is.character(subclass) && length(subclass)==1 && subclass %in% names(data)) {
            subclass <- data[, subclass]
        }
        else stop("The name supplied to subclass is not the name of a variable in data.", call. = FALSE)
        
    }
    
    #Process match.strata
    if (length(match.strata) > 0) {
        if (is.character(match.strata) && length(match.strata)==1) {
            if (match.strata %in% names(data)) {
                match.strata <- data[, match.strata]
            }
            else stop("The name supplied to match.strata is not the name of a variable in data.", call. = FALSE)
        }
        
        weights <- match.strata2weights(covs = covs, 
                                        treat = treat, 
                                        match.strata = match.strata)
    }
    
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
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.numeric(val)) {
                val.df <- setNames(data.frame(val), i)
            }
            else if (is.character(val)) {
                if (length(data) > 0 && any(val %in% names(data))) {
                    val.df <- data[, val[val %in% names(data)], drop = FALSE]
                    val <- val[!val %in% names(data)]
                }
                if (length(val) > 0) {
                    warning(paste("The following variable(s) named in", i, "are not in data and will be ignored: ",
                                  paste(val)))
                }
            }
            else if (is.data.frame(val)) {
                val.df <- val
            }
            else stop(paste("The names supplied to", i, "are not the name of a variable in data."), call. = FALSE)
            
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("treat", "weights", "subclass", "match.strata", "cluster")
    data.frames <- c("covs", "distance", "addl")
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
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), ]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, setNames(temp.merge[order(temp.merge[,3]), -c(1:3)], names(get(i))))
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
    
    #Get s.d.denom
    check.estimand <- check.weights <- bad.s.d.denom <- bad.estimand <- FALSE
    if (length(unique(treat)) <= 2 || !is.numeric(treat)) { #non-continuous
        if (length(s.d.denom) > 0) {
            try.s.d.denom <- tryCatch(match.arg(s.d.denom, c("treated", "control", "pooled")),
                                      error = function(cond) FALSE)
            if (try.s.d.denom == FALSE) {
                check.estimand <- TRUE
                bad.s.d.denom <- TRUE
            }
            else {
                X$s.d.denom <- try.s.d.denom
            }
        }
        else {
            check.estimand <- TRUE
        }
        
        if (check.estimand == TRUE) {
            if (length(estimand) > 0) {
                try.estimand <- tryCatch(match.arg(tolower(estimand), c("att", "atc", "ate")),
                                         error = function(cond) FALSE)
                if (try.estimand == FALSE) {
                    check.weights <- TRUE
                    bad.estimand <- TRUE
                }
                else {
                    X$s.d.denom <- switch(try.estimand, att = "treated", atc = "control", ate = "pooled")
                }
            }
            else {
                check.weights <- TRUE
            }
        }
        
        if (check.weights == TRUE) {
            if (X$method == "weighting") {
                if (max(weights[treat==1 & weights > sqrt(.Machine$double.eps)]) - min(weights[treat==1 & weights > sqrt(.Machine$double.eps)]) < sqrt(.Machine$double.eps) &&
                    max(weights[treat==0 & weights > sqrt(.Machine$double.eps)]) - min(weights[treat==0 & weights > sqrt(.Machine$double.eps)]) >= sqrt(.Machine$double.eps)
                ) { #if treated weights are only all either 0 the same; ATT
                    estimand <- "att"
                    X$s.d.denom <- "treated"
                }
                else if (max(weights[treat==0 & weights > sqrt(.Machine$double.eps)]) - min(weights[treat==0 & weights > sqrt(.Machine$double.eps)]) < sqrt(.Machine$double.eps) &&
                         max(weights[treat==1 & weights > sqrt(.Machine$double.eps)]) - min(weights[treat==1 & weights > sqrt(.Machine$double.eps)]) >= sqrt(.Machine$double.eps)
                ) { #if control weights are only all either 0 the same; ATC
                    estimand <- "atc"
                    X$s.d.denom <- "control"
                }
                else {
                    estimand <- "ate"
                    X$s.d.denom <- "pooled"
                }
            }
            else {
                X$s.d.denom <- "treated"
            }
        }
        
        if (bad.s.d.denom && bad.estimand) {
            message("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\n         Using ", deparse(X$s.d.denom), " instead.")
        }
        else if (bad.estimand) {
            message("Warning: estimand should be one of \"ATT\", \"ATC\", or \"ATE\". Using ", deparse(toupper(estimand)), " instead.")
        }
        else if (check.weights && X$method == "weighting") {
            message("Note: estimand and s.d.denom not specified; assuming ", deparse(toupper(estimand)), " and ", deparse(X$s.d.denom), ".")
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
    X$obj <- data.frame(treat=X$treat, weights=NA)
    X$imp <- factor(imp)
    if (length(weights) > 0) X$obj$weights <- X$weights
    if (length(subclass) > 0) X$obj$subclass <- X$subclass
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
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
              obj=NA,
              cluster=NA)
    #Checks
    
    treat <- cbps.fit$y
    covs <- cbps.fit$data[, !is.na(match(names(cbps.fit$data), attributes(terms(cbps.fit))$term.labels))]
    weights <- cbps.fit$weights
    
    if (!(any(class(cbps.fit) == "CBPSContinuous") || nlevels(as.factor(treat)) > 2)) {
        # if (!std.ok && sum(weights) < 3) {
        #     if ((length(A$estimand > 0) && is.character(A$estimand)) || (length(A$s.d.denom > 0) && is.character(A$s.d.denom))) {
        #         #warning("Standardized weights were used; this may cause reported values to be incorrect. Use unstandardized weights instead.", call. = FALSE)
        #     }
        #     else {
        #         stop("Please specify either the estimand (\"ATT\" or \"ATE\") or an argument to s.d.denom.", call. = FALSE)
        #     }
        # }
        # else {
        #     if (isTRUE(all.equal(weights, treat / cbps.fit$fitted.values + (1-treat) / (1-cbps.fit$fitted.values)))) A$estimand <- "ATE"
        #     else A$estimand <- "ATT"
        # }
        if (length(A$s.d.denom > 0) && is.character(A$s.d.denom)) {
            X$s.d.denom <- tryCatch(match.arg(A$s.d.denom, c("treated", "control", "pooled")),
                                    error = function(cond) {
                                        new.s.d.denom <- switch(tolower(A$estimand), att = "treated", ate = "pooled")
                                        message(paste0("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\nUsing ", deparse(new.s.d.denom), " instead."))
                                        return(new.s.d.denom)})
        }
        else {
            if (abs(max(weights[treat == 1], na.rm = TRUE) - 
                    min(weights[treat == 1], na.rm = TRUE)) < 
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
            if (is.numeric(val)) {
                val.df <- setNames(data.frame(val), i)
            }
            else if (is.character(val)) {
                if (length(data) > 0 && any(val %in% names(data))) {
                    val.df <- data[, val[val %in% names(data)], drop = FALSE]
                    val <- val[!val %in% names(data)]
                    
                }
                if (length(c.data) > 0 && any(val %in% names(c.data))) {
                    if (length(val.df) > 0) val.df <- cbind(val.df, c.data[, val[val %in% names(c.data)], drop = FALSE])
                    else val.df <- c.data[, val[val %in% names(c.data)], drop = FALSE]
                    val <- val[!val %in% names(c.data)]
                    
                }
                if (length(val) > 0) {
                    
                    warning(paste("The following variable(s) named in", i, "are not in any given data set and will be ignored: ",
                                  paste(val)))
                }
            }
            else if (is.data.frame(val)) {
                val.df <- val
            }
            else stop(paste("The names supplied to", i, "are not the name of a variable in any given data set."), call. = FALSE)
            
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
    data.frames <- c("covs", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors[vectors!="weights"], data.frames)) {
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
    X$obj <- list(treat = treat, weights = weights)
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
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
              obj=NA,
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
    weights <- rep(1, length(treat))
    
    if (length(ebalance$w) != sum(treat == 0)) {
        stop("There are more control units in treat than weights in the ebalance object.", call. = FALSE)
    }
    weights[treat == 0] <- ebalance$w
    
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
            if (is.numeric(val)) {
                val.df <- setNames(data.frame(val), i)
            }
            else if (is.character(val)) {
                if (length(data) > 0 && any(val %in% names(data))) {
                    val.df <- data[, val[val %in% names(data)], drop = FALSE]
                    val <- val[!val %in% names(data)]
                }
                if (length(val) > 0) {
                    warning(paste("The following variable(s) named in", i, "are not in data and will be ignored: ",
                                  paste(val)))
                }
            }
            else if (is.data.frame(val)) {
                val.df <- val
            }
            else stop(paste("The names supplied to", i, "are not the name of a variable in data."), call. = FALSE)
            
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("weights", "treat", "cluster")
    data.frames <- c(covs.data, "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors[vectors!="weights"], data.frames)) {
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
    X$obj <- list(treat=treat, weights=weights)
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
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
              obj=NA,
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
    weights <- match.strata2weights(covs = covs, 
                                    treat = treat, 
                                    match.strata = optmatch)
    
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
            if (is.numeric(val)) {
                val.df <- setNames(data.frame(val), i)
            }
            else if (is.character(val)) {
                if (length(data) > 0 && any(val %in% names(data))) {
                    val.df <- data[, val[val %in% names(data)], drop = FALSE]
                    val <- val[!val %in% names(data)]
                }
                if (length(val) > 0) {
                    warning(paste("The following variable(s) named in", i, "are not in data and will be ignored: ",
                                  paste(val)))
                }
            }
            else if (is.data.frame(val)) {
                val.df <- val
            }
            else stop(paste("The names supplied to", i, "are not the name of a variable in data."), call. = FALSE)
            
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    ensure.equal.lengths <- TRUE
    covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
    vectors <- c("weights", "treat", "cluster")
    data.frames <- c(covs.data, "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors[vectors!="weights"], data.frames)) {
            print(i)
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
    X$obj <- list(treat=X$treat, weights=X$weights)
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
    return(X)
    
}