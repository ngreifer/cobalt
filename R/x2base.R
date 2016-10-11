#Functions to convert object to base.bal.tab input

x2base <- function(obj, ...) UseMethod("x2base")

x2base.matchit <- function(m, ...) {
    A <- list(...)
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              subclass=NA,
              method=NA,
              distance=NA,
              obj=NA,
              call=NA,
              cluster=NA)
    
    #Initializing variables
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (!(is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1))) {
            stop("The argument to cluster must be a vector of cluster membership.", call. = FALSE)
        }
        if (length(cluster) != length(m$treat)) {
            stop("cluster must be the same length as the original data set.", call. = FALSE)
        }
    }
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
    if (!all(is.na(m$distance))) X$distance <- m$distance
    else X$distance <- NULL
    if (length(m$model$model) > 0) {
        o.data <- m$model$model #data used in the PS formula, including treatment and covs
        X$covs <- data.frame(o.data[, !is.na(match(names(o.data), attributes(terms(m$model))$term.labels))])
    }
    else {
        X$covs <- data.frame(m$X)
    }

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
              s.d.denom=NA,
              call=NA,
              obj=NA,
              cluster = NA)

    #Initializing variables
    
    if (length(A$full.stop.method) > 1) {
        warning("The argument to fll.stop.method must have length 1; using the first available stop method instead.", call. = FALSE)
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
            warning("full.stop.method should be one of ", paste(deparse(names(ps$w)), collapse = ", "), ".\nUsing ", deparse(names(ps$w)[1]), " instead.", call. = FALSE)
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
    
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (!(is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1))) {
            stop("The argument to cluster must be a vector of cluster membership.", call. = FALSE)
        }
        if (length(cluster) != length(ps$treat)) {
            stop("cluster must be the same length as the original data set.", call. = FALSE)
        }
    }
    X$weights <- as.matrix(ps$w[s])
    X$treat <- ps$treat
    X$distance <- ps$ps[s][, ]
    X$covs <- ps$data[ps$gbm.obj$var.names]
    X$call <- ps$parameters
    X$cluster <- factor(cluster)
    X$obj <- list(treat=X$treat, weights=X$weights)
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
    return(X)
}
x2base.Match <- function(Match, ...) {
    #formula
    #data
    #treat
    #covs
    #s.d.denom
    A <- list(...)
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              method=NA,
              distance=NA,
              call=NA,
              s.d.denom=NA,
              obj=NA,
              cluster = NA)
    #Checks
    if (!is.list(Match) & !is.null(Match)) {
        stop("'Match' object contains no valid matches")}
    useWhich <- function(formula, data, treat, covs) {
        #output: tc, fd, both
        good <- data.frame(formula=0, data=0, covs=0, treat=0)
        if (!is.null(formula) & class(formula)=="formula") good[1] <- 1
        if (!is.null(data) & is.data.frame(data)) good[2] <- 1
        if (!is.null(covs) & is.data.frame(covs)) good[3] <- 1
        if (!is.null(treat) & length(unique(treat))==2) good[4] <- 1

        if (any(c(0,1) == sum(good))) {
            stop("Either formula and data or treat and covs must be specified correctly.", call. = FALSE)}
        else if (sum(good)==2) {
            if (sum(good[1:2])==0) {
                use.which <- "tc"}
            else if (sum(good[1:2])==1) {
                stop("Either formula and data or treat and covs must be specified correctly.", call. = FALSE)}
            else if (sum(good[1:2])==2) {
                use.which <- "fd"}
        }
        else if (sum(good)==3) {
            bad <- names(good)[match(0, good)]
            if (match(0, good) %in% 1:2) {
                warning(paste("Argument to", bad, "is missing; using treat and covs instead."), call. = FALSE, immediate.=TRUE)
                use.which <- "tc"}
            else {
                warning(paste("Argument to", bad, "is missing; using formula and data instead."), call. = FALSE, immediate.=TRUE)
                use.which <- "fd"}
        }
        else if (sum(good)==4) {
            use.which <- "both"
        }
        return(use.which)
    }
    use.which <- useWhich(A$formula, A$data, A$treat, A$covs)

    use.fd <- function(formula, data){
        #outputs a list containing treat [1] and covs [2]
        out.list <- list(treat=NA, covs=data.frame(NA))
        tt <- terms(formula)
        attr(tt, "intercept") <- 0
        mf<- tryCatch(model.frame(tt, data),
                          error = function(cond) stop(paste0(c("All right hand side variables of formula must be variables in data.\nVariables not in data: ",
                                                               paste(attr(tt, "term.labels")[which(!attr(tt, "term.labels") %in% names(data))], collapse=", "))), call. = FALSE))
        out.list$treat <- model.response(mf) #treat
        out.list$covs <- as.data.frame(model.matrix(tt, data=mf)) #covs
        return(out.list)
    }
    use.tc <- function(treat, covs) {
        if (length(treat)!=nrow(covs)) {
            stop("treat must be the same length as covs", call. = FALSE)}
        out.list <- list(treat=treat, covs=covs)
        return(out.list)
    }

    if (use.which == "fd") {
        t.c <- use.fd(A$formula, A$data)}
    else if (use.which == "tc") {
        t.c <- use.tc(A$treat, A$covs)}
    else if (use.which == "both") {
        try.fd <- try({t.c <- use.fd(A$formula, A$data)})
        if (class(try.fd) == "try-error") {
            message("Formula, data, treat, and covs all supplied; ignoring formula and data.")
            t.c <- use.tc(A$treat, A$covs)}
        else {
            message("Formula, data, treat, and covs all supplied; ignoring treat and covs.")}
    }

    if (sum(is.na(t.c$covs))>0)
        stop("Missing values exist in the covariates", call. = FALSE)

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
    covs.list$unmatched <- cbind(covs0[!(1:nobs) %in% c(m$index.treated, m$index.control, m$index.dropped), ], index=as.numeric(row.names(covs0)[!(1:nobs) %in% c(m$index.treated, m$index.control, m$index.dropped)]))
    covs.list$dropped <- cbind(covs0[m$index.dropped, ], index=m$index.dropped)

    treat.list$control <- treat0[m$index.control]
    treat.list$treated <- treat0[m$index.treat]
    treat.list$unmatched <- treat0[!(1:nobs) %in% c(m$index.treated, m$index.control, m$index.dropped)]
    treat.list$dropped <- treat0[m$index.dropped]

    weights.list$control <- weights.list$treated <- m$weights
    weights.list$unmatched <- rep(0, length(treat0[!(1:nobs) %in% c(m$index.treated, m$index.control, m$index.dropped)]))
    weights.list$dropped <- rep(0, length(m$index.dropped))

    data.list <- lapply(1:4, function(x) cbind(data.frame(treat=treat.list[[x]]), data.frame(weights=weights.list[[x]]), covs.list[[x]]))
    o.data <- do.call(rbind, data.list)
    o.data2 <- merge(unique(o.data[, is.na(match(names(o.data), "weights"))]), aggregate(weights~index, data=o.data, FUN=sum), by="index")
    
    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (!(is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1))) {
            stop("The argument to cluster must be a vector of cluster membership.", call. = FALSE)
        }
        if (length(cluster) != length(m$treat)) {
            stop("cluster must be the same length as the original data set.", call. = FALSE)
        }
    }
    
    X$treat <- o.data2$treat
    X$weights <- o.data2$weights
    X$distance <- NULL #NAs in distance bcause of incomplete list in Match object
    X$covs <- o.data2[, is.na(match(names(o.data2), c("treat", "weights", "index")))]
    X$call <- NULL
    X$method <- "matching"
    X$cluster <- factor(cluster)
    X$obj <- list(treat=X$treat, weights=X$weights)
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
    return(X)
}
x2base.formula <- function(formula, ...) {
    #data
    #weights
    #distance
    #subclass
    #addl
    #method
    #cluster
    
    A <- list(...)
    
    X <- list(covs = NA,
              weights = NA,
              treat = NA,
              distance = NA,
              subclass = NA,
              addl = NA,
              method = NA,
              call = NA,
              obj = NA,
              cluster = NA)
    #Checks
    if (is.null(A$data)) {
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
    #addl
    #method
    #cluster
    
    A <- list(...)
    X <- list(covs = NA,
              weights = NA,
              treat = NA,
              distance = NA,
              subclass = NA,
              addl = NA,
              method = NA,
              call = NA,
              obj = NA,
              cluster = NA)

    treat <- A$treat
    data <- A$data
    weights <- A$weights
    distance <- A$distance
    subclass <- A$subclass
    cluster <- A$cluster
    addl <- A$addl
    method <- A$method
    
    #Checks
    if (is.null(covs)) {
        stop("covs dataframe must be specified", call. = FALSE)
    }
    if (!is.data.frame(covs)) {
        stop("covs must be a dataframe", call. = FALSE)
    }
    if (sum(is.na(covs)) > 0) {
        stop("Missing values exist in the covariates", call. = FALSE)
    }
    if (!is.null(data) && !is.data.frame(data)) {
        warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execuction will halt.")
        data <- NULL
    }
    if (!is.null(weights) && !is.character(weights) && !is.numeric(weights)) {
        stop("The argument to weights must be a vector of weights or the (quoted) name of a variable in data that contains weights.", call. = FALSE)
    }
    if (!is.null(distance) && !is.character(distance) && !is.numeric(distance)) {
        stop("The argument to distance must be a vector of distance scores or the (quoted) name of a variable in data that contains distance scores.", call. = FALSE)
    }
    if (!is.null(subclass) && !is.character(subclass) && !is.numeric(subclass)) {
        stop("The argument to subclass must be a vector of subclass membership or the (quoted) name of a variable in data that contains subclass membership.", call. = FALSE)
    }
    if (!is.null(cluster) && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
        stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (is.null(treat)) stop("treat must be specified.", call. = FALSE)
    # else if (!is.character(treat) && !is.numeric(treat)) {
    #     stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status", call. = FALSE)
    # }

    if (is.numeric(treat) || is.factor(treat) || (is.character(treat) && length(treat) > 1)) {
        treat <- treat
    }
    else if (is.character(treat) && length(treat)==1 && treat %in% names(data)) {
        treat <- data[, treat]
    }
    else stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)

    if (length(treat) != nrow(covs)) {
        stop("treat must be the same length as covs", call. = FALSE)}

    if (sum(is.na(treat)) > 0)
        stop("Missing values exist in treat", call. = FALSE)

    if (!is.null(weights)) {
        if (is.numeric(weights)) {
            weights <- weights
        }
        else if (is.character(weights) && length(weights)==1 && weights %in% names(data)) {
            weights <- data[, weights]
        }
        else stop("The name supplied to weights is not the name of a variable in data.", call. = FALSE)

        if (length(weights) != nrow(covs)) {
            stop("weights must be the same length as covs", call. = FALSE)
        }

        if (sum(is.na(weights)) > 0)
            stop("Missing values exist in weights", call. = FALSE)
    }

    if (!is.null(distance)) {
        if (is.numeric(distance)) {
            distance <- distance
        }
        else if (is.character(distance) && length(distance) == 1 && distance %in% names(data)) {
            distance <- data[, distance]
        }
        else stop("The name supplied to distance is not the name of a variable in data.", call. = FALSE)

        if (length(distance) != nrow(covs)) {
            stop("distance must be the same length as covs", call. = FALSE)
        }

        if (sum(is.na(distance)) > 0)
            stop("Missing values exist in distance", call. = FALSE)
    }

    if (!is.null(subclass)) {
        if (is.numeric(subclass)) {
            subclass <- subclass
        }
        else if (is.character(subclass) && length(subclass)==1 && subclass %in% names(data)) {
            subclass <- data[, subclass]
        }
        else stop("The name supplied to subclass is not the name of a variable in data.", call. = FALSE)

        if (length(subclass) != nrow(covs)) {
            stop("subclass must be the same length as covs", call. = FALSE)
        }
    }
    
    if (!is.null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && cluster %in% names(data)) {
            cluster <- data[, cluster]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
        
        if (length(cluster) != nrow(covs)) {
            stop("cluster must be the same length as covs", call. = FALSE)
        }
    }

    if (!is.null(addl)) {
        if (is.character(addl)) {
            if (any(!addl %in% names(data))) {
                warning(paste("The following variable(s) named in addl are not in data and will be ignored: ",
                              paste(addl[which(!addl %in% names(data))], collapse=", ")))
                addl <- data[, addl[which(addl %in% names(data))]]
            }
        }
        else if (is.data.frame(addl)) {
            if (nrow(addl)!=nrow(covs)) {
                stop("If addl is a data.frame, it must have the same number of rows as covs")
            }
        }
        else {
            warning("addl must be a list of names of variables in data or a data.frame containing additional variable(s). addl will be ignored in the following output.")
            addl <- NULL
        }
    }

    if (!is.null(weights)) {
        if (!is.null(method)) {
            X$method <- match.arg(method, c("weighting", "matching", "subclassification"))
        }
        else {
            message("Assuming weights generated through weighting; if not, please specify with argument to method.")
            X$method <- "weighting"
        }
    }
    else if (!is.null(subclass)){
        X$method <- "subclassification"
        weights <- rep(1, length(treat))
    }
    else X$method <- "matching"
    
    X$covs <- covs
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$subclass <- factor(subclass)
    X$cluster <- factor(cluster)
    X$call <- NULL
    X$addl <- addl
    X$obj <- data.frame(treat=X$treat, weights=NA)
    if (!is.null(weights)) X$obj$weights <- X$weights
    if (!is.null(subclass)) X$obj$subclass <- X$subclass
    if (!is.null(cluster)) X$obj$cluster <- X$cluster
    return(X)
}
x2base.CBPS <- function(cbps.fit, std.ok = FALSE, ...) {
    #estimand
    #s.d.denom
    #cluster
    A <- list(...)
    X <- list(covs=NA,
              treat=NA,
              weights=NA,
              distance=NA,
              s.d.denom=NA,
              call=NA,
              obj=NA,
              cluster=NA)
    #Checks

    if (!(any(class(cbps.fit) == "CBPSContinuous") || nlevels(as.factor(cbps.fit$y)) > 2)) {
        if (!std.ok && sum(cbps.fit$weights) < 3) {
            if ((length(A$estimand > 0) && is.character(A$estimand)) || (length(A$s.d.denom > 0) && is.character(A$s.d.denom))) warning("Standardized weights were used; this may cause reported values to be incorrect. Use unstandardized weights instead.", call. = FALSE)
            else stop("Please specify either the estimand (\"ATT\" or \"ATE\") or an argument to s.d.denom.", call. = FALSE)
        }
        else {
            if (isTRUE(all.equal(cbps.fit$weights, cbps.fit$y / cbps.fit$fitted.values + (1-cbps.fit$y) / (1-cbps.fit$fitted.values)))) A$estimand <- "ATE"
            else A$estimand <- "ATT"
        }
        if (length(A$s.d.denom > 0) && is.character(A$s.d.denom)) {
            X$s.d.denom <- tryCatch(match.arg(A$s.d.denom, c("treated", "control", "pooled")),
                                    error = function(cond) {
                                        new.s.d.denom <- switch(tolower(A$estimand), att = "treated", ate = "pooled")
                                        message(paste0("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\nUsing ", deparse(new.s.d.denom), " instead."))
                                        return(new.s.d.denom)})
        }
        else X$s.d.denom <- switch(tolower(A$estimand), att = "treated", ate = "pooled")
    }

    cluster <- A$cluster
    if (length(cluster) > 0) {
        if (!(is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster) > 1))) {
            stop("The argument to cluster must be a vector of cluster membership.", call. = FALSE)
        }
        if (length(cluster) != length(cbps.fit$y)) {
            stop("cluster must be the same length as the original data set.", call. = FALSE)
        }
    }

    if (length(cbps.fit$fitted.values) > 0) X$distance <- cbps.fit$fitted.values
    else X$distance <- NULL
    X$weights <- cbps.fit$weights
    X$treat <- cbps.fit$y
    X$covs <- cbps.fit$data[, !is.na(match(names(cbps.fit$data), attributes(terms(cbps.fit))$term.labels))]
    X$cluster <- factor(cluster)
    X$call <- cbps.fit$call
    X$obj <- list(treat = X$treat, weights = X$weights)
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
    return(X)
}