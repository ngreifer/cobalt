#Functions to convert object to base.bal.tab input
matchit2base <- function(m) {
    X <- list(covs=NA, 
              treat=NA, 
              weights=NA, 
              subclass=NA, 
              method=NA, 
              distance=NA, 
              obj=NA,
              call=NA)
    #Initializing variables
    if (any(class(m)=="matchit.subclass")) {
        m$subclass <- factor(m$subclass)
        X$subclass <- m$subclass
        X$method <- "subclassification"
    }
    else {
        X$subclass <- NULL
        X$method <- "matching"
    }
    X$weights <- m$weights
    X$treat <- m$treat
    if (!all(is.na(m$distance))) X$distance <- m$distance
    else X$distance <- NULL
    o.data <- data.frame(m$model$model) #Just the data used in the PS, including treatment and covs
    X$covs <- o.data[, which(names(o.data) %in% attributes(terms(m$model))$term.labels)]
    #X$covs <- data.frame(m$X)
    X$obj <- m
    
    X$call <- m$call
    return(X)
}
ps2base <- function(ps, full.stop.method, s.d.denom) {
    X <- list(covs=NA, 
              treat=NA, 
              weights=NA, 
              distance=NA, 
              s.d.denom=NA, 
              call=NA, 
              obj=NA)

        if (exists("full.stop.method", mode=c("character"))) {
        rule1 <- tryCatch(match.arg(tolower(full.stop.method), tolower(names(ps$w))),
                                    #error = function() {message(paste0("Warning: '", full.stop.method, "' is not the name of a stop.method in ", deparse(substitute(ps, parent.frame(1))), ". Will use '", names(ps$w)[1], "' instead."));
                                    error = function(cond) {message(paste0("Warning: full.stop.method should be one of ", paste(shQuote(names(ps$w)), collapse = ", "), ".\nUsing  '", names(ps$w)[1], "' instead.")); 
                                                                return(names(ps$w)[1])})
    }
    else {
        rule1 <- names(ps$w)[1]
    }

    #Initializing variables
    s <- names(ps$w)[match(tolower(rule1), tolower(names(ps$w)))]

    attr(ps, "which") <- s
    X$s.d.denom <- ifelse(!exists("s.d.denom", mode=c("character")), switch(substr(tolower(s), nchar(s)-2, nchar(s)), att = "treated", ate = "pooled"), match.arg(s.d.denom, c("treated", "control", "pooled")))
    X$weights <- as.matrix(ps$w[s])
    X$treat <- ps$treat
    X$distance <- ps$ps[s][, ]
    X$covs <- ps$data[ps$gbm.obj$var.names]
    X$call <- ps$parameters
    X$obj <- list(treat=ps$treat, weights=ps$w[s][, ])
    return(X)
}
Match2base <- function(Match, formula=NULL, data=NULL, treat=NULL, covs=NULL, s.d.denom) {
    X <- list(covs=NA, 
              treat=NA, 
              weights=NA, 
              method=NA, 
              distance=NA, 
              call=NA, 
              s.d.denom=NA, 
              obj=NA)
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
        
        if (sum(good) %in% c(0, 1)) {
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
    use.which <- useWhich(formula, data, treat, covs)
    
    use.fd <- function(formula, data){
        #outputs a list containing treat [1] and covs [2]
        out.list <- list(treat=NA, covs=data.frame(NA))
        tt <- terms(formula)
        attr(tt, "intercept") <- 0
        m.try <- try({mf <- model.frame(tt, data)}, TRUE) 
        if (class(m.try) == "try-error") {
            stop(paste0(c("All right hand side variables of formula must be variables in data.\nVariables not in data: ", 
                          paste(attr(tt, "term.labels")[which(!attr(tt, "term.labels") %in% names(data))], collapse=", "))), call. = FALSE)}
        out.list$treat <- model.response(mf) #treat
        out.list$covs <- as.data.frame(model.matrix(tt, data=mf)) #covs
        return(out.list)
    }
    use.tc <- function(treat, covs) {
        if (length(treat)!=nrow(covs)) {
            stop("treat must be same length as covs", call. = FALSE)}
        out.list <- list(treat=treat, covs=covs)
        return(out.list)
    }
    
    if (use.which=="fd") {
        t.c <- use.fd(formula, data)}
    else if (use.which=="tc") {
        t.c <- use.tc(treat, covs)}
    else if (use.which=="both") {
        try.fd <- try({t.c <- use.fd(formula, data)})
        if (class(try.fd)=="try-error") {
            message("Formula, data, treat, and covs all supplied; ignoring formula and data.")
            t.c <- use.tc(treat, covs)}
        else {
            message("Formula, data, treat, and covs all supplied; ignoring treat and covs.")}
    }
    
    if (sum(is.na(t.c$covs))>0)
        stop("Missing values exist in the covariates", call. = FALSE)
    
    #Initializing variables
    m <- Match
    s <- m$estimand
    X$s.d.denom = ifelse(!exists("s.d.denom", mode=c("character")), switch(toupper(s), ATT = "treated", ATE = "treated", ATC = "control"), match.arg(s.d.denom, c("treated", "control", "pooled")))
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
    
    # distance.list$control <- m$mdata$X[2, ]
    # distance.list$treated <- m$mdata$X[1, ]
    # distance.list$unmatched <- rep(NA, length(treat0[!(1:nobs) %in% c(m$index.treated, m$index.control, m$index.dropped)]))
    # distance.list$dropped <- rep(NA, length(m$index.dropped))
    
    #data.list <- lapply(1:4, function(x) cbind(data.frame(treat=treat.list[[x]]), data.frame(weights=weights.list[[x]]), data.frame(distance=distance.list[[x]]), covs.list[[x]]))
    data.list <- lapply(1:4, function(x) cbind(data.frame(treat=treat.list[[x]]), data.frame(weights=weights.list[[x]]), covs.list[[x]]))
    o.data <- do.call(rbind, data.list)
    o.data2 <- merge(unique(o.data[, !names(o.data) %in% "weights"]), aggregate(weights~index, data=o.data, FUN=sum), by="index")
    
    X$treat <- o.data2$treat
    X$weights <- o.data2$weights
    X$distance <- NULL #NAs in distance bcause of incomplete list in Match object
    X$covs <- o.data2[, !names(o.data2) %in% c("treat", "weights", "index")]
    X$call <- NULL
    X$method <- "matching"
    X$obj <- list(treat=X$treat, weights=X$weights)
    return(X)
}
formula2df <- function(formula, data) {
    X <- list(treat=NA, 
              covs=NA)
    #Checks
    if (is.null(data)) {
        stop("Dataframe must be specified", call. = FALSE)}
    if (!is.data.frame(data)) {
        stop("Data must be a dataframe", call. = FALSE)}
    
    #Initializing variables
    tt <- terms(formula)
    attr(tt, "intercept") <- 0
    m.try <- try({mf <- model.frame(tt, data)}, TRUE) 
    if (class(m.try) == "try-error") {
        stop(paste0(c("All right hand side variables of formula must be variables in data.\nVariables not in data: ", 
                      paste(attr(tt, "term.labels")[which(!attr(tt, "term.labels") %in% names(data))], collapse=", "))), call. = FALSE)}
    X$treat <- model.response(mf)
    X$covs <- as.data.frame(model.matrix(tt, data=mf))
    return(X)
}
df2base <- function(covs, treat, data=NULL, weights=NULL, distance=NULL, subclass=NULL, addl=NULL, method) {
    X <- list(covs=NA, 
              weights=NA, 
              treat=NA, 
              distance=NA, 
              subclass=NA, 
              addl=NA, 
              method=NA, 
              call=NA, 
              obj=NA)
    
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
    if (!is.null(data) & !is.data.frame(data)) {
        warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execuction will halt.")
        data <- NULL
    }
    if (!any(is.null(weights), exists("weights", mode=c("character")), exists("weights", mode=c("numeric")))) {
        stop("The argument to weights must be a vector of weights or the (quoted) name of a variable in data that contains weights.", call. = FALSE)
    }
    if (!any(is.null(distance), exists("distance", mode=c("character")), exists("distance", mode=c("numeric")))) {
        stop("The argument to distance must be a vector of distance scores or the (quoted) name of a variable in data that contains distance scores.", call. = FALSE)
    }
    if (!any(is.null(subclass), exists("subclass", mode=c("character")), exists("subclass", mode=c("numeric")))) {
        stop("The argument to subclass must be a vector of subclass membership or the (quoted) name of a variable in data that contains subclass membership.", call. = FALSE)
    }
    if (!any(exists("treat", mode=c("character")), exists("treat", mode=c("numeric")))) {
        stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status", call. = FALSE)
    }
    
    if (is.numeric(get0("treat"))) {
        treat <- treat
    }
    else if (is.character(get0("treat")) & !is.null(data) & get0("treat") %in% names(data)) {
        treat <- data[, get0("treat")]
    }
    else stop("The name supplied to treat is not the name of a variable in data.")
    
    if (length(treat) != nrow(covs)) {
        stop("treat must be same length as covs", call. = FALSE)}
    
    if (sum(is.na(treat)) > 0)
        stop("Missing values exist in treat", call. = FALSE)
    
    if (!is.null(weights)) {
        if (is.numeric(get0("weights"))) {
            weights <- weights
        }
        else if (is.character(get0("weights")) & get0("weights") %in% names(data)) {
            weights <- data[, get0("weights")]
        }
        else stop("The name supplied to weights is not the name of a variable in data.")
        
        if (length(weights) != nrow(covs)) {
            stop("weights must be same length as covs", call. = FALSE)
        }
        
        if (sum(is.na(weights)) > 0)
            stop("Missing values exist in weights", call. = FALSE)
    }
    
    if (!is.null(distance)) {
        if (is.numeric(get0("distance"))) {
            distance <- distance
        }
        else if (is.character(get0("distance")) & get0("distance") %in% names(data)) {
            distance <- data[, get0("distance")]
        }
        else stop("The name supplied to distance is not the name of a variable in data.")
        
        if (length(distance)!=nrow(covs)) {
            stop("distance must be same length as covs", call. = FALSE)
        }
        
        if (sum(is.na(distance))>0)
            stop("Missing values exist in distance", call. = FALSE)
    }
    
    if (!is.null(subclass)) {
        if (is.numeric(get0("subclass"))) {
            subclass <- subclass
        }
        else if (is.character(get0("subclass")) & get0("subclass") %in% names(data)) {
            subclass <- data[, get0("subclass")]
        }
        else stop("The name supplied to subclass is not the name of a variable in data.")
        
        if (length(subclass) != nrow(covs)) {
            stop("subclass must be same length as covs", call. = FALSE)
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
            if (nrow(addl)!=nrow(data)) {
                stop("If addl is a data.frame, it must have the same number of rows as data.")
            }
        }
        else {
            warning("addl must be a list of names of variables in data or a data.frame containing additional variable(s). addl will be ignored in the following output.")
            addl <- NULL
        }
    }
    
    if (!is.null(weights)) {
        if (exists("method", mode="character")) {
            X$method <- match.arg(method, c("weighting", "matching", "subclassification"))
        }
        else {
            message("Assuming weights generated through weighting; if not, please specify with argument to method.")
            X$method <- "weighting"
        }
    }
    else if (!is.null(subclass)){
        X$method <- "subclassification"
    }
    else X$method <- "matching"
    
    X$covs <- covs
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$subclass <- factor(subclass)
    X$call <- NULL
    X$addl <- addl
    X$obj <- data.frame(treat=treat, weights=NA)
    if (!is.null(weights)) X$obj$weights <- weights
    if (!is.null(subclass)) X$obj$subclass <- factor(subclass)
    return(X)
}
CBPS2base <- function(cbps.fit, estimand=NULL, s.d.denom, std.ok = FALSE) {
    X <- list(covs=NA, 
              treat=NA, 
              weights=NA, 
              distance=NA, 
              s.d.denom=NA, 
              call=NA, 
              obj=NA)
    #Checks
    if (!std.ok && sum(cbps.fit$weights) < 3) {
        if (is.null(estimand) && !exists("s.d.denom", mode=c("character"))) stop("Please specify either the estimand (\"ATT\" or \"ATE\") or an argument to s.d.denom.")
        warning("Standardized weights were used; this may cause reported values to be incorrect. Use unstandardized weights instead.", call. = FALSE)
    }
    else {
        if (isTRUE(all.equal(cbps.fit$weights, cbps.fit$y / cbps.fit$fitted.values + (1-cbps.fit$y) / (1-cbps.fit$fitted.values)))) estimand <- "ATE"
        else estimand <- "ATT"
    }
    X$s.d.denom <- ifelse(!exists("s.d.denom", mode=c("character")), switch(tolower(estimand), att = "treated", ate = "pooled"), match.arg(s.d.denom, c("treated", "control", "pooled")))
    if (!is.null(cbps.fit$fitted.values)) X$distance <- cbps.fit$fitted.values
    else X$distance <- NULL
    X$weights <- cbps.fit$weights
    X$treat <- cbps.fit$y
    X$distance <- cbps.fit$fitted.values
    X$covs <- cbps.fit$data[, which(names(cbps.fit$data) %in% attributes(terms(cbps.fit))$term.labels)]
    X$call <- cbps.fit$call
    X$obj <- list(treat = cbps.fit$y, weights = cbps.fit$weights)
    return(X)
}