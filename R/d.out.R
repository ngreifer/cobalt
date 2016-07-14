d.out <- function(obj, ...) {
    #Produces data.frame of observations after adjustment. Similar to MatchIt's match.data
    #obj is the first input to bal.tab
    #
    #ps: full.stop.method (optional, will use first if blank)
    #Match: formula and data or treat and covs
    #formula: data, weights, distance (optional), method (optional, will default to matching)
    #data.frame: treat, data (optional; for naming treat and weights), weights, distance (optional), method (optional, will default to matching)
    
    args <- list(...)
    if (any(class(obj)=="matchit")) X <- matchit2base(obj)
    else if (any(class(obj)=="ps")) X <- ps2base(obj, full.stop.method = args$full.stop.method)
    else if (any(class(obj)=="Match")) X <- Match2base(obj, formula = args$formula, data = args$data, treat = args$treat, covs = args$covs)
    else if (any(class(obj)=="formula")) {
        X0 <- formula2df(obj, data = args$data)
        X <- df2base(X0$covs, X0$treat, data = args$data, weights = args$weights, distance = args$distance, method = args$method)
    }
    else if (is.data.frame(obj)) X <- df2base(obj, treat = args$treat, data = args$data, weights = args$weights, distance = args$distance, method = args$method)
    
    out <- data.frame(X$covs, treat = X$treat, weights = X$weights)
    if (!is.null(X$distance)) out$distance <- X$distance
    if (!is.null(X$subclass)) out$subclass <- X$subclass
    out <- out[, !sapply(out, function(x) all(is.na(x)))]
    return(out)
}