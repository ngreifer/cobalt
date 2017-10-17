d.out <- function(obj, ..., clean = FALSE) {
    #Produces data.frame of observations after adjustment. Similar to MatchIt's match.data
    #obj is the first input to bal.tab
    #
    #ps: full.stop.method (optional, will use first if blank)
    #Match: formula and data or treat and covs
    #formula: data, weights, distance (optional), method (optional, will default to weighting)
    #data.frame: treat, data (optional; for naming treat and weights), weights, distance (optional), method (optional, will default to matching)
    
    args <- list(...)
    X <- x2base(obj, ..., cluster = cluster, std.ok = TRUE)
    if (clean) 
    out <- data.frame(X$covs, treat = X$treat, weights = X$weights)
    if (!is.null(X$distance)) out$distance <- X$distance
    if (!is.null(X$subclass)) out$subclass <- X$subclass
    out <- out[, !sapply(out, function(x) all(is.na(x)))]
    return(out)
}