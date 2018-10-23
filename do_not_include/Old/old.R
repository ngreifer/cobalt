#Retired functions and versions of functions

.diff.selector <- function(x, treat, weights = NULL, subclass = NULL, which.sub = NULL, x.type, continuous, binary, s.d.denom, no.weights = FALSE, s.weights = rep(1, length(treat))) {
    if (no.weights) weights <- rep(1, length(x))
    no.sub <- length(which.sub) == 0
    if (no.sub) ss <- rep(TRUE, length(x))
    else ss <- (subclass == which.sub & weights > 0 & s.weights > 0)
    
    # if (x.type=="Distance")  {
    #     diff <- w.m(x[ss & treat==1], w = weights[ss & treat==1]) - w.m(x[ss & treat==0], w = weights[ss & treat==0])
    # }
    if (x.type=="Binary") {
        if      (binary=="raw") diff <- w.m(x[ss & treat==1], w = weights[ss & treat==1]) - w.m(x[ss & treat==0], w = weights[ss & treat==0])
        else if (binary=="std") {
            if (no.sub) diff <- std.diff(x, treat = treat, weights, denom = s.d.denom, s.weights = s.weights)
            else diff <- std.diff.subclass(x, treat = treat, weights, subclass, which.sub, denom = s.d.denom)
        }
        else diff <- NULL
    }
    else if (any(c("Contin.", "Distance") == x.type)) {
        if      (continuous=="raw") diff <- w.m(x[ss & treat==1], w = weights[ss & treat==1]) - w.m(x[ss & treat==0], w = weights[ss & treat==0])
        else if (continuous=="std") {
            if (no.sub) diff <- std.diff(x, treat = treat, weights, denom = s.d.denom, s.weights = s.weights)
            else diff <- std.diff.subclass(x, treat = treat, weights, subclass, which.sub, denom = s.d.denom)
        }
        else diff <- NULL
    }
    return(diff)
}  

.var.ratio <- function(x, treat, weights, var.type, no.weights = FALSE) {
    if (no.weights) weights <- rep(1, length(x))
    if (var.type != "Binary") {
        ratio <- w.v(x[treat==1], weights[treat==1]) / w.v(x[treat==0], weights[treat==0])
        return(max(ratio, 1 / ratio))
    }
    else return(NA)
}

.std.diff <- function(x, treat, weights, denom, s.weights = rep(1, length(treat))) {
    #Uses variance of original group, as in MatchIt.
    #treated <- as.logical(treat)
    g0 <- x[treat==0];     g1 <- x[treat==1];
    w0 <- weights[treat==0]*s.weights[treat==0]; w1 <- weights[treat==1]*s.weights[treat==1]
    m0 <- w.m(g0, w0);        m1 <- w.m(g1, w1)
    v0 <- w.v(g0, s.weights[treat==0]);           v1 <- w.v(g1, s.weights[treat==1])
    if (!all(is.finite(m0)) || !all(is.finite(m1))) {
        stop("There is an error in your weights. This may also be due to a bug.", call. = FALSE)
    }
    
    s.d <- switch(denom, 
                  control = sqrt(v0), 
                  treated = sqrt(v1), 
                  pooled = sqrt((v0+v1)/2))
    m.dif <- m1 - m0
    
    if (abs(m.dif) < sqrt(.Machine$double.eps)) s.diff <- 0
    else s.diff <- m.dif/s.d
    if (!is.finite(s.diff)) s.diff <- NA
    return(s.diff)
}

.std.diff.subclass <- function(x, treat, weights, subclass, which.sub, denom) {
    #treated <- as.logical(treat)
    if (sum(treat==0 & !is.na(subclass) & subclass==which.sub & weights==1) == 0) {
        warning(paste0("There are no control units in subclass ", which.sub, "."), call. = FALSE)
        return(NA)
    }
    if (sum(treat==1 & !is.na(subclass) & subclass==which.sub & weights==1) == 0) {
        warning(paste0("There are no treated units in subclass ", which.sub, "."), call. = FALSE)
        return(NA)
    }
    g0 <- x[treat==0 & !is.na(subclass) & subclass==which.sub & weights==1]
    g1 <- x[treat==1 & !is.na(subclass) & subclass==which.sub & weights==1]
    m0 <- sum(g0)/length(g0)
    m1 <- sum(g1)/length(g1)
    v0 <- var(x[treat==0])
    v1 <- var(x[treat==1])
    if (!all(is.finite(m0)) || !all(is.finite(m1))) {
        stop("There is an error in your weights. This may also be due to a bug.", call. = FALSE)
    }
    s.d <- switch(denom, 
                  control = sqrt(v0), 
                  treated = sqrt(v1), 
                  pooled = sqrt((v0+v1)/2))
    m.dif <- m1 - m0
    if (abs(m.dif) < sqrt(.Machine$double.eps)) s.diff <- 0
    else s.diff <- m.dif/s.d
    if (!is.finite(s.diff)) s.diff <- NA
    return(s.diff)
}

.col.diff.selector <- function(mat, treat, weights = NULL, subclass = NULL, which.sub = NULL, x.types, continuous, binary, s.d.denom, no.weights = FALSE, s.weights = rep(1, length(treat))) {
    if (no.weights) weights <- rep(1, nrow(mat))
    no.sub <- length(which.sub) == 0
    if (no.sub) ss <- rep(TRUE, nrow(mat))
    else ss <- (subclass == which.sub & weights > 0 & s.weights > 0)
    
    diffs <- rep(NA, ncol(mat))
    if (binary == "raw") {
        diffs[x.types == "Binary"] <- col.w.m(mat[ss & treat == 1, x.types == "Binary"], w = weights[ss & treat==1]) - col.w.m(mat[ss & treat == 0, x.types == "Binary"], w = weights[ss & treat==0])
    }
    else if (binary=="std") {
        if (no.sub) diffs[x.types == "Binary"] <- apply(mat[, x.types == "Binary"], 2, function(x) std.diff(x, treat = treat, weights, denom = s.d.denom, s.weights = s.weights))
        else diffs[x.types == "Binary"] <- apply(mat[, x.types == "Binary"], 2, function(x) std.diff.subclass(x, treat = treat, weights, subclass, which.sub, denom = s.d.denom))
    }
    
    if (continuous == "raw") {
        diffs[x.types != "Binary"] <- col.w.m(mat[ss & treat == 1, x.types != "Binary"], w = weights[ss & treat==1]) - col.w.m(mat[ss & treat == 0, x.types != "Binary"], w = weights[ss & treat==0])
    }
    else if (continuous == "std") {
        if (no.sub) diffs[x.types != "Binary"] <- apply(mat[, x.types != "Binary"], 2, function(x) std.diff(x, treat = treat, weights, denom = s.d.denom, s.weights = s.weights))
        else diffs[x.types != "Binary"] <- apply(mat[, x.types != "Binary"], 2, function(x) std.diff.subclass(x, treat = treat, weights, subclass, which.sub, denom = s.d.denom))
    }
    
    return(diffs)
}  

.ks <- function(x, treat, weights = NULL, var.type, no.weights = FALSE) {
    #Computes ks-statistic
    if (no.weights) weights <- rep(1, length(x))
    if (var.type != "Binary") {
        # if (TRUE) {
        weights[treat == 1] <- weights[treat==1]/sum(weights[treat==1])
        weights[treat == 0] <- -weights[treat==0]/sum(weights[treat==0])
        
        ordered.index <- order(x)
        cumv <- abs(cumsum(weights[ordered.index]))[diff(x[ordered.index]) != 0]
        ks <- ifelse(length(cumv) > 0, max(cumv), 0)
        
        return(ks)
    }
    else return(NA)
}

.bal.tab.optmatch <- function(optmatch, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base", c(list(optmatch), args), quote = TRUE)
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
