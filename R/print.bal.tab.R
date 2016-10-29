print.bal.tab <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", un = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    balance <- x$Balance
    baltal.m <- x$Balanced.Means
    maximbal.m <- x$Max.Imbalance.Means
    baltal.v <- x$Balanced.Variances
    maximbal.v <- x$Max.Imbalance.Variances
    nn <- x$Observations
    p.ops <- x$print.options
    
    #Adjustments to print options
    if (!p.ops$quick) {
        if (!identical(un, "as.is") && p.ops$disp.adj) {
            if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
            p.ops$un <- un
        }
        if (!identical(disp.means, "as.is")) {
            if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
            p.ops$disp.means <- disp.means
        }
        if (!identical(disp.v.ratio, "as.is")) {
            if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
            p.ops$disp.v.ratio <- disp.v.ratio
        }
    }
    if (!identical(disp.m.threshold, "as.is")) {
        if (!is.logical(disp.m.threshold)) stop("disp.m.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
            baltal.m <- NULL
            maximbal.m <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    
    
    keep <- c((1:ncol(balance))*c(TRUE, 
                                  p.ops$un*p.ops$disp.means, 
                                  p.ops$un*p.ops$disp.means, 
                                  p.ops$un, 
                                  p.ops$un*p.ops$disp.v.ratio, 
                                  p.ops$disp.adj*p.ops$disp.means, 
                                  p.ops$disp.adj*p.ops$disp.means, 
                                  p.ops$disp.adj, 
                                  !is.null(p.ops$m.threshold), 
                                  p.ops$disp.adj*p.ops$disp.v.ratio, 
                                  !is.null(p.ops$v.threshold)))
    
    round_df <- function(df, digits) {
        nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
        df[, nums] <- round(df[, nums], digits = digits)
        return(df)
    }
    replaceNA <- function(x) {
        x[is.na(x)] <- ""
        return(x)
    }
    
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n")
    }
    
    cat("\nBalance Measures:\n")
    print.data.frame(replaceNA(round_df(balance[, keep], digits)))
    
    if (!is.null(baltal.m)) {
        cat("\nBalance tally for mean differences:\n")
        print.data.frame(x$Balanced.Means)
    }
    if (!is.null(maximbal.m)) {
        cat("\nVariable with the greatest mean difference:\n")
        print.data.frame(round_df(x$Max.Imbalance.Means, digits))
    }
    if (!is.null(baltal.v)) {
        cat("\nBalance tally for variance ratios:\n")
        print.data.frame(x$Balanced.Variances, digits)
    }
    if (!is.null(maximbal.v)) {
        cat("\nVariable with the greatest variance ratio:\n")
        print.data.frame(round_df(x$Max.Imbalance.Variances, digits))
    }
    if (!is.null(nn)) {
        cat(paste0("\n", attr(x$Observations, "tag"), "\n"))
        print.data.frame(replaceNA(x$Observations), digits = digits)
    }
    invisible(x)
}
print.bal.tab.cont <- function(x, disp.r.threshold = "as.is", un = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    balance <- x$Balance
    baltal.r <- x$Balanced.Corr
    maximbal.r <- x$Max.Imbalance.Corr
    nn <- x$Observations
    p.ops <- x$print.options
    
    #Adjustments to print options
    if (!p.ops$quick) {
        if (!identical(un, "as.is") && p.ops$disp.adj) {
            if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
            p.ops$un <- un
        }
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
            baltal.r <- NULL
            maximbal.r <- NULL
        }
    }
    
    keep <- c((1:ncol(balance))*c(TRUE, 
                                  p.ops$un, 
                                  p.ops$disp.adj, 
                                  !is.null(p.ops$r.threshold)))
    
    round_df <- function(df, digits) {
        nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
        df[, nums] <- round(df[, nums], digits = digits)
        return(df)
    }
    replaceNA <- function(x) {
        x[is.na(x)] <- ""
        return(x)
    }
    
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n")
    }
    
    cat("\nBalance Measures:\n")
    print.data.frame(replaceNA(round_df(balance[, keep], digits)))
    
    if (!is.null(baltal.r)) {
        cat("\nBalance tally for correlations:\n")
        print.data.frame(x$Balanced.Corr)
    }
    if (!is.null(maximbal.r)) {
        cat("\nVariable with the greatest treatment correlation:\n")
        print.data.frame(round_df(x$Max.Imbalance.Corr, digits))
    }
    
    if (!is.null(nn)) {
        cat(paste0("\n", attr(x$Observations, "tag"), "\n"))
        print.data.frame(replaceNA(x$Observations), digits = digits)
    }
    invisible(x)
}
print.bal.tab.subclass <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", un = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", disp.subclass = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    s.balance <- x$Subclass.Balance
    b.a.subclass <- x$Balance.Across.Subclass
    baltal.m.subclass <- x$Balanced.Means.Subclass
    maximbal.m.subclass <- x$Max.Imbalance.Means.Subclass
    baltal.v.subclass <- x$Balanced.Variances.Subclass
    maximbal.v.subclass <- x$Max.Imbalance.Variances.Subclass
    s.nn <- x$Subclass.Observations
    p.ops <- x$print.options
    
    #Adjustments to print options
    if (!p.ops$quick) {
        if (!identical(un, "as.is") && p.ops$disp.adj) {
            if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
            p.ops$un <- un
        }
        if (!identical(disp.means, "as.is")) {
            if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
            p.ops$disp.means <- disp.means
        }
        if (!identical(disp.v.ratio, "as.is")) {
            if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
            p.ops$disp.v.ratio <- disp.v.ratio
        }
    }
    if (!identical(disp.m.threshold, "as.is")) {
        if (!is.logical(disp.m.threshold)) stop("disp.m.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
            baltal.m.subclass <- NULL
            maximbal.m.subclass <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
            baltal.v.subclass <- NULL
            maximbal.v.subclass <- NULL
        }
    }
    if (!identical(disp.subclass, "as.is")) {
        if (!is.logical(disp.subclass)) stop("disp.subclass must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.subclass <- disp.subclass
    }
    
    round_df <- function(df, digits) {
        nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
        df[, nums] <- round(df[, nums], digits = digits)
        return(df)
    }
    replaceNA <- function(x) {
        x[is.na(x)] <- ""
        return(x)
    }
    
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n")
    }
    
    if (p.ops$disp.subclass) {
        s.keep <- c( (1:7) * c(TRUE, 
                               p.ops$disp.means, 
                               p.ops$disp.means, 
                               p.ops$disp.adj, 
                               !is.null(p.ops$m.threshold), 
                               p.ops$disp.adj * p.ops$disp.v.ratio, 
                               !is.null(p.ops$v.threshold)))
        cat("\nBalance by subclass:")
        for (i in names(s.balance)) {
            cat(paste0("\n - - - Subclass ", i, " - - - \n"))
            print.data.frame(replaceNA(round_df(s.balance[[i]][, s.keep], digits)))
        }
    }
    
    a.s.keep <- c( (1:8) * c(TRUE, 
                             p.ops$un*p.ops$disp.means, 
                             p.ops$un*p.ops$disp.means, 
                             p.ops$un, 
                             p.ops$disp.adj*p.ops$disp.means, 
                             p.ops$disp.adj*p.ops$disp.means, 
                             p.ops$disp.adj, 
                             !is.null(p.ops$m.threshold)))
    cat("\nBalance measures across subclasses:\n")
    print.data.frame(replaceNA(round_df(b.a.subclass[, a.s.keep], digits)))
    
    if (!is.null(baltal.m.subclass)) {
        cat("\nBalance tally for mean differences across subclasses:\n")
        print.data.frame(baltal.m.subclass)
    }
    if (!is.null(maximbal.m.subclass)) {
        cat("\nVariable with the greatest mean difference across subclasses:\n")
        print.data.frame(round_df(maximbal.m.subclass, digits))
    }
    if (!is.null(baltal.v.subclass)) {
        cat("\nBalance tally for variance ratios across subclasses:\n")
        print.data.frame(baltal.v.subclass)
    }
    if (!is.null(maximbal.v.subclass)) {
        cat("\nVariable with the greatest variance ratios across subclasses:\n")
        print.data.frame(round_df(maximbal.v.subclass, digits))
    }
    
    if (!is.null(s.nn)) {
        cat(paste0("\n", attr(x$Subclass.Observations, "tag"), "\n"))
        print.data.frame(replaceNA(x$Subclass.Observations), digits = digits)
    }
    
    invisible(x)
}
print.bal.tab.cluster <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", un = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = NULL, digits = max(3, getOption("digits") - 3), ...) {
    #Figure out how to print bal.tab for clusters with subclassification
    call <- x$call
    c.balance <- x$Cluster.Balance
    c.balance.summary <- x$Cluster.Summary
    p.ops <- x$print.options
    
    #Adjustments to print options
    if (!p.ops$quick) {
        if (!identical(un, "as.is") && p.ops$disp.adj) {
            if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
            p.ops$un <- un
        }
        if (!identical(disp.means, "as.is")) {
            if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
            p.ops$disp.means <- disp.means
        }
        if (!identical(disp.v.ratio, "as.is")) {
            if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
            p.ops$disp.v.ratio <- disp.v.ratio
        }
        if (!identical(cluster.summary, "as.is")) {
            if (!is.logical(cluster.summary)) stop("cluster.summary must be TRUE, FALSE, or \"as.is\"")
            p.ops$cluster.summary <- cluster.summary
        }
    }
    if (!identical(disp.m.threshold, "as.is")) {
        if (!is.logical(disp.m.threshold)) stop("disp.m.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
        }
    }
    if (!missing(which.cluster)) {
        p.ops$which.cluster <- which.cluster
    }
    
    cluster.funs <- c("min", "mean", "median", "max")
    if (length(cluster.fun) > 0) {
        cluster.fun <- match.arg(tolower(cluster.fun), cluster.funs, several.ok = TRUE)
        if (!is.na(match("min", cluster.fun)) && p.ops$quick) {
            cluster.fun <- cluster.fun[is.na(match(cluster.fun, "min"))]
            if (length(cluster.fun) == 0) {
                warning("\"min\" was the only valid entry to cluster.fun, but it cannot be used when quick = TRUE in the call to bal.tab(). Using all other cluster.funs instead.", call. = FALSE)
                cluster.fun <- c("mean", "median", "max")
            }
            else warning("\"min\" cannot be requested when quick = TRUE in the call to bal.tab(), so it was ignored.", call. = FALSE)
        }
        else if (length(cluster.fun) == 0) {
            warning("There were no valid entries to cluster.fun. Using all other cluster.funs instead.", call. = FALSE)
            cluster.fun <- c("mean", "median", "max")
        }
    }
    else cluster.fun <- c("mean", "median", "max")
    
    #Checks and Adjustments
    if (length(p.ops$which.cluster) == 0) 
        which.cluster <- seq_along(c.balance)
    else if (is.numeric(p.ops$which.cluster)) {
        which.cluster <- seq_along(c.balance)[seq_along(c.balance) %in% p.ops$which.cluster]
        if (length(which.cluster) == 0) {
            warning("No indices in which.cluster are cluster indices. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(c.balance)
        }
    }
    else if (is.character(p.ops$which.cluster)) {
        which.cluster <- seq_along(c.balance)[names(c.balance) %in% p.ops$which.cluster]
        if (length(which.cluster) == 0) {
            warning("No names in which.cluster are cluster names. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(c.balance)
        }
    }
    else if (is.na(p.ops$which.cluster)) {
        which.cluster <- integer(0)
        p.ops$cluster.summary <- TRUE
    }
    else {
        warning("The argument to which.cluster must be NA, NULL, or a vector of cluster indicies or cluster names. Displaying all clusters instead.", call. = FALSE)
        which.cluster <- seq_along(c.balance)
    }
    
    round_df <- function(df, digits) {
        nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
        df[, nums] <- round(df[, nums], digits = digits)
        return(df)
    }
    replaceNA <- function(x) {
        x[is.na(x)] <- ""
        return(x)
    }
    
    keep <- c((1:ncol(c.balance[[1]]))*c(TRUE, 
                                         p.ops$un*p.ops$disp.means, 
                                         p.ops$un*p.ops$disp.means, 
                                         p.ops$un, 
                                         p.ops$un*p.ops$disp.v.ratio, 
                                         p.ops$disp.adj*p.ops$disp.means, 
                                         p.ops$disp.adj*p.ops$disp.means, 
                                         p.ops$disp.adj, 
                                         !is.null(p.ops$m.threshold), 
                                         p.ops$disp.adj*p.ops$disp.v.ratio, 
                                         !is.null(p.ops$v.threshold)))
    
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n")
    }
    
    if (length(which.cluster)>0) {
        cat("\nBalance by cluster:")
        for (i in which.cluster) {
            cat(paste0("\n - - - Cluster: ", names(c.balance)[i], " - - - \n"))
            print.data.frame(replaceNA(round_df(c.balance[[i]][, keep], digits)))
        }
    }
    
    if (isTRUE(as.logical(p.ops$cluster.summary))) {
        cf <- !is.na(match(cluster.funs, cluster.fun))
        names(cf) <- cluster.funs
        s.keep <- (1:ncol(c.balance.summary))*c(TRUE, 
                                                p.ops$un*cf["min"],
                                                p.ops$un*cf["mean"],
                                                p.ops$un*cf["median"],
                                                p.ops$un*cf["max"],
                                                p.ops$un*p.ops$disp.v.ratio*cf["min"],
                                                p.ops$un*p.ops$disp.v.ratio*cf["mean"],
                                                p.ops$un*p.ops$disp.v.ratio*cf["median"],
                                                p.ops$un*p.ops$disp.v.ratio*cf["max"],
                                                p.ops$disp.adj*cf["min"],
                                                p.ops$disp.adj*cf["mean"],
                                                p.ops$disp.adj*cf["median"],
                                                p.ops$disp.adj*cf["max"],
                                                p.ops$disp.adj*p.ops$disp.v.ratio*cf["min"],
                                                p.ops$disp.adj*p.ops$disp.v.ratio*cf["mean"],
                                                p.ops$disp.adj*p.ops$disp.v.ratio*cf["median"],
                                                p.ops$disp.adj*p.ops$disp.v.ratio*cf["max"])
        cat("\nBalance summary across all clusters:\n")
        print.data.frame(replaceNA(round_df(c.balance.summary[, s.keep], digits)))
    }
    invisible(x)
}
print.bal.tab.cont.cluster <- function(x, disp.r.threshold = "as.is", un = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = NULL, digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    c.balance <- x$Cluster.Balance
    c.balance.summary <- x$Cluster.Summary
    p.ops <- x$print.options
    
    #Adjustments to print options
    if (!p.ops$quick) {
        if (!identical(un, "as.is") && p.ops$disp.adj) {
            if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
            p.ops$un <- un
        }
        if (!identical(cluster.summary, "as.is")) {
            if (!is.logical(cluster.summary)) stop("cluster.summary must be TRUE, FALSE, or \"as.is\"")
            p.ops$cluster.summary <- cluster.summary
        }
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
            baltal.r <- NULL
            maximbal.r <- NULL
        }
    }
    if (!missing(which.cluster)) {
        p.ops$which.cluster <- which.cluster
    }
 
    cluster.funs <- c("min", "mean", "median", "max")
    if (length(cluster.fun) > 0) {
        cluster.fun <- match.arg(tolower(cluster.fun), cluster.funs, several.ok = TRUE)
        if (!is.na(match("min", cluster.fun)) && p.ops$quick) {
            cluster.fun <- cluster.fun[is.na(match(cluster.fun, "min"))]
            if (length(cluster.fun) == 0) {
                warning("\"min\" was the only valid entry to cluster.fun, but it cannot be used when quick = TRUE in the call to bal.tab(). Using all other cluster.funs instead.", call. = FALSE)
                cluster.fun <- c("mean", "median", "max")
            }
            else warning("\"min\" cannot be requested when quick = TRUE in the call to bal.tab(), so it was ignored.", call. = FALSE)
        }
        else if (length(cluster.fun) == 0) {
            warning("There were no valid entries to cluster.fun. Using all other cluster.funs instead.", call. = FALSE)
            cluster.fun <- c("mean", "median", "max")
        }
    }
    else cluster.fun <- c("mean", "median", "max")
    
    #Checks and Adjustments
    if (length(p.ops$which.cluster) == 0) 
        which.cluster <- seq_along(c.balance)
    else if (is.numeric(p.ops$which.cluster)) {
        which.cluster <- seq_along(c.balance)[seq_along(c.balance) %in% p.ops$which.cluster]
        if (length(which.cluster) == 0) {
            warning("No indices in which.cluster are cluster indices. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(c.balance)
        }
    }
    else if (is.character(p.ops$which.cluster)) {
        which.cluster <- seq_along(c.balance)[names(c.balance) %in% p.ops$which.cluster]
        if (length(which.cluster) == 0) {
            warning("No names in which.cluster are cluster names. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(c.balance)
        }
    }
    else if (is.na(p.ops$which.cluster)) {
        which.cluster <- integer(0)
        p.ops$cluster.summary <- TRUE
    }
    else {
        warning("The argument to which.cluster must be NA, NULL, or a vector of cluster indicies or cluster names. Displaying all clusters instead.", call. = FALSE)
        which.cluster <- seq_along(c.balance)
    }
    
    round_df <- function(df, digits) {
        nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
        df[, nums] <- round(df[, nums], digits = digits)
        return(df)
    }
    replaceNA <- function(x) {
        x[is.na(x)] <- ""
        return(x)
    }
    
    keep <- c((1:ncol(c.balance[[1]]))*c(TRUE, 
                                         p.ops$un, 
                                         p.ops$disp.adj, 
                                         !is.null(p.ops$r.threshold)))
    
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n")
    }
    
    if (length(which.cluster)>0) {
        cat("\nBalance by cluster:")
        for (i in which.cluster) {
            cat(paste0("\n - - - Cluster: ", names(c.balance)[i], " - - - \n"))
            print.data.frame(replaceNA(round_df(c.balance[[i]][, keep], digits)))
        }
    }
    
    if (isTRUE(as.logical(p.ops$cluster.summary))) {
        cf <- !is.na(match(cluster.funs, cluster.fun))
        names(cf) <- cluster.funs
        s.keep <- (1:ncol(c.balance.summary))*c(TRUE, 
                                                p.ops$un*cf["min"],
                                                p.ops$un*cf["mean"],
                                                p.ops$un*cf["median"],
                                                p.ops$un*cf["max"],
                                                p.ops$disp.adj*cf["min"],
                                                p.ops$disp.adj*cf["mean"],
                                                p.ops$disp.adj*cf["median"],
                                                p.ops$disp.adj*cf["max"])
        cat("\nBalance summary across all clusters:\n")
        print.data.frame(replaceNA(round_df(c.balance.summary[, s.keep], digits)))
    }
    invisible(x)
}