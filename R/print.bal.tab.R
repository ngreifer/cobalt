print.bal.tab <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.r.threshold = "as.is", un = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    balance <- x$Balance
    baltal.r <- x$Balanced.Corr
    maximbal.r <- x$Max.Imbalance.Corr
    baltal.m <- x$Balanced.Means
    maximbal.m <- x$Max.Imbalance.Means
    baltal.v <- x$Balanced.Variances
    maximbal.v <- x$Max.Imbalance.Variances
    nn <- x$Observations
    p.ops <- x$print.options
    
    #Prevent expnential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
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
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
            baltal.r <- NULL
            maximbal.r <- NULL
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
    
    if (!is.na(match("bal.tab.cont", class(x)))) {
        keep <- as.logical(c(TRUE, 
                             p.ops$un, 
                             p.ops$disp.adj, 
                             !is.null(p.ops$r.threshold)))
    }
    else {
        keep <- as.logical(c(TRUE, 
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
        for (i in rownames(x$Observations)) {
            if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
        }
        if (!is.na(x$Observations["Matched (Unweighted)",]) && all(x$Observations["Matched",] == x$Observations["Matched (Unweighted)",])) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)",]
        cat(paste0("\n", attr(x$Observations, "tag"), ":\n"))
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
    
    #Prevent expnential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
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
    
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n")
    }
    
    if (p.ops$disp.subclass) {
        s.keep <- as.logical(c(TRUE, 
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
    
    a.s.keep <- as.logical(c(TRUE, 
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
        cat(paste0("\n", attr(x$Subclass.Observations, "tag"), ":\n"))
        print.data.frame(replaceNA(x$Subclass.Observations), digits = digits)
    }
    
    invisible(x)
}
print.bal.tab.cluster <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.r.threshold = "as.is", un = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = NULL, digits = max(3, getOption("digits") - 3), ...) {
    #Figure out how to print bal.tab for clusters with subclassification
    call <- x$call
    c.balance <- x$Cluster.Balance
    c.balance.summary <- x$Cluster.Summary
    nn <- x$Observations
    p.ops <- x$print.options
    
    #Prevent expnential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
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
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
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
    
    if (!is.na(match("bal.tab.cont.cluster", class(x)))) {
        keep <- as.logical(c(TRUE, 
                             p.ops$un, 
                             p.ops$disp.adj, 
                             !is.null(p.ops$r.threshold)))
    }
    else {
        keep <- as.logical(c(TRUE, 
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
    }
    
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n")
    }
    
    if (length(which.cluster)>0) {
        cat("Balance by cluster:\n")
        for (i in which.cluster) {
            cat(paste0("\n - - - Cluster: ", names(c.balance)[i], " - - - \n"))
            cat("Balance measures:\n")
            print.data.frame(replaceNA(round_df(c.balance[[i]][["Balance.Table"]][, keep], digits)))
            
            for (j in rownames(c.balance[[i]][["Observations"]])) {
                if (all(c.balance[[i]][["Observations"]][j,] == 0)) c.balance[[i]][["Observations"]] <- c.balance[[i]][["Observations"]][rownames(c.balance[[i]][["Observations"]])!=j,]
            }
            if (!is.na(c.balance[[i]][["Observations"]]["Matched (Unweighted)",]) && all(c.balance[[i]][["Observations"]]["Matched",] == c.balance[[i]][["Observations"]]["Matched (Unweighted)",])) c.balance[[i]][["Observations"]] <- c.balance[[i]][["Observations"]][rownames(c.balance[[i]][["Observations"]])!="Matched (Unweighted)",]
            cat(paste0("\n", attr(c.balance[[i]][["Observations"]], "tag"), ":\n"))
            print.data.frame(replaceNA(round_df(c.balance[[i]][["Observations"]], digits)))
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Cluster: ", names(c.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
    }
    
    if (isTRUE(as.logical(p.ops$cluster.summary))) {
        CF <- !is.na(match(cluster.funs, cluster.fun))
        names(CF) <- cluster.funs
        if (!is.na(match("bal.tab.cont.cluster", class(x)))) {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un*CF["min"],
                                   p.ops$un*CF["mean"],
                                   p.ops$un*CF["median"],
                                   p.ops$un*CF["max"],
                                   p.ops$disp.adj*CF["min"],
                                   p.ops$disp.adj*CF["mean"],
                                   p.ops$disp.adj*CF["median"],
                                   p.ops$disp.adj*CF["max"]))
        }
        else {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un*CF["min"],
                                   p.ops$un*CF["mean"],
                                   p.ops$un*CF["median"],
                                   p.ops$un*CF["max"],
                                   p.ops$un*p.ops$disp.v.ratio*CF["min"],
                                   p.ops$un*p.ops$disp.v.ratio*CF["mean"],
                                   p.ops$un*p.ops$disp.v.ratio*CF["median"],
                                   p.ops$un*p.ops$disp.v.ratio*CF["max"],
                                   p.ops$disp.adj*CF["min"],
                                   p.ops$disp.adj*CF["mean"],
                                   p.ops$disp.adj*CF["median"],
                                   p.ops$disp.adj*CF["max"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*CF["min"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*CF["mean"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*CF["median"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*CF["max"]))
        }
        cat("\nBalance summary across all clusters:\n")
        print.data.frame(replaceNA(round_df(c.balance.summary[, s.keep], digits)))
        
        if (!is.null(nn)) {
            for (i in rownames(x$Observations)) {
                if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
            }
            if (!is.na(x$Observations["Matched (Unweighted)",]) && all(x$Observations["Matched",] == x$Observations["Matched (Unweighted)",])) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)",]
            cat(paste0("\n", attr(x$Observations, "tag"), ":\n"))
            print.data.frame(replaceNA(x$Observations), digits = digits)
        }
    }

    invisible(x)
}
print.bal.tab.imp <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.r.threshold = "as.is", un = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", which.imp, imp.summary = "as.is", imp.fun = NULL, digits = max(3, getOption("digits") - 3), ...) {
    args <- c(as.list(environment()), list(...))[-1]
    
    call <- x$call
    i.balance <- x[["Imputation.Balance"]]
    i.balance.summary <- x[["Balance.Across.Imputations"]]
    nn <- x$Observations
    p.ops <- x$print.options
    
    #Prevent expnential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
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
        if (!identical(imp.summary, "as.is")) {
            if (!is.logical(imp.summary)) stop("imp.summary must be TRUE, FALSE, or \"as.is\"")
            p.ops$imp.summary <- imp.summary
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
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
        }
    }
    if (!missing(which.imp) && which.imp != "as.is") {
        p.ops$which.imp <- which.imp
    }
    
    imp.funs <- c("min", "mean", "median", "max")
    if (length(imp.fun) > 0) {
        imp.fun <- match.arg(tolower(imp.fun), imp.funs, several.ok = TRUE)
        if (!is.na(match("min", imp.fun)) && p.ops$quick) {
            imp.fun <- imp.fun[is.na(match(imp.fun, "min"))]
            if (length(imp.fun) == 0) {
                warning("\"min\" was the only valid entry to imp.fun, but it cannot be used when quick = TRUE in the call to bal.tab(). Using all other imp.funs instead.", call. = FALSE)
                imp.fun <- c("mean", "median", "max")
            }
            else warning("\"min\" cannot be requested when quick = TRUE in the call to bal.tab(), so it was ignored.", call. = FALSE)
        }
        else if (length(imp.fun) == 0) {
            warning("There were no valid entries to imp.fun Using all other imp.funs instead.", call. = FALSE)
            imp.fun <- c("mean", "median", "max")
        }
    }
    else imp.fun <- c("mean", "median", "max")
    
    #Checks and Adjustments
    if (length(p.ops$which.imp) == 0) 
        which.imp <- seq_along(i.balance)
    else if (is.numeric(p.ops$which.imp)) {
        which.imp <- seq_along(i.balance)[seq_along(i.balance) %in% p.ops$which.imp]
        if (length(which.imp) == 0) {
            warning("No numbers in which.imp are imputation numbers. No imputations will be displayed.", call. = FALSE)
            which.imp <- integer(0)
        }
    }
    else if (is.na(p.ops$which.imp)) {
        which.imp <- integer(0)
        p.ops$imp.summary <- TRUE
    }
    else {
        warning("The argument to which.imp must be NA, NULL, or a vector of imputation numbers. No imputations will be displayed.", call. = FALSE)
        which.imp <- integer(0)
    }
    
    #Printing output
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n")
    }
    
    if (length(which.imp) > 0) {
        cat("Balance by imputation:\n")
        for (i in which.imp) {
            cat(paste0("\n - - - Imputation: ", names(i.balance)[i], " - - - "))
            do.call(print, c(list(i.balance[[i]]), args))
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Imputation: ", names(i.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
    }
    
    if (isTRUE(as.logical(p.ops$imp.summary))) {
        IF <- !is.na(match(imp.funs, imp.fun))
        names(IF) <- imp.funs
        if (!is.na(match("bal.tab.cont", class(x)))) { #continuous
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un*IF["min"],
                                   p.ops$un*IF["mean"],
                                   p.ops$un*IF["median"],
                                   p.ops$un*IF["max"],
                                   p.ops$disp.adj*IF["min"],
                                   p.ops$disp.adj*IF["mean"],
                                   p.ops$disp.adj*IF["median"],
                                   p.ops$disp.adj*IF["max"]))
        }
        else { #binary
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un*IF["min"],
                                   p.ops$un*IF["mean"],
                                   p.ops$un*IF["median"],
                                   p.ops$un*IF["max"],
                                   p.ops$un*p.ops$disp.v.ratio*IF["min"],
                                   p.ops$un*p.ops$disp.v.ratio*IF["mean"],
                                   p.ops$un*p.ops$disp.v.ratio*IF["median"],
                                   p.ops$un*p.ops$disp.v.ratio*IF["max"],
                                   p.ops$disp.adj*IF["min"],
                                   p.ops$disp.adj*IF["mean"],
                                   p.ops$disp.adj*IF["median"],
                                   p.ops$disp.adj*IF["max"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*IF["min"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*IF["mean"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*IF["median"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*IF["max"]))
        }
        
        
        cat("\nBalance summary across all imputations:\n")
        print.data.frame(replaceNA(round_df(i.balance.summary[, s.keep], digits)))
        
        if (!is.null(nn)) {
            for (i in rownames(x$Observations)) {
                if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
            }
            if (!is.na(x$Observations["Matched (Unweighted)",]) && all(x$Observations["Matched",] == x$Observations["Matched (Unweighted)",])) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)",]
            cat(paste0("\n", attr(x$Observations, "tag"), ":\n"))
            print.data.frame(replaceNA(x$Observations), digits = digits)
        }
    }

    invisible(x)
    
}
print.bal.tab.imp.cluster <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.r.threshold = "as.is", un = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = NULL, which.imp, imp.summary = "as.is", imp.fun = NULL, digits = max(3, getOption("digits") - 3), ...) {
    args <- c(as.list(environment()), list(...))[-1]
    
    call <- x$call
    i.balance <- x[["Imputation.Balance"]]
    i.balance.c.summary <- x[["Cluster.Balance.Across.Imputations"]]
    i.balance.summary <- x[["Balance.Across.Imputations"]]
    nn <- x$Observations
    p.ops <- x$print.options
    
    #Prevent expnential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
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
        if (!identical(imp.summary, "as.is")) {
            if (!is.logical(imp.summary)) stop("imp.summary must be TRUE, FALSE, or \"as.is\"")
            p.ops$imp.summary <- imp.summary
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
    if (!missing(which.imp) && which.imp != "as.is") {
        p.ops$which.imp <- which.imp
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
    
    if (length(p.ops$which.cluster) == 0) 
        which.cluster <- seq_along(i.balance[[1]][["Cluster.Balance"]])
    else if (is.numeric(p.ops$which.cluster)) {
        which.cluster <- seq_along(i.balance[[1]][["Cluster.Balance"]])[seq_along(i.balance[[1]][["Cluster.Balance"]]) %in% p.ops$which.cluster]
        if (length(which.cluster) == 0) {
            warning("No indices in which.cluster are cluster indices. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(i.balance[[1]][["Cluster.Balance"]])
        }
    }
    else if (is.character(p.ops$which.cluster)) {
        which.cluster <- seq_along(i.balance[[1]][["Cluster.Balance"]])[names(i.balance[[1]][["Cluster.Balance"]]) %in% p.ops$which.cluster]
        if (length(which.cluster) == 0) {
            warning("No names in which.cluster are cluster names. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(i.balance[[1]][["Cluster.Balance"]])
        }
    }
    else if (is.na(p.ops$which.cluster)) {
        which.cluster <- integer(0)
        p.ops$cluster.summary <- TRUE
    }
    else {
        warning("The argument to which.cluster must be NA, NULL, or a vector of cluster indicies or cluster names. Displaying all clusters instead.", call. = FALSE)
        which.cluster <- seq_along(i.balance[[1]][["Cluster.Balance"]])
    }
    
    imp.funs <- c("min", "mean", "median", "max")
    if (length(imp.fun) > 0) {
        imp.fun <- match.arg(tolower(imp.fun), imp.funs, several.ok = TRUE)
        if (!is.na(match("min", imp.fun)) && p.ops$quick) {
            imp.fun <- imp.fun[is.na(match(imp.fun, "min"))]
            if (length(imp.fun) == 0) {
                warning("\"min\" was the only valid entry to imp.fun, but it cannot be used when quick = TRUE in the call to bal.tab(). Using all other imp.funs instead.", call. = FALSE)
                imp.fun <- c("mean", "median", "max")
            }
            else warning("\"min\" cannot be requested when quick = TRUE in the call to bal.tab(), so it was ignored.", call. = FALSE)
        }
        else if (length(imp.fun) == 0) {
            warning("There were no valid entries to imp.fun Using all other imp.funs instead.", call. = FALSE)
            imp.fun <- c("mean", "median", "max")
        }
    }
    else imp.fun <- c("mean", "median", "max")
    
    if (length(p.ops$which.imp) == 0) 
        which.imp <- seq_along(i.balance)
    else if (is.numeric(p.ops$which.imp)) {
        which.imp <- seq_along(i.balance)[seq_along(i.balance) %in% p.ops$which.imp]
        if (length(which.imp) == 0) {
            warning("No numbers in which.imp are imputation numbers. No imputations will be displayed.", call. = FALSE)
            which.imp <- integer(0)
        }
    }
    else if (is.na(p.ops$which.imp)) {
        which.imp <- integer(0)
        p.ops$imp.summary <- TRUE
    }
    else {
        warning("The argument to which.imp must be NA, NULL, or a vector of imputation numbers. No imputations will be displayed.", call. = FALSE)
        which.imp <- integer(0)
    }
    
    #Printing output
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n")
    }
    
    if (length(which.imp) > 0) {
        cat("Balance by imputation:\n")
        for (i in which.imp) {
            cat(paste0("\n - - - - Imputation: ", names(i.balance)[i], " - - - - \n"))
            do.call(print, c(list(i.balance[[i]]), args))
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - - Imputation: ", names(i.balance)[i], " - - - - "))/2)), collapse = ""), " \n"))
    }
    
    if (isTRUE(as.logical(p.ops$imp.summary))) {
        IF <- !is.na(match(imp.funs, imp.fun))
        names(IF) <- imp.funs
        if (!is.na(match("bal.tab.cont", class(x)))) {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un*IF["min"],
                                   p.ops$un*IF["mean"],
                                   p.ops$un*IF["median"],
                                   p.ops$un*IF["max"],
                                   p.ops$disp.adj*IF["min"],
                                   p.ops$disp.adj*IF["mean"],
                                   p.ops$disp.adj*IF["median"],
                                   p.ops$disp.adj*IF["max"]))
        }
        else {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un*IF["min"],
                                   p.ops$un*IF["mean"],
                                   p.ops$un*IF["median"],
                                   p.ops$un*IF["max"],
                                   p.ops$un*p.ops$disp.v.ratio*IF["min"],
                                   p.ops$un*p.ops$disp.v.ratio*IF["mean"],
                                   p.ops$un*p.ops$disp.v.ratio*IF["median"],
                                   p.ops$un*p.ops$disp.v.ratio*IF["max"],
                                   p.ops$disp.adj*IF["min"],
                                   p.ops$disp.adj*IF["mean"],
                                   p.ops$disp.adj*IF["median"],
                                   p.ops$disp.adj*IF["max"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*IF["min"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*IF["mean"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*IF["median"],
                                   p.ops$disp.adj*p.ops$disp.v.ratio*IF["max"]))
        }
        
        
        if (length(which.cluster) > 0) {
            cat("\nCluster balance summary across all imputations:\n")
            for (c in which.cluster) {
                cat(paste0("\n - - - Cluster: ", names(i.balance.c.summary)[c], " - - - \n"))
                cat("Balance summary across imputations:\n")
                print.data.frame(replaceNA(round_df(i.balance.c.summary[[c]][["Cluster.Balance"]][, s.keep], digits)))
                
                for (j in rownames(i.balance.c.summary[[c]][["Cluster.Observations"]])) {
                    if (all(i.balance.c.summary[[c]][["Cluster.Observations"]][j,] == 0)) i.balance.c.summary[[c]][["Cluster.Observations"]] <- i.balance.c.summary[[c]][["Cluster.Observations"]][rownames(i.balance.c.summary[[c]][["Cluster.Observations"]])!=j,]
                }
                if (!is.na(i.balance.c.summary[[c]][["Cluster.Observations"]]["Matched (Unweighted)",]) && all(i.balance.c.summary[[c]][["Cluster.Observations"]]["Matched",] == i.balance.c.summary[[c]][["Cluster.Observations"]]["Matched (Unweighted)",])) i.balance.c.summary[[c]][["Cluster.Observations"]] <- i.balance.c.summary[[c]][["Cluster.Observations"]][rownames(i.balance.c.summary[[c]][["Cluster.Observations"]])!="Matched (Unweighted)",]
                cat(paste0("\n", attr(i.balance.c.summary[[c]][["Cluster.Observations"]], "tag"), ":\n"))
                print.data.frame(replaceNA(round_df(i.balance.c.summary[[c]][["Cluster.Observations"]], digits)))
            }
        }
        cat("\nBalance summary across all imputations and clusters:\n")
        print.data.frame(replaceNA(round_df(i.balance.summary[, s.keep], digits)))
        
        if (!is.null(nn)) {
            for (i in rownames(x$Observations)) {
                if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
            }
            if (!is.na(x$Observations["Matched (Unweighted)",]) && all(x$Observations["Matched",] == x$Observations["Matched (Unweighted)",])) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)",]
            cat(paste0("\n", attr(nn, "tag"), ":\n"))
            print.data.frame(replaceNA(nn), digits = digits)
        }
    }
    invisible(x)
    
}