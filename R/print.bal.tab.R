print.bal.tab <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    balance <- x$Balance
    baltal.r <- x$Balanced.Corr
    maximbal.r <- x$Max.Imbalance.Corr
    baltal.m <- x$Balanced.Means
    maximbal.m <- x$Max.Imbalance.Means
    baltal.v <- x$Balanced.Variances
    maximbal.v <- x$Max.Imbalance.Variances
    baltal.ks <- x$Balanced.KS
    maximbal.ks <- x$Max.Imbalance.KS
    nn <- x$Observations
    p.ops <- x$print.options
  
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("un cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp.means, "as.is")) {
        if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.means == FALSE && disp.means == TRUE) {
            warning("disp.means cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.means <- disp.means
    }
    if (!identical(disp.v.ratio, "as.is")) {
        if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.v.ratio == FALSE && disp.v.ratio == TRUE) {
            warning("disp.v.ratio cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.v.ratio <- disp.v.ratio
    }
    if (!identical(disp.ks, "as.is")) {
        if (!is.logical(disp.ks)) stop("disp.ks must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.ks == FALSE && disp.ks == TRUE) {
            warning("disp.ks cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.ks <- disp.ks
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
    if (length(p.ops$disp.v.ratio) == 0 || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (length(p.ops$disp.ks) == 0 || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (!any(sapply(c(p.ops$m.threshold, 
                              p.ops$v.threshold, 
                              p.ops$ks.threshold, 
                              p.ops$r.threshold), length) > 0)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (!is.na(match("bal.tab.cont", class(x)))) {
        keep <- as.logical(c(TRUE, 
                             p.ops$un, 
                             rep(c(p.ops$disp.adj, 
                                   !is.null(p.ops$r.threshold)), p.ops$nweights + !p.ops$disp.adj)))
        
    }
    else {
        keep <- setNames(as.logical(c(TRUE, 
                                      p.ops$un && p.ops$disp.means, 
                                      p.ops$un && p.ops$disp.means, 
                                      p.ops$un, 
                                      p.ops$un && !p.ops$disp.adj && !is.null(p.ops$m.threshold),
                                      p.ops$un && p.ops$disp.v.ratio, 
                                      p.ops$un && !p.ops$disp.adj && !is.null(p.ops$v.threshold), 
                                      p.ops$un && p.ops$disp.ks, 
                                      p.ops$un && !p.ops$disp.adj && !is.null(p.ops$ks.threshold),
                                      rep(c(p.ops$disp.adj && p.ops$disp.means, 
                                            p.ops$disp.adj && p.ops$disp.means, 
                                            p.ops$disp.adj, 
                                            p.ops$disp.adj && !is.null(p.ops$m.threshold), 
                                            p.ops$disp.adj && p.ops$disp.v.ratio, 
                                            p.ops$disp.adj && !is.null(p.ops$v.threshold), 
                                            p.ops$disp.adj && p.ops$disp.ks, 
                                            p.ops$disp.adj && !is.null(p.ops$ks.threshold)), p.ops$nweights + !p.ops$disp.adj))),
                         names(balance))
    }
    
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n  ")
        cat("\n")
    }
    
    if (p.ops$disp.bal.tab) {
        if (p.ops$imbalanced.only) {
            keep.row <- rowSums(apply(balance[grepl(".Threshold", names(balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        }
        else keep.row <- rep(TRUE, nrow(balance))
        
        cat("Balance Measures:\n")
        print.data.frame(round_df_char(balance[keep.row, keep], digits))
        cat("\n")
    }
    
    if (!is.null(baltal.r)) {
        cat("Balance tally for correlations:\n")
        print.data.frame(x$Balanced.Corr)
        cat("\n")
    }
    if (!is.null(maximbal.r)) {
        cat("Variable with the greatest treatment correlation:\n")
        print.data.frame(round_df(x$Max.Imbalance.Corr, digits), row.names = FALSE)
        cat("\n")
    }
    if (!is.null(baltal.m)) {
        cat("Balance tally for mean differences:\n")
        print.data.frame(x$Balanced.Means)
        cat("\n")
    }
    if (!is.null(maximbal.m)) {
        cat("Variable with the greatest mean difference:\n")
        print.data.frame(round_df(x$Max.Imbalance.Means, digits), row.names = FALSE)
        cat("\n")
    }
    if (!is.null(baltal.v)) {
        cat("Balance tally for variance ratios:\n")
        print.data.frame(x$Balanced.Variances, digits)
        cat("\n")
    }
    if (!is.null(maximbal.v)) {
        cat("Variable with the greatest variance ratio:\n")
        print.data.frame(round_df(x$Max.Imbalance.Variances, digits), row.names = FALSE)
        cat("\n")
    }
    if (!is.null(baltal.ks)) {
        cat("Balance tally for KS statistics:\n")
        print.data.frame(x$Balanced.KS, digits)
        cat("\n")
    }
    if (!is.null(maximbal.ks)) {
        cat("Variable with the greatest KS statistic:\n")
        print.data.frame(round_df(x$Max.Imbalance.KS, digits), row.names = FALSE)
        cat("\n")
    }
    if (!is.null(nn)) {
        for (i in rownames(x$Observations)) {
            if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
        }
        if (!is.na(x$Observations["Matched (Unweighted)",]) && all(x$Observations["Matched",] == x$Observations["Matched (Unweighted)",])) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)",]
        cat(paste0(attr(x$Observations, "tag"), ":\n"))
        print.warning <- FALSE
        if (length(attr(x$Observations, "ss.type")) > 1 && length(unique(attr(x$Observations, "ss.type")[-1])) > 1) {
            ess <- ifelse(attr(x$Observations, "ss.type") == "ess", "*", "")
            x$Observations <- setNames(cbind(x$Observations, ess), c(names(x$Observations), ""))
            print.warning <- TRUE
        }
        print.data.frame(replaceNA(x$Observations), digits = digits)
        if (print.warning) cat("* indicates effective sample size")
    }
    invisible(x)
}
print.bal.tab.subclass <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", disp.subclass = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    s.balance <- x$Subclass.Balance
    b.a.subclass <- x$Balance.Across.Subclass
    baltal.r.subclass <- x$Balanced.Corr.Subclass
    maximbal.r.subclass <- x$Max.Imbalance.Corr.Subclass
    baltal.m.subclass <- x$Balanced.Means.Subclass
    maximbal.m.subclass <- x$Max.Imbalance.Means.Subclass
    baltal.v.subclass <- x$Balanced.Variances.Subclass
    maximbal.v.subclass <- x$Max.Imbalance.Variances.Subclass
    baltal.ks.subclass <- x$Balanced.KS.Subclass
    maximbal.ks.subclass <- x$Max.Imbalance.KS.Subclass
    s.nn <- x$Subclass.Observations
    p.ops <- x$print.options
    
    #Prevent expnential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("un cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp.means, "as.is")) {
        if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.means == FALSE && disp.means == TRUE) {
            warning("disp.means cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.means <- disp.means
    }
    if (!identical(disp.v.ratio, "as.is")) {
        if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.v.ratio == FALSE && disp.v.ratio == TRUE) {
            warning("disp.v.ratio cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.v.ratio <- disp.v.ratio
    }
    if (!identical(disp.ks, "as.is")) {
        if (!is.logical(disp.ks)) stop("disp.ks must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.ks == FALSE && disp.ks == TRUE) {
            warning("disp.ks cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.ks <- disp.ks
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
            baltal.r.subclass <- NULL
            maximbal.r.subclass <- NULL
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
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (length(p.ops$disp.v.ratio) == 0 || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (length(p.ops$disp.ks) == 0 || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (!any(sapply(c(p.ops$m.threshold, 
                              p.ops$v.threshold, 
                              p.ops$ks.threshold, 
                              p.ops$r.threshold), length) > 0)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (!identical(disp.subclass, "as.is")) {
        if (!is.logical(disp.subclass)) stop("disp.subclass must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.subclass <- disp.subclass
    }
    
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n  ")
        cat("\n")
    }
    if (p.ops$disp.bal.tab) {
        if (p.ops$disp.subclass) {
            if (!is.na(match("bal.tab.cont", class(x)))) {
                s.keep <- as.logical(c(TRUE, 
                                       p.ops$disp.adj, 
                                       !is.null(p.ops$r.threshold)))
            }
            else {
                s.keep <- as.logical(c(TRUE, 
                                       p.ops$disp.means, 
                                       p.ops$disp.means, 
                                       p.ops$disp.adj, 
                                       !is.null(p.ops$m.threshold), 
                                       p.ops$disp.adj  &&  p.ops$disp.v.ratio, 
                                       !is.null(p.ops$v.threshold), 
                                       p.ops$disp.adj  &&  p.ops$disp.ks, 
                                       !is.null(p.ops$ks.threshold)))
            }
            cat("Balance by subclass:")
            for (i in names(s.balance)) {
                if (p.ops$imbalanced.only) {
                    keep.row <- rowSums(apply(s.balance[[i]][grepl(".Threshold", names(s.balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
                }
                else keep.row <- rep(TRUE, nrow(s.balance[[i]]))
                cat(paste0("\n - - - Subclass ", i, " - - - \n"))
                print.data.frame(replaceNA(round_df(s.balance[[i]][keep.row, s.keep, drop = FALSE], digits)))
            }
            cat("\n")
        }
        
        if (!is.null(b.a.subclass)) {
            if (p.ops$imbalanced.only) {
                keep.row <- rowSums(apply(b.a.subclass[grepl(".Threshold", names(b.a.subclass), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
            }
            else keep.row <- rep(TRUE, nrow(b.a.subclass))
            a.s.keep <- as.logical(c(TRUE, 
                                     p.ops$un && p.ops$disp.means, 
                                     p.ops$un && p.ops$disp.means, 
                                     p.ops$un, 
                                     p.ops$disp.adj && p.ops$disp.means, 
                                     p.ops$disp.adj && p.ops$disp.means, 
                                     p.ops$disp.adj, 
                                     !is.null(p.ops$m.threshold)))
            cat("Balance measures across subclasses:\n")
            print.data.frame(replaceNA(round_df(b.a.subclass[keep.row, a.s.keep, drop = FALSE], digits)))
            cat("\n")
        }
    }
    if (!is.null(baltal.r.subclass)) {
        cat("Balance tally for correlations across subclasses:\n")
        print.data.frame(baltal.r.subclass)
        cat("\n")
    }
    if (!is.null(maximbal.r.subclass)) {
        cat("Variable with the greatest treatment correlation across subclasses:\n")
        print.data.frame(round_df(maximbal.r.subclass, digits), row.names = TRUE)
        cat("\n")
    }
    if (!is.null(baltal.m.subclass)) {
        cat("Balance tally for mean differences across subclasses:\n")
        print.data.frame(baltal.m.subclass)
        cat("\n")
    }
    if (!is.null(maximbal.m.subclass)) {
        cat("Variable with the greatest mean difference across subclasses:\n")
        print.data.frame(round_df(maximbal.m.subclass, digits), row.names = TRUE)
        cat("\n")
    }
    if (!is.null(baltal.v.subclass)) {
        cat("Balance tally for variance ratios across subclasses:\n")
        print.data.frame(baltal.v.subclass)
        cat("\n")
    }
    if (!is.null(maximbal.v.subclass)) {
        cat("Variable with the greatest variance ratios across subclasses:\n")
        print.data.frame(round_df(maximbal.v.subclass, digits), row.names = TRUE)
        cat("\n")
    }
    if (!is.null(baltal.ks.subclass)) {
        cat("Balance tally for KS statistics across subclasses:\n")
        print.data.frame(baltal.ks.subclass)
        cat("\n")
    }
    if (!is.null(maximbal.ks.subclass)) {
        cat("Variable with the greatest KS statistc across subclasses:\n")
        print.data.frame(round_df(maximbal.ks.subclass, digits), row.names = TRUE)
        cat("\n")
    }
    
    if (!is.null(s.nn)) {
        cat(paste0(attr(x$Subclass.Observations, "tag"), ":\n"))
        print.data.frame(replaceNA(x$Subclass.Observations), digits = digits)
    }
    
    invisible(x)
}
print.bal.tab.cluster <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = NULL, digits = max(3, getOption("digits") - 3), ...) {
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
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("un cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp.means, "as.is")) {
        if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.means == FALSE && disp.means == TRUE) {
            warning("disp.means cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.means <- disp.means
    }
    if (!identical(disp.v.ratio, "as.is")) {
        if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.v.ratio == FALSE && disp.v.ratio == TRUE) {
            warning("disp.v.ratio cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.v.ratio <- disp.v.ratio
    }
    if (!identical(disp.ks, "as.is")) {
        if (!is.logical(disp.ks)) stop("disp.ks must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.ks == FALSE && disp.ks == TRUE) {
            warning("disp.ks cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.ks <- disp.ks
    }
    if (!identical(cluster.summary, "as.is")) {
        if (!is.logical(cluster.summary)) stop("cluster.summary must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$cluster.summary == FALSE && cluster.summary == TRUE) {
            warning("cluster.summary cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$cluster.summary <- cluster.summary
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
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (length(p.ops$disp.v.ratio) == 0 || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (length(p.ops$disp.ks) == 0 || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
        }
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (!any(sapply(c(p.ops$m.threshold, 
                              p.ops$v.threshold, 
                              p.ops$ks.threshold, 
                              p.ops$r.threshold), length) > 0)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
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
                             rep(c(p.ops$disp.adj, 
                                   !is.null(p.ops$r.threshold)), p.ops$nweights)))
    }
    else {
        keep <- as.logical(c(TRUE, 
                             p.ops$un && p.ops$disp.means, 
                             p.ops$un && p.ops$disp.means, 
                             p.ops$un, 
                             p.ops$un && !p.ops$disp.adj && !is.null(p.ops$m.threshold),
                             p.ops$un && p.ops$disp.v.ratio, 
                             p.ops$un && !p.ops$disp.adj && !is.null(p.ops$v.threshold), 
                             p.ops$un && p.ops$disp.ks, 
                             p.ops$un && !p.ops$disp.adj && !is.null(p.ops$ks.threshold), 
                             rep(c(p.ops$disp.adj && p.ops$disp.means, 
                                   p.ops$disp.adj && p.ops$disp.means, 
                                   p.ops$disp.adj, 
                                   p.ops$disp.adj && !is.null(p.ops$m.threshold), 
                                   p.ops$disp.adj && p.ops$disp.v.ratio, 
                                   p.ops$disp.adj && !is.null(p.ops$v.threshold),
                                   p.ops$disp.adj && p.ops$disp.ks, 
                                   p.ops$disp.adj && !is.null(p.ops$ks.threshold)), 
                                 p.ops$nweights + !p.ops$disp.adj)))
    }
    
    #Printing
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n  ")
        cat("\n")
    }
    
    if (length(which.cluster)>0) {
        cat("Balance by cluster:\n")
        for (i in which.cluster) {
            if (p.ops$imbalanced.only) {
                keep.row <- rowSums(apply(c.balance[[i]][["Balance"]][grepl(".Threshold", names(c.balance[[i]][["Balance"]]), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
            }
            else keep.row <- rep(TRUE, nrow(c.balance[[i]][["Balance"]]))
            cat(paste0("\n - - - Cluster: ", names(c.balance)[i], " - - - \n"))
            if (p.ops$disp.bal.tab) {
                cat("Balance measures:\n")
                print.data.frame(round_df_char(c.balance[[i]][["Balance"]][keep.row, keep], digits))
            }
            for (j in rownames(c.balance[[i]][["Observations"]])) {
                if (all(c.balance[[i]][["Observations"]][j,] == 0)) c.balance[[i]][["Observations"]] <- c.balance[[i]][["Observations"]][rownames(c.balance[[i]][["Observations"]])!=j,]
            }
            if (!is.na(c.balance[[i]][["Observations"]]["Matched (Unweighted)",]) && all(c.balance[[i]][["Observations"]]["Matched",] == c.balance[[i]][["Observations"]]["Matched (Unweighted)",])) c.balance[[i]][["Observations"]] <- c.balance[[i]][["Observations"]][rownames(c.balance[[i]][["Observations"]])!="Matched (Unweighted)",]
            cat(paste0("\n", attr(c.balance[[i]][["Observations"]], "tag"), ":\n"))
            print.warning <- FALSE
            if (length(attr(c.balance[[i]][["Observations"]], "ss.type")) > 1 && length(unique(attr(c.balance[[i]][["Observations"]], "ss.type")[-1])) > 1) {
                ess <- ifelse(attr(c.balance[[i]][["Observations"]], "ss.type") == "ess", "*", "")
                xc.balance[[i]][["Observations"]] <- setNames(cbind(c.balance[[i]][["Observations"]], ess), c(names(c.balance[[i]][["Observations"]]), ""))
                print.warning <- TRUE
            }
            print.data.frame(round_df_char(c.balance[[i]][["Observations"]], digits))
            if (print.warning) cat("* indicates effective sample size")
            
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Cluster: ", names(c.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$cluster.summary))) {
        CF <- !is.na(match(cluster.funs, cluster.fun))
        names(CF) <- cluster.funs
        if (!is.na(match("bal.tab.cont.cluster", class(x)))) {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && CF["min"],
                                   p.ops$un && CF["mean"],
                                   p.ops$un && CF["median"],
                                   p.ops$un && CF["max"],
                                   rep(c(p.ops$disp.adj && CF["min"],
                                         p.ops$disp.adj && CF["mean"],
                                         p.ops$disp.adj && CF["median"],
                                         p.ops$disp.adj && CF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        else {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && CF["min"],
                                   p.ops$un && CF["mean"],
                                   p.ops$un && CF["median"],
                                   p.ops$un && CF["max"],
                                   p.ops$un && p.ops$disp.v.ratio && CF["min"],
                                   p.ops$un && p.ops$disp.v.ratio && CF["mean"],
                                   p.ops$un && p.ops$disp.v.ratio && CF["median"],
                                   p.ops$un && p.ops$disp.v.ratio && CF["max"],
                                   p.ops$un && p.ops$disp.ks && CF["min"],
                                   p.ops$un && p.ops$disp.ks && CF["mean"],
                                   p.ops$un && p.ops$disp.ks && CF["median"],
                                   p.ops$un && p.ops$disp.ks && CF["max"],
                                   rep(c(p.ops$disp.adj && CF["min"],
                                         p.ops$disp.adj && CF["mean"],
                                         p.ops$disp.adj && CF["median"],
                                         p.ops$disp.adj && CF["max"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && CF["min"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && CF["mean"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && CF["median"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && CF["max"],
                                         p.ops$disp.adj && p.ops$disp.ks && CF["min"],
                                         p.ops$disp.adj && p.ops$disp.ks && CF["mean"],
                                         p.ops$disp.adj && p.ops$disp.ks && CF["median"],
                                         p.ops$disp.adj && p.ops$disp.ks && CF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        if (p.ops$disp.bal.tab) {
            cat("Balance summary across all clusters:\n")
            print.data.frame(round_df_char(c.balance.summary[, s.keep], digits))
            cat("\n")
        }
        
        if (!is.null(nn)) {
            for (i in rownames(x$Observations)) {
                if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
            }
            if (!is.na(x$Observations["Matched (Unweighted)",]) && all(x$Observations["Matched",] == x$Observations["Matched (Unweighted)",])) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)",]
            cat(paste0(attr(x$Observations, "tag"), ":\n"))
            print.warning <- FALSE
            if (length(attr(x$Observations, "ss.type")) > 1 && length(unique(attr(x$Observations, "ss.type")[-1])) > 1) {
                ess <- ifelse(attr(x$Observations, "ss.type") == "ess", "*", "")
                x$Observations <- setNames(cbind(x$Observations, ess), c(names(x$Observations), ""))
                print.warning <- TRUE
            }
            print.data.frame(replaceNA(x$Observations), digits = digits)
            if (print.warning) cat("* indicates effective sample size")
        }
    }
    
    invisible(x)
}
print.bal.tab.imp <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.imp, imp.summary = "as.is", imp.fun = NULL, digits = max(3, getOption("digits") - 3), ...) {
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
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("un cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp.means, "as.is")) {
        if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.means == FALSE && disp.means == TRUE) {
            warning("disp.means cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.means <- disp.means
    }
    if (!identical(disp.v.ratio, "as.is")) {
        if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.v.ratio == FALSE && disp.v.ratio == TRUE) {
            warning("disp.v.ratio cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.v.ratio <- disp.v.ratio
    }
    if (!identical(disp.ks, "as.is")) {
        if (!is.logical(disp.ks)) stop("disp.ks must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.ks == FALSE && disp.ks == TRUE) {
            warning("disp.ks cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.ks <- disp.ks
    }
    if (!identical(imp.summary, "as.is")) {
        if (!is.logical(imp.summary)) stop("imp.summary must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$imp.summary == FALSE && imp.summary == TRUE) {
            warning("imp.summary cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$imp.summary <- imp.summary
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
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (length(p.ops$disp.v.ratio) == 0 || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (length(p.ops$disp.ks) == 0 || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
        }
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (!any(sapply(c(p.ops$m.threshold, 
                              p.ops$v.threshold, 
                              p.ops$ks.threshold, 
                              p.ops$r.threshold), length) > 0)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    # if (p.ops$imbalanced.only) {
    #     keep.row <- rowSums(apply(balance[grepl(".Threshold", names(balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
    # }
    # else keep.row <- rep(TRUE, nrow(balance))
    
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
        p.ops$imp.summary <- TRUE
    }
    
    #Printing output
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n  ")
        cat("\n")
    }
    
    if (length(which.imp) > 0) {
        cat("Balance by imputation:\n")
        for (i in which.imp) {
            cat(paste0("\n - - - Imputation: ", names(i.balance)[i], " - - - \n"))
            do.call(print, c(list(i.balance[[i]]), args))
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Imputation: ", names(i.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$imp.summary))) {
        IF <- !is.na(match(imp.funs, imp.fun))
        names(IF) <- imp.funs
        if (!is.na(match("bal.tab.cont", class(x)))) { #continuous
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && IF["min"],
                                   p.ops$un && IF["mean"],
                                   p.ops$un && IF["median"],
                                   p.ops$un && IF["max"],
                                   rep(c(p.ops$disp.adj && IF["min"],
                                         p.ops$disp.adj && IF["mean"],
                                         p.ops$disp.adj && IF["median"],
                                         p.ops$disp.adj && IF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        else { #binary
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && IF["min"],
                                   p.ops$un && IF["mean"],
                                   p.ops$un && IF["median"],
                                   p.ops$un && IF["max"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["min"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["mean"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["median"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["max"],
                                   p.ops$un && p.ops$disp.ks && IF["min"],
                                   p.ops$un && p.ops$disp.ks && IF["mean"],
                                   p.ops$un && p.ops$disp.ks && IF["median"],
                                   p.ops$un && p.ops$disp.ks && IF["max"],
                                   rep(c(p.ops$disp.adj && IF["min"],
                                         p.ops$disp.adj && IF["mean"],
                                         p.ops$disp.adj && IF["median"],
                                         p.ops$disp.adj && IF["max"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["min"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["mean"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["median"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["max"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["min"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["mean"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["median"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        
        if (p.ops$disp.bal.tab) {
            cat("Balance summary across all imputations:\n")
            print.data.frame(round_df_char(i.balance.summary[, s.keep], digits))
            cat("\n")
        }
        
        if (!is.null(nn)) {
            for (i in rownames(x$Observations)) {
                if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
            }
            if (!is.na(x$Observations["Matched (Unweighted)",]) && all(x$Observations["Matched",] == x$Observations["Matched (Unweighted)",])) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)",]
            cat(paste0(attr(x$Observations, "tag"), ":\n"))
            print.warning <- FALSE
            if (length(attr(x$Observations, "ss.type")) > 1 && length(unique(attr(x$Observations, "ss.type")[-1])) > 1) {
                ess <- ifelse(attr(x$Observations, "ss.type") == "ess", "*", "")
                x$Observations <- setNames(cbind(x$Observations, ess), c(names(x$Observations), ""))
                print.warning <- TRUE
            }
            print.data.frame(replaceNA(x$Observations), digits = digits)
            if (print.warning) cat("* indicates effective sample size")
        }
    }
    
    invisible(x)
    
}
print.bal.tab.imp.cluster <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = NULL, which.imp, imp.summary = "as.is", imp.fun = NULL, digits = max(3, getOption("digits") - 3), ...) {
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
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("un cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp.means, "as.is")) {
        if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.means == FALSE && disp.means == TRUE) {
            warning("disp.means cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.means <- disp.means
    }
    if (!identical(disp.v.ratio, "as.is")) {
        if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.v.ratio == FALSE && disp.v.ratio == TRUE) {
            warning("disp.v.ratio cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.v.ratio <- disp.v.ratio
    }
    if (!identical(disp.ks, "as.is")) {
        if (!is.logical(disp.ks)) stop("disp.ks must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.ks == FALSE && disp.ks == TRUE) {
            warning("disp.ks cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.ks <- disp.ks
    }
    if (!identical(cluster.summary, "as.is")) {
        if (!is.logical(cluster.summary)) stop("cluster.summary must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$cluster.summary == FALSE && cluster.summary == TRUE) {
            warning("cluster.summary cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$cluster.summary <- cluster.summary
    }
    if (!identical(imp.summary, "as.is")) {
        if (!is.logical(imp.summary)) stop("imp.summary must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$imp.summary == FALSE && imp.summary == TRUE) {
            warning("imp.summary cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$imp.summary <- imp.summary
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
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (length(p.ops$disp.v.ratio) == 0 || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (length(p.ops$disp.ks) == 0 || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (!any(sapply(c(p.ops$m.threshold, 
                              p.ops$v.threshold, 
                              p.ops$ks.threshold, 
                              p.ops$r.threshold), length) > 0)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
     
    # if (p.ops$imbalanced.only) {
    #     keep.row <- rowSums(apply(balance[grepl(".Threshold", names(balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
    # }
    # else keep.row <- rep(TRUE, nrow(balance))
    
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
        cat("\nCall:", deparse(call), sep = "\n  ")
        cat("\n")
    }
    
    if (length(which.imp) > 0) {
        cat("Balance by imputation:\n")
        for (i in which.imp) {
            cat(paste0("\n - - - - Imputation: ", names(i.balance)[i], " - - - - \n"))
            do.call(print, c(list(i.balance[[i]]), args))
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - - Imputation: ", names(i.balance)[i], " - - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$imp.summary))) {
        IF <- !is.na(match(imp.funs, imp.fun))
        names(IF) <- imp.funs
        if (!is.na(match("bal.tab.cont", class(x)))) {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && IF["min"],
                                   p.ops$un && IF["mean"],
                                   p.ops$un && IF["median"],
                                   p.ops$un && IF["max"],
                                   rep(c(p.ops$disp.adj && IF["min"],
                                         p.ops$disp.adj && IF["mean"],
                                         p.ops$disp.adj && IF["median"],
                                         p.ops$disp.adj && IF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        else {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && IF["min"],
                                   p.ops$un && IF["mean"],
                                   p.ops$un && IF["median"],
                                   p.ops$un && IF["max"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["min"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["mean"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["median"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["max"],
                                   p.ops$un && p.ops$disp.ks && IF["min"],
                                   p.ops$un && p.ops$disp.ks && IF["mean"],
                                   p.ops$un && p.ops$disp.ks && IF["median"],
                                   p.ops$un && p.ops$disp.ks && IF["max"],
                                   rep(c(p.ops$disp.adj && IF["min"],
                                         p.ops$disp.adj && IF["mean"],
                                         p.ops$disp.adj && IF["median"],
                                         p.ops$disp.adj && IF["max"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["min"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["mean"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["median"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["max"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["min"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["mean"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["median"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        
        
        if (length(which.cluster) > 0) {
            cat("Cluster balance summary across all imputations:\n")
            for (c in which.cluster) {
                cat(paste0("\n - - - Cluster: ", names(i.balance.c.summary)[c], " - - - \n"))
                if (p.ops$disp.bal.tab) {
                    cat("Balance summary across imputations:\n")
                    print.data.frame(round_df_char(i.balance.c.summary[[c]][["Cluster.Balance"]][, s.keep], digits))
                }
                for (j in rownames(i.balance.c.summary[[c]][["Cluster.Observations"]])) {
                    if (all(i.balance.c.summary[[c]][["Cluster.Observations"]][j,] == 0)) i.balance.c.summary[[c]][["Cluster.Observations"]] <- i.balance.c.summary[[c]][["Cluster.Observations"]][rownames(i.balance.c.summary[[c]][["Cluster.Observations"]])!=j,]
                }
                if (!is.na(i.balance.c.summary[[c]][["Cluster.Observations"]]["Matched (Unweighted)",]) && all(i.balance.c.summary[[c]][["Cluster.Observations"]]["Matched",] == i.balance.c.summary[[c]][["Cluster.Observations"]]["Matched (Unweighted)",])) i.balance.c.summary[[c]][["Cluster.Observations"]] <- i.balance.c.summary[[c]][["Cluster.Observations"]][rownames(i.balance.c.summary[[c]][["Cluster.Observations"]])!="Matched (Unweighted)",]
                cat(paste0("\n", attr(i.balance.c.summary[[c]][["Cluster.Observations"]], "tag"), ":\n"))
                print.warning <- FALSE
                if (length(attr(i.balance.c.summary[[c]][["Cluster.Observations"]], "ss.type")) > 1 && length(unique(attr(i.balance.c.summary[[c]][["Cluster.Observations"]], "ss.type")[-1])) > 1) {
                    ess <- ifelse(attr(i.balance.c.summary[[c]][["Cluster.Observations"]], "ss.type") == "ess", "*", "")
                    i.balance.c.summary[[c]][["Cluster.Observations"]] <- setNames(cbind(i.balance.c.summary[[c]][["Cluster.Observations"]], ess), c(names(i.balance.c.summary[[c]][["Cluster.Observations"]]), ""))
                    print.warning <- TRUE
                }
                print.data.frame(round_df_char(i.balance.c.summary[[c]][["Cluster.Observations"]], digits))
                if (print.warning) cat("* indicates effective sample size")
            }
            cat("\n")
        }
        if (p.ops$disp.bal.tab) {
        cat("Balance summary across all imputations and clusters:\n")
        print.data.frame(round_df_char(i.balance.summary[, s.keep], digits))
        cat("\n")
        }
        if (!is.null(nn)) {
            for (i in rownames(x$Observations)) {
                if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
            }
            if (!is.na(x$Observations["Matched (Unweighted)",]) && all(x$Observations["Matched",] == x$Observations["Matched (Unweighted)",])) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)",]
            cat(paste0(attr(x$Observations, "tag"), ":\n"))
            print.warning <- FALSE
            if (length(attr(x$Observations, "ss.type")) > 1 && length(unique(attr(x$Observations, "ss.type")[-1])) > 1) {
                ess <- ifelse(attr(x$Observations, "ss.type") == "ess", "*", "")
                x$Observations <- setNames(cbind(x$Observations, ess), c(names(x$Observations), ""))
                print.warning <- TRUE
            }
            print.data.frame(replaceNA(x$Observations), digits = digits)
            if (print.warning) cat("* indicates effective sample size")
        }
    }
    invisible(x)
    
}
print.bal.tab.multi <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.treat, multi.summary = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    call <- x$call
    m.balance <- x[["Pair.Balance"]]
    m.balance.summary <- x[["Balance.Across.Pairs"]]
    nn <- x$Observations
    p.ops <- x$print.options
    
    #Prevent expnential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("un cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp.means, "as.is")) {
        if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.means == FALSE && disp.means == TRUE) {
            warning("disp.means cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.means <- disp.means
    }
    if (!identical(disp.v.ratio, "as.is")) {
        if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.v.ratio == FALSE && disp.v.ratio == TRUE) {
            warning("disp.v.ratio cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.v.ratio <- disp.v.ratio
    }
    if (!identical(disp.ks, "as.is")) {
        if (!is.logical(disp.ks)) stop("disp.ks must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.ks == FALSE && disp.ks == TRUE) {
            warning("disp.ks cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.ks <- disp.ks
    }
    if (!identical(multi.summary, "as.is")) {
        if (!is.logical(multi.summary)) stop("multi.summary must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$multi.summary == FALSE && multi.summary == TRUE) {
            warning("multi.summary cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$multi.summary <- multi.summary
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
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (length(p.ops$disp.v.ratio) == 0 || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (length(p.ops$disp.ks) == 0 || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (!any(sapply(c(p.ops$m.threshold, 
                              p.ops$v.threshold, 
                              p.ops$ks.threshold, 
                              p.ops$r.threshold), length) > 0)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (p.ops$imbalanced.only) {
        keep.row <- rowSums(apply(m.balance.summary[grepl(".Threshold", names(m.balance.summary), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
    }
    else keep.row <- rep(TRUE, nrow(m.balance.summary))
    
    if (!missing(which.treat)) {
        p.ops$which.treat <- which.treat
    }
    
    #Checks and Adjustments
    if (length(p.ops$which.treat) == 0) 
        which.treat <- p.ops$treat.names
    else if (is.numeric(p.ops$which.treat)) {
        which.treat <- p.ops$treat.names[seq_along(p.ops$treat.names) %in% p.ops$which.treat]
        if (length(which.treat) == 0) {
            warning("No numbers in which.treat correspond to treatment values. No treatment pairs will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
    }
    else if (is.character(p.ops$which.treat)) {
        which.treat <- p.ops$treat.names[p.ops$treat.names %in% p.ops$which.treat]
        if (length(which.treat) == 0) {
            warning("No names in which.treat correspond to treatment values. No treatment pairs will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
    }
    else if (is.na(p.ops$which.treat)) {
        which.treat <- character(0)
        p.ops$multi.summary <- TRUE
    }
    else {
        warning("The argument to which.treat must be NA, NULL, or a vector of treatment names or indices. No treatment pairs will be displayed.", call. = FALSE)
        which.treat <- character(0)
        p.ops$multi.summary <- TRUE
    }
    
    if (length(which.treat) == 0) {
        disp.treat.pairs <- character(0)
    }
    else {
        if (p.ops$pairwise) {
            if (length(which.treat) == 1) {
                disp.treat.pairs <- names(m.balance)[sapply(names(m.balance), function(x) any(m.balance[[x]]$print.options$treat.names == which.treat))]
            }
            else {
                disp.treat.pairs <- names(m.balance)[sapply(names(m.balance), function(x) all(m.balance[[x]]$print.options$treat.names %in% which.treat))]
            }
        }
        else {
            if (length(which.treat) == 1) {
                disp.treat.pairs <- names(m.balance)[sapply(names(m.balance), function(x) {
                    treat.names <- m.balance[[x]]$print.options$treat.names
                    any(treat.names[treat.names != "Others"] == which.treat)})]
            }
            else {
                disp.treat.pairs <- names(m.balance)[sapply(names(m.balance), function(x) {
                    treat.names <- m.balance[[x]]$print.options$treat.names
                    all(treat.names[treat.names != "Others"] %in% which.treat)})]
            }
        }
    }
    
    #Printing output
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n  ")
        cat("\n")
    }
    
    if (length(disp.treat.pairs) > 0) {
        headings <- setNames(character(length(disp.treat.pairs)), disp.treat.pairs)
        if (p.ops$pairwise) cat("Balance by treatment pair:\n")
        else cat("Balance by treatment group:\n")
        for (i in disp.treat.pairs) {
            headings[i] <- paste0("\n - - - ", m.balance[[i]]$print.options$treat.names[1]," (0) vs. ",
                                  m.balance[[i]]$print.options$treat.names[2]," (1) - - - \n")
            cat(headings[i])
            do.call(print, c(list(m.balance[[i]]), args))
        }
        cat(paste0(paste(rep(" -", round(max(nchar(headings))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$multi.summary))) {
        s.keep <- as.logical(c(TRUE, 
                               p.ops$un,
                               p.ops$un && !p.ops$disp.adj && !is.null(p.ops$m.threshold),
                               p.ops$un && p.ops$disp.v.ratio, 
                               p.ops$un && !p.ops$disp.adj && !is.null(p.ops$v.threshold), 
                               p.ops$un && p.ops$disp.ks, 
                               p.ops$un && !p.ops$disp.adj && !is.null(p.ops$ks.threshold),
                               rep(c(p.ops$disp.adj, 
                                     p.ops$disp.adj && !is.null(p.ops$m.threshold), 
                                     p.ops$disp.adj && p.ops$disp.v.ratio, 
                                     p.ops$disp.adj && !is.null(p.ops$v.threshold), 
                                     p.ops$disp.adj && p.ops$disp.ks, 
                                     p.ops$disp.adj && !is.null(p.ops$ks.threshold)), p.ops$nweights + !p.ops$disp.adj)))
        
        if (p.ops$disp.bal.tab) {
        cat("Balance summary across all treatment pairs:\n")
        print.data.frame(round_df_char(m.balance.summary[keep.row, s.keep], digits))
        cat("\n")
        }
        
        if (!is.null(nn)) {
            tag <- attr(x$Observations, "tag")
            ss.type <- attr(x$Observations, "ss.type")
            for (i in rownames(x$Observations)) {
                if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
            }
            if (!is.na(x$Observations["Matched (Unweighted)",]) && all(x$Observations["Matched",] == x$Observations["Matched (Unweighted)",])) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)",]
            cat(paste0(tag, ":\n"))
            print.warning <- FALSE
            if (length(ss.type) > 1 && length(unique(ss.type[-1])) > 1) {
                ess <- ifelse(ss.type == "ess", "*", "")
                x$Observations <- setNames(cbind(x$Observations, ess), c(names(x$Observations), ""))
                print.warning <- TRUE
            }
            print.data.frame(replaceNA(x$Observations), digits = digits)
            if (print.warning) cat("* indicates effective sample size")
        }
    }
    
    invisible(x)
    
}
print.bal.tab.msm <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = NULL, which.time, msm.summary = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    args <- c(as.list(environment()), list(...))[-1]
    args <- args[!sapply(args, function(x) identical(x, quote(expr =)))]
    
    call <- x$call
    msm.balance <- x[["Time.Balance"]]
    msm.balance.summary <- x[["Balance.Across.Times"]]
    nn <- x$Observations
    p.ops <- x$print.options
    
    #Prevent expnential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("un cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp.means, "as.is")) {
        if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.means == FALSE && disp.means == TRUE) {
            warning("disp.means cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.means <- disp.means
    }
    if (!identical(disp.v.ratio, "as.is")) {
        if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.v.ratio == FALSE && disp.v.ratio == TRUE) {
            warning("disp.v.ratio cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.v.ratio <- disp.v.ratio
    }
    if (!identical(disp.ks, "as.is")) {
        if (!is.logical(disp.ks)) stop("disp.ks must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.ks == FALSE && disp.ks == TRUE) {
            warning("disp.ks cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.ks <- disp.ks
    }
    if (!identical(msm.summary, "as.is")) {
        if (!is.logical(msm.summary)) stop("msm.summary must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$msm.summary == FALSE && msm.summary == TRUE) {
            warning("msm.summary cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$msm.summary <- msm.summary
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
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (length(p.ops$disp.v.ratio) == 0 || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (length(p.ops$disp.ks) == 0 || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (!is.null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
        }
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (!any(sapply(c(p.ops$m.threshold, 
                              p.ops$v.threshold, 
                              p.ops$ks.threshold, 
                              p.ops$r.threshold), length) > 0)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (length(msm.balance.summary) > 0) {
        if (p.ops$imbalanced.only) {
            keep.row <- rowSums(apply(msm.balance.summary[grepl(".Threshold", names(msm.balance.summary), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        }
        else keep.row <- rep(TRUE, nrow(msm.balance.summary))
    }
    
    if (!missing(which.time) && which.time != "as.is") {
        p.ops$which.time <- which.time
    }
    
    #Checks and Adjustments
    if (length(p.ops$which.time) == 0) 
        which.time <- seq_along(msm.balance)
    else if (is.numeric(p.ops$which.time)) {
        which.time <- seq_along(msm.balance)[seq_along(msm.balance) %in% p.ops$which.time]
        if (length(which.time) == 0) {
            warning("No numbers in which.time are treatment time points. No time points will be displayed.", call. = FALSE)
            which.time <- integer(0)
        }
    }
    else if (is.character(p.ops$which.time)) {
        which.time <- seq_along(msm.balance)[names(msm.balance) %in% p.ops$which.time]
        if (length(which.time) == 0) {
            warning("No names in which.time are treatment names. Displaying all time points instead.", call. = FALSE)
            which.time <- seq_along(msm.balance)
        }
    }
    else if (is.na(p.ops$which.time)) {
        which.time <- integer(0)
        p.ops$msm.summary <- TRUE
    }
    else {
        warning("The argument to which.time must be NA, NULL, or a vector of time point numbers. No time points will be displayed.", call. = FALSE)
        which.time <- integer(0)
        p.ops$msm.summary <- TRUE
    }
    
    #Printing output
    if (!is.null(call)) {
        cat("\nCall:", deparse(call), sep = "\n  ")
        cat("\n")
    }
    
    if (length(which.time) > 0) {
        cat("Balance by Time Point:\n")
        for (i in which.time) {
            cat(paste0("\n - - - Time: ", i, " - - - \n"))
            do.call(print, c(list(x = msm.balance[[i]]), args))
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Time: ", i, " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$msm.summary))) {
        if (!is.na(match("bal.tab.cont", class(x)))) { #continuous
            s.keep <- as.logical(c(TRUE, 
                                   TRUE,
                                   p.ops$un,
                                   p.ops$un && !p.ops$disp.adj && !is.null(p.ops$r.threshold),
                                   rep(c(p.ops$disp.adj, 
                                         p.ops$disp.adj && !is.null(p.ops$r.threshold) 
                                         ), p.ops$nweights + !p.ops$disp.adj)))
        }
        else { #binary
            s.keep <- as.logical(c(TRUE, 
                                   TRUE,
                                   p.ops$un,
                                   p.ops$un && !p.ops$disp.adj && !is.null(p.ops$m.threshold),
                                   p.ops$un && p.ops$disp.v.ratio, 
                                   p.ops$un && !p.ops$disp.adj && !is.null(p.ops$v.threshold), 
                                   p.ops$un && p.ops$disp.ks, 
                                   p.ops$un && !p.ops$disp.adj && !is.null(p.ops$ks.threshold),
                                   rep(c(p.ops$disp.adj, 
                                         p.ops$disp.adj && !is.null(p.ops$m.threshold), 
                                         p.ops$disp.adj && p.ops$disp.v.ratio, 
                                         p.ops$disp.adj && !is.null(p.ops$v.threshold), 
                                         p.ops$disp.adj && p.ops$disp.ks, 
                                         p.ops$disp.adj && !is.null(p.ops$ks.threshold)), p.ops$nweights + !p.ops$disp.adj)))
        }
        
        if (p.ops$disp.bal.tab && length(msm.balance.summary) > 0) {
            cat("Balance summary across all time points:\n")
            print.data.frame(round_df_char(msm.balance.summary[keep.row, s.keep, drop = FALSE], digits))
            cat("\n")
        }
        
        if (!is.null(nn)) {
            print.warning <- FALSE
            cat(paste0(attr(x$Observations[[1]], "tag"), ":\n"))
            
            for (ti in seq_along(x$Observations)) {
                cat(paste0(" - Time ", ti, ":\n"))
                for (i in rownames(x$Observations[[ti]])) {
                    if (all(x$Observations[[ti]][i,] == 0)) x$Observations[[ti]] <- x$Observations[[ti]][rownames(x$Observations[[ti]])!=i,]
                }
                if (!is.na(x$Observations[[ti]]["Matched (Unweighted)",]) && all(x$Observations[[ti]]["Matched",] == x$Observations[[ti]]["Matched (Unweighted)",])) x$Observations[[ti]] <- x$Observations[[ti]][rownames(x$Observations[[ti]])!="Matched (Unweighted)",]
                if (length(attr(x$Observations[[ti]], "ss.type")) > 1 && length(unique(attr(x$Observations[[ti]], "ss.type")[-1])) > 1) {
                    ess <- ifelse(attr(x$Observations[[ti]], "ss.type") == "ess", "*", "")
                    x$Observations[[ti]] <- setNames(cbind(x$Observations[[ti]], ess), c(names(x$Observations[[ti]]), ""))
                    print.warning <- TRUE
                }
                print.data.frame(replaceNA(x$Observations[[ti]]), digits = digits)
            }
            
            if (print.warning) cat("* indicates effective sample size")
        }
    }
    
    invisible(x)
}