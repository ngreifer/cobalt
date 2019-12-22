print.bal.tab <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.sds = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", digits = max(3, getOption("digits") - 3), ...) {

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
    p.ops <- attr(x, "print.options")
    
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
    if (!identical(disp.sds, "as.is")) {
        if (!is.logical(disp.sds)) stop("disp.sds must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.sds == FALSE && disp.sds == TRUE) {
            warning("disp.sds cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.sds <- disp.sds
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
        if (is_not_null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
            baltal.r <- NULL
            maximbal.r <- NULL
        }
    }
    if (!identical(disp.m.threshold, "as.is")) {
        if (!is.logical(disp.m.threshold)) stop("disp.m.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
            baltal.m <- NULL
            maximbal.m <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (is_null(p.ops$disp.v.ratio)|| !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (is_null(p.ops$disp.ks) || !p.ops$disp.ks) {
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
            if (all(sapply(c(p.ops$m.threshold, 
                             p.ops$v.threshold, 
                             p.ops$ks.threshold, 
                             p.ops$r.threshold), is_null))) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    
    if ("bal.tab.cont" %in% class(x)) {
        keep <- setNames(as.logical(c(TRUE, 
                                      p.ops$un && p.ops$disp.means,
                                      p.ops$un && p.ops$disp.sds,
                                      p.ops$un, 
                                      p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$r.threshold), 
                                      rep(c(p.ops$disp.adj && p.ops$disp.means,
                                            p.ops$disp.adj && p.ops$disp.sds,
                                            p.ops$disp.adj, 
                                            p.ops$disp.adj && is_not_null(p.ops$r.threshold)), 
                                          p.ops$nweights + !p.ops$disp.adj))),
                         names(balance))
        
    }
    else {
        keep <- setNames(as.logical(c(TRUE, 
                                      p.ops$un && p.ops$disp.means, 
                                      p.ops$un && p.ops$disp.sds, 
                                      p.ops$un && p.ops$disp.means, 
                                      p.ops$un && p.ops$disp.sds, 
                                      p.ops$un, 
                                      p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$m.threshold),
                                      p.ops$un && p.ops$disp.v.ratio, 
                                      p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                                      p.ops$un && p.ops$disp.ks, 
                                      p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$ks.threshold),
                                      rep(c(p.ops$disp.adj && p.ops$disp.means, 
                                            p.ops$disp.adj && p.ops$disp.sds, 
                                            p.ops$disp.adj && p.ops$disp.means, 
                                            p.ops$disp.adj && p.ops$disp.sds, 
                                            p.ops$disp.adj, 
                                            p.ops$disp.adj && is_not_null(p.ops$m.threshold), 
                                            p.ops$disp.adj && p.ops$disp.v.ratio, 
                                            p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                                            p.ops$disp.adj && p.ops$disp.ks, 
                                            p.ops$disp.adj && is_not_null(p.ops$ks.threshold)), 
                                          p.ops$nweights + !p.ops$disp.adj))),
                         names(balance))
    }
    
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (p.ops$disp.bal.tab) {
        if (p.ops$imbalanced.only) {
            keep.row <- rowSums(apply(balance[grepl(".Threshold", names(balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        }
        else keep.row <- rep(TRUE, nrow(balance))
        
        cat(underline("Balance Measures") %+% "\n")
        if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
        else print.data.frame_(round_df_char(balance[keep.row, keep], digits))
        cat("\n")
    }
    
    if (is_not_null(baltal.r)) {
        cat(underline("Balance tally for correlations") %+% "\n")
        print.data.frame_(baltal.r)
        cat("\n")
    }
    if (is_not_null(maximbal.r)) {
        cat(underline("Variable with the greatest treatment correlation") %+% "\n")
        print.data.frame_(round_df_char(maximbal.r, digits), row.names = FALSE)
        cat("\n")
    }
    if (is_not_null(baltal.m)) {
        cat(underline("Balance tally for mean differences") %+% "\n")
        print.data.frame_(baltal.m)
        cat("\n")
    }
    if (is_not_null(maximbal.m)) {
        cat(underline("Variable with the greatest mean difference") %+% "\n")
        print.data.frame_(round_df_char(maximbal.m, digits), row.names = FALSE)
        cat("\n")
    }
    if (is_not_null(baltal.v)) {
        cat(underline("Balance tally for variance ratios") %+% "\n")
        print.data.frame_(baltal.v, digits)
        cat("\n")
    }
    if (is_not_null(maximbal.v)) {
        cat(underline("Variable with the greatest variance ratio") %+% "\n")
        print.data.frame_(round_df_char(maximbal.v, digits), row.names = FALSE)
        cat("\n")
    }
    if (is_not_null(baltal.ks)) {
        cat(underline("Balance tally for KS statistics") %+% "\n")
        print.data.frame_(baltal.ks, digits)
        cat("\n")
    }
    if (is_not_null(maximbal.ks)) {
        cat(underline("Variable with the greatest KS statistic") %+% "\n")
        print.data.frame_(round_df_char(maximbal.ks, digits), row.names = FALSE)
        cat("\n")
    }
    if (is_not_null(nn)) {
        for (i in seq_len(NROW(nn))) {
            if (all(nn[i,] == 0)) {
                nn <- nn[-i, , drop = FALSE]
                attr(nn, "ss.type") <- attr(nn, "ss.type")[-i]
            }
        }
        if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
            all(check_if_zero(nn["Matched (ESS)",] - nn["Matched (Unweighted)",]))) {
            nn <- nn[rownames(nn)!="Matched (Unweighted)", , drop = FALSE]
            rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
        }
        cat(underline(attr(nn, "tag")) %+% "\n")
        print.warning <- FALSE
        if (length(attr(nn, "ss.type")) > 1 && nunique.gt(attr(nn, "ss.type")[-1], 1)) {
            ess <- ifelse(attr(nn, "ss.type") == "ess", "*", "")
            nn <- setNames(cbind(nn, ess), c(names(nn), ""))
            print.warning <- TRUE
        }
        print.data.frame_(round_df_char(nn, digits = max(0, digits-1)))
        if (print.warning) cat(italic("* indicates effective sample size"))
    }
    invisible(x)
}
print.bal.tab.cluster <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.sds = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = "as.is", digits = max(3, getOption("digits") - 3), ...) {

    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval.parent(.call))
    }
    
    args <- c(as.list(environment()), list(...))[-1]
    
    call <- x$call
    c.balance <- x$Cluster.Balance
    c.balance.summary <- x$Cluster.Summary
    nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
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
    if (!identical(disp.sds, "as.is")) {
        if (!is.logical(disp.sds)) stop("disp.sds must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.sds == FALSE && disp.sds == TRUE) {
            warning("disp.sds cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.sds <- disp.sds
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
        if (is_not_null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
        }
    }
    if (is_null(p.ops$disp.v.ratio) || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
        }
    }
    if (is_null(p.ops$disp.ks) || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$r.threshold) && !disp.r.threshold) {
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
            if (all(sapply(c(p.ops$m.threshold, 
                             p.ops$v.threshold, 
                             p.ops$ks.threshold, 
                             p.ops$r.threshold), is_null))) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (!missing(which.cluster)) {
        if (paste(deparse(substitute(which.cluster)), collapse = "") == ".none") which.cluster <- NA
        else if (paste(deparse(substitute(which.cluster)), collapse = "") == ".all") which.cluster <- NULL
        if (!identical(which.cluster, "as.is")) {
            p.ops$which.cluster <- which.cluster
        }
    }
    
    cluster.funs <- c("min", "mean", "max")
    if (!identical(cluster.fun, "as.is")) {
        p.ops$cluster.fun <- cluster.fun
    }
    if (is_not_null(p.ops$cluster.fun)) {
        if (!is.character(p.ops$cluster.fun)) stop(paste0("cluster.fun must be ", word_list(c(cluster.funs, "as.is"), and.or = "or", quotes = TRUE)))
        p.ops$cluster.fun <- match_arg(tolower(p.ops$cluster.fun), cluster.funs, several.ok = TRUE)
        if (is_null(p.ops$cluster.fun)) {
            warning("There were no valid entries to cluster.fun. Using the default cluster.funs instead.", call. = FALSE)
            if (p.ops$abs) p.ops$cluster.fun <- c("mean", "max")
            else p.ops$cluster.fun <- c("min", "mean", "max")
        }
    }
    else {
        if (p.ops$abs) p.ops$cluster.fun <- c("mean", "max")
        else p.ops$cluster.fun <- c("min", "mean", "max")
    }
    
    #Checks and Adjustments
    if (is_null(p.ops$which.cluster)) 
        which.cluster <- seq_along(c.balance)
    else if (anyNA(p.ops$which.cluster)) {
        which.cluster <- integer(0)
        if (!p.ops$cluster.summary) message("No clusters will be displayed; displaying the summary across clusters.")
        p.ops$cluster.summary <- TRUE
    }
    else if (is.numeric(p.ops$which.cluster)) {
        which.cluster <- seq_along(c.balance)[seq_along(c.balance) %in% p.ops$which.cluster]
        if (is_null(which.cluster)) {
            warning("No indices in which.cluster are cluster indices. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(c.balance)
        }
    }
    else if (is.character(p.ops$which.cluster)) {
        which.cluster <- seq_along(c.balance)[names(c.balance) %in% p.ops$which.cluster]
        if (is_null(which.cluster)) {
            warning("No names in which.cluster are cluster names. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(c.balance)
        }
    }
    else {
        warning("The argument to which.cluster must be NA, NULL, or a vector of cluster indicies or cluster names. Displaying all clusters instead.", call. = FALSE)
        which.cluster <- seq_along(c.balance)
    }
    
    # if (p.ops$cluster.summary && is_null(c.balance.summary)) {
    #     warning("No summary across clusters was produced. This can occur if cluster.summary is FALSE and quick is TRUE.", call. = FALSE)
    # }
    
    #Printing
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(which.cluster)) {
        cat(underline("Balance by cluster") %+% "\n")
        for (i in which.cluster) {

            cat("\n - - - " %+% italic("Cluster: " %+% names(c.balance)[i]) %+% " - - - \n")
            do.call(print, c(list(c.balance[[i]]), args), quote = TRUE)
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Cluster: ", names(c.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$cluster.summary)) && is_not_null(c.balance.summary)) {
        CF <- setNames(cluster.funs %in% p.ops$cluster.fun, cluster.funs)
        if ("bal.tab.cont" %in% class(c.balance[[1]])) {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && CF["min"],
                                   p.ops$un && CF["mean"],
                                   p.ops$un && CF["max"],
                                   rep(c(p.ops$disp.adj && CF["min"],
                                         p.ops$disp.adj && CF["mean"],
                                         p.ops$disp.adj && CF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        else {
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && CF["min"],
                                   p.ops$un && CF["mean"],
                                   p.ops$un && CF["max"],
                                   p.ops$un && p.ops$disp.v.ratio && CF["min"],
                                   p.ops$un && p.ops$disp.v.ratio && CF["mean"],
                                   p.ops$un && p.ops$disp.v.ratio && CF["max"],
                                   p.ops$un && p.ops$disp.ks && CF["min"],
                                   p.ops$un && p.ops$disp.ks && CF["mean"],
                                   p.ops$un && p.ops$disp.ks && CF["max"],
                                   rep(c(p.ops$disp.adj && CF["min"],
                                         p.ops$disp.adj && CF["mean"],
                                         p.ops$disp.adj && CF["max"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && CF["min"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && CF["mean"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && CF["max"],
                                         p.ops$disp.adj && p.ops$disp.ks && CF["min"],
                                         p.ops$disp.adj && p.ops$disp.ks && CF["mean"],
                                         p.ops$disp.adj && p.ops$disp.ks && CF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all clusters") %+% "\n")
            print.data.frame_(round_df_char(c.balance.summary[, s.keep], digits))
            cat("\n")
        }
        
        if (is_not_null(nn)) {
            for (i in rownames(nn)) {
                if (all(nn[i,] == 0)) nn <- nn[rownames(nn)!=i,]
            }
            if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
                all(check_if_zero(nn["Matched (ESS)",] - nn["Matched (Unweighted)",]))) {
                nn <- nn[rownames(nn)!="Matched (Unweighted)", , drop = FALSE]
                rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
            }
            cat(underline(attr(nn, "tag")) %+% "\n")
            print.warning <- FALSE
            if (length(attr(nn, "ss.type")) > 1 && nunique.gt(attr(nn, "ss.type")[-1], 1)) {
                ess <- ifelse(attr(nn, "ss.type") == "ess", "*", "")
                nn <- setNames(cbind(nn, ess), c(names(nn), ""))
                print.warning <- TRUE
            }
            print.data.frame_(round_df_char(nn, digits = max(0, digits-1)))
            if (print.warning) cat(italic("* indicates effective sample size"))
        }
    }
    
    invisible(x)
}
print.bal.tab.imp <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.sds = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.imp, imp.summary = "as.is", imp.fun = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval.parent(.call))
    }
    
    args <- c(as.list(environment()), list(...))[-1]
    
    call <- x$call
    i.balance <- x[["Imputation.Balance"]]
    i.balance.summary <- x[["Balance.Across.Imputations"]]
    nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
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
    if (!identical(disp.sds, "as.is")) {
        if (!is.logical(disp.sds)) stop("disp.sds must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.sds == FALSE && disp.sds == TRUE) {
            warning("disp.sds cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.sds <- disp.sds
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
        if (is_not_null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (is_null(p.ops$disp.v.ratio) || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (is_null(p.ops$disp.ks) || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$r.threshold) && !disp.r.threshold) {
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
            if (all(sapply(c(p.ops$m.threshold, 
                             p.ops$v.threshold, 
                             p.ops$ks.threshold, 
                             p.ops$r.threshold), is_null))) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (!missing(which.imp)) {
        if (paste(deparse(substitute(which.imp)), collapse = "") == ".none") which.imp <- NA
        else if (paste(deparse(substitute(which.imp)), collapse = "") == ".all") which.imp <- NULL
        if (!identical(which.imp, "as.is")) {
            p.ops$which.imp <- which.imp
        }
    }
    
    imp.funs <- c("min", "mean", "max")
    if (!identical(imp.fun, "as.is")) {
        p.ops$imp.fun <- imp.fun
    }
    if (is_not_null(p.ops$imp.fun)) {
        if (!is.character(p.ops$imp.fun)) stop(paste0("imp.fun must be ", word_list(c(imp.funs, "as.is"), and.or = "or", quotes = TRUE)))
        p.ops$imp.fun <- match_arg(tolower(p.ops$imp.fun), imp.funs, several.ok = TRUE)
        if (is_null(p.ops$imp.fun)) {
            warning("There were no valid entries to imp.fun. Using the default imp.funs instead.", call. = FALSE)
            if (p.ops$abs) p.ops$imp.fun <- c("mean", "max")
            else p.ops$imp.fun <- c("min", "mean", "max")
        }
    }
    else {
        if (p.ops$abs) p.ops$imp.fun <- c("mean", "max")
        else p.ops$imp.fun <- c("min", "mean", "max")
    }
    
    #Checks and Adjustments
    if (is_null(p.ops$which.imp)) 
        which.imp <- seq_along(i.balance)
    else if (anyNA(p.ops$which.imp)) {
        which.imp <- integer(0)
        if (!p.ops$imp.summary) message("No imputations will be displayed; displaying the summary across imputations.")
        p.ops$imp.summary <- TRUE
    }
    else if (is.numeric(p.ops$which.imp)) {
        which.imp <- seq_along(i.balance)[seq_along(i.balance) %in% p.ops$which.imp]
        if (is_null(which.imp)) {
            warning("No numbers in which.imp are imputation numbers. No imputations will be displayed.", call. = FALSE)
            which.imp <- integer(0)
            if (!p.ops$imp.summary) message("No imputations will be displayed; displaying the summary across imputations.")
            p.ops$imp.summary <- TRUE
        }
    }
    else {
        warning("The argument to which.imp must be NA, NULL, or a vector of imputation numbers. No imputations will be displayed. Displaying the summary across imputations", call. = FALSE)
        which.imp <- integer(0)
        p.ops$imp.summary <- TRUE
    }
    
    # if (p.ops$imp.summary && is_null(i.balance.summary)) {
    #     warning("No summary across imputations was produced. This can occur if imp.summary is FALSE and quick is TRUE.", call. = FALSE)
    # }
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(which.imp)) {
        cat(underline("Balance by imputation") %+% "\n")
        for (i in which.imp) {
            cat("\n - - - " %+% italic("Imputation " %+% names(i.balance)[i]) %+% " - - - \n")
            do.call(print, c(list(i.balance[[i]]), args), quote = TRUE)
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Imputation: ", names(i.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$imp.summary)) && is_not_null(i.balance.summary)) {
        IF <- setNames(imp.funs %in% p.ops$imp.fun, imp.funs)
        if ("bal.tab.cont" %in% class(i.balance[[1]])) { #continuous
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && IF["min"],
                                   p.ops$un && IF["mean"],
                                   p.ops$un && IF["max"],
                                   rep(c(p.ops$disp.adj && IF["min"],
                                         p.ops$disp.adj && IF["mean"],
                                         p.ops$disp.adj && IF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        else { #binary
            s.keep <- as.logical(c(TRUE, 
                                   p.ops$un && IF["min"],
                                   p.ops$un && IF["mean"],
                                   p.ops$un && IF["max"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["min"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["mean"],
                                   p.ops$un && p.ops$disp.v.ratio && IF["max"],
                                   p.ops$un && p.ops$disp.ks && IF["min"],
                                   p.ops$un && p.ops$disp.ks && IF["mean"],
                                   p.ops$un && p.ops$disp.ks && IF["max"],
                                   rep(c(p.ops$disp.adj && IF["min"],
                                         p.ops$disp.adj && IF["mean"],
                                         p.ops$disp.adj && IF["max"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["min"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["mean"],
                                         p.ops$disp.adj && p.ops$disp.v.ratio && IF["max"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["min"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["mean"],
                                         p.ops$disp.adj && p.ops$disp.ks && IF["max"]), p.ops$nweights + !p.ops$disp.adj)))
        }
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all imputations") %+% "\n")
            print.data.frame_(round_df_char(i.balance.summary[, s.keep, drop = FALSE], digits))
            cat("\n")
        }
        
        if (is_not_null(nn)) {
            for (i in rownames(nn)) {
                if (all(nn[i,] == 0)) nn <- nn[rownames(nn)!=i,]
            }
            if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
                all(check_if_zero(nn["Matched (ESS)",] - nn["Matched (Unweighted)",]))) {
                nn <- nn[rownames(nn)!="Matched (Unweighted)", , drop = FALSE]
                rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
            }
            cat(underline(attr(nn, "tag")) %+% "\n")
            print.warning <- FALSE
            if (length(attr(nn, "ss.type")) > 1 && nunique.gt(attr(nn, "ss.type")[-1], 1)) {
                ess <- ifelse(attr(nn, "ss.type") == "ess", "*", "")
                nn <- setNames(cbind(nn, ess), c(names(nn), ""))
                print.warning <- TRUE
            }
            print.data.frame_(round_df_char(nn, digits = max(0, digits-1)))
            if (print.warning) cat(italic("* indicates effective sample size"))
        }
    }
    
    invisible(x)
    
}
print.bal.tab.multi <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.sds = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.treat, multi.summary = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval.parent(.call))
    }
    
    args <- c(as.list(environment()), list(...))[-1]
    
    call <- x$call
    m.balance <- x[["Pair.Balance"]]
    m.balance.summary <- x[["Balance.Across.Pairs"]]
    nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
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
    if (!identical(disp.sds, "as.is")) {
        if (!is.logical(disp.sds)) stop("disp.sds must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.sds == FALSE && disp.sds == TRUE) {
            warning("disp.sds cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.sds <- disp.sds
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
        if (is_not_null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (is_null(p.ops$disp.v.ratio) || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (is_null(p.ops$disp.ks) || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (is_not_null(m.balance.summary)) {
        if (p.ops$disp.bal.tab) {
            if (!identical(imbalanced.only, "as.is")) {
                if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
                p.ops$imbalanced.only <- imbalanced.only
            }
            if (p.ops$imbalanced.only) {
                if (all(sapply(c(p.ops$m.threshold, 
                                 p.ops$v.threshold, 
                                 p.ops$ks.threshold, 
                                 p.ops$r.threshold), is_null))) {
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
    }
    
    if (!missing(which.treat)) {
        if (paste(deparse(substitute(which.treat)), collapse = "") == ".none") which.treat <- NA
        else if (paste(deparse(substitute(which.treat)), collapse = "") == ".all") which.treat <- NULL
        if (!identical(which.treat, "as.is")) {
            p.ops$which.treat <- which.treat
        }
    }
    
    #Checks and Adjustments
    if (is_null(p.ops$which.treat)) 
        which.treat <- p.ops$treat_names_multi
    else if (anyNA(p.ops$which.treat)) {
        which.treat <- character(0)
        if (!p.ops$multi.summary) message("No treatment pairs will be displayed; displaying the summary across treatments.")
        p.ops$multi.summary <- TRUE
    }
    else if (is.numeric(p.ops$which.treat)) {
        which.treat <- p.ops$treat_names_multi[seq_along(p.ops$treat_names_multi) %in% p.ops$which.treat]
        if (is_null(which.treat)) {
            warning("No numbers in which.treat correspond to treatment values. No treatment pairs will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
    }
    else if (is.character(p.ops$which.treat)) {
        which.treat <- p.ops$treat_names_multi[p.ops$treat_names_multi %in% p.ops$which.treat]
        if (is_null(which.treat)) {
            warning("No names in which.treat correspond to treatment values. No treatment pairs will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
    }
    else {
        warning("The argument to which.treat must be NA, NULL, or a vector of treatment names or indices. No treatment pairs will be displayed.", call. = FALSE)
        which.treat <- character(0)
        if (!p.ops$multi.summary) message("No treatment pairs will be displayed; displaying the summary across treatments.")
        p.ops$multi.summary <- TRUE
    }
    
    # if (p.ops$multi.summary && is_null(m.balance.summary)) {
    #     warning("No summary across treatment pairs was produced. This can occur if multi.summary is FALSE and quick is TRUE.", call. = FALSE)
    # }
    
    if (is_null(which.treat)) {
        disp.treat.pairs <- character(0)
    }
    else {
        if (p.ops$pairwise) {
            if (length(which.treat) == 1) {
                disp.treat.pairs <- names(m.balance)[sapply(names(m.balance), function(x) any(attr(m.balance[[x]], "print.options")$treat_names == which.treat))]
            }
            else {
                disp.treat.pairs <- names(m.balance)[sapply(names(m.balance), function(x) all(attr(m.balance[[x]], "print.options")$treat_names %in% which.treat))]
            }
        }
        else {
            if (length(which.treat) == 1) {
                disp.treat.pairs <- names(m.balance)[sapply(names(m.balance), function(x) {
                    treat_names <- attr(m.balance[[x]], "print.options")$treat_names
                    any(treat_names[treat_names != "Others"] == which.treat)})]
            }
            else {
                disp.treat.pairs <- names(m.balance)[sapply(names(m.balance), function(x) {
                    treat_names <- attr(m.balance[[x]], "print.options")$treat_names
                    all(treat_names[treat_names != "Others"] %in% which.treat)})]
            }
        }
    }
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(disp.treat.pairs)) {
        headings <- setNames(character(length(disp.treat.pairs)), disp.treat.pairs)
        if (p.ops$pairwise) cat(underline("Balance by treatment pair") %+% "\n")
        else cat(underline("Balance by treatment group") %+% "\n")
        for (i in disp.treat.pairs) {
            headings[i] <- "\n - - - " %+% italic(attr(m.balance[[i]], "print.options")$treat_names[1] %+% " (0) vs. " %+%
                                                      attr(m.balance[[i]], "print.options")$treat_names[2] %+% " (1)") %+% " - - - \n"
            cat(headings[i])
            do.call(print, c(list(m.balance[[i]]), args), quote = TRUE)
        }
        cat(paste0(paste(rep(" -", round(max(nchar(headings))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$multi.summary)) && is_not_null(m.balance.summary)) {
        s.keep <- as.logical(c(TRUE, 
                               p.ops$un,
                               p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$m.threshold),
                               p.ops$un && p.ops$disp.v.ratio, 
                               p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                               p.ops$un && p.ops$disp.ks, 
                               p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$ks.threshold),
                               rep(c(p.ops$disp.adj, 
                                     p.ops$disp.adj && is_not_null(p.ops$m.threshold), 
                                     p.ops$disp.adj && p.ops$disp.v.ratio, 
                                     p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                                     p.ops$disp.adj && p.ops$disp.ks, 
                                     p.ops$disp.adj && is_not_null(p.ops$ks.threshold)), p.ops$nweights + !p.ops$disp.adj)))
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all treatment pairs") %+% "\n")
            if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(m.balance.summary[keep.row, s.keep], digits))
            cat("\n")
        }
        
        if (is_not_null(nn)) {
            tag <- attr(nn, "tag")
            ss.type <- attr(nn, "ss.type")
            for (i in rownames(nn)) {
                if (all(nn[i,] == 0)) nn <- nn[rownames(nn)!=i,]
            }
            if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
                all(check_if_zero(nn["Matched (ESS)",] - nn["Matched (Unweighted)",]))) {
                nn <- nn[rownames(nn)!="Matched (Unweighted)", , drop = FALSE]
                rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
            }
            cat(underline(tag) %+% "\n")
            print.warning <- FALSE
            if (length(ss.type) > 1 && nunique.gt(ss.type[-1], 1)) {
                ess <- ifelse(ss.type == "ess", "*", "")
                nn <- setNames(cbind(nn, ess), c(names(nn), ""))
                print.warning <- TRUE
            }
            print.data.frame_(round_df_char(nn, digits = max(0, digits-1)))
            if (print.warning) cat(italic("* indicates effective sample size"))
        }
    }
    
    invisible(x)
    
}
print.bal.tab.msm <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.sds = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = "as.is", which.time, msm.summary = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval.parent(.call))
    }
    
    args <- c(as.list(environment()), list(...))[-1]
    args <- args[!sapply(args, function(x) identical(x, quote(expr =)))]
    
    call <- x$call
    msm.balance <- x[["Time.Balance"]]
    msm.balance.summary <- x[["Balance.Across.Times"]]
    nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
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
    if (!identical(disp.sds, "as.is")) {
        if (!is.logical(disp.sds)) stop("disp.sds must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.sds == FALSE && disp.sds == TRUE) {
            warning("disp.sds cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.sds <- disp.sds
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
        if (is_not_null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (is_null(p.ops$disp.v.ratio) || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (is_null(p.ops$disp.ks) || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$r.threshold) && !disp.r.threshold) {
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
            if (all(sapply(c(p.ops$m.threshold, 
                             p.ops$v.threshold, 
                             p.ops$ks.threshold, 
                             p.ops$r.threshold), is_null))) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (is_not_null(msm.balance.summary)) {
        if (p.ops$imbalanced.only) {
            keep.row <- rowSums(apply(msm.balance.summary[grepl(".Threshold", names(msm.balance.summary), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        }
        else keep.row <- rep(TRUE, nrow(msm.balance.summary))
    }
    
    if (!missing(which.time)) {
        if (paste(deparse(substitute(which.time)), collapse = "") == ".none") which.time <- NA
        else if (paste(deparse(substitute(which.time)), collapse = "") == ".all") which.time <- NULL
        if (!identical(which.time, "as.is")) {
            p.ops$which.time <- which.time
        }
    }
    
    #Checks and Adjustments
    if (is_null(p.ops$which.time)) 
        which.time <- seq_along(msm.balance)
    else if (anyNA(p.ops$which.time)) {
        which.time <- integer(0)
        if (!p.ops$msm.summary) message("No time points will be displayed; displaying the summary across time points.")
        p.ops$msm.summary <- TRUE
    }
    else if (is.numeric(p.ops$which.time)) {
        which.time <- seq_along(msm.balance)[seq_along(msm.balance) %in% p.ops$which.time]
        if (is_null(which.time)) {
            warning("No numbers in which.time are treatment time points. No time points will be displayed.", call. = FALSE)
            which.time <- integer(0)
        }
    }
    else if (is.character(p.ops$which.time)) {
        which.time <- seq_along(msm.balance)[names(msm.balance) %in% p.ops$which.time]
        if (is_null(which.time)) {
            warning("No names in which.time are treatment names. Displaying all time points instead.", call. = FALSE)
            which.time <- seq_along(msm.balance)
        }
    }
    else {
        warning("The argument to which.time must be NA, NULL, or a vector of time point numbers. No time points will be displayed.", call. = FALSE)
        which.time <- integer(0)
        if (!p.ops$msm.summary) message("No time points will be displayed; displaying the summary across time points.")
        p.ops$msm.summary <- TRUE
    }
    
    # if (p.ops$msm.summary && is_null(msm.balance.summary)) {
    #     warning("No summary across time points was produced. This can occur if msm.summary is FALSE and quick is TRUE or if the treatments are not all binary or all continuous.", call. = FALSE)
    # }
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(which.time)) {
        cat(underline("Balance by Time Point") %+% "\n")
        for (i in which.time) {
            cat("\n - - - " %+% italic("Time: " %+% as.character(i)) %+% " - - - \n")
            do.call(print, c(list(x = msm.balance[[i]]), args), quote = TRUE)
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Time: ", i, " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$msm.summary)) && is_not_null(msm.balance.summary)) {
        if (all(sapply(msm.balance, function(x) "bal.tab.cont" %in% class(x)))) { #continuous
            s.keep <- as.logical(c(TRUE, 
                                   TRUE,
                                   p.ops$un,
                                   p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$r.threshold),
                                   rep(c(p.ops$disp.adj, 
                                         p.ops$disp.adj && is_not_null(p.ops$r.threshold) 
                                   ), p.ops$nweights + !p.ops$disp.adj)))
            
        }
        else { #binary
            s.keep <- as.logical(c(TRUE, 
                                   TRUE,
                                   p.ops$un,
                                   p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$m.threshold),
                                   p.ops$un && p.ops$disp.v.ratio, 
                                   p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                                   p.ops$un && p.ops$disp.ks, 
                                   p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$ks.threshold),
                                   rep(c(p.ops$disp.adj, 
                                         p.ops$disp.adj && is_not_null(p.ops$m.threshold), 
                                         p.ops$disp.adj && p.ops$disp.v.ratio, 
                                         p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                                         p.ops$disp.adj && p.ops$disp.ks, 
                                         p.ops$disp.adj && is_not_null(p.ops$ks.threshold)), p.ops$nweights + !p.ops$disp.adj)))
        }
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all time points") %+% "\n")
            if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(msm.balance.summary[keep.row, s.keep, drop = FALSE], digits))
            cat("\n")
        }
        
        if (is_not_null(nn)) {
            print.warning <- FALSE
            cat(underline(attr(nn[[1]], "tag")) %+% "\n")
            
            for (ti in seq_along(nn)) {
                cat(" - " %+% italic("Time " %+% as.character(ti)) %+% "\n")
                for (i in rownames(nn[[ti]])) {
                    if (all(nn[[ti]][i,] == 0)) nn[[ti]] <- nn[[ti]][rownames(nn[[ti]])!=i,]
                }
                if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn[[ti]])) && 
                    all(check_if_zero(nn[[ti]]["Matched (ESS)",] - nn[[ti]]["Matched (Unweighted)",]))) {
                    nn[[ti]] <- nn[[ti]][rownames(nn[[ti]])!="Matched (Unweighted)", , drop = FALSE]
                    rownames(nn[[ti]])[rownames(nn[[ti]]) == "Matched (ESS)"] <- "Matched"
                }
                if (length(attr(nn[[ti]], "ss.type")) > 1 && nunique.gt(attr(nn[[ti]], "ss.type")[-1], 1)) {
                    ess <- ifelse(attr(nn[[ti]], "ss.type") == "ess", "*", "")
                    nn[[ti]] <- setNames(cbind(nn[[ti]], ess), c(names(nn[[ti]]), ""))
                    print.warning <- TRUE
                }
                print.data.frame_(round_df_char(nn[[ti]], digits = max(0, digits-1)))
            }
            
            if (print.warning) cat(italic("* indicates effective sample size"))
        }
    }
    
    invisible(x)
}

print.bal.tab.subclass <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.sds = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", disp.subclass = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
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
    p.ops <- attr(x, "print.options")
    
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
    if (!identical(disp.sds, "as.is")) {
        if (!is.logical(disp.sds)) stop("disp.sds must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.sds == FALSE && disp.sds == TRUE) {
            warning("disp.sds cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.sds <- disp.sds
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
        if (is_not_null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
            baltal.r.subclass <- NULL
            maximbal.r.subclass <- NULL
        }
    }
    if (!identical(disp.m.threshold, "as.is")) {
        if (!is.logical(disp.m.threshold)) stop("disp.m.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
            baltal.m.subclass <- NULL
            maximbal.m.subclass <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (is_null(p.ops$disp.v.ratio) || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (is_null(p.ops$disp.ks) || !p.ops$disp.ks) {
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
            if (all(sapply(c(p.ops$m.threshold, 
                             p.ops$v.threshold, 
                             p.ops$ks.threshold, 
                             p.ops$r.threshold), is_null))) {
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
    
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (p.ops$disp.bal.tab) {
        if (p.ops$disp.subclass) {
            if ("bal.tab.cont" %in% class(x)) {
                s.keep <- as.logical(c(TRUE, 
                                       p.ops$disp.means, 
                                       p.ops$disp.sds, 
                                       p.ops$disp.adj, 
                                       is_not_null(p.ops$r.threshold)))
            }
            else {
                s.keep <- as.logical(c(TRUE, 
                                       p.ops$disp.means, 
                                       p.ops$disp.sds, 
                                       p.ops$disp.means, 
                                       p.ops$disp.sds, 
                                       p.ops$disp.adj, 
                                       is_not_null(p.ops$m.threshold), 
                                       p.ops$disp.adj  &&  p.ops$disp.v.ratio, 
                                       is_not_null(p.ops$v.threshold), 
                                       p.ops$disp.adj  &&  p.ops$disp.ks, 
                                       is_not_null(p.ops$ks.threshold)))
            }
            cat(underline("Balance by subclass"))
            for (i in names(s.balance)) {
                if (p.ops$imbalanced.only) {
                    keep.row <- rowSums(apply(s.balance[[i]][grepl(".Threshold", names(s.balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
                }
                else keep.row <- rep(TRUE, nrow(s.balance[[i]]))
                
                cat("\n - - - " %+% italic("Subclass " %+% as.character(i)) %+% " - - - \n")
                if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
                else print.data.frame_(round_df_char(s.balance[[i]][keep.row, s.keep, drop = FALSE], digits))
            }
            cat("\n")
        }
        
        if (is_not_null(b.a.subclass)) {
            if (p.ops$imbalanced.only) {
                keep.row <- rowSums(apply(b.a.subclass[grepl(".Threshold", names(b.a.subclass), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
            }
            else keep.row <- rep(TRUE, nrow(b.a.subclass))
            if ("bal.tab.cont" %in% class(x)) {
                a.s.keep <- setNames(as.logical(c(TRUE, 
                                                  p.ops$un && p.ops$disp.means,
                                                  p.ops$un && p.ops$disp.sds,
                                                  p.ops$un, 
                                                  p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$r.threshold), 
                                                  rep(c(p.ops$disp.adj && p.ops$disp.means,
                                                        p.ops$disp.adj && p.ops$disp.sds,
                                                        p.ops$disp.adj, 
                                                        p.ops$disp.adj && is_not_null(p.ops$r.threshold)), 
                                                      p.ops$disp.adj))),
                                     names(b.a.subclass))
            }
            else {
            a.s.keep <- setNames(as.logical(c(TRUE, 
                                          p.ops$un && p.ops$disp.means, 
                                          p.ops$un && p.ops$disp.sds, 
                                          p.ops$un && p.ops$disp.means, 
                                          p.ops$un && p.ops$disp.sds, 
                                          p.ops$un, 
                                          p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$m.threshold),
                                          p.ops$un && p.ops$disp.v.ratio, 
                                          p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                                          p.ops$un && p.ops$disp.ks, 
                                          p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$ks.threshold),
                                          rep(c(p.ops$disp.adj && p.ops$disp.means, 
                                                p.ops$disp.adj && p.ops$disp.sds, 
                                                p.ops$disp.adj && p.ops$disp.means, 
                                                p.ops$disp.adj && p.ops$disp.sds, 
                                                p.ops$disp.adj, 
                                                p.ops$disp.adj && is_not_null(p.ops$m.threshold), 
                                                p.ops$disp.adj && p.ops$disp.v.ratio, 
                                                p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                                                p.ops$disp.adj && p.ops$disp.ks, 
                                                p.ops$disp.adj && is_not_null(p.ops$ks.threshold)), 
                                              p.ops$disp.adj))),
                             names(b.a.subclass))
            }
            
            cat(underline("Balance measures across subclasses") %+% "\n")
            if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(b.a.subclass[keep.row, a.s.keep, drop = FALSE], digits))
            cat("\n")
        }
    }
    if (is_not_null(baltal.r.subclass)) {
        cat(underline("Balance tally for correlations across subclasses") %+% "\n")
        print.data.frame_(baltal.r.subclass)
        cat("\n")
    }
    if (is_not_null(maximbal.r.subclass)) {
        cat(underline("Variable with the greatest treatment correlation across subclasses") %+% "\n")
        print.data.frame_(round_df_char(maximbal.r.subclass, digits), row.names = TRUE)
        cat("\n")
    }
    if (is_not_null(baltal.m.subclass)) {
        cat(underline("Balance tally for mean differences across subclasses") %+% "\n")
        print.data.frame_(baltal.m.subclass)
        cat("\n")
    }
    if (is_not_null(maximbal.m.subclass)) {
        cat(underline("Variable with the greatest mean difference across subclasses") %+% "\n")
        print.data.frame_(round_df_char(maximbal.m.subclass, digits), row.names = TRUE)
        cat("\n")
    }
    if (is_not_null(baltal.v.subclass)) {
        cat(underline("Balance tally for variance ratios across subclasses") %+% "\n")
        print.data.frame_(baltal.v.subclass)
        cat("\n")
    }
    if (is_not_null(maximbal.v.subclass)) {
        cat(underline("Variable with the greatest variance ratios across subclasses") %+% "\n")
        print.data.frame_(round_df_char(maximbal.v.subclass, digits), row.names = TRUE)
        cat("\n")
    }
    if (is_not_null(baltal.ks.subclass)) {
        cat(underline("Balance tally for KS statistics across subclasses") %+% "\n")
        print.data.frame_(baltal.ks.subclass)
        cat("\n")
    }
    if (is_not_null(maximbal.ks.subclass)) {
        cat(underline("Variable with the greatest KS statistc across subclasses") %+% "\n")
        print.data.frame_(round_df_char(maximbal.ks.subclass, digits), row.names = TRUE)
        cat("\n")
    }
    
    if (is_not_null(s.nn)) {
        cat(underline(attr(s.nn, "tag")) %+% "\n")
        print.data.frame_(round_df_char(s.nn, digits = max(0, digits-1)))
    }
    
    invisible(x)
}
