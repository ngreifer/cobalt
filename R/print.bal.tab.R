print.bal.tab <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    balance <- x$Balance
    baltal.m <- x$Balanced.Means
    maximbal.m <- x$Max.Imbalance.Means
    baltal.v <- x$Balanced.Variances
    maximbal.v <- x$Max.Imbalance.Variances
    nn <- x$Observations
    p.ops <- x$print.options
    
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
    cat("\n")
    
    if (!is.null(baltal.m)) {
        cat("\nBalance tally for mean differences:\n")
        print.data.frame(x$Balanced.Means)
        cat("\n")
    }
    if (!is.null(maximbal.m)) {
        cat("\nVariable with the greatest mean difference:\n")
        print.data.frame(round_df(x$Max.Imbalance.Means, digits))
        cat("\n")
    }
    if (!is.null(baltal.v)) {
        cat("\nBalance tally for variance ratios:\n")
        print.data.frame(x$Balanced.Variances, digits)
        cat("\n")
    }
    if (!is.null(maximbal.v)) {
        cat("\nVariable with the greatest variance ratio:\n")
        print.data.frame(round_df(x$Max.Imbalance.Variances, digits))
        cat("\n")
    }
    if (!is.null(nn)) {
        cat(paste0("\n", attr(x$Observations, "tag"), "\n"))
        print.data.frame(replaceNA(x$Observations), digits = digits)
        cat("\n")
    }
    invisible(x)
}
print.bal.tab.cont <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    balance <- x$Balance
    baltal.r <- x$Balanced.Corr
    maximbal.r <- x$Max.Imbalance.Corr
    nn <- x$Observations
    p.ops <- x$print.options
    
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
    cat("\n")
    
    if (!is.null(baltal.r)) {
        cat("\nBalance tally for correlations:\n")
        print.data.frame(x$Balanced.Corr)
        cat("\n")
    }
    if (!is.null(maximbal.r)) {
        cat("\nVariable with the greatest treatment correlation:\n")
        print.data.frame(round_df(x$Max.Imbalance.Corr, digits))
        cat("\n")
    }

    if (!is.null(nn)) {
        cat(paste0("\n", attr(x$Observations, "tag"), "\n"))
        print.data.frame(replaceNA(x$Observations), digits = digits)
        cat("\n")
    }
    invisible(x)
}
print.bal.tab.subclass <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    s.balance <- x$Subclass.Balance
    b.a.subclass <- x$Balance.Across.Subclass
    baltal.m.subclass <- x$Balanced.Means.Subclass
    maximbal.m.subclass <- x$Max.Imbalance.Means.Subclass
    baltal.v.subclass <- x$Balanced.Variances.Subclass
    maximbal.v.subclass <- x$Max.Imbalance.Variances.Subclass
    s.nn <- x$Subclass.Observations
    p.ops <- x$print.options
    
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
    cat("\n")
    
    if (!is.null(baltal.m.subclass)) {
        cat("\nBalance tally for mean differences across subclasses:\n")
        print.data.frame(baltal.m.subclass)
        cat("\n")
    }
    if (!is.null(maximbal.m.subclass)) {
        cat("\nVariable with the greatest mean difference across subclasses:\n")
        print.data.frame(round_df(maximbal.m.subclass, digits))
        cat("\n")
    }
    if (!is.null(baltal.v.subclass)) {
        cat("\nBalance tally for variance ratios across subclasses:\n")
        print.data.frame(baltal.v.subclass)
        cat("\n")
    }
    if (!is.null(maximbal.v.subclass)) {
        cat("\nVariable with the greatest variance ratios across subclasses:\n")
        print.data.frame(round_df(maximbal.v.subclass, digits))
        cat("\n")
    }
    
    if (!is.null(s.nn)) {
        cat(paste0("\n", attr(x$Subclass.Observations, "tag"), "\n"))
        print.data.frame(replaceNA(x$Subclass.Observations), digits = digits)
        cat("\n")
    }
    
    invisible(x)
}
print.bal.tab.cluster <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    #Figure out how to print bal.tab for clusters with subclassification
    call <- x$call
    c.balance <- x$Cluster.Balance
    c.balance.summary <- x$Cluster.Summary
    p.ops <- x$print.options
    
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
        cat("\nBalance summary across all clusters after adjustment:\n")
        print.data.frame(replaceNA(round_df(c.balance.summary, digits)))
        cat("\n")
    }
    invisible(x)
}
print.bal.tab.cont.cluster <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    call <- x$call
    c.balance <- x$Cluster.Balance
    c.balance.summary <- x$Cluster.Summary
    p.ops <- x$print.options
    
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
        cat("\nBalance summary across all clusters after adjustment:\n")
        print.data.frame(replaceNA(round_df(c.balance.summary, digits)))
        cat("\n")
    }
    invisible(x)
}