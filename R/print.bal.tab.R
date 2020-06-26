print.bal.tab <- function(x, imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", stats = "as.is", disp.thresholds = "as.is", disp = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    A <- list(...)
    call <- x$call
    p.ops <- attr(x, "print.options")
    balance <- x$Balance
    
    baltal <- maximbal <- list()
    for (s in p.ops$compute) {
        baltal[[s]] <- x[[paste.("Balanced", s)]]
        maximbal[[s]] <- x[[paste.("Max.Imbalance", s)]]
    }
    nn <- x$Observations
    
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!rlang::is_bool(un)) stop("'un' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("'un' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp, "as.is")) {
        if (!is.character(disp)) stop("'disp' must be a character vector.")
        allowable.disp <- c("means", "sds", all_STATS(p.ops$type))
        if (any(disp %nin% allowable.disp)) {
            stop(paste(word_list(disp[disp %nin% allowable.disp], and.or = "and", quotes = 2, is.are = TRUE),
                       "not allowed in 'disp'."), call. = FALSE)
        }
        if (any(disp %nin% p.ops$compute)) {
            warning(paste("'disp' cannot include", word_list(disp[disp %nin% p.ops$compute], and.or = "or", quotes = 2), "if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
        }
        else p.ops$disp <- disp
    }
    if (is_not_null(A[["disp.means"]]) && !identical(A[["disp.means"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.means"]])) stop("'disp.means' must be TRUE, FALSE, or \"as.is\".")
        if ("means" %nin% p.ops$compute && A[["disp.means"]] == TRUE) {
            warning("'disp.means' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "means"[A[["disp.means"]]]))
    }
    if (is_not_null(A[["disp.sds"]]) && !identical(A[["disp.sds"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.sds"]])) stop("'disp.sds' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if ("sds" %nin% p.ops$compute && A[["disp.sds"]] == TRUE) {
            warning("'disp.sds' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "sds"[A[["disp.sds"]]]))
    }
    if (!identical(stats, "as.is")) {
        if (!is_(stats, "character")) stop("'stats' must be a string.")
        stats <- match_arg(stats, all_STATS(p.ops$type), several.ok = TRUE)
        stats_in_p.ops <- stats %in% p.ops$compute
        if (any(!stats_in_p.ops)) {
            stop(paste0("'stats' cannot contain ", word_list(stats[!stats_in_p.ops], and.or = "or", quotes = 2), " when ", 
                        if (sum(!stats_in_p.ops) > 1) "they were " else "it was ", 
                        "not requested in the original call to bal.tab()."), call. = TRUE)
        }
        else p.ops$disp <- unique(c(p.ops$disp[p.ops$disp %nin% all_STATS()], stats))
    }
    for (s in all_STATS(p.ops$type)) {
        if (is_not_null(A[[STATS[[s]]$disp_stat]]) && !identical(A[[STATS[[s]]$disp_stat]], "as.is")) {
            if (!rlang::is_bool(A[[STATS[[s]]$disp_stat]])) {
                stop(paste0("'", STATS[[s]]$disp_stat, "' must be TRUE, FALSE, or \"as.is\"."), call. = FALSE)
            }
            if (s %nin% p.ops$compute && isTRUE(A[[STATS[[s]]$disp_stat]])) {
                warning(paste0("'", STATS[[s]]$disp_stat, "' cannot be set to TRUE if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            }
            else p.ops$disp <- unique(c(p.ops$disp, s))
        }
    }
    
    for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
        if (STATS[[s]]$threshold %in% names(A) && !identical(temp.thresh <- A[[STATS[[s]]$threshold]], "as.is")) {
            if (is_not_null(temp.thresh) &&
                (!is.numeric(temp.thresh) || length(temp.thresh) != 1 ||
                 is_null(p.ops[["thresholds"]][[s]]) ||
                 p.ops[["thresholds"]][[s]] != temp.thresh))
                stop(paste0("'", STATS[[s]]$threshold, "' must be NULL or \"as.is\"."))
            if (is_null(temp.thresh)) {
                p.ops[["thresholds"]][[s]] <- NULL
                baltal[[s]] <- NULL
                maximbal[[s]] <- NULL
            }
        }
        if (s %nin% p.ops$disp) {
            p.ops[["thresholds"]][[s]] <- NULL
            baltal[[s]] <- NULL
            maximbal[[s]] <- NULL
        }
    }
    if (!identical(disp.thresholds, "as.is")) {
        if (!is.logical(disp.thresholds) || anyNA(disp.thresholds)) stop("'disp.thresholds' must only contain TRUE or FALSE.", call. = FALSE)
        if (is_null(names(disp.thresholds))) {
            if (length(disp.thresholds) <= length(p.ops[["thresholds"]])) {
                names(disp.thresholds) <- names(p.ops[["thresholds"]])[seq_along(disp.thresholds)]
            }
            else {
                stop("More entries were given to 'disp.thresholds' than there are thresholds in the bal.tab object.", call. = FALSE)
            }
        }
        
        if (!all(names(disp.thresholds) %pin% names(p.ops[["thresholds"]]))) {
            warning(paste0(word_list(names(disp.thresholds)[!names(disp.thresholds) %pin% names(p.ops[["thresholds"]])],
                                     quotes = 2, is.are = TRUE), " not available in thresholds and will be ignored."), call. = FALSE)
            disp.thresholds <- disp.thresholds[names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]
        }
        names(disp.thresholds) <- match_arg(names(disp.thresholds), names(p.ops[["thresholds"]]), several.ok = TRUE)
        for (x in names(disp.thresholds)) {
            if (!disp.thresholds[x]) {
                p.ops[["thresholds"]][[x]] <- NULL
                baltal[[x]] <- NULL
                maximbal[[x]] <- NULL
            }
        }
    }
    
    if (!identical(disp.bal.tab, "as.is")) {
        if (!rlang::is_bool(disp.bal.tab)) stop("'disp.bal.tab' must be TRUE, FALSE, or \"as.is\".")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!rlang::is_bool(imbalanced.only)) stop("'imbalanced.only' must be TRUE, FALSE, or \"as.is\".")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (is_null(p.ops$thresholds)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse1(call), collapse = "\n") %+% "\n\n")
    }
    
    if (p.ops$disp.bal.tab) {
        if (p.ops$imbalanced.only) {
            keep.row <- rowSums(apply(balance[grepl(".Threshold", names(balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        }
        else keep.row <- rep(TRUE, nrow(balance))
        
        keep.col <- setNames(as.logical(c(TRUE, 
                                          rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                              p.ops$un && s %in% p.ops$disp
                                          })), switch(p.ops$type, bin = 2, cont = 1)),
                                          unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                              c(p.ops$un && s %in% p.ops$disp,
                                                p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]]))
                                          })),
                                          rep(c(rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                              p.ops$disp.adj && s %in% p.ops$disp
                                          })), switch(p.ops$type, bin = 2, cont = 1)),
                                          unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                              c(p.ops$disp.adj && s %in% p.ops$disp,
                                                p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]]))
                                          }))
                                          ), 
                                          p.ops$nweights + !p.ops$disp.adj))),
                             names(balance))
        
        cat(underline("Balance Measures") %+% "\n")
        if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
        else print.data.frame_(round_df_char(balance[keep.row, keep.col, drop = FALSE], digits))
        cat("\n")
    }
    
    for (s in p.ops$compute) {
        if (is_not_null(baltal[[s]])) {
            cat(underline(paste("Balance tally for", STATS[[s]]$balance_tally_for)) %+% "\n")
            print.data.frame_(baltal[[s]])
            cat("\n")
        }
        if (is_not_null(maximbal[[s]])) {
            cat(underline(paste("Variable with the greatest", STATS[[s]]$variable_with_the_greatest)) %+% "\n")
            print.data.frame_(round_df_char(maximbal[[s]], digits), row.names = FALSE)
            cat("\n")
        }
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
print.bal.tab.cluster <- function(x, imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", stats = "as.is", disp.thresholds = "as.is", disp = "as.is", which.cluster, cluster.summary = "as.is", cluster.fun = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval.parent(.call))
    }
    
    A <- list(...)
    call <- x$call
    c.balance <- x$Cluster.Balance
    c.balance.summary <- x$Balance.Across.Clusters
    nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
    baltal <- maximbal <- list()
    for (s in p.ops$stats) {
        baltal[[s]] <- x[[paste.("Balanced", s)]]
        maximbal[[s]] <- x[[paste.("Max.Imbalance", s)]]
    }
    
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!rlang::is_bool(un)) stop("'un' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("'un' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp, "as.is")) {
        if (!is.character(disp)) stop("'disp.means' must be a character vector.")
        allowable.disp <- c("means", "sds", all_STATS(p.ops$type))
        if (any(disp %nin% allowable.disp)) {
            stop(paste(word_list(disp[disp %nin% allowable.disp], and.or = "and", quotes = 2, is.are = TRUE),
                       "not allowed in 'disp'."), call. = FALSE)
        }
        if (any(disp %nin% p.ops$compute)) {
            warning(paste("'disp' cannot include", word_list(disp[disp %nin% p.ops$compute], and.or = "or", quotes = 2), "if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
        }
        else p.ops$disp <- disp
    }
    if (is_not_null(A[["disp.means"]]) && !identical(A[["disp.means"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.means"]])) stop("'disp.means' must be TRUE, FALSE, or \"as.is\".")
        if ("means" %nin% p.ops$compute && A[["disp.means"]] == TRUE) {
            warning("'disp.means' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "means"[A[["disp.means"]]]))
        
        A[["disp.means"]] <- NULL
    }
    if (is_not_null(A[["disp.sds"]]) && !identical(A[["disp.sds"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.sds"]])) stop("'disp.sds' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if ("sds" %nin% p.ops$compute && A[["disp.sds"]] == TRUE) {
            warning("'disp.sds' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "sds"[A[["disp.sds"]]]))
        
        A[["disp.sds"]] <- NULL
    }
    if (!identical(stats, "as.is")) {
        if (!is_(stats, "character")) stop("'stats' must be a string.")
        stats <- match_arg(stats, all_STATS(p.ops$type), several.ok = TRUE)
        stats_in_p.ops <- stats %in% p.ops$compute
        if (any(!stats_in_p.ops)) {
            stop(paste0("'stats' cannot contain ", word_list(stats[!stats_in_p.ops], and.or = "or", quotes = 2), " if quick = TRUE in the original call to bal.tab()."), call. = TRUE)
        }
        else p.ops$disp <- unique(c(p.ops$disp[p.ops$disp %nin% all_STATS()], stats))
    }
    for (s in all_STATS(p.ops$type)) {
        if (is_not_null(A[[STATS[[s]]$disp_stat]]) && !identical(A[[STATS[[s]]$disp_stat]], "as.is")) {
            if (!rlang::is_bool(A[[STATS[[s]]$disp_stat]])) {
                stop(paste0("'", STATS[[s]]$disp_stat, "' must be TRUE, FALSE, or \"as.is\"."), call. = FALSE)
            }
            if (s %nin% p.ops$compute && isTRUE(A[[STATS[[s]]$disp_stat]])) {
                warning(paste0("'", STATS[[s]]$disp_stat, "' cannot be set to TRUE if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            }
            else p.ops$disp <- unique(c(p.ops$disp, s))
            
            A[[STATS[[s]]$disp_stat]] <- NULL
        }
    }
    
    for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
        if (STATS[[s]]$threshold %in% names(A) && !identical(temp.thresh <- A[[STATS[[s]]$threshold]], "as.is")) {
            if (is_not_null(temp.thresh) &&
                (!is.numeric(temp.thresh) || length(temp.thresh) != 1 ||
                 is_null(p.ops[["thresholds"]][[s]]) ||
                 p.ops[["thresholds"]][[s]] != temp.thresh))
                stop(paste0("'", STATS[[s]]$threshold, "' must be NULL or \"as.is\"."))
            if (is_null(temp.thresh)) {
                p.ops[["thresholds"]][[s]] <- NULL
                baltal[[s]] <- NULL
                maximbal[[s]] <- NULL
            }
            A[[STATS[[s]]$threshold]] <- NULL
        }
        if (s %nin% p.ops$disp) {
            p.ops[["thresholds"]][[s]] <- NULL
            baltal[[s]] <- NULL
            maximbal[[s]] <- NULL
        }
    }
    if (!identical(disp.thresholds, "as.is")) {
        if (!is.logical(disp.thresholds) || anyNA(disp.thresholds)) stop("'disp.thresholds' must only contain TRUE or FALSE.", call. = FALSE)
        if (is_null(names(disp.thresholds))) {
            if (length(disp.thresholds) <= length(p.ops[["thresholds"]])) {
                names(disp.thresholds) <- names(p.ops[["thresholds"]])[seq_along(disp.thresholds)]
            }
            else {
                stop("More entries were given to 'disp.thresholds' than there are thresholds in the bal.tab object.", call. = FALSE)
            }
        }
        
        if (!all(names(disp.thresholds) %pin% names(p.ops[["thresholds"]]))) {
            warning(paste0(word_list(names(disp.thresholds)[!names(disp.thresholds) %pin% names(p.ops[["thresholds"]])],
                                     quotes = 2, is.are = TRUE), " not available in thresholds and will be ignored."), call. = FALSE)
            disp.thresholds <- disp.thresholds[names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]
        }
        names(disp.thresholds) <- match_arg(names(disp.thresholds), names(p.ops[["thresholds"]]), several.ok = TRUE)
        for (x in names(disp.thresholds)) {
            if (!disp.thresholds[x]) {
                p.ops[["thresholds"]][[x]] <- NULL
                baltal[[x]] <- NULL
                maximbal[[x]] <- NULL
            }
        }
    }
    if (!identical(cluster.summary, "as.is")) {
        if (!rlang::is_bool(cluster.summary)) stop("'cluster.summary' must be TRUE, FALSE, or \"as.is\".")
        if (p.ops$quick && p.ops$cluster.summary == FALSE && cluster.summary == TRUE) {
            warning("'cluster.summary' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$cluster.summary <- cluster.summary
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!rlang::is_bool(imbalanced.only)) stop("'imbalanced.only' must be TRUE, FALSE, or \"as.is\".")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (is_null(p.ops$thresholds)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (!missing(which.cluster)) {
        if (paste(deparse1(substitute(which.cluster)), collapse = "") == ".none") which.cluster <- NA
        else if (paste(deparse1(substitute(which.cluster)), collapse = "") == ".all") which.cluster <- NULL
        if (!identical(which.cluster, "as.is")) {
            p.ops$which.cluster <- which.cluster
        }
    }
    
    if (!p.ops$quick || is_null(p.ops$cluster.fun)) computed.cluster.funs <- c("min", "mean", "max")
    else computed.cluster.funs <- p.ops$cluster.fun
    if (is_not_null(cluster.fun) && !identical(cluster.fun, "as.is")) {
        if (!is.character(cluster.fun) || !all(cluster.fun %pin% computed.cluster.funs)) stop(paste0("'cluster.fun' must be ", word_list(c(computed.cluster.funs, "as.is"), and.or = "or", quotes = 2)), call. = FALSE)
    }
    else {
        if (p.ops$abs) cluster.fun <- c("mean", "max")
        else cluster.fun <- c("min", "mean", "max")
    }
    cluster.fun <- match_arg(tolower(cluster.fun), computed.cluster.funs, several.ok = TRUE)
    
    #Checks and Adjustments
    if (is_null(p.ops$which.cluster)) 
        which.cluster <- seq_along(c.balance)
    else if (anyNA(p.ops$which.cluster)) {
        which.cluster <- integer(0)
    }
    else if (is.numeric(p.ops$which.cluster)) {
        which.cluster <- intersect(seq_along(c.balance), p.ops$which.cluster)
        if (is_null(which.cluster)) {
            warning("No indices in 'which.cluster' are cluster indices. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(c.balance)
        }
    }
    else if (is.character(p.ops$which.cluster)) {
        which.cluster <- intersect(names(c.balance), p.ops$which.cluster)
        if (is_null(which.cluster)) {
            warning("No names in 'which.cluster' are cluster names. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(c.balance)
        }
    }
    else {
        warning("The argument to 'which.cluster' must be .all, .none, or a vector of cluster indices or cluster names. Displaying all clusters instead.", call. = FALSE)
        which.cluster <- seq_along(c.balance)
    }
    
    #Printing
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse1(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(which.cluster)) {
        cat(underline("Balance by cluster") %+% "\n")
        for (i in which.cluster) {
            
            cat("\n - - - " %+% italic("Cluster: " %+% names(c.balance)[i]) %+% " - - - \n")
            do.call(print, c(list(c.balance[[i]]), p.ops[names(p.ops) %nin% names(A)], A), quote = TRUE)
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Cluster: ", names(c.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$cluster.summary)) && is_not_null(c.balance.summary)) {
        s.keep.col <- as.logical(c(TRUE, 
                                   unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)], function(s) {
                                       c(unlist(lapply(computed.cluster.funs, function(af) {
                                           p.ops$un && s %in% p.ops$disp && af %in% cluster.fun
                                       })), 
                                       p.ops$un && !p.ops$disp.adj && length(cluster.fun) == 1 && is_not_null(p.ops$thresholds[[s]]))
                                   })),
                                   rep(
                                       unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)], function(s) {
                                           c(unlist(lapply(computed.cluster.funs, function(af) {
                                               p.ops$disp.adj && s %in% p.ops$disp && af %in% cluster.fun
                                           })), 
                                           p.ops$disp.adj && length(cluster.fun) == 1 && is_not_null(p.ops$thresholds[[s]]))
                                       })),
                                       p.ops$nweights + !p.ops$disp.adj)
        ))
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all clusters") %+% "\n")
            print.data.frame_(round_df_char(c.balance.summary[, s.keep.col, drop = FALSE], digits))
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
print.bal.tab.imp <- function(x, imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", stats = "as.is", disp.thresholds = "as.is", disp = "as.is", which.imp, imp.summary = "as.is", imp.fun = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval.parent(.call))
    }
    
    A <- list(...) 
    call <- x$call
    i.balance <- x[["Imputation.Balance"]]
    i.balance.summary <- x[["Balance.Across.Imputations"]]
    nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
    baltal <- maximbal <- list()
    for (s in p.ops$stats) {
        baltal[[s]] <- x[[paste.("Balanced", s)]]
        maximbal[[s]] <- x[[paste.("Max.Imbalance", s)]]
    }
    
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!rlang::is_bool(un)) stop("'un' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("'un' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp, "as.is")) {
        if (!is.character(disp)) stop("'disp.means' must be a character vector.")
        allowable.disp <- c("means", "sds", all_STATS(p.ops$type))
        if (any(disp %nin% allowable.disp)) {
            stop(paste(word_list(disp[disp %nin% allowable.disp], and.or = "and", quotes = 2, is.are = TRUE),
                       "not allowed in 'disp'."), call. = FALSE)
        }
        if (any(disp %nin% p.ops$compute)) {
            warning(paste("'disp' cannot include", word_list(disp[disp %nin% p.ops$compute], and.or = "or", quotes = 2), "if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
        }
        else p.ops$disp <- disp
    }
    if (is_not_null(A[["disp.means"]]) && !identical(A[["disp.means"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.means"]])) stop("'disp.means' must be TRUE, FALSE, or \"as.is\".")
        if ("means" %nin% p.ops$compute && A[["disp.means"]] == TRUE) {
            warning("'disp.means' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "means"[A[["disp.means"]]]))
        A[["disp.means"]] <- NULL
    }
    if (is_not_null(A[["disp.sds"]]) && !identical(A[["disp.sds"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.sds"]])) stop("'disp.sds' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if ("sds" %nin% p.ops$compute && A[["disp.sds"]] == TRUE) {
            warning("'disp.sds' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "sds"[A[["disp.sds"]]]))
        A[["disp.sds"]] <- NULL
    }
    if (!identical(stats, "as.is")) {
        if (!is_(stats, "character")) stop("'stats' must be a string.")
        stats <- match_arg(stats, all_STATS(p.ops$type), several.ok = TRUE)
        stats_in_p.ops <- stats %in% p.ops$compute
        if (any(!stats_in_p.ops)) {
            stop(paste0("'stats' cannot contain ", word_list(stats[!stats_in_p.ops], and.or = "or", quotes = 2), " if quick = TRUE in the original call to bal.tab()."), call. = TRUE)
        }
        else p.ops$disp <- unique(c(p.ops$disp[p.ops$disp %nin% all_STATS()], stats))
    }
    for (s in all_STATS(p.ops$type)) {
        if (is_not_null(A[[STATS[[s]]$disp_stat]]) && !identical(A[[STATS[[s]]$disp_stat]], "as.is")) {
            if (!rlang::is_bool(A[[STATS[[s]]$disp_stat]])) {
                stop(paste0("'", STATS[[s]]$disp_stat, "' must be TRUE, FALSE, or \"as.is\"."), call. = FALSE)
            }
            if (s %nin% p.ops$compute && isTRUE(A[[STATS[[s]]$disp_stat]])) {
                warning(paste0("'", STATS[[s]]$disp_stat, "' cannot be set to TRUE if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            }
            else p.ops$disp <- unique(c(p.ops$disp, s))
            
            A[[STATS[[s]]$disp_stat]] <- NULL
        }
    }
    
    for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
        if (STATS[[s]]$threshold %in% names(A) && !identical(temp.thresh <- A[[STATS[[s]]$threshold]], "as.is")) {
            if (is_not_null(temp.thresh) &&
                (!is.numeric(temp.thresh) || length(temp.thresh) != 1 ||
                 is_null(p.ops[["thresholds"]][[s]]) ||
                 p.ops[["thresholds"]][[s]] != temp.thresh))
                stop(paste0("'", STATS[[s]]$threshold, "' must be NULL or \"as.is\"."))
            if (is_null(temp.thresh)) {
                p.ops[["thresholds"]][[s]] <- NULL
                baltal[[s]] <- NULL
                maximbal[[s]] <- NULL
            }
            A[[STATS[[s]]$threshold]] <- NULL
        }
        if (s %nin% p.ops$disp) {
            p.ops[["thresholds"]][[s]] <- NULL
            baltal[[s]] <- NULL
            maximbal[[s]] <- NULL
        }
    }
    if (!identical(disp.thresholds, "as.is")) {
        if (!is.logical(disp.thresholds) || anyNA(disp.thresholds)) stop("'disp.thresholds' must only contain TRUE or FALSE.", call. = FALSE)
        if (is_null(names(disp.thresholds))) {
            if (length(disp.thresholds) <= length(p.ops[["thresholds"]])) {
                names(disp.thresholds) <- names(p.ops[["thresholds"]])[seq_along(disp.thresholds)]
            }
            else {
                stop("More entries were given to 'disp.thresholds' than there are thresholds in the bal.tab object.", call. = FALSE)
            }
        }
        
        if (!all(names(disp.thresholds) %pin% names(p.ops[["thresholds"]]))) {
            warning(paste0(word_list(names(disp.thresholds)[!names(disp.thresholds) %pin% names(p.ops[["thresholds"]])],
                                     quotes = 2, is.are = TRUE), " not available in thresholds and will be ignored."), call. = FALSE)
            disp.thresholds <- disp.thresholds[names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]
        }
        names(disp.thresholds) <- match_arg(names(disp.thresholds), names(p.ops[["thresholds"]]), several.ok = TRUE)
        for (x in names(disp.thresholds)) {
            if (!disp.thresholds[x]) {
                p.ops[["thresholds"]][[x]] <- NULL
                baltal[[x]] <- NULL
                maximbal[[x]] <- NULL
            }
        }
    }
    if (!identical(imp.summary, "as.is")) {
        if (!rlang::is_bool(imp.summary)) stop("'imp.summary' must be TRUE, FALSE, or \"as.is\".")
        if (p.ops$quick && p.ops$imp.summary == FALSE && imp.summary == TRUE) {
            warning("'imp.summary' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$imp.summary <- imp.summary
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!rlang::is_bool(disp.bal.tab)) stop("'disp.bal.tab' must be TRUE, FALSE, or \"as.is\".")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!rlang::is_bool(imbalanced.only)) stop("'imbalanced.only' must be TRUE, FALSE, or \"as.is\".")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (is_null(p.ops$thresholds)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (!missing(which.imp)) {
        if (paste(deparse1(substitute(which.imp)), collapse = "") == ".none") which.imp <- NA
        else if (paste(deparse1(substitute(which.imp)), collapse = "") == ".all") which.imp <- NULL
        if (!identical(which.imp, "as.is")) {
            p.ops$which.imp <- which.imp
        }
    }
    
    if (!p.ops$quick || is_null(p.ops$imp.fun)) computed.imp.funs <- c("min", "mean", "max")
    else computed.imp.funs <- p.ops$imp.fun
    if (is_not_null(imp.fun) && !identical(imp.fun, "as.is")) {
        if (!is.character(imp.fun) || !all(imp.fun %pin% computed.imp.funs)) stop(paste0("'imp.fun' must be ", word_list(c(computed.imp.funs, "as.is"), and.or = "or", quotes = 2)), call. = FALSE)
    }
    else {
        if (p.ops$abs) imp.fun <- c("mean", "max")
        else imp.fun <- c("min", "mean", "max")
    }
    imp.fun <- match_arg(tolower(imp.fun), computed.imp.funs, several.ok = TRUE)
    
    #Checks and Adjustments
    if (is_null(p.ops$which.imp)) 
        which.imp <- seq_along(i.balance)
    else if (anyNA(p.ops$which.imp)) {
        which.imp <- integer(0)
    }
    else if (is.numeric(p.ops$which.imp)) {
        which.imp <- intersect(seq_along(i.balance), p.ops$which.imp)
        if (is_null(which.imp)) {
            warning("No numbers in 'which.imp' are imputation numbers. No imputations will be displayed.", call. = FALSE)
            which.imp <- integer(0)
        }
    }
    else {
        warning("The argument to 'which.imp' must be .all, .none, or a vector of imputation numbers.", call. = FALSE)
        which.imp <- integer(0)
    }
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse1(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(which.imp)) {
        cat(underline("Balance by imputation") %+% "\n")
        for (i in which.imp) {
            cat("\n - - - " %+% italic("Imputation " %+% names(i.balance)[i]) %+% " - - - \n")
            do.call(print, c(list(i.balance[[i]]), p.ops, A), quote = TRUE)
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Imputation: ", names(i.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$imp.summary)) && is_not_null(i.balance.summary)) {
        s.keep.col <- as.logical(c(TRUE, 
                                   unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)], function(s) {
                                       c(unlist(lapply(computed.imp.funs, function(af) {
                                           p.ops$un && s %in% p.ops$disp && af %in% imp.fun
                                       })), 
                                       p.ops$un && !p.ops$disp.adj && length(imp.fun) == 1 && is_not_null(p.ops$thresholds[[s]]))
                                   })),
                                   rep(
                                       unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)], function(s) {
                                           c(unlist(lapply(computed.imp.funs, function(af) {
                                               p.ops$disp.adj && s %in% p.ops$disp && af %in% imp.fun
                                           })), 
                                           p.ops$disp.adj && length(imp.fun) == 1 && is_not_null(p.ops$thresholds[[s]]))
                                       })),
                                       p.ops$nweights + !p.ops$disp.adj)
        ))
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all imputations") %+% "\n")
            print.data.frame_(round_df_char(i.balance.summary[, s.keep.col, drop = FALSE], digits))
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
print.bal.tab.multi <- function(x, imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", stats = "as.is", disp.thresholds = "as.is", disp = "as.is", which.treat, multi.summary = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval.parent(.call))
    }
    
    A <- list(...)    
    call <- x$call
    m.balance <- x[["Pair.Balance"]]
    m.balance.summary <- x[["Balance.Across.Pairs"]]
    nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
    baltal <- maximbal <- list()
    for (s in p.ops$stats) {
        baltal[[s]] <- x[[paste.("Balanced", s)]]
        maximbal[[s]] <- x[[paste.("Max.Imbalance", s)]]
    }
    
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!rlang::is_bool(un)) stop("'un' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("'un' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp, "as.is")) {
        if (!is.character(disp)) stop("'disp.means' must be a character vector.")
        allowable.disp <- c("means", "sds", all_STATS(p.ops$type))
        if (any(disp %nin% allowable.disp)) {
            stop(paste(word_list(disp[disp %nin% allowable.disp], and.or = "and", quotes = 2, is.are = TRUE),
                       "not allowed in 'disp'."), call. = FALSE)
        }
        if (any(disp %nin% p.ops$compute)) {
            warning(paste("'disp' cannot include", word_list(disp[disp %nin% p.ops$compute], and.or = "or", quotes = 2), "if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
        }
        else p.ops$disp <- disp
    }
    if (is_not_null(A[["disp.means"]]) && !identical(A[["disp.means"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.means"]])) stop("'disp.means' must be TRUE, FALSE, or \"as.is\".")
        if ("means" %nin% p.ops$compute && A[["disp.means"]] == TRUE) {
            warning("'disp.means' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "means"[A[["disp.means"]]]))
    }
    if (is_not_null(A[["disp.sds"]]) && !identical(A[["disp.sds"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.sds"]])) stop("'disp.sds' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if ("sds" %nin% p.ops$compute && A[["disp.sds"]] == TRUE) {
            warning("'disp.sds' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "sds"[A[["disp.sds"]]]))
    }
    if (!identical(stats, "as.is")) {
        if (!is_(stats, "character")) stop("'stats' must be a string.")
        stats <- match_arg(stats, all_STATS(p.ops$type), several.ok = TRUE)
        stats_in_p.ops <- stats %in% p.ops$compute
        if (any(!stats_in_p.ops)) {
            stop(paste0("'stats' cannot contain ", word_list(stats[!stats_in_p.ops], and.or = "or", quotes = 2), " if quick = TRUE in the original call to bal.tab()."), call. = TRUE)
        }
        else p.ops$disp <- unique(c(p.ops$disp[p.ops$disp %nin% all_STATS()], stats))
    }
    for (s in all_STATS(p.ops$type)) {
        if (is_not_null(A[[STATS[[s]]$disp_stat]]) && !identical(A[[STATS[[s]]$disp_stat]], "as.is")) {
            if (!rlang::is_bool(A[[STATS[[s]]$disp_stat]])) {
                stop(paste0("'", STATS[[s]]$disp_stat, "' must be TRUE, FALSE, or \"as.is\"."), call. = FALSE)
            }
            if (s %nin% p.ops$compute && isTRUE(A[[STATS[[s]]$disp_stat]])) {
                warning(paste0("'", STATS[[s]]$disp_stat, "' cannot be set to TRUE if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            }
            else p.ops$disp <- unique(c(p.ops$disp, s))
        }
    }
    
    for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
        if (STATS[[s]]$threshold %in% names(A) && !identical(temp.thresh <- A[[STATS[[s]]$threshold]], "as.is")) {
            if (is_not_null(temp.thresh) &&
                (!is.numeric(temp.thresh) || length(temp.thresh) != 1 ||
                 is_null(p.ops[["thresholds"]][[s]]) ||
                 p.ops[["thresholds"]][[s]] != temp.thresh))
                stop(paste0("'", STATS[[s]]$threshold, "' must be NULL or \"as.is\"."))
            if (is_null(temp.thresh)) {
                p.ops[["thresholds"]][[s]] <- NULL
                baltal[[s]] <- NULL
                maximbal[[s]] <- NULL
            }
        }
        if (s %nin% p.ops$disp) {
            p.ops[["thresholds"]][[s]] <- NULL
            baltal[[s]] <- NULL
            maximbal[[s]] <- NULL
        }
    }
    if (!identical(disp.thresholds, "as.is")) {
        if (!is.logical(disp.thresholds) || anyNA(disp.thresholds)) stop("'disp.thresholds' must only contain TRUE or FALSE.", call. = FALSE)
        if (is_null(names(disp.thresholds))) {
            if (length(disp.thresholds) <= length(p.ops[["thresholds"]])) {
                names(disp.thresholds) <- names(p.ops[["thresholds"]])[seq_along(disp.thresholds)]
            }
            else {
                stop("More entries were given to 'disp.thresholds' than there are thresholds in the bal.tab object.", call. = FALSE)
            }
        }
        
        if (!all(names(disp.thresholds) %pin% names(p.ops[["thresholds"]]))) {
            warning(paste0(word_list(names(disp.thresholds)[!names(disp.thresholds) %pin% names(p.ops[["thresholds"]])],
                                     quotes = 2, is.are = TRUE), " not available in thresholds and will be ignored."), call. = FALSE)
            disp.thresholds <- disp.thresholds[names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]
        }
        names(disp.thresholds) <- match_arg(names(disp.thresholds), names(p.ops[["thresholds"]]), several.ok = TRUE)
        for (x in names(disp.thresholds)) {
            if (!disp.thresholds[x]) {
                p.ops[["thresholds"]][[x]] <- NULL
                baltal[[x]] <- NULL
                maximbal[[x]] <- NULL
            }
        }
    }
    if (!identical(multi.summary, "as.is")) {
        if (!rlang::is_bool(multi.summary)) stop("'multi.summary' must be TRUE, FALSE, or \"as.is\".")
        if (p.ops$quick && p.ops$multi.summary == FALSE && multi.summary == TRUE) {
            warning("'multi.summary' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$multi.summary <- multi.summary
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!rlang::is_bool(disp.bal.tab)) stop("'disp.bal.tab' must be TRUE, FALSE, or \"as.is\".")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (is_not_null(m.balance.summary)) {
        if (p.ops$disp.bal.tab) {
            if (!identical(imbalanced.only, "as.is")) {
                if (!rlang::is_bool(imbalanced.only)) stop("'imbalanced.only' must be TRUE, FALSE, or \"as.is\".")
                p.ops$imbalanced.only <- imbalanced.only
            }
            if (p.ops$imbalanced.only) {
                if (is_null(p.ops$thresholds)) {
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
        if (paste(deparse1(substitute(which.treat)), collapse = "") == ".none") which.treat <- NA
        else if (paste(deparse1(substitute(which.treat)), collapse = "") == ".all") which.treat <- NULL
        if (!identical(which.treat, "as.is")) {
            p.ops$which.treat <- which.treat
        }
    }
    
    #Checks and Adjustments
    if (is_null(p.ops$which.treat)) 
        which.treat <- p.ops$treat_names_multi
    else if (anyNA(p.ops$which.treat)) {
        which.treat <- character(0)
    }
    else if (is.numeric(p.ops$which.treat)) {
        which.treat <- p.ops$treat_names_multi[seq_along(p.ops$treat_names_multi) %in% p.ops$which.treat]
        if (is_null(which.treat)) {
            warning("No numbers in 'which.treat' correspond to treatment values. No treatment pairs will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
    }
    else if (is.character(p.ops$which.treat)) {
        which.treat <- p.ops$treat_names_multi[p.ops$treat_names_multi %in% p.ops$which.treat]
        if (is_null(which.treat)) {
            warning("No names in 'which.treat' correspond to treatment values. No treatment pairs will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
    }
    else {
        warning("The argument to 'which.treat' must be .all, .none, or a vector of treatment names or indices. No treatment pairs will be displayed.", call. = FALSE)
        which.treat <- character(0)
    }
    
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
                    any(treat_names[treat_names != "All"] == which.treat)})]
            }
            else {
                disp.treat.pairs <- names(m.balance)[sapply(names(m.balance), function(x) {
                    treat_names <- attr(m.balance[[x]], "print.options")$treat_names
                    all(treat_names[treat_names != "All"] %in% which.treat)})]
            }
        }
    }
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse1(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(disp.treat.pairs)) {
        headings <- setNames(character(length(disp.treat.pairs)), disp.treat.pairs)
        if (p.ops$pairwise) cat(underline("Balance by treatment pair") %+% "\n")
        else cat(underline("Balance by treatment group") %+% "\n")
        for (i in disp.treat.pairs) {
            headings[i] <- "\n - - - " %+% italic(attr(m.balance[[i]], "print.options")$treat_names[1] %+% " (0) vs. " %+%
                                                      attr(m.balance[[i]], "print.options")$treat_names[2] %+% " (1)") %+% " - - - \n"
            cat(headings[i])
            do.call(print, c(list(m.balance[[i]]), p.ops[names(p.ops) %nin% names(A)], A), quote = TRUE)
        }
        cat(paste0(paste(rep(" -", round(max(nchar(headings))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$multi.summary)) && is_not_null(m.balance.summary)) {
        computed.agg.funs <- "max"
        s.keep.col <- as.logical(c(TRUE, 
                                   unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS("bin")], function(s) {
                                       c(unlist(lapply(computed.agg.funs, function(af) {
                                           p.ops$un && s %in% p.ops$disp && af %in% "max"
                                       })), 
                                       p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]]))
                                   })),
                                   rep(
                                       unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS("bin")], function(s) {
                                           c(unlist(lapply(computed.agg.funs, function(af) {
                                               p.ops$disp.adj && s %in% p.ops$disp && af %in% "max"
                                           })), 
                                           p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]]))
                                       })),
                                       p.ops$nweights + !p.ops$disp.adj)
        ))
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all treatment pairs") %+% "\n")
            if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(m.balance.summary[keep.row, s.keep.col, drop = FALSE], digits))
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
print.bal.tab.msm <- function(x, imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", stats = "as.is", disp.thresholds = "as.is", disp = "as.is", which.time, msm.summary = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval.parent(.call))
    }
    
    A <- list(...)
    A <- clear_null(A[!vapply(A, function(x) identical(x, quote(expr =)), logical(1L))])
    
    call <- x$call
    msm.balance <- x[["Time.Balance"]]
    msm.balance.summary <- x[["Balance.Across.Times"]]
    nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
    baltal <- maximbal <- list()
    for (s in p.ops$stats) {
        baltal[[s]] <- x[[paste.("Balanced", s)]]
        maximbal[[s]] <- x[[paste.("Max.Imbalance", s)]]
    }
    
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!rlang::is_bool(un)) stop("'un' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("'un' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp, "as.is")) {
        if (!is.character(disp)) stop("'disp.means' must be a character vector.")
        allowable.disp <- c("means", "sds", all_STATS(p.ops$type))
        if (any(disp %nin% allowable.disp)) {
            stop(paste(word_list(disp[disp %nin% allowable.disp], and.or = "and", quotes = 2, is.are = TRUE),
                       "not allowed in 'disp'."), call. = FALSE)
        }
        if (any(disp %nin% p.ops$compute)) {
            warning(paste("'disp' cannot include", word_list(disp[disp %nin% p.ops$compute], and.or = "or", quotes = 2), "if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
        }
        else p.ops$disp <- disp
    }
    if (is_not_null(A[["disp.means"]]) && !identical(A[["disp.means"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.means"]])) stop("'disp.means' must be TRUE, FALSE, or \"as.is\".")
        if ("means" %nin% p.ops$compute && A[["disp.means"]] == TRUE) {
            warning("'disp.means' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "means"[A[["disp.means"]]]))
    }
    if (is_not_null(A[["disp.sds"]]) && !identical(A[["disp.sds"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.sds"]])) stop("'disp.sds' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if ("sds" %nin% p.ops$compute && A[["disp.sds"]] == TRUE) {
            warning("'disp.sds' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "sds"[A[["disp.sds"]]]))
    }
    if (!identical(stats, "as.is")) {
        if (!is_(stats, "character")) stop("'stats' must be a string.")
        stats <- match_arg(stats, all_STATS(p.ops$type), several.ok = TRUE)
        stats_in_p.ops <- stats %in% p.ops$compute
        if (any(!stats_in_p.ops)) {
            stop(paste0("'stats' cannot contain ", word_list(stats[!stats_in_p.ops], and.or = "or", quotes = 2), " if quick = TRUE in the original call to bal.tab()."), call. = TRUE)
        }
        else p.ops$disp <- unique(c(p.ops$disp[p.ops$disp %nin% all_STATS()], stats))
    }
    for (s in all_STATS(p.ops$type)) {
        if (is_not_null(A[[STATS[[s]]$disp_stat]]) && !identical(A[[STATS[[s]]$disp_stat]], "as.is")) {
            if (!rlang::is_bool(A[[STATS[[s]]$disp_stat]])) {
                stop(paste0("'", STATS[[s]]$disp_stat, "' must be TRUE, FALSE, or \"as.is\"."), call. = FALSE)
            }
            if (s %nin% p.ops$compute && isTRUE(A[[STATS[[s]]$disp_stat]])) {
                warning(paste0("'", STATS[[s]]$disp_stat, "' cannot be set to TRUE if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            }
            else p.ops$disp <- unique(c(p.ops$disp, s))
        }
    }
    
    for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
        if (STATS[[s]]$threshold %in% names(A) && !identical(temp.thresh <- A[[STATS[[s]]$threshold]], "as.is")) {
            if (is_not_null(temp.thresh) &&
                (!is.numeric(temp.thresh) || length(temp.thresh) != 1 ||
                 is_null(p.ops[["thresholds"]][[s]]) ||
                 p.ops[["thresholds"]][[s]] != temp.thresh))
                stop(paste0("'", STATS[[s]]$threshold, "' must be NULL or \"as.is\"."))
            if (is_null(temp.thresh)) {
                p.ops[["thresholds"]][[s]] <- NULL
                baltal[[s]] <- NULL
                maximbal[[s]] <- NULL
            }
        }
        if (s %nin% p.ops$disp) {
            p.ops[["thresholds"]][[s]] <- NULL
            baltal[[s]] <- NULL
            maximbal[[s]] <- NULL
        }
    }
    if (!identical(disp.thresholds, "as.is")) {
        if (!is.logical(disp.thresholds) || anyNA(disp.thresholds)) stop("'disp.thresholds' must only contain TRUE or FALSE.", call. = FALSE)
        if (is_null(names(disp.thresholds))) {
            if (length(disp.thresholds) <= length(p.ops[["thresholds"]])) {
                names(disp.thresholds) <- names(p.ops[["thresholds"]])[seq_along(disp.thresholds)]
            }
            else {
                stop("More entries were given to 'disp.thresholds' than there are thresholds in the bal.tab object.", call. = FALSE)
            }
        }
        
        if (!all(names(disp.thresholds) %pin% names(p.ops[["thresholds"]]))) {
            warning(paste0(word_list(names(disp.thresholds)[!names(disp.thresholds) %pin% names(p.ops[["thresholds"]])],
                                     quotes = 2, is.are = TRUE), " not available in thresholds and will be ignored."), call. = FALSE)
            disp.thresholds <- disp.thresholds[names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]
        }
        names(disp.thresholds) <- match_arg(names(disp.thresholds), names(p.ops[["thresholds"]]), several.ok = TRUE)
        for (x in names(disp.thresholds)) {
            if (!disp.thresholds[x]) {
                p.ops[["thresholds"]][[x]] <- NULL
                baltal[[x]] <- NULL
                maximbal[[x]] <- NULL
            }
        }
    }
    if (!identical(msm.summary, "as.is")) {
        if (!rlang::is_bool(msm.summary)) stop("'msm.summary' must be TRUE, FALSE, or \"as.is\".")
        if (p.ops$quick && p.ops$msm.summary == FALSE && msm.summary == TRUE) {
            warning("'msm.summary' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$msm.summary <- msm.summary
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!rlang::is_bool(disp.bal.tab)) stop("'disp.bal.tab' must be TRUE, FALSE, or \"as.is\".")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!rlang::is_bool(imbalanced.only)) stop("'imbalanced.only' must be TRUE, FALSE, or \"as.is\".")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (is_null(p.ops$thresholds)) {
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
        if (paste(deparse1(substitute(which.time)), collapse = "") == ".none") which.time <- NA
        else if (paste(deparse1(substitute(which.time)), collapse = "") == ".all") which.time <- NULL
        if (!identical(which.time, "as.is")) {
            p.ops$which.time <- which.time
        }
    }
    
    #Checks and Adjustments
    if (is_null(p.ops$which.time)) 
        which.time <- seq_along(msm.balance)
    else if (anyNA(p.ops$which.time)) {
        which.time <- integer(0)
    }
    else if (is.numeric(p.ops$which.time)) {
        which.time <- seq_along(msm.balance)[seq_along(msm.balance) %in% p.ops$which.time]
        if (is_null(which.time)) {
            warning("No numbers in 'which.time' are treatment time points. No time points will be displayed.", call. = FALSE)
            which.time <- integer(0)
        }
    }
    else if (is.character(p.ops$which.time)) {
        which.time <- seq_along(msm.balance)[names(msm.balance) %in% p.ops$which.time]
        if (is_null(which.time)) {
            warning("No names in 'which.time' are treatment names. No time points will be displayed.", call. = FALSE)
            which.time <- integer(0)
        }
    }
    else {
        warning("The argument to 'which.time' must be .all, .none, or a vector of time point numbers. No time points will be displayed.", call. = FALSE)
        which.time <- integer(0)
    }
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse1(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(which.time)) {
        cat(underline("Balance by Time Point") %+% "\n")
        for (i in which.time) {
            cat("\n - - - " %+% italic("Time: " %+% as.character(i)) %+% " - - - \n")
            do.call(print, c(list(x = msm.balance[[i]]), p.ops[names(p.ops) %nin% names(A)], A), quote = TRUE)
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Time: ", i, " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$msm.summary)) && is_not_null(msm.balance.summary)) {
        computed.agg.funs <- "max"
        s.keep.col <- as.logical(c(TRUE, 
                                   TRUE,
                                   unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)], function(s) {
                                       c(unlist(lapply(computed.agg.funs, function(af) {
                                           p.ops$un && s %in% p.ops$disp && af %in% "max"
                                       })), 
                                       p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]]))
                                   })),
                                   rep(
                                       unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)], function(s) {
                                           c(unlist(lapply(computed.agg.funs, function(af) {
                                               p.ops$disp.adj && s %in% p.ops$disp && af %in% "max"
                                           })), 
                                           p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]]))
                                       })),
                                       p.ops$nweights + !p.ops$disp.adj)
        ))
        
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all time points") %+% "\n")
            if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(msm.balance.summary[keep.row, s.keep.col, drop = FALSE], digits))
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

print.bal.tab.subclass <- function(x, imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", stats = "as.is", disp.thresholds = "as.is", disp = "as.is", disp.subclass = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    A <- list(...)
    call <- x$call
    s.balance <- x$Subclass.Balance
    b.a.subclass <- x$Balance.Across.Subclass
    s.nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
    baltal <- maximbal <- list()
    for (s in p.ops$compute) {
        baltal[[s]] <- x[[paste.("Balanced", s, "Subclass")]]
        maximbal[[s]] <- x[[paste.("Max.Imbalance", s, "Subclass")]]
    }
    
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!rlang::is_bool(un)) stop("'un' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("'un' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp, "as.is")) {
        if (!is.character(disp)) stop("'disp.means' must be a character vector.")
        allowable.disp <- c("means", "sds", all_STATS(p.ops$type))
        if (any(disp %nin% allowable.disp)) {
            stop(paste(word_list(disp[disp %nin% allowable.disp], and.or = "and", quotes = 2, is.are = TRUE),
                       "not allowed in 'disp'."), call. = FALSE)
        }
        if (any(disp %nin% p.ops$compute)) {
            warning(paste("'disp' cannot include", word_list(disp[disp %nin% p.ops$compute], and.or = "or", quotes = 2), "if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
        }
        else p.ops$disp <- disp
    }
    if (is_not_null(A[["disp.means"]]) && !identical(A[["disp.means"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.means"]])) stop("'disp.means' must be TRUE, FALSE, or \"as.is\".")
        if ("means" %nin% p.ops$compute && A[["disp.means"]] == TRUE) {
            warning("'disp.means' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "means"[A[["disp.means"]]]))
    }
    if (is_not_null(A[["disp.sds"]]) && !identical(A[["disp.sds"]], "as.is")) {
        if (!rlang::is_bool(A[["disp.sds"]])) stop("'disp.sds' must be TRUE, FALSE, or \"as.is\".", call. = FALSE)
        if ("sds" %nin% p.ops$compute && A[["disp.sds"]] == TRUE) {
            warning("'disp.sds' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "sds"[A[["disp.sds"]]]))
    }
    if (!identical(stats, "as.is")) {
        if (!is_(stats, "character")) stop("'stats' must be a string.")
        stats <- match_arg(stats, all_STATS(p.ops$type), several.ok = TRUE)
        stats_in_p.ops <- stats %in% p.ops$compute
        if (any(!stats_in_p.ops)) {
            stop(paste0("'stats' cannot contain ", word_list(stats[!stats_in_p.ops], and.or = "or", quotes = 2), " if quick = TRUE in the original call to bal.tab()."), call. = TRUE)
        }
        else p.ops$disp <- unique(c(p.ops$disp[p.ops$disp %nin% all_STATS()], stats))
    }
    for (s in all_STATS(p.ops$type)) {
        if (is_not_null(A[[STATS[[s]]$disp_stat]]) && !identical(A[[STATS[[s]]$disp_stat]], "as.is")) {
            if (!rlang::is_bool(A[[STATS[[s]]$disp_stat]])) {
                stop(paste0("'", STATS[[s]]$disp_stat, "' must be TRUE, FALSE, or \"as.is\"."), call. = FALSE)
            }
            if (s %nin% p.ops$compute && isTRUE(A[[STATS[[s]]$disp_stat]])) {
                warning(paste0("'", STATS[[s]]$disp_stat, "' cannot be set to TRUE if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            }
            else p.ops$disp <- unique(c(p.ops$disp, s))
        }
    }
    
    for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
        if (STATS[[s]]$threshold %in% names(A) && !identical(temp.thresh <- A[[STATS[[s]]$threshold]], "as.is")) {
            if (is_not_null(temp.thresh) &&
                (!is.numeric(temp.thresh) || length(temp.thresh) != 1 ||
                 is_null(p.ops[["thresholds"]][[s]]) ||
                 p.ops[["thresholds"]][[s]] != temp.thresh))
                stop(paste0("'", STATS[[s]]$threshold, "' must be NULL or \"as.is\"."))
            if (is_null(temp.thresh)) {
                p.ops[["thresholds"]][[s]] <- NULL
                baltal[[s]] <- NULL
                maximbal[[s]] <- NULL
            }
        }
        if (s %nin% p.ops$disp) {
            p.ops[["thresholds"]][[s]] <- NULL
            baltal[[s]] <- NULL
            maximbal[[s]] <- NULL
        }
    }
    if (!identical(disp.thresholds, "as.is")) {
        if (!is.logical(disp.thresholds) || anyNA(disp.thresholds)) stop("'disp.thresholds' must only contain TRUE or FALSE.", call. = FALSE)
        if (is_null(names(disp.thresholds))) {
            if (length(disp.thresholds) <= length(p.ops[["thresholds"]])) {
                names(disp.thresholds) <- names(p.ops[["thresholds"]])[seq_along(disp.thresholds)]
            }
            else {
                stop("More entries were given to 'disp.thresholds' than there are thresholds in the bal.tab object.", call. = FALSE)
            }
        }
        
        if (!all(names(disp.thresholds) %pin% names(p.ops[["thresholds"]]))) {
            warning(paste0(word_list(names(disp.thresholds)[!names(disp.thresholds) %pin% names(p.ops[["thresholds"]])],
                                     quotes = 2, is.are = TRUE), " not available in thresholds and will be ignored."), call. = FALSE)
            disp.thresholds <- disp.thresholds[names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]
        }
        names(disp.thresholds) <- match_arg(names(disp.thresholds), names(p.ops[["thresholds"]]), several.ok = TRUE)
        for (x in names(disp.thresholds)) {
            if (!disp.thresholds[x]) {
                p.ops[["thresholds"]][[x]] <- NULL
                baltal[[x]] <- NULL
                maximbal[[x]] <- NULL
            }
        }
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!rlang::is_bool(disp.bal.tab)) stop("'disp.bal.tab' must be TRUE, FALSE, or \"as.is\".")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!rlang::is_bool(imbalanced.only)) stop("'imbalanced.only' must be TRUE, FALSE, or \"as.is\".")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (is_null(p.ops$thresholds)) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    if (!identical(disp.subclass, "as.is")) {
        if (!rlang::is_bool(disp.subclass)) stop("'disp.subclass' must be TRUE, FALSE, or \"as.is\".")
        p.ops$disp.subclass <- disp.subclass
    }
    
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse1(call), collapse = "\n") %+% "\n\n")
    }
    
    if (p.ops$disp.bal.tab) {
        if (p.ops$disp.subclass) {
            s.keep.col <- setNames(c(TRUE,
                                     rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                         s %in% p.ops$disp
                                     })), switch(p.ops$type, bin = 2, cont = 1)),
                                     unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                         c(s %in% p.ops$disp,
                                           is_not_null(p.ops$thresholds[[s]]))
                                     }))),
                                   names(s.balance[[1]]))
            
            cat(underline("Balance by subclass"))
            for (i in names(s.balance)) {
                if (p.ops$imbalanced.only) {
                    s.keep.row <- rowSums(apply(s.balance[[i]][grepl(".Threshold", names(s.balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
                }
                else s.keep.row <- rep(TRUE, nrow(s.balance[[i]]))
                
                cat("\n - - - " %+% italic("Subclass " %+% as.character(i)) %+% " - - - \n")
                if (all(!s.keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
                else print.data.frame_(round_df_char(s.balance[[i]][s.keep.row, s.keep.col, drop = FALSE], digits))
            }
            cat("\n")
        }
        
        if (is_not_null(b.a.subclass)) {
            if (p.ops$imbalanced.only) {
                a.s.keep.row <- rowSums(apply(b.a.subclass[grepl(".Threshold", names(b.a.subclass), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
            }
            else a.s.keep.row <- rep(TRUE, nrow(b.a.subclass))
            
            a.s.keep.col <- setNames(as.logical(c(TRUE, 
                                                  rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                                      p.ops$un && s %in% p.ops$disp
                                                  })), switch(p.ops$type, bin = 2, cont = 1)),
                                                  unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                                      c(p.ops$un && s %in% p.ops$disp,
                                                        p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]]))
                                                  })),
                                                  rep(c(rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                                      p.ops$disp.adj && s %in% p.ops$disp
                                                  })), 2),
                                                  unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                                      c(p.ops$disp.adj && s %in% p.ops$disp,
                                                        p.ops$disp.adj && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]]))
                                                  }))
                                                  ), 
                                                  p.ops$disp.adj))),
                                     names(b.a.subclass))
            
            cat(underline("Balance measures across subclasses") %+% "\n")
            if (all(!a.s.keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(b.a.subclass[a.s.keep.row, a.s.keep.col, drop = FALSE], digits))
            cat("\n")
        }
    }
    
    for (s in p.ops$stats) {
        if (is_not_null(baltal[[s]])) {
            cat(underline(paste("Balance tally for", STATS[[s]]$balance_tally_for, "across subclasses")) %+% "\n")
            print.data.frame_(baltal[[s]])
            cat("\n")
        }
        if (is_not_null(maximbal[[s]])) {
            cat(underline(paste("Variable with the greatest", STATS[[s]]$variable_with_the_greatest, "across subclasses")) %+% "\n")
            print.data.frame_(round_df_char(maximbal[[s]], digits), row.names = FALSE)
            cat("\n")
        }
    }
    
    if (is_not_null(s.nn)) {
        cat(underline(attr(s.nn, "tag")) %+% "\n")
        print.data.frame_(round_df_char(s.nn, digits = max(0, digits-1)))
    }
    
    invisible(x)
}
