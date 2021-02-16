print.bal.tab <- function(x, imbalanced.only, un, disp.bal.tab, disp.call,
                          stats, disp.thresholds, disp, 
                          which.subclass, subclass.summary,
                          which.imp, imp.summary, imp.fun, 
                          which.treat, multi.summary, 
                          which.time, msm.summary,
                          which.cluster, cluster.summary, cluster.fun, 
                          digits = max(3, getOption("digits") - 3), ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval.parent(.call))
    }
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    unpack_p.ops <- function(b) {
        out <- do.call("print_process", c(list(b), args), quote = TRUE)
        
        if (is_(b, c("bal.tab.bin", "bal.tab.cont"))) return(out)
        else {
            b_ <- b[[which(endsWith(names(b), ".Balance"))]][[1]]
            if (is_(b_, "bal.tab")) out <- c(out, unpack_p.ops(b_))
            else return(out)
        }
    }
    
    p.ops <- unpack_p.ops(x)
    
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))

    bal.tab_print(x, p.ops)
}

bal.tab_print <- function(x, p.ops) {
    UseMethod("bal.tab_print")
}
bal.tab_print.bal.tab <- function(x, p.ops) {
    
    call <- if (p.ops$disp.call) x$call else NULL
    
    balance <- x$Balance
    
    thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
    baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
    maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)

    nn <- x$Observations
    
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
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
                                          unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                              c(p.ops$un && s %in% p.ops$disp,
                                                if (p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                          })),
                                          rep(c(rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                              p.ops$disp.adj && s %in% p.ops$disp
                                          })), switch(p.ops$type, bin = 2, cont = 1)),
                                          unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                              c(p.ops$disp.adj && s %in% p.ops$disp,
                                                if (p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                          }))
                                          ), 
                                          p.ops$nweights + !p.ops$disp.adj))),
                             names(balance))

        cat(underline("Balance Measures") %+% "\n")
        if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
        else print.data.frame_(round_df_char(balance[keep.row, keep.col, drop = FALSE], p.ops$digits, na_vals = "."))
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
            print.data.frame_(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
            cat("\n")
        }
    }
    
    if (is_not_null(nn)) {
        
        drop.nn <- rowSums(nn) == 0
        ss.type <- attr(nn, "ss.type")[!drop.nn]
        nn <- nn[!drop.nn, , drop = FALSE]
        
        if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
            all(check_if_zero(nn["Matched (ESS)",] - nn["Matched (Unweighted)",]))) {
            nn <- nn[rownames(nn)!="Matched (Unweighted)", , drop = FALSE]
            rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
        }
        cat(underline(attr(nn, "tag")) %+% "\n")
        print.warning <- FALSE
        if (length(ss.type) > 1 && nunique.gt(ss.type[-1], 1)) {
            ess <- ifelse(ss.type == "ess", "*", "")
            nn <- setNames(cbind(nn, ess), c(names(nn), ""))
            print.warning <- TRUE
        }
        print.data.frame_(round_df_char(nn, digits = min(2, p.ops$digits), pad = " "))
        if (print.warning) cat(italic("* indicates effective sample size"))
    }
    invisible(x)
}
bal.tab_print.bal.tab.cluster <- function(x, p.ops) {
    
    call <- if (p.ops$disp.call) x$call else NULL
    c.balance <- x$Cluster.Balance
    c.balance.summary <- x$Balance.Across.Clusters
    
    thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
    baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
    maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)
    
    nn <- x$Observations
    
    #Printing
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(p.ops$which.cluster)) {
        cat(underline("Balance by cluster") %+% "\n")
        for (i in p.ops$which.cluster) {
            cat("\n - - - " %+% italic("Cluster: " %+% names(c.balance)[i]) %+% " - - - \n")
            bal.tab_print(c.balance[[i]], p.ops)
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Cluster: ", names(c.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$cluster.summary)) && is_not_null(c.balance.summary)) {
        s.keep.col <- setNames(as.logical(c(TRUE, 
                                   unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                       c(unlist(lapply(p.ops$computed.cluster.funs, function(af) {
                                           p.ops$un && s %in% p.ops$disp && af %in% p.ops$cluster.fun
                                       })), 
                                       if (p.ops$un && !p.ops$disp.adj && length(p.ops$cluster.fun) == 1 && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                   })),
                                   rep(
                                       unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                           c(unlist(lapply(p.ops$computed.cluster.funs, function(af) {
                                               p.ops$disp.adj && s %in% p.ops$disp && af %in% p.ops$cluster.fun
                                           })), 
                                           if (p.ops$disp.adj && length(p.ops$cluster.fun) == 1 && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                       })),
                                       p.ops$nweights + !p.ops$disp.adj)
        )), names(c.balance.summary))
    
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all clusters") %+% "\n")
            print.data.frame_(round_df_char(c.balance.summary[, s.keep.col, drop = FALSE], p.ops$digits, na_vals = "."))
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
                print.data.frame_(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
                cat("\n")
            }
        }
        
        if (is_not_null(nn)) {
            drop.nn <- rowSums(nn) == 0
            ss.type <- attr(nn, "ss.type")[!drop.nn]
            nn <- nn[!drop.nn, , drop = FALSE]
            if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
                all(check_if_zero(nn["Matched (ESS)",] - nn["Matched (Unweighted)",]))) {
                nn <- nn[rownames(nn)!="Matched (Unweighted)", , drop = FALSE]
                rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
            }
            cat(underline(attr(nn, "tag")) %+% "\n")
            print.warning <- FALSE
            if (length(ss.type) > 1 && nunique.gt(ss.type[-1], 1)) {
                ess <- ifelse(ss.type == "ess", "*", "")
                nn <- setNames(cbind(nn, ess), c(names(nn), ""))
                print.warning <- TRUE
            }
            print.data.frame_(round_df_char(nn, digits = min(2, p.ops$digits), pad = " "))
            if (print.warning) cat(italic("* indicates effective sample size"))
        }
    }
    
    invisible(x)
}
bal.tab_print.bal.tab.imp <- function(x, p.ops) {
    
    call <- if (p.ops$disp.call) x$call else NULL
    i.balance <- x[["Imputation.Balance"]]
    i.balance.summary <- x[["Balance.Across.Imputations"]]
    
    thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
    baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
    maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)
    
    nn <- x$Observations
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(p.ops$which.imp)) {
        cat(underline("Balance by imputation") %+% "\n")
        for (i in p.ops$which.imp) {
            cat("\n - - - " %+% italic("Imputation " %+% names(i.balance)[i]) %+% " - - - \n")
            bal.tab_print(i.balance[[i]], p.ops)
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Imputation: ", names(i.balance)[i], " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$imp.summary)) && is_not_null(i.balance.summary)) {
        s.keep.col <- as.logical(c(TRUE, 
                                   unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                       c(unlist(lapply(p.ops$computed.imp.funs, function(af) {
                                           p.ops$un && s %in% p.ops$disp && af %in% p.ops$imp.fun
                                       })), 
                                       if (p.ops$un && !p.ops$disp.adj && length(p.ops$imp.fun) == 1 && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                   })),
                                   rep(
                                       unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                           c(unlist(lapply(p.ops$computed.imp.funs, function(af) {
                                               p.ops$disp.adj && s %in% p.ops$disp && af %in% p.ops$imp.fun
                                           })), 
                                           if (p.ops$disp.adj && length(p.ops$imp.fun) == 1 && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                       })),
                                       p.ops$nweights + !p.ops$disp.adj)
        ))
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all imputations") %+% "\n")
            print.data.frame_(round_df_char(i.balance.summary[, s.keep.col, drop = FALSE], p.ops$digits, na_vals = "."))
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
                print.data.frame_(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
                cat("\n")
            }
        }
        
        if (is_not_null(nn)) {
            drop.nn <- rowSums(nn) == 0
            ss.type <- attr(nn, "ss.type")[!drop.nn]
            nn <- nn[!drop.nn, , drop = FALSE]
            if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
                all(check_if_zero(nn["Matched (ESS)",] - nn["Matched (Unweighted)",]))) {
                nn <- nn[rownames(nn)!="Matched (Unweighted)", , drop = FALSE]
                rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
            }
            cat(underline(attr(nn, "tag")) %+% "\n")
            print.warning <- FALSE
            if (length(ss.type) > 1 && nunique.gt(ss.type[-1], 1)) {
                ess <- ifelse(ss.type == "ess", "*", "")
                nn <- setNames(cbind(nn, ess), c(names(nn), ""))
                print.warning <- TRUE
            }
            print.data.frame_(round_df_char(nn, digits = min(2, p.ops$digits), pad = " "))
            if (print.warning) cat(italic("* indicates effective sample size"))
        }
    }
    
    invisible(x)
    
}
bal.tab_print.bal.tab.multi <- function(x, p.ops) {

    call <- if (p.ops$disp.call) x$call else NULL
    m.balance <- x[["Pair.Balance"]]
    m.balance.summary <- x[["Balance.Across.Pairs"]]
    
    thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
    baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
    maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)
    
    nn <- x$Observations
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(p.ops$disp.treat.pairs)) {
        headings <- setNames(character(length(p.ops$disp.treat.pairs)), p.ops$disp.treat.pairs)
        if (p.ops$pairwise) cat(underline("Balance by treatment pair") %+% "\n")
        else cat(underline("Balance by treatment group") %+% "\n")
        for (i in p.ops$disp.treat.pairs) {
            headings[i] <- "\n - - - " %+% italic(attr(m.balance[[i]], "print.options")$treat_names[1] %+% " (0) vs. " %+%
                                                      attr(m.balance[[i]], "print.options")$treat_names[2] %+% " (1)") %+% " - - - \n"
            cat(headings[i])
            bal.tab_print(m.balance[[i]], p.ops)
        }
        cat(paste0(paste(rep(" -", round(max(nchar(headings))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$multi.summary)) && is_not_null(m.balance.summary)) {
        
        if (p.ops$imbalanced.only) {
            keep.row <- rowSums(apply(m.balance.summary[grepl(".Threshold", names(m.balance.summary), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        }
        else keep.row <- rep(TRUE, nrow(m.balance.summary))
        
        computed.agg.funs <- "max"
        s.keep.col <- as.logical(c(TRUE, 
                                   unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS("bin")[!get_from_STATS("adj_only")]], function(s) {
                                       c(unlist(lapply(computed.agg.funs, function(af) {
                                           p.ops$un && s %in% p.ops$disp && af %in% "max"
                                       })), 
                                       if (p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                   })),
                                   rep(
                                       unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS("bin")], function(s) {
                                           c(unlist(lapply(computed.agg.funs, function(af) {
                                               p.ops$disp.adj && s %in% p.ops$disp && af %in% "max"
                                           })), 
                                           if (p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                       })),
                                       p.ops$nweights + !p.ops$disp.adj)
        ))
        names(s.keep.col) <- names(m.balance.summary)
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all treatment pairs") %+% "\n")
            if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(m.balance.summary[keep.row, s.keep.col, drop = FALSE], p.ops$digits, na_vals = "."))
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
                print.data.frame_(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
                cat("\n")
            }
        }
        
        if (is_not_null(nn)) {
            tag <- attr(nn, "tag")
            drop.nn <- rowSums(nn) == 0
            ss.type <- attr(nn, "ss.type")[!drop.nn]
            nn <- nn[!drop.nn, , drop = FALSE]
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
            print.data.frame_(round_df_char(nn, digits = min(2, p.ops$digits), pad = " "))
            if (print.warning) cat(italic("* indicates effective sample size"))
        }
    }
    
    invisible(x)
    
}
bal.tab_print.bal.tab.msm <- function(x, p.ops){
 
    call <- if (p.ops$disp.call) x$call else NULL
    msm.balance <- x[["Time.Balance"]]
    msm.balance.summary <- x[["Balance.Across.Times"]]
    
    thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
    baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
    maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)
    
    nn <- x$Observations
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(p.ops$which.time)) {
        cat(underline("Balance by Time Point") %+% "\n")
        for (i in p.ops$which.time) {
            cat("\n - - - " %+% italic("Time: " %+% as.character(i)) %+% " - - - \n")
            bal.tab_print(msm.balance[[i]], p.ops)
        }
        cat(paste0(paste(rep(" -", round(nchar(paste0("\n - - - Time: ", i, " - - - "))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$msm.summary)) && is_not_null(msm.balance.summary)) {
        
        if (p.ops$imbalanced.only) {
            keep.row <- rowSums(apply(msm.balance.summary[grepl(".Threshold", names(msm.balance.summary), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        }
        else keep.row <- rep(TRUE, nrow(msm.balance.summary))
        
        computed.agg.funs <- "max"
        s.keep.col <- as.logical(c(TRUE, 
                                   TRUE,
                                   unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                       c(unlist(lapply(computed.agg.funs, function(af) {
                                           p.ops$un && s %in% p.ops$disp && af %in% "max"
                                       })), 
                                       if (p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                   })),
                                   rep(
                                       unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                           c(unlist(lapply(computed.agg.funs, function(af) {
                                               p.ops$disp.adj && s %in% p.ops$disp && af %in% "max"
                                           })), 
                                           if (p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                       })),
                                       p.ops$nweights + !p.ops$disp.adj)
        ))
        
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Balance summary across all time points") %+% "\n")
            if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(msm.balance.summary[keep.row, s.keep.col, drop = FALSE], p.ops$digits, na_vals = "."))
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
                print.data.frame_(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
                cat("\n")
            }
        }
        
        if (is_not_null(nn)) {
            print.warning <- FALSE
            cat(underline(attr(nn[[1]], "tag")) %+% "\n")
            
            for (ti in seq_along(nn)) {
                cat(" - " %+% italic("Time " %+% as.character(ti)) %+% "\n")
                drop.nn <- rowSums(nn[[ti]]) == 0
                ss.type <- attr(nn[[ti]], "ss.type")[!drop.nn]
                nn[[ti]] <- nn[[ti]][!drop.nn, , drop = FALSE]
                if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn[[ti]])) && 
                    all(check_if_zero(nn[[ti]]["Matched (ESS)",] - nn[[ti]]["Matched (Unweighted)",]))) {
                    nn[[ti]] <- nn[[ti]][rownames(nn[[ti]])!="Matched (Unweighted)", , drop = FALSE]
                    rownames(nn[[ti]])[rownames(nn[[ti]]) == "Matched (ESS)"] <- "Matched"
                }
                if (length(ss.type) > 1 && nunique.gt(ss.type[-1], 1)) {
                    ess <- ifelse(ss.type == "ess", "*", "")
                    nn[[ti]] <- setNames(cbind(nn[[ti]], ess), c(names(nn[[ti]]), ""))
                    print.warning <- TRUE
                }
                print.data.frame_(round_df_char(nn[[ti]], digits = min(2, p.ops$digits), pad = " "))
            }
            
            if (print.warning) cat(italic("* indicates effective sample size"))
        }
    }
    
    invisible(x)
}

bal.tab_print.bal.tab.subclass <- function(x, p.ops) {
    
    call <- if (p.ops$disp.call) x$call else NULL
    s.balance <- x$Subclass.Balance
    b.a.subclass <- x$Balance.Across.Subclass
    s.nn <- x$Observations

    thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
    baltal <- setNames(x[paste.("Balanced", thresholds, "Subclass")], thresholds)
    maximbal <- setNames(x[paste.("Max.Imbalance", thresholds, "Subclass")], thresholds)
    
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    #Print subclass balance
    if (p.ops$disp.bal.tab) {
        if (is_not_null(p.ops$which.subclass)) {
            s.keep.col <- setNames(c(TRUE,
                                     rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                         s %in% p.ops$disp
                                     })), switch(p.ops$type, bin = 2, cont = 1)),
                                     unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                         c(s %in% p.ops$disp,
                                           if (is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                     }))),
                                   names(s.balance[[1]]))
            
            cat(underline("Balance by subclass"))
            for (i in p.ops$which.subclass) {
                if (p.ops$imbalanced.only) {
                    s.keep.row <- rowSums(apply(s.balance[[i]][grepl(".Threshold", names(s.balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
                }
                else s.keep.row <- rep(TRUE, nrow(s.balance[[i]]))
                
                cat("\n - - - " %+% italic("Subclass " %+% as.character(i)) %+% " - - - \n")
                if (all(!s.keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
                else print.data.frame_(round_df_char(s.balance[[i]][s.keep.row, s.keep.col, drop = FALSE], p.ops$digits, na_vals = "."))
            }
            cat("\n")
        }
    }
    
    #Print balance across subclasses
    if (p.ops$subclass.summary && is_not_null(b.a.subclass)) {
        
        if (p.ops$disp.bal.tab) {
            if (p.ops$imbalanced.only) {
                a.s.keep.row <- rowSums(apply(b.a.subclass[grepl(".Threshold", names(b.a.subclass), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
            }
            else a.s.keep.row <- rep(TRUE, nrow(b.a.subclass))
            
            a.s.keep.col <- setNames(as.logical(c(TRUE, 
                                                  rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                                      p.ops$un && s %in% p.ops$disp
                                                  })), switch(p.ops$type, bin = 2, cont = 1)),
                                                  unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                                      c(p.ops$un && s %in% p.ops$disp,
                                                        if (p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                                  })),
                                                  rep(c(rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                                      p.ops$disp.adj && s %in% p.ops$disp
                                                  })), 2),
                                                  unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                                      c(p.ops$disp.adj && s %in% p.ops$disp,
                                                        if (p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                                  }))
                                                  ), 
                                                  p.ops$disp.adj)
                                                  )),
                                     names(b.a.subclass))
            
            cat(underline("Balance measures across subclasses") %+% "\n")
            if (all(!a.s.keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(b.a.subclass[a.s.keep.row, a.s.keep.col, drop = FALSE], p.ops$digits, na_vals = "."))
            cat("\n")
        }
        
        for (s in p.ops$compute) {
            if (is_not_null(baltal[[s]])) {
                cat(underline(paste("Balance tally for", STATS[[s]]$balance_tally_for, "across subclasses")) %+% "\n")
                print.data.frame_(baltal[[s]])
                cat("\n")
            }
            if (is_not_null(maximbal[[s]])) {
                cat(underline(paste("Variable with the greatest", STATS[[s]]$variable_with_the_greatest, "across subclasses")) %+% "\n")
                print.data.frame_(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
                cat("\n")
            }
        }
        
        if (is_not_null(s.nn)) {
            cat(underline(attr(s.nn, "tag")) %+% "\n")
            print.data.frame_(round_df_char(s.nn, digits = min(2, p.ops$digits), pad = " "))
        }
    }
    
    invisible(x)
}


#Process arguments
print_process <- function(x, ...) {
    UseMethod("print_process")
}
print_process.bal.tab.cluster <- function(x, which.cluster, cluster.summary, cluster.fun, ...) {
    A <- list(...)
    
    c.balance <- x$Cluster.Balance
    p.ops <- attr(x, "print.options")
    
    if (!missing(cluster.summary)) {
        if (!rlang::is_bool(cluster.summary)) stop("'cluster.summary' must be TRUE or FALSE.")
        if (p.ops$quick && p.ops$cluster.summary == FALSE && cluster.summary == TRUE) {
            warning("'cluster.summary' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$cluster.summary <- cluster.summary
    }
    
    if (!missing(which.cluster)) {
        if (paste(deparse1(substitute(which.cluster)), collapse = "") == ".none") which.cluster <- NA
        else if (paste(deparse1(substitute(which.cluster)), collapse = "") == ".all") which.cluster <- NULL
        p.ops$which.cluster <- which.cluster
    }
    
    if (!p.ops$quick || is_null(p.ops$cluster.fun)) computed.cluster.funs <- c("min", "mean", "max")
    else computed.cluster.funs <- p.ops$cluster.fun
    if (!missing(cluster.fun) && is_not_null(cluster.fun)) {
        if (!is.character(cluster.fun) || !all(cluster.fun %pin% computed.cluster.funs)) stop(paste0("'cluster.fun' must be ", word_list(computed.cluster.funs, and.or = "or", quotes = 2)), call. = FALSE)
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
        which.cluster <- seq_along(c.balance)[names(c.balance) %in% p.ops$which.cluster]
        if (is_null(which.cluster)) {
            warning("No names in 'which.cluster' are cluster names. Displaying all clusters instead.", call. = FALSE)
            which.cluster <- seq_along(c.balance)
        }
    }
    else {
        warning("The argument to 'which.cluster' must be .all, .none, or a vector of cluster indices or cluster names. Displaying all clusters instead.", call. = FALSE)
        which.cluster <- seq_along(c.balance)
    }
    
    out <- list(cluster.summary = p.ops$cluster.summary,
                cluster.fun = cluster.fun,
                which.cluster = which.cluster,
                computed.cluster.funs = computed.cluster.funs)
    
    out
}
print_process.bal.tab.imp <- function(x, which.imp, imp.summary, imp.fun, ...) {
    A <- list(...) 
    i.balance <- x[["Imputation.Balance"]]
    p.ops <- attr(x, "print.options")

    if (!missing(imp.summary)) {
        if (!rlang::is_bool(imp.summary)) stop("'imp.summary' must be TRUE or FALSE.")
        if (p.ops$quick && p.ops$imp.summary == FALSE && imp.summary == TRUE) {
            warning("'imp.summary' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$imp.summary <- imp.summary
    }
    
    if (!missing(which.imp)) {
        if (paste(deparse1(substitute(which.imp)), collapse = "") == ".none") which.imp <- NA
        else if (paste(deparse1(substitute(which.imp)), collapse = "") == ".all") which.imp <- NULL
        p.ops$which.imp <- which.imp
    }
    
    if (!p.ops$quick || is_null(p.ops$imp.fun)) computed.imp.funs <- c("min", "mean", "max")
    else computed.imp.funs <- p.ops$imp.fun
    if (!missing(imp.fun) && is_not_null(imp.fun)) {
        if (!is.character(imp.fun) || !all(imp.fun %pin% computed.imp.funs)) stop(paste0("'imp.fun' must be ", word_list(computed.imp.funs, and.or = "or", quotes = 2)), call. = FALSE)
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
    
    out <- list(imp.summary = p.ops$imp.summary,
                imp.fun = imp.fun,
                which.imp = which.imp,
                computed.imp.funs = computed.imp.funs)
    
    out
}
print_process.bal.tab.multi <- function(x, which.treat, multi.summary, ...) {
    A <- list(...)    

    m.balance <- x[["Pair.Balance"]]
    m.balance.summary <- x[["Balance.Across.Pairs"]]
    
    p.ops <- attr(x, "print.options")
    
    if (!missing(multi.summary)) {
        if (!rlang::is_bool(multi.summary)) stop("'multi.summary' must be TRUE or FALSE.")
        if (p.ops$quick && p.ops$multi.summary == FALSE && multi.summary == TRUE) {
            warning("'multi.summary' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$multi.summary <- multi.summary
    }
    
    if (!missing(which.treat)) {
        if (paste(deparse1(substitute(which.treat)), collapse = "") == ".none") which.treat <- NA
        else if (paste(deparse1(substitute(which.treat)), collapse = "") == ".all") which.treat <- NULL
        p.ops$which.treat <- which.treat
    }
    
    #Checks and Adjustments
    if (is_null(p.ops$which.treat)) 
        which.treat <- p.ops$treat_vals_multi
    else if (anyNA(p.ops$which.treat)) {
        which.treat <- character(0)
    }
    else if (!is_(p.ops$which.treat, c("numeric", "character"))) {
        warning("The argument to 'which.treat' must be .all, .none, or a vector of treatment names or indices. No treatment pairs will be displayed.", call. = FALSE)
        which.treat <- character(0)
    }
    else {
        if (length(p.ops$treat_vals_multi) == 2) p.ops$which.treat <- as.character(p.ops$which.treat)
        
        if (is.numeric(p.ops$which.treat)) {
            which.treat <- p.ops$treat_vals_multi[seq_along(p.ops$treat_vals_multi) %in% p.ops$which.treat]
            if (is_null(which.treat)) {
                warning("No numbers in 'which.treat' correspond to treatment values. No treatment pairs will be displayed.", call. = FALSE)
                which.treat <- character(0)
            }
        }
        else if (is.character(p.ops$which.treat)) {
            which.treat <- p.ops$treat_vals_multi[p.ops$treat_vals_multi %in% p.ops$which.treat]
            if (is_null(which.treat)) {
                warning("No names in 'which.treat' correspond to treatment values. No treatment pairs will be displayed.", call. = FALSE)
                which.treat <- character(0)
            }
        }
    }
    
    if (is_null(which.treat)) {
        disp.treat.pairs <- character(0)
    }
    else {
        if (p.ops$pairwise) {
            if (length(which.treat) == 1) {
                disp.treat.pairs <- names(m.balance)[sapply(m.balance, function(z) {
                    treat_names <- attr(z, "print.options")$treat_names
                    any(p.ops$treat_vals_multi[treat_names] == which.treat)
                })]
            }
            else {
                disp.treat.pairs <- names(m.balance)[sapply(m.balance, function(z) {
                    treat_names <- attr(z, "print.options")$treat_names
                    all(p.ops$treat_vals_multi[treat_names] %in% which.treat)
                })]
            }
        }
        else {
            if (length(which.treat) == 1) {
                disp.treat.pairs <- names(m.balance)[sapply(m.balance, function(z) {
                    treat_names <- attr(z, "print.options")$treat_names
                    any(p.ops$treat_vals_multi[treat_names[treat_names != "All"]] == which.treat)
                })]
            }
            else {
                disp.treat.pairs <- names(m.balance)[sapply(m.balance, function(z) {
                    treat_names <- attr(z, "print.options")$treat_names
                    all(p.ops$treat_vals_multi[treat_names[treat_names != "All"]] %in% which.treat)
                })]
            }
        }
    }
    
    out <- list(disp.treat.pairs = disp.treat.pairs,
                multi.summary = p.ops$multi.summary,
                pairwise = p.ops$pairwise)
    
    out
}
print_process.bal.tab.msm <- function(x, which.time, msm.summary, ...) {

    A <- list(...)
    A <- clear_null(A[!vapply(A, function(x) identical(x, quote(expr =)), logical(1L))])
    
    msm.balance <- x[["Time.Balance"]]
    
    p.ops <- attr(x, "print.options")
    
    if (!missing(msm.summary)) {
        if (!rlang::is_bool(msm.summary)) stop("'msm.summary' must be TRUE or FALSE.")
        if (p.ops$quick && p.ops$msm.summary == FALSE && msm.summary == TRUE) {
            warning("'msm.summary' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$msm.summary <- msm.summary
    }
    
    if (!missing(which.time)) {
        if (paste(deparse1(substitute(which.time)), collapse = "") == ".none") which.time <- NA
        else if (paste(deparse1(substitute(which.time)), collapse = "") == ".all") which.time <- NULL
        p.ops$which.time <- which.time
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
    
    out <- list(msm.summary = p.ops$msm.summary,
                which.time = which.time)
    
    out
}
print_process.bal.tab <- function(x, imbalanced.only, un, disp.bal.tab, disp.call, stats, disp.thresholds, disp, digits = max(3, getOption("digits") - 3), ...) {
    
    A <- list(...)
    p.ops <- attr(x, "print.options")

    drop.thresholds <- c()
    
    #Adjustments to print options
    if (!missing(un) && p.ops$disp.adj) {
        if (!rlang::is_bool(un)) stop("'un' must be TRUE or FALSE.", call. = FALSE)
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("'un' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!missing(disp)) {
        if (!rlang::is_character(disp)) stop("'disp' must be a character vector.", call. = FALSE)
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
    if (is_not_null(A[["disp.means"]])) {
        if (!rlang::is_bool(A[["disp.means"]])) stop("'disp.means' must be TRUE or FALSE.")
        if ("means" %nin% p.ops$compute && A[["disp.means"]] == TRUE) {
            warning("'disp.means' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "means"[A[["disp.means"]]]))
    }
    if (is_not_null(A[["disp.sds"]])) {
        if (!rlang::is_bool(A[["disp.sds"]])) stop("'disp.sds' must be TRUE or FALSE.", call. = FALSE)
        if ("sds" %nin% p.ops$compute && A[["disp.sds"]] == TRUE) {
            warning("'disp.sds' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "sds"[A[["disp.sds"]]]))
    }
    if (!missing(stats)) {
        if (!rlang::is_character(stats)) stop("'stats' must be a string.", call. = FALSE)
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
        if (is_not_null(A[[STATS[[s]]$disp_stat]])) {
            if (!rlang::is_bool(A[[STATS[[s]]$disp_stat]])) {
                stop(paste0("'", STATS[[s]]$disp_stat, "' must be TRUE or FALSE."), call. = FALSE)
            }
            if (s %nin% p.ops$compute && isTRUE(A[[STATS[[s]]$disp_stat]])) {
                warning(paste0("'", STATS[[s]]$disp_stat, "' cannot be set to TRUE if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            }
            else p.ops$disp <- unique(c(p.ops$disp, s))
        }
    }
    
    for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
        if (STATS[[s]]$threshold %in% names(A)) {
            temp.thresh <- A[[STATS[[s]]$threshold]]
            if (is_not_null(temp.thresh) &&
                (!is.numeric(temp.thresh) || length(temp.thresh) != 1 ||
                 is_null(p.ops[["thresholds"]][[s]]) ||
                 p.ops[["thresholds"]][[s]] != temp.thresh))
                stop(paste0("'", STATS[[s]]$threshold, "' must be NULL or left unspecified."))
            if (is_null(temp.thresh)) {
                drop.thresholds <- c(drop.thresholds, s)
            }
        }
        if (s %nin% p.ops$disp) {
            drop.thresholds <- c(drop.thresholds, s)
        }
    }
    if (!missing(disp.thresholds)) {
        if (!rlang::is_logical(disp.thresholds) || anyNA(disp.thresholds)) stop("'disp.thresholds' must only contain TRUE or FALSE.", call. = FALSE)
        if (is_null(names(disp.thresholds))) {
            if (length(disp.thresholds) <= length(p.ops[["thresholds"]])) {
                if (length(disp.thresholds) == 1) disp.thresholds <- rep(disp.thresholds, length(p.ops[["thresholds"]]))
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
        for (i in names(disp.thresholds)) {
            if (!disp.thresholds[i]) {
                drop.thresholds <- c(drop.thresholds, i)
            }
        }
    }
    
    if (!missing(disp.bal.tab)) {
        if (!rlang::is_bool(disp.bal.tab)) stop("'disp.bal.tab' must be TRUE or FALSE.")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!missing(imbalanced.only)) {
            if (!rlang::is_bool(imbalanced.only)) stop("'imbalanced.only' must be TRUE or FALSE.")
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
    
    if (!missing(disp.call)) {
        if (!rlang::is_bool(disp.call)) stop("'disp.call' must be TRUE or FALSE.", call. = FALSE)
        if (disp.call && is_null(x$call)) {
            warning("'disp.call' cannot be set to TRUE if the input object does not have a call component.", call. = FALSE)
        }
        else p.ops$disp.call <- disp.call
    }
    
    out <- list(un = p.ops$un,
                disp = p.ops$disp,
                compute = p.ops$compute,
                drop.thresholds = drop.thresholds,
                disp.bal.tab = p.ops$disp.bal.tab,
                imbalanced.only = p.ops$imbalanced.only,
                digits = digits,
                disp.adj = p.ops$disp.adj,
                thresholds = p.ops$thresholds,
                type = p.ops$type,
                nweights = p.ops$nweights,
                disp.call = p.ops$disp.call)
    
    out
    
}
print_process.bal.tab.subclass <- function(x, imbalanced.only, un, disp.bal.tab, disp.call, stats, disp.thresholds, disp, digits = max(3, getOption("digits") - 3), which.subclass, subclass.summary, ...) {
    A <- list(...)
    
    s.balance <- x$Subclass.Balance
    p.ops <- attr(x, "print.options")
    
    drop.thresholds <- c()
    
    #Adjustments to print options
    if (!missing(un) && p.ops$disp.adj) {
        if (!rlang::is_bool(un)) stop("'un' must be TRUE or FALSE.", call. = FALSE)
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("'un' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!missing(disp)) {
        if (!rlang::is_character(disp)) stop("'disp' must be a character vector.", call. = FALSE)
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
    if (is_not_null(A[["disp.means"]])) {
        if (!rlang::is_bool(A[["disp.means"]])) stop("'disp.means' must be TRUE or FALSE.")
        if ("means" %nin% p.ops$compute && A[["disp.means"]] == TRUE) {
            warning("'disp.means' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "means"[A[["disp.means"]]]))
    }
    if (is_not_null(A[["disp.sds"]])) {
        if (!rlang::is_bool(A[["disp.sds"]])) stop("'disp.sds' must be TRUE or FALSE.", call. = FALSE)
        if ("sds" %nin% p.ops$compute && A[["disp.sds"]] == TRUE) {
            warning("'disp.sds' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$disp <- unique(c(p.ops$disp, "sds"[A[["disp.sds"]]]))
    }
    if (!missing(stats)) {
        if (!rlang::is_character(stats)) stop("'stats' must be a string.", call. = FALSE)
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
        if (is_not_null(A[[STATS[[s]]$disp_stat]])) {
            if (!rlang::is_bool(A[[STATS[[s]]$disp_stat]])) {
                stop(paste0("'", STATS[[s]]$disp_stat, "' must be TRUE or FALSE."), call. = FALSE)
            }
            if (s %nin% p.ops$compute && isTRUE(A[[STATS[[s]]$disp_stat]])) {
                warning(paste0("'", STATS[[s]]$disp_stat, "' cannot be set to TRUE if quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            }
            else p.ops$disp <- unique(c(p.ops$disp, s))
        }
    }
    
    for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
        if (STATS[[s]]$threshold %in% names(A)) {
            temp.thresh <- A[[STATS[[s]]$threshold]]
            if (is_not_null(temp.thresh) &&
                (!is.numeric(temp.thresh) || length(temp.thresh) != 1 ||
                 is_null(p.ops[["thresholds"]][[s]]) ||
                 p.ops[["thresholds"]][[s]] != temp.thresh))
                stop(paste0("'", STATS[[s]]$threshold, "' must be NULL or left unspecified."))
            if (is_null(temp.thresh)) {
                drop.thresholds <- c(drop.thresholds, s)
            }
        }
        if (s %nin% p.ops$disp) {
            drop.thresholds <- c(drop.thresholds, s)
        }
    }
    if (!missing(disp.thresholds)) {
        if (!rlang::is_logical(disp.thresholds) || anyNA(disp.thresholds)) stop("'disp.thresholds' must only contain TRUE or FALSE.", call. = FALSE)
        if (is_null(names(disp.thresholds))) {
            if (length(disp.thresholds) <= length(p.ops[["thresholds"]])) {
                if (length(disp.thresholds) == 1) disp.thresholds <- rep(disp.thresholds, length(p.ops[["thresholds"]]))
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
                drop.thresholds <- c(drop.thresholds, x)
            }
        }
    }
    
    if (!missing(disp.bal.tab)) {
        if (!rlang::is_bool(disp.bal.tab)) stop("'disp.bal.tab' must be TRUE or FALSE.")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!missing(imbalanced.only)) {
            if (!rlang::is_bool(imbalanced.only)) stop("'imbalanced.only' must be TRUE or FALSE.")
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
    
    if (!missing(disp.call)) {
        if (!rlang::is_bool(disp.call)) stop("'disp.call' must be TRUE or FALSE.", call. = FALSE)
        if (disp.call && is_null(x$call)) {
            warning("'disp.call' cannot be set to TRUE if the input object does not have a call component.", call. = FALSE)
        }
        else p.ops$disp.call <- disp.call
    }
    
    if (!missing(subclass.summary)) {
        if (!rlang::is_bool(subclass.summary)) stop("'subclass.summary' must be TRUE or FALSE.")
        if (p.ops$quick && p.ops$subclass.summary == FALSE && subclass.summary == TRUE) {
            warning("'subclass.summary' cannot be set to TRUE if quick = TRUE in the original call to bal.tab().", call. = FALSE)
        }
        else p.ops$subclass.summary <- subclass.summary
    }
    
    if (!missing(which.subclass)) {
        p.ops$which.subclass <- which.subclass
    }
    else if (is_not_null(A[["disp.subclass"]])) {
        p.ops$which.subclass <- if (isTRUE(A[["disp.subclass"]])) NULL else NA
    }
    
    #Checks and Adjustments
    if (is_null(p.ops$which.subclass)) 
        which.subclass <- seq_along(s.balance)
    else if (anyNA(p.ops$which.subclass)) {
        which.subclass <- integer(0)
    }
    else if (is.numeric(p.ops$which.subclass)) {
        which.subclass <- intersect(seq_along(s.balance), p.ops$which.subclass)
        if (is_null(which.subclass)) {
            warning("No indices in 'which.subclass' are subclass indices. No subclasses will be displayed.", call. = FALSE)
            which.subclass <- NA
        }
    }
    else {
        warning("The argument to 'which.subclass' must be .all, .none, or a vector of subclass indices. No subclasses will be displayed.", call. = FALSE)
        which.subclass <- NA
    }
    
    out <- list(un = p.ops$un,
                disp = p.ops$disp,
                compute = p.ops$compute,
                drop.thresholds = drop.thresholds,
                disp.bal.tab = p.ops$disp.bal.tab,
                imbalanced.only = p.ops$imbalanced.only,
                digits = digits,
                disp.adj = p.ops$disp.adj,
                thresholds = p.ops$thresholds,
                type = p.ops$type,
                subclass.summary = p.ops$subclass.summary,
                which.subclass = which.subclass,
                disp.call = p.ops$disp.call)
    
    out
}
