.base.bal.tab.binary <- function(X, int = FALSE, poly = 1, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), thresholds = list(), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), abs = FALSE, quick = TRUE, ...) {
    #Preparations
    A <- clear_null(list(...))
    
    X$treat <- process_treat(X$treat) 
    
    type = "bin"
    
    if (get.treat.type(X$treat) != "binary") {
        stop("Treatment indicator must be a binary variable---e.g., treatment (1) or control (0)", call. = FALSE)
    }
    
    # thresholds <- list()
    
    for (s in all_STATS[get_from_STATS("type") == type]) {
        #If disp.stat is TRUE, add stat to stats
        if (isTRUE(A[[STATS[[s]]$disp_stat]])) {
            X$stats <- unique(c(X$stats, s))
        }
        
        #Process and check thresholds
        if (is_not_null(temp.thresh <- get0(STATS[[s]]$threshold))) {
            thresholds[[s]] <- STATS[[s]]$abs(temp.thresh)
            if (!between(thresholds[[s]], STATS[[s]]$threshold_range)) {
                thresholds[[s]] <- NULL
                warning(paste0(STATS[[s]]$threshold, " must be between ", word_list(STATS[[s]]$threshold_range),
                               "; ignoring ", STATS[[s]]$threshold, "."), call. = FALSE)
            }
            else X$stats <- unique(c(X$stats, s))
        }
    }
    
    # if (isTRUE(A[["disp.diff"]])) X$stats <- unique(c(X$stats, "mean.diffs"))
    # if (isTRUE(A[["disp.v.ratio"]])) X$stats <- unique(c(X$stats, "variance.ratios"))
    # if (isTRUE(A[["disp.ks"]])) X$stats <- unique(c(X$stats, "ks.statistics"))
    
    # if (is_not_null(m.threshold)) {
    #     m.threshold <- abs(m.threshold)
    #     X$stats <- unique(c(X$stats, "mean.diffs"))
    # }
    # if (is_not_null(v.threshold)) {
    #     v.threshold <- abs_(v.threshold, ratio = TRUE)
    #     X$stats <- unique(c(X$stats, "variance.ratios"))
    # }
    # if (is_null(ks.threshold) && is_null(A[["k.threshold"]])) {
    #     ks.threshold <- A[["k.threshold"]]
    # }
    # if (is_not_null(ks.threshold)) {
    #     if (ks.threshold > 1) {
    #         warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
    #         ks.threshold <- NULL
    #     }
    #     else X$stats <- unique(c(X$stats, "ks.statistics"))
    # }
    
    if (is_null(X$weights)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        check_if_zero_weights(X$weights, X$treat)
        if(ncol(X$weights) == 1) names(X$weights) <- "Adj"
    }
    if (is_null(X$s.weights)) {
        X$s.weights <- rep(1, length(X$treat))
    }
    
    #Actions
    out.names <- c("Balance", 
                   expand.grid_string(c("Balanced", "Max.Imbalance"),
                                      X$stats,
                                      # c("Means", "Variances", "KS"),
                                      collapse = "."), 
                   "Observations", "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    C <- get.C(covs = X$covs, addl = X$addl, distance = X$distance, int = int, poly = poly, ...)
    co.names <- attr(C, "co.names")
    
    out[["Balance"]] <- balance.table.bin(C, weights = X$weights, treat = X$treat, 
                                          s.d.denom = X$s.d.denom, s.weights = X$s.weights, 
                                          continuous, binary, 
                                          thresholds = thresholds,
                                          # m.threshold = thresholds[["mean.diffs"]], v.threshold = thresholds[["variance.ratios"]], ks.threshold = thresholds[["ks.statistics"]], 
                                          # m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold,
                                          un = un, disp.means = disp.means, 
                                          disp.sds = disp.sds, stats = X$stats, 
                                          abs = abs, no.adj = no.adj, quick = quick, s.d.denom.list = X$s.d.denom.list)
    
    #Reassign disp... and ...threshold based on balance table output
    stats <- attr(out[["Balance"]], "stats")
    thresholds <- attr(out[["Balance"]], "thresholds")
    
    for (s in stats) {
        if (is_not_null(thresholds[[s]])) {
            if (no.adj) {
                out[[paste.("Balanced", s)]] <- baltal(out[["Balance"]][[paste.(STATS[[s]]$Threshold, "Un")]])
                out[[paste.("Max.Imbalance", s)]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], 
                                                               paste.(STATS[[s]]$bal.tab._column_prefix, "Un"), 
                                                               paste.(STATS[[s]]$Threshold, "Un"), 
                                                               ratio = s == "variance.ratios")
            }
            else if (ncol(X$weights) == 1) {
                out[[paste.("Balanced", s)]] <- baltal(out[["Balance"]][[STATS[[s]]$Threshold]])
                out[[paste.("Max.Imbalance", s)]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], 
                                                               paste.(STATS[[s]]$bal.tab._column_prefix, "Adj"), 
                                                               STATS[[s]]$Threshold, 
                                                               ratio = s == "variance.ratios")
            }
            else if (ncol(X$weights) > 1) {
                out[[paste.("Balanced", s)]] <- setNames(do.call("cbind", lapply(names(X$weights), function(x) baltal(out[["Balance"]][[paste.(STATS[[s]]$Threshold, x)]]))),
                                                         names(X$weights))
                out[[paste.("Max.Imbalance", s)]] <- cbind(Weights = names(X$weights),
                                                           do.call("rbind", lapply(names(X$weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], 
                                                                                                                                    paste.(STATS[[s]]$bal.tab._column_prefix, x), 
                                                                                                                                    paste.(STATS[[s]]$Threshold, x), 
                                                                                                                                    ratio = s == "variance.ratios"),
                                                                                                                          c("Variable", 
                                                                                                                            STATS[[s]]$bal.tab._column_prefix, 
                                                                                                                            STATS[[s]]$Threshold)))),
                                                           stringsAsFactors = FALSE)
            }
        }
        else {
            out[[paste.("Balanced", s)]] <- NULL
            out[[paste.("Max.Imbalance", s)]] <- NULL
        }
    }
    
    # for (i in names(attr(out[["Balance"]], "threshold"))) {
    #     assign(paste0(i, ".threshold"), attr(out[["Balance"]], "threshold")[i])
    # }
    # 
    # S <- list(diff = list(threshold = m.threshold,
    #                       Names = "Means",
    #                       Threshold = "M.Threshold",
    #                       Stat = "Diff"),
    #           v.ratio = list(threshold = v.threshold,
    #                          Names = "Variances",
    #                          Threshold = "V.Threshold",
    #                          Stat = "V.Ratio"),
    #           ks = list(threshold = ks.threshold,
    #                     Names = "KS",
    #                     Threshold = "KS.Threshold",
    #                     Stat = "KS"))
    # 
    # for (s in S) {
    #     if (is_not_null(s[["threshold"]])) {
    #         if (no.adj) {
    #             out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[paste.(s[["Threshold"]], "Un")]])
    #             out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Un"), paste.(s[["Threshold"]], "Un"), ratio = s$Stat == "V.Ratio")
    #         }
    #         else if (ncol(X$weights) == 1) {
    #             out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[s[["Threshold"]]]])
    #             out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio")
    #         }
    #         else if (ncol(X$weights) > 1) {
    #             out[[paste.("Balanced", s[["Names"]])]] <- setNames(do.call("cbind", lapply(names(X$weights), function(x) baltal(out[["Balance"]][[paste.(s[["Threshold"]], x)]]))),
    #                                                                 names(X$weights))
    #             out[[paste.("Max.Imbalance", s[["Names"]])]] <- cbind(Weights = names(X$weights),
    #                                                                   do.call("rbind", lapply(names(X$weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], x), paste.(s[["Threshold"]], x), ratio = s$Stat == "V.Ratio"),
    #                                                                                                                                  c("Variable", s[["Stat"]], s[["Threshold"]])))),
    #                                                                   stringsAsFactors = FALSE)
    #         }
    #     }
    #     else {
    #         out[[paste.("Balanced", s[["Names"]])]] <- NULL
    #         out[[paste.("Max.Imbalance", s[["Names"]])]] <- NULL
    #     }
    # }
    
    out[["Observations"]] <- samplesize(treat = X$treat, weights = X$weights, s.weights = X$s.weights, method = X$method, discarded = X$discarded)
    out[["call"]] <- X$call
    attr(out, "print.options") <- list(thresholds = thresholds,
                                       # m.threshold=thresholds[["mean.diffs"]], 
                                       # v.threshold=thresholds[["variance.ratios"]], 
                                       # ks.threshold=thresholds[["ks.statistics"]], 
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.means=disp.means, 
                                       disp.sds=disp.sds,
                                       # disp.diff=disp.diff, 
                                       # disp.v.ratio=disp.v.ratio, 
                                       # disp.ks=disp.ks, 
                                       stats = stats, 
                                       disp.adj=!no.adj,
                                       disp.bal.tab = disp.bal.tab, 
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       nweights = ifelse(no.adj, 0, ncol(X$weights)),
                                       weight.names = names(X$weights),
                                       treat_names = treat_vals(X$treat),
                                       type = type,
                                       co.names = co.names)
    class(out) <- c("bal.tab.bin", "bal.tab")
    
    return(out)
}
.base.bal.tab.cont <- function(X, int = FALSE, poly = 1, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "std"), thresholds = list(), r.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), abs = FALSE, quick = TRUE, ...) {
    
    #Preparations
    A <- clear_null(list(...))
    
    X$treat <- process_treat(X$treat) 
    
    type <- "cont"
    
    for (s in all_STATS[get_from_STATS("type") == type]) {
        #If disp.stat is TRUE, add stat to stats
        if (isTRUE(A[[STATS[[s]]$disp_stat]])) {
            X$stats <- unique(c(X$stats, s))
        }
        
        #Process and check thresholds
        if (is_not_null(temp.thresh <- get0(STATS[[s]]$threshold))) {
            thresholds[[s]] <- STATS[[s]]$abs(temp.thresh)
            if (!between(thresholds[[s]], STATS[[s]]$threshold_range)) {
                thresholds[[s]] <- NULL
                warning(paste0(STATS[[s]]$threshold, " must be between ", word_list(STATS[[s]]$threshold_range),
                               "; ignoring ", STATS[[s]]$threshold, "."), call. = FALSE)
            }
            else X$stats <- unique(c(X$stats, s))
        }
    }
    
    # if (isTRUE(A[["disp.corr"]])) X$stats <- unique(c(X$stats, "correlations"))
    
    # if (is_not_null(r.threshold)) {
    #     r.threshold <- abs(r.threshold)
    #     if (r.threshold > 1) {
    #         warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
    #         r.threshold <- NULL
    #     }
    #     else X$stats <- unique(c(X$stats, "correlations"))
    # }
    if (is_null(X$weights)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        check_if_zero_weights(X$weights)
        if (ncol(X$weights) == 1) names(X$weights) <- "Adj"
    }
    if (is_null(X$s.weights)) {
        X$s.weights <- rep(1, length(X$treat))
    }    
    
    #Actions
    
    out.names <- c("Balance", 
                   expand.grid_string(c("Balanced", "Max.Imbalance"),
                                      X$stats,
                                      # c("Means", "Variances", "KS"),
                                      collapse = "."), 
                   "Observations", "call")
    
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    C <- get.C(covs = X$covs, int = int, poly = poly, addl = X$addl, distance = X$distance, ...)
    co.names <- attr(C, "co.names")
    
    out[["Balance"]] <- balance.table.cont(C, weights = X$weights, treat = X$treat, 
                                           s.d.denom = X$s.d.denom, s.weights = X$s.weights, 
                                           continuous = continuous, binary = binary, 
                                           # r.threshold = r.threshold, 
                                           # r.threshold = thresholds[["correlations"]],
                                           thresholds = thresholds,
                                           un = un, 
                                           disp.means = disp.means, disp.sds = disp.sds, 
                                           stats = X$stats, abs = abs, 
                                           no.adj = no.adj, quick = quick)
    
    #Reassign disp... and ...threshold based on balance table output
    stats <- attr(out[["Balance"]], "stats")
    thresholds <- attr(out[["Balance"]], "thresholds")
    
    for (s in stats) {
        if (is_not_null(thresholds[[s]])) {
            if (no.adj) {
                out[[paste.("Balanced", s)]] <- baltal(out[["Balance"]][[paste.(STATS[[s]]$Threshold, "Un")]])
                out[[paste.("Max.Imbalance", s)]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], 
                                                               paste.(STATS[[s]]$bal.tab._column_prefix, "Un"), 
                                                               paste.(STATS[[s]]$Threshold, "Un"), 
                                                               ratio = s == "variance.ratios")
            }
            else if (ncol(X$weights) == 1) {
                out[[paste.("Balanced", s)]] <- baltal(out[["Balance"]][[STATS[[s]]$Threshold]])
                out[[paste.("Max.Imbalance", s)]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], 
                                                               paste.(STATS[[s]]$bal.tab._column_prefix, "Adj"), 
                                                               STATS[[s]]$Threshold, 
                                                               ratio = s == "variance.ratios")
            }
            else if (ncol(X$weights) > 1) {
                out[[paste.("Balanced", s)]] <- setNames(do.call("cbind", lapply(names(X$weights), function(x) baltal(out[["Balance"]][[paste.(STATS[[s]]$Threshold, x)]]))),
                                                         names(X$weights))
                out[[paste.("Max.Imbalance", s)]] <- cbind(Weights = names(X$weights),
                                                           do.call("rbind", lapply(names(X$weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], 
                                                                                                                                    paste.(STATS[[s]]$bal.tab._column_prefix, x), 
                                                                                                                                    paste.(STATS[[s]]$Threshold, x), 
                                                                                                                                    ratio = s == "variance.ratios"),
                                                                                                                          c("Variable", 
                                                                                                                            STATS[[s]]$bal.tab._column_prefix, 
                                                                                                                            STATS[[s]]$Threshold)))),
                                                           stringsAsFactors = FALSE)
            }
        }
        else {
            out[[paste.("Balanced", s)]] <- NULL
            out[[paste.("Max.Imbalance", s)]] <- NULL
        }
    }
    
    # for (i in names(attr(out[["Balance"]], "threshold"))) {
    #     assign(paste0(i, ".threshold"), attr(out[["Balance"]], "threshold")[i])
    # }
    # 
    # S <- list(corr = list(threshold = r.threshold,
    #                       Names = "Corr",
    #                       Threshold = "R.Threshold",
    #                       Stat = "Corr"))
    # 
    # for (s in S) {
    #     if (is_not_null(s[["threshold"]])) {
    #         if (no.adj) {
    #             out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[paste.(s[["Threshold"]], "Un")]])
    #             out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Un"), paste.(s[["Threshold"]], "Un"), ratio = s$Stat == "V.Ratio")
    #         }
    #         else if (ncol(X$weights) == 1) {
    #             out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[s[["Threshold"]]]])
    #             out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio")
    #         }
    #         else if (ncol(X$weights) > 1) {
    #             out[[paste.("Balanced", s[["Names"]])]] <- setNames(do.call("cbind", lapply(names(X$weights), function(x) baltal(out[["Balance"]][[paste.(s[["Threshold"]], x)]]))),
    #                                                                 names(X$weights))
    #             out[[paste.("Max.Imbalance", s[["Names"]])]] <- cbind(Weights = names(X$weights),
    #                                                                   do.call("rbind", lapply(names(X$weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], x), paste.(s[["Threshold"]], x), ratio = s$Stat == "V.Ratio"),
    #                                                                                                                                  c("Variable", s[["Stat"]], s[["Threshold"]])))),
    #                                                                   stringsAsFactors = FALSE)
    #         }
    #     }
    # }
    
    out[["Observations"]] <- samplesize.cont(treat = X$treat, weights = X$weights, s.weights = X$s.weights, method = X$method, discarded = X$discarded)
    out[["call"]] <- X$call
    attr(out, "print.options") <- list(thresholds = thresholds,
                                       # r.threshold = thresholds[["correlations"]],
                                       # r.threshold=r.threshold, 
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.means=disp.means, 
                                       disp.sds = disp.sds,
                                       # disp.corr = disp.corr, 
                                       stats = stats,
                                       disp.adj=!no.adj,
                                       disp.bal.tab = disp.bal.tab,
                                       continuous = continuous,
                                       binary = binary,
                                       abs = abs,
                                       quick = quick,
                                       nweights = if (no.adj) 0 else ncol(X$weights),
                                       weight.names = names(X$weights),
                                       type = type,
                                       co.names = co.names)
    class(out) <- c("bal.tab.cont", "bal.tab")
    
    return(out)
}
