base.bal.tab <- function(X, ...) {
    UseMethod("base.bal.tab")
}

base.bal.tab.binary <- function(X, int = FALSE, poly = 1, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), abs = FALSE, quick = TRUE, ...) {
    #Preparations
    A <- clear_null(list(...))
    
    X$treat <- process_treat(X$treat) 
    
    if (get.treat.type(X$treat) != "binary") {
        stop("Treatment indicator must be a binary variable---e.g., treatment (1) or control (0)", call. = FALSE)
    }
    
    if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
    if (is_not_null(v.threshold)) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (is_null(ks.threshold) && is_null(A[["k.threshold"]])) {
        ks.threshold <- A[["k.threshold"]]
    }
    if (is_not_null(ks.threshold)) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
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
                                      c("Means", "Variances", "KS"),
                                      collapse = "."), 
                   "Observations", "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    C <- get.C(covs = X$covs, addl = X$addl, distance = X$distance, int = int, poly = poly, ...)
    co.names <- attr(C, "co.names")
    
    out[["Balance"]] <- balance.table.bin(C, weights = X$weights, treat = X$treat, 
                                          s.d.denom = X$s.d.denom, s.weights = X$s.weights, 
                                          continuous, binary, m.threshold = m.threshold, v.threshold = v.threshold, 
                                          ks.threshold = ks.threshold, un = un, disp.means = disp.means, 
                                          disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, 
                                          abs = abs, no.adj = no.adj, quick = quick, s.d.denom.list = X$s.d.denom.list)
    
    #Reassign disp... and ...threshold based on balance table output
    for (i in names(attr(out[["Balance"]], "disp"))) {
        assign(paste0("disp.", i), attr(out[["Balance"]], "disp")[i])
    }
    for (i in names(attr(out[["Balance"]], "threshold"))) {
        assign(paste0(i, ".threshold"), attr(out[["Balance"]], "threshold")[i])
    }
    
    S <- list(diff = list(threshold = m.threshold,
                          Names = "Means",
                          Threshold = "M.Threshold",
                          Stat = "Diff"),
              v.ratio = list(threshold = v.threshold,
                             Names = "Variances",
                             Threshold = "V.Threshold",
                             Stat = "V.Ratio"),
              ks = list(threshold = ks.threshold,
                        Names = "KS",
                        Threshold = "KS.Threshold",
                        Stat = "KS"))
    
    for (s in S) {
        if (is_not_null(s[["threshold"]])) {
            if (no.adj) {
                out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[paste.(s[["Threshold"]], "Un")]])
                out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Un"), paste.(s[["Threshold"]], "Un"), ratio = s$Stat == "V.Ratio")
            }
            else if (ncol(X$weights) == 1) {
                out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[s[["Threshold"]]]])
                out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio")
            }
            else if (ncol(X$weights) > 1) {
                out[[paste.("Balanced", s[["Names"]])]] <- setNames(do.call("cbind", lapply(names(X$weights), function(x) baltal(out[["Balance"]][[paste.(s[["Threshold"]], x)]]))),
                                                                    names(X$weights))
                out[[paste.("Max.Imbalance", s[["Names"]])]] <- cbind(Weights = names(X$weights),
                                                                      do.call("rbind", lapply(names(X$weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], x), paste.(s[["Threshold"]], x), ratio = s$Stat == "V.Ratio"),
                                                                                                                                     c("Variable", s[["Stat"]], s[["Threshold"]])))),
                                                                      stringsAsFactors = FALSE)
            }
        }
        else {
            out[[paste.("Balanced", s[["Names"]])]] <- NULL
            out[[paste.("Max.Imbalance", s[["Names"]])]] <- NULL
        }
    }
    
    out[["Observations"]] <- samplesize(treat = X$treat, weights = X$weights, s.weights = X$s.weights, method = X$method, discarded = X$discarded)
    out[["call"]] <- X$call
    attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                       v.threshold=v.threshold, 
                                       ks.threshold=ks.threshold, 
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.means=disp.means, 
                                       disp.sds=disp.sds,
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks, 
                                       disp.adj=!no.adj,
                                       disp.bal.tab = disp.bal.tab, 
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       nweights = ifelse(no.adj, 0, ncol(X$weights)),
                                       weight.names = names(X$weights),
                                       treat_names = treat_vals(X$treat),
                                       co.names = co.names)
    class(out) <- c("bal.tab.bin", "bal.tab")
    
    return(out)
}
base.bal.tab.cont <- function(X, int = FALSE, poly = 1, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "std"), r.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), abs = FALSE, quick = TRUE, ...) {
    
    #Preparations
    A <- clear_null(list(...))
    
    X$treat <- process_treat(X$treat) 
    
    if (is_not_null(r.threshold)) {
        r.threshold <- abs(r.threshold)
        if (r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            r.threshold <- NULL
        }
    }
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
    
    out.names <- c("Balance", "Balanced.Corr", 
                   "Max.Imbalance.Corr", 
                   "Observations", 
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    C <- get.C(covs = X$covs, int = int, poly = poly, addl = X$addl, distance = X$distance, ...)
    co.names <- attr(C, "co.names")
    
    out[["Balance"]] <- balance.table.cont(C, weights = X$weights, treat = X$treat, 
                                           s.d.denom = X$s.d.denom, s.weights = X$s.weights, 
                                           continuous = continuous, binary = binary, 
                                           r.threshold = r.threshold, un = un, 
                                           disp.means = disp.means, disp.sds = disp.sds, abs = abs, 
                                           no.adj = no.adj, quick = quick)
    
    #Reassign disp... and ...threshold based on balance table output
    for (i in names(attr(out[["Balance"]], "disp"))) {
        assign(paste0("disp.", i), attr(out[["Balance"]], "disp")[i])
    }
    for (i in names(attr(out[["Balance"]], "threshold"))) {
        assign(paste0(i, ".threshold"), attr(out[["Balance"]], "threshold")[i])
    }
    
    S <- list(corr = list(threshold = r.threshold,
                          Names = "Corr",
                          Threshold = "R.Threshold",
                          Stat = "Corr"))
    
    for (s in S) {
        if (is_not_null(s[["threshold"]])) {
            if (no.adj) {
                out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[paste.(s[["Threshold"]], "Un")]])
                out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Un"), paste.(s[["Threshold"]], "Un"), ratio = s$Stat == "V.Ratio")
            }
            else if (ncol(X$weights) == 1) {
                out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[s[["Threshold"]]]])
                out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio")
            }
            else if (ncol(X$weights) > 1) {
                out[[paste.("Balanced", s[["Names"]])]] <- setNames(do.call("cbind", lapply(names(X$weights), function(x) baltal(out[["Balance"]][[paste.(s[["Threshold"]], x)]]))),
                                                                    names(X$weights))
                out[[paste.("Max.Imbalance", s[["Names"]])]] <- cbind(Weights = names(X$weights),
                                                                      do.call("rbind", lapply(names(X$weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], x), paste.(s[["Threshold"]], x), ratio = s$Stat == "V.Ratio"),
                                                                                                                                     c("Variable", s[["Stat"]], s[["Threshold"]])))),
                                                                      stringsAsFactors = FALSE)
            }
        }
    }
    
    out[["Observations"]] <- samplesize.cont(treat = X$treat, weights = X$weights, s.weights = X$s.weights, method = X$method, discarded = X$discarded)
    out[["call"]] <- X$call
    attr(out, "print.options") <- list(r.threshold=r.threshold, 
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.means=disp.means, 
                                       disp.sds = disp.sds,
                                       disp.adj=!no.adj,
                                       disp.bal.tab = disp.bal.tab,
                                       continuous = continuous,
                                       binary = binary,
                                       abs = abs,
                                       quick = quick,
                                       nweights = if (no.adj) 0 else ncol(X$weights),
                                       weight.names = names(X$weights),
                                       co.names = co.names)
    class(out) <- c("bal.tab.cont", "bal.tab")
    
    return(out)
}
base.bal.tab.imp <- function(X, which.imp = NA, imp.summary = getOption("cobalt_imp.summary", TRUE), imp.fun = getOption("cobalt_imp.fun", NULL), ...) {
    A <- clear_null(list(...))
    
    X$treat <- process_treat(X$treat)
    
    #Preparations
    if (is_not_null(A$m.threshold)) A$m.threshold <- abs(A$m.threshold)
    if (is_not_null(A$v.threshold)) {
        A$v.threshold <- max(A$v.threshold, 1/A$v.threshold)
        A$disp.v.ratio <- TRUE
    }
    if (is_not_null(A$ks.threshold)) {
        if (A$ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            A$ks.threshold <- NULL
        }
        else A$disp.ks <- TRUE
    }
    if (is_not_null(A$r.threshold)) {
        A$r.threshold <- abs(A$r.threshold)
        if (A$r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            A$r.threshold <- NULL
        }
    }
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    imp <- factor(X$imp)
    
    if (is_not_null(A[["agg.fun"]])) A$imp.fun <- A[["agg.fun"]]
    
    #Setup output object
    out.names <- c("Imputation.Balance", 
                   "Balance.Across.Imputations", 
                   "Observations", 
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    #Get list of bal.tabs for each imputation
    
    out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) {
        X_i <- assign.X.class(subset_X(X, imp == i)) 
        X_i$call <- NULL
        do.call(base.bal.tab, c(list(X_i), A[names(A) %nin% names(X_i)]), quote = TRUE)
    })
    
    names(out[["Imputation.Balance"]]) <- levels(imp)
    
    #Create summary of lists
    
    if (imp.summary || !A$quick) {
        out[["Balance.Across.Imputations"]] <- balance.summary(out[["Imputation.Balance"]], 
                                                               weight.names = names(X$weights),
                                                               Agg.Fun = if_null_then(A$imp.fun, c("Min", "Mean", "Max")),
                                                               no.adj = is_null(X$weights),
                                                               abs = attr(out[["Imputation.Balance"]][[1]], "print.options")[["abs"]],
                                                               quick = A$quick,
                                                               types = NULL)
        
        observations <- lapply(out[["Imputation.Balance"]], function(x) x[["Observations"]])
        
        out[["Observations"]] <- samplesize.across.imps(observations)
    }
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- c(attr(out[["Imputation.Balance"]][[1]], "print.options"),
                                    list(which.imp = which.imp,
                                         imp.summary = imp.summary,
                                         imp.fun = imp.fun))
    class(out) <- c("bal.tab.imp", "bal.tab")
    
    return(out)
}
base.bal.tab.multi <- function(X, pairwise = TRUE, which.treat, multi.summary = getOption("cobalt_multi.summary", TRUE), ...) {
    A <- clear_null(list(...))
    
    X$treat <- process_treat(X$treat)
    
    #Preparations
    if (is_not_null(A$m.threshold)) A$m.threshold <- abs(A$m.threshold)
    if (is_not_null(A$v.threshold)) {
        A$v.threshold <- max(A$v.threshold, 1/A$v.threshold)
        A$disp.v.ratio <- TRUE
    }
    if (is_not_null(A$ks.threshold)) {
        if (A$ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            A$ks.threshold <- NULL
        }
        else A$disp.ks <- TRUE
    }
    if (is_not_null(A$r.threshold)) {
        A$r.threshold <- abs(A$r.threshold)
        if (A$r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            A$r.threshold <- NULL
        }
    }
    if (is_not_null(X$weights))  {
        check_if_zero_weights(X$weights, X$treat)
        if (ncol(X$weights) == 1) names(X$weights) <- "Adj"
    }
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    #Treat is a factor variable of 3+ levels
    if (is_null(X$focal)) {
        if (pairwise) treat.combinations <- combn(levels(X$treat), 2, list)
        else treat.combinations <- lapply(levels(X$treat), function(x) c(x, "Others"))
    }
    else {
        if (length(X$focal) > 1) stop("focal must be a vector of length 1 containing the name or index of the focal treatment group.", call. = FALSE)
        
        if (is.numeric(X$focal)) {
            X$focal <- levels(X$treat)[X$focal]
        }
        if (!is.character(X$focal)) {
            stop("focal must be the name or index of the focal treatment group.", call. = FALSE)
        }
        
        treat.combinations <- lapply(levels(X$treat)[levels(X$treat) != X$focal], function(x) rev(c(X$focal, x)))
        pairwise <- TRUE
    }
    
    #Setup output object
    out.names <- c("Pair.Balance", 
                   "Balance.Across.Pairs", 
                   "Observations", 
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    C <- do.call(get.C, c(X, A), quote = TRUE)
    bin.vars <- apply(C, 2, is_binary)
    if (is_null(X$weights)) {
        X$s.d.denom.list <- list(compute_s.d.denom(C, X$treat, s.d.denom = X$s.d.denom, s.weights = X$s.weights, bin.vars = bin.vars))
    }
    else {
        X$s.d.denom.list <- setNames(lapply(seq_along(X$s.d.denom), function(i) compute_s.d.denom(C, X$treat,
                                                                                                  s.d.denom = X$s.d.denom[i], s.weights = X$s.weights, 
                                                                                                  bin.vars = bin.vars, weighted.weights = X$weights[[i]])),
                                     names(X$s.d.denom))
    }
    
    if (pairwise || is_not_null(X$focal)) {
        balance.tables <- lapply(treat.combinations, function(t) {
            X_t <- assign.X.class(subset_X(X, X$treat %in% t))
            X_t$call <- NULL
            do.call(base.bal.tab, c(list(X_t), A[names(A) %nin% names(X_t)]), quote = TRUE)
        })
    }
    else {
        if (any(treat_vals(X$treat) == "Others")) stop ("\"Others\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
        balance.tables <- lapply(treat.combinations, function(t) {
            treat_ <- factor(X$treat, levels = c(levels(X$treat), "Others"))
            treat_[treat_ != t[1]] <- "Others"
            treat_ <- factor(treat_, rev(t))
            X_t <- X
            X_t$treat <- treat_
            X_t <- assign.X.class(X_t)
            X_t$call <- NULL
            do.call(base.bal.tab, c(list(X_t), A[names(A) %nin% names(X_t)]), quote = TRUE)
        })
    }
    
    for (i in seq_along(balance.tables)) {
        names(balance.tables)[i] <- paste(rev(treat.combinations[[i]]), collapse = " vs. ")
    }
    
    out[["Pair.Balance"]] <- balance.tables
    
    if (missing(which.treat)) {
        if (is_null(X$imp)) {
            which.treat <- NA
            multi.summary <- TRUE
        }
        else which.treat <- NULL
    }
    
    if ((multi.summary || !A$quick) && is_null(X$imp)) {
        out[["Balance.Across.Pairs"]] <- balance.summary(balance.tables, 
                                                         weight.names = names(X$weights),
                                                         Agg.Fun = "Max",
                                                         m.threshold = A$m.threshold,
                                                         v.threshold = A$v.threshold,
                                                         ks.threshold = A$ks.threshold,
                                                         no.adj = is_null(X$weights),
                                                         quick = A$quick,
                                                         types = NULL)
        out[["Observations"]] <- samplesize.multi(balance.tables, treat_names(X$treat), X$focal)
    }
    
    out[["call"]] <- X$call
    
    attr(out, "print.options") <- c(attr(out[["Pair.Balance"]][[1]], "print.options"),
                                    list(treat_names_multi = treat_vals(X$treat),
                                         which.treat = which.treat,
                                         multi.summary = multi.summary,
                                         pairwise = pairwise))
    
    attr(out, "print.options")[["treat_names"]] <- NULL
    
    class(out) <- c("bal.tab.multi", "bal.tab")
    
    return(out)
    
}
base.bal.tab.msm <- function(X, which.time = NULL, msm.summary = getOption("cobalt_msm.summary", TRUE), ...) {
    #One vector of weights
    #treat.list should be a df/list of treatment vectors, one for each time period
    #cov.list should be a list of covariate data.frames, one for each time period; 
    #   should include all covs from previous time points, but no treatment statuses
    
    A <- clear_null(list(...))
    
    X$treat.list <- process_treat.list(X$treat)
    
    #Preparations
    if (is_not_null(A$m.threshold)) A$m.threshold <- abs(A$m.threshold)
    if (is_not_null(A$v.threshold)) {
        A$v.threshold <- max(A$v.threshold, 1/A$v.threshold)
        A$disp.v.ratio <- TRUE
    }
    if (is_not_null(A$ks.threshold)) {
        if (A$ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            A$ks.threshold <- NULL
        }
        else A$disp.ks <- TRUE
    }
    if (is_not_null(A$r.threshold)) {
        A$r.threshold <- abs(A$r.threshold)
        if (A$r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            A$r.threshold <- NULL
        }
    }
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    #Setup output object
    out.names <- c("Time.Balance", 
                   "Balance.Across.Times", 
                   "Observations", 
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    out[["Time.Balance"]] <- vector("list", length(X$covs.list))
    
    treat.types <- vapply(X$treat.list, function(x) get.treat.type(x), character(1L))
    
    #Get list of bal.tabs for each time period
    out[["Time.Balance"]] <- lapply(seq_along(X$covs.list), function(ti) {
        X_ti <- X
        X_ti <- c(X_ti, list(
            covs = X_ti$covs.list[[ti]], 
            treat = X_ti$treat.list[[ti]], 
            addl = X_ti$addl.list[[ti]], 
            distance = X_ti$distance.list[[ti]]
        ))
        X_ti[c("covs.list", "treat.list", "addl.list", "distance.list")] <- NULL
        X_ti$call <- NULL
        X_ti <- assign.X.class(X_ti)
        
        do.call(base.bal.tab, c(list(X_ti), A[names(A) %nin% names(X_ti)]), quote = TRUE)
    })
    
    if (length(names(X$treat.list)) == length(X$treat.list)) {
        names(out[["Time.Balance"]]) <- names(X$treat.list)
    }
    else names(out[["Time.Balance"]]) <- seq_along(X$treat.list)
    
    if (!(A$quick && !msm.summary) && all_the_same(treat.types) && "multinomial" %nin% treat.types && is_null(X$imp)) {
        out[["Balance.Across.Times"]] <- balance.table.msm.summary(out[["Time.Balance"]],
                                                                   weight.names = names(X$weights),
                                                                   no.adj = is_null(X$weights),
                                                                   m.threshold = A$m.threshold, 
                                                                   v.threshold = A$v.threshold, 
                                                                   ks.threshold = A$ks.threshold, 
                                                                   r.threshold = A$r.threshold, 
                                                                   quick = A$quick, 
                                                                   types = NULL)
        out[["Observations"]] <- lapply(out[["Time.Balance"]], function(x) x$Observations)
    }
    
    out[["call"]] <- X$call
    
    attr(out, "print.options") <- c(attr(out[["Time.Balance"]][[1]], "print.options"),
                                    list(which.time = which.time,
                                         msm.summary = msm.summary))
    
    class(out) <- c("bal.tab.msm", "bal.tab")
    
    return(out)
}
base.bal.tab.cluster <- function(X, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), ...) {
    A <- clear_null(list(...))
    
    #Preparations
    if (is_not_null(A$m.threshold)) A$m.threshold <- abs(A$m.threshold)
    if (is_not_null(A$v.threshold)) {
        A$v.threshold <- max(A$v.threshold, 1/A$v.threshold)
        A$disp.v.ratio <- TRUE
    }
    if (is_not_null(A$ks.threshold)) {
        if (A$ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            A$ks.threshold <- NULL
        }
        else A$disp.ks <- TRUE
    }
    if (is_not_null(A$r.threshold)) {
        A$r.threshold <- abs(A$r.threshold)
        if (A$r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            A$r.threshold <- NULL
        }
    }
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    cluster <- factor(X$cluster)
    
    if (is_not_null(A[["agg.fun"]])) A$cluster.fun <- A[["agg.fun"]]
    
    #Setup output object
    out.names <- c("Cluster.Balance", 
                   "Cluster.Summary", 
                   "Observations", 
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    #Get list of bal.tabs for each imputation
    
    out[["Cluster.Balance"]] <- lapply(levels(cluster), function(cl) {
        X_cl <- assign.X.class(subset_X(X, cluster == cl)) 
        X_cl$call <- NULL
        do.call(base.bal.tab, c(list(X_cl), A[names(A) %nin% names(X_cl)]), quote = TRUE)
    })
    
    names(out[["Cluster.Balance"]]) <- levels(cluster)
    
    #Create summary of lists
    
    if ((cluster.summary || !A$quick) && is_null(X$covs.list) && get.treat.type(X$treat) != "multinomial" && is_null(X$imp)) {
        out[["Cluster.Summary"]] <- balance.summary(out[["Cluster.Balance"]], 
                                                    weight.names = names(X$weights),
                                                    Agg.Fun = if_null_then(A$cluster.fun, c("Min", "Mean", "Max")),
                                                    no.adj = is_null(X$weights),
                                                    abs = A$abs,
                                                    quick = A$quick,
                                                    types = NULL)
        observations <- lapply(out[["Cluster.Balance"]], function(x) x[["Observations"]])
        
        out[["Observations"]] <- samplesize.across.clusters(observations)
    }
    
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- c(attr(out[["Cluster.Balance"]][[1]], "print.options"),
                                    list(which.cluster = which.cluster,
                                         cluster.summary = cluster.summary,
                                         cluster.fun = cluster.fun))
    class(out) <- c("bal.tab.cluster", "bal.tab")
    
    return(out)
}

base.bal.tab.subclass.bin <- function(X, int = FALSE, poly = 1, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), abs = FALSE, quick = TRUE, ...) {
    #Preparations
    A <- clear_null(list(...))
    
    if (!is_(X$treat, "processed.treat")) X$treat <- process_treat(X$treat) 
    
    if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
    if (is_not_null(v.threshold)) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (is_null(ks.threshold) && is_null(A$k.threshold)) {
        ks.threshold <- A$k.threshold
    }
    if (is_not_null(ks.threshold)) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    
    no.adj <- FALSE
    
    if (is_null(X$s.weights)) {
        X$s.weights <- rep(1, length(X$treat))
    }
    
    #Actions
    
    out.names <- c("Subclass.Balance", "Balance.Across.Subclass", 
                   expand.grid_string(c("Balanced", "Max.Imbalance"),
                                      c("Means", "Variances", "KS"),
                                      "Subclass", collapse = "."), 
                   "Subclass.Observations", "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    C <- get.C(covs = X$covs, int = int, poly = poly, addl = X$addl, distance = X$distance, ...)
    co.names <- attr(C, "co.names")
    
    out[["Subclass.Balance"]] <- balance.table.subclass.bin(C, weights=NULL, treat=X$treat, subclass=X$subclass, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom[1], m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold = ks.threshold, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, abs = abs, quick = quick)
    out[["Subclass.Observations"]] <- samplesize(treat = X$treat, weights = NULL, subclass = X$subclass, s.weights = X$s.weights, method = X$method, discarded = X$discarded)
    out[["Balance.Across.Subclass"]] <- balance.table.bin(C, 
                                                          weights = data.frame(Adj = strata2weights(X$subclass, X$treat)),
                                                          treat = X$treat,
                                                          continuous = continuous,
                                                          binary = binary,
                                                          s.d.denom = X$s.d.denom[1],
                                                          m.threshold = m.threshold,
                                                          v.threshold = v.threshold,
                                                          ks.threshold = ks.threshold,
                                                          un = un,
                                                          disp.means = disp.means,
                                                          disp.sds = disp.sds,
                                                          disp.v.ratio = disp.v.ratio,
                                                          disp.ks = disp.ks, no.adj = FALSE,
                                                          abs = abs, quick = quick)
    
    #Reassign disp... and ...threshold based on balance table output
    for (i in names(attr(out[["Subclass.Balance"]], "disp"))) {
        assign(paste.("disp", i), attr(out[["Subclass.Balance"]], "disp")[i])
    }
    for (i in names(attr(out[["Subclass.Balance"]], "threshold"))) {
        assign(paste.(i, "threshold"), attr(out[["Subclass.Balance"]], "threshold")[i])
    }
    
    S <- list(diff = list(threshold = m.threshold,
                          Names = "Means",
                          Threshold = "M.Threshold",
                          Stat = "Diff"),
              v.ratio = list(threshold = v.threshold,
                             Names = "Variances",
                             Threshold = "V.Threshold",
                             Stat = "V.Ratio"),
              ks = list(threshold = ks.threshold,
                        Names = "KS",
                        Threshold = "KS.Threshold",
                        Stat = "KS"))
    
    for (s in S) {
        if (is_not_null(s[["threshold"]])) {
            out[[paste.("Balanced", s[["Names"]], "Subclass")]] <- as.data.frame(lapply(levels(X$subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][[s[["Threshold"]]]])))
            names(out[[paste.("Balanced", s[["Names"]], "Subclass")]]) <- paste("Subclass", levels(X$subclass))
            max.imbal.list <- lapply(levels(X$subclass), function(x) {
                return(max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][["Type"]] != "Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio"))
            } )
            out[[paste.("Max.Imbalance", s[["Names"]], "Subclass")]] <- data.frame(do.call("rbind", max.imbal.list), 
                                                                                   row.names = paste("Subclass", levels(X$subclass)))
        }
        else {
            out[[paste.("Balanced", s[["Names"]], "Subclass")]] <- NULL
            out[[paste.("Max.Imbalance", s[["Names"]], "Subclass")]] <- NULL
        }
    }
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                       v.threshold=v.threshold, 
                                       ks.threshold=ks.threshold, 
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.means=disp.means, 
                                       disp.sds=disp.sds,
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks, 
                                       disp.adj=!no.adj, 
                                       disp.subclass=disp.subclass,
                                       disp.bal.tab = disp.bal.tab, 
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       treat_names = treat_vals(X$treat),
                                       co.names = co.names)
    class(out) <- c("bal.tab.subclass", "bal.tab.bin", "bal.tab")
    
    
    return(out)
}
base.bal.tab.subclass.cont <- function(X, int = FALSE, poly = 1, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "std"), r.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), abs = FALSE, quick = TRUE, ...) {
    #Preparations
    A <- clear_null(list(...))
    
    if (!is_(X$treat, "processed.treat")) X$treat <- process_treat(X$treat) 
    
    if (is_not_null(r.threshold)) r.threshold <- abs(r.threshold)
    
    no.adj <- FALSE
    
    if (is_null(X$s.weights)) {
        X$s.weights <- rep(1, length(X$treat))
    }
    
    #Actions
    
    out.names <- c("Subclass.Balance", "Balance.Across.Subclass", 
                   paste.(c("Balanced", "Max.Imbalance"),
                          "Corr",
                          "Subclass"), 
                   "Subclass.Observations", "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    C <- get.C(covs = X$covs, int = int, poly = poly, addl = X$addl, distance = X$distance, ...)
    co.names <- attr(C, "co.names")
    
    out[["Subclass.Balance"]] <- balance.table.subclass.cont(C, weights=NULL, treat=X$treat, subclass=X$subclass, continuous=continuous, binary=binary, r.threshold=r.threshold, disp.means = disp.means, disp.sds = disp.sds, abs = abs, quick = quick)
    out[["Subclass.Observations"]] <- samplesize.cont(treat = X$treat, weights = NULL, subclass = X$subclass, s.weights = X$s.weights, method = X$method, discarded = X$discarded)
    out[["Balance.Across.Subclass"]] <- balance.table.across.subclass.cont(balance.table.cont(C, weights = NULL, treat = X$treat, 
                                                                                              s.d.denom = X$s.d.denom, s.weights = X$s.weights, 
                                                                                              continuous = continuous, binary = binary, 
                                                                                              r.threshold = r.threshold, un = un, 
                                                                                              disp.means = disp.means, disp.sds = disp.sds, abs = abs, 
                                                                                              no.adj = TRUE, quick = quick), 
                                                                           balance.table.subclass.list = out[["Subclass.Balance"]], 
                                                                           subclass.obs = out[["Subclass.Observations"]], 
                                                                           r.threshold = r.threshold)
    
    #Reassign disp... and ...threshold based on balance table output
    for (i in names(attr(out[["Subclass.Balance"]], "disp"))) {
        assign(paste.("disp", i), attr(out[["Subclass.Balance"]], "disp")[i])
    }
    for (i in names(attr(out[["Subclass.Balance"]], "threshold"))) {
        assign(paste.(i, "threshold"), attr(out[["Subclass.Balance"]], "threshold")[i])
    }
    
    S <- list(corr = list(threshold = r.threshold,
                          Names = "Corr",
                          Threshold = "R.Threshold",
                          Stat = "Corr"))
    
    for (s in S) {
        if (is_not_null(s[["threshold"]])) {
            out[[paste.("Balanced", s[["Names"]], "Subclass")]] <- as.data.frame(lapply(levels(X$subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][[s[["Threshold"]]]])))
            names(out[[paste.("Balanced", s[["Names"]], "Subclass")]]) <- paste("Subclass", levels(X$subclass))
            max.imbal.list <- lapply(levels(X$subclass), function(x) {
                return(max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][["Type"]] != "Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio"))
            } )
            out[[paste.("Max.Imbalance", s[["Names"]], "Subclass")]] <- data.frame(do.call("rbind", max.imbal.list), 
                                                                                   row.names = paste("Subclass", levels(X$subclass)))
        }
        else {
            out[[paste.("Balanced", s[["Names"]], "Subclass")]] <- NULL
            out[[paste.("Max.Imbalance", s[["Names"]], "Subclass")]] <- NULL
        }
    }
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- list(r.threshold=r.threshold, 
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.means=disp.means, 
                                       disp.sds=disp.sds,
                                       disp.adj=!no.adj, 
                                       disp.subclass=disp.subclass,
                                       disp.bal.tab = disp.bal.tab, 
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       co.names = co.names)
    class(out) <- c("bal.tab.subclass", "bal.tab.cont", "bal.tab")
    
    
    return(out)
}
