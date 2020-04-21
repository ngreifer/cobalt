base.bal.tab <- function(X, ...) {
    UseMethod("base.bal.tab")
}
base.bal.tab.base <- function(X, type, int = FALSE, poly = 1, continuous, binary, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp = NULL, disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), abs = FALSE, quick = TRUE, ...) {
    #Preparations

    X$treat <- process_treat(X$treat) 
    
    if (type == "bin") {
        if (get.treat.type(X$treat) != "binary") {
            stop("Treatment indicator must be a binary variable---e.g., treatment (1) or control (0)", call. = FALSE)}
        if (missing(continuous)) continuous <- getOption("cobalt_continuous", "std")
        if (missing(binary)) binary <- getOption("cobalt_binary", "raw")
    }
    else if (type == "cont"){
        if (missing(continuous)) continuous <- getOption("cobalt_continuous", "std")
        if (missing(binary)) binary <- getOption("cobalt_binary", "std")
    }
    
    if (is_null(X$weights)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (type == "bin") check_if_zero_weights(X$weights, X$treat)
        else if (type == "cont") check_if_zero_weights(X$weights)
        if (ncol(X$weights) == 1) names(X$weights) <- "Adj"
    }
    if (is_null(X$s.weights)) {
        X$s.weights <- rep(1, length(X$treat))
    }
    
    disp <- process_disp(disp, ...)
    
    #Actions
    out.names <- c("Balance", 
                   expand.grid_string(c("Balanced", "Max.Imbalance"),
                                      X$stats,
                                      collapse = "."), 
                   "Observations", "call")
    out <- make_list(out.names)
    
    C <- get.C(covs = X$covs, addl = X$addl, distance = X$distance, int = int, poly = poly, ...)
    co.names <- attr(C, "co.names")
    
    out[["Balance"]] <- balance.table(C, type = type, weights = X$weights, treat = X$treat, 
                                      s.d.denom = X$s.d.denom, s.weights = X$s.weights, 
                                      continuous = continuous, binary = binary, 
                                      thresholds = X$thresholds,
                                      un = un, disp = disp, 
                                      stats = X$stats, abs = abs, 
                                      no.adj = no.adj, quick = quick, 
                                      s.d.denom.list = X$s.d.denom.list, ...)
    
    #Reassign disp... and ...threshold based on balance table output
    compute <- attr(out[["Balance"]], "compute")
    thresholds <- attr(out[["Balance"]], "thresholds")
    disp <- attr(out[["Balance"]], "disp")
    
    for (s in compute) {
        if (is_not_null(thresholds[[s]])) {
            if (no.adj) {
                out[[paste.("Balanced", s)]] <- baltal(out[["Balance"]][[paste.(STATS[[s]]$Threshold, "Un")]])
                out[[paste.("Max.Imbalance", s)]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], 
                                                               col.name = paste.(STATS[[s]]$bal.tab_column_prefix, "Un"), 
                                                               thresh.col.name = paste.(STATS[[s]]$Threshold, "Un"), 
                                                               abs_stat = STATS[[s]]$abs)
            }
            else if (ncol(X$weights) == 1) {
                out[[paste.("Balanced", s)]] <- baltal(out[["Balance"]][[STATS[[s]]$Threshold]])
                out[[paste.("Max.Imbalance", s)]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], 
                                                               col.name = paste.(STATS[[s]]$bal.tab_column_prefix, "Adj"), 
                                                               thresh.col.name = STATS[[s]]$Threshold, 
                                                               abs_stat = STATS[[s]]$abs)
            }
            else if (ncol(X$weights) > 1) {
                out[[paste.("Balanced", s)]] <- setNames(do.call("cbind", lapply(names(X$weights), function(x) baltal(out[["Balance"]][[paste.(STATS[[s]]$Threshold, x)]]))),
                                                         names(X$weights))
                out[[paste.("Max.Imbalance", s)]] <- cbind(Weights = names(X$weights),
                                                           do.call("rbind", lapply(names(X$weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], 
                                                                                                                                    col.name = paste.(STATS[[s]]$bal.tab_column_prefix, x), 
                                                                                                                                    thresh.col.name = paste.(STATS[[s]]$Threshold, x), 
                                                                                                                                    abs_stat = STATS[[s]]$abs),
                                                                                                                          c("Variable", 
                                                                                                                            STATS[[s]]$bal.tab_column_prefix, 
                                                                                                                            STATS[[s]]$Threshold)))),
                                                           stringsAsFactors = FALSE)
            }
        }
        else {
            out[[paste.("Balanced", s)]] <- NULL
            out[[paste.("Max.Imbalance", s)]] <- NULL
        }
    }
    
    out[["Observations"]] <- samplesize(treat = X$treat, type = type, weights = X$weights, s.weights = X$s.weights, method = X$method, discarded = X$discarded)

    out[["call"]] <- X$call
    attr(out, "print.options") <- list(thresholds = thresholds,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       compute = compute, 
                                       disp = disp,
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
    class(out) <- c(paste.("bal.tab", type), "bal.tab")
    
    return(out)
}
base.bal.tab.binary <- function(X, ...) {
    base.bal.tab.base(X, type = "bin", ...)
}
base.bal.tab.cont <- function(X, ...) {
    base.bal.tab.base(X, type = "cont", ...)
}

base.bal.tab.imp <- function(X, which.imp = NA, imp.summary = getOption("cobalt_imp.summary", TRUE), imp.fun = getOption("cobalt_imp.fun", NULL), ...) {
    A <- clear_null(list(...))
    
    X$treat <- process_treat(X$treat)
    
    #Preparations
    
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    imp <- factor(X$imp)
    
    if (is_not_null(A[["agg.fun"]])) imp.fun <- A[["agg.fun"]]
    
    #Setup output object
    out.names <- c("Imputation.Balance", 
                   "Balance.Across.Imputations", 
                   "Observations", 
                   "call")
    out <- make_list(out.names)
    
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
                                                               agg.funs = if_null_then(imp.fun, c("min", "mean", "max")))
        
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
    
    if (is_not_null(X$weights))  {
        check_if_zero_weights(X$weights, X$treat)
        if (ncol(X$weights) == 1) names(X$weights) <- "Adj"
    }
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    #Treat is a factor variable of 3+ levels
    if (is_null(X$focal)) {
        if (pairwise) treat.combinations <- combn(levels(X$treat), 2, list)
        else treat.combinations <- lapply(levels(X$treat), function(x) c(x, "All"))
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
    out <- make_list(out.names)
    
    C <- do.call(get.C, c(X, A), quote = TRUE)
    bin.vars <- is_binary_col(C)
    if ("mean.diffs" %in% X$stats) {
        if (is_null(X$weights)) {
            X$s.d.denom.list <- list(compute_s.d.denom(C, X$treat, s.d.denom = X$s.d.denom, s.weights = X$s.weights, bin.vars = bin.vars))
        }
        else {
            X$s.d.denom.list <- setNames(lapply(seq_along(X$s.d.denom), function(i) compute_s.d.denom(C, X$treat,
                                                                                                      s.d.denom = X$s.d.denom[i], s.weights = X$s.weights, 
                                                                                                      bin.vars = bin.vars, weighted.weights = X$weights[[i]])),
                                         names(X$s.d.denom))
        }
    }
    
    if (pairwise || is_not_null(X$focal)) {
        balance.tables <- lapply(treat.combinations, function(t) {
            X_t <- assign.X.class(subset_X(X, X$treat %in% t))
            X_t$call <- NULL
            do.call(base.bal.tab, c(list(X_t), A[names(A) %nin% names(X_t)]), quote = TRUE)
        })
    }
    else {
        if (any(treat_vals(X$treat) == "All")) stop ("\"All\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
        balance.tables <- lapply(treat.combinations, function(t) {
            n <- length(X$treat)
            X_t <- X
            X_t$call <- NULL
            X_t <- subset_X(X_t, c(seq_len(n), which(X$treat == t[1])))
            X_t$treat <- factor(c(rep("All", n), rep(t[1], sum(X$treat == t[1]))), nmax = 2,
                                levels = c("All", t[1]))
            X_t <- assign.X.class(X_t)
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
                                                         agg.funs = "max")
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
    
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    #Setup output object
    out.names <- c("Time.Balance", 
                   "Balance.Across.Times", 
                   "Observations", 
                   "call")
    out <- make_list(out.names)
    
    out[["Time.Balance"]] <- make_list(length(X$covs.list))
    
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
        out[["Balance.Across.Times"]] <- balance.summary(out[["Time.Balance"]],
                                                         agg.funs = "max",
                                                         include.times = TRUE)
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
    
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    if (is_null(A[["abs"]])) A[["abs"]] <- FALSE
    
    cluster <- factor(X$cluster)
    
    if (is_not_null(A[["agg.fun"]])) cluster.fun <- A[["agg.fun"]]
    
    #Setup output object
    out.names <- c("Cluster.Balance", 
                   "Balance.Across.Clusters", 
                   "Observations", 
                   "call")
    out <- make_list(out.names)
    
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
                                                    agg.funs = if_null_then(cluster.fun, c("min", "mean", "max")))
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

#NEEDS UPDATING with STATS
base.bal.tab.subclass <- function(X, type, int = FALSE, poly = 1, continuous, binary, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp = NULL, disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), abs = FALSE, quick = TRUE, ...) {
    #Preparations
    if (type == "bin") {
        if(get.treat.type(X$treat) != "binary") 
            stop("Treatment indicator must be a binary variable---e.g., treatment (1) or control (0)", call. = FALSE)
        if (missing(continuous)) continuous <- getOption("cobalt_continuous", "std")
        if (missing(binary)) binary <- getOption("cobalt_binary", "raw")
    }
    else {
        if (missing(continuous)) continuous <- getOption("cobalt_continuous", "std")
        if (missing(binary)) binary <- getOption("cobalt_binary", "std")
    }
    
    no.adj <- FALSE

    if (is_null(X$s.weights)) {
        X$s.weights <- rep(1, length(X$treat))
    }
    
    disp <- process_disp(disp, ...)
    
    #Actions
    
    out.names <- c("Subclass.Balance", "Balance.Across.Subclass", 
                   expand.grid_string(c("Balanced", "Max.Imbalance"),
                                      X$stats,
                                      "Subclass",
                                      collapse = "."), 
                   "Observations", "call")
    out <- make_list(out.names)
    
    C <- get.C(covs = X$covs, addl = X$addl, distance = X$distance, int = int, poly = poly, ...)
    co.names <- attr(C, "co.names")
    
    out[["Observations"]] <- samplesize(treat = X$treat, 
                                        type = type,
                                        weights = NULL, 
                                        subclass = X$subclass,
                                        s.weights = X$s.weights, 
                                        method = X$method, 
                                        discarded = X$discarded)
    out[["Subclass.Balance"]] <- balance.table.subclass(C, type = type, 
                                                        weights = NULL, 
                                                        treat = X$treat, 
                                                        subclass = X$subclass,
                                                        continuous=continuous, binary=binary, 
                                                        s.d.denom=X$s.d.denom[1], 
                                                        thresholds = X$thresholds,
                                                        disp = disp,
                                                        stats = X$stats, 
                                                        s.weights = X$s.weights, 
                                                        abs = abs, 
                                                        quick = quick)
    if (type == "bin") {
        out[["Balance.Across.Subclass"]] <- balance.table(C, 
                                                          type = type, 
                                                          weights = data.frame(Adj = strata2weights(X$subclass, X$treat)), 
                                                          treat = X$treat, 
                                                          s.d.denom = X$s.d.denom[1], 
                                                          s.weights = X$s.weights, 
                                                          continuous = continuous, 
                                                          binary = binary, 
                                                          thresholds = X$thresholds,
                                                          un = un, 
                                                          disp = disp,
                                                          stats = X$stats, 
                                                          abs = abs, 
                                                          no.adj = FALSE, quick = quick, 
                                                          s.d.denom.list = X$s.d.denom.list)
    }
    else if (type == "cont") {
        out[["Balance.Across.Subclass"]] <- balance.table.across.subclass.cont(balance.table(C, 
                                                                                             type = type, 
                                                                                             weights = NULL,
                                                                                             treat = X$treat, 
                                                                                             s.d.denom = X$s.d.denom[1], 
                                                                                             s.weights = X$s.weights, 
                                                                                             continuous = continuous, 
                                                                                             binary = binary, 
                                                                                             thresholds = X$thresholds,
                                                                                             un = un, 
                                                                                             disp = disp,
                                                                                             stats = X$stats, 
                                                                                             abs = abs, 
                                                                                             no.adj = TRUE, 
                                                                                             quick = quick, 
                                                                                             s.d.denom.list = X$s.d.denom.list), 
                                                                               balance.table.subclass.list = out[["Subclass.Balance"]], 
                                                                               subclass.obs = out[["Observations"]], 
                                                                               r.threshold = X$thresholds[["correlations"]])
    }
    
    #Reassign disp... and ...threshold based on balance table output
    compute <- attr(out[["Subclass.Balance"]], "compute")
    thresholds <- attr(out[["Subclass.Balance"]], "thresholds")
    disp <- attr(out[["Subclass.Balance"]], "disp")
    
    for (s in compute) {
        if (is_not_null(thresholds[[s]])) {
            out[[paste.("Balanced", s, "Subclass")]] <- setNames(do.call(data.frame, lapply(out[["Subclass.Balance"]], function(x) baltal(x[[STATS[[s]]$Threshold]]))),
                                                                 paste("Subclass", levels(X$subclass)))
            max.imbal.list <- lapply(out[["Subclass.Balance"]], function(x) {
                return(max.imbal(x[x[["Type"]] != "Distance", , drop = FALSE], 
                                 col.name = paste.(STATS[[s]]$bal.tab_column_prefix, "Adj"), 
                                 thresh.col.name = STATS[[s]]$Threshold, 
                                 abs_stat = STATS[[s]]$abs))
            } )
            out[[paste.("Max.Imbalance", s, "Subclass")]] <- data.frame(do.call(rbind, max.imbal.list), 
                                                                        row.names = paste("Subclass", levels(X$subclass)))
        }
        else {
            out[[paste.("Balanced", s, "Subclass")]] <- NULL
            out[[paste.("Max.Imbalance", s, "Subclass")]] <- NULL
        }
    }
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- list(thresholds = thresholds,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp=disp,
                                       compute = compute, 
                                       disp.adj = !no.adj, 
                                       disp.subclass = disp.subclass,
                                       disp.bal.tab = disp.bal.tab, 
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       treat_names = treat_vals(X$treat),
                                       type = type,
                                       co.names = co.names)
    class(out) <- c("bal.tab.subclass", paste.("bal.tab", type), "bal.tab")
    
    return(out)
}
base.bal.tab.subclass.binary <- function(X, ...) {
    base.bal.tab.subclass(X, type = "bin", ...)
}
base.bal.tab.subclass.cont <- function(X, ...) {
    base.bal.tab.subclass(X, type = "cont", ...)
}