base.bal.tab <- function(X, ...) {
    UseMethod("base.bal.tab")
}
base.bal.tab.base <- function(X, type, int = FALSE, poly = 1, continuous, binary, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp = NULL, disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), disp.call = getOption("cobalt_disp.call", TRUE), abs = FALSE, quick = TRUE, ...) {
    #Preparations

    A <- clear_null(list(...))
    A$subset <- NULL
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
    
    disp <- do.call("process_disp", c(list(disp), A), quote = TRUE)
    
    #Actions
    out <- list()
    
    C <- do.call("get.C2", c(X, A[names(A) %nin% names(X)], list(int = int, poly = poly)), quote = TRUE)
    
    co.names <- attr(C, "co.names")
    
    out[["Balance"]] <- do.call("balance.table", c(list(C, type = type, weights = X$weights, treat = X$treat, 
                                                        s.d.denom = X$s.d.denom, s.weights = X$s.weights, 
                                                        continuous = continuous, binary = binary, 
                                                        thresholds = X$thresholds,
                                                        un = un, disp = disp, 
                                                        stats = X$stats, abs = abs, 
                                                        no.adj = no.adj, quick = quick, 
                                                        s.d.denom.list = X$s.d.denom.list), A), quote = TRUE)
    
    #Reassign disp... and ...threshold based on balance table output
    compute <- attr(out[["Balance"]], "compute")
    thresholds <- attr(out[["Balance"]], "thresholds")
    disp <- attr(out[["Balance"]], "disp")
    
    out <- c(out, threshold.summary(compute = compute,
                                    thresholds = thresholds,
                                    no.adj = no.adj,
                                    balance.table = out[["Balance"]],
                                    weight.names = names(X$weights)))
    
    out[["Observations"]] <- samplesize(treat = X$treat, type = type, weights = X$weights, s.weights = X$s.weights, method = X$method, discarded = X$discarded)
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- list(thresholds = thresholds,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       compute = compute, 
                                       disp = disp,
                                       disp.adj=!no.adj,
                                       disp.bal.tab = disp.bal.tab,
                                       disp.call = disp.call, 
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       nweights = ifelse(no.adj, 0, ncol(X$weights)),
                                       weight.names = names(X$weights),
                                       treat_names = treat_names(X$treat),
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

base.bal.tab.imp <- function(X, which.imp = NA, imp.summary = getOption("cobalt_imp.summary"), imp.fun = getOption("cobalt_imp.fun", NULL), ...) {
    A <- list(...)
    
    # X$treat <- process_treat(X$treat)
    
    #Preparations
    
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    imp <- factor(X$imp)
    
    if (missing(which.imp)) {
        which.imp <- NA
    }
    
    if (is_null(imp.summary)) {
        imp.summary <- is_not_null(which.imp) && 
            (anyNA(which.imp) || !is.numeric(which.imp) || 
                 (is.numeric(which.imp) && !any(which.imp %in% seq_len(nlevels(imp)))))
    }
    
    all.agg.funs <- c("min", "mean", "max")
    agg.fun <- tolower(as.character(if_null_then(imp.fun, A[["agg.fun"]], all.agg.funs)))
    agg.fun <- match_arg(agg.fun, all.agg.funs, several.ok = TRUE)
    
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
        do.call("base.bal.tab", c(list(X_i), A[names(A) %nin% names(X_i)]), quote = TRUE)
    })
    
    names(out[["Imputation.Balance"]]) <- levels(imp)
    
    #Create summary of lists
    
    if (imp.summary || !A$quick) {
        out[["Balance.Across.Imputations"]] <- balance.summary(out[["Imputation.Balance"]], 
                                                               agg.funs = agg.fun)
        
        if (length(agg.fun) == 1) {
            out <- c(out, threshold.summary(compute = attr(out[["Imputation.Balance"]][[1]][["Balance"]], "compute"),
                                            thresholds = attr(out[["Imputation.Balance"]][[1]][["Balance"]], "thresholds"),
                                            no.adj = !attr(out[["Imputation.Balance"]][[1]], "print.options")$disp.adj,
                                            balance.table = out[["Balance.Across.Imputations"]],
                                            weight.names = attr(out[["Imputation.Balance"]][[1]], "print.options")$weight.names,
                                            agg.fun = agg.fun))
        }
        
        observations <- lapply(out[["Imputation.Balance"]], function(x) x[["Observations"]])
        
        out[["Observations"]] <- samplesize.across.imps(observations)
    }
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- c(attr(out[["Imputation.Balance"]][[1]], "print.options"),
                                    list(which.imp = which.imp,
                                         imp.summary = imp.summary,
                                         imp.fun = agg.fun))
    class(out) <- c("bal.tab.imp", "bal.tab")
    
    return(out)
}
base.bal.tab.multi <- function(X, pairwise = TRUE, which.treat, multi.summary = getOption("cobalt_multi.summary"), ...) {
    A <- list(...)
    
    # X$treat <- process_treat(X$treat)
    
    #Preparations
    
    if (is_not_null(X$weights))  {
        check_if_zero_weights(X$weights, X$treat)
        if (ncol(X$weights) == 1) names(X$weights) <- "Adj"
    }
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    if (missing(which.treat)) {
        if (is_null(X$imp)) which.treat <- NA
        else which.treat <- NULL
    }
    
    if (is_null(multi.summary)) {
        multi.summary <- is_not_null(which.treat) && anyNA(which.treat)
    }
    
    #Treat is a factor variable
    if (is_null(X$focal)) {
        if (pairwise) treat.combinations <- combn(treat_names(X$treat), 2, simplify = FALSE)
        else treat.combinations <- lapply(treat_names(X$treat), function(x) c(x, "All"))
    }
    else {
        if (length(X$focal) > 1) stop("'focal' must be a vector of length 1 containing the name or index of the focal treatment group.", call. = FALSE)
        
        if (is.numeric(X$focal)) {
            X$focal <- levels(X$treat)[X$focal]
        }
        if (!is.character(X$focal)) {
            stop("'focal' must be the name or index of the focal treatment group.", call. = FALSE)
        }
        
        treat.combinations <- lapply(levels(X$treat)[levels(X$treat) != X$focal], function(x) rev(c(X$focal, x)))
        pairwise <- TRUE
    }
    
    #Setup output object
    # out.names <- c("Pair.Balance", 
    #                "Balance.Across.Pairs", 
    #                "Observations", 
    #                "call")
    # out <- make_list(out.names)
    out <- list()
    
    if ("mean.diffs" %in% X$stats) {
        C <- do.call("get.C2", c(X, A[names(A) %nin% names(X)]), quote = TRUE)
        bin.vars <- is_binary_col(C)
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
            X_t <- subset_X(X, X$treat %in% treat_vals(X$treat)[t])
            X_t$treat <- process_treat(X_t$treat)
            X_t <- assign.X.class(X_t)
            X_t$call <- NULL
            do.call("base.bal.tab", c(list(X_t), A[names(A) %nin% names(X_t)]), quote = TRUE)
        })
    }
    else {
        if (any(treat_vals(X$treat) == "All")) stop ("\"All\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
        balance.tables <- lapply(treat.combinations, function(t) {
            n <- length(X$treat)
            X_t <- X
            X_t$call <- NULL
            X_t <- subset_X(X_t, c(seq_len(n), which(X$treat == treat_vals(X$treat)[t[1]])))
            X_t$treat <- factor(rep(0:1, times = c(n, sum(X$treat == treat_vals(X$treat)[t[1]]))),
                                levels = c(0, 1), labels = c("All", t[1]))
            # if (is_not_null(X_t$weights)) X_t$weights[X_t$treat == "All",] <- 1 #Uncomment to compare each group to unweighted dist.
            X_t <- assign.X.class(X_t)
            do.call("base.bal.tab", c(list(X_t), A[names(A) %nin% names(X_t)]), quote = TRUE)
        })
    }
    
    for (i in seq_along(balance.tables)) {
        names(balance.tables)[i] <- paste(rev(treat.combinations[[i]]), collapse = " vs. ")
    }
    
    out[["Pair.Balance"]] <- balance.tables
    
    if ((multi.summary || !A$quick) && is_null(X$imp)) {
        out[["Balance.Across.Pairs"]] <- balance.summary(balance.tables, 
                                                         agg.funs = "max")
        
        out <- c(out, threshold.summary(compute = attr(out[["Pair.Balance"]][[1]][["Balance"]], "compute"),
                                        thresholds = attr(out[["Pair.Balance"]][[1]][["Balance"]], "thresholds"),
                                        no.adj = !attr(out[["Pair.Balance"]][[1]], "print.options")$disp.adj,
                                        balance.table = out[["Balance.Across.Pairs"]],
                                        weight.names = attr(out[["Pair.Balance"]][[1]], "print.options")$weight.names,
                                        agg.fun = "max"))
        
        out[["Observations"]] <- samplesize.multi(balance.tables, treat_names(X$treat), X$focal)
    }
    
    out[["call"]] <- X$call
    
    attr(out, "print.options") <- c(attr(out[["Pair.Balance"]][[1]], "print.options"),
                                    list(treat_vals_multi = treat_vals(X$treat),
                                         which.treat = which.treat,
                                         multi.summary = multi.summary,
                                         pairwise = pairwise))
    
    attr(out, "print.options")[["treat_names"]] <- NULL
    
    class(out) <- c("bal.tab.multi", "bal.tab")
    
    return(out)
    
}
base.bal.tab.msm <- function(X, which.time, msm.summary = getOption("cobalt_msm.summary"), ...) {
    #One vector of weights
    #treat.list should be a df/list of treatment vectors, one for each time period
    #covs.list should be a list of covariate data.frames, one for each time period; 
    #   should include all covs from previous time points, but no treatment statuses
    
    A <- list(...)
    
    # X$treat.list <- process_treat.list(X$treat)
    
    #Preparations
    
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    treat.types <- vapply(X$treat.list, get.treat.type, character(1L))
    
    if (missing(which.time)) {
        if (all_the_same(treat.types) && "multinomial" %nin% treat.types && is_null(X$imp)) which.time <- NA
        else which.time <- NULL
    }
    
    if (is_null(msm.summary)) {
        msm.summary <- is_not_null(which.time) && (anyNA(which.time) || !is_(which.time, c("character", "numeric")) ||
                                                       (is.numeric(which.time) && !any(which.time %in% seq_along(X$treat.list))) ||
                                                       (is.character(which.time) && !any(which.time %in% names(X$treat.list))))
    }
    
    #Setup output object
    # out.names <- c("Time.Balance", 
    #                "Balance.Across.Times", 
    #                "Observations", 
    #                "call")
    # out <- make_list(out.names)
    out <- list()
    
    out[["Time.Balance"]] <- make_list(length(X$covs.list))
    
    
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
        
        do.call("base.bal.tab", c(list(X_ti), A[names(A) %nin% names(X_ti)]), quote = TRUE)
    })
    
    if (length(names(X$treat.list)) == length(X$treat.list)) {
        names(out[["Time.Balance"]]) <- names(X$treat.list)
    }
    else names(out[["Time.Balance"]]) <- seq_along(X$treat.list)
    
    if (!(A$quick && !msm.summary) && all_the_same(treat.types) && "multinomial" %nin% treat.types && is_null(X$imp)) {
        out[["Balance.Across.Times"]] <- balance.summary(out[["Time.Balance"]],
                                                         agg.funs = "max",
                                                         include.times = TRUE)
        
        out <- c(out, threshold.summary(compute = attr(out[["Time.Balance"]][[1]][["Balance"]], "compute"),
                                        thresholds = attr(out[["Time.Balance"]][[1]][["Balance"]], "thresholds"),
                                        no.adj = !attr(out[["Time.Balance"]][[1]], "print.options")$disp.adj,
                                        balance.table = out[["Balance.Across.Times"]],
                                        weight.names = attr(out[["Time.Balance"]][[1]], "print.options")$weight.names,
                                        agg.fun = "max"))
        
        out[["Observations"]] <- lapply(out[["Time.Balance"]], function(x) x$Observations)
    }
    
    out[["call"]] <- X$call
    
    attr(out, "print.options") <- c(attr(out[["Time.Balance"]][[1]], "print.options"),
                                    list(which.time = which.time,
                                         msm.summary = msm.summary))
    
    class(out) <- c("bal.tab.msm", "bal.tab")
    
    return(out)
}
base.bal.tab.cluster <- function(X, which.cluster, cluster.summary = getOption("cobalt_cluster.summary"), cluster.fun = getOption("cobalt_cluster.fun", NULL), ...) {
    A <- list(...)
    
    #Preparations
    
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    if (is_null(A[["abs"]])) A[["abs"]] <- FALSE
    
    cluster <- factor(X$cluster)
    
    #Process cluster.summary
    if (missing(which.cluster)) {
        which.cluster <- NULL
    }
    
    if (is_null(cluster.summary)) {
        cluster.summary <- is_not_null(which.cluster) && anyNA(which.cluster)
    }
    
    all.agg.funs <- c("min", "mean", "max")
    agg.fun <- tolower(as.character(if_null_then(cluster.fun, A[["agg.fun"]], all.agg.funs)))
    agg.fun <- match_arg(agg.fun, all.agg.funs, several.ok = TRUE)
    
    #Setup output object
    # out.names <- c("Cluster.Balance", 
    #                "Balance.Across.Clusters", 
    #                "Observations", 
    #                "call")
    # out <- make_list(out.names)
    out <- list()
    
    #Get list of bal.tabs for each imputation
    out[["Cluster.Balance"]] <- lapply(levels(cluster), function(cl) {
        X_cl <- assign.X.class(subset_X(X, cluster == cl)) 
        X_cl$call <- NULL
        do.call("base.bal.tab", c(list(X_cl), A[names(A) %nin% names(X_cl)]), quote = TRUE)
    })
    
    names(out[["Cluster.Balance"]]) <- levels(cluster)
    
    #Create summary of lists
    
    if ((cluster.summary || !A$quick) && is_null(X$covs.list) && get.treat.type(X$treat) != "multinomial" && is_null(X$imp)) {
        out[["Balance.Across.Clusters"]] <- balance.summary(out[["Cluster.Balance"]], 
                                                            agg.funs = if_null_then(agg.fun, c("min", "mean", "max")))
        
        if (length(agg.fun) == 1) {
            out <- c(out, threshold.summary(compute = attr(out[["Cluster.Balance"]][[1]][["Balance"]], "compute"),
                                            thresholds = attr(out[["Cluster.Balance"]][[1]][["Balance"]], "thresholds"),
                                            no.adj = !attr(out[["Cluster.Balance"]][[1]], "print.options")$disp.adj,
                                            balance.table = out[["Balance.Across.Clusters"]],
                                            weight.names = attr(out[["Cluster.Balance"]][[1]], "print.options")$weight.names,
                                            agg.fun = agg.fun))
        }
        
        observations <- lapply(out[["Cluster.Balance"]], function(x) x[["Observations"]])
        
        out[["Observations"]] <- samplesize.across.clusters(observations)
    }
    
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- c(attr(out[["Cluster.Balance"]][[1]], "print.options"),
                                    list(which.cluster = which.cluster,
                                         cluster.summary = cluster.summary,
                                         cluster.fun = agg.fun))
    class(out) <- c("bal.tab.cluster", "bal.tab")
    
    return(out)
}

#NEEDS UPDATING with STATS
base.bal.tab.subclass <- function(X, type, int = FALSE, poly = 1, continuous, binary, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp = NULL, which.subclass = NA, subclass.summary = getOption("cobalt_subclass.summary"), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), disp.call = getOption("cobalt_disp.call", TRUE), abs = FALSE, quick = TRUE, ...) {
    #Preparations
    
    A <- clear_null(list(...))
    A$subset <- NULL
    
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
    
    subclass <- factor(X$subclass)
    
    if (missing(which.subclass)) {
        if (isTRUE(A[["disp.subclass"]])) which.subclass <- seq_len(nlevels(subclass))
        else which.subclass <- NA
    }
    
    if (is_null(subclass.summary)) {
        subclass.summary <- is_not_null(which.subclass) && 
            (anyNA(which.subclass) || !is.numeric(which.subclass) || 
                 (is.numeric(which.subclass) && !any(which.subclass %in% seq_len(nlevels(subclass)))))
    }
    
    no.adj <- FALSE
    
    if (is_null(X$s.weights)) {
        X$s.weights <- rep(1, length(X$treat))
    }
    
    disp <- process_disp(disp, ...)
    
    #Actions
    
    out <- list()
    
    C <- do.call("get.C2", c(X, A[names(A) %nin% names(X)], list(int = int, poly = poly)), quote = TRUE)
    co.names <- attr(C, "co.names")
    
    out[["Subclass.Balance"]] <- do.call("balance.table.subclass", 
                                         c(list(C, type = type, 
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
                                                quick = quick), A), quote = TRUE)
    
    if (subclass.summary || !quick) {
        if (type == "bin") {
            out[["Balance.Across.Subclass"]] <- do.call("balance.table", 
                                                        c(list(C, 
                                                               type = type, 
                                                               weights = data.frame(Adj = strata2weights(X$subclass, X$treat, X$estimand)), 
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
                                                               s.d.denom.list = X$s.d.denom.list), A), quote = TRUE)
        }
        else if (type == "cont") {
            out[["Balance.Across.Subclass"]] <- do.call("balance.table.across.subclass.cont", 
                                                        c(list(do.call("balance.table", c(list(C, 
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
                                                                                               s.d.denom.list = X$s.d.denom.list), A), quote = TRUE), 
                                                               balance.table.subclass.list = out[["Subclass.Balance"]], 
                                                               subclass.obs = out[["Observations"]], 
                                                               r.threshold = X$thresholds[["correlations"]]), A), quote = TRUE)
        }
    }
    
    #Reassign disp... and ...threshold based on balance table output
    compute <- attr(out[["Subclass.Balance"]], "compute")
    thresholds <- attr(out[["Subclass.Balance"]], "thresholds")
    disp <- attr(out[["Subclass.Balance"]], "disp")
    
    for (s in compute) {
        if (is_not_null(thresholds[[s]])) {
            out[[paste.("Balanced", s, "Subclass")]] <- setNames(do.call("data.frame", lapply(out[["Subclass.Balance"]], function(x) baltal(x[[STATS[[s]]$Threshold]]))),
                                                                 paste("Subclass", levels(X$subclass)))
            max.imbal.list <- lapply(out[["Subclass.Balance"]], function(x) {
                return(max.imbal(x[x[["Type"]] != "Distance", , drop = FALSE], 
                                 col.name = paste.(STATS[[s]]$bal.tab_column_prefix, "Adj"), 
                                 thresh.col.name = STATS[[s]]$Threshold, 
                                 abs_stat = STATS[[s]]$abs))
            } )
            out[[paste.("Max.Imbalance", s, "Subclass")]] <- data.frame(do.call("rbind", max.imbal.list), 
                                                                        row.names = paste("Subclass", levels(X$subclass)))
        }
    }
    
    out[["Observations"]] <- samplesize(treat = X$treat, 
                                        type = type,
                                        weights = NULL, 
                                        subclass = X$subclass,
                                        s.weights = X$s.weights, 
                                        method = X$method, 
                                        discarded = X$discarded)
    
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- list(thresholds = thresholds,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp=disp,
                                       compute = compute, 
                                       disp.adj = !no.adj, 
                                       which.subclass = which.subclass,
                                       subclass.summary = subclass.summary,
                                       disp.bal.tab = disp.bal.tab, 
                                       disp.call = disp.call,
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       treat_names = treat_vals(X$treat),
                                       type = type,
                                       co.names = co.names)
    class(out) <- c("bal.tab.subclass", "bal.tab")
    
    return(out)
}
base.bal.tab.subclass.binary <- function(X, ...) {
    base.bal.tab.subclass(X, type = "bin", ...)
}
base.bal.tab.subclass.cont <- function(X, ...) {
    stop("Subclasses are not yet compatible with continuous treatments.", call. = FALSE)
    base.bal.tab.subclass(X, type = "cont", ...)
}
