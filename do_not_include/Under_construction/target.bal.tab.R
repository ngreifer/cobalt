target.bal.tab <- function(obj, ..., method, int = FALSE, addl = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, which.treat = NA, target.summary = TRUE, s.weights = NULL, estimand = NULL, abs = FALSE, subset = NULL, quick = FALSE) {
    tryCatch(identity(obj), error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals())[vapply(formals(), function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    obj <- is.designmatch(obj)
    obj <- is.time.list(obj)
    if (any(class(obj) == "time.list")) {
        if (all(sapply(obj, is.formula))) class(obj) <- "formula.list"
        else if (all(sapply(obj, is.data.frame))) class(obj) <- "data.frame.list"
        else stop("If obj is a list, it must be a list of formulas specifying the treatment/covariate relationships at each time point or a list of data frames containing covariates to be assessed at each time point.", call. = FALSE)
    }
    
    #Initializing variables
    X <- do.call("x2base", c(list(obj), args), quote = TRUE)
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    if (is_not_null(X$s.d.denom) && X$s.d.denom %in% c("treated", "control")) {
        stop("When the estimand is the ATT or ATC, use bal.tab to assess balance between non-focal groups and the focal (target) group.", call. = FALSE)
    }
    if (is_not_null(X$subclass)) {
        stop("Target balance assessment is not yet supported with subclassification.", call. = FALSE)
    }
    if (is_not_null(X$cluster)) {
        stop("Target balance assessment is not yet supported with clustered data.", call. = FALSE)
    }
    if (is_not_null(X$imp)) {
        stop("Target balance assessment is not yet supported with multiply imputed data.", call. = FALSE)
    }
    if (is_not_null(X$treat.list) && is_not_null(X$covs.list)) {
        warning("Target balance assessment will only occur for the baseline covariates.", call. = FALSE)
        X$treat <- X$treat.list[[1]]
        X$covs <- X$covs.list[[1]]
        X$distance <- X$distance.list[[1]]
        X$addl <- X$addl.list[[1]]
    }
    
    out <- do.call("base.bal.tab.target", c(X, args),
                   quote = TRUE)
    return(out)
}

base.bal.tab.target <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), which.treat = NA, target.summary = getOption("cobalt_target.summary", TRUE), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
    #Preparations
    args <- list(...)
    
    if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
    if (is_not_null(v.threshold)) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (is_not_null(ks.threshold)) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (is_null(weights) && is_null(subclass)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (is_null(s.weights)) {
        s.weights <- rep(1, length(treat))
    }
    if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
    
    #Create new Target group
    target.name <- "Target"
    n <- length(treat)
    
    if (isTRUE(get.treat.type(treat) == "continuous") || (is.numeric(treat) && !is_binary(treat))) {#if continuous treatment
        covs <- data.frame(treat = treat, covs)
        treat <- factor(rep(c("All", target.name), each = n))
        target.summary <- FALSE
        which.treat <- NULL
        needs.summary <- FALSE
        treat_names <- unique.treat <- "All"
    }
    else {
        if (is.factor(treat) || is.character(treat)) {
            if (is.factor(treat)) treat_names <- unique.treat <- levels(treat)
            else treat_names <- unique.treat <- unique(treat, nmax = n - 1)
        }
        else {
            treat_names <- c("Control", "Treated")
            unique.treat <- sort(unique(treat, nmax = 2))
        }
        names(treat_names) <- unique.treat
        
        treat <- factor(c(treat_names[as.character(treat)], rep(target.name, n)))
        needs.summary <- TRUE
    }
    
    covs <- rbind(covs, covs)
    if (is_not_null(weights)) weights <- rbind(weights, as.data.frame(array(1, dim = dim(weights), 
                                                                            dimnames = dimnames(weights))))
    distance <- rbind(distance, distance)
    addl <- rbind(addl, addl)
    s.weights <- c(s.weights, s.weights)
    if (is_not_null(discarded)) discarded <- c(discarded, rep(FALSE, length(discarded)))
    s.d.denom <- "treated"
    
    treat.target.combinations <- lapply(treat_names, function(x) c(x, target.name))
    
    if (is_not_null(cluster)) {
        stop("Clusters are not yet supported with target balance assessment.", call. = FALSE)
    }
    else if (is_not_null(subclass)) {
        stop("Subclassification is not yet supported with target balance assessment.", call. = FALSE)
    }
    else {
        #Setup output object
        out.names <- c("Target.Balance", 
                       "Balance.Across.Treatments", 
                       "Observations", 
                       "call")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        
        if (any(treat_names == "Target")) stop ("\"Target\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
        args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
        balance.tables <- lapply(treat.target.combinations, function(t) do.call("base.bal.tab.binary", c(list(weights = weights[treat %in% t, , drop = FALSE], treat = factor(treat[treat %in% t], t), distance = distance[treat %in% t, , drop = FALSE], subclass = subclass[treat %in% t], covs = covs[treat %in% t, , drop = FALSE], call = NULL, int = int, poly = poly, addl = addl[treat %in% t, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster[treat %in% t], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[treat %in% t], discarded = discarded[treat %in% t], quick = quick), args), quote = TRUE))
        
        for (i in seq_along(balance.tables)) {
            names(balance.tables)[i] <- paste(treat.target.combinations[[i]], collapse = " vs. ")
            balance.tables[[i]][["Observations"]][[2]] <- NULL
        }
        
        out[["Target.Balance"]] <- balance.tables
        
        out[["Observations"]] <- samplesize.target(balance.tables, treat_names, target.name) 
        
        if (needs.summary && (target.summary || !quick)) {
            out[["Balance.Across.Treatments"]] <- balance.table.target.summary(balance.tables, 
                                                                               weight.names = names(weights),
                                                                               m.threshold = m.threshold,
                                                                               v.threshold = v.threshold,
                                                                               ks.threshold = ks.threshold,
                                                                               no.adj = no.adj,
                                                                               quick = quick,
                                                                               types = NULL)
        }
        
        out[["call"]] <- call
        
        attr(out, "print.options") <- list(m.threshold=m.threshold,
                                           v.threshold=v.threshold,
                                           ks.threshold=ks.threshold,
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.adj=!no.adj, 
                                           which.cluster=which.cluster,
                                           cluster.summary=cluster.summary,
                                           cluster.fun = cluster.fun,
                                           abs = abs,
                                           continuous = continuous,
                                           binary = binary,
                                           quick = quick,
                                           disp.means=disp.means, 
                                           disp.sds = disp.sds,
                                           disp.v.ratio=disp.v.ratio, 
                                           disp.ks=disp.ks,
                                           disp.bal.tab = disp.bal.tab,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights),
                                           treat_names = treat_names,
                                           target.name = target.name,
                                           which.treat = which.treat,
                                           target.summary = target.summary,
                                           co.names = attr(out[["Target.Balance"]][[1]], "print.options")[["co.names"]])
        
        class(out) <- c("bal.tab.target", "bal.tab")
    }
    return(out)
    
}
