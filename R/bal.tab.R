bal.tab <- function(...) {
    A <- list(...)
    if (is_null(A)) stop("No arguments were supplied.", call. = FALSE)
    A[[1]] <- is.designmatch(A[[1]])
    A[[1]] <- is.time.list(A[[1]])
    UseMethod("bal.tab", A[[1]])
}

#Point treatments
bal.tab.matchit <- function(m, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.matchit", c(list(m), args), quote = TRUE)
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab", c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.ps <- function(ps, stop.method, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    args <- as.list(environment())[-(1:2)]
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.ps", c(list(ps), args), quote = TRUE)
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab", c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.mnps <- function(mnps, stop.method, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    args <- as.list(environment())[-(1:2)]
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.mnps", c(list(mnps), args), quote = TRUE)
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab.multi", c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.ps.cont <- function(ps.cont, stop.method, int = FALSE, distance = NULL, addl = NULL, data = NULL, r.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    args <- as.list(environment())[-(1:2)]
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.ps.cont", c(list(ps.cont), args), quote = TRUE)
    
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab.cont", c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.Match <- function(M, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.Match", c(list(M), args), quote = TRUE) 
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab", c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.formula <- function(formula, data = NULL, ...) {
    
    args <- c(as.list(environment()), list(...))[-(1:2)]
    args[["covs"]] <- NULL
    args[["treat"]] <- NULL
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    t.c <- get.covs.and.treat.from.formula(formula, data)
    covs <- t.c[["reported.covs"]]
    treat <- t.c[["treat"]]
    
    out <- do.call("bal.tab", c(list(covs = covs, treat = treat), args), quote = TRUE)
    
    return(out)
}
bal.tab.data.frame <- function(covs, treat, data = NULL, weights = NULL, distance = NULL, subclass = NULL, match.strata = NULL, method, int = FALSE, addl = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = TRUE, s.weights = NULL, estimand = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    X <- do.call("x2base.data.frame", c(list(covs), args), quote = TRUE)
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    if (is_not_null(X$imp)) {
        
        out <- do.call("base.bal.tab.imp", c(X, args),
                       quote = TRUE)
    }
    else if (is_binary(X$treat)) {
        out <- do.call("base.bal.tab", c(X, args),
                       quote = TRUE)
    }
    else if (is.factor(X$treat) || is.character(X$treat)) {
        
        out <- do.call("base.bal.tab.multi", c(X, args),
                       quote = TRUE)
    }
    else if (is.numeric(X$treat)) {
        
        out <- do.call("base.bal.tab.cont", c(X, args),
                       quote = TRUE)
    }
    else stop("Something went wrong. Contact the maintainer.", call. = FALSE)
    
    return(out)
}
bal.tab.CBPS <- function(cbps, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    if (any(class(cbps) == "CBMSM")) {
        if (is_not_null(cluster)) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
        
        #Initializing variables
        X <- do.call("x2base.CBMSM", c(list(cbps), args), quote = TRUE)
        
        
        args <- args[names(args) %nin% attr(X, "X.names")]
        
        X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
        
        out <- do.call("base.bal.tab.msm", c(X, args),
                       quote = TRUE)
    }
    else {
        #Initializing variables
        X <- do.call("x2base.CBPS", c(list(cbps), args), quote = TRUE)
        
        
        args <- args[names(args) %nin% attr(X, "X.names")]
        
        X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
        
        if (any(class(cbps) == "CBPSContinuous")) {
            
            out <- do.call("base.bal.tab.cont", c(X, args),
                           quote = TRUE)
        }
        else if (!is_binary(X$treat)) {
            
            out <- do.call("base.bal.tab.multi", c(X, args),
                           quote = TRUE)
        }
        else {
            out <- do.call("base.bal.tab", c(X, args),
                           quote = TRUE)
        }
    }
    return(out)
}
bal.tab.ebalance <- function(ebal, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.ebalance", c(list(ebal), args), quote = TRUE)
    
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab", c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.ebalance.trim <- bal.tab.ebalance
bal.tab.optmatch <- function(optmatch, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.optmatch", c(list(optmatch), args), quote = TRUE)
    
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab", c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.weightit <- function(weightit, int = FALSE, distance = NULL, addl = NULL, data = NULL,  continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, which.treat = NA, pairwise = TRUE, focal = NULL, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ... ) {
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    if (any(class(weightit) == "weightitMSM")) {
        if (is_not_null(cluster)) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
        if (is_not_null(imp)) stop("Multiply imputed data is not yet supported with longitudinal treatments.", call. = FALSE)
        if (is_not_null(args$addl.list)) addl <- args$addl.list
        
        #Initializing variables
        X <- do.call("x2base.weightitMSM", c(list(weightit), args), quote = TRUE)
        
        args <- args[names(args) %nin% attr(X, "X.names")]
        
        X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
        
        out <- do.call("base.bal.tab.msm", c(X, args),
                       quote = TRUE)
    }
    else {
        #Initializing variables
        X <- do.call("x2base.weightit", c(list(weightit), args), quote = TRUE)
        
        args <- args[names(args) %nin% attr(X, "X.names")]
        
        X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
        
        if (is_not_null(X$imp)) {
            
            out <- do.call("base.bal.tab.imp", c(X, args),
                           quote = TRUE)
        }
        else if (isTRUE(weightit$treat.type == "binary") || isTRUE(attr(weightit$treat, "treat.type") == "binary")) {
            out <- do.call("base.bal.tab", c(X, args),
                           quote = TRUE)
        }
        else if (isTRUE(weightit$treat.type == "multinomial") || isTRUE(attr(weightit$treat, "treat.type") == "multinomial")) {
            out <- do.call("base.bal.tab.multi", c(X, args),
                           quote = TRUE)
        }
        else if (isTRUE(weightit$treat.type == "continuous") || isTRUE(attr(weightit$treat, "treat.type") == "continuous")) {
            
            out <- do.call("base.bal.tab.cont", c(X, args),
                           quote = TRUE)
        }
        else stop("Something went wrong. Contact the maintainer.", call. = FALSE)
    }
    
    return(out)
}
bal.tab.designmatch <- function(dm, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.designmatch", c(list(dm), args), quote = TRUE)
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab", c(X, args),
                   quote = TRUE)
    return(out)
}

#MSMs wth multiple time points
bal.tab.time.list <- function(time.list, data, treat.list = NULL, weights = NULL, int = FALSE, distance.list = NULL, addl.list = NULL, method, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, pairwise = TRUE, which.treat = NA, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, s.weights = NULL, estimand = "ATE", abs = FALSE, subset = NULL, quick = FALSE, ...) {
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    if (is_not_null(args$cluster)) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
    if (is_not_null(args$imp)) stop("Multiply imputed data is not yet supported with longitudinal treatments.", call. = FALSE)
    
    if (all(vapply(time.list, is.formula, logical(1L)))) {
        X <- do.call("x2base.formula.list", c(list(formula.list = time.list), args), quote = TRUE)
    }
    else if (all(vapply(time.list, is.data.frame, logical(1L)))) {
        X <- do.call("x2base.data.frame.list", c(list(covs.list = time.list), args), quote = TRUE)
    }
    else {
        stop("If the first argument is a list, it must be a list of formulas specifying the treatment/covariate relationships at each time point or a list of data frames containing covariates to be assessed at each time point.", call. = FALSE)
    }
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab.msm", c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.iptw <- function(iptw, stop.method, int = FALSE, distance.list = NULL, addl.list = NULL, data = NULL, continuous = getOption("cobalt_cont", "std"), binary = getOption("cobalt_bin", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, pairwise = TRUE, which.treat = NA, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    args <- as.list(environment())[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    if (is_not_null(args$cluster)) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
    if (is_not_null(args$imp)) stop("Multiply imputed data sets are not yet supported with longitudinal treatments.", call. = FALSE)
    
    #Initializing variables
    X <- do.call("x2base.iptw", c(list(iptw), args), quote = TRUE)
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    out <- do.call("base.bal.tab.msm", c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.CBMSM <- bal.tab.CBPS

#default method
bal.tab.default <- function(obj, ...) {
    args <- c(as.list(environment()), list(...))[-(1:2)]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match.arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.default", c(list(obj = obj),
                                     args),
                 quote = TRUE)
    
    args <- args[names(args) %nin% attr(X, "X.names")]
    
    X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
    
    if (is_not_null(X$treat.list) && is_not_null(X$covs.list)) { #MSM
        if (is_not_null(X$cluster)) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
        if (is_not_null(X$imp)) stop("Multiply imputed data sets are not yet supported with longitudinal treatments.", call. = FALSE)
        
        out <- do.call("base.bal.tab.msm", c(X, args),
                       quote = TRUE)
    }
    else if (is_not_null(X$imp)) {
        out <- do.call("base.bal.tab.imp", c(X, args),
                       quote = TRUE)
    }
    else if (!is_binary(X$treat) && is.numeric(X$treat)) {
        out <- do.call("base.bal.tab.cont", c(X, args),
                       quote = TRUE)
    }
    else if (!is_binary(X$treat) && (is.factor(X$treat) || is.character(X$treat))) {
        out <- do.call("base.bal.tab.multi", c(X, args),
                       quote = TRUE)
    }
    else {
        out <- do.call("base.bal.tab", c(X, args),
                       quote = TRUE)
    }
    
    return(out)
}
