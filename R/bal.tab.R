bal.tab <- function(...) {
    if (...length() == 0L) stop("No arguments were supplied.", call. = FALSE)
    .obj <- ...elt(1)
    .obj <- is.designmatch(.obj)
    .obj <- is.time.list(.obj)
    UseMethod("bal.tab", .obj)
}

#Point treatments
bal.tab.matchit <- function(m, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
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
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.ps <- function(ps, stop.method, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
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
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.mnps <- function(mnps, stop.method, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, pairwise = TRUE, focal = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
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
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.ps.cont <- function(ps.cont, stop.method, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, r.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
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
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.Match <- function(M, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
  
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base", c(list(M), args), quote = TRUE) 
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.formula <- function(formula, data = NULL, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.formula", c(list(formula = formula), args), quote = TRUE)
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.data.frame <- function(covs, treat, data = NULL, weights = NULL, distance = NULL, subclass = NULL, match.strata = NULL, method, int = FALSE, poly = 1, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, focal = NULL, s.weights = NULL, estimand = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-c(1, length(formals()))])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    X <- do.call("x2base.data.frame", c(list(covs), args), quote = TRUE)
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.CBPS <- function(cbps, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, pairwise = TRUE, focal = NULL, s.weights = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base", c(list(cbps), args), quote = TRUE)
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.weightit <- function(weightit, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL,  continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, focal = NULL, abs = FALSE, subset = NULL, quick = FALSE, ... ) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base", c(list(weightit), args), quote = TRUE)
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.ebalance <- function(ebal, ...) {
    bal.tab.Match(ebal, ...)
}
bal.tab.ebalance.trim <- bal.tab.ebalance
bal.tab.optmatch <- function(optmatch, ...) {
    bal.tab.Match(optmatch, ...)
}
bal.tab.designmatch <- function(dm, ...) {
    class(dm) <- "designmatch"
    bal.tab.Match(dm, ...)
}

#MSMs wth multiple time points
bal.tab.formula.list <- function(formula.list, data = NULL, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    X <- do.call("x2base.formula.list", c(list(formula.list = formula.list), args), quote = TRUE)
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.data.frame.list <- function(covs.list, treat.list = NULL, data = NULL, weights = NULL, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, method, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, pairwise = TRUE, s.weights = NULL, estimand = "ATE", abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }

    X <- do.call("x2base.data.frame.list", c(list(covs.list = covs.list), args), quote = TRUE)
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.iptw <- function(iptw, stop.method, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
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
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    return(out)
}
bal.tab.CBMSM <- bal.tab.CBPS

#default method
bal.tab.default <- function(obj, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
    for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
    
    blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                args[[arg.name]] <- NULL
            }
        }
    }
    
    #Initializing variables
    X <- do.call("x2base.default", c(list(obj = obj), args),
                 quote = TRUE)
    
    args <- args[names(args) %nin% names(X)]
    
    out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                   quote = TRUE)
    
    return(out)
}
