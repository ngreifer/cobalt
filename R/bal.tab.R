bal.tab <- function(...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    .alls <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.all)), logical(1L))
    .nones <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.none)), logical(1L))
    if (any(c(.alls, .nones))) {
        .call[.alls] <- expression(NULL)
        .call[.nones] <- expression(NA)
        return(eval.parent(.call))
    }
    
    if (...length() == 0L) stop("No arguments were supplied.", call. = FALSE)
    .obj <- ...elt(1)

    .obj <- process_designmatch(.obj)
    .obj <- process_time.list(.obj)
    .obj <- process_cem.match.list(.obj)
    
    UseMethod("bal.tab", .obj)
}

#Point treatments
bal.tab.matchit <- function(m, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.ps <- function(ps, stop.method, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.mnps <- function(mnps, stop.method, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, focal = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.ps.cont <- function(ps.cont, stop.method, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.Match <- function(M, formula = NULL, data = NULL, treat = NULL, covs = NULL, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.formula <- function(formula, data = NULL, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.data.frame <- function(covs, treat, data = NULL, weights = NULL, subclass = NULL, match.strata = NULL, method, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, focal = NULL, s.weights = NULL, estimand = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.numeric <- function(treat, covs, ...) {
    bal.tab.data.frame(as.data.frame(covs), treat = treat, ...)
}
bal.tab.factor <- bal.tab.character <- bal.tab.logical <- bal.tab.numeric
bal.tab.matrix <- function(covs, treat, ...) {
    bal.tab.data.frame(as.data.frame(covs), treat, ...)
}
bal.tab.CBPS <- function(cbps, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, focal = NULL, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.weightit <- function(weightit, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, focal = NULL, abs = FALSE, subset = NULL, quick = TRUE, ... ) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
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
bal.tab.cem.match <- function(cem.match, data, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    X <- do.call("x2base.cem.match", c(list(cem.match), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.designmatch <- function(dm, ...) {
    class(dm) <- "designmatch"
    bal.tab.Match(dm, ...)
}
bal.tab.mimids <- function(mimids, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    X <- do.call("x2base.mimids", c(list(mimids), args), quote = TRUE)

    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.wimids <- function(wimids, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, pairwise = TRUE, focal = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    X <- do.call("x2base.wimids", c(list(wimids), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.sbwcau <- function(sbwcau, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    X <- do.call("x2base.sbwcau", c(list(sbwcau), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}

#MSMs wth multiple time points
bal.tab.formula.list <- function(formula.list, data = NULL, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.data.frame.list <- function(covs.list, treat.list = NULL, data = NULL, weights = NULL, stats, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, method, continuous, binary, s.d.denom, thresholds = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, estimand = "ATE", abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.iptw <- function(iptw, stop.method, stats, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    X <- do.call("x2base.iptw", c(list(iptw), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.CBMSM <- bal.tab.CBPS

#default method
bal.tab.default <- function(obj, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
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
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call(base.bal.tab, c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
