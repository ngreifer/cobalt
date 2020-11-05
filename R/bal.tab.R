bal.tab <- function(x, ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    .alls <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.all)), logical(1L))
    .nones <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.none)), logical(1L))
    if (any(c(.alls, .nones))) {
        .call[.alls] <- expression(NULL)
        .call[.nones] <- expression(NA)
        return(eval.parent(.call))
    }

    tryCatch(force(x), error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    if (!is_(x, "cobalt.processed.obj")) {
        x <- process_obj(x)
        bal.tab(x, ...)
    }
    else {
        UseMethod("bal.tab")
    }
}

#Point treatments
bal.tab.matchit <- function(x, method, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL

    #Initializing variables
    X <- do.call("x2base.matchit", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.ps <- function(x, stop.method, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.ps", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.mnps <- function(x, stop.method, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.mnps", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.ps.cont <- function(x, stop.method, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.ps.cont", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.Match <- function(x, formula = NULL, data = NULL, treat = NULL, covs = NULL, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base", c(list(x), args), quote = TRUE) 
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.formula <- function(x, data = NULL, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.formula", c(list(formula = x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.data.frame <- function(x, treat, data = NULL, weights = NULL, subclass = NULL, match.strata = NULL, method, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, focal = NULL, s.weights = NULL, estimand = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments

    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.data.frame", c(covs = list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.matrix <- function(x, treat, ...) {
    bal.tab.data.frame(as.data.frame.matrix(x), treat, ...)
}
bal.tab.CBPS <- function(x, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.weightit <- function(x, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = TRUE, ... ) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.ebalance <- function(x, ...) {
    bal.tab.Match(x, ...)
}
bal.tab.optmatch <- function(x, ...) {
    bal.tab.Match(x, ...)
}
bal.tab.cem.match <- function(x, data, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.cem.match", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.designmatch <- function(x, ...) {
    class(x) <- "designmatch"
    bal.tab.Match(x, ...)
}
bal.tab.mimids <- function(x, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.mimids", c(list(x), args), quote = TRUE)

    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.wimids <- function(x, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.wimids", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.sbwcau <- function(x, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.sbwcau", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}

#MSMs wth multiple time points
bal.tab.formula.list <- function(x, data = NULL, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    X <- do.call("x2base.formula.list", c(list(formula.list = x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.data.frame.list <- function(x, treat.list = NULL, data = NULL, weights = NULL, stats, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, method, continuous, binary, s.d.denom, thresholds = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, estimand = "ATE", abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    X <- do.call("x2base.data.frame.list", c(list(covs.list = x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.iptw <- function(x, stop.method, stats, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.iptw", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.CBMSM <- function(x, stats, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.CBMSM", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.weightitMSM <- function(x, stats, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = TRUE, ... ) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.weightitMSM", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}

#default method
bal.tab.default <- function(x, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.default", c(list(obj = x), args),
                 quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
