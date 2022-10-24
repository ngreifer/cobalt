bal.tab <- function(x, ...) {
    
    .call <- match.call()
    
    tryCatch(force(x), error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Replace .all and .none with NULL and NA respectively
    if (!is_(x, "cobalt.processed.obj")) {
        
        .alls <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.all)), logical(1L))
        .nones <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.none)), logical(1L))
        if (any(c(.alls, .nones))) {
            .call[.alls] <- expression(NULL)
            .call[.nones] <- expression(NA)
        }
        .call[["x"]] <- process_obj(x)
        
        return(eval.parent(.call))
    }
    else {
        UseMethod("bal.tab")
    }
}

#Point treatments
bal.tab.matchit <-    function(x,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               method, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.matchit", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.ps <-         function(x, stop.method,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.ps", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.mnps <-       function(x, stop.method,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.mnps", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.ps.cont <-    function(x,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.ps.cont", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.Match <-      function(x, formula = NULL, data = NULL, treat = NULL, covs = NULL,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base", c(list(x), args), quote = TRUE) 
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.formula <-    function(x, data = NULL,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               subclass = NULL, match.strata = NULL, method, estimand = NULL, focal = NULL, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.formula", c(list(formula = x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.data.frame <- function(x, treat,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               subclass = NULL, match.strata = NULL, method, estimand = NULL, focal = NULL, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments

    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.data.frame", c(covs = list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
bal.tab.matrix <- bal.tab.data.frame
bal.tab.CBPS <-       function(x,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.weightit <-   function(x,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.ebalance <- bal.tab.Match
bal.tab.optmatch <- bal.tab.Match
bal.tab.cem.match <-  function(x, data,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.cem.match", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.designmatch <- bal.tab.Match
bal.tab.mimids <-     function(x,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
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
bal.tab.wimids <-     function(x,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.wimids", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.sbwcau <-     function(x,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.sbwcau", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}

#MSMs wth multiple time points
bal.tab.formula.list <- function(x,
                                 stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                                 ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    X <- do.call("x2base.formula.list", c(list(formula.list = x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.data.frame.list <- function(x, treat.list, 
                                    stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                                    ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    X <- do.call("x2base.data.frame.list", c(list(covs.list = x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.iptw <- function(x,
                         stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                         stop.method, ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.iptw", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.CBMSM <- function(x,
                          stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                          ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.CBMSM", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}
bal.tab.weightitMSM <- function(x, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                                ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.weightitMSM", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}

#default method
bal.tab.default <- function(x, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                            ...) {
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Adjustments to arguments
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.default", c(list(obj = x), args),
                 quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    
    return(out)
}
