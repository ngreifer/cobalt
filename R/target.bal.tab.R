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
    
    if (X$s.d.denom %in% c("treated", "control")) {
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
        X$treat <- X$treat.list[[1]]
        X$covs <- X$covs.list[[1]]
        X$distance <- X$distance.list[[1]]
        X$addl <- X$addl.list[[1]]
    }
    
    out <- do.call("base.bal.tab.target", c(X, args),
                   quote = TRUE)
    return(out)
}