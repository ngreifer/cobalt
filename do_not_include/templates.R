#Templates
.x2base.template.point <- function(OBJ, ...) {
    A <- list(...)
    
    #Process OBJ
    
    #Process data and get imp
    OBJ.data <- OBJ$data
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", data = data, missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat
    treat <- process_treat(OBJ$treat, data = list(data, OBJ.data))
    
    #Process covs
    covs <- OBJ$covs
    
    #Get estimand
    estimand <- ps$estimand
    
    #Get method
    method <- "method"
    
    #Process addl 
    addl <- data.frame.process("addl", A[["addl"]], treat, covs, data, OBJ.data)
    
    #Process distance
    distance <- data.frame.process("distance", A[["distance"]], treat, covs, data, OBJ.data)
    if (any(is.finite(OBJ$distance))) {
        if (is_not_null(distance)) distance <- cbind(distance, distance = OBJ$distance)
        else distance <- data.frame(distance = OBJ$distance)
    }
    
    #Process subclass
    if (is_not_null(subclass <- OBJ$subclass)) {
        subclass <- factor(subclass)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- OBJ$match.strata)) {
        match.strata <- vector.process(match.strata, 
                                       data = list(data, OBJ.data),
                                       name = "match.strata", 
                                       which = "matching strata",
                                       missing.okay = FALSE)
        weights <- data.frame(weights = strata2weights(match.strata,
                                                       treat = treat))
    }
    
    #Process weights
    if (is_not_null(weights <- OBJ$weights)) {
        weights <- data.frame(weights = OBJ$weights)
        if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) warning("Negative weights found.", call. = FALSE)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        s.weights <- vector.process(s.weights, 
                                    data = list(data, OBJ.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, OBJ.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("covs", "weights", "distance", "addl"),
                       imp = imp,
                       original.call.to = "FUN()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        focal <- focal
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat, focal = focal)
    
    #Missing values warning
    if (any(c(anyNA(covs), anyNA(addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- ps$parameters
    
    #Process output
    X <- initialize_X()
    
    for (i in names(X)) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[names(X)], names(X))
    
    class(X) <- "binary"
    
    return(X)
}
.x2base.template.msm <- function(OBJ, ...) {
    A <- list(...)
    
    #Process OBJ
    
    #Process data and get imp
    OBJ.data <- OBJ$data
    imp <- A$imp
    if (is_not_null(data <- A$data)) {
        if (is_(data, "mids")) {
            data <- imp.complete(data)
            if (is_null(imp)) imp <- data[[".imp"]]
        }
        else if (!is_(data, "data.frame"))
        {
            # warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
            data <- NULL
        }
    }
    
    #Process imp
    if (is_not_null(imp)) {
        imp <- vector.process(imp, "imp", "imputation identifiers", data = data, missing.okay = FALSE)
        imp <- factor(imp)
    }
    
    #Process treat.list
    treat.list <- process_treat.list(OBJ$treat.list, data = list(data, OBJ.data))
    
    #Process covs.list
    all.covs <- unique(unlist(lapply(OBJ$covs.list, names)))
    covs.list <- lapply(OBJ$covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
    
    #Get estimand
    estimand <- A$estimand
    
    #Get method
    method <- "weighting"
    
    #Process addl.list 
    ntimes <- iptw$nFits
    addl.list <- list.process("addl", A[["addl"]], ntimes, 
                              "the original call to iptw()",
                              treat.list,
                              covs.list,
                              data,
                              OBJ.data)
    
    #Process distance
    distance.list <- list.process("distance", A[["distance"]], ntimes, 
                              "the original call to iptw()",
                              treat.list,
                              covs.list,
                              data,
                              OBJ.data)
    if (any(is.finite(OBJ$distance))) {
        if (is_not_null(distance)) distance <- cbind(distance, distance = OBJ$distance)
        else distance <- data.frame(distance = OBJ$distance)
    }
    
    #Process subclass
    if (is_not_null(subclass <- OBJ$subclass)) {
        subclass <- factor(subclass)
    }
    
    #Process match.strata
    if (is_not_null(match.strata <- OBJ$match.strata)) {
        match.strata <- vector.process(match.strata, 
                                       data = list(data, OBJ.data),
                                       name = "match.strata", 
                                       which = "matching strata",
                                       missing.okay = FALSE)
        weights <- data.frame(weights = strata2weights(match.strata,
                                                       treat = treat))
    }
    
    #Process weights
    if (is_not_null(weights <- OBJ$weights)) {
        weights <- data.frame(weights = OBJ$weights)
        if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) warning("Negative weights found.", call. = FALSE)
    }
    
    #Process s.weights
    if (is_not_null(s.weights <- A$s.weights)) {
        s.weights <- vector.process(s.weights, 
                                    data = list(data, OBJ.data),
                                    name = "s.weights", 
                                    which = "sampling weights",
                                    missing.okay = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster <- A$cluster)) {
        cluster <- vector.process(cluster, 
                                  data = list(data, OBJ.data),
                                  name = "cluster", 
                                  which = "cluster membership",
                                  missing.okay = FALSE)
        cluster <- factor(cluster)
    }
    
    #Process subset
    if (is_not_null(subset <- A$subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (anyNA(subset)) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    #Process discarded
    
    #Process length
    length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                       data.frames = c("weights"),
                       lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                       imp = imp,
                       original.call.to = "FUN()")
    
    #Process focal
    if (is_not_null(focal <- A$focal)) {
        focal <- focal
    }
    
    #Get s.d.denom
    s.d.denom <- get.s.d.denom(A$s.d.denom, estimand = estimand, weights = weights, treat = treat.list[[1]], focal = focal)
    
    #Missing values warning
    if (any(c(any(vapply(covs.list, anyNA, logical(1L))), any(vapply(addl.list, anyNA, logical(1L)))))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    #Get call
    call <- OBJ$call
    
    #Process output
    X <- initialize_X_msm()
    
    for (i in names(X)) {
        X[[i]] <- get0(i, inherits = FALSE)
    }
    
    X <- subset_X(X, subset)
    X <- setNames(X[names(X)], names(X))
    
    class(X) <- "binary"
    
    return(X)
}
.bal.tab.template <- function(OBJ, int = FALSE, distance = NULL, addl = NULL, data = NULL,  continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, which.treat = NA, pairwise = TRUE, focal = NULL, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, abs = FALSE, subset = NULL, quick = FALSE, ... ) {
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
    
    if (any(class(OBJ) == "MSM")) {
        if (is_not_null(cluster)) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
        if (is_not_null(imp)) stop("Multiply imputed data is not yet supported with longitudinal treatments.", call. = FALSE)
        if (is_not_null(args$addl.list)) addl <- args$addl.list
        
        #Initializing variables
        X <- do.call("x2base.OBJ_MSM", c(list(OBJ), args), quote = TRUE)
        
        args <- args[names(args) %nin% attr(X, "X.names")]
        
        X <- setNames(X[attr(X, "X.names")], attr(X, "X.names"))
        
        out <- do.call("base.bal.tab.msm", c(X, args),
                       quote = TRUE)
    }
    else {
        #Initializing variables
        X <- do.call("x2base.OBJ", c(list(OBJ), args), quote = TRUE)
        
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
        else if (is.factor(X[["treat"]]) || is.character(X[["treat"]])) {
            
            out <- do.call("base.bal.tab.multi", c(X, args),
                           quote = TRUE)
        }
        else if (is.numeric(X[["treat"]])) {
            
            out <- do.call("base.bal.tab.cont", c(X, args),
                           quote = TRUE)
        }
        else stop("Something went wrong. Contact the maintainer.", call. = FALSE)
    }
    
    return(out)
}
