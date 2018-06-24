bal.tab.default <- function(obj, formula = NULL, data = NULL, treat = NULL, covs = NULL, weights = NULL, distance = NULL, subclass = NULL, match.strata = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = TRUE, s.weights = NULL, estimand = NULL, abs = FALSE, subset = NULL, quick = FALSE, ...) {
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    blank.args <- sapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)))
    if (any(blank.args)) {
        for (arg.name in names(blank.args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base.default(obj, covs = covs,
                        treat = treat,
                        formula = formula,
                        data = data,
                        weights = weights,
                        distance = distance,
                        subclass = subclass,
                        match.strata = match.strata,
                        addl = addl,
                        s.d.denom = s.d.denom,
                        method = method,
                        cluster = cluster,
                        estimand = estimand,
                        imp = imp,
                        s.weights = s.weights,
                        focal = focal,
                        subset = subset)
    
    out <- base.bal.tab(weights=X$weights, 
                        treat=X$treat, 
                        distance=X$distance, 
                        subclass=X$subclass, 
                        covs=X$covs, 
                        call=X$call, 
                        int=int, 
                        addl=X$addl, 
                        continuous=continuous, 
                        binary=binary, 
                        s.d.denom=X$s.d.denom, 
                        m.threshold=m.threshold, 
                        v.threshold=v.threshold, 
                        ks.threshold=ks.threshold, 
                        imbalanced.only = imbalanced.only,
                        un=un, 
                        disp.means=disp.means, 
                        disp.v.ratio=disp.v.ratio, 
                        disp.ks=disp.ks, 
                        disp.subclass=disp.subclass,                                 
                        disp.bal.tab = disp.bal.tab,
                        method=X$method, 
                        cluster = X$cluster, 
                        which.cluster = which.cluster, 
                        cluster.summary = cluster.summary, 
                        abs = abs,
                        discarded = X$discarded, 
                        quick = quick)
    return(out)
}

x2base.default <- function(obj, ...) {
    #treat
    #covs
    #data
    #formula
    #weights
    #distance
    #subclass
    #match.strata
    #addl
    #s.d.denom
    #method
    #cluster
    #estimand
    #imp
    
    A <- list(...)
    X <- list(covs = NA,
              weights = NA,
              treat = NA,
              distance = NA,
              subclass = NA,
              match.strata = NA,
              addl = NA,
              method = NA,
              call = NA,
              cluster = NA,
              imp = NA,
              s.weights = NA)
    in.obj <- setNames(vector("list", 10), 
                       c("treat", "covs", "formula", "data",
                       ))
    
    if (is_not_null(obj[["treat"]])) {
        
        in.obj[["treat"]] <- TRUE
    }
    else in.obj[["treat"]] <- FALSE
    
    if (is_not_null(obj[["covs"]])) {
        
        in.obj[["covs"]] <- TRUE
    }
    else in.obj[["covs"]] <- FALSE
    
    if (is_not_null(obj[["formula"]])) {
        
        in.obj[["formula"]] <- TRUE
    }
    else in.obj[["formula"]] <- FALSE
    
    if (is_not_null(obj[["data"]])) {
        
        in.obj[["data"]] <- TRUE
    }
    else in.obj[["data"]] <- FALSE
    
    if (is_not_null(obj[["weights"]])) {
        
        in.obj[["weights"]] <- TRUE
    }
    else in.obj[["weights"]] <- FALSE
    
    if (is_not_null(obj[["distance"]]) || is_not_null(obj[["ps"]])) {
        
        in.obj[["distance"]] <- TRUE
    }
    else in.obj[["distance"]] <- FALSE
    
    if (is_not_null(obj[["subclass"]])) {
        
        in.obj[["subclass"]] <- TRUE
    }
    else in.obj[["subclass"]] <- TRUE
    
    if (is_not_null(obj[["match.strata"]])) {
        
        in.obj[["match.strata"]] <- TRUE
    }
    else in.obj[["match.strata"]] <- TRUE
    
    if (is_not_null(obj[["cluster"]]) || is_not_null(obj[["exact"]])) {
        
        in.obj[["cluster"]] <- TRUE
    }
    else in.obj[["cluster"]] <- TRUE
    
    if (is_not_null(obj[["estimand"]]) || is_not_null(obj[["ATT"]]) || is_not_null(obj[["ATT"]])  || 
        is_not_null(obj[["target"]])) {
        
        in.obj[["estimand"]] <- TRUE
    }
    else in.obj[["estimand"]] <- TRUE
    
    if (is_not_null(obj[["cluster"]]) || is_not_null(obj[["exact"]])) {
        
        in.obj[["cluster"]] <- TRUE
    }
    else in.obj[["cluster"]] <- TRUE
    
    if (is_not_null(obj[["s.weights"]]) || is_not_null(obj[["sweights"]]) || is_not_null(obj[["sampw"]])) {
        
        in.obj[["s.weights"]] <- TRUE
    }
    else in.obj[["s.weights"]] <- TRUE
    
    if (is_not_null(obj[["focal"]])) {
        
        in.obj[["focal"]] <- TRUE
    }
    else in.obj[["focal"]] <- TRUE
    
    treat <- A$treat
    data <- A$data
    weights <- A$weights
    distance <- A$distance
    subclass <- A$subclass
    match.strata <- A$match.strata
    cluster <- A$cluster
    addl <- A$addl
    s.d.denom <- A$s.d.denom
    method <- A$method
    estimand <- A$estimand
    imp <- A$imp
    s.weights <- A$s.weights
    subset <- A$subset
    focal <- A$focal
    
    if (is_not_null(weights) && any(sapply(weights, function(x) any(is.na(x))))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (is_not_null(s.weights) && any(sapply(s.weights, function(x) any(is.na(x))))) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    #Checks
    if (is_null(covs)) {
        stop("covs dataframe must be specified.", call. = FALSE)
    }
    if (!is.data.frame(covs)) {
        stop("covs must be a data.frame.", call. = FALSE)
    }
    if (is_not_null(data) && !is.data.frame(data)) {
        warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
        data <- NULL
    }
    # if (length(distance) > 0 && !is.character(distance) && !is.numeric(distance) && !is.data.frame(distance)) {
    #     stop("The argument to distance must be a vector of distance scores or the (quoted) name of a variable in data that contains distance scores.", call. = FALSE)
    # }
    
    specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
    if (is_not_null(weights)) {
        if (!is.character(weights) && !is.numeric(weights) && !is.data.frame(weights) && !is.list(weights)) {
            stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
        }
        specified["weights"] <- TRUE
    }
    if (is_not_null(subclass)){
        if (!is.character(subclass) && !is.numeric(subclass)) {
            stop("The argument to subclass must be a vector of subclass membership or the (quoted) name of a variable in data that contains subclass membership.", call. = FALSE)
        }
        specified["subclass"] <- TRUE
    }
    if (is_not_null(match.strata)) {
        if (!is.character(match.strata) && !is.numeric(match.strata) && !is.factor(match.strata)) {
            stop("The argument to match.strata must be a vector of match stratum membership or the (quoted) name of a variable in data that contains match stratum membership.", call. = FALSE)
        }
        specified["match.strata"] <- TRUE
    }
    
    #Getting method
    if (is_null(method)) {
        if (specified["match.strata"]) {
            if (sum(specified) > 1) {
                message(word.list(names(specified)[specified]), " are specified. Assuming \"matching\" and using match.strata and ignoring ", word.list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                weights <- subclass <- NULL
            }
            X$method <- "matching"
        }
        else if (specified["subclass"]) {
            if (sum(specified) > 1) {
                message(word.list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word.list(names(specified)[specified & names(specified)!="subclass"]), ".")
                weights <- match.strata <- NULL
            }
            X$method <- "subclassification"
            #weights <- rep(1, nrow(covs))
        }
        else if (specified["weights"]) {
            if (sum(specified) > 1) {
                message(word.list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word.list(names(specified)[specified & names(specified)!="subclass"]), ".")
                match.strata <- subclass <- NULL
            }
            else {
                message("Assuming \"weighting\". If not, specify with an argument to method.")
            }
            X$method <- "weighting"
        }
        else {
            X$method <- "matching"
        }
    }
    else if (length(method) == 1) {
        specified.method <- match.arg(method, c("weighting", "matching", "subclassification"))
        if (specified.method == "weighting") {
            if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word.list(names(specified)[specified]), " are specified. Using weights and ignoring ", word.list(names(specified)[specified & names(specified)!="weights"]), ".")
                    match.strata <- subclass <- NULL
                }
                X$method <- "weighting"
            }
            else if (specified["match.strata"]) {
                message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
                subclass <- NULL
                X$method <- "matching"
            }
            else if (specified["subclass"]) {
                message("method = \"weighting\" is specified, but no weights are present. Assuming \"subclassification\" and using subclass instead.")
                X$method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else {
                X$method <- "matching"
            }
        }
        else if (specified.method == "matching") {
            if (specified["match.strata"]) {
                if (sum(specified) > 1) {
                    message(word.list(names(specified)[specified]), " are specified. Using match.strata and ignoring ", word.list(names(specified)[specified & names(specified)!="match.strata"]), ".")
                    weights <- subclass <- NULL
                }
                X$method <- "matching"
            }
            else if (specified["weights"]) {
                if (sum(specified) > 1) {
                    message(word.list(names(specified)[specified]), " are specified. Using weights and ignoring ", word.list(names(specified)[specified & names(specified)!="weights"]), ".")
                    match.strata <- subclass <- NULL
                }
                X$method <- "matching"
            }
            else if (specified["subclass"]) {
                message("method = \"matching\" is specified, but no weights or match.strata are present. Assuming \"subclassification\" and using subclass instead.")
                X$method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else {
                X$method <- "matching"
            }
        }
        else if (specified.method == "subclassification") {
            if (specified["subclass"]) {
                if (sum(specified) > 1) {
                    message(word.list(names(specified)[specified]), " are specified. Using subclass and ignoring ", word.list(names(specified)[specified & names(specified)!="subclass"]), ".")
                    weights <- match.strata <- NULL
                }
                X$method <- "subclassification"
                #weights <- rep(1, nrow(covs))
            }
            else if (specified["match.strata"]) {
                message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
                weights <- NULL
                X$method <- "matching"
            }
            else if (specified["weights"]) {
                message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"weighting\" and using weights instead.")
                X$method <- "weighting"
            }
        }
    }
    else {
        specified.method <- match.arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
        if (any(specified.method == "subclassification") || specified["subclass"]) {
            stop("Subclassification cannot be specified along with other methods.", call. = FALSE)
        }
        else if (specified["match.strata"]) {
            stop("Only weights can be specified with mutiple methods.", call. = FALSE)
        }
        else if (!specified["weights"]) {
            warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
            X$method <- "matching"
        }
        else {
            #Matching and/or weighting with various weights
            X$method <- specified.method
            match.strata <- subclass <- NULL
        }
    }
    
    if (is_not_null(cluster) && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
        stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (is_not_null(imp) && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
        stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
    }
    
    #Process treat
    if (is_null(treat)) stop("treat must be specified.", call. = FALSE)
    
    if (is.numeric(treat) || is.factor(treat) || (is.character(treat) && length(treat) > 1)) {
        treat <- treat
    }
    else if (is.character(treat) && length(treat)==1 && any(names(data) == treat)) {
        treat <- data[[treat]]
    }
    else stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
    
    if (sum(is.na(treat)) > 0)
        stop("Missing values exist in treat.", call. = FALSE)
    
    if (nunique(treat) == 2) {
        treat <- binarize(treat)
    }
    else if (is.character(treat)) {
        treat <- factor(treat)
    }
    
    #Process weights, addl, and distance
    for (i in c("weights", "addl", "distance")) {
        assign(i, data.frame.process(i, A[[i]], treat, covs, data))
    }
    
    #Process subclass
    if (is_not_null(subclass)) {
        if (is.numeric(subclass) || is.factor(subclass) || (is.character(subclass) && length(subclass) > 1)) {
            subclass <- factor(subclass)
        }
        else if (is.character(subclass) && length(subclass)==1 && any(names(data) == subclass)) {
            subclass <- factor(data[[subclass]])
        }
        else stop("The name supplied to subclass is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process match.strata
    if (is_not_null(match.strata)) {
        if (is.character(match.strata) && length(match.strata)==1) {
            if (any(names(data) == match.strata)) {
                match.strata <- data[[match.strata]]
            }
            else stop("The name supplied to match.strata is not the name of a variable in data.", call. = FALSE)
        }
    }
    
    #Process sampling weights
    if (is_not_null(s.weights)) {
        if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
            stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
        }
        if (is.character(s.weights) && length(s.weights)==1) {
            if (any(names(data) == s.weights)) {
                s.weights <- data[[s.weights]]
            }
            else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
        }
    }
    
    #Process cluster
    if (is_not_null(cluster)) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
            cluster <- data[[cluster]]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
        if (!is.logical(subset)) {
            stop("The argument to subset must be a logical vector.", call. = FALSE)
        }
        if (any(is.na(subset))) {
            warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
            subset[is.na(subset)] <- FALSE
        }
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    #Process imp
    if (is_not_null(imp)) {
        if (is.numeric(imp) || is.factor(imp) || (is.character(imp) && length(imp)>1)) {
            imp <- imp
        }
        else if (is.character(imp) && length(imp)==1 && any(names(data) == imp)) {
            imp <- data[[imp]]
        }
        else stop("The name supplied to imp is not the name of a variable in data.", call. = FALSE)
        
        imp.lengths <- sapply(unique(imp), function(i) sum(imp == i))
        
        if (all_the_same(imp.lengths)) { #all the same
            for (i in vectors) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) { 
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        
                        temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, temp.merge[[-c(1:3)]][order(temp.merge[[3]])])
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
            for (i in data.frames) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) {
                    if (lengths[i] == imp.lengths[1]) {
                        temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                               order2 = seq_along(imp))
                        temp.var <- data.frame(sort(imp),rep(seq_len(lengths[i]), length(imp.lengths)),
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, setNames(temp.merge[[-c(1:3)]][order(temp.merge[[3]])], names(get(i))))
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
        }
        else {
            problematic <- lengths > 0 & lengths != length(imp)
        }
        if (any(problematic)) {
            stop(paste0(word.list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
        }
        else ensure.equal.lengths <- FALSE
    }
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames[data.frames!="covs"])) {
            if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        stop(paste0(word.list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
    }
    
    #Turn match.strata into weights
    if (is_not_null(match.strata)) {
        weights <- data.frame(weights = match.strata2weights(match.strata = match.strata,
                                                             treat = treat,
                                                             covs = covs
        ))
    }
    
    if (is_not_null(weights)) {
        if (any(sapply(weights, function(x) !is.numeric(x)))) {
            stop("All weights must be numeric.", call. = FALSE)
        }
        if (length(X$method) == 1) {
            X$method <- rep(X$method, ncol(weights))
        }
        else if (length(X$method) != ncol(weights)) {
            stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
        }
        
    }
    
    #Check focal
    if (is_not_null(focal) && is.factor(treat)) {
        if (is.numeric(focal)) {
            if (focal <= nunique(treat)) focal <- levels(treat)[focal]
            else 
                stop(paste0("focal was specified as ", focal, 
                            ", but there are only ", levels(treat), " treatment groups."), call. = FALSE)
        }
        else {
            if (!any(levels(treat) == focal)) 
                stop(paste0("The name specified to focal is not the name of any treatment group."), call. = FALSE)
        }
    }
    
    #Get s.d.denom
    if (nunique(treat) <= 2 || !is.numeric(treat)) { #non-continuous
        check.estimand <- check.weights <- check.focal <- bad.s.d.denom <- bad.estimand <- FALSE
        if (is_not_null(s.d.denom)) {
            try.s.d.denom <- tryCatch(match.arg(s.d.denom, c("treated", "control", "pooled"), several.ok = TRUE),
                                      error = function(cond) FALSE)
            if (any(try.s.d.denom == FALSE)) {
                check.estimand <- TRUE
                bad.s.d.denom <- TRUE
            }
            else {
                if (length(try.s.d.denom) > 1 && length(try.s.d.denom) != ncol(weights)) {
                    stop("s.d.denom must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
                }
                else X$s.d.denom <- try.s.d.denom
            }
        }
        else {
            check.estimand <- TRUE
        }
        
        if (check.estimand == TRUE) {
            if (is_not_null(estimand)) {
                try.estimand <- tryCatch(match.arg(tolower(estimand), c("att", "atc", "ate"), several.ok = TRUE),
                                         error = function(cond) FALSE)
                if (any(try.estimand == FALSE)) {
                    check.focal <- TRUE
                    bad.estimand <- TRUE
                }
                else {
                    if (length(try.estimand) > 1 && length(try.estimand) != ncol(weights)) {
                        stop("estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
                    }
                    else X$s.d.denom <- sapply(try.estimand, function(x) switch(x, att = "treated", atc = "control", ate = "pooled"))
                }
            }
            else {
                check.focal <- TRUE
            }
        }
        if (check.focal == TRUE) {
            if (is_not_null(focal)) {
                X$s.d.denom <- "treated"
                estimand <- "att"
            }
            else check.weights <- TRUE
        }
        if (check.weights == TRUE) {
            if (is_null(weights)) {
                X$s.d.denom <- "pooled"
                estimand <- "ate"
            }
            else {
                X$s.d.denom <- estimand <- character(ncol(weights))
                for (i in seq_len(ncol(weights))) {
                    if (X$method[i] == "weighting") {
                        if (nunique(treat) <= 2) {
                            if (all_the_same(weights[[i]][treat==1 & !check_if_zero(weights[[i]])]) &&
                                !all_the_same(weights[[i]][treat==0 & !check_if_zero(weights[[i]])])
                            ) { #if treated weights are the same and control weights differ; ATT
                                estimand[i] <- "att"
                                X$s.d.denom[i] <- "treated"
                            }
                            else if (all_the_same(weights[[i]][treat==0 & !check_if_zero(weights[[i]])]) &&
                                     !all_the_same(weights[[i]][treat==1 & !check_if_zero(weights[[i]])])
                            ) { #if control weights are the same and treated weights differ; ATC
                                estimand[i] <- "atc"
                                X$s.d.denom[i] <- "control"
                            }
                            else {
                                estimand[i] <- "ate"
                                X$s.d.denom[i] <- "pooled"
                            }
                        }
                        else {
                            if (length(focal) == 1) {
                                estimand[i] <- "att"
                                X$s.d.denom[i] <- "treated"
                            }
                            else {
                                estimand[i] <- "ate"
                                X$s.d.denom[i] <- "pooled"
                            }
                        }
                    }
                    else {
                        estimand[i] <- "att"
                        X$s.d.denom[i] <- "treated"
                    }
                }
            }
        }
        if (is_not_null(weights) && length(X$s.d.denom) == 1) X$s.d.denom <- rep(X$s.d.denom, ncol(weights))
        
        if (bad.s.d.denom && bad.estimand) {
            message("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\n         Using \"", word.list(X$s.d.denom), "\" instead.")
        }
        else if (bad.estimand) {
            message("Warning: estimand should be one of \"ATT\", \"ATC\", or \"ATE\". Using \"", ifelse(all_the_same(estimand), toupper(estimand)[1], word.list(toupper(estimand))), "\" instead.")
        }
        else if (check.focal || check.weights) {
            message("Note: estimand and s.d.denom not specified; assuming ", ifelse(nunique(toupper(estimand)) == 1, toupper(unique(estimand)), word.list(toupper(estimand))), " and ", ifelse(nunique(X$s.d.denom) == 1, unique(X$s.d.denom), word.list(X$s.d.denom)), ".")
        }
        
        if (all(X$method %in% c("weighting", "matching"))) {
            if (is_not_null(weights) && length(X$s.d.denom) != ncol(weights)) {
                stop("Valid inputs to s.d.denom or estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
            }
        }
    }
    
    if (any(is.na(c(covs, addl)))) {
        warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    if (is_null(subset)) subset <- rep(TRUE, length(treat))
    
    X$covs <- covs[subset, , drop = FALSE]
    X$weights <- weights[subset, , drop = FALSE]
    X$treat <- treat[subset]
    X$distance <- distance[subset, , drop = FALSE]
    X$subclass <- subclass[subset]
    X$cluster <- factor(cluster[subset])
    X$call <- NULL
    X$addl <- addl[subset, , drop = FALSE]
    X$imp <- factor(imp[subset])
    X$s.weights <- s.weights[subset]
    
    return(X)
}

skew.diff <- function(x, group, weights, var.type) {
    #Calculates weighted skew and skew difference for groups; uses formulas at http://www.nematrian.com/WeightedMomentsAndCumulants
    skew.cols <- data.frame(Skew.C = NA, Skew.T = NA, Skew.Diff = NA)
    if (is.null(weights)) weights <- rep(1, length(x))
    if (var.type=="Contin.") {
        # skew for each group
        x0 <- x[group==0];       x1 <- x[group==1];
        w0 <- weights[group==0]; w1 <- weights[group==1]
        m0 <- w.m(x0, w0);       m1 <- w.m(x1, w1)
        wvp <- function(x, w, m) return(sum(w*(x-m)^2, na.rm=TRUE)/sum(w, na.rm=TRUE)) #weighted variance population
        Fv <- function(w) return((sum(w0, na.rm=TRUE)^2-sum(w0^2, na.rm=TRUE))/(sum(w0, na.rm=TRUE)^2)) #adjustment for sample variance
        wskp <- function(x, w, m, sp) return(sum(w*((x-m)/sp)^3, na.rm=TRUE)/sum(w, na.rm=TRUE)) #weighted skew population
        Fsk <- function(w, sp, s) return(((sum(w0, na.rm=TRUE)^3 - 3*sum(w0^2, na.rm=TRUE)*sum(w0, na.rm=TRUE) + 2*sum(w0^3, na.rm=TRUE))/(sum(w0, na.rm=TRUE)^3))*(sp/s)^3 ) #adjustment for sample skew
        #group==0
        sp0 <- sqrt(wvp(x0, w0, m0))
        s0 <- sqrt(wvp(x0, w0, m0)/Fv(w0))
        skp0 <- wskp(x0, w0, m0, sp0)
        sk0 <- skp0/Fsk(w0, sp0, s0)
        #group==1
        sp1 <- sqrt(wvp(x1, w1, m1))
        s1 <- sqrt(wvp(x1, w1, m1)/Fv(w1))
        skp1 <- wskp(x1, w1, m1, sp1)
        sk1 <- skp1/Fsk(w1, sp1, s1)
        
        skew.cols[, ] <- c(sk0, sk1, sk1-sk0)
    }
    return(skew.cols)
}

