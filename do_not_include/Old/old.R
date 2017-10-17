#Retired functions and versions of functions

x2base.data.frame0 <- function(covs, ...) {
    #treat
    #data
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
              obj = NA,
              cluster = NA,
              imp = NA)
    
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
    
    #Checks
    if (length(covs) == 0) {
        stop("covs dataframe must be specified.", call. = FALSE)
    }
    if (!is.data.frame(covs)) {
        stop("covs must be a dataframe.", call. = FALSE)
    }
    if (sum(is.na(covs)) > 0) {
        stop("Missing values exist in the covariates.", call. = FALSE)
    }
    if (length(data) > 0 && !is.data.frame(data)) {
        warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execuction will halt.")
        data <- NULL
    }
    if (length(distance) > 0 && !is.character(distance) && !is.numeric(distance) && !is.data.frame(distance)) {
        stop("The argument to distance must be a vector of distance scores or the (quoted) name of a variable in data that contains distance scores.", call. = FALSE)
    }
    
    specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
    if (length(weights) > 0) {
        if (!is.character(weights) && !is.numeric(weights)) {
            stop("The argument to weights must be a vector of weights or the (quoted) name of a variable in data that contains weights.", call. = FALSE)
        }
        specified["weights"] <- TRUE
    }
    if (length(subclass) > 0){
        if (!is.character(subclass) && !is.numeric(subclass)) {
            stop("The argument to subclass must be a vector of subclass membership or the (quoted) name of a variable in data that contains subclass membership.", call. = FALSE)
        }
        specified["subclass"] <- TRUE
    }
    if (length(match.strata) > 0) {
        if (!is.character(match.strata) && !is.numeric(match.strata) && !is.factor(match.strata)) {
            stop("The argument to match.strata must be a vector of match stratum membership or the (quoted) name of a variable in data that contains match stratum membership.", call. = FALSE)
        }
        specified["match.strata"] <- TRUE
    }
    
    #Getting method
    if (length(method) == 0) {
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
            weights <- rep(1, nrow(covs))
        }
        else if (specified["weights"]) {
            if (sum(specified) > 1) {
                message(word.list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word.list(names(specified)[specified & names(specified)!="subclass"]), ".")
                match.strata <- match.strata <- NULL
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
    else {
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
                weights <- rep(1, nrow(covs))
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
                weights <- rep(1, nrow(covs))
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
                weights <- rep(1, nrow(covs))
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
    
    if (length(cluster) > 0 && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
        stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (length(imp) > 0 && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
        stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
    }
    if (length(treat) == 0) stop("treat must be specified.", call. = FALSE)
    
    #Process treat
    if (is.numeric(treat) || is.factor(treat) || (is.character(treat) && length(treat) > 1)) {
        treat <- treat
    }
    else if (is.character(treat) && length(treat)==1 && treat %in% names(data)) {
        treat <- data[, treat]
    }
    else stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
    
    if (sum(is.na(treat)) > 0)
        stop("Missing values exist in treat.", call. = FALSE)
    
    #Process weights
    if (length(weights) > 0) {
        if (is.numeric(weights)) {
            weights <- weights
        }
        else if (is.character(weights) && length(weights)==1 && weights %in% names(data)) {
            weights <- data[, weights]
        }
        else stop("The name supplied to weights is not the name of a variable in data.", call. = FALSE)
        
        if (sum(is.na(weights)) > 0)
            stop("Missing values exist in weights.", call. = FALSE)
    }
    
    #Process subclass
    if (length(subclass) > 0) {
        if (is.numeric(subclass) || is.factor(subclass) || (is.character(subclass) && length(subclass)>1)) {
            subclass <- subclass
        }
        else if (is.character(subclass) && length(subclass)==1 && subclass %in% names(data)) {
            subclass <- data[, subclass]
        }
        else stop("The name supplied to subclass is not the name of a variable in data.", call. = FALSE)
        
    }
    
    #Process match.strata
    if (length(match.strata) > 0) {
        if (is.character(match.strata) && length(match.strata)==1) {
            if (match.strata %in% names(data)) {
                match.strata <- data[, match.strata]
            }
            else stop("The name supplied to match.strata is not the name of a variable in data.", call. = FALSE)
        }
        
        weights <- match.strata2weights(covs = covs, 
                                        treat = treat, 
                                        match.strata = match.strata)
    }
    
    #Process cluster
    if (length(cluster) > 0) {
        if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
            cluster <- cluster
        }
        else if (is.character(cluster) && length(cluster)==1 && cluster %in% names(data)) {
            cluster <- data[, cluster]
        }
        else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process addl and distance
    for (i in c("addl", "distance")) {
        val <- A[[i]]
        val.df <- NULL
        if (length(val) > 0) {
            if (is.numeric(val)) {
                val.df <- setNames(data.frame(val), i)
            }
            else if (is.character(val)) {
                if (length(data) > 0 && any(val %in% names(data))) {
                    val.df <- data[, val[val %in% names(data)], drop = FALSE]
                    val <- val[!val %in% names(data)]
                }
                if (length(val) > 0) {
                    warning(paste("The following variable(s) named in", i, "are not in data and will be ignored: ",
                                  paste(val)))
                }
            }
            else if (is.data.frame(val)) {
                val.df <- val
            }
            else stop(paste("The names supplied to", i, "are not the name of a variable in data."), call. = FALSE)
            
            if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
            }
        }
        assign(i, val.df)
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("treat", "weights", "subclass", "match.strata", "cluster")
    data.frames <- c("covs", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(sapply(vectors, 
                                 function(x) length(get(x))), 
                          sapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 })), c(vectors, data.frames))
    #Process imp
    if (length(imp) > 0) {
        if (is.numeric(imp) || is.factor(imp) || (is.character(imp) && length(imp)>1)) {
            imp <- imp
        }
        else if (is.character(imp) && length(imp)==1 && imp %in% names(data)) {
            imp <- data[, imp]
        }
        else stop("The name supplied to imp is not the name of a variable in data.", call. = FALSE)
        
        imp.lengths <- sapply(unique(imp), function(i) sum(imp == i))
        
        if (abs(max(imp.lengths) - min(imp.lengths)) < sqrt(.Machine$double.eps)) { #all the same
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
                        assign(i, temp.merge[order(temp.merge[,3]), -c(1:3), drop = TRUE])
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
                                               get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), ]
                        )
                        temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                            by.y = 1:2, sort = FALSE)
                        assign(i, setNames(temp.merge[order(temp.merge[,3]), -c(1:3)], names(get(i))))
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
    
    #Get s.d.denom
    check.estimand <- check.weights <- bad.s.d.denom <- bad.estimand <- FALSE
    if (length(unique(treat)) <= 2 || !is.numeric(treat)) { #non-continuous
        if (length(s.d.denom) > 0) {
            try.s.d.denom <- tryCatch(match.arg(s.d.denom, c("treated", "control", "pooled")),
                                      error = function(cond) FALSE)
            if (try.s.d.denom == FALSE) {
                check.estimand <- TRUE
                bad.s.d.denom <- TRUE
            }
            else {
                X$s.d.denom <- try.s.d.denom
            }
        }
        else {
            check.estimand <- TRUE
        }
        
        if (check.estimand == TRUE) {
            if (length(estimand) > 0) {
                try.estimand <- tryCatch(match.arg(tolower(estimand), c("att", "atc", "ate")),
                                         error = function(cond) FALSE)
                if (try.estimand == FALSE) {
                    check.weights <- TRUE
                    bad.estimand <- TRUE
                }
                else {
                    X$s.d.denom <- switch(try.estimand, att = "treated", atc = "control", ate = "pooled")
                }
            }
            else {
                check.weights <- TRUE
            }
        }
        
        if (check.weights == TRUE) {
            if (X$method == "weighting") {
                if (max(weights[treat==1 & weights > sqrt(.Machine$double.eps)]) - min(weights[treat==1 & weights > sqrt(.Machine$double.eps)]) < sqrt(.Machine$double.eps) &&
                    max(weights[treat==0 & weights > sqrt(.Machine$double.eps)]) - min(weights[treat==0 & weights > sqrt(.Machine$double.eps)]) >= sqrt(.Machine$double.eps)
                ) { #if treated weights are only all either 0 the same; ATT
                    estimand <- "att"
                    X$s.d.denom <- "treated"
                }
                else if (max(weights[treat==0 & weights > sqrt(.Machine$double.eps)]) - min(weights[treat==0 & weights > sqrt(.Machine$double.eps)]) < sqrt(.Machine$double.eps) &&
                         max(weights[treat==1 & weights > sqrt(.Machine$double.eps)]) - min(weights[treat==1 & weights > sqrt(.Machine$double.eps)]) >= sqrt(.Machine$double.eps)
                ) { #if control weights are only all either 0 the same; ATC
                    estimand <- "atc"
                    X$s.d.denom <- "control"
                }
                else {
                    estimand <- "ate"
                    X$s.d.denom <- "pooled"
                }
            }
            else {
                X$s.d.denom <- "treated"
            }
        }
        
        if (bad.s.d.denom && bad.estimand) {
            message("Warning: s.d.denom should be one of \"treated\", \"control\", or \"pooled\".\n         Using ", deparse(X$s.d.denom), " instead.")
        }
        else if (bad.estimand) {
            message("Warning: estimand should be one of \"ATT\", \"ATC\", or \"ATE\". Using ", deparse(toupper(estimand)), " instead.")
        }
        else if (check.weights && X$method == "weighting") {
            message("Note: estimand and s.d.denom not specified; assuming ", deparse(toupper(estimand)), " and ", deparse(X$s.d.denom), ".")
        }
    }
    
    X$covs <- covs
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$subclass <- factor(subclass)
    X$cluster <- factor(cluster)
    X$call <- NULL
    X$addl <- addl
    X$obj <- data.frame(treat=X$treat, weights=NA)
    X$imp <- factor(imp)
    if (length(weights) > 0) X$obj$weights <- X$weights
    if (length(subclass) > 0) X$obj$subclass <- X$subclass
    if (length(cluster) > 0) X$obj$cluster <- X$cluster
    return(X)
}
