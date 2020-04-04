#bal.tab
process_designmatch <- function(x) {
    dm.b.names <- c("obj_total", "obj_dist_mat", "t_id", 
                    "c_id", "group_id", "time")
    dm.n.names <- c("obj_total", "obj_dist_mat", "id_1", 
                    "id_2", "group_id", "time")
    if (length(x) >= min(length(dm.b.names), length(dm.n.names)) && 
        (all(dm.b.names %in% names(x)) || all(dm.n.names %in% names(x)))) {
        class(x) <- c("designmatch")
    }
    return(x)
}
process_time.list <- function(x) {
    if (is_(x, "list")) {
        if (all(vapply(x, is.formula, logical(1)))) {
            class(x) <- c("formula.list", "time.list", class(x))
        }
        else if (all(vapply(x, is.data.frame, logical(1)))) {
            class(x) <- c("data.frame.list", "time.list", class(x))
        }
    }
    return(x)
}

#x2base
process_treat <- function(treat, data = NULL) {
    
    if (missing(treat)) stop("treat must be specified.", call. = FALSE)
    
    if (inherits(treat, "unprocessed.treat")) {
        attrs <- attributes(treat)
        renamed_original <- setNames(names(treat_vals(treat)), treat_vals(treat))
        treat <- factor(renamed_original[as.character(treat)], levels = renamed_original)
        for (at in c("treat_names", "treat_vals", "treat.type", "names"))
            attr(treat, at) <- attrs[[at]]
    }
    else {
        treat <- vector.process(treat, name = "treat", 
                                which = "treatment statuses", 
                                data = data, missing.okay = FALSE)
        
        treat <- assign.treat.type(treat)
        treat.type <- get.treat.type(treat)
        
        if (treat.type == "binary") {
            treat <- factor(treat, levels = sort(unique(treat, nmax = 2)))
            original_values <- levels(treat)
            if (!anyNA(suppressWarnings(as.numeric(as.character(treat))))) {
                treat_names(treat) <- setNames(c("Control", "Treated"), c("control", "treated"))
            }
            else treat_names(treat) <- setNames(original_values, c("control", "treated"))
            
            treat_vals(treat) <- setNames(original_values, treat_names(treat))
        }
        else if (treat.type == "multinomial") {
            treat <- factor(treat)
            treat_names(treat) <- setNames(levels(treat), levels(treat))
            treat_vals(treat) <- setNames(levels(treat), treat_names(treat))
        }
        attr(treat, "treat.type") <- treat.type
    }
    class(treat) <- c("processed.treat", class(treat))
    return(treat)
}
unprocess_treat <- function(treat) {
    if (inherits(treat, "processed.treat")) {
        t.names <- names(treat)
        attrs <- attributes(treat)
        treat <- treat_vals(treat)[as.character(treat)]
        attributes(treat) <- attrs
        class(treat) <- c("unprocessed.treat", class(treat_vals(treat)))
    }
    return(treat)
}
process_treat.list <- function(treat.list, data = NULL) {
    
    if (is_null(treat.list)) stop("treat.list must be specified.", call. = FALSE)
    if (!is_(treat.list, "list")) {
        treat.list <- as.list(treat.list)
    }
    
    treat.list.names <- vapply(seq_along(treat.list), function(ti) {
        if (is.character(treat.list[[ti]]) && length(treat.list[[ti]])==1L && is_not_null(data)) {
            treat.list[[ti]]
        }
        else if (is_not_null(names(treat.list))) names(treat.list)[ti]
        else as.character(ti)
    }, character(1L))
    treat.list <- lapply(treat.list, process_treat, data = data)
    names(treat.list) <- treat.list.names
    
    return(treat.list)
}
`treat_names<-` <- function(treat, value) {
    `attr<-`(treat, "treat_names", value)
}
treat_names <- function(treat) {
    attr(treat, "treat_names")
}
`treat_vals<-` <- function(treat, value) {
    `attr<-`(treat, "treat_vals", value)
}
treat_vals <- function(treat) {
    attr(treat, "treat_vals")
}
`[.processed.treat` <- function(x, ...) {
    y <- NextMethod("[")
    treat_names(y) <- treat_names(x)
    treat_vals(y) <- treat_vals(x)
    attr(y, "treat.type") <- attr(x, "treat.type")
    class(y) <- class(x)
    y
}
`[.unprocessed.treat` <- `[.processed.treat`

initialize_X <- function() {
    X.names <- c("covs",
                 "treat",
                 "weights",
                 "distance",
                 "addl",
                 "s.d.denom",
                 "call",
                 "cluster",
                 "imp",
                 "s.weights",
                 "focal",
                 "discarded",
                 "method",
                 "subclass",
                 "stats",
                 "thresholds")
    X <- make_list(X.names)
    return(X)
}
initialize_X_msm <- function() {
    X.names <- c("covs.list",
                 "treat.list",
                 "weights",
                 "distance.list",
                 "addl.list",
                 "s.d.denom",
                 "call",
                 "cluster",
                 "imp",
                 "s.weights",
                 "focal",
                 "discarded",
                 "method",
                 "subclass",
                 "stats",
                 "thresholds")
    X <- make_list(X.names)
    return(X)
}
weight.check <- function(w) {
    wname <- deparse(substitute(w))
    if (!is.list(w)) w <- list(w)
    if (any(vapply(w, anyNA, logical(1L)))) stop(paste0("NAs are not allowed in the ", wname, "."), call. = FALSE)
    if (any(vapply(w, function(x) any(!is.numeric(x)), logical(1L)))) stop(paste0("All ", wname, " must be numeric."), call. = FALSE)
    if (any(vapply(w, function(x) any(!is.finite(x)), logical(1L)))) stop(paste0("Infinite ", wname, " are not allowed."), call. = FALSE)
    if (any(vapply(w, function(x) any(x < 0), logical(1L)))) warning(paste0("Negative ", wname, " found."), call. = FALSE)
}
strata2weights <- function(strata, treat) {
    #Process strata into weights (similar to weight.subclass from MatchIt)
    
    #Checks
    if (!is_(strata, c("factor", "atomic"))) {
        stop("strata must be an atomic vector or factor.", call. = FALSE)
    }
    #Process treat
    treat <- process_treat(treat)
    
    matched <- !is.na(strata)
    strata <- factor(strata)
    strata.matched <- factor(strata[matched])
    
    weights <- rep(0, length(treat))
    
    if (get.treat.type(treat) == "continuous") {
        
        stop("strata cannot be turned into weights for continuous treatments.", call. = FALSE)
    }
    else {
        sub.tab <- table(treat[matched], strata.matched)[,levels(strata.matched)]
        totals <- colSums(sub.tab)
        
        # estimand <- get.estimand(subclass = strata, treat = treat)
        s.d.denom <- get.s.d.denom(NULL, subclass = strata, treat = treat, quietly = TRUE)
        
        if (s.d.denom %in% treat_vals(treat)) {
            sub.weights <- setNames(sub.tab[s.d.denom, levels(strata.matched)] / sub.tab[treat_vals(treat) != s.d.denom, levels(strata.matched)], 
                                    levels(strata.matched))
            weights[matched & treat == s.d.denom] <- 1
            weights[matched & treat != s.d.denom] <- sub.weights[strata[matched & treat != s.d.denom]]
        }
        
        else {
            for (tn in treat_vals(treat)) {
                sub.weights <- setNames(totals[levels(strata.matched)]/sub.tab[tn, levels(strata.matched)], 
                                        levels(strata.matched))
                weights[matched & treat == tn] <- sub.weights[strata[matched & treat == tn]]
                
            }
        }
        
        if (all(check_if_zero(weights))) 
            stop("No units were stratified", call. = FALSE)
        else {
            for (tnn in names(treat_names(treat))) {
                if (all(check_if_zero(weights[treat == treat_vals(treat)[treat_names(treat)[tnn]]])))
                    stop(paste("No", tnn, "units were stratified."), call. = FALSE)
            }
        }
    }
    return(weights)
}
use.tc.fd <- function(formula = NULL, data = NULL, treat = NULL, covs = NULL, needs.treat = TRUE, needs.covs = TRUE) {
    if (is_not_null(formula) && class(formula) == "formula") {
        D <- NULL
        if (is_not_null(data)) D <- data
        if (is_not_null(covs) && is_(covs, c("data.frame", "matrix"))) {
          if (is_not_null(D)) D <- data.frame(D, covs) 
          else D <- as.data.frame(covs)
        }
        t.c <- get.covs.and.treat.from.formula(formula, D, treat = treat)
        t.c <- list(treat = t.c[["treat"]], covs = t.c[["reported.covs"]], treat.name = t.c[["treat.name"]])
        attr(t.c, "which") <- "fd"
    }
    else {
        if (is.matrix(covs)) covs <- as.data.frame(covs)
        else if (is.character(covs)) {
          if (is_not_null(data) && is_(data, c("data.frame", "matrix"))) {
            if (any(covs %nin% colnames(data))) {
              stop("All entries in covs must be names of variables in data.", call. = FALSE)
            }
            covs <- as.data.frame(data)[covs]
          }
          else {
            stop("If covs is a character vector, data must be specified as a data.frame.", call. = FALSE)
          }
        }
        else if (!is.data.frame(covs)) stop("covs must be a data.frame of covariates.", call. = FALSE)
        if (!is.atomic(treat)) stop("treat must be an atomic vector of treatment statuses.", call. = FALSE)
        t.c <- list(treat = treat, covs = covs)
        attr(t.c, "which") <- "tc"
    }
    
    if (needs.covs && is_null(t.c[["covs"]])) stop("No covariates were specified.", call. = FALSE)
    if (needs.treat && is_null(t.c[["treat"]])) stop("No treatment variable was specified.", call. = FALSE)
    
    return(t.c)
}
process.val <- function(val, i, treat, covs, ...) {
    if (is.numeric(val)) {
        val.df <- setNames(data.frame(val), i)
    }
    else if (is.character(val)) {
        data.sets <- list(...)
        data.sets <- data.sets[!vapply(data.sets, is_null, logical(1))]
        if ((is_not_null(data.sets) && length(val) > max(vapply(data.sets, ncol, numeric(1)))) || length(val) == NROW(covs) || length(val) == length(treat)){
            val.df <- setNames(data.frame(val), i)
        }
        else {
            if (is_not_null(data.sets)) {
                val <- unique(val)
                val.df <- setNames(as.data.frame(matrix(NA, ncol = length(val), nrow = max(vapply(data.sets, nrow, numeric(1))))),
                                   val)
                not.found <- setNames(rep(FALSE, length(val)), val)
                for (v in val) {
                    found <- FALSE
                    k <- 1
                    while (found == FALSE && k <= length(data.sets)) {
                        if (v %in% names(data.sets[[k]])) {
                            val.df[[v]] <- data.sets[[k]][[v]]
                            found <- TRUE
                        }
                        else k <- k + 1
                    }
                    if (!found) not.found[v] <- TRUE
                }
                if (any(not.found)) {
                    warning(paste("The following variable(s) named in", i, "are not in any available data sets and will be ignored: ",
                                  paste(val[not.found])), call. = FALSE)
                    val.df <- val.df[!not.found]
                }
            }
            else {
                val.df <- NULL
                warning(paste0("Names were provided to ", i, ", but no argument to data was provided. Ignoring ", i,"."), 
                        call. = FALSE)
            }
        }
    }
    else if (is.data.frame(val)) {
        val.df <- val
    }
    else stop(paste("The argument supplied to", i, "must be a vector, a data.frame, or the names of variables in an available data set."), call. = FALSE)
    
    return(val.df)
}
data.frame.process <- function(i, df, treat, covs, ...) {
    val <- df
    val.df <- NULL
    if (is_not_null(val)) {
        if (is_(val, "list")) {
            val.list <- lapply(val, function(x) process.val(x, i, treat, covs, ...))
            if (is_null(names(val.list)) || "" %in% names(val.list)) {
                stop(paste("All entries in", i, "must have names."), call. = FALSE)
            }
            val.list <- lapply(seq_along(val.list), function(x) {
                if (NCOL(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                return(val.list[[x]])})
            if (!all_the_same(vapply(val.list, nrow, numeric(1)))) {
                stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
            }
            
            val.df <- setNames(do.call("cbind", val.list),
                               c(sapply(val.list, names)))
        }
        else {
            val.df <- process.val(val, i, treat, covs, ...)
        }
        if (is_not_null(val.df)) { if (anyNA(val.df)) {
            stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
        }
    }
    return(val.df)
}
list.process <- function(i, List, ntimes, call.phrase, treat.list, covs.list, ...) {
    val.List <- List
    if (is_not_null(val.List)) {
        if (class(val.List)[1] != "list") {
            val.List <- list(val.List)
        }
        if (length(val.List) == 1) {
            val.List <- replicate(ntimes, val.List)
        }
        else if (length(val.List) == ntimes) {
            
        }
        else {
            stop(paste0("The argument to ", i, " must be a list of the same length as the number of time points in ",  call.phrase, "."), call. = FALSE)
        }
        for (ti in seq_along(val.List)) {
            val <- val.List[[ti]]
            val.df <- NULL
            if (is_not_null(val)) {
                if (is_(val, "list")) {
                    val.list <- lapply(val, function(x) process.val(x, strsplit(i, ".list", fixed = TRUE)[[1]], treat.list[[ti]], covs.list[[ti]], ...))
                    val.list <- lapply(seq_along(val.list), function(x) {
                        if (NCOL(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                        val.list[[x]]})
                    if (!all_the_same(vapply(val.list, nrow, numeric(1)))) {
                        stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                    }
                    
                    val.df <- setNames(do.call("cbind", val.list),
                                       c(vapply(val.list, names, character(1))))
                }
                else {
                    val.df <- process.val(val, strsplit(i, ".list", fixed = TRUE)[[1]], treat.list[[ti]], covs.list[[ti]], ...)
                }
                if (is_not_null(val.df)) { if (anyNA(val.df)) {
                    stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
                }
                val.List[[ti]] <- val.df
            }
            
        }
        val.df.lengths <- vapply(val.List[lengths(val.List) > 0], nrow, numeric(1))
        if (max(val.df.lengths) != min(val.df.lengths)) {
            stop(paste("All columns in", i, "need to have the same number of rows."), call. = FALSE)
        }
    }
    return(val.List)
}
vector.process <- function(vec, name = deparse(substitute(vec)), which = name, data = NULL, missing.okay = FALSE) {
    bad.vec <- FALSE
    if (is.character(vec) && length(vec)==1L && is_not_null(data)) {
        if (is_(data, "list")) {
            for (i in seq_along(data)) {
                if (is_(data[[i]], "matrix") && vec %in% colnames(data[[i]])) {
                    vec <- data[[i]][,vec]
                    break
                }
                else if (is_(data[[i]], "data.frame") && vec %in% names(data[[i]])) {
                    vec <- data[[i]][[vec]]
                    break
                }
                else if (i == length(data)) bad.vec <- TRUE
            }
        }
        else if (is_(data, "matrix") && vec %in% colnames(data)) {
            vec <- data[,vec]
        }
        else if (is_(data, "data.frame") && vec %in% names(data)) {
            vec <- data[[vec]]
        }
        else bad.vec <- TRUE
    }
    else if (is_(vec, c("atomic", "factor")) && length(vec) > 1L) {
        vec <- vec
    }
    else {
        bad.vec <- TRUE
    }
    
    if (bad.vec) stop(paste0("The argument to ", name, " must be a vector of ", which, " or the (quoted) name of a variable in data that contains ", which, "."), call. = FALSE)
    
    if (!missing.okay && anyNA(vec)) stop(paste0("Missing values exist in ", name, "."), call. = FALSE)
    
    return(vec)
} 
get.s.d.denom <- function(s.d.denom, estimand = NULL, weights = NULL, subclass = NULL, treat = NULL, focal = NULL, quietly = FALSE) {
    check.estimand <- check.weights <- check.focal <- bad.s.d.denom <- bad.estimand <- FALSE
    s.d.denom.specified <- is_not_null(s.d.denom)
    estimand.specified <- is_not_null(estimand)
    treat.is.processed <- is_(treat, "processed.treat")
    
    if (s.d.denom.specified) {
        if (!treat.is.processed) {
            treat <- process_treat(treat)
        }
        unique.treats <- as.character(treat_vals(treat))
        allowable.s.d.denoms <- c("pooled", "all", "weighted", "hedges")
        if (length(treat_names(treat)) == 2 && all(c("treated", "control") %in% names(treat_names(treat))))
            allowable.s.d.denoms <- c(allowable.s.d.denoms, "treated", "control")
        if (is_not_null(focal)) allowable.s.d.denoms <- c(allowable.s.d.denoms, "focal")
        
        if (length(s.d.denom) > 1 && length(s.d.denom) != NCOL(weights)) {
            stop("s.d.denom must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
        }
        
        s.d.denom <- match_arg(s.d.denom, unique(c(unique.treats, allowable.s.d.denoms)), 
                               several.ok = TRUE)
        
        if (any(s.d.denom %in% c("treated", "control")) && length(treat_names(treat) == 2)) {
            s.d.denom[s.d.denom %in% c("treated", "control")] <- treat_vals(treat)[treat_names(treat)[s.d.denom[s.d.denom %in% c("treated", "control")]]]
        }
        else if (any(s.d.denom %in% "focal")) check.focal <- TRUE
        
    }
    else {
        check.estimand <- TRUE
    }
    
    if (check.estimand) {
        if (estimand.specified) {
            if (!treat.is.processed) treat <- process_treat(treat)
            try.estimand <- tryCatch(match_arg(toupper(estimand), c("ATT", "ATC", "ATE", "ATO", "ATM"), several.ok = TRUE),
                                     error = function(cond) NA_character_)
            if (anyNA(try.estimand)) {
                check.focal <- TRUE
                # bad.estimand <- TRUE
            }
            else if (any(try.estimand %in% c("ATC", "ATT")) && length(treat_vals(treat)) > 2) {
                check.focal <- TRUE
            }
            else {
                if (length(try.estimand) > 1 && length(try.estimand) != NCOL(weights)) {
                    stop("estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
                }
                else s.d.denom <- vapply(try.estimand, function(x) switch(x, 
                                                                          ATT = treat_vals(treat)[treat_names(treat)["treated"]], 
                                                                          ATC = treat_vals(treat)[treat_names(treat)["control"]], 
                                                                          ATO = "weighted",
                                                                          ATM = "weighted",
                                                                          "pooled"), FUN.VALUE = character(1L))
            }
        }
        else {
            check.focal <- TRUE
        }
    }
    if (check.focal) {
        if (is_not_null(focal) && NCOL(weights) == 1) {
            s.d.denom <- focal
        }
        else check.weights <- TRUE
    }
    if (check.weights) {
        if (!treat.is.processed) treat <- process_treat(treat)
        if (is_null(weights) && is_null(subclass)) {
            s.d.denom <- "pooled"
        }
        else if (is_not_null(subclass)) {
            sub.tab <- table(treat, subclass)[treat_vals(treat), ]
            sub.tab <- rbind(sub.tab, table(subclass)[colnames(sub.tab)])
            dimnames(sub.tab) <- list(c(treat_vals(treat), "pooled"), colnames(sub.tab))
            
            ranges <- apply(sub.tab, 1, function(x) mean.abs.dev(x)/sum(x))
            s.d.denom <- rownames(sub.tab)[which.min(ranges)]
        }
        else {
            s.d.denom <- vapply(weights, function(w) {
                for (tv in treat_vals(treat)) {
                    if (all_the_same(w[treat == tv]) &&
                        !all_the_same(w[treat != tv])) {
                        return(tv)
                    }
                }
                return(NA_character_)
            }, character(1L))
            s.d.denom[is.na(s.d.denom)] <- "pooled"
        }
    }
    if (is_not_null(weights) && length(s.d.denom) == 1) s.d.denom <- rep.int(s.d.denom, NCOL(weights))
    
    if (s.d.denom.specified && bad.s.d.denom && (!estimand.specified || bad.estimand)) {
        attr(s.d.denom, "note") <- paste0("Warning: s.d.denom should be one of ", word_list(unique(c(unique.treats, allowable.s.d.denoms)), "or", quotes = TRUE), 
                                          ".\n         Using ", word_list(s.d.denom, quotes = TRUE), " instead.")
    }
    else if (estimand.specified && bad.estimand) {
        attr(s.d.denom, "note") <- paste0("Warning: estimand should be one of ", word_list(c("ATT", "ATC", "ATE"), "or", quotes = TRUE), 
                                          ". Ignoring estimand.")
    }
    else if ((check.focal || check.weights) && any(s.d.denom %nin% treat_vals(treat))) {
        attr(s.d.denom, "note") <- paste0("Note: s.d.denom not specified; assuming ", 
                                          if (all_the_same(s.d.denom)) s.d.denom[1] 
                                          else word_list(paste0("\"", vapply(s.d.denom, function(s) {
                                              if (s %in% treat_vals(treat) && all(treat_vals(treat) %in% c("0", "1"))) names(treat_names(treat))[treat_names(treat) == names(treat_vals(treat))[treat_vals(treat) == s]] else s
                                              }, character(1L)), "\" for ", names(weights))), ".")
    }
    
    if (is_not_null(weights) && length(s.d.denom) != NCOL(weights)) {
        stop("Valid inputs to s.d.denom or estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
    }
    
    if (!quietly && is_not_null(attr(s.d.denom, "note"))) message(attr(s.d.denom, "note"))
    
    return(s.d.denom)
}
get.s.d.denom.cont <- function(s.d.denom, weights = NULL, subclass = NULL, quietly = FALSE) {
    bad.s.d.denom <- FALSE
    s.d.denom.specified <- is_not_null(s.d.denom)

    if (is_not_null(subclass)) {
        s.d.denom <- "all"
    }
    else if (s.d.denom.specified) {
        allowable.s.d.denoms <- c("all", "weighted")

        if (length(s.d.denom) > 1 && length(s.d.denom) != NCOL(weights)) {
            stop("s.d.denom must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
        }
        
        s.d.denom <- match_arg(s.d.denom, unique(allowable.s.d.denoms), 
                               several.ok = TRUE)
    }
    else {
        s.d.denom <- "all"
    }

    if (is_not_null(weights) && length(s.d.denom) == 1) s.d.denom <- rep.int(s.d.denom, ncol(weights))
    
    if (!quietly) {
        if (s.d.denom.specified && bad.s.d.denom) {
            message(paste0("Warning: s.d.denom should be ", word_list(unique(allowable.s.d.denoms), "or", quotes = TRUE), 
                           ".\n         Using ", word_list(s.d.denom, quotes = TRUE), " instead."))
        }
    }
    
    if (is_not_null(weights) && length(s.d.denom) != NCOL(weights)) {
        stop("Valid inputs to s.d.denom or estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
    }
    
    if (is_not_null(weights)) names(s.d.denom) <- names(weights)
    
    return(s.d.denom)
}
get.estimand <- function(estimand = NULL, weights = NULL, subclass = NULL, treat = NULL, focal = NULL, quietly = TRUE) {
    check.weights <- check.focal <- FALSE
    
    if (is_not_null(estimand)) {
        try.estimand <- tryCatch(match_arg(toupper(estimand), c("ATT", "ATC", "ATE"), several.ok = TRUE),
                                 error = function(cond) NA_character_)
        if (anyNA(try.estimand)) {
            check.focal <- TRUE
            # bad.estimand <- TRUE
        }
        else {
            if (length(try.estimand) > 1 && length(try.estimand) != NCOL(weights)) {
                stop("estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
            }
            else estimand <- try.estimand
        }
    }
    else {
        check.focal <- TRUE
    }
    
    if (check.focal) {
        if (is_not_null(focal)) {
            estimand <- "ATT"
        }
        else check.weights <- TRUE
    }
    if (check.weights) {
        if (is_null(weights) && is_null(subclass)) {
            estimand <- "ATE"
        }
        else if (is_not_null(subclass)) {
            sub.tab <- table(treat, subclass)[treat_vals(treat), ]
            sub.tab <- rbind(sub.tab, table(subclass)[colnames(sub.tab)])
            dimnames(sub.tab) <- list(c(treat_vals(treat), "Total"), colnames(sub.tab))
            
            ranges <- apply(sub.tab, 1, function(x) mean.abs.dev(x)/sum(x))
            min.range <- which.min(ranges)
            if (rownames(sub.tab)[min.range] == treat_vals(treat)[treat_names(treat)["control"]]) estimand <- "ATC"
            else if (rownames(sub.tab)[min.range] == treat_vals(treat)[treat_names(treat)["treated"]]) estimand <- "ATT"
            else estimand <- "ATE"
            
        }
        else {
            estimand <- vapply(weights, function(w) {
                for (tnn in names(treat_names(treat))) {
                    if (all_the_same(w[treat == treat_vals(treat)[treat_names(treat)[tnn]]]) &&
                        !all_the_same(w[treat != treat_vals(treat)[treat_names(treat)[tnn]]])) {
                        return(switch(tnn, "control" = "ATC", "treated" = "ATT", NA_character_))
                    }
                }
                return(NA_character_)
            }, character(1L))
            estimand[is.na(estimand)] <- "ATE"
            
        }
    }
    if (is_not_null(weights) && length(estimand) == 1) estimand <- rep.int(estimand, ncol(weights))
    
    if (!quietly && (check.focal || check.weights)) {
        message("Note: estimand not specified; assuming ", ifelse(all_the_same(toupper(estimand)), toupper(estimand[1]), word_list(paste("\"", toupper(estimand), "\" for ", names(weights)))), ".")
    }
    
    if (is_not_null(weights) && length(estimand) != ncol(weights)) {
        stop("Valid inputs to estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
    }
    return(estimand)
}
assign.X.class <- function(X) {
    X <- clear_null(X)
    
    if (is_not_null(X[["treat"]]) && !has.treat.type(X[["treat"]])) X[["treat"]] <- assign.treat.type(X[["treat"]])
    
    if (is_not_null(X[["subclass"]])) {
        if (get.treat.type(X[["treat"]]) == "binary") X.class <- "subclass.binary"
        else if (get.treat.type(X[["treat"]]) == "continuous") X.class <- "subclass.cont"
        else stop("Multi-category treatments are not currently compatible with subclasses.", call. = FALSE)
    }
    else if (is_not_null(X[["cluster"]]) && nlevels(X[["cluster"]]) > 1) X.class <- "cluster"
    else if (is_not_null(X[["covs.list"]])) X.class <- "msm"
    else if (get.treat.type(X[["treat"]]) == "multinomial") X.class <- "multi"
    else if (is_not_null(X[["imp"]]) && nlevels(X[["imp"]]) > 1) X.class <- "imp"
    else if (get.treat.type(X[["treat"]]) == "binary") X.class <- "binary"
    else if (get.treat.type(X[["treat"]]) == "continuous") X.class <- "cont"
    else probably.a.bug()
    
    class(X) <- X.class
    
    return(X)
}
subsettable <- function() {
    c("covs",
      "treat",
      "weights",
      "distance",
      "addl",
      "cluster",
      "imp",
      "s.weights",
      "discarded",
      "subclass",
      "covs.list",
      "treat.list",
      "distance.list",
      "addl.list")
}
subset_X <- function(X, subset = NULL) {
    if (is_not_null(subset)) {
        if (is.logical(subset)) {
            if (!any(subset)) stop("All subset set to FALSE.", call. = FALSE)
            to_be_subset <- !all(subset)
            subset <- which(subset)
        }
        else if (is.numeric(subset)) to_be_subset <- TRUE
        else stop("subset must be logical or numeric.")
        
        if (to_be_subset) {
            subset_X_internal <- function(x, subset) {
                if (is_not_null(x)) {
                    if (is.factor(x)) factor(x[subset])
                    else if (is.atomic(x)) x[subset]
                    else if ((is.matrix(x) || is.data.frame(x))) x[subset, , drop = FALSE]
                    else if (is.list(x)) lapply(x, subset_X_internal, subset = subset)
                    else x
                }
                else x
            }
            X[names(X) %in% subsettable()] <- lapply(X[names(X) %in% subsettable()], subset_X_internal, subset)
        }
    }
    return(X)
}
imp.complete <- function(data) {
    if (!is_(data, "mids")) stop("'data' not of class 'mids'")
    
    single.complete <- function (data, where, imp, ell) {
        if (is.null(where)) where <- is.na(data)
        idx <- seq_len(ncol(data))[apply(where, 2, any)]
        for (j in idx) {
            if (is_null(imp[[j]])) data[where[, j], j] <- NA
            else data[where[, j], j] <- imp[[j]][, ell]
        }
        return(data)
    }
    
    m <- as.integer(data$m)
    idx <- 1L:m
    
    mylist <- lapply(idx, function(i) single.complete(data$data, data$where, data$imp, i))
    
    cmp <- data.frame(.imp = rep(idx, each = nrow(data$data)), 
                      .id = rep.int(1L:nrow(data$data), length(idx)), 
                      do.call(rbind, mylist))
    
    if (is.integer(attr(data$data, "row.names"))) 
        row.names(cmp) <- seq_len(nrow(cmp))
    else row.names(cmp) <- as.character(seq_len(nrow(cmp)))
    
    return(cmp)
}
length_imp_process <- function(vectors = NULL, data.frames = NULL, lists = NULL, imp = NULL, data = NULL, original.call.to = NULL, env = sys.frame(-1)) {
    #Processes imp and assigns it to imp in env
    #Processes vectors, data.frames, and lists for length
    #If object correspond to one imputation, assigns expanded object in env 
    
    ensure.equal.lengths <- TRUE
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
    lengths <- setNames(c(vapply(vectors, 
                                 function(x) {len(get0(x, envir = env, inherits = FALSE))
                                 }, numeric(1L)), 
                          vapply(data.frames, 
                                 function(x) {len(get0(x, envir = env, inherits = FALSE))
                                 }, numeric(1L)),
                          vapply(lists, function(x) {
                              if (is_null(get0(x, envir = env, inherits = FALSE))) 0 
                              else max(vapply(get(x, envir = env, inherits = FALSE), len, numeric(1L)))
                          }, numeric(1L))), 
                        c(vectors, data.frames, lists))
    
    #Process imp further
    if (is_not_null(imp)) {
        
        imp.lengths <- vapply(levels(imp), function(i) sum(imp == i), numeric(1L))
        
        if (all_the_same(imp.lengths)) { #all the same
            unsorted.imp <- is.unsorted(imp)
            for (i in vectors) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) { 
                    if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                    if (lengths[i] == imp.lengths[1]) {
                        i_obj <- get(i, envir = env, inherits = FALSE)
                        new_i <- i_obj[rep(seq_along(i_obj), length(imp.lengths))]
                        if (unsorted.imp) {for (i_ in levels(imp)) new_i[imp == i_] <- i_obj}
                        assign(i, new_i, pos = env)
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
            for (i in data.frames) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) {
                    if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                    if (lengths[i] == imp.lengths[1]) {
                        i_obj <- get(i, envir = env, inherits = FALSE)
                        new_i <- i_obj[rep(seq_len(nrow(i_obj)), length(imp.lengths)), , drop = FALSE]
                        if (unsorted.imp) {for (i_ in levels(imp)) new_i[imp == i_,] <- i_obj}
                        assign(i, new_i, pos = env)
                    }
                    else {
                        problematic[i] <- TRUE
                    }
                }
            }
            for (i in lists) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) {
                    if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
                    if (lengths[i] == imp.lengths[1]) {
                        assign(i, lapply(get(i, envir = env, inherits = FALSE), function(j) {
                            if (is_(j, c("atomic", "factor"))) {
                                newj <- j[rep(seq_along(j), length(imp.lengths))]
                                if (unsorted.imp) {for (i_ in levels(imp)) newj[imp == i_] <- j}
                                return(newj)
                            }
                            else if (is_(j, c("data.frame", "matrix"))) {
                                newj <- j[rep(seq_len(nrow(j)), length(imp.lengths)), , drop = FALSE]
                                if (unsorted.imp) {for (i_ in levels(imp)) newj[imp == i_,] <- j}
                                return(newj)
                            }
                            else {
                                stop(paste(i, "can only contain vectors or data.frames."), call. = FALSE)
                            }
                        }), pos = env)
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
            stop(paste0(word_list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
        }
        else ensure.equal.lengths <- FALSE
        
        assign("imp", imp, pos = env)
    }
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
        for (i in c(vectors, data.frames, lists)) {
            if (lengths[i] > 0 && lengths[i] != lengths[c(lists, data.frames, vectors)[1]]) {
                problematic[i] <- TRUE
            }
        }
    }
    if (any(problematic)) {
        if (is_null(original.call.to)) anchor <- c(lists, data.frames, vectors)[1]
        else anchor <- paste("in the original call to", original.call.to)
        
        stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as ", anchor, "."), call. = FALSE)
    }
}
process_stats <- function(stats = NULL, treat) {
    if (is_(treat, "list")) {
        return(unique(unlist(lapply(treat, function(x) process_stats(stats, x)))))
    }
    else {
        if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
        treat.type <- get.treat.type(treat)
        
        if (treat.type %in% c("binary", "multinomial")) {
            if (is_null(stats)) stats <- getOption("cobalt_stats", "mean.diffs")
            stats <- match_arg(stats, all_STATS("bin"), several.ok = TRUE)
        }
        else if (treat.type %in% c("continuous")) {
            if (is_null(stats)) stats <- getOption("cobalt_stats", "correlations")
            stats <- match_arg(stats, all_STATS("cont"), several.ok = TRUE)
        }
        return(unique(stats))
    }
}
process_thresholds <- function(thresholds, stats) {
    if (is_not_null(thresholds)) {
        if (is.list(thresholds)) thresholds <- unlist(thresholds)
        if (!all(is.na(thresholds)) && !is.numeric(thresholds)) {
            stop("thresholds must be numeric.", call. = FALSE)
        }
        
        if (is_not_null(names(thresholds))) {
            names(thresholds) <- stats[pmatch(names(thresholds), stats, duplicates.ok = TRUE)]
            thresholds <- thresholds[!is.na(names(thresholds))]
        }
        else {
            names(thresholds) <- stats[1:min(length(stats), length(thresholds))]
        }
        
        thresholds[names(thresholds)] <- as.numeric(thresholds)
    }
    thresholds <- as.list(thresholds[!is.na(thresholds)])
    return(thresholds)
}
process_subset <- function(subset, n) {
    if (!is_(subset, c("logical", "numeric"))) {
        stop("The argument to subset must be a logical or numeric vector.", call. = FALSE)
    }
    else if (is.numeric(subset)) {
        if (any(abs(subset) > n)) stop("Numeric values for subset cannot be larger than the number of units.", call. = FALSE)
        subset <- subset[!is.na(subset) & subset != 0]
        if (any(subset < 0) && any(subset > 0)) stop("Positive and negative indices cannot be mixed with subset.")
        if (any(abs(subset) > n)) stop("If subset is numeric, none of its values can exceed the number of units.")
        logical.subset <- rep(any(subset < 0), n)
        logical.subset[abs(subset)] <- !logical.subset[abs(subset)]
        subset <- logical.subset
    }
    if (anyNA(subset)) {
        warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
        subset[is.na(subset)] <- FALSE
    }
    return(subset)
}

#get.C
#Functions to turn input covariates into usable form
#int.poly.f creates interactions and polynomials
#splitfactor splits factor variable into indicators (now in utilities)
#binarize transforms 2-value variable into binary (0,1)
#get.C controls flow and handles redunancy
#get.types gets variables types (contin./binary)

int.poly.f <- function(mat, ex=NULL, int=FALSE, poly=1, center = FALSE, sep, co.names) {
    #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
    #Only to be used in base.bal.tab; for general use see int.poly()
    #mat=matrix input
    #ex=matrix of variables to exclude in interactions and polynomials; a subset of df
    #int=whether to include interactions or not; currently only 2-way are supported
    #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
    #nunder=number of underscores between variables
    
    if (is_not_null(ex)) d <- mat[, colnames(mat) %nin% colnames(ex), drop = FALSE]
    else d <- mat
    binary.vars <- apply(d, 2, is_binary)
    if (center) {
        d[,!binary.vars] <- center(d[, !binary.vars, drop = FALSE])
    }
    nd <- NCOL(d)
    nrd <- NROW(d)
    no.poly <- binary.vars
    npol <- nd - sum(no.poly)
    new <- matrix(0, ncol = (poly-1)*npol + int*(.5*(nd)*(nd-1)), nrow = nrd)
    nc <- NCOL(new)
    new.co.names <- make_list(nc)
    if (poly > 1 && npol != 0) {
        for (i in 2:poly) {
            new[, (1 + npol*(i - 2)):(npol*(i - 1))] <- apply(d[, !no.poly, drop = FALSE], 2, function(x) x^i)
            new.co.names[(1 + npol*(i - 2)):(npol*(i - 1))] <- lapply(colnames(d)[!no.poly], function(x) setNames(list(c(co.names[[x]][["component"]], num_to_superscript(i)), c(co.names[[x]][["type"]], "power")), c("component", "type")))
            
        }
    }
    if (int && nd > 1) {
        new[,(nc - .5*nd*(nd-1) + 1):nc] <- matrix(t(apply(d, 1, combn, 2, prod)), nrow = nrd)
        new.co.names[(nc - .5*nd*(nd-1) + 1):nc] <- lapply(as.data.frame(combn(colnames(d), 2), stringsAsFactors = FALSE), 
                                                           function(x) setNames(list(c(co.names[[x[1]]][["component"]], sep, co.names[[x[2]]][["component"]]),
                                                                                     c(co.names[[x[1]]][["type"]], "isep", co.names[[x[2]]][["type"]])),
                                                                                c("component", "type")))
    }
    
    colnames(new) <- vapply(new.co.names, function(x) paste0(x[["component"]], collapse = ""), character(1))
    names(new.co.names) <- colnames(new)
    new <- new[, !apply(new, 2, all_the_same), drop = FALSE]
    attr(new, "co.names") <- new.co.names
    return(new)
}
get.C <- function(covs, int = FALSE, poly = 1, addl = NULL, distance = NULL, cluster = NULL, ...) {
    #gets C data.frame, which contains all variables for which balance is to be assessed. Used in balance.table.
    A <- list(...)
    if (is_null(A[["int_sep"]])) A[["int_sep"]] <- getOption("cobalt_int_sep", default = " * ")
    if (is_null(A[["factor_sep"]])) A[["factor_sep"]] <- getOption("cobalt_factor_sep", default = "_")
    if (is_null(A[["center"]]) || A[["center"]] %nin% c(TRUE, FALSE)) A[["center"]] <- getOption("cobalt_center", default = FALSE)
    
    C <- covs
    if (!is.null(addl)) {
        if (!is.data.frame(addl)) {
            if (is.character(addl)) stop("The argument to addl must be a data.frame containing the the values of the additional variables you want to include in the balance assessment.", call. = FALSE)
            else stop("The argument to addl must be a data.frame. Wrap data.frame() around the argument if it is a matrix or vector.", call. = FALSE)
        }
        else {
            repeat.name.indices <- vapply(names(addl), function(x) x %in% names(C), logical(1))
            if (any(repeat.name.indices)) {
                warning(paste("The following variables in addl have the same name as covariates and will be ignored:\n",
                              paste(names(addl)[repeat.name.indices], collapse = " ")), call. = FALSE)
                addl <- addl[!repeat.name.indices]
            }
            C <- cbind(C, addl)
        }
    } 
    
    covs.with.inf <- vapply(C, function(x) !is.character(x) && any(!is.finite(x) & !is.na(x)), logical(1L))
    if (any(covs.with.inf)) {
        s <- if (sum(covs.with.inf) == 1) c("", "s") else c("s", "")
        stop(paste0("The variable", s[1], " ", word_list(names(C)[covs.with.inf], quotes = TRUE), 
                    " contain", s[2], " non-finite values, which are not allowed."), call. = FALSE)
    }
    
    vars.w.missing <- data.frame(placed.after = names(C),
                                 has.missing = FALSE, 
                                 has.Inf = FALSE,
                                 row.names = names(C),
                                 stringsAsFactors = FALSE)
    co.names <- setNames(lapply(names(C), function(x) setNames(list(x, "base"), c("component", "type"))), names(C))
    #component types: base, fsep, isep, power, na, level
    is.0.1.cov <- setNames(rep(FALSE, ncol(C)), names(C))
    for (i in names(C)) {
        if (nunique(C[[i]]) == 2) {
            #if (is.logical(C[[i]])) C[[i]] <- as.numeric(C[[i]])
            if (((is.numeric(C[[i]]) || is.logical(C[[i]])) && 
                 all(check_if_zero(as.numeric(C[[i]]) - binarize(C[[i]])), na.rm = TRUE)) ||
                all(as.character(C[[i]]) %in% c("0", "1"))) {
                is.0.1.cov[i] <- TRUE
            }
            
            C[[i]] <- factor(C[[i]])
            C[[i]] <- relevel(C[[i]], levels(C[[i]])[2])
        }
        else if (is.character(C[[i]]) || is.factor(C[[i]])) C[[i]] <- factor(C[[i]])
        if (is_not_null(cluster) && !nunique.gt(C[[i]], nunique(cluster)) && 
            equivalent.factors(C[[i]], cluster)) {
            C <- C[names(C) != i] #Remove variable if it is the same (linear combo) as cluster variable
        }
        else {
            if (anyNA(C[[i]])) vars.w.missing[i, "has.missing"] <- TRUE
            if (!is.numeric(C[[i]])) {
                old.C.names <- names(C)
                C <- splitfactor(C, i, replace = TRUE, sep = A[["factor_sep"]], drop.first = FALSE, 
                                 drop.singleton = FALSE)
                newly.added.names <- names(C)[names(C) %nin% old.C.names]
                vars.w.missing[i, "placed.after"] <- last(newly.added.names)
                co.names <- c(co.names, setNames(lapply(newly.added.names, function(x) {
                    split.points <- c(nchar(i), nchar(i) + nchar(A[["factor_sep"]]))
                    split.names <- substring(x,
                                             c(1, split.points[1] + 1, split.points[2] + 1),
                                             c(split.points[1], split.points[2], nchar(x))
                    )
                    setNames(list(split.names, c("base", "fsep", "level")), 
                             c("component", "type"))
                }), newly.added.names))
            }
        }
    }
    
    if (NCOL(C) == 0) stop("There are no variables for which to display balance.", call. = FALSE)
    
    #Make sure categorical variable have missingness indicators done correctly
    
    C <- C[!vapply(C, all_the_same, logical(1L))]
    C <- as.matrix(C)
    
    #Process int and poly
    if (length(int) != 1L || !is.finite(int) || !(is.logical(int) || is.numeric(int))) {
        stop("int must be TRUE, FALSE, or a numeric value of length 1.", call. = FALSE)
    }
    if (int < 0 || !check_if_zero(abs(int - round(int)))) {
        stop("int must be TRUE, FALSE, or a numeric (integer) value greater than 1.", call. = FALSE)
    }
    int <- as.integer(round(int))
    
    if (length(poly) != 1L || !is.finite(poly) || !is.numeric(poly)) {
        stop("poly must be a numeric value of length 1.", call. = FALSE)
    }
    if (poly < 0 || !check_if_zero(abs(poly - round(poly)))) {
        stop("poly must be a numeric (integer) value greater than 1.", call. = FALSE)
    }
    poly <- as.integer(round(poly))
    
    if (int || poly) {
        if (int) { 
            #Prevent duplicate var names with `sep`s
            nsep <- 1
            repeat {
                all.possible.names <- outer(colnames(C), colnames(C), paste, sep = paste0(rep.int(A[["int_sep"]], nsep), collapse = ""))
                if (!any(colnames(C) %in% all.possible.names)) break
                else nsep <- nsep + 1
            }
            
            if (poly < int) poly <- int
            
            int <- TRUE
        }
        
        new <- int.poly.f(C, int = int, poly = poly, center = A[["center"]], 
                          sep = rep.int(A[["int_sep"]], nsep), co.names = co.names)
        C <- cbind(C, new)
        co.names <- c(co.names, attr(new, "co.names"))
    }
    
    #Add missingness indicators
    vars.w.missing <- vars.w.missing[vars.w.missing$placed.after %in% colnames(C) & vars.w.missing$has.missing, , drop = FALSE]
    if (NROW(vars.w.missing) > 0) {
        missing.ind <- apply(C[,colnames(C) %in% vars.w.missing$placed.after, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
        colnames(missing.ind) <- rownames(vars.w.missing)
        vars.w.missing <- vars.w.missing[colnames(missing.ind), , drop = FALSE]
        colnames(missing.ind) <- paste0(colnames(missing.ind), ":<NA>")
        original.var.order <- setNames(seq_len(NCOL(C)), colnames(C))
        new.var.order <- original.var.order + cumsum(c(0,(colnames(C) %in% vars.w.missing$placed.after)[-NCOL(C)]))
        new.C <- matrix(NA, nrow = NROW(C), ncol = NCOL(C) + NCOL(missing.ind),
                        dimnames = list(rownames(C), seq_len(NCOL(C) + NCOL(missing.ind))))
        new.C[, new.var.order] <- C
        new.C[, -new.var.order] <- missing.ind
        colnames(new.C)[new.var.order] <- colnames(C)
        colnames(new.C)[-new.var.order] <- colnames(missing.ind)
        miss.co.names <- setNames(lapply(rownames(vars.w.missing), function(x) setNames(list(c(x, ":<NA>"),
                                                                                             c("base", "na")), c("component", "type"))),
                                  colnames(missing.ind))
        C <- new.C
        missing.ind.C <- colnames(missing.ind)
        if (int) {
            new <- int.poly.f(missing.ind, int = TRUE, poly = 1, sep = rep.int(A[["int_sep"]], nsep), co.names = miss.co.names)
            C <- cbind(C, new)
            missing.ind.C <- c(missing.ind.C, colnames(new))
            miss.co.names <- c(miss.co.names, attr(new, "co.names"))
        }
        co.names <- c(co.names, miss.co.names)
        
    }

    #Remove duplicate & redundant variables
    C <- remove.perfect.col(C) 
        
    if (is_not_null(distance)) {
        if (any(names(distance) %in% colnames(C))) stop("distance variable(s) share the same name as a covariate. Please ensure each variable name is unique.", call. = FALSE)
        if (any(apply(distance, 2, function(x) anyNA(x)))) stop("Missing values are not allowed in the distance measure.", call. = FALSE)
        C <- cbind(as.matrix(distance), C, row.names = NULL)
        dist.co.names <- setNames(lapply(names(distance), function(x) setNames(list(x, "base"), c("component", "type"))), names(distance))
        co.names <- c(co.names, dist.co.names)
    }
    
    co.names <- co.names[colnames(C)]
    
    #Get rid of _1 for binary covs
    for (i in colnames(C)) {
        in.is.0.1.cov <- vapply(names(is.0.1.cov)[is.0.1.cov], 
                                function(i1) length(co.names[[i]][["component"]]) == 3 && 
                                    co.names[[i]][["component"]][1] == i1 && 
                                    co.names[[i]][["component"]][2] == A[["factor_sep"]] &&
                                    co.names[[i]][["component"]][3] %in% c("1", "TRUE"), 
                                logical(1L))
        
        if (any(in.is.0.1.cov)) {
            name.index <- which(colnames(C) == i)
            new.name <- names(is.0.1.cov)[is.0.1.cov][in.is.0.1.cov][1]
            colnames(C)[name.index] <- new.name
            names(co.names)[name.index] <- new.name
            co.names[[name.index]][["component"]] <- new.name
            co.names[[name.index]][["type"]] <- "base"
        }
    }
    
    attr(co.names, "seps") <- c(factor = A[["factor_sep"]], int = A[["int_sep"]])
    attr(C, "co.names") <- co.names
    if (is_not_null(distance)) attr(C, "distance.names") <- names(distance)
    if (NROW(vars.w.missing) > 0) attr(C, "missing.ind") <- missing.ind.C
    
    return(C)
    
}
get.types <- function(C) {
    vapply(colnames(C), function(x) {
        if (any(attr(C, "distance.names") == x)) "Distance"
        else if (is_binary(C[,x]))  "Binary"
        else "Contin."
    }, character(1))
}
remove.perfect.col <- function(C) {
    C.no.miss <- C[,colnames(C) %nin% attr(C, "missing.ind"), drop = FALSE]
    #If many rows, select subset to test redundancy
    if (NROW(C.no.miss) > 1500) {
        repeat {
            mini.C.no.miss <- C.no.miss[sample(seq_len(NROW(C.no.miss)), 1000),,drop=FALSE]
            single.value <- apply(mini.C.no.miss, 2, all_the_same)
            if (all(!single.value)) break
        }
        suppressWarnings(C.cor <- cor(mini.C.no.miss, use = "pairwise.complete.obs"))
    }
    else suppressWarnings(C.cor <- cor(C.no.miss, use = "pairwise.complete.obs"))
    
    s <- !lower.tri(C.cor, diag=TRUE) & !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
    redundant.vars <- colnames(C.no.miss)[apply(s, 2, any)]
    C <- C[, colnames(C) %nin% redundant.vars, drop = FALSE] 
    return(C)
}
# remove.perfect.col2 <- function(C) {
#   C.no.miss <- C[,colnames(C) %nin% attr(C, "missing.ind"), drop = FALSE]
#   cor_ <- function(x, y, mean_x, sd_x, mean_y, sd_y, n) {
#     cov_ <- sum((x - mean_x) * (y-mean_y))/(n-1)
#     return(cov_/(sd_x*sd_y))
#   }
#   n <- nrow(C.no.miss)
#   redundant.vars <- rep(FALSE, ncol(C.no.miss))
#   means <- colMeans(C.no.miss)
#   sds <- sqrt(col.w.v(C.no.miss))
#   nas <- apply(C.no.miss, 2, anyNA)
#   for (i in ncol(C.no.miss):2) {
#     j <- i-1
#     mean_i <- mean(C.no.miss[,i])
#     sd_i <- sd(C.no.miss[,i])
#     while(!redundant.vars[i] && j > 0) {
#       if (nas[i] && nas[j] && sum(is.na(C.no.miss[,i])) != sum(is.na(C.no.miss[,j]))) {} 
#       else if (check_if_zero(1-abs(suppressWarnings(cor_(C.no.miss[,j], C.no.miss[,i], means[j], sds[j], means[i], sds[i], n))))) {
#         redundant.vars[i] <- TRUE
#       }
#       j <- j - 1
#     }
#   }
#   C <- C[, colnames(C) %nin% colnames(C.no.miss)[redundant.vars], drop = FALSE] 
#   return(C)
# }
# remove.perfect.col3 <- function(C) {
#   C.no.miss <- C[,colnames(C) %nin% attr(C, "missing.ind"), drop = FALSE]
#   #If many rows, select subset to test redundancy
#   if (NROW(C.no.miss) > 1500) {
#     repeat {
#       mini.C.no.miss <- C.no.miss[sample(seq_len(NROW(C.no.miss)), 1000),,drop=FALSE]
#       single.value <- apply(mini.C.no.miss, 2, all_the_same)
#       if (all(!single.value)) break
#     }
#     suppressWarnings(C.cor <- coop::pcor(mini.C.no.miss, use = "pairwise.complete.obs"))
#   }
#   else suppressWarnings(C.cor <- coop::pcor(C.no.miss, use = "pairwise.complete.obs"))
#   
#   s <- !lower.tri(C.cor, diag=TRUE) & !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
#   redundant.vars <- colnames(C.no.miss)[apply(s, 2, any)]
#   C <- C[, colnames(C) %nin% redundant.vars, drop = FALSE] 
#   return(C)
# }
# remove.perfect.col4 <- function(C) {
#   C.no.miss <- C[,colnames(C) %nin% attr(C, "missing.ind"), drop = FALSE]
#   #If many rows, select subset to test redundancy
#   if (NROW(C.no.miss) > 1500) {
#     repeat {
#       mini.C.no.miss <- C.no.miss[sample(seq_len(NROW(C.no.miss)), 1000),,drop=FALSE]
#       single.value <- apply(mini.C.no.miss, 2, all_the_same)
#       if (all(!single.value)) break
#     }
#     suppressWarnings(C.cor <- Rfast::cora(mini.C.no.miss))
#   }
#   else suppressWarnings(C.cor <- Rfast::cora(C.no.miss))
#   
#   s <- !lower.tri(C.cor, diag=TRUE) & !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
#   redundant.vars <- colnames(C.no.miss)[apply(s, 2, any)]
#   C <- C[, colnames(C) %nin% redundant.vars, drop = FALSE] 
#   return(C)
# }

#base.bal.tab
check_if_zero_weights <- function(weights.df, treat = NULL) {
    #Checks if all weights are zero in each treat group for each set of weights
    if (is_not_null(treat)) {
        w.t.mat <- expand.grid(colnames(weights.df), treat_vals(treat), stringsAsFactors = FALSE)
        colnames(w.t.mat) <- c("weight_names", "treat_vals")
        
        if (NROW(w.t.mat) > 0) {
            problems <- vapply(seq_len(NROW(w.t.mat)), function(x) all(check_if_zero(weights.df[treat == w.t.mat[x, "treat_vals"], w.t.mat[x, "weight_names"]])), logical(1L))
            if (any(problems)) {
                prob.w.t.mat <- droplevels(w.t.mat[problems,])
                if (NCOL(weights.df) == 1) {
                    error <- paste0("All weights are zero when treat is ", word_list(prob.w.t.mat[, "treat_vals"], "or"), ".")
                }
                else {
                    errors <- setNames(character(nlevels(prob.w.t.mat[,"weight_names"])), levels(prob.w.t.mat[,"weight_names"]))
                    
                    for (i in levels(prob.w.t.mat[,"weight_names"])) {
                        errors[i] <- paste0("\"", i, "\" weights are zero when treat is ", word_list(prob.w.t.mat[prob.w.t.mat[,"weight_names"] == i, "treat_vals"], "or"))
                    }
                    errors <- paste(c("All", rep("all", length(errors)-1)), errors)
                    error <- paste0(word_list(errors, "and"), ".")
                }
                stop(error, call. = FALSE)
            }
        }
    }
    else {
        if (length(colnames(weights.df)) > 0) {
            problems <- vapply(colnames(weights.df), function(wn) all(check_if_zero(weights.df[, wn])), logical(1L))
            if (any(problems)) {
                prob.wts <- colnames(weights.df)[problems]
                if (NCOL(weights.df) == 1) {
                    error <- paste0("All weights are zero.")
                }
                else {
                    errors <- setNames(character(length(prob.wts)), prob.wts)
                    
                    for (i in prob.wts) {
                        errors[i] <- paste0("\"", i, "\" weights are zero")
                    }
                    errors <- paste(c("All", rep("all", length(errors)-1)), errors)
                    error <- paste0(word_list(errors, "and"), ".")
                }
                stop(error, call. = FALSE)
            }
        }
    }
}
compute_s.d.denom <- function(mat, treat, s.d.denom = "pooled", s.weights = NULL, bin.vars = NULL, subset = NULL, weighted.weights = NULL, to.sd = rep(TRUE, ncol(mat)), na.rm = TRUE) {
    denoms <- rep(1, ncol(mat))
    if (is.character(s.d.denom) && length(s.d.denom) == 1L) {
        
        if (is_null(bin.vars)) {
            bin.vars <- rep(FALSE, ncol(mat))
            bin.vars[to.sd] <- apply(mat[subset, to.sd,drop = FALSE], 2, is_binary)
        }
        else if (!is.atomic(bin.vars) || anyNA(as.logical(bin.vars)) ||
                 length(bin.vars) != ncol(mat)) {
            stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
        }
        
        possibly.supplied <- c("mat", "treat", "weighted.weights", "s.weights", "subset")
        lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                            possibly.supplied)
        supplied <- lengths > 0
        if (!all_the_same(lengths[supplied])) {
            stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
        }
        
        if (lengths["weighted.weights"] == 0) weighted.weights <- rep(1, NROW(mat))
        if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
        if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
        else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
        
        unique.treats <- if (is_(treat, "processed.treat") && all(subset)) as.character(treat_vals(treat)) else as.character(unique(treat[subset]))
        
        s.d.denom <- get.s.d.denom(as.character(s.d.denom), weights = weighted.weights[subset], treat = treat[subset])
        
        if (s.d.denom %in% c("treated", "control")) s.d.denom <- treat_vals(treat)[treat_names(treat)[s.d.denom]]
        
        treat <- as.character(treat)
        
        if (s.d.denom %in% unique.treats) denoms[to.sd] <- sqrt(col.w.v(mat[treat == s.d.denom, to.sd, drop = FALSE], 
                                                                    w = s.weights[treat == s.d.denom],
                                                                    bin.vars = bin.vars[to.sd], na.rm = na.rm))

        else if (s.d.denom == "pooled") denoms[to.sd] <- sqrt(Reduce("+", lapply(unique.treats, 
                                                                        function(t) col.w.v(mat[treat == t, to.sd, drop = FALSE], 
                                                                                            w = s.weights[treat == t],
                                                                                            bin.vars = bin.vars[to.sd], na.rm = na.rm))) / length(unique.treats))
        else if (s.d.denom == "all") denoms[to.sd] <- sqrt(col.w.v(mat[, to.sd, drop = FALSE], 
                                                               w = s.weights,
                                                               bin.vars = bin.vars[to.sd], na.rm = na.rm))
        else if (s.d.denom == "weighted") denoms[to.sd] <- sqrt(col.w.v(mat[, to.sd, drop = FALSE], 
                                                                    w = weighted.weights * s.weights,
                                                                    bin.vars = bin.vars[to.sd], na.rm = na.rm))
        else if (s.d.denom == "hedges") denoms[to.sd] <- (1 - 3/(4*length(treat) - 9))^-1 * sqrt(Reduce("+", lapply(unique.treats, 
                                                                                 function(t) (sum(treat == t) - 1) * col.w.v(mat[treat == t, to.sd, drop = FALSE], 
                                                                                                     w = s.weights[treat == t],
                                                                                                     bin.vars = bin.vars[to.sd], na.rm = na.rm))) / (length(treat) - 2))
        else stop("s.d.denom is not an allowed value.")
        
        if (any(zero_sds <- check_if_zero(denoms[to.sd]))) {
            denoms[to.sd][zero_sds] <- sqrt(col.w.v(mat[, to.sd, drop = FALSE][, zero_sds, drop = FALSE], 
                                                w = s.weights,
                                                bin.vars = bin.vars[to.sd][zero_sds], na.rm = na.rm))
        }
    }
    else {
        if (!is.numeric(s.d.denom) || length(s.d.denom) %nin% c(sum(to.sd), ncol(mat))) {
            stop("s.d.denom must be an allowable value or a numeric vector of with length equal to the number of columns of mat. See ?col_w_smd for allowable values.")
        }
        else if (length(s.d.denom) == sum(to.sd)) denoms[to.sd] <- s.d.denom
        else denoms <- s.d.denom
    }
    return(denoms)
}
baltal <- function(threshold) {
    #threshold: vector of threshold values (i.e., "Balanced"/"Not Balanced")
    threshnames <- names(table(threshold))
    balstring <- threshnames[nchar(threshnames) > 0][1]
    thresh.val <- substring(balstring, 1 + regexpr("[><]", balstring), nchar(balstring))
    b <- data.frame(count=c(sum(threshold==paste0("Balanced, <", thresh.val)), 
                            sum(threshold==paste0("Not Balanced, >", thresh.val))))
    rownames(b) <- c(paste0("Balanced, <", thresh.val), paste0("Not Balanced, >", thresh.val))
    return(b)
}
max.imbal <- function(balance.table, col.name, thresh.col.name, abs_stat) {
    balance.table.clean <- balance.table[balance.table$Type != "Distance" & is.finite(balance.table[, col.name]),]
    maxed <- balance.table.clean[which.max(abs_stat(balance.table.clean[, col.name])), match(c(col.name, thresh.col.name), names(balance.table.clean))]
    maxed <- data.frame(Variable = rownames(maxed), maxed)
    return(maxed)
}

balance.table <- function(C, type, weights = NULL, treat, continuous, binary, s.d.denom, 
                          thresholds = list(),
                          un = FALSE, disp.means = FALSE, disp.sds = FALSE, stats = NULL, 
                          s.weights = rep(1, length(treat)), abs = FALSE, no.adj = FALSE, 
                          var_types = NULL, s.d.denom.list = NULL, quick = TRUE) {
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    if (no.adj) weight.names <- "Adj"
    else weight.names <- names(weights)
    
    if (is_not_null(s.d.denom.list)) names(s.d.denom.list) <- weight.names
    if (is_not_null(s.d.denom)) names(s.d.denom) <- weight.names
    
    disp <- c("means"[disp.means], "sds"[disp.sds], all_STATS(type)[all_STATS(type) %in% stats])
    compute <- if (quick) disp else c("means", "sds", all_STATS(type))
    
    #B=Balance frame
    Bnames <- c("Type", 
                expand.grid_string(c(if (type == "bin") expand.grid_string(c("M"["means" %in% compute], "SD"["sds" %in% compute]), c("0", "1"), collapse = ".")
                                     else if (type == "cont") c("M"["means" %in% compute], "SD"["sds" %in% compute]), 
                                     unlist(lapply(compute[compute %in% all_STATS(type)], function(s) {
                                         c(STATS[[s]]$bal.tab_column_prefix,
                                           STATS[[s]]$Threshold)
                                     }))
                ),
                c("Un", weight.names), collapse = "."))
    B <- as.data.frame(matrix(nrow = NCOL(C), ncol = length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- colnames(C)
    
    #Set var type (binary/continuous)
    B[["Type"]] <- if_null_then(var_types, get.types(C))
    bin.vars <- B[["Type"]] == "Binary"
    
    #Means for each group
    if ("means" %in% compute) {
        if (type == "bin") {
            tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
            
            if (un || !quick) {
                for (t in c("0", "1")) {
                    B[[paste.("M", t, "Un")]] <- col_w_mean(C[treat == tn01[t], , drop = FALSE], weights = NULL,
                                                            s.weights = s.weights[treat==tn01[t]])
                }
            }
            
            if (!no.adj) {
                for (i in weight.names) {
                    for (t in c("0", "1")) {
                        B[[paste.("M", t, i)]] <- col_w_mean(C[treat == tn01[t], , drop = FALSE], weights = weights[[i]][treat==tn01[t]],
                                                             s.weights = s.weights[treat==tn01[t]])
                    }
                }
            }
        }
        else if (type == "cont") {
            if (un || !quick) {
                B[["M.Un"]] <- col_w_mean(C, weights = NULL, s.weights = s.weights)
            }
            if (!no.adj && !quick) {
                for (i in weight.names) {
                    B[[paste.("M", i)]] <- col_w_mean(C, weights = weights[[i]], s.weights = s.weights)
                }
            }
        }
    }
    
    #SDs for each group
    binary <- match_arg(binary, c("raw", "std"))
    if ("sds" %in% compute) {
        sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else !bin.vars
        if (type == "bin") {
            if (un || !quick) {
                for (t in c("0", "1")) {
                    sds <- rep(NA_real_, NCOL(C))
                    if (any(sd.computable)) {
                        sds[sd.computable] <- col_w_sd(C[treat == tn01[t], sd.computable, drop = FALSE],
                                                       weights = NULL, s.weights = s.weights[treat==tn01[t]],
                                                       bin.vars = bin.vars[sd.computable])
                    }
                    B[[paste.("SD", t, "Un")]] <- sds
                }
            }
            
            if (!no.adj && !quick) {
                for (i in weight.names) {
                    for (t in c("0", "1")) {
                        sds <- rep(NA_real_, NCOL(C))
                        if (any(sd.computable)) {
                            sds[sd.computable] <- col_w_sd(C[treat == tn01[t], sd.computable, drop = FALSE],
                                                           weights = weights[[i]][treat==tn01[t]], s.weights = s.weights[treat==tn01[t]],
                                                           bin.vars = bin.vars[sd.computable])
                        }
                        B[[paste.("SD", t, i)]] <- sds
                    }
                }
            }
        }
        else if (type == "cont") {
            if (un || !quick) {
                sds <- rep(NA_real_, NCOL(C))
                if (any(sd.computable)) {
                    sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE],
                                                   weights = NULL, s.weights = s.weights,
                                                   bin.vars = bin.vars[sd.computable])
                }
                B[["SD.Un"]] <- sds
            }
            if (!no.adj && !quick) {
                for (i in weight.names) {
                    sds <- rep(NA_real_, NCOL(C))
                    if (any(sd.computable)) {
                        sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE],
                                                       weights = weights[[i]], s.weights = s.weights,
                                                       bin.vars = bin.vars[sd.computable])
                    }
                    B[[paste.("SD", i)]] <- sds
                }
            }
        }
        if (all(vapply(B[startsWith(names(B), "SD.")], function(x) all(!is.finite(x)), logical(1L)))) {
            disp <- disp[disp != "sds"]
        }
    }
    
    for (s in all_STATS(type)) {
        if (s %in% compute) {
            if (!quick || un) {
                B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]] <- STATS[[s]]$fun(C, treat = treat, weights = NULL, 
                                                                                      std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                                      s.d.denom = if_null_then(s.d.denom.list[[1]], s.d.denom[1]),
                                                                                      abs = abs, s.weights = s.weights, bin.vars = bin.vars,
                                                                                      weighted.weights = weights[[1]])
            }
            
            if (!no.adj) {
                for (i in weight.names) {
                    B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]] <- STATS[[s]]$fun(C, treat = treat, weights = weights[[i]],
                                                                                       std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                                       s.d.denom = if_null_then(s.d.denom.list[[i]], s.d.denom[i]),
                                                                                       abs = abs, s.weights = s.weights, bin.vars = bin.vars)
                }
            }
            
            
            if (all(vapply(B[paste.(STATS[[s]]$bal.tab_column_prefix, c("Un", weight.names))], 
                           function(x) all(!is.finite(x)), logical(1L)))) {
                disp <- disp[disp != s] 
                thresholds[[s]] <- NULL
            }
            
            if (is_not_null(thresholds[[s]])) {
                if (no.adj) {
                    B[[paste.(STATS[[s]]$Threshold, "Un")]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]]), "", 
                                                                       STATS[[s]]$abs(B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                                       paste0("Not Balanced, >", round(thresholds[[s]], 3)))
                }
                else {
                    for (i in weight.names) {
                        B[[paste.(STATS[[s]]$Threshold, i)]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]]), "", 
                                                                        STATS[[s]]$abs(B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                                        paste0("Not Balanced, >", round(thresholds[[s]], 3)))
                    }
                }
            }
            if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
        }
    }
    
    attr(B, "thresholds") <- thresholds
    attr(B, "disp") <- disp
    attr(B, "compute") <- compute
    
    return(B)
}

samplesize <- function(treat, type, weights = NULL, subclass = NULL, s.weights = NULL, method=c("matching", "weighting", "subclassification"), discarded = NULL) {
    #Computes sample size info. for unadjusted and adjusted samples.
    # method is what method the weights are to be used for. 
    # method="subclassification" is for subclass sample sizes only.
    
    if (is_null(s.weights)) s.weights <- rep(1, length(treat))
    if (is_null(discarded)) discarded <- rep(FALSE, length(treat))
    
    if (type == "bin") {
        if (length(method) == 1 && method == "subclassification") {
            if (is_null(subclass)) stop("subclass must be a vector of subclasses.")
            
            nn <- matrix(0, nrow = length(treat_names(treat)) + 1, ncol = 1 + nlevels(subclass))
            
            nn[, 1 + nlevels(subclass)] <- c(vapply(treat_vals(treat), function(tn) sum(treat==tn), numeric(1L)), length(treat))
            
            matched <- !is.na(subclass)
            k <- 1
            for (i in levels(subclass)) {
                qt <- treat[matched & subclass == i]
                for (tnn in names(treat_names(treat))) {
                    if (sum(qt==treat_vals(treat)[treat_names(treat)[tnn]]) < 2)
                        warning(paste0("Not enough ", tnn, " units in subclass ", i, "."), call. = FALSE)
                }
                nn[, k] <- c(vapply(treat_vals(treat), function(tn) sum(qt==tn), numeric(1L)), length(qt))
                k <- k + 1 #Use a counter because subclass names may not be numbers
            }
            nn <- as.data.frame.matrix(nn, optional = TRUE)
            rownames(nn) <- c(treat_names(treat), "Total")
            names(nn) <- c(levels(subclass), "All")
            attr(nn, "tag") <- "Sample sizes by subclass"
        }
        else {
            if (is_null(weights)) {
                
                nn <- as.data.frame(matrix(0, ncol = length(treat_vals(treat)), nrow = 1))
                nn[1, ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
                dimnames(nn) <- list(c("All"), 
                                     c(treat_names(treat)))
                if (nunique.gt(s.weights, 2) || !any(s.weights==1) || any(s.weights %nin% c(0,1))) {
                    attr(nn, "ss.type") <- c("ess")
                }
                else {
                    attr(nn, "ss.type") <- c("ss")
                }
                
            }
            else if (NCOL(weights) == 1) {
                if (method=="matching") {
                    nn <- as.data.frame(matrix(0, ncol=length(treat_vals(treat)), nrow=5))
                    nn[1, ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn), numeric(1L))
                    nn[2, ] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, 1]), numeric(1L))
                    nn[3, ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & weights[,1] > 0), numeric(1L))
                    nn[4, ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & weights[,1]==0 & !discarded), numeric(1L))
                    nn[5, ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & weights[,1]==0 & discarded), numeric(1L))
                    dimnames(nn) <- list(c("All", "Matched (ESS)", "Matched (Unweighted)", "Unmatched", "Discarded"), 
                                         c(treat_names(treat)))
                    
                    attr(nn, "ss.type") <- rep("ss", NROW(nn))
                    
                    if (!any(discarded)) {
                        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                        nn <- nn[rownames(nn) != "Discarded", ,drop = FALSE]
                    }
                }
                else if (method == "weighting") {
                    nn <- as.data.frame(matrix(0, ncol = length(treat_vals(treat)), nrow = 3))
                    nn[1, ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
                    nn[2, ] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, 1]*s.weights[treat==tn]), numeric(1L))
                    nn[3, ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & discarded), numeric(1L))
                    dimnames(nn) <- list(c("Unadjusted", "Adjusted", "Discarded"), 
                                         c(treat_names(treat)))
                    attr(nn, "ss.type") <- c("ss", "ess", "ss")
                    
                    if (!any(discarded)) {
                        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                        nn <- nn[rownames(nn) != "Discarded", ,drop = FALSE]
                    }
                }
            }
            else {
                nn <- as.data.frame(matrix(0, ncol = length(treat_vals(treat)), nrow = 1 + NCOL(weights)))
                nn[1, ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
                for (i in seq_len(NCOL(weights))) {
                    if (method[i] == "matching") {
                        nn[1+i,] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, i]), numeric(1L))
                    }
                    else if (method[i] == "weighting") {
                        nn[1+i,] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, i]*s.weights[treat==tn]), numeric(1L))
                    }
                    
                }
                dimnames(nn) <- list(c("All", names(weights)), 
                                     treat_names(treat))
                attr(nn, "ss.type") <- c("ss", rep("ess", length(method)))
                
            }
            if (length(attr(nn, "ss.type")) > 1 && all(attr(nn, "ss.type")[-1] == "ess")) {
                attr(nn, "tag") <- "Effective sample sizes"
            }
            else attr(nn, "tag") <- "Sample sizes"
        }
    }
    else if (type == "cont") {
        if (length(method) == 1 && method == "subclassification") {
            if (is_null(subclass)) stop("subclass must be a vector of subclasses.")
            
            nn <- matrix(0, nrow = 1, ncol = 1 + nlevels(subclass))
            
            nn[, 1 + nlevels(subclass)] <- length(treat)
            
            matched <- !is.na(subclass)
            k <- 1
            for (i in levels(subclass)) {
                qt <- treat[matched & subclass == i]
                if (length(qt) < 2)
                    warning(paste0("Not enough units in subclass ", i, "."), call. = FALSE)
                
                nn[, k] <- length(qt)
                k <- k + 1 #Use a counter because subclass names may not be numbers
            }
            nn <- as.data.frame.matrix(nn, optional = TRUE)
            rownames(nn) <- c("Total")
            names(nn) <- c(levels(subclass), "All")
            attr(nn, "tag") <- "Sample sizes by subclass"
        }
        else {
            if (is_null(weights)) {
                
                nn <- as.data.frame(matrix(0, ncol = 1, nrow = 1))
                nn[1, ] <- ESS(s.weights)
                dimnames(nn) <- list(c("All"), 
                                     c("Total"))
                if (nunique.gt(s.weights, 2) || !any(s.weights==1) || any(s.weights %nin% c(0,1))) {
                    attr(nn, "ss.type") <- c("ess")
                }
                else {
                    attr(nn, "ss.type") <- c("ss")
                }
                
            }
            else if (NCOL(weights) == 1) {
                if (method=="matching") {
                    nn <- as.data.frame(matrix(0, ncol=1, nrow=5))
                    nn[1, ] <- length(treat)
                    nn[2, ] <- ESS(weights[, 1])
                    nn[3, ] <- sum(weights[,1] > 0 & !discarded)
                    nn[4, ] <- sum(weights[,1]==0 & !discarded)
                    nn[5, ] <- sum(discarded)
                    dimnames(nn) <- list(c("All", "Matched (ESS)", "Matched (Unweighted)", "Unmatched", "Discarded"), 
                                         c("Total"))
                    
                    attr(nn, "ss.type") <- rep("ss", NROW(nn))
                    
                    if (!any(discarded)) {
                        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                        nn <- nn[rownames(nn) != "Discarded", ,drop = FALSE]
                    }
                }
                else if (method == "weighting") {
                    nn <- as.data.frame(matrix(0, ncol = 1, nrow = 3))
                    nn[1, ] <- ESS(s.weights)
                    nn[2, ] <- ESS(weights[!discarded, 1]*s.weights[!discarded])
                    nn[3, ] <- sum(discarded)
                    dimnames(nn) <- list(c("Unadjusted", "Adjusted", "Discarded"), 
                                         c("Total"))
                    attr(nn, "ss.type") <- c("ss", "ess", "ss")
                    
                    if (!any(discarded)) {
                        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                        nn <- nn[rownames(nn) != "Discarded", ,drop = FALSE]
                    }
                }
            }
            else {
                nn <- as.data.frame(matrix(0, ncol = 1, nrow = 1 + NCOL(weights)))
                nn[1, ] <- ESS(s.weights)
                for (i in seq_len(NCOL(weights))) {
                    if (method[i] == "matching") {
                        nn[1+i,] <- ESS(weights[!discarded, i])
                    }
                    else if (method[i] == "weighting") {
                        nn[1+i,] <- ESS(weights[!discarded, i]*s.weights[!discarded])
                    }
                    
                }
                dimnames(nn) <- list(c("All", names(weights)), 
                                     "Total")
                attr(nn, "ss.type") <- c("ss", rep("ess", length(method)))
                
            }
            if (length(attr(nn, "ss.type")) > 1 && all(attr(nn, "ss.type")[-1] == "ess")) {
                attr(nn, "tag") <- "Effective sample sizes"
            }
            else attr(nn, "tag") <- "Sample sizes"
        }
    }
    return(nn)
}

balance.summary <- function(bal.tab.list, agg.funs, include.times = FALSE) {
    type <- attr(bal.tab.list[[1]], "print.options")[["type"]]
    disp <- attr(bal.tab.list[[1]], "print.options")[["disp"]]
    compute <- attr(bal.tab.list[[1]], "print.options")[["compute"]]
    thresholds <- attr(bal.tab.list[[1]], "print.options")[["thresholds"]]
    quick <- attr(bal.tab.list[[1]], "print.options")[["quick"]]
    weight.names <- attr(bal.tab.list[[1]], "print.options")[["weight.names"]]
    abs <- attr(bal.tab.list[[1]], "print.options")[["abs"]]
    no.adj <- attr(bal.tab.list[[1]], "print.options")[["nweights"]] == 0
    
    balance.list <- clear_null(lapply(bal.tab.list, function(x) x[["Balance"]]))
    
    Brownames <- unique(unlist(lapply(balance.list, rownames), use.names = FALSE))
    
    agg.funs <- tolower(agg.funs)
    Agg.Funs <- firstup(match_arg(agg.funs, c("min", "mean", "max"), several.ok = TRUE))
    all.Agg.Funs <- firstup(c("min", "mean", "max"))
    
    if (length(Agg.Funs) == 1 && Agg.Funs == "Max") abs <- TRUE
    
    Bcolnames <- c("Times", "Type", 
                   expand.grid_string(
                       c(unlist(lapply(compute, function(s) {
                           c(paste.(all.Agg.Funs, STATS[[s]]$bal.tab_column_prefix),
                             STATS[[s]]$Threshold)
                       }))), c("Un", weight.names), collapse = "."))
    
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
    names(B) <- Bcolnames
    
    B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(vapply(balance.list, function(y) if (x %in% rownames(y)) y[[x, "Type"]] else NA_character_, character(1))))), use.names = FALSE)
    
    if (include.times) B[["Times"]] <- vapply(Brownames, function(x) paste(seq_along(balance.list)[vapply(balance.list, function(y) x %in% rownames(y), logical(1L))], collapse = ", "), character(1))[Brownames]
    else B[["Times"]] <- NULL
    
    for (Agg.Fun in all.Agg.Funs) {
        if (Agg.Fun %in% Agg.Funs || !quick) {
            for (s in compute) {
                abs0 <- function(x) {if (is_null(x)) NA_real_ else if (abs) STATS[[s]]$abs(x) else (x)}
                agg <- function(x, ...) {
                    if (!any(is.finite(x))) NA_real_
                    else if (s == "variance.ratios" && tolower(Agg.Fun) == "mean") geom.mean(x)
                    else if (tolower(Agg.Fun) == "rms") sqrt(mean_fast(STATS[[s]]$abs(x)^2, TRUE))
                    else get(tolower(Agg.Fun))(x, ...)
                }
                for (sample in c("Un", weight.names)) {
                    if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
                        B[[paste.(Agg.Fun, STATS[[s]]$bal.tab_column_prefix, sample)]] <- vapply(Brownames, function(x) agg(unlist(lapply(balance.list, function(y) if (x %in% rownames(y)) abs0(y[[x, paste.(STATS[[s]]$bal.tab_column_prefix, sample)]]))), na.rm = TRUE), numeric(1))
                    }
                }
            }
        }
    }
    
    if (length(Agg.Funs) == 1) {
        #Assign X.Threshold values
        for (s in compute) {
            if (is_not_null(thresholds[[s]])) {
                if (no.adj) {
                    B[[paste.(STATS[[s]]$Threshold, "Un")]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix, "Un")]]), "", 
                                                                       STATS[[s]]$abs(B[[paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix, "Un")]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                                       paste0("Not Balanced, >", round(thresholds[[s]], 3)))
                }
                else {
                    for (i in weight.names) {
                        B[[paste.(STATS[[s]]$Threshold, i)]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix, i)]]), "", 
                                                                        STATS[[s]]$abs(B[[paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix, i)]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                                        paste0("Not Balanced, >", round(thresholds[[s]], 3)))
                    }
                }
            }
            if (no.adj || length(weight.names) <= 1) names(B)[names(B) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
        }
    }
    
    return(B)
}

#base.bal.tab.imp
samplesize.across.imps <- function(obs.list) {
    obs.list <- clear_null(obs.list)
    
    obs <- Reduce("+", obs.list)/length(obs.list)
    attr(obs, "tag") <- paste0("Average ", tolower(attr(obs.list[[1]], "tag")), " across imputations")
    return(obs)
}

#base.bal.tab.multi
samplesize.multi <- function(bal.tab.multi.list, treat_names, focal = NULL) {
    if (is_not_null(focal)) which <- c(treat_names[treat_names != focal], focal)
    else which <- treat_names
    bal.tab.multi.list <- clear_null(bal.tab.multi.list)
    obs <- do.call("cbind", unname(lapply(bal.tab.multi.list, function(x) x[["Observations"]])))[, which]
    attr(obs, "tag") <- attr(bal.tab.multi.list[[1]][["Observations"]], "tag")
    attr(obs, "ss.type") <- attr(bal.tab.multi.list[[1]][["Observations"]], "ss.type")
    return(obs)
}

#base.bal.tab.msm
samplesize.msm <- function(bal.tab.msm.list) {
    obs <- do.call("cbind", lapply(bal.tab.msm.list, function(x) x[["Observations"]]))
    attr(obs, "tag") <- attr(bal.tab.msm.list[[1]][["Observations"]], "tag")
    attr(obs, "ss.type") <- attr(bal.tab.msm.list[[1]][["Observations"]], "ss.type")
    return(obs)
}

#base.bal.tab.cluster
samplesize.across.clusters <- function(obs.list) {
    obs.list <- clear_null(obs.list)
    obs <- Reduce("+", obs.list)
    attr(obs, "tag") <- paste0("Total ", tolower(attr(obs.list[[1]], "tag")), " across clusters")
    return(obs)
}

#base.bal.tab.subclass
balance.table.subclass <- function(C, type, weights = NULL, treat, subclass,
                                       continuous, binary, s.d.denom, 
                                       thresholds = list(),
                                       un = FALSE, disp.means = FALSE, disp.sds = FALSE, stats = NULL, 
                                       s.weights = rep(1, length(treat)), abs = FALSE, var_types = NULL, quick = TRUE) {
    #Creates list SB of balance tables for each subclass
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    disp <- c("means"[disp.means], "sds"[disp.sds], all_STATS(type)[all_STATS(type) %in% stats])
    compute <- if (quick) disp else c("means", "sds", all_STATS(type))
    
    #B=Balance frame
    Bnames <- c("Type", 
                expand.grid_string(c(if (type == "bin") expand.grid_string(c("M"["means" %in% compute], "SD"["sds" %in% compute]), c("0", "1"), collapse = ".")
                                     else if (type == "cont") c("M"["means" %in% compute], "SD"["sds" %in% compute]), 
                                     unlist(lapply(compute[compute %in% all_STATS(type)], function(s) {
                                         c(STATS[[s]]$bal.tab_column_prefix,
                                           STATS[[s]]$Threshold)
                                     }))
                ),
                c("Adj"), collapse = "."))
    B <- as.data.frame(matrix(nrow = NCOL(C), ncol = length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- colnames(C)
    
    #Set var type (binary/continuous)
    B[["Type"]] <- if_null_then(var_types, get.types(C))
    bin.vars <- B[["Type"]] == "Binary"
    
    SB <- setNames(lapply(levels(subclass), function(i) B), levels(subclass))
    
    binary <- match_arg(binary, c("raw", "std"))
    sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else !bin.vars
    
    #-------------------------------------
    for (i in levels(subclass)) {
        
        in.subclass <- !is.na(subclass) & subclass==i
        
        #Means for each group
        if ("means" %in% compute) {
            if (type == "bin") {
                tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
                for (t in c("0", "1")) {
                    SB[[i]][[paste.("M", t, "Adj")]] <- col_w_mean(C, subset = treat==tn01[t] & in.subclass)
                }
            }
            else if (type == "cont") {
                SB[[i]][["M.Adj"]] <- col_w_mean(C, subset = in.subclass)
            }
        }
        
        #SDs for each group
        if ("sds" %in% compute) {
            if (type == "bin") {
                for (t in c("0", "1")) {
                    sds <- rep(NA_real_, NCOL(C))
                    sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == tn01[t] & in.subclass)
                    SB[[i]][[paste.("SD", t, "Adj")]] <- sds
                }
            }
            else if (type == "cont") {
                sds <- rep(NA_real_, NCOL(C))
                sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == in.subclass)
                SB[[i]][["SD.Adj"]] <- sds
            }
        }
        
        for (s in all_STATS(type)) {
            if (s %in% compute) {
                SB[[i]][[paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")]] <- STATS[[s]]$fun(C, treat = treat, weights = NULL, 
                                                                                             std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                                             s.d.denom = s.d.denom,
                                                                                             abs = FALSE, s.weights = s.weights, 
                                                                                             bin.vars = bin.vars, subset = in.subclass)
                
                if (all(vapply(SB[[i]][paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")], 
                               function(x) all(!is.finite(x)), logical(1L)))) {
                    disp <- disp[disp != s] 
                    thresholds[[s]] <- NULL
                }
                
                
                if (is_not_null(thresholds[[s]])) {
                    
                    SB[[i]][[paste.(STATS[[s]]$Threshold, "Adj")]] <- ifelse(SB[[i]][["Type"]]!="Distance" & is.finite(SB[[i]][[paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")]]), 
                                                                             paste0(ifelse(abs_(SB[[i]][[paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")]]) < thresholds[[s]], "Balanced, <", "Not Balanced, >"), round(thresholds[[s]], 3)), "")
                }
                names(SB[[i]])[names(SB[[i]]) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
            }
            
        }

    }
    
    if (all(vapply(SB, function(sb) all(vapply(sb[startsWith(names(sb), "SD.")], function(x) all(!is.finite(x)), logical(1L))), logical(1L)))) {
        disp <- disp[disp != "sds"]
    }

    attr(SB, "thresholds") <- thresholds
    attr(SB, "disp") <- disp
    attr(SB, "compute") <- compute
    
    return(SB)
}
balance.table.across.subclass.cont <- function(balance.table, balance.table.subclass.list, subclass.obs, r.threshold = NULL) {
    #Not specified
    
    B.A <- balance.table.subclass.list[[1]][c("M.Adj", "SD.Adj", "Corr.Adj")]
    
    for(i in rownames(B.A)) {
        for(j in colnames(B.A)) {
            if (startsWith(j, "SD.")) {
                B.A[[i, j]] <- sqrt(sum(vapply(seq_along(balance.table.subclass.list),
                                               function(s) subclass.obs[[s]]/sum(subclass.obs) * (balance.table.subclass.list[[s]][[i, j]]^2), numeric(1))))
            }
            else {
                B.A[[i, j]] <- sum(vapply(seq_along(balance.table.subclass.list),
                                          function(s) subclass.obs[[s]]/sum(subclass.obs) * (balance.table.subclass.list[[s]][[i, j]]), numeric(1)))
                
            }
        }
    }
    B.A.df <- data.frame(balance.table[c("Type", "M.Un", "SD.Un", "Corr.Un", "R.Threshold.Un")], 
                         B.A, R.Threshold = NA_character_)
    if (is_not_null(r.threshold)) B.A.df[["R.Threshold"]] <- ifelse(B.A.df[["Type"]]=="Distance", "", paste0(ifelse(is.finite(B.A.df[["Corr.Adj"]]) & abs_(B.A.df[["Corr.Adj"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold))
    return(B.A.df)
}

#base.bal.tab.target
samplesize.target <- function(bal.tab.target.list, treat_names, target.name) {
    which <- treat_names[treat_names != target.name]
    obs <- do.call("cbind", unname(lapply(bal.tab.target.list, function(x) x[["Observations"]])))[, which]
    attr(obs, "tag") <- attr(bal.tab.target.list[[1]][["Observations"]], "tag")
    attr(obs, "ss.type") <- attr(bal.tab.target.list[[1]][["Observations"]], "ss.type")
    return(obs)
}

#love.plot
isColor <- function(x) {
    tryCatch(is.matrix(col2rgb(x)), 
             error = function(e) FALSE)
}
f.recode <- function(f, ...) {
    #Simplified version of forcats::fct_recode
    f <- factor(f)
    new_levels <- unlist(list(...), use.names = TRUE)
    old_levels <- levels(f)
    idx <- match(new_levels, old_levels)
    
    old_levels[idx] <- names(new_levels)
    
    levels(f) <- old_levels
    return(f)
}
seq_int_cycle <- function(begin, end, max) {
    seq(begin, end, by = 1) - max*(seq(begin-1, end-1, by = 1) %/% max)
}
assign.shapes <- function(colors, default.shape = "circle") {
    if (nunique(colors) < length(colors)) {
        shapes <- seq_int_cycle(19, 19 + length(colors) - 1, max = 25)
    }
    else shapes <- rep.int(default.shape, length(colors))
    return(shapes)
}
shapes.ok <- function(shapes, nshapes) {
    shape_names <- c(
        "circle", paste("circle", c("open", "filled", "cross", "plus", "small")), "bullet",
        "square", paste("square", c("open", "filled", "cross", "plus", "triangle")),
        "diamond", paste("diamond", c("open", "filled", "plus")),
        "triangle", paste("triangle", c("open", "filled", "square")),
        paste("triangle down", c("open", "filled")),
        "plus", "cross", "asterisk"
    )
    shape_nums <- 1:25
    return((length(shapes) == 1 || length(shapes) == nshapes) && ((is.numeric(shapes) && all(shapes %in% shape_nums)) || (is.character(shapes) && all(shapes %in% shape_names))))
}
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
ggarrange_simple <- function (plots, nrow = NULL, ncol = NULL) {
    #A thin version of egg:ggarrange
    
    gtable_frame <- function (g, width = unit(1, "null"), height = unit(1, "null")) {
        panels <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]),]
        pargins <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]),]
        ll <- unique(panels$l)
        margins <- if (length(ll) == 1) unit(0, "pt") else g$widths[ll[-length(ll)] + 2]
        tt <- unique(panels$t)
        fixed_ar <- g$respect
        if (fixed_ar) {
            ar <- as.numeric(g$heights[tt[1]])/as.numeric(g$widths[ll[1]])
            height <- width * (ar/length(ll))
            g$respect <- FALSE
        }
        core <- g[seq(min(tt), max(tt)), seq(min(ll), max(ll))]
        top <- g[seq(1, min(tt) - 1), seq(min(ll), max(ll))]
        bottom <- g[seq(max(tt) + 1, nrow(g)), seq(min(ll), max(ll))]
        left <- g[seq(min(tt), max(tt)), seq(1, min(ll) - 1)]
        right <- g[seq(min(tt), max(tt)), seq(max(ll) + 1, ncol(g))]
        fg <- grid::nullGrob()
        if (length(left)) {
            lg <- gtable::gtable_add_cols(left, unit(1, "null"), 0)
            lg <- gtable::gtable_add_grob(lg, fg, 1, l = 1)
        }
        else {
            lg <- fg
        }
        if (length(right)) {
            rg <- gtable::gtable_add_cols(right, unit(1, "null"))
            rg <- gtable::gtable_add_grob(rg, fg, 1, l = ncol(rg))
        }
        else {
            rg <- fg
        }
        if (length(top)) {
            tg <- gtable::gtable_add_rows(top, unit(1, "null"), 0)
            tg <- gtable::gtable_add_grob(tg, fg, t = 1, l = 1)
        }
        else {
            tg <- fg
        }
        if (length(bottom)) {
            bg <- gtable::gtable_add_rows(bottom, unit(1, "null"), 
                                          -1)
            bg <- gtable::gtable_add_grob(bg, fg, t = nrow(bg), l = 1)
        }
        else {
            bg <- fg
        }
        grobs <- list(fg, tg, fg, lg, core, rg, fg, bg, fg)
        widths <- grid::unit.c(sum(left$widths), width, sum(right$widths))
        heights <- grid::unit.c(sum(top$heights), height, sum(bottom$heights))
        all <- gtable::gtable_matrix("all", grobs = matrix(grobs, ncol = 3, nrow = 3, byrow = TRUE), 
                                     widths = widths, heights = heights)
        
        all[["layout"]][5, "name"] <- "panel"
        if (fixed_ar) 
            all$respect <- TRUE
        all
    }
    
    n <- length(plots)
    
    grobs <- lapply(plots, ggplot2::ggplotGrob)
    
    if (is_null(nrow) && is_null(ncol)) {
        nm <- grDevices::n2mfrow(n)
        nrow <- nm[1]
        ncol <- nm[2]
    }
    
    hw <- lapply(rep(1, n), unit, "null")
    
    fg <- lapply(seq_along(plots), function(i) gtable_frame(g = grobs[[i]], 
                                                            width = hw[[i]], height = hw[[i]]))
    
    spl <- split(fg, rep(1, n))
    
    rows <- lapply(spl, function(r) do.call(gridExtra::gtable_cbind, r))
    
    gt <- do.call(gridExtra::gtable_rbind, rows)
    
    invisible(gt)
}
unpack_bal.tab <- function(b) {
    unpack_bal.tab_internal <- function(b) {
        if (is_(b, c("bal.tab.bin", "bal.tab.cont"))) return(b[["Balance"]])
        else {
            b_ <- b[[which(endsWith(names(b), ".Balance"))]]
            
            b_list <- lapply(b_, function(i) {
                if (is_(i, c("bal.tab.bin", "bal.tab.cont"))) return(i[["Balance"]])
                else return(unpack_bal.tab_internal(i))
            })
            return(b_list)
        }
    }
    LinearizeNestedList <- function(NList, NameSep) {
        # LinearizeNestedList:
        #
        # https://sites.google.com/site/akhilsbehl/geekspace/
        #         articles/r/linearize_nested_lists_in_r
        #
        # Akhil S Bhel
        # 
        
        if (is.data.frame(NList)) return(NList)
        
        A <- 1
        B <- length(NList)
        
        while (A <= B) {
            Element <- NList[[A]]
            EName <- names(NList)[A]
            if (is.list(Element)) {
                
                if (A == 1) {
                    Before <- NULL
                } else {
                    Before <- NList[1:(A - 1)]
                }
                if (A == B) {
                    After <- NULL
                } else {
                    After <- NList[(A + 1):B]
                }
                
                if (is.data.frame(Element)) {
                    Jump <- 1
                } else {
                    NList[[A]] <- NULL
                    
                    Element <- LinearizeNestedList(Element, NameSep)
                    names(Element) <- paste(EName, names(Element), sep=NameSep)
                    Jump <- length(Element)
                    NList <- c(Before, Element, After)
                }
            } else {
                Jump <- 1
            }
            A <- A + Jump
            B <- length(NList)
        }
        return(NList)
    }
    class_sequence <- function(b) {
        if (is_(b, c("bal.tab.bin", "bal.tab.cont"))) return(NULL)
        else {
            b_ <- b[[which(endsWith(names(b), ".Balance"))]][[1]]
            return(c(class(b)[1], class_sequence(b_)))
        }
    }
    
    namesep <- paste(c("|", do.call(c, lapply(1:20, function(i) sample(LETTERS, 1))), "|"), collapse = "")
    
    out_ <- unpack_bal.tab_internal(b)
    out <- LinearizeNestedList(out_, NameSep = namesep)
    
    attr(out, "namesep") <- namesep
    attr(out, "class_sequence") <- class_sequence(b)
    
    return(out)
}

#bal.plot
get.var.from.list.with.time <- function(var.name, covs.list) {
    var.name.in.covs <- vapply(covs.list, function(x) var.name %in% names(x), logical(1))
    var <- unlist(lapply(covs.list[var.name.in.covs], function(x) x[[var.name]]))
    times <- rep(var.name.in.covs, each = NCOL(covs.list[[1]]))
    return(list(var = var, times = times))
}

#print.bal.tab
print.data.frame_ <- function(x, ...) {
    if (is_not_null(x) && NROW(x) > 0 && NCOL(x) > 0) {
        print.data.frame(x, ...)
    }
}

#Balance summary
process.bin.vars <- function(bin.vars, mat) {
    if (missing(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
    else if (is_null(bin.vars)) bin.vars <- rep(FALSE, ncol(mat))
    else {
        if (is.logical(bin.vars)) {
            bin.vars[is.na(bin.vars)] <- FALSE
            if (length(bin.vars) != ncol(mat)) stop("If bin.vars is logical, it must have length equal to the number of columns of mat.")
        }
        else if (is.numeric(bin.vars)) {
            bin.vars <- bin.vars[!is.na(bin.vars) & bin.vars != 0]
            if (any(bin.vars < 0) && any(bin.vars > 0)) stop("Positive and negative indices cannot be mixed with bin.vars.")
            if (any(abs(bin.vars) > ncol(mat))) stop("If bin.vars is numeric, none of its values can exceed the number of columns of mat.")
            logical.bin.vars <- rep(any(bin.vars < 0), ncol(mat))
            logical.bin.vars[abs(bin.vars)] <- !logical.bin.vars[abs(bin.vars)]
            bin.vars <- logical.bin.vars
        }
        else if (is.character(bin.vars)) {
            bin.vars <- bin.vars[!is.na(bin.vars) & bin.vars != ""]
            if (is_null(colnames(mat))) stop("If bin.vars is character, mat must have column names.")
            if (any(bin.vars %nin% colnames(mat))) stop("If bin.vars is character, all its values must be column names of mat.")
            bin.vars <- colnames(mat) %in% bin.vars
        }
        else stop("bin.vars must be a logical, numeric, or character vector.")
    }
    return(bin.vars)
}

#set.cobalt.options
acceptable.options <- function() {
    TF <- c(TRUE, FALSE)
    return(list(stats = c("mean.diffs"),
                un = TF,
                continuous = c("raw", "std"),
                binary = c("raw", "std"),
                imbalanced.only = TF,
                disp.means = TF,
                disp.sds = TF,
                disp.v.ratio = TF,
                disp.ks = TF,
                disp.subclass = TF,
                disp.bal.tab = TF,
                cluster.summary = TF,
                cluster.fun = c("min", "mean", "max"),
                imp.summary = TF,
                imp.fun = c("min", "mean", "max"),
                multi.summary = TF,
                msm.summary = TF,
                target.summary = TF,
                int_sep = " * ",
                factor_sep = "_",
                center = TF))
}

#On attach
.onAttach <- function(...) {
    cobaltLib <- dirname(system.file(package = "cobalt"))
    version <- packageDescription("cobalt", lib.loc = cobaltLib)$Version
    BuildDate <- packageDescription("cobalt", lib.loc = cobaltLib)$Date
    
    foo <- paste0(" cobalt (Version ", version, ", Build Date: ", BuildDate, ")\n", 
                  "   Please read the documentation at ?bal.tab to understand the default outputs.\n",
                  "   Submit bug reports and feature requests to https://github.com/ngreifer/cobalt/issues\n",
                  "   Install the development version (not guaranteed to be stable) with:\n",
                  "     devtools::install_github(\"ngreifer/cobalt\")\n",
                  "   Thank you for using cobalt!")
    packageStartupMessage(foo)
}

.onLoad <- function(libname, pkgname) {
    backports::import(pkgname)
}

#To pass CRAN checks:
utils::globalVariables(c("distance", "addl", "addl.list", "distance.list",
                         "quick", "treat", "Sample", "min.stat",
                         "max.stat", "mean.stat", "count",
                         "cum.pt", "treat.mean", "var.mean"))