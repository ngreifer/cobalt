#bal.tab
is.designmatch <- function(x) {
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
is.time.list <- function(x) {
    if (is.vector(x, mode = "list")) {
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
    
    if (is_null(treat)) stop("treat must be specified.", call. = FALSE)
    
    if (inherits(treat, "unprocessed.treat")) {
        attrs <- attributes(treat)
        renamed_original <- setNames(names(treat_vals(treat)), treat_vals(treat))
        treat <- factor(renamed_original[as.character(treat)], levels = renamed_original)
        for (at in c("treat_names", "treat_vals", "treat.type", "names"))
            attr(treat, at) <- attrs[[at]]
    }
    else {
        bad.treat <- FALSE
        if (is.character(treat) && length(treat)==1L && is_not_null(data)) {
            if (is_(data, "list")) {
                for (i in seq_along(data)) {
                    if (is_(data[[i]], "matrix") && treat %in% colnames(data[[i]])) {
                        treat <- data[[i]][,treat]
                        break
                    }
                    else if (is_(data[[i]], "data.frame") && treat %in% names(data[[i]])) {
                        treat <- data[[i]][[treat]]
                        break
                    }
                    else if (i == length(data)) bad.treat <- TRUE
                }
            }
            else if (is_(data, "matrix") && treat %in% colnames(data)) {
                treat <- data[,treat]
            }
            else if (is_(data, "data.frame") && treat %in% names(data)) {
                treat <- data[[treat]]
            }
            else bad.treat <- TRUE
        }
        else if (is_(treat, c("numeric", "logical", "character", "factor")) && length(treat) > 1L) {
            treat <- treat
        }
        else {
            bad.treat <- TRUE
        }
        
        if (bad.treat) stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
        
        if (anyNA(treat)) stop("Missing values exist in treat.", call. = FALSE)
        
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
    if (!is.vector(treat.list, "list")) {
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
        if (is_not_null(covs)) {if (is_not_null(D)) D <- cbind(D, covs) else D <- covs}
        t.c <- get.covs.and.treat.from.formula(formula, D, treat = treat)
        t.c <- list(treat = t.c[["treat"]], covs = t.c[["reported.covs"]], treat.name = t.c[["treat.name"]])
        attr(t.c, "which") <- "fd"
    }
    else {
        if (is.matrix(covs)) covs <- as.data.frame(covs)
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
        if (is.vector(val, mode = "list")) {
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
                if (is.vector(val, mode = "list")) {
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
        allowable.s.d.denoms <- c("pooled", "all", "weighted")
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
            try.estimand <- tryCatch(match_arg(toupper(estimand), c("ATT", "ATC", "ATE"), several.ok = TRUE),
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
    
    if (!quietly) {
        if (s.d.denom.specified && bad.s.d.denom && (!estimand.specified || bad.estimand)) {
            message(paste0("Warning: s.d.denom should be one of ", word_list(unique(c(unique.treats, allowable.s.d.denoms)), "or", quotes = TRUE), 
                           ".\n         Using ", word_list(s.d.denom, quotes = TRUE), " instead."))
        }
        else if (estimand.specified && bad.estimand) {
            message(paste0("Warning: estimand should be one of ", word_list(c("ATT", "ATC", "ATE"), "or", quotes = TRUE), 
                           ". Ignoring estimand."))
        }
        else if ((check.focal || check.weights) && any(s.d.denom %nin% treat_names(treat))) {
            message("Note: s.d.denom not specified; assuming ", ifelse(all_the_same(s.d.denom), unique(s.d.denom), word_list(s.d.denom)), ".")
        }
    }
    
    if (is_not_null(weights) && length(s.d.denom) != NCOL(weights)) {
        stop("Valid inputs to s.d.denom or estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
    }
    
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
        message("Note: estimand not specified; assuming ", ifelse(all_the_same(toupper(estimand)), toupper(unique(estimand)), word_list(toupper(estimand))), ".")
    }
    
    if (is_not_null(weights) && length(estimand) != ncol(weights)) {
        stop("Valid inputs to estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
    }
    return(estimand)
}
get.X.class <- function(X) {
    if (is_not_null(X[["imp"]])) {
        if (is_not_null(X[["treat.list"]])) stop("Multiply imputed data is not yet supported with longitudinal treatments.", call. = FALSE)
        else if (get.treat.type(X[["treat"]]) == "multinomial") stop("Multiply imputed data is not yet supported with multinomial treatments.", call. = FALSE)
        else X.class <- "imp"
    }
    else if (is_not_null(X[["treat.list"]])) X.class <- "msm"
    else if (get.treat.type(X[["treat"]]) == "binary") X.class <- "binary"
    else if (get.treat.type(X[["treat"]]) == "multinomial") X.class <- "multi"
    else if (get.treat.type(X[["treat"]]) == "continuous") X.class <- "cont"
    else probably.a.bug()
    
    return(X.class)
}
subset_X <- function(X, subset = NULL) {
    if (is_not_null(subset)) {
        if (!any(subset)) {
            stop("All subset set to FALSE.", call. = FALSE)
        }
        else if (!all(subset)) {
            n <- length(subset)
            subset_X_internal <- function(x, subset) {
                if (is_not_null(x)) {
                    if (is.factor(x) && length(x) == n) factor(x[subset])
                    else if (is.atomic(x) && length(x) == n) x[subset]
                    else if ((is.matrix(x) || is.data.frame(x)) && nrow(x) == n) x[subset, , drop = FALSE]
                    else if (is.list(x)) lapply(x, subset_X_internal, subset = subset)
                    else x
                }
                else x
            }
            lapply(X, subset_X_internal, subset)
        }
        else X
    }
    else X
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
                      do.call("rbind", mylist))
    
    if (is.integer(attr(data$data, "row.names"))) 
        row.names(cmp) <- seq_len(nrow(cmp))
    else row.names(cmp) <- as.character(seq_len(nrow(cmp)))
    
    return(cmp)
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
    new.co.names <- vector("list", (nc))
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
                vars.w.missing[i, "placed.after"] <- newly.added.names[length(newly.added.names)]
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

#base.bal.tab
check_if_zero_weights <- function(weights.df, treat, treat.type = "cat") {
    #Checks if all weights are zero in each treat group for each set of weights
    if (treat.type == "cat") {
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
    else if (treat.type == "cont") {
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
    else stop("treat.type must be either \"cat\" or \"cont\".")
    
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
samplesize <- function(treat, weights = NULL, subclass = NULL, s.weights = NULL, method=c("matching", "weighting", "subclassification"), cluster = NULL, which.cluster = NULL, discarded = NULL) {
    #Computes sample size info. for unadjusted and adjusted samples.
    # method is what method the weights are to be used for. 
    # method="subclassification" is for subclass sample sizes only.
    
    if (is_not_null(cluster) && is_not_null(which.cluster)) in.cluster <- cluster == which.cluster
    else in.cluster <- rep(TRUE, length(treat))
    if (is_null(s.weights)) s.weights <- rep(1, length(treat))
    if (is_null(discarded)) discarded <- rep(FALSE, length(treat))
    
    if (length(method) == 1 && method == "subclassification") {
        if (is_null(subclass)) stop("subclass must be a vector of subclasses.")
        qbins <- nlevels(subclass)
        
        nn <- as.data.frame(matrix(0, nrow = length(treat_names(treat)) + 1, ncol = qbins))
        
        dimnames(nn) <- list(c(treat_names(treat), "Total"), 
                             paste("Subclass", levels(subclass)))
        
        matched <- !is.na(subclass)
        k <- 0
        for (i in levels(subclass)) {
            qi <- subclass[matched] == i
            qt <- treat[matched][qi]
            for (tnn in names(treat_names(treat))) {
                if (sum(qt==treat_vals(treat)[treat_names(treat)[tnn]]) < 2)
                    warning(paste0("Not enough ", tnn, " units in subclass ", i, "."), call. = FALSE)
            }
            k <- k + 1 #Use a counter because subclass names may not be numbers
            nn[, k] <- c(vapply(treat_vals(treat), function(tn) sum(qt==tn), numeric(1L)), length(qt))
        }
        attr(nn, "tag") <- "Sample sizes by subclass"
    }
    else if (is_null(weights)) {
        
        t <- treat[in.cluster]
        sw <- s.weights[in.cluster]
        
        nn <- as.data.frame(matrix(0, ncol = length(treat_vals(treat)), nrow = 1))
        nn[1, ] <- vapply(treat_vals(treat), function(tn) ESS(sw[t==tn]), numeric(1L))
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
            nn[1, ] <- vapply(treat_vals(treat), function(tn) sum(in.cluster & treat==tn), numeric(1L))
            nn[2, ] <- vapply(treat_vals(treat), function(tn) ESS(weights[in.cluster & treat==tn, 1]), numeric(1L))
            nn[3, ] <- vapply(treat_vals(treat), function(tn) sum(in.cluster & treat==tn & weights[,1] > 0), numeric(1L))
            nn[4, ] <- vapply(treat_vals(treat), function(tn) sum(in.cluster & treat==tn & weights[,1]==0 & discarded==0), numeric(1L))
            nn[5, ] <- vapply(treat_vals(treat), function(tn) sum(in.cluster & treat==tn & weights[,1]==0 & discarded==1), numeric(1L))
            dimnames(nn) <- list(c("All", "Matched (ESS)", "Matched (Unweighted)", "Unmatched", "Discarded"), 
                                 c(treat_names(treat)))
            
            attr(nn, "ss.type") <- rep("ss", NROW(nn))
            
            if (!any(discarded)) {
                attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                nn <- nn[rownames(nn) != "Discarded", ,drop = FALSE]
            }
        }
        else if (method == "weighting") {
            
            t <- treat[in.cluster]
            w <- weights[in.cluster, 1]
            sw <- s.weights[in.cluster]
            dc <- discarded[in.cluster]
            
            nn <- as.data.frame(matrix(0, ncol = length(treat_vals(treat)), nrow = 3))
            nn[1, ] <- vapply(treat_vals(treat), function(tn) ESS(sw[t==tn]), numeric(1L))
            nn[2, ] <- vapply(treat_vals(treat), function(tn) ESS(w[t==tn]*sw[t==tn]), numeric(1L))
            nn[3, ] <- vapply(treat_vals(treat), function(tn) sum(t==tn & dc==1), numeric(1L))
            dimnames(nn) <- list(c("Unadjusted", "Adjusted", "Discarded"), 
                                 c(treat_names(treat)))
            attr(nn, "ss.type") <- c("ss", "ess", "ss")
            
            if (!any(dc)) {
                attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                nn <- nn[rownames(nn) != "Discarded", ,drop = FALSE]
            }
        }
    }
    else {
        t <- treat[in.cluster]
        sw <- s.weights[in.cluster]
        
        nn <- as.data.frame(matrix(0, ncol=length(treat_vals(treat)), nrow=1+NCOL(weights)))
        nn[1, ] <- vapply(treat_vals(treat), function(tn) ESS(sw[t==tn]), numeric(1L))
        for (i in seq_len(NCOL(weights))) {
            if (method[i] == "matching") {
                w <- weights[in.cluster, i]
                nn[1+i,] <- vapply(treat_vals(treat), function(tn) ESS(w[t==tn]), numeric(1L))
            }
            else if (method[i] == "weighting") {
                w <- weights[in.cluster, i]
                nn[1+i,] <- vapply(treat_vals(treat), function(tn) ESS(w[t==tn]*sw[t==tn]), numeric(1L))
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
    return(nn)
}
samplesize.across.clusters <- function(samplesize.list) {
    samplesize.list <- clear_null(samplesize.list)
    obs <- Reduce("+", samplesize.list)
    attr(obs, "tag") <- paste0("Total ", tolower(attr(samplesize.list[[1]], "tag")), " across clusters")
    return(obs)
}
max.imbal <- function(balance.table, col.name, thresh.col.name, ratio = FALSE) {
    balance.table.clean <- balance.table[balance.table$Type != "Distance" & is.finite(balance.table[, col.name]),]
    maxed <- balance.table.clean[which.max(abs_(balance.table.clean[, col.name], ratio = ratio)), match(c(col.name, thresh.col.name), names(balance.table.clean))]
    maxed <- data.frame(Variable = rownames(maxed), maxed)
    return(maxed)
}
balance.table <- function(C, weights, treat, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.sds = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, 
                          s.weights = rep(1, length(treat)), abs = FALSE, no.adj = FALSE, types = NULL, addl.sds = NULL, disp.pop = FALSE, pop.means = NULL, pop.sds = NULL, quick = TRUE) {
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    if (no.adj) weight.names <- "Adj"
    else weight.names <- names(weights)
    names(s.d.denom) <- weight.names
    
    #B=Balance frame
    Bnames <- c("Type", 
                apply(expand.grid(c("M.0", "SD.0", "M.1", "SD.1", 
                                    # "M.Pop", "SD.Pop", 
                                    "Diff", "M.Threshold", 
                                    # "Diff.0.Pop", "Diff.1.Pop", 
                                    "V.Ratio", "V.Threshold", "KS", "KS.Threshold"),
                                  c("Un", weight.names)), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow = NCOL(C), ncol = length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- colnames(C)
    
    #Set var type (binary/continuous)
    if (is_not_null(types)) B[["Type"]] <- types
    else B[["Type"]] <- get.types(C)
    bin.vars <- B[["Type"]] == "Binary"
    tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
    
    #Means for each group
    # if (!((!un || !disp.means) && quick)) {
    for (t in c("0", "1")) {
        B[[paste.("M", t, "Un")]] <- col_w_mean(C[treat == tn01[t], , drop = FALSE], weights = NULL,
                                                s.weights = s.weights[treat==tn01[t]])
    }
    # if (!(!disp.pop && quick)) {
    #     B[["M.Pop.Un"]] <- col.w.m(C, w = s.weights)
    # }
    # }
    if (!no.adj) {
        # if (!no.adj && !(!disp.means && quick)) {
        for (i in weight.names) {
            for (t in c("0", "1")) {
                B[[paste.("M", t, i)]] <- col_w_mean(C[treat == tn01[t], , drop = FALSE], weights = weights[[i]][treat==tn01[t]],
                                                     s.weights = s.weights[treat==tn01[t]])
            }
            # if (!(!disp.pop && quick)) {
            #     B[[paste.("M.Pop.Un", i)]] <- col.w.m(C, w = weights[[i]]*s.weights)
            # }
        }
    }
    
    #SDs for each group
    if (missing(binary) || is_null(binary)) {
        binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
    }
    else binary <- match_arg(binary, c("raw", "std"))
    
    sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else !bin.vars
    # if (!((!un || !disp.sds) && quick)) {
    for (t in c("0", "1")) {
        sds <- rep(NA_real_, NCOL(C))
        if (any(sd.computable)) {
            sds[sd.computable] <- col_w_sd(C[treat == tn01[t], sd.computable, drop = FALSE],
                                           weights = NULL, s.weights = s.weights[treat==tn01[t]],
                                           bin.vars = bin.vars[sd.computable])
        }
        B[[paste.("SD", t, "Un")]] <- sds
    }
    # if (!(!disp.pop && quick)) {
    #     sds <- rep(NA_real_, NCOL(C))
    #     sds[non.binary] <- sqrt(col.w.v(C[, non.binary, drop = FALSE], w = s.weights))
    #     B[["SD.Pop.Un"]] <- sds
    # }
    # }
    # if (!no.adj && !(!disp.sds && quick)) {
    if (!no.adj) {
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
            # if (!(!disp.pop && quick)) {
            #     sds <- rep(NA_real_, NCOL(C))
            #     sds[non.binary] <- sqrt(col.w.v(C[, non.binary, drop = FALSE], w = weights[[i]]*s.weights))
            #     B[[paste.("SD.Pop", i)]] <- sds
            # }
        }
    }
    if (!any(sapply(B[startsWith(names(B), "SD.")], is.finite))) {disp.sds <- FALSE}
    
    #Mean differences
    # if (!(!un && quick)) #Always compute unadjusted diffs
    
    B[["Diff.Un"]] <- col_w_smd(C, treat = treat, weights = NULL,
                                std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                s.d.denom = if (is_null(addl.sds)) s.d.denom[1] else addl.sds[[s.d.denom[1]]],
                                abs = abs, s.weights = s.weights, bin.vars = bin.vars)
    
    if (!no.adj) {
        for (i in weight.names) {
            B[[paste.("Diff", i)]] <- col_w_smd(C, treat = treat, weights = weights[[i]],
                                                std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                s.d.denom = if (is_null(addl.sds[[s.d.denom[i]]])) s.d.denom[i] 
                                                else addl.sds[[s.d.denom[i]]],
                                                abs = abs, s.weights = s.weights, bin.vars = bin.vars)
        }
    }
    
    #Variance ratios
    if (!(!disp.v.ratio && quick)) {
        if (!(!un && quick)) {
            vrs <- rep(NA_real_, NCOL(C))
            if (any(!bin.vars)) {
                vrs[!bin.vars] <- col_w_vr(C[, !bin.vars, drop = FALSE], treat, weights = NULL, abs = abs, 
                                           s.weights = s.weights, bin.vars = bin.vars[!bin.vars])
            }
            B[["V.Ratio.Un"]] <- vrs
        }
        if (!no.adj) {
            for (i in weight.names) {
                vrs <- rep(NA_real_, NCOL(C))
                if (any(!bin.vars)) {
                    vrs[!bin.vars] <- col_w_vr(C[, !bin.vars, drop = FALSE], treat, weights = weights[[i]], abs = abs, 
                                               s.weights = s.weights, bin.vars = bin.vars[!bin.vars])
                }
                B[[paste.("V.Ratio", i)]] <- vrs
                
            }
        }
    }
    if (!any(sapply(B[startsWith(names(B), "V.Ratio.")], is.finite))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
    
    #KS Statistics
    if (!(!disp.ks && quick)) {
        if (!(!un && quick)) {
            B[["KS.Un"]] <- col_w_ks(C, treat = treat, weights = NULL, s.weights = s.weights, bin.vars = bin.vars)
        }
        if (!no.adj) {
            for (i in weight.names) {
                B[[paste.("KS", i)]] <- col_w_ks(C, treat = treat, weights = weights[[i]], s.weights = s.weights, bin.vars = bin.vars)
            }
        }
    }
    if (!any(sapply(B[startsWith(names(B), "KS.")], is.finite))) {disp.ks <- FALSE; ks.threshold <- NULL}
    
    
    if (is_not_null(m.threshold)) {
        if (no.adj) {
            B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Diff.Un"]]), paste0(ifelse(abs_(B[["Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), round(m.threshold, 3)), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("M.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Diff", i)]]), paste0(ifelse(abs_(B[[paste.("Diff", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), round(m.threshold, 3)), "")
            }
        }
        
    }
    if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
    
    if (is_not_null(v.threshold)) {
        if (no.adj) {
            B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["V.Ratio.Un"]]), paste0(ifelse(abs_(B[["V.Ratio.Un"]], ratio = TRUE) < v.threshold, "Balanced, <", "Not Balanced, >"), round(v.threshold, 3)), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("V.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("V.Ratio", i)]]), paste0(ifelse(abs_(B[[paste.("V.Ratio", i)]], ratio = TRUE) < v.threshold, "Balanced, <", "Not Balanced, >"), round(v.threshold, 3)), "")
            }
        }
        
    }
    if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == "V.Threshold.Adj"] <- "V.Threshold"
    
    if (is_not_null(ks.threshold)) {
        if (no.adj) {
            B[["KS.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["KS.Un"]]), paste0(ifelse(B[["KS.Un"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), round(ks.threshold, 3)), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("KS.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("KS", i)]]), paste0(ifelse(B[[paste.("KS", i)]] < ks.threshold, "Balanced, <", "Not Balanced, >"), round(ks.threshold, 3)), "")
            }
        }
        
    }
    if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == "KS.Threshold.Adj"] <- "KS.Threshold"
    
    attr(B, "thresholds") <- c(m = m.threshold,
                               v = v.threshold,
                               ks = ks.threshold)
    attr(B, "disp") <- c(means = disp.means,
                         sds = disp.sds,
                         v.ratio = disp.v.ratio,
                         ks = disp.ks)
    
    return(B)
}
balance.table.subclass <- function(C, weights = NULL, treat, subclass, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, disp.means = FALSE, disp.sds = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, s.weights = rep(1, length(treat)), types = NULL, abs = FALSE, quick = TRUE) {
    #Creates list SB of balance tables for each subclass
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    #B=Balance frame
    Bnames <- c("Type", "M.0.Adj", "SD.0.Adj", "M.1.Adj", "SD.1.Adj", "Diff.Adj", "M.Threshold", "V.Ratio.Adj", "V.Threshold", "KS.Adj", "KS.Threshold")
    B <- as.data.frame(matrix(NA_real_, nrow = NCOL(C), ncol = length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- colnames(C)
    #Set var type (binary/continuous)
    if (is_not_null(types)) B[["Type"]] <- types
    else B[["Type"]] <- get.types(C)
    bin.vars <- B[["Type"]] == "Binary"
    tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
    
    SB <- vector("list", nlevels(subclass))
    names(SB) <- levels(subclass)
    
    if (missing(binary) || is_null(binary)) {
        binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
    }
    else binary <- match_arg(binary, c("raw", "std"))
    
    #-------------------------------------
    for (i in levels(subclass)) {
        
        SB[[i]] <- B
        in.subclass <- !is.na(subclass) & subclass==i
        
        #Means for each group
        # if (!(!disp.means && quick)) {
        for (t in c("0", "1")) {
            SB[[i]][[paste.("M", t, "Adj")]] <- col_w_mean(C, subset = treat==tn01[t] & in.subclass)
        }
        # }
        
        #SDs for each group
        # if (!(!disp.sds && quick)) {
        sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else !bin.vars
        for (t in c("0", "1")) {
            sds <- rep(NA_real_, NCOL(C))
            sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == tn01[t] & in.subclass)
            SB[[i]][[paste.("SD", t, "Adj")]] <- sds
        }
        
        # }
        
        #Mean differences
        SB[[i]][["Diff.Adj"]] <- col_w_smd(C, treat = treat, weights = NULL,
                                           std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                           s.d.denom = s.d.denom, abs = FALSE, s.weights = NULL, 
                                           bin.vars = bin.vars, subset = in.subclass)
        
        #Variance ratios
        if (!(!disp.v.ratio && quick)) {
            vrs <- rep(NA_real_, NCOL(C))
            if (any(!bin.vars)) {
                vrs[!bin.vars] <- col_w_vr(C[, !bin.vars, drop = FALSE], treat, weights = NULL, abs = abs, 
                                           s.weights = NULL, bin.vars = bin.vars[!bin.vars],
                                           subset = in.subclass)
            }
            SB[[i]][["V.Ratio.Adj"]] <- vrs
        }
        
        #KS Statistics
        if (!(!disp.ks && quick)) {
            SB[[i]][["KS.Adj"]] <- col_w_ks(C, treat = treat, weights = NULL, s.weights = NULL, bin.vars = bin.vars,
                                            subset = in.subclass)
        }
    }
    
    if (all(sapply(SB, function(x) !any(is.finite(c(x[["SD.0.Adj"]], x[["SD.1.Adj"]])))))) {
        attr(SB, "dont.disp.sds") <- TRUE
        disp.sds <- FALSE
    }
    
    if (is_not_null(m.threshold)) {
        for (i in levels(subclass)) {
            SB[[i]][["M.Threshold"]] <- ifelse(SB[[i]][["Type"]]=="Distance", "", 
                                               paste0(ifelse(is.finite(SB[[i]][["Diff.Adj"]]) & abs_(SB[[i]][["Diff.Adj"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), round(m.threshold, 3)))
        }
    }
    
    if (all(sapply(SB, function(x) !any(is.finite(x[["V.Ratio.Adj"]]))))) {
        attr(SB, "dont.disp.v.ratio") <- TRUE; v.threshold <- NULL
        disp.v.ratio <- FALSE
    }
    if (is_not_null(v.threshold)) {
        for (i in levels(subclass)) {
            SB[[i]][["V.Threshold"]] <- ifelse(SB[[i]][["Type"]]!="Distance" & is.finite(SB[[i]][["V.Ratio.Adj"]]), 
                                               paste0(ifelse(abs_(SB[[i]][["V.Ratio.Adj"]], ratio = TRUE) < v.threshold, "Balanced, <", "Not Balanced, >"), round(v.threshold, 3)), "")
        }
    }
    if (all(sapply(SB, function(x) !any(is.finite(x[["KS.Adj"]]))))) {
        attr(SB, "dont.disp.ks") <- TRUE
        disp.ks <- FALSE
    }
    if (is_not_null(ks.threshold)) {
        for (i in levels(subclass)) {
            SB[[i]][["KS.Threshold"]] <- ifelse(SB[[i]][["Type"]]!="Distance" & is.finite(SB[[i]][["KS.Adj"]]), 
                                                paste0(ifelse(SB[[i]][["KS.Adj"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), round(ks.threshold, 3)), "")
        }
    }
    
    attr(SB, "thresholds") <- c(m = m.threshold,
                                v = v.threshold,
                                ks = ks.threshold)
    attr(SB, "disp") <- c(means = disp.means,
                          sds = disp.sds,
                          v.ratio = disp.v.ratio,
                          ks = disp.ks)
    
    return(SB)
}
balance.table.across.subclass <- function(balance.table, balance.table.subclass.list, subclass.obs, sub.by = NULL, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, s.d.denom = NULL) {
    #Variance ratio, v.threshold, and KS not yet supported
    if (is_not_null(s.d.denom)){
        sub.by <- switch(s.d.denom, treated = "treat",
                         control = "control", "all")
    }
    wsub <- switch(sub.by, treat = "Treated",
                   control = "Control",
                   "Total")
    
    B.A <- balance.table.subclass.list[[1]][c("M.0.Adj", "M.1.Adj", "Diff.Adj")]
    
    for(i in rownames(B.A)) {
        for(j in colnames(B.A)) {
            B.A[i, j] <- sum(vapply(seq_along(balance.table.subclass.list),
                                    function(s) subclass.obs[wsub, s]/sum(subclass.obs[wsub, ]) * (balance.table.subclass.list[[s]][i, j]), numeric(1)))
        }
    }
    B.A.df <- data.frame(balance.table[c("Type", "M.0.Un", "SD.0.Un", "M.1.Un", "SD.1.Un", "Diff.Un", "V.Ratio.Un", "KS.Un")], 
                         B.A, M.Threshold = NA_character_)
    if (is_not_null(m.threshold)) B.A.df[["M.Threshold"]] <- ifelse(B.A.df[["Type"]]=="Distance", "", paste0(ifelse(is.finite(B.A.df[["Diff.Adj"]]) & abs_(B.A.df[["Diff.Adj"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
    return(B.A.df)
}
balance.table.cluster.summary <- function(balance.table.clusters.list, weight.names = NULL, no.adj = FALSE, abs = FALSE, quick = TRUE, types = NULL) {
    
    balance.table.clusters.list <- clear_null(balance.table.clusters.list)
    cont.treat <- "Corr.Un" %in% unique(do.call("c", lapply(balance.table.clusters.list, names)))
    if (no.adj) weight.names <- "Adj"
    
    Brownames <- unique(do.call("c", lapply(balance.table.clusters.list, rownames)))
    #cluster.functions <- c("Min", "Mean", "Median", "Max")
    cluster.functions <- c("Min", "Mean", "Max")
    stats <- if (cont.treat) "Corr" else c("Diff", "V.Ratio", "KS")
    Bcolnames <- c("Type", apply(expand.grid(cluster.functions, stats, c("Un", weight.names)), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
    names(B) <- Bcolnames
    
    if (is_not_null(types)) B[["Type"]] <- types
    else B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(vapply(balance.table.clusters.list, function(y) y[[x, "Type"]], character(1))))), use.names = FALSE)
    
    abs0 <- function(x) {if (abs) abs(x) else (x)}
    funs <- vfuns <- structure(vector("list", length(cluster.functions)), names = cluster.functions)
    for (Fun in cluster.functions) {
        funs[[Fun]] <- function(x, ...) {
            if (!any(is.finite(x))) NA_real_
            else get(tolower(Fun))(x, ...)
        }
        vfuns[[Fun]] <- function(x, ...) {
            if (!any(is.finite(x))) NA_real_
            else if (Fun == "Mean") geom.mean(x, ...)
            else get(tolower(Fun))(x, ...)
        }
        for (sample in c("Un", weight.names)) {
            if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
                if (cont.treat) {
                    B[[paste.(Fun, "Corr", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) abs0(y[[x, paste.("Corr", sample)]])), na.rm = TRUE), numeric(1))
                }
                else {
                    B[[paste.(Fun, "Diff", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) abs0(y[[x, paste.("Diff", sample)]])), na.rm = TRUE), numeric(1))
                    B[[paste.(Fun, "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else vfuns[[Fun]](sapply(balance.table.clusters.list, function(y) y[[x, paste.("V.Ratio", sample)]]), na.rm = TRUE), numeric(1))
                    B[[paste.(Fun, "KS", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) y[[x, paste.("KS", sample)]]), na.rm = TRUE), numeric(1))
                }            
            }
        }
    }
    
    return(B)
}

#base.bal.tab.cont
samplesize.cont <- function(treat, weights = NULL, subclass = NULL, s.weights = NULL, method=c("matching", "weighting", "subclassification"), cluster = NULL, which.cluster = NULL, discarded = NULL) {
    #Computes sample size info. for unadjusted and adjusted samples.
    # method is what method the weights are to be used for. 
    # method="subclassification" is for subclass sample sizes only.
    #method <- match_arg(method)
    if (nlevels(cluster) > 0 && is_not_null(which.cluster)) in.cluster <- cluster == which.cluster
    else in.cluster <- rep(TRUE, length(treat))
    if (is_null(discarded)) discarded <- rep(0, length(treat))
    
    if (length(method) == 1 && method == "subclassification") {
        #stop("Subclassification is not yet surpported with continuous treatments.", call. = FALSE)
        if (is_null(subclass)) stop("subclass must be a vector of subclasses.")
        qbins <- nlevels(subclass)
        
        nn <- as.data.frame(matrix(0, nrow = 1, ncol = qbins))
        
        dimnames(nn) <- list(c("Total"), 
                             paste("Subclass", levels(subclass)))
        
        matched <- !is.na(subclass)
        k <- 0
        for (i in levels(subclass)) {
            qi <- subclass[matched]==i
            qt <- treat[matched][qi]
            if (length(qt)<2){
                if (sum(qt==1)<2)
                    warning("Not enough units in subclass ", i, call. = FALSE)
            }
            k <- k + 1
            nn[, k] <- c(length(qt))
        }
        attr(nn, "tag") <- "Sample sizes by subclass"
    }
    else if (is_null(weights)) {
        nn <- as.data.frame(matrix(0, ncol = 1, nrow = 1))
        if (nunique.gt(s.weights, 2) || !any(s.weights==1) || !all(s.weights %in% c(0,1))) {
            sw <- s.weights[in.cluster]
            
            nn[1, ] <- ESS(sw)
        }
        else {
            nn[1, ] <- sum(in.cluster)
            
        }
        dimnames(nn) <- list(c("All"), 
                             c("Total"))
        attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
    }
    else if (length(weights) == 1) {
        if (method=="matching") {
            
            nn <- as.data.frame(matrix(0, ncol = 1, nrow = 3))
            nn[1, ] <- c(length(treat[in.cluster]))
            nn[2, ] <- c(sum(in.cluster & weights[,1] > 0))
            nn[3, ] <- c(sum(in.cluster & weights[,1] == 0))
            dimnames(nn) <- list(c("All", "Matched", "Unmatched"), 
                                 c("Total"))
            attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
            
            #attr(nn, "tag") <- "Sample sizes"
        }
        else if (method == "weighting") {
            w <- weights[in.cluster, 1]
            sw <- s.weights[in.cluster]
            
            nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
            nn[1, ] <- ESS(sw)
            nn[2, ] <- ESS(w*sw)
            dimnames(nn) <- list(c("Unadjusted", "Adjusted"), 
                                 c("Total"))
            attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
            #attr(nn, "tag") <- "Effective sample sizes"
        }
    }
    else {
        #t <- treat[in.cluster]
        sw <- s.weights[in.cluster]
        nn <- as.data.frame(matrix(0, ncol=1, nrow=1+NCOL(weights)))
        nn[1, ] <- ESS(sw)
        for (i in seq_len(NCOL(weights))) {
            if (method[i] == "matching") {
                nn[1+i,] <- c(sum(in.cluster & weights[,i] > 0))
            }
            else if (method[i] == "weighting") {
                w <- weights[in.cluster, i]
                nn[1+i,] <- ESS(w*sw)
            }
            
        }
        dimnames(nn) <- list(c("Unadjusted", names(weights)), 
                             c("Total"))
        attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
        # if (all(obs$ss.type == "ess")) attr(obs, "tag") <- "Effective sample sizes"
        # else attr(obs, "tag") <- "Sample sizes"
        
    }
    if (length(attr(nn, "ss.type")) > 1 && all(attr(nn, "ss.type")[-1] == "ess")) {
        attr(nn, "tag") <- "Effective sample sizes"
    }
    else attr(nn, "tag") <- "Sample sizes"
    
    return(nn)
}
balance.table.cont <- function(C, weights, treat, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.sds = FALSE, s.weights = rep(1, length(treat)), abs = FALSE, no.adj = FALSE, types = NULL, target.means = NULL, target.sds = NULL, quick = TRUE) {
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    if (no.adj) weight.names <- "Adj"
    else weight.names <- names(weights)
    
    #B=Balance frame
    Bnames <- c("Type", 
                apply(expand.grid(c("M", "SD", "Corr", "R.Threshold"),
                                  c("Un", weight.names)), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow = NCOL(C), ncol = length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- colnames(C)
    
    #Set var type (binary/continuous)
    if (is_not_null(types)) B[,"Type"] <- types
    else B[["Type"]] <- get.types(C)
    bin.vars <- B[["Type"]] == "Binary"
    
    #Means
    if (!((!un || !disp.means) && quick)) {
        B[["M.Un"]] <- col_w_mean(C, weights = NULL, s.weights = s.weights)
    }
    if (!no.adj && !(!disp.means && quick)) {
        for (i in weight.names) {
            B[[paste.("M", i)]] <- col_w_mean(C, weights = weights[[i]], s.weights = s.weights)
        }
    }
    
    #SDs
    sd.computable <- !bin.vars
    if (!((!un || !disp.sds) && quick)) {
        sds <- rep(NA_real_, NCOL(C))
        if (any(sd.computable)) {
            sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE],
                                           weights = NULL, s.weights = s.weights,
                                           bin.vars = bin.vars[sd.computable])
        }
        B[["SD.Un"]] <- sds
    }
    if (!no.adj && !(!disp.sds && quick)) {
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
    if (!any(sapply(B[startsWith(names(B), "SD.")], is.finite))) {disp.sds <- FALSE}
    
    #Correlations
    # if (!(!un && quick)) #Always calculate unadjusted corrs
    B[["Corr.Un"]] <- col_w_corr(C, treat, weights = NULL, abs = abs, s.weights = s.weights, 
                                 bin.vars = bin.vars, na.rm = TRUE)
    if (!no.adj) {
        for (i in weight.names) {
            B[[paste.("Corr", i)]]  <- col_w_corr(C, treat, weights = weights[[i]], abs = abs, s.weights = s.weights, 
                                                  bin.vars = bin.vars, na.rm = TRUE)
        }
    }
    
    if (is_not_null(r.threshold)) {
        if (no.adj) {
            B[["R.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Corr.Un"]]), paste0(ifelse(abs_(B[["Corr.Un"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), round(r.threshold, 3)), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("R.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Corr", i)]]), paste0(ifelse(abs_(B[[paste.("Corr", i)]]) < r.threshold, "Balanced, <", "Not Balanced, >"), round(r.threshold, 3)), "")
            }
        }
    }
    if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == "R.Threshold.Adj"] <- "R.Threshold"
    
    
    attr(B, "thresholds") <- c(r = r.threshold)
    attr(B, "disp") <- c(means = disp.means,
                         sds = disp.sds)
    return(B)
    
}
balance.table.subclass.cont <- function(C, weights = NULL, treat, subclass, r.threshold = NULL, disp.means = FALSE, disp.sds = FALSE, s.weights = rep(1, length(treat)), types = NULL, quick = TRUE) {
    #Creates list SB of balance tables for each subclass
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    #B=Balance frame
    Bnames <- c("Type", "M.Adj", "SD.Adj", "Corr.Adj", "R.Threshold")
    B <- as.data.frame(matrix(nrow=NCOL(C), ncol=length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- colnames(C)
    #Set var type (binary/continuous)
    if (is_not_null(types)) B[["Type"]] <- types
    else B[["Type"]] <- get.types(C)
    
    SB <- vector("list", nlevels(subclass))
    names(SB) <- levels(subclass)
    
    #-------------------------------------
    for (i in levels(subclass)) {
        
        SB[[i]] <- B
        in.subclass <- !is.na(subclass) & subclass==i
        
        if (!(!disp.means && quick)) {
            SB[[i]][["M.Adj"]] <- colMeans(C[in.subclass, , drop = FALSE])
        }
        if (!(!disp.sds && quick)) {
            non.binary <- B[["Type"]] != "Binary"
            sds <- rep(NA_real_, NCOL(C))
            sds[non.binary] <- apply(C[in.subclass, non.binary, drop = FALSE], 2, sd)
            SB[[i]][["SD.Adj"]] <- sds
        }
        
        #Correlations
        # SB[[i]][["Corr.Adj"]] <- apply(C, 2, function(x) w.r(x[in.subclass], y = treat[in.subclass]))
        SB[[i]][["Corr.Adj"]] <- col.w.r(C[in.subclass,], y = treat[in.subclass])
        
    }
    
    if (all(sapply(SB, function(x) !any(is.finite(x[["SD.Adj"]]))))) {
        attr(SB, "dont.disp.sds") <- TRUE
        disp.sds <- FALSE
    }
    
    if (is_not_null(r.threshold)) {
        for (i in levels(subclass)) {
            SB[[i]][["R.Threshold"]] <- ifelse(SB[[i]][["Type"]]=="Distance", "", 
                                               paste0(ifelse(is.finite(SB[[i]][["Corr.Adj"]]) & abs_(SB[[i]][["Corr.Adj"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), round(r.threshold, 3)))
        }
    }
    
    attr(SB, "thresholds") <- c(r = r.threshold)
    attr(SB, "disp") <- c(means = disp.means,
                          sds = disp.sds)
    
    return(SB)
}
balance.table.across.subclass.cont <- function(balance.table, balance.table.subclass.list, subclass.obs, sub.by = NULL, r.threshold = NULL) {
    #Not specified
}

#base.bal.tab.imp
balance.table.imp.summary <- function(bal.tab.imp.list, weight.names = NULL, no.adj = FALSE, abs = FALSE, quick = TRUE, types = NULL) {
    if ("bal.tab" %in% unique(do.call("c", lapply(bal.tab.imp.list, class)))) {
        bal.tab.imp.list <- lapply(bal.tab.imp.list, function(x) x[["Balance"]])}
    cont.treat <- "Corr.Un" %in% unique(do.call("c", lapply(bal.tab.imp.list, names)))
    if (length(weight.names) <= 1) weight.names <- "Adj"
    bal.tab.imp.list <- clear_null(bal.tab.imp.list)
    
    Brownames <- unique(do.call("c", lapply(bal.tab.imp.list, rownames)))
    #imp.functions <- c("Min", "Mean", "Median", "Max")
    imp.functions <- c("Min", "Mean", "Max")
    stats <- if (cont.treat) "Corr" else c("Diff", "V.Ratio", "KS")
    Bcolnames <- c("Type", apply(expand.grid(imp.functions, stats, c("Un", weight.names)), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
    names(B) <- Bcolnames
    
    if (is_not_null(types)) B[["Type"]] <- types
    else B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(sapply(bal.tab.imp.list, function(y) y[[x, "Type"]])))), use.names = FALSE)
    
    abs0 <- function(x) {if (is_null(x)) NA_real_ else if (abs) abs(x) else (x)}
    funs <- vfuns <- structure(vector("list", length(imp.functions)), names = imp.functions)
    for (Fun in imp.functions) {
        funs[[Fun]] <- function(x, ...) {
            if (!any(is.finite(x))) NA_real_
            else get(tolower(Fun))(x, ...)
        }
        vfuns[[Fun]] <- function(x, ...) {
            if (!any(is.finite(x))) NA_real_
            else if (Fun == "Mean") geom.mean(x, ...)
            else get(tolower(Fun))(x, ...)
        }
        for (sample in c("Un", weight.names)) {
            if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
                if (cont.treat) {
                    B[[paste.(Fun, "Corr", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(bal.tab.imp.list, function(y) abs0(y[x, paste.("Corr", sample)])), na.rm = TRUE), numeric(1))
                }
                else {
                    B[[paste.(Fun, "Diff", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(bal.tab.imp.list, function(y) abs0(y[[x, paste.("Diff", sample)]])), na.rm = TRUE), numeric(1))
                    B[[paste.(Fun, "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else vfuns[[Fun]](sapply(bal.tab.imp.list, function(y) y[[x, paste.("V.Ratio", sample)]]), na.rm = TRUE), numeric(1))
                    B[[paste.(Fun, "KS", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(bal.tab.imp.list, function(y) y[[x, paste.("KS", sample)]]), na.rm = TRUE), numeric(1))
                }
            }
        }
    }
    return(B)
}
balance.table.clust.imp.summary <- function(summary.tables, weight.names = NULL, no.adj = FALSE, abs = FALSE, quick = TRUE, types = NULL) {
    #cont.treat <- !is.na(match("bal.tab.cont", unique(do.call("c", lapply(bal.tab.imp.list, class)))))
    #clusters <- unique(do.call("c", lapply(bal.tab.imp.list, function(x) names(x[["Cluster.Balance"]]))))
    #cluster.tables <- lapply(clusters, function(x) lapply(bal.tab.imp.list, function(y) y[["Cluster.Balance"]][[x]]))
    #cluster.balance.across.imps <- lapply(cluster.tables, balance.table.imp.summary, no.adj, quick, types)
    #names(cluster.balance.across.imps) <- clusters
    
    if (!all(vapply(summary.tables, is_null, logical(1)))) {
        Brownames <- unique(do.call("c", lapply(summary.tables, rownames)))
        Bcolnames <- unique(do.call("c", lapply(summary.tables, colnames)))
        cont.treat <- !anyNA(charmatch("Mean.Corr.Un", Bcolnames))
        if (length(weight.names) <= 1) weight.names <- "Adj"
        #imp.functions <- c("Min", "Mean", "Median", "Max")
        imp.functions <- c("Min", "Mean", "Max")
        
        B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)))
        dimnames(B) <- list(Brownames, Bcolnames)
        
        if (is_not_null(types)) B[["Type"]] <- types
        else B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(vapply(summary.tables, function(y) y[[x, "Type"]], character(1))))), use.names = FALSE)
        
        abs0 <- function(x) {if (abs) abs(x) else (x)}
        funs <- vfuns <- structure(vector("list", length(imp.functions)), names = imp.functions)
        for (Fun in imp.functions) {
            funs[[Fun]] <- function(x, ...) {
                if (!any(is.finite(x))) NA_real_
                else get(tolower(Fun))(x, ...)
            }
            vfuns[[Fun]] <- function(x, ...) {
                if (!any(is.finite(x))) NA_real_
                else if (Fun == "Mean") geom.mean(x, ...)
                else get(tolower(Fun))(x, ...)
            }
            for (sample in c("Un", weight.names)) {
                if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
                    if (cont.treat) {
                        B[[paste.(Fun, "Corr", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(summary.tables, function(y) abs0(y[[x, paste.(Fun, "Corr", sample)]])), na.rm = TRUE), numeric(1))
                    }
                    else {
                        B[[paste.(Fun, "Diff", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(summary.tables, function(y) abs0(y[[x, paste.(Fun, "Diff", sample)]])), na.rm = TRUE), numeric(1))
                        B[[paste.(Fun, "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else vfuns[[Fun]](sapply(summary.tables, function(y) y[[x, paste.(Fun, "V.Ratio", sample)]]), na.rm = TRUE), numeric(1))
                        B[[paste.(Fun, "KS", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(summary.tables, function(y) y[[x, paste.(Fun, "KS", sample)]]), na.rm = TRUE), numeric(1))
                    }
                }
            }
        }
    }
    else B <- NULL
    
    return(B)
}
samplesize.across.imps <- function(obs.list) {
    #obs.list <- lapply(bal.tab.imp.list, function(x) x[["Observations"]])
    obs.list <- clear_null(obs.list)
    
    obs <- Reduce("+", obs.list)/length(obs.list)
    attr(obs, "tag") <- paste0("Average ", tolower(attr(obs.list[[1]], "tag")), " across imputations")
    return(obs)
}

#base.bal.tab.multi
balance.table.multi.summary <- function(bal.tab.multi.list, weight.names = NULL, no.adj = FALSE, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, quick = TRUE, types = NULL) {
    if ("bal.tab" %in% unique(do.call("c", lapply(bal.tab.multi.list, class)))) {
        bal.tab.multi.list <- lapply(bal.tab.multi.list, function(x) x[["Balance"]])}
    if (length(weight.names) <= 1) weight.names <- "Adj"
    bal.tab.multi.list <- clear_null(bal.tab.multi.list)
    
    Brownames <- unique(do.call("c", lapply(bal.tab.multi.list, rownames)))
    Bcolnames <- c("Type", expand.grid_string(c("Max.Diff", "M.Threshold", "Max.V.Ratio", "V.Threshold", "Max.KS", "KS.Threshold"), 
                                              c("Un", weight.names), collapse = "."))
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
    names(B) <- Bcolnames
    
    if (is_not_null(types)) B[["Type"]] <- types
    else B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(vapply(bal.tab.multi.list, function(y) y[[x, "Type"]], character(1))))), use.names = FALSE)
    
    for (sample in c("Un", weight.names)) {
        if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
            B[[paste.("Max", "Diff", sample)]] <- vapply(Brownames, function(x) max_(sapply(bal.tab.multi.list, function(y) abs(y[[x, paste.("Diff", sample)]])), na.rm = TRUE), numeric(1))
            B[[paste.("Max", "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else max_(sapply(bal.tab.multi.list, function(y) y[[x, paste.("V.Ratio", sample)]]), na.rm = TRUE), numeric(1))
            B[[paste.("Max", "KS", sample)]] <- vapply(Brownames, function(x) max_(sapply(bal.tab.multi.list, function(y) y[[x, paste.("KS", sample)]]), na.rm = TRUE), numeric(1))
        }
    }
    
    if (is_not_null(m.threshold)) {
        if (no.adj) {
            B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Diff.Un"]]), paste0(ifelse(abs_(B[["Max.Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("M.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.Diff", i)]]), paste0(ifelse(abs_(B[[paste.("Max.Diff", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
    
    if (is_not_null(v.threshold)) {
        if (no.adj) {
            B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.V.Ratio.Un"]]), paste0(ifelse(B[, "Max.V.Ratio.Un"] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("V.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.V.Ratio", i)]]), paste0(ifelse(B[[paste.("Max.V.Ratio", i)]] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "V.Threshold.Adj"] <- "V.Threshold"
    
    if (is_not_null(ks.threshold)) {
        if (no.adj) {
            B[["KS.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.KS.Un"]]), paste0(ifelse(B[["Max.KS.Un"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("KS.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.KS", i)]]), paste0(ifelse(B[[paste.("Max.KS", i)]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "KS.Threshold.Adj"] <- "KS.Threshold"
    
    return(B)
}
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
balance.table.msm.summary <- function(bal.tab.msm.list, weight.names = NULL, no.adj = FALSE, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, quick = TRUE, types = NULL) {
    if ("bal.tab" %in% unique(do.call("c", lapply(bal.tab.msm.list, class)))) {
        bal.tab.msm.list <- lapply(bal.tab.msm.list, function(x) x[["Balance"]])}
    cont.treat <- "Corr.Un" %in% unique(do.call("c", lapply(bal.tab.msm.list, names)))
    if (length(weight.names) <= 1) weight.names <- "Adj"
    
    Brownames <- unique(do.call("c", lapply(bal.tab.msm.list, rownames)))
    Brownames.appear <- vapply(Brownames, function(x) paste(seq_along(bal.tab.msm.list)[sapply(bal.tab.msm.list, function(y) x %in% rownames(y))], collapse = ", "), character(1))
    if (cont.treat) {
        Bcolnames <- c("Type", expand.grid_string(c("Max.Corr", "R.Threshold"), 
                                                  c("Un", weight.names), collapse = "."))
    }
    else {
        Bcolnames <- c("Type", expand.grid_string(c("Max.Diff", "M.Threshold", "Max.V.Ratio", "V.Threshold", "Max.KS", "KS.Threshold"), 
                                                  c("Un", weight.names), collapse = "."))
    }
    
    B <- as.data.frame(matrix(NA, nrow = length(Brownames), ncol = 1 + length(Bcolnames)), row.names = Brownames)
    names(B) <- c("Times", Bcolnames)
    
    if (is_not_null(types)) B[["Type"]] <- types
    else B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, "Type"]] else NA_character_, character(1))))), use.names = FALSE)
    
    B[["Times"]] <- Brownames.appear[Brownames]
    
    for (sample in c("Un", weight.names)) {
        if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
            if (cont.treat) {
                B[[paste.("Max", "Corr", sample)]] <- vapply(Brownames, function(x) max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) abs(y[[x, paste.("Corr", sample)]]) else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
            }
            else {
                B[[paste.("Max", "Diff", sample)]] <- vapply(Brownames, function(x) max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) abs(y[[x, paste.("Diff", sample)]]) else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
                B[[paste.("Max", "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, paste.("V.Ratio", sample)]] else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
                B[[paste.("Max", "KS", sample)]] <- vapply(Brownames, function(x) max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, paste.("KS", sample)]] else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
            }
        }
    }
    
    if (is_not_null(m.threshold)) {
        if (no.adj) {
            B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Diff.Un"]]), paste0(ifelse(abs_(B[["Max.Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("M.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.Diff", i)]]), paste0(ifelse(abs_(B[[paste.("Max.Diff", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
    
    if (is_not_null(v.threshold)) {
        if (no.adj) {
            B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.V.Ratio.Un"]]), paste0(ifelse(B[, "Max.V.Ratio.Un"] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("V.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.V.Ratio", i)]]), paste0(ifelse(B[[paste.("Max.V.Ratio", i)]] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "V.Threshold.Adj"] <- "V.Threshold"
    
    if (is_not_null(ks.threshold)) {
        if (no.adj) {
            B[["KS.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.KS.Un"]]), paste0(ifelse(B[["Max.KS.Un"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("KS.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.KS", i)]]), paste0(ifelse(B[[paste.("Max.KS", i)]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "KS.Threshold.Adj"] <- "KS.Threshold"
    
    if (is_not_null(r.threshold)) {
        if (no.adj) {
            B[["R.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Corr.Un"]]), paste0(ifelse(B[["Max.Corr.Un"]] < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("R.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.Corr", i)]]), paste0(ifelse(B[[paste.("Max.Corr", i)]] < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "R.Threshold.Adj"] <- "R.Threshold"
    
    
    return(B)
}
samplesize.msm <- function(bal.tab.msm.list) {
    obs <- do.call("cbind", lapply(bal.tab.msm.list, function(x) x[["Observations"]]))
    attr(obs, "tag") <- attr(bal.tab.msm.list[[1]][["Observations"]], "tag")
    attr(obs, "ss.type") <- attr(bal.tab.msm.list[[1]][["Observations"]], "ss.type")
    return(obs)
}

#base.bal.tab.target
balance.table.target.summary <- balance.table.multi.summary
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

#set.cobalt.options
acceptable.options <- function() {
    TF <- c(TRUE, FALSE)
    return(list(un = TF,
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