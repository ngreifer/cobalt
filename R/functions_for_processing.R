#bal.tab
process_obj <- function(obj) {
    if (is.null(obj)) obj <- list()
    else {
        if (isS4(obj)) obj <- asS3(obj)
        
        #npCBPS
        if (inherits(obj, "npCBPS")) {
            class(obj) <- c("CBPS", "npCBPS")
        }
        #ebalance.trim
        else if (inherits(obj, "ebalance.trim")) {
            class(obj) <- "ebalance"
        }
        #Matchby
        else if (inherits(obj, "Matchby")) {
            class(obj) <- "Match"
        }
        #cem.match.list
        else if (inherits(obj, "cem.match.list")) {
            class(obj) <- c("cem.match", "cem.match.list")
        }
        #time.list
        else if (is.list(obj) && !is.data.frame(obj) &&
                 all(vapply(obj, rlang::is_formula, logical(1L)))) {
            class(obj) <- c("formula.list", "time.list", class(obj))
        }
        else if (is.list(obj) && !is.data.frame(obj) &&
                 all(vapply(obj, is.data.frame, logical(1L)))) {
            class(obj) <- c("data.frame.list", "time.list", class(obj))
        }
        #designmatch
        else if (is.list(obj) && length(obj) >= 6) {
            dm.b.names <- c("obj_total", "obj_dist_mat", "t_id", 
                            "c_id", "group_id", "time")
            dm.n.names <- c("obj_total", "obj_dist_mat", "id_1", 
                            "id_2", "group_id", "time")
            if (all(dm.b.names %in% names(obj)) || all(dm.n.names %in% names(obj))) {
                class(obj) <- c("designmatch")
            }
        }
    }
    class(obj) <- c(class(obj), "cobalt.processed.obj")
    
    obj
}

#x2base
process_treat <- function(treat, datalist = list(), keep_values = FALSE) {
    
    .chk_not_missing(treat, "`treat`")
    
    if (inherits(treat, "unprocessed.treat")) {
        attrs <- attributes(treat)
        renamed_original <- setNames(names(treat_vals(treat)), treat_vals(treat))
        treat <- factor(renamed_original[as.character(treat)], levels = renamed_original)
        for (at in c("treat_names", "treat_vals", "keep_values", "treat.type", "names"))
            attr(treat, at) <- attrs[[at]]
    }
    else {
        # keep_values <- isTRUE(attr(treat, "keep_values")) || 
        #     (has.treat.type(treat) && get.treat.type(treat) == "multinomial")
        
        treat <- .process_vector(treat, name = "treat", 
                                 which = "treatment statuses", 
                                 datalist = datalist, missing.okay = FALSE)
        
        treat <- assign.treat.type(treat)
        treat.type <- get.treat.type(treat)
        
        if (treat.type == "binary") {
            if (!is.factor(treat)) treat <- factor(treat, levels = sort(unique(treat, nmax = 2)))
            original_values <- levels(treat)[levels(treat) %in% unique(treat, nmax = 2)]
            treat_names(treat) <- {
                if (!keep_values && can_str2num(as.character(treat)) && all(original_values %in% c("0", "1"))) {
                    setNames(c("Control", "Treated"), c("control", "treated"))
                }
                else {
                    setNames(original_values, c("control", "treated"))
                }
            }
            
            treat_vals(treat) <- setNames(original_values, treat_names(treat))
        }
        else if (treat.type == "multinomial") {
            treat <- factor(treat, ordered = FALSE)
            treat_names(treat) <- setNames(levels(treat), levels(treat))
            treat_vals(treat) <- setNames(levels(treat), treat_names(treat))
        }
        attr(treat, "treat.type") <- treat.type
        # attr(treat, "keep_values") <- keep_values
    }
    class(treat) <- c("processed.treat", class(treat))
    
    treat
}
unprocess_treat <- function(treat) {
    if (!inherits(treat, "processed.treat")) return(treat)
    
    attrs <- attributes(treat)
    treat <- treat_vals(treat)[as.character(treat)]
    attributes(treat) <- attrs[setdiff(names(attrs), "class")]
    class(treat) <- c("unprocessed.treat", class(treat_vals(treat)))
    
    treat
}
process_treat.list <- function(treat.list, datalist = list()) {
    .chk_not_missing(treat.list, "`treat.list`")
    
    if (!is.list(treat.list)) {
        treat.list <- as.list(treat.list)
    }
    
    treat.list.names <- vapply(seq_along(treat.list), function(ti) {
        if (is.character(treat.list[[ti]]) && length(treat.list[[ti]])==1L && is_not_null(datalist)) {
            treat.list[[ti]]
        }
        else if (rlang::is_named(treat.list)) names(treat.list)[ti]
        else as.character(ti)
    }, character(1L))
    treat.list <- lapply(treat.list, process_treat, datalist = datalist)
    
    setNames(treat.list, treat.list.names)
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
subset_processed.treat <- function(x, index) {
    y <- x[index]
    treat_names(y) <- treat_names(x)[treat_vals(x) %in% unique(y)]
    treat_vals(y) <- treat_vals(x)[treat_vals(x) %in% unique(y)]
    y <- assign.treat.type(y)
    class(y) <- class(x)
    y
}
# `[.processed.treat` <- function(x, ...) {
#     y <- NextMethod("[")
#     # if (is.factor(y)) y <- droplevels(y)
#     treat_names(y) <- treat_names(x)[treat_vals(x) %in% unique(y)]
#     treat_vals(y) <- treat_vals(x)[treat_vals(x) %in% unique(y)]
#     y <- assign.treat.type(y)
#     class(y) <- class(x)
#     y
# }
# `[.unprocessed.treat` <- `[.processed.treat`

initialize_X <- function() {
    X.names <- c("covs",
                 "treat",
                 "weights",
                 "distance",
                 "addl",
                 "s.d.denom",
                 "estimand",
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
    make_list(X.names)
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
    make_list(X.names)
}
.weight_check <- function(w) {
    wname <- deparse1(substitute(w))
    if (!is.list(w)) w <- list(w)
    if (anyNA(w, recursive = TRUE)) .err(sprintf("`NA`s are not allowed in the %s", wname))
    for (x in w) {if (!all(is.numeric(x))) .err(sprintf("all %s must be numeric", wname))}
    for (x in w) {if (!all(is.finite(x))) .err(sprintf("infinite %s are not allowed", wname))}
    # for (x in w) {if (any(x < 0)) {
    #   .wrn(paste0("Negative ", wname, " found. This may yield nonsensical results and errors as the square root of negative weights is encountered."))
    #   break
    # }}
}
.cluster_check <- function(cluster, treat) {
    if (!is.list(treat)) treat <- list(treat)
    
    stop_warn <- c(cw = FALSE, bw = FALSE, bs = FALSE)
    for (t in treat) {
        if (!has.treat.type(t)) t <- assign.treat.type(t)
        treat.type <- get.treat.type(t)
        
        if (treat.type == "continuous" && !stop_warn["cw"]) {
            tabu <- tabulate(cluster, nbins = nlevels(cluster))
            if (any(tabu == 1)) stop_warn["cw"] <- TRUE
        }
        else if (treat.type != "continuous" && !all(stop_warn[c("bw", "bs")])) {
            tab <- table(cluster, t)
            if (any(tab == 0)) stop_warn["bs"] <- TRUE
            else if (any(tab == 1)) stop_warn["bw"] <- TRUE
        }
    }
    if (stop_warn["cw"]) .wrn("Some clusters have only one unit in them, which may yield unxpected results")
    if (stop_warn["bw"]) .wrn("Some clusters have only one member of a treatment group in them, which may yield unxpected results")
    if (stop_warn["bs"]) .err("not all treatment levels are present in all clusters")
    
}
strata2weights <- function(strata, treat, estimand = NULL, focal = NULL) {
    #Process strata into weights (similar to weight.subclass from MatchIt)
    
    #Checks
    if (!is.atomic(strata) || is_not_null(dim(strata))) {
        .err("`strata` must be an atomic vector or factor")
    }
    
    #Process treat
    treat <- process_treat(treat)
    
    if (get.treat.type(treat) == "continuous") {
        .err("`strata` cannot be turned into weights for continuous treatments")
    }
    
    s.d.denom <- .get_s.d.denom(NULL, estimand = estimand, subclass = strata, treat = treat, focal = focal, quietly = TRUE)
    if (s.d.denom %in% treat_vals(treat)) focal <- process_focal(s.d.denom, treat)
    else focal <- NULL
    
    NAsub <- is.na(strata)
    imat <- do.call("cbind", lapply(treat_vals(treat), function(t) treat == t & !NAsub))
    colnames(imat) <- treat_vals(treat)
    
    weights <- rep(0, length(treat))
    
    if (!is.factor(strata)) {
        strata <- factor(strata, nmax = min(colSums(imat)))
        levels(strata) <- seq_len(nlevels(strata))
    }
    
    t_by_sub <- do.call("rbind", lapply(treat_vals(treat), function(t) tabulate(strata[imat[,t]], nlevels(strata))))
    dimnames(t_by_sub) <- list(treat_vals(treat), levels(strata))
    
    total_by_sub <- colSums(t_by_sub)
    
    strata.c <- as.character(strata)
    
    if (is_not_null(focal)) {
        focal <- process_focal(focal, treat)
        for (t in treat_vals(treat)) {
            weights[imat[,t]] <- {
                if (t == focal)  1
                else (t_by_sub[focal,]/t_by_sub[t,])[strata.c[imat[,t]]]
            }
        }
    }
    else {
        for (t in treat_vals(treat)) {
            weights[imat[,t]] <- (total_by_sub/t_by_sub[t,])[strata.c[imat[,t]]]
        }
    }
    
    if (any(na.w <- !is.finite(weights))) {
        weights[na.w] <- 0
        .wrn("Some units were given weights of zero due to zeros in stratum membership")
    }
    
    if (all(check_if_zero(weights))) 
        .err("no units were stratified")
    
    for (tnn in names(treat_names(treat))) {
        if (all(check_if_zero(weights[treat == treat_vals(treat)[treat_names(treat)[tnn]]])))
            .err(sprintf("No %s units were stratified", tnn))
    }
    
    attr(weights, "match.strata") <- strata
    weights
}

.use_tc_fd <- function(formula = NULL, data = NULL, treat = NULL, covs = NULL, needs.treat = TRUE, needs.covs = TRUE) {
    
    treat_f <- treat_c <- covs_f <- covs_c <- NULL
    
    if (is_not_null(formula) && 
        (rlang::is_formula(formula <- try(as.formula(formula), silent = TRUE)))) {
        D <- NULL
        if (is_not_null(data)) D <- data
        if (is_mat_like(covs)) {
            D <- {
                if (is_not_null(D)) cbind(D, as.data.frame(covs)) 
                else as.data.frame(covs)
            }
        }
        treat_f <- try(get_treat_from_formula(formula, D, treat = treat), silent = TRUE)
        covs_f <- try(get_covs_from_formula(formula, D), silent = TRUE)
    }
    
    if (is_not_null(covs)) {
        covs_c <- {
            if (is_mat_like(covs)) try(get_covs_from_formula(data = covs), silent = TRUE)
            else if (is.character(covs)) {
                if (!is_mat_like(data)) {
                    .err("if `covs` is a character vector, `data` must be specified as a data.frame")
                }
                
                if (!all(covs %in% colnames(data))) {
                    .err("all entries in `covs` must be names of variables in `data`")
                }
                try(get_covs_from_formula(f.build(covs), data = as.data.frame(data)), silent = TRUE)
                
            }
            else .err("`covs` must be a data.frame of covariates")
        }
    }
    
    if (is_not_null(treat)) {
        treat_c <- try({
            if (!is.atomic(treat)) .err("`treat` must be an vector of treatment statuses")
            if (is.character(treat) && length(treat) == 1) get_treat_from_formula(reformulate(".", treat), data = data)
            else get_treat_from_formula(treat ~ ., treat = treat)
        }, silent = TRUE)
    }
    
    covs_to_use <- treat_to_use <- "c"
    if (is_error(covs_c)) {
        if (is_error(covs_f)) {
            .err(attr(covs_c, "condition")$message, tidy = FALSE)
        }
        covs_to_use <- "f"
    }
    else if (is_null(covs_c)) {
        if (is_error(covs_f)) {
            .err(attr(covs_f, "condition")$message, tidy = FALSE)
        }
        if (is_null(covs_f) && needs.covs) {
            .err("no covariates were specified")
        }
        covs_to_use <- "f"
    }
    
    if (is_error(treat_c)) {
        if (is_error(treat_f)) {
            .err(attr(treat_c, "condition")$message, tidy = FALSE)
        }
        treat_to_use <- "f"
    }
    else if (is_null(treat_c)) {
        if (is_error(treat_f)) {
            .err(attr(treat_f, "condition")$message, tidy = FALSE)
        }
        if (is_null(treat_f) && needs.treat) {
            .err("no treatment variable was specified")
        }
        treat_to_use <- "f"
    }
    
    t.c <- list(treat = switch(treat_to_use, "f" = treat_f, "c" = treat_c), 
                covs = switch(covs_to_use, "f" = covs_f, "c" = covs_c))
    t.c$treat.name <- attr(t.c$treat, "treat.name")
    
    t.c
}
.process_val <- function(val, i, treat = NULL, covs = NULL, addl.data = list(), ...) {
    
    if (is.data.frame(val)) {
        return(val)
    }
    
    if (is.numeric(val)) {
        return(setNames(data.frame(val), i))
    }
    
    if (!is.character(val)) {
        if (i == "weights") {
            .err("the argument supplied to `weights` must be a named list of weights, names of variables containing weights in an available data set, or objects with a `get.w()` method")
        }
        else {
            .err(sprintf("the argument supplied to %s must be a vector, a data.frame, or the names of variables in an available data set", add_quotes(i, "`")))
        }
    }
    
    addl.data <- addl.data[lengths(addl.data) > 0]
    if ((is_null(covs) && is_null(treat)) ||
        length(val) == NROW(covs) ||
        length(val) == length(treat) ||
        (is_not_null(addl.data) && length(val) > max(vapply(addl.data, ncol, numeric(1))))) {
        return(setNames(data.frame(val), i))
    }
    
    if (is_null(addl.data)) {
        .wrn(sprintf("names were provided to %s, but no argument to `data` was provided. Ignoring %s",
                     add_quotes(i, "`"), add_quotes(i, "`")))
        return(NULL)
    }
    
    val <- unique(val)
    val.df <- make_df(val, nrow = max(vapply(addl.data, nrow, numeric(1))))
    not.found <- rlang::rep_named(val, FALSE)
    for (v in val) {
        found <- FALSE
        k <- 1
        while (!found && k <= length(addl.data)) {
            if (v %in% names(addl.data[[k]])) {
                val.df[[v]] <- addl.data[[k]][[v]]
                found <- TRUE
            }
            else k <- k + 1
        }
        if (!found) not.found[v] <- TRUE
    }
    
    if (any(not.found)) {
        .wrn(sprintf("the following variable(s) named in %s are not in any available data sets and will be ignored: %s",
                     add_quotes(i, "`"), paste(val[not.found])))
        val.df <- val.df[!not.found]
    }
    
    val.df
}
.process_data_frame <- function(i, df, treat = NULL, covs = NULL, addl.data = list(), ...) {
    if (is_null(df)) return(NULL)
    
    val <- df
    val.df <- NULL
    
    if (i == "weights" && any(has_method(class(val), "get.w"))) {
        val <- list(val)
    }
    
    if (!is.list(val) || is.data.frame(val)) {
        return(.process_val(val, i, treat, covs, addl.data = addl.data))
    }
   
    if (i == "weights") {
        #Use get.w() on inputs
        for (x in seq_along(val)) {
            val[[x]] <- process_obj(val[[x]])
            if (any(has_method(class(val[[x]]), "get.w"))) {
                get.w_class <- class(val[[x]])[has_method(class(val[[x]]), "get.w")][1]
                val[[x]] <- get.w(val[[x]], treat = treat, ...)
                if (rlang::names2(val)[x] == "") names(val)[x] <- get.w_class
            }
        }
    }
    
    val.list <- lapply(val, function(x) .process_val(x, i, treat, covs, addl.data = addl.data))
    
    if (!rlang::is_named(val.list)) {
        .err(sprintf("all entries in `%s` must have names", i))
    }
    
    for (x in seq_along(val.list)) {
        if (NCOL(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
    }
    
    if (!all_the_same(vapply(val.list, nrow, numeric(1)))) {
        .err(sprintf("not all items in `%s` have the same length", i))
    }
    
    setNames(do.call("cbind", val.list),
             unlist(lapply(val.list, names)))

}
.process_list <- function(i, List, ntimes, call.phrase, treat.list = list(), covs.list = list(), addl.data = list(), ...) {
    if (is_null(List)) return(NULL)
    
    val.List <- List
    
    if (!is.list(val.List)) {
        val.List <- list(val.List)
    }
    
    if (length(val.List) == 1) {
        val.List <- replicate(ntimes, val.List)
    }
    else if (length(val.List) != ntimes) {
        .err(sprintf("the argument to %s must be a list of the same length as the number of time points in %s", add_quotes(i, "`"), call.phrase))
    }
    
    for (ti in seq_along(val.List)) {
        val <- val.List[[ti]]
        if (is_null(val)) next
        
        val.df <- NULL
        
        if (is.list(val) && !is.data.frame(val)) {
            val.list <- lapply(val, function(x) .process_val(x, strsplit(i, ".list", fixed = TRUE)[[1]], treat.list[[ti]], covs.list[[ti]], addl.data = addl.data))
            
            if (!all_the_same(vapply(val.list, nrow, numeric(1)))) {
                .err(sprintf("not all items in `%s` have the same length", i))
            }
            
            for (x in val.list) {
                if (NCOL(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
            }
            
            val.df <- setNames(do.call("cbind", val.list),
                               vapply(val.list, names, character(1)))
        }
        else {
            val.df <- .process_val(val, strsplit(i, ".list", fixed = TRUE)[[1]], treat.list[[ti]], covs.list[[ti]], addl.data = addl.data)
        }
        
        if (is_not_null(val.df) && anyNA(val.df)) {
            .err(sprintf("missing values exist in %s", add_quotes(i, "`")))
        }
        
        val.List[[ti]] <- val.df
    }
    
    val.df.lengths <- vapply(val.List[lengths(val.List) > 0], nrow, numeric(1))
    
    if (max(val.df.lengths) != min(val.df.lengths)) {
        .err(sprintf("all columns in `%s` need to have the same number of rows", i))
    }
    
    val.List
}
.process_vector <- function(vec, name = deparse1(substitute(vec)), which = name, datalist = list(), missing.okay = FALSE) {
    bad.vec <- FALSE
    if (is.character(vec) && length(vec) == 1L && is_not_null(datalist)) {
        for (i in seq_along(datalist)) {
            if (is.matrix(datalist[[i]]) && vec %in% colnames(datalist[[i]])) {
                vec <- datalist[[i]][,vec]
                break
            }
            else if (is.data.frame(datalist[[i]]) && vec %in% names(datalist[[i]])) {
                vec <- datalist[[i]][[vec]]
                break
            }
            else if (i == length(datalist)) bad.vec <- TRUE
        }
    }
    else if (is.atomic(vec) && length(vec) > 1L) {
        vec <- vec
    }
    else {
        bad.vec <- TRUE
    }
    
    if (bad.vec) {
        .err(sprintf("The argument to %s must be a vector of %s or the (quoted) name of a variable in `data` that contains %s",
                     add_quotes(name, "`"), which, which))
    }
    
    if (!missing.okay && anyNA(vec)) {
        .err(sprintf("Missing values exist in %s", add_quotes(name, "`")))
    }
    
    vec
} 
.get_s.d.denom <- function(s.d.denom = NULL, estimand = NULL, weights = NULL,
                           subclass = NULL, treat = NULL, focal = NULL, quietly = FALSE) {
    check.estimand <- check.weights <- check.focal <- bad.s.d.denom <- bad.estimand <- FALSE
    s.d.denom.specified <- !missing(s.d.denom) && is_not_null(s.d.denom)
    estimand.specified <- is_not_null(estimand)
    
    if (s.d.denom.specified) {
        
        treat <- process_treat(treat)
        
        unique.treats <- as.character(treat_vals(treat))
        allowable.s.d.denoms <- c("pooled", "all", "weighted", "hedges")
        if (length(treat_names(treat)) == 2 && all(c("treated", "control") %in% names(treat_names(treat))))
            allowable.s.d.denoms <- c(allowable.s.d.denoms, "treated", "control")
        if (is_not_null(focal)) allowable.s.d.denoms <- c(allowable.s.d.denoms, "focal")
        
        if (length(s.d.denom) > 1 && length(s.d.denom) != NCOL(weights)) {
            .err(sprintf("`s.d.denom` must have length 1 or equal to the number of valid sets of weights, which is %s",
                         NCOL(weights)))
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
            if (!inherits(treat, "processed.treat")) treat <- process_treat(treat)
            try.estimand <- tryCatch(match_arg(toupper(estimand), c("ATT", "ATC", "ATE", "ATO", "ATM"), several.ok = TRUE),
                                     error = function(cond) NA_character_)
            if (anyNA(try.estimand) || any(try.estimand %in% c("ATC", "ATT")) && get.treat.type(treat) != "binary") {
                check.focal <- TRUE
            }
            else {
                if (length(try.estimand) > 1 && length(try.estimand) != NCOL(weights)) {
                    .err(sprintf("`estimand` must have length 1 or equal to the number of valid sets of weights, which is %s",
                                 NCOL(weights)))
                }
                
                s.d.denom <- vapply(try.estimand, switch, character(1L),
                                         ATT = treat_vals(treat)[treat_names(treat)["treated"]], 
                                         ATC = treat_vals(treat)[treat_names(treat)["control"]], 
                                         ATO = "weighted",
                                         ATM = "weighted",
                                         "pooled")
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
        if (!inherits(treat, "processed.treat")) treat <- process_treat(treat)
        if (is_null(weights) && is_null(subclass)) {
            s.d.denom <- "pooled"
        }
        else if (is_not_null(subclass)) {
            sub.tab <- table(treat, subclass)[treat_vals(treat), ]
            sub.tab <- rbind(sub.tab, table(subclass)[colnames(sub.tab)])
            dimnames(sub.tab) <- list(c(treat_vals(treat), "pooled"), colnames(sub.tab))
            
            ranges <- apply(sub.tab, 1, function(x) .mean_abs_dev(x)/sum(x))
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
                NA_character_
            }, character(1L))
            s.d.denom[is.na(s.d.denom)] <- "pooled"
        }
    }
    if (is_not_null(weights) && length(s.d.denom) == 1) s.d.denom <- rep.int(s.d.denom, NCOL(weights))
    
    if (s.d.denom.specified && bad.s.d.denom && (!estimand.specified || bad.estimand)) {
        attr(s.d.denom, "note") <- sprintf("warning: `s.d.denom` should be one of %s.\n         Using %s instead",
                                           word_list(unique(c(unique.treats, allowable.s.d.denoms)), "or", quotes = 2),
                                           word_list(s.d.denom, quotes = 2))
    }
    else if (estimand.specified && bad.estimand) {
        attr(s.d.denom, "note") <- sprintf("warning: `estimand` should be one of %s. Ignoring `estimand`",
                                           word_list(c("ATT", "ATC", "ATE"), "or", quotes = 2))
    }
    else if ((check.focal || check.weights) && !all(s.d.denom %in% treat_vals(treat))) {
        attr(s.d.denom, "note") <- sprintf("note: `s.d.denom` not specified; assuming %s", 
                                           if (all_the_same(s.d.denom)) s.d.denom[1] 
                                           else word_list(paste0(add_quotes(vapply(s.d.denom, function(s) {
                                               if (s %in% treat_vals(treat) && all(treat_vals(treat) %in% c("0", "1"))) {
                                                   names(treat_names(treat))[treat_names(treat) == names(treat_vals(treat))[treat_vals(treat) == s]]
                                               }
                                               else s
                                           }, character(1L))), " for ", names(weights))))
    }
    
    if (is_not_null(weights) && length(s.d.denom) != NCOL(weights)) {
        .err(sprintf("Valid inputs to `s.d.denom` or `estimand` must have length 1 or equal to the number of valid sets of weights, which is %s",
                     NCOL(weights)))
    }
    
    if (!quietly && is_not_null(attr(s.d.denom, "note"))) .msg(attr(s.d.denom, "note"))
    
    s.d.denom
}
.get_s.d.denom.cont <- function(s.d.denom, weights = NULL, subclass = NULL, quietly = FALSE) {
    bad.s.d.denom <- FALSE
    s.d.denom.specified <- !missing(s.d.denom) && is_not_null(s.d.denom)
    
    if (is_not_null(subclass)) {
        s.d.denom <- "all"
    }
    else if (s.d.denom.specified) {
        allowable.s.d.denoms <- c("all", "weighted")
        
        if (length(s.d.denom) > 1 && length(s.d.denom) != NCOL(weights)) {
            .err("`s.d.denom` must have length 1 or equal to the number of valid sets of weights")
        }
        
        s.d.denom <- match_arg(s.d.denom, unique(allowable.s.d.denoms), 
                               several.ok = TRUE)
    }
    else {
        s.d.denom <- "all"
    }
    
    if (is_not_null(weights) && NCOL(weights) > 1 && length(s.d.denom) == 1) s.d.denom <- rep.int(s.d.denom, NCOL(weights))
    
    if (!quietly) {
        if (s.d.denom.specified && bad.s.d.denom) {
            .msg(sprintf("warning: `s.d.denom` should be %s.\n         Using %s instead",
                         word_list(unique(allowable.s.d.denoms), "or", quotes = 2),
                         word_list(s.d.denom, quotes = 2)))
        }
    }
    
    if (is_not_null(weights) && length(s.d.denom) != NCOL(weights)) {
        .err("valid inputs to `s.d.denom` or `estimand` must have length 1 or equal to the number of valid sets of weights")
    }
    
    if (is_not_null(weights)) names(s.d.denom) <- names(weights)
    
    s.d.denom
}
.compute_s.d.denom <- function(mat, treat, s.d.denom = "pooled", s.weights = NULL,
                               bin.vars = NULL, subset = NULL, weighted.weights = NULL,
                               to.sd = rep(TRUE, ncol(mat)), na.rm = TRUE) {
    denoms <- setNames(rep(1, ncol(mat)), colnames(mat))
    if (is.character(s.d.denom) && length(s.d.denom) == 1L) {
        if (is_null(bin.vars)) {
            bin.vars <- rep(FALSE, ncol(mat))
            bin.vars[to.sd] <- is_binary_col(mat[subset, to.sd,drop = FALSE])
        }
        else if (!is.atomic(bin.vars) || length(bin.vars) != ncol(mat) ||
                 anyNA(as.logical(bin.vars))) {
            .err("`bin.vars` must be a logical vector with length equal to the number of columns of `mat`")
        }
        
        possibly.supplied <- c("mat", "treat", "weighted.weights", "s.weights", "subset")
        lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                            possibly.supplied)
        supplied <- lengths > 0
        if (!all_the_same(lengths[supplied])) {
            .err(sprintf("%s must have the same number of units",
                         word_list(possibly.supplied[supplied], quotes = "`")))
        }
        
        if (lengths["weighted.weights"] == 0) weighted.weights <- rep(1, NROW(mat))
        if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
        
        if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
        else if (anyNA(as.logical(subset))) .err("`subset` must be a logical vector")
        
        if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
        cont.treat <- get.treat.type(treat) == "continuous"
        
        if (cont.treat) {
            unique.treats <- NULL
            s.d.denom <- .get_s.d.denom.cont(as.character(s.d.denom), weights = weighted.weights[subset])
        }
        else {
            unique.treats <- {
                if (inherits(treat, "processed.treat") && all(subset)) as.character(treat_vals(treat))
                else as.character(unique(treat[subset]))
            }
            s.d.denom <- .get_s.d.denom(as.character(s.d.denom), weights = weighted.weights[subset], treat = treat[subset])
            if (s.d.denom %in% c("treated", "control") && s.d.denom %nin% unique.treats) {
                s.d.denom <- treat_vals(treat)[treat_names(treat)[s.d.denom]]
            }
            treat <- as.character(treat)
        }
        
        if (s.d.denom %in% unique.treats)
            denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                                  unique.treats, na.rm) {
                sqrt(col.w.v(mat[treat == s.d.denom, , drop = FALSE],
                             w = s.weights[treat == s.d.denom],
                             bin.vars = bin.vars, na.rm = na.rm))
            }
        else if (s.d.denom == "pooled")
            denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                                  unique.treats, na.rm) {
                sqrt(Reduce("+", lapply(unique.treats,
                                        function(t) col.w.v(mat[treat == t, , drop = FALSE],
                                                            w = s.weights[treat == t],
                                                            bin.vars = bin.vars, na.rm = na.rm))) / length(unique.treats))
            }
        else if (s.d.denom == "all")
            denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                                  unique.treats, na.rm) {
                sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm))
            }
        else if (s.d.denom == "weighted")
            denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                                  unique.treats, na.rm) {
                sqrt(col.w.v(mat, w = weighted.weights * s.weights, bin.vars = bin.vars, na.rm = na.rm))
            }
        else if (s.d.denom == "hedges")
            denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                                  unique.treats, na.rm) {
                (1 - 3/(4*length(treat) - 9))^-1 *
                    sqrt(Reduce("+", lapply(unique.treats,
                                            function(t) (sum(treat == t) - 1) * col.w.v(mat[treat == t, , drop = FALSE],
                                                                                        w = s.weights[treat == t],
                                                                                        bin.vars = bin.vars, na.rm = na.rm))) / (length(treat) - 2))
            }
        else .err("`s.d.denom` is not an allowed value")
        
        denoms[to.sd] <- denom.fun(mat = mat[, to.sd, drop = FALSE], treat = treat, s.weights = s.weights,
                                   weighted.weights = weighted.weights, bin.vars = bin.vars[to.sd],
                                   unique.treats = unique.treats, na.rm = na.rm)
        
        if (any(zero_sds <- !is.finite(denoms[to.sd]) | check_if_zero(denoms[to.sd]))) {
            denoms[to.sd][zero_sds] <- sqrt(col.w.v(mat[, to.sd, drop = FALSE][, zero_sds, drop = FALSE],
                                                    w = s.weights,
                                                    bin.vars = bin.vars[to.sd][zero_sds], na.rm = na.rm))
        }
        
        if (cont.treat) {
            treat.sd <- denom.fun(mat = treat, s.weights = s.weights,
                                  weighted.weights = weighted.weights, bin.vars = FALSE,
                                  na.rm = na.rm)
            denoms[to.sd] <- denoms[to.sd]*treat.sd
        }
    }
    else if (is.numeric(s.d.denom)) {
        if (is_not_null(names(s.d.denom)) && any(colnames(mat) %in% names(s.d.denom))) {
            denoms[colnames(mat)[colnames(mat) %in% names(s.d.denom)]] <- s.d.denom[names(s.d.denom)[names(s.d.denom) %in% colnames(mat)]]
        }
        else if (length(s.d.denom) == sum(to.sd)) {
            denoms[to.sd] <- s.d.denom
        }
        else if (length(s.d.denom) == ncol(mat)) {
            denoms[] <- s.d.denom
        }
        else {
            .err("`s.d.denom` must be an allowable value or a numeric vector of with length equal to the number of columns of `mat`. See ?col_w_smd for allowable values")
        }
    }
    else {
        .err("`s.d.denom` must be an allowable value or a numeric vector of with length equal to the number of columns of `mat`. See ?col_w_smd for allowable values")
    }
    
    denoms
}
.assign_X_class <- function(X) {
    X <- clear_null(X)
    
    if (is_not_null(X[["treat"]]) && !has.treat.type(X[["treat"]])) X[["treat"]] <- assign.treat.type(X[["treat"]])
    
    if (is_not_null(X[["subclass"]])) {
        if (get.treat.type(X[["treat"]]) == "binary") X.class <- "subclass.binary"
        else if (get.treat.type(X[["treat"]]) == "continuous") X.class <- "subclass.cont"
        else .err("Multi-category treatments are not currently compatible with subclasses")
    }
    else if (is_not_null(X[["cluster"]]) && nlevels(X[["cluster"]]) > 1) X.class <- "cluster"
    else if (is_not_null(X[["covs.list"]])) X.class <- "msm"
    else if (get.treat.type(X[["treat"]]) == "multinomial") X.class <- "multi"
    else if (is_not_null(X[["imp"]]) && nlevels(X[["imp"]]) > 1) X.class <- "imp"
    else if (get.treat.type(X[["treat"]]) == "binary") X.class <- "binary"
    else if (get.treat.type(X[["treat"]]) == "continuous") X.class <- "cont"
    else probably.a.bug()
    
    attr(X, "X.class") <- X.class
    
    X
}
get_length_X <- function(X) {
    if (is_not_null(X[["treat"]])) length(X[["treat"]])
    else if (is_not_null(X[["covs"]])) nrow(X[["covs"]])
    else if (is_not_null(X[["treat.list"]])) length(X[["treat.list"]][[1]])
    else if (is_not_null(X[["covs.list"]])) nrow(X[["covs.list"]][[1]])
    else .err("couldn't determine length of `X` components")
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
    if (is_not_null(subset) && any(names(X) %in% subsettable())) {
        n <- get_length_X(X)
        if (is.logical(subset)) {
            if (length(subset) != n) .err("`subset` must have the same length as the other entries")
            if (!any(subset)) .err("All `subset` set to FALSE")
            to_be_subset <- !all(subset)
            subset <- which(subset)
        }
        else if (is.numeric(subset)) {
            if (max(subset) > n) .err("Subset indices cannot be higher than the length of the other entries")
            to_be_subset <- TRUE
        }
        else .err("`subset` must be logical or numeric")
        
        if (to_be_subset) {
            subset_X_internal <- function(x, subset) {
                attrs <- attributes(x)
                attrs_to_subset <- names(attrs)[vapply(attrs, function(a) all(len(a) == n), logical(1L))]
                if (is_not_null(attrs_to_subset)) {
                    subsetted_attrs <- lapply(attrs[attrs_to_subset],
                                              subset_X_internal, subset = subset)
                }
                
                if (is_null(x)) out <- x
                else if (inherits(x, "processed.treat")) out <- subset_processed.treat(x, subset)
                else if ((is.matrix(x) || is.data.frame(x))) out <- x[subset, , drop = FALSE]
                else if (is.factor(x)) out <- factor(x[subset], nmax = nlevels(x))
                else if (is.atomic(x)) out <- x[subset]
                else if (is.list(x)) out <- lapply(x, subset_X_internal, subset = subset)
                else out <- x
                
                if (is_not_null(attrs)) {
                    if (inherits(x, "processed.treat")) {
                        out <- process_treat(out, keep_values = TRUE)
                    }
                    else {
                        for (i in setdiff(names(attrs), names(attributes(out)))) {
                            if (i %in% attrs_to_subset) attr(out, i) <- subsetted_attrs[[i]]
                            else attr(out, i) <- attrs[[i]]
                        }
                    }
                }
                
                out
            }
            
            for (i in which(names(X) %in% subsettable())) {
                X[[i]] <- subset_X_internal(X[[i]], subset)
            }
        }
    }
    
    X
}
.mids_complete <- function(data) {
    if (!inherits(data, "mids")) .err("`data` not of class `mids`")
    
    single.complete <- function(data, where = NULL, imp, ell) {
        if (is_null(where)) where <- is.na(data)
        idx <- seq_len(ncol(data))[which(colSums(where) > 0)]
        for (j in idx) {
            if (is_null(imp[[j]])) data[where[, j], j] <- NA
            else data[where[, j], j] <- imp[[j]][, ell]
        }
        
        data
    }
    
    m <- as.integer(data$m)
    idx <- seq_len(m)
    
    mylist <- lapply(idx, function(i) single.complete(data$data, data$where, data$imp, i))
    
    cmp <- data.frame(.imp = rep(idx, each = nrow(data$data)), 
                      .id = rep.int(seq_len(nrow(data$data)), length(idx)), 
                      do.call("rbind", mylist))
    
    row.names(cmp) <- {
        if (is.integer(attr(data$data, "row.names"))) seq_len(nrow(cmp))
        else as.character(seq_len(nrow(cmp)))
    }
    
    cmp
}
length_imp_process <- function(vectors = NULL, data.frames = NULL, lists = NULL,
                               imp = NULL, data = NULL, original.call.to = NULL,
                               env = sys.frame(-1)) {
    #Processes imp and assigns it to imp in env
    #Processes vectors, data.frames, and lists for length
    #If object correspond to one imputation, assigns expanded object in env 
    
    all.objects <- c(vectors, data.frames, lists)
    ensure.equal.lengths <- TRUE
    problematic <- rlang::rep_named(all.objects, FALSE)
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
                        all.objects)
    
    #Process imp further
    if (is_not_null(imp)) {
        
        imp.lengths <- vapply(levels(imp), function(i) sum(imp == i), numeric(1L))
        
        if (all_the_same(imp.lengths)) { #all the same
            unsorted.imp <- is.unsorted(imp)
            for (i in vectors) {
                if (lengths[i] > 0 && lengths[i] != length(imp)) { 
                    if (!all_the_same(imp.lengths)) .err("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation")
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
                    if (!all_the_same(imp.lengths)) .err("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation")
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
                    if (!all_the_same(imp.lengths)) .err("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation")
                    if (lengths[i] == imp.lengths[1]) {
                        assign(i, lapply(get(i, envir = env, inherits = FALSE), function(j) {
                            if (is.factor(j)) {
                                newj <- j[rep(seq_along(j), length(imp.lengths))]
                                if (unsorted.imp) {for (i_ in levels(imp)) newj[imp == i_] <- j}
                                return(newj)
                            }
                            if (is_mat_like(j)) {
                                newj <- j[rep(seq_len(nrow(j)), length(imp.lengths)), , drop = FALSE]
                                if (unsorted.imp) {for (i_ in levels(imp)) newj[imp == i_,] <- j}
                                return(newj)
                            }
                            .err(sprintf("% can only contain vectors or data frames",
                                         add_quotes(i, "`")))
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
            .err(sprintf("%s must have the same number of observations as `imp`",
                         word_list(names(problematic)[problematic], quotes = "`")))
        }
        
        ensure.equal.lengths <- FALSE
        
        assign("imp", imp, pos = env)
    }
    
    #Ensure all input lengths are the same.
    anchor <- {
        if ("treat" %in% all.objects) "treat"
        else if ("treat.list" %in% all.objects) "treat.list"
        else all.objects[which(lengths[all.objects] != 0)[1]]
    }
    if (ensure.equal.lengths) {
        problematic[lengths %nin% c(0, lengths[anchor])] <- TRUE
    }
    if (any(problematic)) {
        if (is_not_null(original.call.to)) anchor <- paste("in the original call to", add_quotes(original.call.to, "`"))
        
        .err("%s must have the same number of observations as %s",
             word_list(names(problematic)[problematic], quotes = "`"),
             anchor)
    }
}
process_stats <- function(stats = NULL, treat) {
    if (is.list(treat) && !is.data.frame(treat)) {
        stats.list <- lapply(treat, function(x) process_stats(stats, x))
        if (all_the_same(vapply(stats.list, attr, character(1L), "type")))
            type <- attr(stats.list[[1]], "type")
        else type <- NULL
        stats <- unique(unlist(stats.list))
        attr(stats, "type") <- type
    }
    else {
        if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
        treat.type <- get.treat.type(treat)
        
        if (treat.type %in% c("binary", "multinomial")) {
            if (is_null(stats)) stats <- getOption("cobalt_stats", "mean.diffs")
            stats <- unique(match_arg(stats, all_STATS("bin"), several.ok = TRUE))
            attr(stats, "type") <- "bin"
        }
        else if (treat.type %in% c("continuous")) {
            if (is_null(stats)) stats <- getOption("cobalt_stats", "correlations")
            stats <- unique(match_arg(stats, all_STATS("cont"), several.ok = TRUE))
            attr(stats, "type") <- "cont"
        }
    }
    
    stats
}
process_thresholds <- function(thresholds, stats) {
    if (is_not_null(thresholds)) {
        if (is.list(thresholds)) thresholds <- unlist(thresholds)
        if (!all(is.na(thresholds)) && !is.numeric(thresholds)) {
            .err("`thresholds` must be numeric")
        }
        
        if (rlang::is_named(thresholds)) {
            names(thresholds) <- stats[pmatch(names(thresholds), stats, duplicates.ok = TRUE)]
            thresholds <- thresholds[!is.na(names(thresholds))]
        }
        else {
            names(thresholds) <- stats[1:min(length(stats), length(thresholds))]
        }
        
        thresholds[names(thresholds)] <- as.numeric(thresholds)
    }
    
    as.list(na.rem(thresholds))
}
process_subset <- function(subset, n) {
    if (!is.logical(subset) && !is.numeric(subset)) {
        .err("the argument to `subset` must be a logical or numeric vector")
    }
    if (is.numeric(subset)) {
        if (any(abs(subset) > n)) .err("numeric values for `subset` cannot be larger than the number of units")
        subset <- subset[!is.na(subset) & subset != 0]
        if (any(subset < 0) && any(subset > 0)) .err("positive and negative indices cannot be mixed with `subset`")
        if (any(abs(subset) > n)) .err("if `subset` is numeric, none of its values can exceed the number of units")
        logical.subset <- rep(any(subset < 0), n)
        logical.subset[abs(subset)] <- !logical.subset[abs(subset)]
        subset <- logical.subset
    }
    if (anyNA(subset)) {
        .wrn("NAs were present in `subset`. Treating them like FALSE")
        subset[is.na(subset)] <- FALSE
    }
    
    subset
}
process_focal <- function(focal, treat) {
    if (is.numeric(focal)) {
        if (can_str2num(treat) && focal %in% str2num(treat)) {focal <- as.character(focal)}
        else if (focal <= length(treat_vals(treat))) focal <- treat_vals(treat)[focal]
        else {
            .err(sprintf("`focal` was specified as %s, but there are only %s treatment groups",
                         focal, length(treat_vals(treat))))
        }
    }
    else if (focal %nin% treat_vals(treat)) {
        .err("the name specified to `focal` is not the name of any treatment group")
    }
    
    focal
}
process_weights <- function(obj = NULL, A = NULL, treat = NULL, covs = NULL,
                            method = character(0), addl.data = list(), ...) {
    A[["x"]] <- NULL
    A[["treat"]] <- NULL
    
    weights <- list()
    if (is_not_null(obj)) {
        weights <- do.call("get.w", c(list(obj, treat = treat), A))
        
        if (is_not_null(weights)) {
            weights <- {
                if (is_mat_like(weights)) as.data.frame(weights)
                else setNames(data.frame(weights), class(obj)[has_method(class(obj), "get.w")][1])
            }
        }
        else {
            weights <- list()
        }
    }
    
    addl.weights <- .process_data_frame("weights", A[["weights"]], treat, covs, addl.data = addl.data, ...)
    if (is_not_null(addl.weights)) {
        if (is_null(A[["method"]])) addl.methods <- rep.int("weighting", ncol(addl.weights))
        else if (length(A[["method"]]) == 1) {
            addl.methods <- rep.int(match_arg(A[["method"]], c("weighting", "matching")), ncol(addl.weights))
        }
        else {
            addl.methods <- match_arg(A[["method"]], c("weighting", "matching"), several.ok = TRUE)
            if (length(addl.methods) != ncol(addl.weights)) {
                .err("Valid inputs to `method` must have length 1 or equal to the number of valid sets of additional weights")
            }
        }
        
        w.names <- c(names(weights), names(addl.weights))
        unique.names <- unique(w.names)
        if (length(unique.names) != length(w.names)) {
            for (i in unique.names) {
                if (sum(w.names == i) > 1) w.names[w.names == i] <- make.unique(c(i, w.names[w.names == i]))[-1]
            } 
        }
        if (is_not_null(weights)) {
            weights <- setNames(cbind(weights, addl.weights), w.names)
            method <- setNames(c(method, addl.methods), w.names)
        }
        else {
            weights <- setNames(addl.weights, w.names)
            method <- setNames(addl.methods, w.names)
        }
    }
    .weight_check(weights)
    
    attr(weights, "method") <- method
    
    weights
}
process_disp <- function(disp = NULL, ...) {
    A <- list(...)
    
    .chk_null_or(disp, .chk_character)
    disp <- {
        if (is_null(disp)) getOption("cobalt_disp")
        else match_arg(disp, acceptable.options()[["disp"]], several.ok = TRUE)
    }
    
    for (d in c("means", "sds")) {
        if (getOption(paste.("cobalt_disp", d), FALSE)) disp <- unique(c(disp, d))
        if (is_not_null(A[[paste.("disp", d)]])) {
            if (!rlang::is_bool(A[[paste.("disp", d)]])) .err(sprintf("`disp.%s` must be `TRUE` or `FALSE`", d))
            disp <- unique(c(disp, d[A[[paste.("disp", d)]]]))
            if (A[[paste.("disp", d)]]) disp <- unique(c(disp, d))
            else disp <- unique(disp[disp != d])
        }
    }
    
    disp
}
process_addl <- function(addl = NULL, datalist = list()) {
    if (is_not_null(addl) && !is.atomic(addl) && !rlang::is_formula(addl) &&
        !is.matrix(addl) && !is.data.frame(addl)) {
        .err("`addl` must be a formula or variable containing the distance values")
    }
    
    data <- do.call("data.frame", unname(clear_null(datalist)))
    
    if (is.atomic(addl) && 
        (!is.character(addl) || is_null(datalist) ||
         length(addl) == nrow(data))) {
        addl <- data.frame(addl = addl)
    }
    else if (is.character(addl)) {
        addl <- reformulate(addl)
    }
    
    addl_t.c <- {
        if (is.data.frame(addl)) {
            .use_tc_fd(data = data, covs = addl, 
                       needs.treat = FALSE, needs.covs = FALSE)
        }
        else {
            .use_tc_fd(formula = addl, data = data, 
                       needs.treat = FALSE, needs.covs = FALSE)
        }
    }
    
    addl_t.c[["covs"]]
}
process_addl.list <- function(addl.list = NULL, datalist = list(), covs.list = list()) {
    datalist <- clear_null(c(datalist, covs.list))
    
    if (is.list(addl.list) && !is.data.frame(addl.list)) {
        if (length(addl.list) != length(covs.list)) .err("`addl` must have an entry for each time point")
        addl.list.out <- lapply(addl.list, process_addl, datalist = datalist)
    }
    else {
        addl <- process_addl(addl.list, datalist = datalist)
        addl.list.out <- lapply(seq_along(covs.list), function(x) addl)
    }
    
    addl.list.out
}
process_distance <- function(distance = NULL, datalist = list(), obj.distance = NULL,
                             obj.distance.name = "distance") {
    
    if (is_not_null(distance) && !is.atomic(distance) && !rlang::is_formula(distance) &&
        !is.matrix(distance) && !is.data.frame(distance)) {
        .err("`distance` must be a formula or variable containing the distance values")
    }
    
    data <- do.call("data.frame", unname(clear_null(datalist)))
    
    if (is.atomic(distance) && 
        (!is.character(distance) || is_null(datalist) ||
         length(distance) == nrow(data))) {
        distance <- data.frame(distance = distance)
    }
    else if (is.character(distance)) {
        distance <- reformulate(distance)
    }
    
    distance_t.c <- {
        if (is.data.frame(distance)) {
            .use_tc_fd(data = data, covs = distance, 
                       needs.treat = FALSE, needs.covs = FALSE)
        }
        else {
            .use_tc_fd(formula = distance, data = data, 
                       needs.treat = FALSE, needs.covs = FALSE)
        }
    }
    
    distance <- distance_t.c[["covs"]]
    
    if (is_not_null(obj.distance) && !all(is.na(obj.distance))) {
        obj.distance <- setNames(data.frame(obj.distance), obj.distance.name)
        obj.distance <- get_covs_from_formula(data = obj.distance)
        distance <- co.cbind(if_null_then(distance, NULL), obj.distance)
    }
    
    distance
}
process_distance.list <- function(distance.list = NULL, datalist = list(),
                                  covs.list = list(), obj.distance = NULL, obj.distance.name = "distance") {
    datalist <- clear_null(c(datalist, covs.list))
    
    if (is_null(obj.distance)) {
        obj.distance <- lapply(seq_along(covs.list), function(x) NULL)
    }
    else if (!is.list(obj.distance) || is.data.frame(obj.distance)) {
        obj.distance <- lapply(seq_along(covs.list), function(x) obj.distance)
    }
    
    if (is_null(distance.list)) {
        distance.list.out <- lapply(seq_along(covs.list), function(x) process_distance(NULL, datalist = datalist, 
                                                                                       obj.distance = obj.distance[[x]], obj.distance.name = obj.distance.name))
    }
    else if (is.list(distance.list) && !is.data.frame(distance.list)) {
        if (length(distance.list) != length(covs.list)) .err("`distance` must have an entry for each time point")
        distance.list.out <- lapply(seq_along(distance.list), function(x) process_distance(distance.list[[x]], datalist = datalist, 
                                                                                           obj.distance = obj.distance[[x]], obj.distance.name = obj.distance.name))
    }
    else {
        distance.list.out <- lapply(seq_along(covs.list), function(x) process_distance(distance.list, datalist = datalist, 
                                                                                       obj.distance = obj.distance[[x]], obj.distance.name = obj.distance.name))
    }
    
    distance.list.out
}
process_focal_and_estimand <- function(focal, estimand, treat, treated = NULL) {
    reported.estimand <- estimand
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    unique.treat <- unique(treat, nmax = switch(treat.type, "binary" = 2, "multinomial" = length(treat)/4))
    
    #Check focal
    if (is_not_null(focal) && (length(focal) > 1L || focal %nin% unique.treat)) {
        .err("The argument supplied to `focal` must be the name of a level of treatment")
    }
    
    if (treat.type == "multinomial") {
        
        if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
            .wrn(sprintf("%s is not compatible with `focal`. Setting `estimand` to \"ATT\"",
                         add_quotes(estimand)))
            reported.estimand <- estimand <- "ATT"
        }
        
        if (estimand == "ATT") {
            if (is_null(focal)) {
                if (is_null(treated) || treated %nin% unique.treat) {
                    .err("When estimand = \"ATT\" for multinomial treatments, an argument must be supplied to `focal`")
                }
                focal <- treated
            }
        }
        else if (estimand == "ATC") {
            if (is_null(focal)) {
                .err("When estimand = \"ATC\" for multinomial treatments, an argument must be supplied to `focal`")
            }
        }
    }
    else if (treat.type == "binary") {
        unique.treat.bin <- unique(binarize(treat), nmax = 2)
        if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
            .wrn(sprintf("%s is not compatible with `focal`. Setting `estimand` to \"ATT\"",
                         add_quotes(estimand)))
            reported.estimand <- estimand <- "ATT"
        }
        
        if (is_null(treated) || treated %nin% unique.treat) {
            if (is_null(focal)) {
                if (all(as.character(unique.treat.bin) == as.character(unique.treat))) {
                    treated <- unique.treat[unique.treat.bin == 1]
                }
                else {
                    if (is.factor(treat)) treated <- levels(treat)[2]
                    else treated <- unique.treat[unique.treat.bin == 1]
                    
                    if (estimand == "ATT") {
                        .msg(sprintf("assuming %s the treated level. If not, supply an argument to `focal`",
                                     word_list(treated, quotes = !is.numeric(treat), is.are = TRUE)))
                        
                    }
                    else if (estimand == "ATC") {
                        .msg(sprintf("assuming %s the control level. If not, supply an argument to `focal`",
                                     word_list(setdiff(unique.treat, treated),
                                               quotes = !is.numeric(treat), is.are = TRUE)))
                    }
                    
                }
                if (estimand == "ATT")
                    focal <- treated
                else if (estimand == "ATC")
                    focal <- setdiff(unique.treat, treated)
            }
            else {
                if (estimand == "ATT")
                    treated <- focal
                else if (estimand == "ATC")
                    treated <- setdiff(unique.treat, focal)
            }
            if (estimand == "ATC") estimand <- "ATT"
        }
        else {
            if (is_null(focal)) {
                if (estimand == "ATT")
                    focal <- treated
                else if (estimand == "ATC")
                    focal <- setdiff(unique.treat, treated)
            }
            if (estimand == "ATC") estimand <- "ATT"
        }
    }
    
    list(focal = as.character(focal),
         estimand = estimand,
         reported.estimand = reported.estimand,
         treated = if (is.factor(treated)) as.character(treated) else treated)
}

#.get_C2
get_ints_from_co.names <- function(co.names) {
    if (is_null(co.names)) return(list())
    
    clear_null(lapply(co.names, function(co) {
        if ("isep" %in% co[["type"]]) {
            co[["type"]] <- c(co[["type"]], "isep")
            which_isep <- which(co[["type"]] == "isep")
            vapply(seq_along(which_isep), function(x) {
                if (x == 1) paste0(co[["component"]][1:(which_isep[1]-1)], collapse = "")
                else paste0(co[["component"]][(which_isep[x-1]+1):(which_isep[x]-1)], collapse = "")
            }, character(1L))
        }
        else NULL
    }))
}
get_treat_from_formula <- function(f, data = NULL, treat = NULL) {
    
    if (is.character(f)) {
        f <- try(as.formula(f), silent = TRUE)
    }
    
    .chk_formula(f)
    
    env <- rlang::f_env(f)
    
    f <- update(f, ~ 0)
    
    #Check if data exists
    if (is_not_null(data)) {
        if (is.data.frame(data)) {
            data.specified <- TRUE
        }
        else {
            .wrn("the argument supplied to `data` is not a data.frame object. Ignoring `data`")
            data <- env
            data.specified <- FALSE
        }
    }
    else {
        data <- env
        data.specified <- FALSE
    }
    
    tryCatch(tt <- terms(f, data = data),
             error = function(e) {
                 .err(conditionMessage(e))
             })
    
    if (rlang::is_formula(tt, lhs = TRUE)) {
        resp.vars.mentioned <- as.character(tt)[2]
        resp.vars.failed <- vapply(resp.vars.mentioned, function(v) {
            test <- tryCatch(eval(str2expression(v), data, env), error = function(e) e)
            if (inherits(test, "simpleError")) {
                if (conditionMessage(test) == sprintf("object '%s' not found", v)) return(TRUE)
                .err(conditionMessage(test), tidy = FALSE)
            }
            if (is.function(test)) .err(sprintf("invalid type (function) for variable '%s'", v))
            
            is_null(test)
        }, logical(1L))
        
        if (any(resp.vars.failed)) {
            if (is_null(treat)) .err(sprintf("the given response variable, %s, is not a variable in %s",
                                             add_quotes(as.character(tt)[2]),
                                             word_list(c("`data`", "the global environment")[c(data.specified, TRUE)], "or")))
            tt <- delete.response(tt)
        }
    }
    else resp.vars.failed <- TRUE
    
    if (any(!resp.vars.failed)) {
        treat.name <- resp.vars.mentioned[!resp.vars.failed][1]
        treat <- eval(str2expression(treat.name), data, env)
    }
    else {
        treat.name <- NULL
    }
    
    attr(treat, "treat.name") <- treat.name
    
    treat
}
get_covs_from_formula <- function(f, data = NULL, factor_sep = "_", int_sep = " * ") {
    
    rebuild_f <- function(ttfactors, tics = FALSE) {
        #Set tics = TRUE if returned formula is used with tmpcovs,
        #i.e., a model.frame. If used with data, leave FALSE
        if (tics) rownames(ttfactors)[!startsWith(rownames(ttfactors), "`")] <- add_quotes(rownames(ttfactors)[!startsWith(rownames(ttfactors), "`")], "`")
        as.formula(paste("~ 0 +", paste(vapply(seq_len(ncol(ttfactors)), 
                                               function(x) paste0(rownames(ttfactors)[ttfactors[,x] > 0], collapse = ":"),
                                               character(1L)), 
                                        collapse = "+")))
    }
    append.ttfactor <- function(ttfactor, term, after) {
        addrow <- matrix(0, nrow = length(term), ncol = ncol(ttfactor), 
                         dimnames = list(term, colnames(ttfactor)))
        addcol <- matrix(0, nrow = nrow(ttfactor) + length(term), ncol = length(term),
                         dimnames = list(c(rownames(ttfactor), term), term))
        addcol[-seq_len(nrow(ttfactor)), ] <- diag(length(term))
        
        ttfactor <- rbind(ttfactor, addrow)
        if (after == 0) {
            return(cbind(addcol, ttfactor))
        }
        if (after == ncol(ttfactor)) {
            return(cbind(ttfactor, addcol))
        }
        
        cbind(ttfactor[,seq_len(after), drop = FALSE], 
              addcol, 
              ttfactor[,-seq_len(after), drop = FALSE])
        
    }
    
    #Check if data exists
    data.specified <- FALSE
    if (!is.null(data)) {
        if (is.matrix(data)) data <- as.data.frame.matrix(data)
        
        if (!is.data.frame(data)) {
            .err("the argument supplied to `data` must be a data.frame object")
        }
        
        data.specified <- TRUE
    }
    
    if (missing(f)) f <- f.build(names(data))
    else {
        if (is.character(f)) {
            f <- try(as.formula(f), silent = TRUE)
        }
        .chk_formula(f)
    }
    
    env <- rlang::f_env(f)
    if (!data.specified) data <- env
    
    # rlang::f_lhs(f) <- NULL
    
    tryCatch(tt <- terms(f, data = data),
             error = function(e) {
                 if (conditionMessage(e) == "'.' in formula and no 'data' argument") {
                     .err("'.' is not allowed in formulas")
                 }
                 else .err(conditionMessage(e))
             })
    
    #Process RHS 
    tt.covs <- delete.response(tt)
    attr(tt.covs, "intercept") <- 0
    
    ttfactors <- attr(tt.covs, "factors")
    ttvars <- setNames(vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1], rownames(ttfactors))
    
    rhs.df.type <- setNames(vapply(ttvars, function(v) {
        if (length(dim(try(eval(str2expression(add_quotes(v, "`")), data, env), silent = TRUE))) == 2) "lit"
        else if (length(dim(try(eval(str2expression(v), data, env), silent = TRUE))) == 2) "exp"
        else "not.a.df"
    }, character(1L)), ttvars)
    
    rhs.df <- setNames(rhs.df.type != "not.a.df", ttvars)
    
    if (any(rhs.df)) {
        term_is_interaction <- apply(ttfactors, 2, function(x) sum(x != 0) > 1)
        if (any(vapply(seq_along(ttvars)[rhs.df], function(x) any(ttfactors[x,] != 0 & term_is_interaction), logical(1L)))) {
            .err("interactions with data.frames are not allowed in the input formula")
        }
        addl.dfs <- setNames(lapply(ttvars[rhs.df], function(v) {
            if (rhs.df.type[v] == "lit") df <- eval(str2expression(add_quotes(v, "`")), data, env)
            else df <- eval(str2expression(v), data, env)
            
            if (inherits(df, "rms")) {
                df <- setNames(as.data.frame.matrix(as.matrix(df)), colnames(df))
                return(df)
            }
            if (is.data.frame(df)) {
                #Deal with the fact that data.frames may contain matrices and data.frames, which
                #may contain data/frames, and so on
                non.vec.col <- which(vapply(df, function(x) is_not_null(dim(x)), logical(1L)))
                while (is_not_null(non.vec.col)) {
                    for (i in non.vec.col) {
                        if (NCOL(df[[i]]) == 1 && is_null(colnames(df[[i]]))) colnames(df[[i]]) <- names(df)[i]
                        else if (can_str2num(colnames(df[[i]]))) colnames(df[[i]]) <- paste(names(df)[i], colnames(df[[i]]), sep = "_")
                    }
                    names(df)[non.vec.col] <- ""
                    
                    df <- as.data.frame(do.call("cbind", df))
                    non.vec.col <- which(vapply(df, function(x) is_not_null(dim(x)), logical(1L)))
                }
                
                if (ncol(df) == 1 && is_null(colnames(df))) colnames(df) <- v
                else if (can_str2num(colnames(df))) colnames(df) <- paste(v, colnames(df), sep = "_")
                return(df)
            }
            
            if (ncol(df) == 1 && is_null(colnames(df))) colnames(df) <- v
            else if (can_str2num(colnames(df))) colnames(df) <- paste(v, colnames(df), sep = "_")
            
            as.data.frame(df)
        }),
        ttvars[rhs.df])
        
        for (i in colnames(ttfactors)[colnames(ttfactors) %in% names(ttvars)[rhs.df]]) {
            for (j in seq_len(ncol(addl.dfs[[ttvars[i]]]))) {
                if (names(addl.dfs[[ttvars[i]]])[j] %in% c(ttvars[!rhs.df], unlist(lapply(addl.dfs[seq_len(which(names(addl.dfs) == ttvars[i]) - 1)], names)))) {
                    names(addl.dfs[[ttvars[i]]])[j] <- paste0(ttvars[i], "_", names(addl.dfs[[ttvars[i]]])[j])
                }
            }
            ind <- which(colnames(ttfactors) == i)
            ttfactors <- append.ttfactor(ttfactors,
                                         add_quotes(names(addl.dfs[[ttvars[i]]]), "`"),
                                         ind)[,-ind, drop = FALSE]
        }
        
        if (data.specified) {
            data <- do.call("cbind", unname(c(addl.dfs, list(data))))
        }
        else {
            data <- do.call("cbind", unname(addl.dfs))
            data.specified <- TRUE
        }
        
        new.form <- rebuild_f(ttfactors)
        tt.covs <- terms(new.form, data = data)
        
        ttfactors <- attr(tt.covs, "factors")
        ttvars <- setNames(vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1], rownames(ttfactors))
    }
    
    #Check to make sure variables are valid
    original_ttvars <- rownames(ttfactors)
    for (i in seq_along(rownames(ttfactors))) {
        #Check if evaluable
        #If not, check if evaluable after changing to literal using ``
        #If not, stop()
        #If eventually evaluable, check if function
        
        evaled.var <- try(eval(str2expression(rownames(ttfactors)[i]), data, env), silent = TRUE)
        if (null_or_error(evaled.var)) {
            evaled.var <- try(eval(str2expression(add_quotes(rownames(ttfactors)[i], "`")), data, env), silent = TRUE)
            if (null_or_error(evaled.var)) {
                ee <- conditionMessage(attr(evaled.var, "condition"))
                if (startsWith(ee, "object '") && endsWith(ee, "' not found")) {
                    v <- sub("object '([^']+)' not found", "\\1", ee)
                    .err(sprintf("the variable \"%s\" cannot be found. Be sure it is entered correctly or supply a dataset that contains this varialble to `data`", v))
                }
                
                .err(ee)
            }
            rownames(ttfactors)[i] <- add_quotes(rownames(ttfactors)[i], "`")
        }
        
        if (is.function(evaled.var)) {
            .err(sprintf("invalid type (function) for variable '%s'", rownames(ttfactors)[i]))
        }
    }
    
    if (!identical(original_ttvars, rownames(ttfactors))) {
        new.form <- rebuild_f(ttfactors)
        tt.covs <- terms(new.form, data = data)
        
        ttfactors <- attr(tt.covs, "factors")
        ttvars <- vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1]
    }
    
    tryCatch({tmpcovs <- model.frame2(tt.covs, data)},
             error = function(e) .err(conditionMessage(e)))
    
    for (i in ttvars) {
        if (is_binary(tmpcovs[[i]])) tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = 2)
        else {
            if (is.character(tmpcovs[[i]])) tmpcovs[[i]] <- factor(tmpcovs[[i]])
            if (is.factor(tmpcovs[[i]])) {
                if (nlevels(tmpcovs[[i]]) == 1) tmpcovs[[i]] <- factor(tmpcovs[[i]], levels = c(paste0(".", levels(tmpcovs[[i]])), levels(tmpcovs[[i]])))
                else tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = nlevels(tmpcovs[[i]]))
            }
        }
    }
    
    #Process NAs: make NA variables
    if (anyNA(tmpcovs)) {
        vars_with_NA <- colnames(tmpcovs)[anyNA_col(tmpcovs)]
        for (i in rev(vars_with_NA)) {
            #Find which of ttlabels i first appears, and put `i: <NA>` after it
            for (x in seq_along(colnames(ttfactors))) {
                if (i %in% c(colnames(ttfactors)[x], all.vars(str2expression(colnames(ttfactors)[x])))) {
                    ind <- x
                    break
                }
            }
            ttfactors <- append.ttfactor(ttfactors, 
                                         paste0("`", i, ":<NA>`"),
                                         ind)
            
            tmpcovs[[paste0(i, ":<NA>")]] <- as.numeric(is.na(tmpcovs[[i]]))
        }
        new.form <- rebuild_f(ttfactors, tics = TRUE)
        
        tt.covs <- terms(new.form, data = tmpcovs)
        ttfactors <- attr(tt.covs, "factors")
        ttvars <- vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1]
        
        na_vars <- paste0(vars_with_NA, ":<NA>")
        
        tryCatch({tmpcovs <- model.frame2(tt.covs, tmpcovs)},
                 error = function(e) .err(conditionMessage(e)))
        
        for (i in setdiff(ttvars, na_vars)) {
            if (is_binary(tmpcovs[[i]])) tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = 2)
            else {
                if (is.character(tmpcovs[[i]])) tmpcovs[[i]] <- factor(tmpcovs[[i]])
                if (is.factor(tmpcovs[[i]])) {
                    if (nlevels(tmpcovs[[i]]) == 1) tmpcovs[[i]] <- factor(tmpcovs[[i]], levels = c(paste0(".",levels(tmpcovs[[i]])), levels(tmpcovs[[i]])))
                    else tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = nlevels(tmpcovs[[i]]))
                }
            }
        }
    }
    else {
        na_vars <- character(0)
    }
    
    #Re-check ttfactors
    original_ttvars <- rownames(ttfactors)
    for (i in seq_along(rownames(ttfactors))) {
        #Check if evaluable in tmpcovs
        #If not, check if evaluable i tmpcovs after changing to literal using ``
        #If not, stop() (shouldn't occur)
        
        evaled.var <- try(eval(str2expression(rownames(ttfactors)[i]), tmpcovs), silent = TRUE)
        if (null_or_error(evaled.var)) {
            evaled.var <- try(eval(str2expression(add_quotes(rownames(ttfactors)[i], "`")), tmpcovs), silent = TRUE)
            if (null_or_error(evaled.var)) {
                .err(conditionMessage(attr(evaled.var, "condition")))
            }
            rownames(ttfactors)[i] <- add_quotes(rownames(ttfactors)[i], "`")
        }
    }
    
    if (!identical(original_ttvars, rownames(ttfactors))) {
        new.form <- rebuild_f(ttfactors)
        tt.covs <- terms(new.form, data = data)
        
        ttfactors <- attr(tt.covs, "factors")
        ttvars <- vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1]
    }
    
    tmpcovs <- model.frame2(tt.covs, data = tmpcovs, drop.unused.levels = TRUE)
    
    #Check for infinite values
    covs.with.inf <- vapply(tmpcovs, function(x) is.numeric(x) && any(!is.na(x) & !is.finite(x)), logical(1L))
    if (any(covs.with.inf)) {
        s <- if (sum(covs.with.inf) == 1) c("", "s") else c("s", "")
        .err(sprintf("the variable%s %s contain%s non-finite values, which are not allowed",
                     s[1],
                     word_list(names(tmpcovs)[covs.with.inf], quotes = 1),
                     s[2]))
    }
    
    attr(tt.covs, "intercept") <- 1 #Add intercept to correctly process single-level factors
    mm <- model.matrix(tt.covs, data = tmpcovs, 
                       contrasts.arg = lapply(Filter(is.factor, tmpcovs),
                                              function(x) contrasts(x, contrasts = all_the_same(x))))
    
    rownames(ttfactors) <- trim_string(rownames(ttfactors), "`")
    
    mmassign <- attr(mm, "assign")[-1]
    mmassign2 <- setNames(factor(mmassign, levels = sort(unique(mmassign), na.last = TRUE),
                                 labels = colnames(ttfactors)), colnames(mm)[-1])
    
    vars_in_each_term <- setNames(lapply(colnames(ttfactors), function(x) {
        rownames(ttfactors)[ttfactors[,x] != 0]
    }), colnames(ttfactors))
    all_factor_levels <- lapply(vars_in_each_term, function(v) {
        do.call("expand.grid", c(clear_null(setNames(lapply(v, function(fa) colnames(attr(mm, "contrasts")[[fa]])), v)),
                                 list(stringsAsFactors = FALSE)))
    })
    expanded <- setNames(lapply(seq_along(mmassign2), function(x) {
        terms <- vars_in_each_term[[mmassign2[x]]]
        k <- sum(seq_along(mmassign2) <= x & mmassign2 == mmassign2[x])
        setNames(lapply(terms, function(t) {
            if (t %in% names(all_factor_levels[[mmassign2[x]]])) {
                all_factor_levels[[mmassign2[x]]][[t]][k]
            }
            else character(0)
        }), terms)
    }), names(mmassign2))
    
    #component types: base, fsep, isep, power, na, level
    co.names <- lapply(expanded, function(x) {
        Reduce(function(x1, x2 = NULL) {
            list(component = c(x1[["component"]], int_sep, x2[["component"]]),
                 type = c(x1[["type"]], "isep", x2[["type"]]))
        },
        lapply(seq_along(x), function(i) {
            base <- gsub("`", "", names(x)[i], fixed = TRUE)
            if (base %in% na_vars) {
                base <- substr(base, 1, nchar(base) - 5)
                list(component = c(base, ":<NA>"),
                     type = c("base", "na"))
            }
            else if (is_null(x[[i]])) {
                list(component = base,
                     type = "base")
            }
            else {
                list(component = c(base, factor_sep, x[[i]]),
                     type = c("base", "fsep", "level"))
            }
        }))
    })
    
    names(co.names) <- vapply(co.names, function(x) paste0(x[["component"]], collapse = ""), character(1L))
    covs <- clear_attr(mm)[,-1, drop = FALSE] #drop the intercept
    
    attr(co.names, "seps") <- c(factor = factor_sep, int = int_sep)
    attr(covs, "co.names") <- co.names
    
    colnames(covs) <- names(co.names)
    
    covs
}
.get_C2 <- function(covs = NULL, int = FALSE, poly = 1, addl = NULL, distance = NULL,
                    treat = NULL, cluster = NULL, drop = TRUE, factor_sep = "_",
                    int_sep = " * ", ...) {
    #gets C data.frame, which contains all variables for which balance is to be assessed. Used in balance.table.
    if (inherits(covs, "processed_C")) return(covs)
    if (is_null(covs)) drop <- FALSE
    
    .chk_string(factor_sep)
    .chk_string(int_sep)
    
    #Process int and poly
    .chk_whole_number(poly)
    .chk_gte(poly, 1)
    poly <- round(poly)
    
    if (is.numeric(int)) {
        if (!chk::vld_whole_number(int) ||
            !chk::vld_gt(int, 1)) {
            .err("`int` must be TRUE, FALSE, or a numeric (integer) value greater than 1")
        }
        if (int > poly) poly <- int
        int <- TRUE
    }
    .chk_flag(int)
    
    A <- list(...)
    if (!chk::vld_flag(A[["center"]])) A[["center"]] <- getOption("cobalt_center", default = FALSE)
    if (!chk::vld_flag(A[["orth"]])) A[["orth"]] <- getOption("cobalt_orth", default = FALSE)
    
    co.names <- attr(covs, "co.names")
    seps <- attr(co.names, "seps")
    
    if (is_not_null(addl)) {
        addl.co.names <- attr(addl, "co.names")
        
        same.name <- names(addl.co.names) %in% names(co.names)
        addl <- addl[,!same.name, drop = FALSE]
        addl.co.names[same.name] <- NULL
        
        #Remove variables in addl that are redundant with covs
        if (drop && getOption("cobalt_remove_perfect_col", max(ncol(addl), ncol(covs)) <= 900)) {
            redundant.var.indices <- find_perfect_col(addl, covs)
            if (is_not_null(redundant.var.indices)) {
                addl <- addl[,-redundant.var.indices, drop = FALSE]
                addl.co.names[redundant.var.indices] <- NULL
            }
        }
        
        covs <- cbind(covs, addl)
        co.names <- c(co.names, addl.co.names)
    } 
    
    #Drop colinear with treat
    if (drop) {
        test.treat <- is_not_null(treat) && get.treat.type(treat) != "continuous"
        
        # test.cluster <- is_not_null(cluster) && !all_the_same(cluster, na.rm = FALSE)
        
        drop_vars <- vapply(seq_len(ncol(covs)), 
                            function(i) {
                                # if (all_the_same(covs[,i], na.rm = FALSE)) return(TRUE)
                                if (anyNA(covs[,i])) return(FALSE)
                                if (test.treat && equivalent.factors2(covs[,i], treat)) return(TRUE)
                                FALSE
                            }, logical(1L))
        
        if (any(drop_vars)) {
            covs <- covs[,!drop_vars, drop = FALSE]
            co.names[drop_vars] <- NULL
        }
    }
    
    C_list <- list(C = covs)
    co_list <- list(C = co.names)
    rm(co.names)
    
    if (int || (poly > 1)) {
        nsep <- 1
        
        #Exclude NA and ints from interactions and poly
        exclude <- vapply(co_list[["C"]], function(x) any(c("na", "isep") %in% x[["type"]]), logical(1L))
        
        new <- .int_poly_f2(C_list[["C"]], ex = exclude, int = int, poly = poly, center = A[["center"]], 
                            orth = A[["orth"]], sep = rep.int(seps["int"], nsep), co.names = co_list[["C"]])
        
        C_list[["int.poly"]] <- new
        co_list[["int.poly"]] <- attr(new, "co.names")
        names(co_list[["int.poly"]]) <- vapply(co_list[["int.poly"]], 
                                               function(x) paste0(x[["component"]], collapse = ""), character(1L))
    }
    
    #Drop 0 category of 0/1 variables and rename 1 category
    if (drop) {
        drop_0_1 <- rep(NA, length(co_list[["C"]]))
        for (i in seq_along(co_list[["C"]])) {
            if (!is.na(drop_0_1[i])) next
            
            if ("isep" %in% co_list[["C"]][[i]][["type"]] || "fsep" %nin% co_list[["C"]][[i]][["type"]]) {
                drop_0_1[i] <- FALSE
                next
            }
            
            which_are_buddies <- which(vapply(co_list[["C"]], function(j) "isep" %nin% j[["type"]] && 
                                                  "fsep" %in% j[["type"]] &&
                                                  j[["component"]][j[["type"]] == "base"][1] == co_list[["C"]][[i]][["component"]][co_list[["C"]][[i]][["type"]] == "base"][1], 
                                              logical(1L)))
            
            buddies <- co_list[["C"]][which_are_buddies]
            
            if (!is.null(cluster)) {
                #Remove variables perfectly redundant with cluster
                unsplit_var <- unsplitfactor(as.data.frame(C_list[["C"]][, which_are_buddies, drop = FALSE]),
                                             buddies[[1]][["component"]][buddies[[1]][["type"]] == "base"],
                                             sep = attr(co_list[["C"]], "seps")["factor"])[[1]]
                tab <- table(cluster, unsplit_var)
                tab <- tab[rowSums(tab == 0) < ncol(tab),,
                           drop = FALSE]
                
                if (all(rowSums(tab > 0) == 1)) {
                    drop_0_1[which_are_buddies] <- TRUE
                    next
                }
            }
            
            if (length(buddies) > 2) {
                drop_0_1[which_are_buddies] <- FALSE
                next
            }
            
            buddy_is_0 <- vapply(buddies, function(x) x[["component"]][x[["type"]] == "level"] %in% c("0", "FALSE"), logical(1L))
            buddy_is_1 <- vapply(buddies, function(x) x[["component"]][x[["type"]] == "level"] %in% c("1", "TRUE"), logical(1L))
            if (!all(buddy_is_0 | buddy_is_1)) {
                drop_0_1[which_are_buddies] <- c(TRUE, FALSE)
                next
            }
            
            drop_0_1[which_are_buddies[buddy_is_0]] <- TRUE
            drop_0_1[which_are_buddies[buddy_is_1]] <- FALSE
            
            buddy_1 <- which_are_buddies[buddy_is_1]
            co_list[["C"]][[buddy_1]][["component"]] <- co_list[["C"]][[buddy_1]][["component"]][co_list[["C"]][[buddy_1]][["type"]] == "base"][1]
            co_list[["C"]][[buddy_1]][["type"]] <- "base"
        }
        
        if (any(drop_0_1)) {
            C_list[["C"]] <- C_list[["C"]][,!drop_0_1, drop = FALSE]
            co_list[["C"]][drop_0_1] <- NULL
        }
    }
    
    names(co_list[["C"]]) <- vapply(co_list[["C"]], function(x) paste0(x[["component"]], collapse = ""), character(1L))
    
    if (is_not_null(distance)) {
        if (anyNA(distance, recursive = TRUE)) .err("missing values are not allowed in the distance measure")
        
        distance.co.names <- attr(distance, "co.names")
        
        same.name <- names(distance.co.names) %in% unlist(lapply(co_list, names))
        if (any(same.name)) {
            distance <- distance[,!same.name, drop = FALSE]
            distance.co.names[same.name] <- NULL
        }
        
        unique.distance.names <- unique(names(distance.co.names))
        distance <- distance[,unique.distance.names, drop = FALSE]
        distance.co.names <- distance.co.names[unique.distance.names]
        
        C_list[["distance"]] <- distance
        co_list[["distance"]] <- distance.co.names
    }
    
    # C_list <- clear_null(C_list)
    # co_list <- clear_null(co_list)
    
    #Remove duplicate & redundant variables
    if (drop) {
        for (x in setdiff(names(C_list), "distance")) {
            
            #Remove self-redundant variables
            # if (getOption("cobalt_remove_perfect_col", ncol(C_list[[x]]) <= 900)) {
            #     
            #     redundant.var.indices <- find_perfect_col(C_list[[x]])
            #     if (is_not_null(redundant.var.indices)) {
            #         C_list[[x]] <- C_list[[x]][,-redundant.var.indices, drop = FALSE]
            #         co_list[[x]][redundant.var.indices] <- NULL
            #     }
            # }
            if (x != "C") {
                #Remove variables in C that have same name as other variables
                if (any(dups <- names(co_list[["C"]]) %in% co_list[[x]])) {
                    C_list[["C"]] <- C_list[["C"]][,!dups, drop = FALSE]
                    co_list[["C"]][dups] <- NULL 
                }
                #Remove variables in C that are redundant with current piece
                # if (getOption("cobalt_remove_perfect_col", max(ncol(C_list[[x]]), ncol(C_list[["C"]])) <= 900)) {
                #     redundant.var.indices <- find_perfect_col(C_list[["C"]], C_list[[x]])
                #     if (is_not_null(redundant.var.indices)) {
                #         C_list[["C"]] <- C_list[["C"]][,-redundant.var.indices, drop = FALSE]
                #         co_list[["C"]][redundant.var.indices] <- NULL
                #     }
                # }
            }
        }
    }
    
    C <- do.call("cbind", clear_null(C_list[c("distance", "C", "int.poly")]))
    co.names <- do.call("c", co_list[c("distance", "C", "int.poly")])
    
    for (i in seq_along(co.names)) {
        co.names[[i]]$component[co.names[[i]]$type == "fsep"] <- factor_sep
        co.names[[i]]$component[co.names[[i]]$type == "isep"] <- int_sep
    }
    seps["factor"] <- factor_sep
    seps["int"] <- int_sep
    
    colnames(C) <- names(co.names) <- vapply(co.names, function(x) paste0(x[["component"]], collapse = ""), character(1L))
    
    attr(co.names, "seps") <- seps
    
    attr(C, "co.names") <- co.names
    
    attr(C, "missing.ind") <- colnames(C)[vapply(co.names, function(x) "na" %in% x[["type"]], logical(1L))]
    if ("distance" %in% names(C_list)) attr(C, "distance.names") <- names(co_list[["distance"]])
    
    attr(C, "var_types") <- .get_types(C)
    class(C) <- c(class(C), "processed_C")
    
    C
}
.int_poly_f2 <- function(mat, ex = NULL, int = FALSE, poly = 1, center = FALSE,
                         orth = FALSE, sep = " * ", co.names = NULL) {
    #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
    #Only to be used in base.bal.tab; for general use see int.poly()
    #mat=matrix input
    #ex=names of variables to exclude in interactions and polynomials; a subset of df
    #int=whether to include interactions or not; currently only 2-way are supported
    #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
    #orth=use orthogonal polynomials
    #nunder=number of underscores between variables
    
    cn <- is_not_null(co.names)
    if (is_null(ex)) ex <- rep(FALSE, ncol(mat))
    d <- mat
    
    binary.vars <- is_binary_col(d)
    interaction.vars <- {
        if (cn) vapply(colnames(d), function(x) "isep" %in% co.names[[x]][["type"]], logical(1L))
        else rep(FALSE, ncol(d))
    }
    
    if (center) {
        d[,!binary.vars] <- center(d[, !binary.vars, drop = FALSE])
    }
    nd <- NCOL(d)
    
    if (poly > 1) {
        poly_terms <- poly_co.names <- make_list(poly-1)
        no.poly <- binary.vars | interaction.vars | ex
        npol <- nd - sum(no.poly)
        if (npol > 0) {
            for (i in 2:poly) {
                poly_terms[[i - 1]] <- apply(d[, !no.poly, drop = FALSE], 2, function(x) {
                    if (orth) poly(x, degree = poly)[,i] else x^i
                })
                poly_co.names[[i - 1]] <- {
                    if (cn) {
                        lapply(colnames(d)[!no.poly], function(x) {
                            list(component = c(co.names[[x]][["component"]], num_to_superscript(i)), 
                                 type = c(co.names[[x]][["type"]], "power"))
                        })
                    }
                    else {
                        paste0(colnames(d)[!no.poly], num_to_superscript(i))
                    }
                }
            }
        }
    }
    else {
        poly_terms <- poly_co.names <- list()
    }
    
    if (int && nd > 1) {
        int_terms <- int_co.names <- make_list(1)
        ints_to_make <- combn(colnames(d)[!ex], 2, simplify = FALSE)
        
        #Don't make ints out of ints that already exist
        # ints_that_already_exist <- get_ints_from_co.names(co.names[interaction.vars])
        # ints_to_make[vapply(ints_to_make, function(x) {
        #   any(vapply(ints_that_already_exist, function(y) {
        #     identical(sort(x), sort(y))
        #   }, logical(1L)))
        # }, logical(1L))] <- NULL
        
        #Don't make ints out of multiple members of the same categorical variable
        ints_to_make[vapply(ints_to_make, function(x) {
            "fsep" %in% co.names[[x[1]]][["type"]] && 
                "fsep" %in% co.names[[x[2]]][["type"]] &&
                identical(co.names[[x[1]]][["component"]][co.names[[x[1]]][["type"]] == "base"],
                          co.names[[x[2]]][["component"]][co.names[[x[2]]][["type"]] == "base"])
        }, logical(1L))] <- NULL
        
        int_terms[[1]] <- do.call("cbind", lapply(ints_to_make, function(i) d[,i[1]] * d[,i[2]]))
        
        if (cn) int_co.names[[1]] <- lapply(ints_to_make, function(x) list(component = c(co.names[[x[1]]][["component"]], sep, co.names[[x[2]]][["component"]]),
                                                                           type = c(co.names[[x[1]]][["type"]], "isep", co.names[[x[2]]][["type"]])))
        else int_co.names[[1]] <- vapply(ints_to_make, paste, character(1L), collapse = sep)
    }
    else {
        int_terms <- int_co.names <- list()
    }
    
    out <- do.call("cbind", c(poly_terms, int_terms))
    out_co.names <- c(do.call("c", poly_co.names), do.call("c", int_co.names))
    
    if (cn) {
        names(out_co.names) <- vapply(out_co.names, 
                                      function(x) paste0(x[["component"]], collapse = ""), character(1))
        colnames(out) <- names(out_co.names)
    }
    else {
        colnames(out) <- unlist(out_co.names)
    }
    
    #Remove single values
    single_value <- apply(out, 2, all_the_same)
    out <- out[, !single_value, drop = FALSE]
    if (cn && is_not_null(out)) {
        attr(out, "co.names") <- out_co.names[!single_value]
    }
    
    out
}
co.cbind <- function(..., deparse.level = 1) {
    args <- clear_null(list(...))
    if (length(args) <= 1) return(args[[1]])
    
    co.names.list <- lapply(args, attr, "co.names")
    
    seps <- attr(co.names.list[[1]], "seps")
    
    out <- do.call("cbind", args)
    
    attr(out, "co.names") <- do.call("c", co.names.list)
    attr(attr(out, "co.names"), "seps") <- seps
    colnames(out) <- names(attr(out, "co.names")) <- vapply(attr(out, "co.names"), function(x) paste0(x[["component"]], collapse = ""), character(1L))
    
    out
}
co.rbind <- function(..., deparse.level = 1) {
    args <- clear_null(list(...))
    if (length(args) <= 1) return(args[[1]])
    
    co.names <- attr(args[[1]], "co.names")
    
    out <- do.call("rbind", args)
    
    attr(out, "co.names") <- co.names
    colnames(out) <- names(attr(out, "co.names")) <- vapply(attr(out, "co.names"), function(x) paste0(x[["component"]], collapse = ""), character(1L))
    
    out
}

df_clean <- function(df) {
    if (!is.list(df)) df <- as.data.frame(df)
    #If a data.frame has matrix columns, clean those into vector columns and rename
    if (any(vapply(df, function(x) is_not_null(dim(x)), logical(1L)))) {
        rn <- rownames(df)
        
        df <- do.call("cbind", lapply(seq_along(df), function(x) {
            if (is_null(dim(df[[x]]))) {
                setNames(data.frame(df[[x]]), names(df)[x])
            }
            else {
                if (inherits(df[[x]], "rms")) {
                    setNames(as.data.frame.matrix(as.matrix(df[[x]])), colnames(df[[x]]))
                }
                else if (can_str2num(colnames(df[[x]]))) {
                    setNames(as.data.frame(df[[x]]), paste(names(df)[x], colnames(df[[x]]), sep = "_"))
                }
                else {
                    setNames(as.data.frame(df[[x]]), colnames(df[[x]]))
                }
            }
        }))
        
        rownames(df) <- rn
    }
    
    df
}
.get_types <- function(C) {
    vapply(colnames(C), function(x) {
        if (any(attr(C, "distance.names") == x)) "Distance"
        else if (all_the_same(C[,x]) || is_binary(C[,x]))  "Binary"
        else "Contin."
    }, character(1))
}
find_perfect_col <- function(C1, C2 = NULL, fun = stats::cor) {
    
    #Finds indices of redundant vars in C1.
    C1.no.miss <- C1[,colnames(C1) %nin% attr(C1, "missing.ind"), drop = FALSE]
    if (is_null(C2)) {
        use <- if (anyNA(C1)) "pairwise.complete.obs" else "everything"
        suppressWarnings(C.cor <- fun(C1.no.miss, use = use))
        s <- !lower.tri(C.cor, diag = TRUE) & !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
    }
    else {
        C2.no.miss <- C2[,colnames(C2) %nin% attr(C2, "missing.ind"), drop = FALSE]
        use <- if (anyNA(C1) || anyNA(C2)) "pairwise.complete.obs" else "everything"
        suppressWarnings(C.cor <- fun(C2.no.miss, y = C1.no.miss, use = use))
        s <- !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
    }
    
    which(colSums(s) > 0)
}

model.frame2 <- function(formula, data = NULL, na.action = "na.pass", ...) {
    withCallingHandlers(force(data),
                        error = function(e) .err(conditionMessage(e)),
                        warning = function(w) .wrn(conditionMessage(w)))
    
    tryCatch({
        mf  <- stats::model.frame(formula, data = data, na.action = na.action, ...)
    },
    error = function(e) {
        ee <- conditionMessage(e)
        if (startsWith(ee, "object '") && endsWith(ee, "' not found")) {
            v <- sub("object '([^']+)' not found", "\\1", ee)
            .err(sprintf("the variable \"%s\" cannot be found. Be sure it is entered correctly or supply a dataset that contains this varialble to `data`", v))
        }
        
        .err(ee)
    })
    
    mf
}

#base.bal.tab
check_if_zero_weights <- function(weights.df, treat = NULL) {
    #Checks if all weights are zero in each treat group for each set of weights
    if (is_not_null(treat)) {
        w.t.mat <- expand.grid(weight_names = colnames(weights.df), 
                               treat_vals = treat_vals(treat), 
                               stringsAsFactors = FALSE)
        if (NROW(w.t.mat) > 0) {
            problems <- vapply(seq_len(NROW(w.t.mat)), function(x) all(check_if_zero(weights.df[treat == w.t.mat[x, "treat_vals"], w.t.mat[x, "weight_names"]])), logical(1L))
            if (any(problems)) {
                prob.w.t.mat <- w.t.mat[problems,]
                if (NCOL(weights.df) == 1) {
                    error <- sprintf("All weights are zero when treat is %s.",
                                     word_list(prob.w.t.mat[, "treat_vals"], "or", quotes = 2))
                }
                else {
                    errors <- vapply(unique(prob.w.t.mat[,"weight_names"]), function(i) {
                        sprintf("%s weights are zero when treat is %s",
                                add_quotes(i),
                                word_list(prob.w.t.mat[prob.w.t.mat[,"weight_names"] == i, "treat_vals"], "or",
                                          quotes = 2))
                    }, character(1L))
                    errors <- paste(c("All", rep("all", length(errors)-1)), errors)
                    error <- paste0(word_list(errors, "and"), ".")
                }
                .err(error, tidy = FALSE)
            }
        }
    }
    else {
        if (length(colnames(weights.df)) > 0) {
            problems <- vapply(colnames(weights.df), function(wn) all(check_if_zero(weights.df[, wn])), logical(1L))
            if (any(problems)) {
                prob.wts <- colnames(weights.df)[problems]
                if (NCOL(weights.df) == 1) {
                    error <- "All weights are zero."
                }
                else {
                    errors <- vapply(prob.wts, function(i) {
                        sprintf("%s weights are zero", add_quotes(i))
                    }, character(1L))
                    errors <- paste(c("All", rep("all", length(errors)-1)), errors)
                    error <- paste0(word_list(errors, "and"), ".")
                }
                .err(error, tidy = FALSE)
            }
        }
    }
}
.baltal <- function(threshold) {
    #threshold: vector of threshold values (i.e., "Balanced"/"Not Balanced")
    threshnames <- names(table(threshold))
    balstring <- threshnames[nchar(threshnames) > 0][1]
    thresh.val <- substring(balstring, 1 + regexpr("[><]", balstring), nchar(balstring))
    b <- data.frame(count=c(sum(threshold == paste0("Balanced, <", thresh.val)), 
                            sum(threshold == paste0("Not Balanced, >", thresh.val))))
    rownames(b) <- c(paste0("Balanced, <", thresh.val), paste0("Not Balanced, >", thresh.val))
    
    b
}
.max_imbal <- function(balance.table, col.name, thresh.col.name, abs_stat) {
    balance.table.clean <- balance.table[balance.table$Type != "Distance" & is.finite(balance.table[, col.name]),]
    maxed <- balance.table.clean[which.max(abs_stat(balance.table.clean[, col.name])), match(c(col.name, thresh.col.name), names(balance.table.clean))]
    
    data.frame(Variable = rownames(maxed), maxed)
}
threshold.summary <- function(compute, thresholds, no.adj, balance.table, weight.names = NULL, agg.fun = NULL) {
    out <- do.call("c", lapply(compute, function(s) make_list(paste.(c("Balanced", "Max.Imbalance"), s))))
    
    for (s in compute) {
        
        if (is_not_null(thresholds[[s]])) {
            
            if (no.adj) {
                out[[paste.("Balanced", s)]] <- .baltal(balance.table[[paste.(STATS[[s]]$Threshold, "Un")]])
                out[[paste.("Max.Imbalance", s)]] <- .max_imbal(balance.table[balance.table[["Type"]]!="Distance", , drop = FALSE], 
                                                                col.name = if (is_null(agg.fun)) paste.(STATS[[s]]$bal.tab_column_prefix, "Un")
                                                                else paste.(firstup(agg.fun), STATS[[s]]$bal.tab_column_prefix, "Un"), 
                                                                thresh.col.name = paste.(STATS[[s]]$Threshold, "Un"), 
                                                                abs_stat = STATS[[s]]$abs)
            }
            else if (length(weight.names) == 1) {
                out[[paste.("Balanced", s)]] <- .baltal(balance.table[[STATS[[s]]$Threshold]])
                out[[paste.("Max.Imbalance", s)]] <- .max_imbal(balance.table[balance.table[["Type"]]!="Distance", , drop = FALSE], 
                                                                col.name = if (is_null(agg.fun)) paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")
                                                                else paste.(firstup(agg.fun), STATS[[s]]$bal.tab_column_prefix, "Adj"), 
                                                                thresh.col.name = STATS[[s]]$Threshold, 
                                                                abs_stat = STATS[[s]]$abs)
            }
            else if (length(weight.names) > 1) {
                out[[paste.("Balanced", s)]] <- setNames(do.call("cbind", lapply(weight.names, function(x) .baltal(balance.table[[paste.(STATS[[s]]$Threshold, x)]]))),
                                                         weight.names)
                out[[paste.("Max.Imbalance", s)]] <- cbind(Weights = weight.names,
                                                           do.call("rbind", lapply(weight.names, function(x) setNames(.max_imbal(balance.table[balance.table[["Type"]]!="Distance", , drop = FALSE], 
                                                                                                                                 col.name = if (is_null(agg.fun)) paste.(STATS[[s]]$bal.tab_column_prefix, x)
                                                                                                                                 else paste.(firstup(agg.fun), STATS[[s]]$bal.tab_column_prefix, x),  
                                                                                                                                 thresh.col.name = paste.(STATS[[s]]$Threshold, x), 
                                                                                                                                 abs_stat = STATS[[s]]$abs),
                                                                                                                      c("Variable", 
                                                                                                                        STATS[[s]]$bal.tab_column_prefix, 
                                                                                                                        STATS[[s]]$Threshold)))),
                                                           stringsAsFactors = FALSE)
            }
        }
        else {
            out[[paste.("Balanced", s)]] <- NULL
            out[[paste.("Max.Imbalance", s)]] <- NULL
        }
    }
    
    out
}

balance.table <- function(C, type, weights = NULL, treat, continuous, binary, s.d.denom, 
                          thresholds = list(), un = FALSE, disp = NULL, stats = NULL, 
                          s.weights = rep(1, length(treat)), abs = FALSE, no.adj = FALSE, 
                          var_types = NULL, s.d.denom.list = NULL, quick = TRUE, ...) {
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    weight.names <- if (no.adj) "Adj" else names(weights)
    
    if (is_not_null(s.d.denom.list)) names(s.d.denom.list) <- weight.names
    if (is_not_null(s.d.denom)) names(s.d.denom) <- weight.names
    
    disp <- c(disp, all_STATS(type)[all_STATS(type) %in% stats])
    compute <- if (quick) disp else c("means", "sds", all_STATS(type)[all_STATS(type) %in% stats])
    
    #B=Balance frame
    Bnames <- c("Type", 
                expand.grid_string(c(if (type == "bin") expand.grid_string(c("M"["means" %in% compute],
                                                                             "SD"["sds" %in% compute]),
                                                                           c("0", "1"), collapse = ".")
                                     else if (type == "cont") c("M"["means" %in% compute], "SD"["sds" %in% compute]), 
                                     unlist(lapply(intersect(compute, all_STATS(type)[!get_from_STATS("adj_only")[get_from_STATS("type") == type]]), function(s) {
                                         c(STATS[[s]]$bal.tab_column_prefix,
                                           if (no.adj && is_not_null(thresholds[[s]])) STATS[[s]]$Threshold)
                                     }))
                ),
                "Un", collapse = "."),
                expand.grid_string(c(if (type == "bin") expand.grid_string(c("M"["means" %in% compute], "SD"["sds" %in% compute]), c("0", "1"), collapse = ".")
                                     else if (type == "cont") c("M"["means" %in% compute], "SD"["sds" %in% compute]), 
                                     unlist(lapply(intersect(compute, all_STATS(type)), function(s) {
                                         c(STATS[[s]]$bal.tab_column_prefix,
                                           if (!no.adj && is_not_null(thresholds[[s]])) STATS[[s]]$Threshold)
                                     }))
                ),
                weight.names, collapse = "."))
    B <- make_df(Bnames, NCOL(C))
    rownames(B) <- colnames(C)
    
    #Set var type (binary/continuous)
    B[["Type"]] <- if_null_then(var_types, .get_types(C))
    bin.vars <- B[["Type"]] == "Binary"
    
    #Means for each group
    if ("means" %in% compute) {
        if (type == "bin") {
            tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
            
            if (un || !quick) {
                for (t in c("0", "1")) {
                    B[[paste.("M", t, "Un")]] <- col_w_mean(C, weights = NULL, s.weights = s.weights, subset = treat == tn01[t])
                }
            }
            
            if (!no.adj && (!quick || "means" %in% disp)) {
                for (i in weight.names) {
                    for (t in c("0", "1")) {
                        B[[paste.("M", t, i)]] <- col_w_mean(C, weights = weights[[i]], s.weights = s.weights, subset = treat == tn01[t])
                    }
                }
            }
        }
        else if (type == "cont") {
            if (un || !quick) {
                B[["M.Un"]] <- col_w_mean(C, weights = NULL, s.weights = s.weights)
            }
            if (!no.adj && (!quick || "means" %in% disp)) {
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
            tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
            
            if (un || !quick) {
                for (t in c("0", "1")) {
                    sds <- rep(NA_real_, NCOL(C))
                    if (any(sd.computable)) {
                        sds[sd.computable] <- col_w_sd(C[, sd.computable,drop = FALSE], weights = NULL, s.weights = s.weights,
                                                       bin.vars = bin.vars[sd.computable], subset = treat == tn01[t])
                    }
                    B[[paste.("SD", t, "Un")]] <- sds
                }
            }
            
            if (!no.adj && (!quick || "sds" %in% disp)) {
                for (i in weight.names) {
                    for (t in c("0", "1")) {
                        sds <- rep(NA_real_, NCOL(C))
                        if (any(sd.computable)) {
                            sds[sd.computable] <- col_w_sd(C[, sd.computable,drop = FALSE], weights = weights[[i]], s.weights = s.weights,
                                                           bin.vars = bin.vars[sd.computable], subset = treat == tn01[t])
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
            if (!no.adj && (!quick || "sds" %in% disp)) {
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
            if (!get_from_STATS("adj_only")[s] && (!quick || un)) {
                B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]] <- STATS[[s]]$fun(C, treat = treat, weights = NULL, 
                                                                                      std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                                      s.d.denom = if_null_then(s.d.denom.list[[1]], s.d.denom[1]),
                                                                                      abs = abs, s.weights = s.weights, bin.vars = bin.vars,
                                                                                      weighted.weights = weights[[1]], ...)
            }
            
            if (!no.adj && (!quick || s %in%  disp)) {
                for (i in weight.names) {
                    B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]] <- STATS[[s]]$fun(C, treat = treat, weights = weights[[i]],
                                                                                       std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                                       s.d.denom = if_null_then(s.d.denom.list[[i]], s.d.denom[i]),
                                                                                       abs = abs, s.weights = s.weights, bin.vars = bin.vars, ...)
                }
            }
            
            
            if (all(vapply(intersect(names(B), paste.(STATS[[s]]$bal.tab_column_prefix, c("Un", weight.names))), 
                           function(x) all(!is.finite(B[[x]])), logical(1L)))) {
                disp <- disp[disp != s] 
                thresholds[[s]] <- NULL
            }
            
            if (is_not_null(thresholds[[s]])) {
                if (!get_from_STATS("adj_only")[s] && no.adj) {
                    B[[paste.(STATS[[s]]$Threshold, "Un")]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]]), "", 
                                                                       STATS[[s]]$abs(B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                                       paste0("Not Balanced, >", round(thresholds[[s]], 3)))
                }
                else if (!no.adj) {
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
    
    B
}

samplesize <- function(treat, type, weights = NULL, subclass = NULL, s.weights = NULL,
                       method = c("matching", "weighting", "subclassification"), discarded = NULL) {
    #Computes sample size info. for unadjusted and adjusted samples.
    # method is what method the weights are to be used for. 
    # method="subclassification" is for subclass sample sizes only.
    
    if (is_null(s.weights)) s.weights <- rep(1, length(treat))
    if (is_null(discarded)) discarded <- rep(FALSE, length(treat))
    
    if (type == "bin") {
        if (length(method) == 1 && method == "subclassification") {
            if (is_null(subclass)) .err("`subclass` must be a vector of subclasses")
            
            nn <- make_df(c(levels(subclass), "Discarded", "All"), c(treat_names(treat), "Total"))
            
            nn[["All"]] <- c(vapply(treat_vals(treat), function(tn) sum(treat==tn), numeric(1L)), length(treat))
            nn[["Discarded"]] <- {
                if (any(discarded)) c(vapply(treat_vals(treat), function(tn) sum(treat[discarded]==tn), numeric(1L)), length(treat))
                else NULL
            }
            
            matched <- !is.na(subclass)
            for (k in levels(subclass)) {
                qt <- treat[matched & subclass == k]
                nn[[k]] <- c(vapply(treat_vals(treat), function(tn) sum(qt==tn), numeric(1L)), length(qt))
            }
            for (tnn in names(treat_names(treat))) {
                small.subclass <- nn[treat_names(treat)[tnn], levels(subclass)] <= 1
                if (any(small.subclass))
                    .wrn(sprintf("Not enough %s units in subclass%s %s",
                                 tnn,
                                 ngettext(sum(small.subclass), "", "es"), 
                                 word_list(levels(subclass)[small.subclass])))
            }
            attr(nn, "tag") <- "Sample sizes by subclass"
        }
        else {
            if (is_null(weights)) {
                nn <- make_df(treat_names(treat), "All")
                nn["All", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
                if (nunique.gt(s.weights, 2) || !any(s.weights == 1) || !all(s.weights %in% c(0,1))) {
                    attr(nn, "ss.type") <- c("ess")
                }
                else {
                    attr(nn, "ss.type") <- c("ss")
                }
                
            }
            else if (NCOL(weights) == 1) {
                if (method == "matching") {
                    nn <- make_df(treat_names(treat), c("All (ESS)", "All (Unweighted)", "Matched (ESS)",
                                                        "Matched (Unweighted)", "Unmatched", "Discarded"))
                    nn["All (ESS)", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
                    nn["All (Unweighted)", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & s.weights > 0), numeric(1L))
                    nn["Matched (ESS)", ] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, 1]*s.weights[treat==tn]), numeric(1L))
                    nn["Matched (Unweighted)", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & weights[,1] > 0 & s.weights > 0), numeric(1L))
                    nn["Unmatched", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & weights[,1]==0 & !discarded), numeric(1L))
                    nn["Discarded", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & discarded), numeric(1L))
                    
                    attr(nn, "ss.type") <- rep("ss", NROW(nn))
                    
                    if (!any(discarded)) {
                        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                        nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
                    }
                }
                else if (method == "weighting") {
                    nn <- make_df(treat_names(treat), c("Unadjusted", "Adjusted", "Discarded"))
                    nn["Unadjusted", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
                    nn["Adjusted", ] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, 1]*s.weights[treat==tn]), numeric(1L))
                    nn["Discarded", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & discarded), numeric(1L))
                    attr(nn, "ss.type") <- c("ss", "ess", "ss")
                    
                    if (!any(discarded)) {
                        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                        nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
                    }
                }
            }
            else {
                nn <- make_df(treat_names(treat), c("All", names(weights)))
                nn["All", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
                for (i in seq_len(NCOL(weights))) {
                    nn[1+i,] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, i] * s.weights[treat==tn]), numeric(1L))
                }
                attr(nn, "ss.type") <- c("ss", rep("ess", length(method)))
            }
            
            attr(nn, "tag") <- {
                if (length(attr(nn, "ss.type")) > 1 && all(attr(nn, "ss.type")[-1] == "ess")) {
                    "Effective sample sizes"
                }
                else "Sample sizes"
            }
        }
    }
    else if (type == "cont") {
        if (length(method) == 1 && method == "subclassification") {
            if (is_null(subclass)) .err("`subclass` must be a vector of subclasses")
            
            nn <- make_df(c(levels(subclass), "All"), c("Total"))
            
            nn[, "All"] <- length(treat)
            
            matched <- !is.na(subclass)
            for (k in levels(subclass)) {
                nn[[k]] <- sum(matched & subclass == k)
            }
            small.subclass <- nn[, levels(subclass)] <= 1
            if (any(small.subclass))
                .wrn(sprintf("not enough units in subclass%s %s",
                             ngettext(sum(small.subclass), "", "es"), 
                             word_list(levels(subclass)[small.subclass])))
            attr(nn, "tag") <- "Sample sizes by subclass"
        }
        else {
            if (is_null(weights)) {
                nn <- make_df("Total", "All")
                nn["All", ] <- ESS(s.weights)
                if (nunique.gt(s.weights, 2) || !any(s.weights == 1) || !all(s.weights %in% c(0,1))) {
                    attr(nn, "ss.type") <- c("ess")
                }
                else {
                    attr(nn, "ss.type") <- c("ss")
                }
                
            }
            else if (NCOL(weights) == 1) {
                if (method == "matching") {
                    nn <- make_df("Total", c("All (ESS)", "All (Unweighted)", "Matched (ESS)", "Matched (Unweighted)", "Unmatched", "Discarded"))
                    nn["All (ESS)", ] <- ESS(s.weights)
                    nn["All (Unweighted)", ] <- sum(s.weights > 0)
                    nn["Matched (ESS)", ] <- ESS(weights[, 1] * s.weights)
                    nn["Matched (Unweighted)", ] <- sum(weights[,1] > 0 & s.weights > 0 & !discarded)
                    nn["Unmatched", ] <- sum(weights[,1] == 0 & !discarded)
                    nn["Discarded", ] <- sum(discarded)
                    
                    attr(nn, "ss.type") <- rep("ss", NROW(nn))
                    
                    if (!any(discarded)) {
                        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                        nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
                    }
                }
                else if (method == "weighting") {
                    nn <- make_df("Total", c("Unadjusted", "Adjusted", "Discarded"))
                    nn["Unadjusted", ] <- ESS(s.weights)
                    nn["Adjusted", ] <- ESS(weights[!discarded, 1]*s.weights[!discarded])
                    nn["Discarded", ] <- sum(discarded)
                    attr(nn, "ss.type") <- c("ss", "ess", "ss")
                    
                    if (!any(discarded)) {
                        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
                        nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
                    }
                }
            }
            else {
                nn <- make_df("Total", c("All", names(weights)))
                nn["All", ] <- ESS(s.weights)
                for (i in seq_len(NCOL(weights))) {
                    if (method[i] == "matching") {
                        nn[1+i,] <- ESS(weights[!discarded, i])
                    }
                    else if (method[i] == "weighting") {
                        nn[1+i,] <- ESS(weights[!discarded, i]*s.weights[!discarded])
                    }
                    
                }
                attr(nn, "ss.type") <- c("ss", rep("ess", length(method)))
                
            }
            
            attr(nn, "tag") <- {
                if (length(attr(nn, "ss.type")) > 1 && all(attr(nn, "ss.type")[-1] == "ess")) {
                    "Effective sample sizes"
                }
                else "Sample sizes"
            }
        }
    }
    
    nn
}

balance.summary <- function(bal.tab.list, agg.funs, include.times = FALSE) {
    type <- attr(bal.tab.list[[1]], "print.options")[["type"]]
    disp <- attr(bal.tab.list[[1]], "print.options")[["disp"]]
    compute <- attr(bal.tab.list[[1]], "print.options")[["compute"]]
    thresholds <- attr(bal.tab.list[[1]], "print.options")[["thresholds"]]
    quick <- attr(bal.tab.list[[1]], "print.options")[["quick"]]
    weight.names <- if_null_then(attr(bal.tab.list[[1]], "print.options")[["weight.names"]], "Adj")
    abs <- attr(bal.tab.list[[1]], "print.options")[["abs"]]
    no.adj <- attr(bal.tab.list[[1]], "print.options")[["nweights"]] == 0
    
    balance.list <- clear_null(grab(bal.tab.list, "Balance"))
    
    Brownames <- unique(unlist(lapply(balance.list, rownames), use.names = FALSE))
    
    Agg.Funs <- firstup(if (quick) agg.funs else c("min", "mean", "max"))
    Agg.Funs.Given <- firstup(agg.funs)
    
    if (length(Agg.Funs) == 1 && Agg.Funs == "Max") abs <- TRUE
    
    Bcolnames <- c("Times", "Type", 
                   paste.(c(unlist(lapply(compute[compute %in% all_STATS(type)[!get_from_STATS("adj_only")[get_from_STATS("type") == type]]], function(s) {
                       c(paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix),
                         if (no.adj && is_not_null(thresholds[[s]]) && length(Agg.Funs.Given) == 1) STATS[[s]]$Threshold)
                   }))), "Un"),
                   unlist(lapply(weight.names, function(wn) {
                       paste.(c(unlist(lapply(compute[compute %in% all_STATS(type)], function(s) {
                           c(paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix),
                             if (!no.adj && is_not_null(thresholds[[s]]) && length(Agg.Funs.Given) == 1) STATS[[s]]$Threshold)
                       }))), wn)
                   })))
    
    B <- make_df(Bcolnames, Brownames)
    
    B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(vapply(balance.list, function(y) if (x %in% rownames(y)) y[[x, "Type"]] else NA_character_, character(1))))), use.names = FALSE)
    
    B[["Times"]] <- {
        if (include.times) vapply(Brownames, function(x) paste(seq_along(balance.list)[vapply(balance.list, function(y) x %in% rownames(y), logical(1L))], collapse = ", "), character(1))[Brownames]
        else NULL
    }
    
    for (Agg.Fun in Agg.Funs) {
        for (s in compute[compute %in% all_STATS(type)]) {
            abs0 <- function(x) {if (is_null(x)) NA_real_ else if (abs) STATS[[s]]$abs(x) else (x)}
            agg <- function(x, ...) {
                if (!any(is.finite(x))) NA_real_
                else if (s == "variance.ratios" && tolower(Agg.Fun) == "mean") .geam_mean(x)
                else if (tolower(Agg.Fun) == "rms") sqrt(mean_fast(STATS[[s]]$abs(x)^2, TRUE))
                else get(tolower(Agg.Fun))(x, ...)
            }
            for (sample in c("Un", weight.names)) {
                if ((sample == "Un" || !no.adj) && (sample != "Un" || !get_from_STATS("adj_only")[s])) {
                    B[[paste.(Agg.Fun, STATS[[s]]$bal.tab_column_prefix, sample)]] <- vapply(Brownames, function(x) agg(unlist(lapply(balance.list, function(y) if (x %in% rownames(y)) abs0(y[[x, paste.(STATS[[s]]$bal.tab_column_prefix, sample)]]))), na.rm = TRUE), numeric(1))
                }
            }
        }
    }
    
    if (length(Agg.Funs.Given) == 1) {
        #Assign X.Threshold values
        for (s in compute[compute %in% all_STATS(type)]) {
            if (is_not_null(thresholds[[s]])) {
                if (!get_from_STATS("adj_only")[s] && no.adj) {
                    B[[paste.(STATS[[s]]$Threshold, "Un")]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, "Un")]]), "", 
                                                                       STATS[[s]]$abs(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, "Un")]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                                       paste0("Not Balanced, >", round(thresholds[[s]], 3)))
                }
                else {
                    for (i in weight.names) {
                        B[[paste.(STATS[[s]]$Threshold, i)]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, i)]]), "", 
                                                                        STATS[[s]]$abs(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, i)]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                                        paste0("Not Balanced, >", round(thresholds[[s]], 3)))
                    }
                }
            }
            
            if (no.adj || length(weight.names) <= 1) {
                names(B)[names(B) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
            }
        }
    }
    
    B
}


#base.bal.tab.imp
samplesize.across.imps <- function(obs.list) {
    obs.list <- clear_null(obs.list)
    
    obs <- Reduce("+", obs.list)/length(obs.list)
    attr(obs, "tag") <- sprintf("Average %s across imputations",
                                tolower(attr(obs.list[[1]], "tag")))
    obs
}

#base.bal.tab.multi
samplesize.multi <- function(bal.tab.multi.list, treat_names, focal = NULL) {
    which <- {
        if (is_null(focal)) treat_names
        else c(treat_names[treat_names != focal], focal)
    }
    
    bal.tab.multi.list <- clear_null(bal.tab.multi.list)
    obs <- do.call("cbind", unname(grab(bal.tab.multi.list, "Observations")))[, which]
    attr(obs, "tag") <- attr(bal.tab.multi.list[[1]][["Observations"]], "tag")
    attr(obs, "ss.type") <- attr(bal.tab.multi.list[[1]][["Observations"]], "ss.type")
    
    obs
}

#base.bal.tab.msm
samplesize.msm <- function(bal.tab.msm.list) {
    obs <- do.call("cbind", grab(bal.tab.msm.list, "Observations"))
    attr(obs, "tag") <- attr(bal.tab.msm.list[[1]][["Observations"]], "tag")
    attr(obs, "ss.type") <- attr(bal.tab.msm.list[[1]][["Observations"]], "ss.type")
    
    obs
}

#base.bal.tab.cluster
samplesize.across.clusters <- function(obs.list) {
    obs.list <- clear_null(obs.list)
    obs <- Reduce("+", obs.list)
    attr(obs, "tag") <- sprintf("Total %s across clusters",
                                tolower(attr(obs.list[[1]], "tag")))
    
    obs
}

#base.bal.tab.subclass
balance.table.subclass <- function(C, type, weights = NULL, treat, subclass,
                                   continuous, binary, s.d.denom, 
                                   thresholds = list(), un = FALSE, disp = NULL, stats = NULL, 
                                   s.weights = rep(1, length(treat)), abs = FALSE, 
                                   var_types = NULL, quick = TRUE, ...) {
    #Creates list SB of balance tables for each subclass
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    disp <- unique(c(disp, all_STATS(type)[all_STATS(type) %in% stats]))
    compute <- if (quick) disp else c("means", "sds", all_STATS(type))
    
    #B=Balance frame
    Bnames <- c("Type", 
                expand.grid_string(c(if (type == "bin") expand.grid_string(c("M"["means" %in% compute], "SD"["sds" %in% compute]), c("0", "1"), collapse = ".")
                                     else if (type == "cont") c("M"["means" %in% compute], "SD"["sds" %in% compute]), 
                                     unlist(lapply(compute[compute %in% all_STATS(type)], function(s) {
                                         c(STATS[[s]]$bal.tab_column_prefix,
                                           if (is_not_null(thresholds[[s]])) STATS[[s]]$Threshold)
                                     }))
                ),
                c("Adj"), collapse = "."))
    B <- make_df(Bnames, colnames(C))
    
    #Set var type (binary/continuous)
    B[["Type"]] <- if_null_then(var_types, .get_types(C))
    bin.vars <- B[["Type"]] == "Binary"
    
    SB <- setNames(lapply(levels(subclass), function(i) B), levels(subclass))
    
    binary <- match_arg(binary, c("raw", "std"))
    sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else !bin.vars
    
    if (type == "bin") 
        subclass_w_empty <- vapply(levels(subclass), function(i) {
            any(vapply(treat_vals(treat), function(t) !any(treat == t & subclass == i), logical(1L)))
        }, logical(1L))
    else {
        subclass_w_empty <- vapply(levels(subclass), function(i) !any(subclass == i), logical(1L))
    }
    
    for (i in levels(subclass)) {
        
        in.subclass <- !is.na(subclass) & subclass==i
        
        #Means for each group
        if ("means" %in% compute) {
            if (type == "bin") {
                tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
                for (t in c("0", "1")) {
                    SB[[i]][[paste.("M", t, "Adj")]] <- col_w_mean(C, subset = treat==tn01[t] & in.subclass, s.weights = s.weights)
                }
            }
            else if (type == "cont") {
                SB[[i]][["M.Adj"]] <- col_w_mean(C, subset = in.subclass, s.weights = s.weights)
            }
        }
        
        #SDs for each group
        if ("sds" %in% compute) {
            if (type == "bin") {
                for (t in c("0", "1")) {
                    sds <- rep(NA_real_, NCOL(C))
                    sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == tn01[t] & in.subclass, s.weights = s.weights)
                    SB[[i]][[paste.("SD", t, "Adj")]] <- sds
                }
            }
            else if (type == "cont") {
                sds <- rep(NA_real_, NCOL(C))
                sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == in.subclass, s.weights = s.weights)
                SB[[i]][["SD.Adj"]] <- sds
            }
        }
        
        for (s in all_STATS(type)) {
            if (s %in% compute && !subclass_w_empty[i]) {
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
    
    SB
}

# !!! NEEDS TO BE UPDATED !!!
balance.table.across.subclass.cont <- function(balance.table, balance.table.subclass.list, subclass.obs, r.threshold = NULL) {
    #NEEDS TO BE UPDATED
    
    B.A <- balance.table.subclass.list[[1]][names(balance.table.subclass.list[[1]]) %in% c("M.Adj", "SD.Adj", "Corr.Adj")]
    
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
    B.A.df <- cbind(balance.table[c("Type", "M.Un", "SD.Un", "Corr.Un", "R.Threshold.Un")], 
                    B.A, R.Threshold = NA_character_)
    if (is_not_null(r.threshold)) {
        B.A.df[["R.Threshold"]] <- ifelse(B.A.df[["Type"]]=="Distance", "", paste0(ifelse(is.finite(B.A.df[["Corr.Adj"]]) & abs_(B.A.df[["Corr.Adj"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold))
    }
    
    B.A.df
}

#Misc
`%+%` <- function(...) {
    if (is.atomic(..1) && is.atomic(..2)) crayon::`%+%`(as.character(..1), as.character(..2))
    else ggplot2::`%+%`(...)
}

check_arg_lengths <- function(...) {
    dots_names <- vapply(match.call(expand.dots = FALSE)$..., deparse1,
                         character(1L))
    lengths <- setNames(vapply(list(...), len, integer(1L)),
                        dots_names)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        .err(sprintf("%s must have the same number of units",
                     word_list(dots_names[supplied], quotes = "`")))
    }
}

intapprox <- function(f, from, to, steps, method = "midpoint") {
    method <- match_arg(method, c("midpoint", "trapezoidal", "simpsons"))
    if (method == "midpoint") {
        seg <- seq(from, to, length = steps)
        delta <- seg[2] - seg[1]
        mids <- .5 * (seg[-1] + seg[-steps])
        s <- sum(f(mids)) * delta
    }
    else if (method == "trapezoidal") {
        seg <- seq(from, to, length = steps)
        delta <- seg[2] - seg[1]
        s <- (f(from) + 2 * sum(f(seg[-c(1, steps)])) + f(to)) * delta / 2
    }
    else if (method == "simpsons") {
        s <- (2 * intapprox(f, from, to, steps, "midpoint") + intapprox(f, from, to, steps, "trapezoidal")) / 3
    }
    
    s
}
