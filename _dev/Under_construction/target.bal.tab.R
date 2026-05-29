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
    
    if (is_not_null(X$s.d.denom) && X$s.d.denom %in% c("treated", "control")) {
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
        warning("Target balance assessment will only occur for the baseline covariates.", call. = FALSE)
        X$treat <- X$treat.list[[1]]
        X$covs <- X$covs.list[[1]]
        X$distance <- X$distance.list[[1]]
        X$addl <- X$addl.list[[1]]
    }
    
    out <- do.call("base.bal.tab.target", c(X, args),
                   quote = TRUE)
    return(out)
}

base.bal.tab.target <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), which.treat = NA, target.summary = getOption("cobalt_target.summary", TRUE), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
    #Preparations
    args <- list(...)
    
    if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
    if (is_not_null(v.threshold)) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (is_not_null(ks.threshold)) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (is_null(weights) && is_null(subclass)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (is_null(s.weights)) {
        s.weights <- rep(1, length(treat))
    }
    if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
    
    #Create new Target group
    target.name <- "Target"
    n <- length(treat)
    
    if (isTRUE(get.treat.type(treat) == "continuous") || (is.numeric(treat) && !is_binary(treat))) {#if continuous treatment
        covs <- data.frame(treat = treat, covs)
        treat <- factor(rep(c("All", target.name), each = n))
        target.summary <- FALSE
        which.treat <- NULL
        needs.summary <- FALSE
        treat_names <- unique.treat <- "All"
    }
    else {
        if (is.factor(treat) || is.character(treat)) {
            if (is.factor(treat)) treat_names <- unique.treat <- levels(treat)
            else treat_names <- unique.treat <- unique(treat, nmax = n - 1)
        }
        else {
            treat_names <- c("Control", "Treated")
            unique.treat <- sort(unique(treat, nmax = 2))
        }
        names(treat_names) <- unique.treat
        
        treat <- factor(c(treat_names[as.character(treat)], rep(target.name, n)))
        needs.summary <- TRUE
    }
    
    covs <- rbind(covs, covs)
    if (is_not_null(weights)) weights <- rbind(weights, as.data.frame(array(1, dim = dim(weights), 
                                                                            dimnames = dimnames(weights))))
    distance <- rbind(distance, distance)
    addl <- rbind(addl, addl)
    s.weights <- c(s.weights, s.weights)
    if (is_not_null(discarded)) discarded <- c(discarded, rep(FALSE, length(discarded)))
    s.d.denom <- "treated"
    
    treat.target.combinations <- lapply(treat_names, function(x) c(x, target.name))
    
    if (is_not_null(cluster)) {
        stop("Clusters are not yet supported with target balance assessment.", call. = FALSE)
    }
    else if (is_not_null(subclass)) {
        stop("Subclassification is not yet supported with target balance assessment.", call. = FALSE)
    }
    else {
        #Setup output object
        out.names <- c("Target.Balance", 
                       "Balance.Across.Treatments", 
                       "Observations", 
                       "call")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        
        if (any(treat_names == "Target")) stop ("\"Target\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
        args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
        balance.tables <- lapply(treat.target.combinations, function(t) do.call("base.bal.tab.binary", c(list(weights = weights[treat %in% t, , drop = FALSE], treat = factor(treat[treat %in% t], t), distance = distance[treat %in% t, , drop = FALSE], subclass = subclass[treat %in% t], covs = covs[treat %in% t, , drop = FALSE], call = NULL, int = int, poly = poly, addl = addl[treat %in% t, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster[treat %in% t], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[treat %in% t], discarded = discarded[treat %in% t], quick = quick), args), quote = TRUE))
        
        for (i in seq_along(balance.tables)) {
            names(balance.tables)[i] <- paste(treat.target.combinations[[i]], collapse = " vs. ")
            balance.tables[[i]][["Observations"]][[2]] <- NULL
        }
        
        out[["Target.Balance"]] <- balance.tables
        
        out[["Observations"]] <- samplesize.target(balance.tables, treat_names, target.name) 
        
        if (needs.summary && (target.summary || !quick)) {
            out[["Balance.Across.Treatments"]] <- balance.table.target.summary(balance.tables, 
                                                                               weight.names = names(weights),
                                                                               m.threshold = m.threshold,
                                                                               v.threshold = v.threshold,
                                                                               ks.threshold = ks.threshold,
                                                                               no.adj = no.adj,
                                                                               quick = quick,
                                                                               types = NULL)
        }
        
        out[["call"]] <- call
        
        attr(out, "print.options") <- list(m.threshold=m.threshold,
                                           v.threshold=v.threshold,
                                           ks.threshold=ks.threshold,
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.adj=!no.adj, 
                                           which.cluster=which.cluster,
                                           cluster.summary=cluster.summary,
                                           cluster.fun = cluster.fun,
                                           abs = abs,
                                           continuous = continuous,
                                           binary = binary,
                                           quick = quick,
                                           disp.means=disp.means, 
                                           disp.sds = disp.sds,
                                           disp.v.ratio=disp.v.ratio, 
                                           disp.ks=disp.ks,
                                           disp.bal.tab = disp.bal.tab,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights),
                                           treat_names = treat_names,
                                           target.name = target.name,
                                           which.treat = which.treat,
                                           target.summary = target.summary,
                                           co.names = attr(out[["Target.Balance"]][[1]], "print.options")[["co.names"]])
        
        class(out) <- c("bal.tab.target", "bal.tab")
    }
    return(out)
    
}

print.bal.tab.target <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.sds = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", which.treat, target.summary = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval(.call))
    }
    
    args <- c(as.list(environment()), list(...))[-1]
    
    call <- x$call
    t.balance <- x[["Target.Balance"]]
    t.balance.summary <- x[["Balance.Across.Treatments"]]
    nn <- x$Observations
    p.ops <- attr(x, "print.options")
    
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("un cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp.means, "as.is")) {
        if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.means == FALSE && disp.means == TRUE) {
            warning("disp.means cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.means <- disp.means
    }
    if (!identical(disp.sds, "as.is")) {
        if (!is.logical(disp.sds)) stop("disp.sds must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.sds == FALSE && disp.sds == TRUE) {
            warning("disp.sds cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.sds <- disp.sds
    }
    if (!identical(disp.v.ratio, "as.is")) {
        if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.v.ratio == FALSE && disp.v.ratio == TRUE) {
            warning("disp.v.ratio cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.v.ratio <- disp.v.ratio
    }
    if (!identical(disp.ks, "as.is")) {
        if (!is.logical(disp.ks)) stop("disp.ks must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.ks == FALSE && disp.ks == TRUE) {
            warning("disp.ks cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.ks <- disp.ks
    }
    if (!identical(target.summary, "as.is")) {
        if (!is.logical(target.summary)) stop("target.summary must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$target.summary == FALSE && target.summary == TRUE) {
            warning("target.summary cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$target.summary <- target.summary
    }
    if (!identical(disp.m.threshold, "as.is")) {
        if (!is.logical(disp.m.threshold)) stop("disp.m.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (is_null(p.ops$disp.v.ratio) || !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (is_null(p.ops$disp.ks) || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (is_not_null(t.balance.summary)) {
        if (p.ops$disp.bal.tab) {
            if (!identical(imbalanced.only, "as.is")) {
                if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
                p.ops$imbalanced.only <- imbalanced.only
            }
            if (p.ops$imbalanced.only) {
                if (all(sapply(c(p.ops$m.threshold, 
                                 p.ops$v.threshold, 
                                 p.ops$ks.threshold, 
                                 p.ops$r.threshold), is_null))) {
                    warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                    p.ops$imbalanced.only <- FALSE
                }
            }
        }
        else p.ops$imbalanced.only <- FALSE
        
        if (p.ops$imbalanced.only) {
            keep.row <- rowSums(apply(t.balance.summary[grepl(".Threshold", names(t.balance.summary), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        }
        else keep.row <- rep(TRUE, nrow(t.balance.summary))
    }
    
    if (!missing(which.treat)) {
        if (paste(deparse(substitute(which.treat)), collapse = "") == ".none") which.treat <- NA
        else if (paste(deparse(substitute(which.treat)), collapse = "") == ".all") which.treat <- NULL
        if (!identical(which.treat, "as.is")) {
            p.ops$which.treat <- which.treat
        }
    }
    
    #Checks and Adjustments
    if (is_null(p.ops$which.treat)) 
        which.treat <- p.ops$treat_names
    else if (anyNA(p.ops$which.treat)) {
        which.treat <- character(0)
        if (!p.ops$target.summary) message("No treatments will be displayed; displaying the summary across treatments.")
        p.ops$target.summary <- TRUE
    }
    else if (is.numeric(p.ops$which.treat)) {
        which.treat <- p.ops$treat_names[seq_along(p.ops$treat_names) %in% p.ops$which.treat]
        if (is_null(which.treat)) {
            warning("No numbers in which.treat correspond to treatment values. No treatments will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
    }
    else if (is.character(p.ops$which.treat)) {
        which.treat <- p.ops$treat_names[p.ops$treat_names %in% p.ops$which.treat]
        if (is_null(which.treat)) {
            warning("No names in which.treat correspond to treatment values. No treatments will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
    }
    else {
        warning("The argument to which.treat must be NA, NULL, or a vector of treatment names or indices. No treatments will be displayed.", call. = FALSE)
        which.treat <- character(0)
        if (!p.ops$target.summary) message("No treatments will be displayed; displaying the summary across treatments.")
        p.ops$target.summary <- TRUE
    }
    
    if (p.ops$target.summary && is_null(t.balance.summary)) {
        warning("No summary across treatments was produced. This can occur if target.summary is FALSE and quick is TRUE.", call. = FALSE)
    }
    
    if (is_null(which.treat)) {
        disp.treat.pairs <- character(0)
    }
    else {
        disp.treat.pairs <- names(t.balance)[sapply(names(t.balance), function(x) any(attr(t.balance[[x]], "print.options")$treat_names %in% which.treat))]
    }
    
    #Printing output
    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (is_not_null(disp.treat.pairs)) {
        headings <- setNames(character(length(disp.treat.pairs)), disp.treat.pairs)
        cat(underline("Target balance by treatment group") %+% "\n")
        for (i in disp.treat.pairs) {
            headings[i] <- "\n - - - " %+% italic(attr(t.balance[[i]], "print.options")$treat_names[1] %+% " (0) vs. " %+%
                                                      p.ops$target.name %+% " (1)") %+% " - - - \n"
            cat(headings[i])
            do.call(print, c(list(t.balance[[i]]), args))
        }
        cat(paste0(paste(rep(" -", round(max(nchar(headings))/2)), collapse = ""), " \n"))
        cat("\n")
    }
    
    if (isTRUE(as.logical(p.ops$target.summary)) && is_not_null(t.balance.summary)) {
        s.keep <- as.logical(c(TRUE, 
                               p.ops$un,
                               p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$m.threshold),
                               p.ops$un && p.ops$disp.v.ratio, 
                               p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                               p.ops$un && p.ops$disp.ks, 
                               p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$ks.threshold),
                               rep(c(p.ops$disp.adj, 
                                     p.ops$disp.adj && is_not_null(p.ops$m.threshold), 
                                     p.ops$disp.adj && p.ops$disp.v.ratio, 
                                     p.ops$disp.adj && is_not_null(p.ops$v.threshold), 
                                     p.ops$disp.adj && p.ops$disp.ks, 
                                     p.ops$disp.adj && is_not_null(p.ops$ks.threshold)), p.ops$nweights + !p.ops$disp.adj)))
        
        if (p.ops$disp.bal.tab) {
            cat(underline("Target balance summary across all treatment groups") %+% "\n")
            if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
            else print.data.frame_(round_df_char(t.balance.summary[keep.row, s.keep], digits))
            cat("\n")
        }
        
        if (is_not_null(nn)) {
            tag <- attr(x$Observations, "tag")
            ss.type <- attr(x$Observations, "ss.type")
            for (i in rownames(x$Observations)) {
                if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i,]
            }
            if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(x$Observations)) && 
                all(check_if_zero(x$Observations["Matched (ESS)",] - x$Observations["Matched (Unweighted)",]))) {
                x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)", , drop = FALSE]
                rownames(x$Observations)[rownames(x$Observations) == "Matched (ESS)"] <- "Matched"
            }
            cat(underline(tag) %+% "\n")
            print.warning <- FALSE
            if (length(ss.type) > 1 && nunique.gt(ss.type[-1], 1)) {
                ess <- ifelse(ss.type == "ess", "*", "")
                x$Observations <- setNames(cbind(x$Observations, ess), c(names(x$Observations), ""))
                print.warning <- TRUE
            }
            print.data.frame_(round_df_char(x$Observations, digits = max(0, digits-1)))
            if (print.warning) cat(italic("* indicates effective sample size"))
        }
    }
    
    invisible(x)
    
}

#base.bal.tab.target
samplesize.target <- function(bal.tab.target.list, treat_names, target.name) {
    which <- treat_names[treat_names != target.name]
    obs <- do.call("cbind", unname(lapply(bal.tab.target.list, function(x) x[["Observations"]])))[, which]
    attr(obs, "tag") <- attr(bal.tab.target.list[[1]][["Observations"]], "tag")
    attr(obs, "ss.type") <- attr(bal.tab.target.list[[1]][["Observations"]], "ss.type")
    return(obs)
}