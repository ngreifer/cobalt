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
match.strata2weights <- function(match.strata, treat, covs = NULL) {
    #Process match.strata into weights (similar to weight.subclass from MatchIt)
    if (is_null(covs)) names(treat) <- seq_along(treat)
    else names(treat) <- row.names(covs)
    matched <- !is.na(match.strata); unmatched <- !matched
    treat.matched <- treat[matched]
    match.strata.matched <- match.strata[matched]
    
    labels <- names(treat.matched)
    tlabels <- labels[treat.matched == 1]
    clabels <- labels[treat.matched == 0]
    
    weights.matched <- rep(0, length(treat.matched))
    names(weights.matched) <- labels
    weights.matched[tlabels] <- 1
    
    for (j in unique(match.strata.matched)){
        qn0 <- sum(treat.matched==0 & match.strata.matched==j)
        qn1 <- sum(treat.matched==1 & match.strata.matched==j)
        weights.matched[treat.matched==0 & match.strata.matched==j] <- qn1/qn0
    }
    if (all(check_if_zero(weights.matched[clabels]))) #all control weights are 0
        weights.matched[clabels] <- rep(0, length(weights.matched[clabels]))
    else {
        ## Number of C units that were matched to at least 1 T
        num.cs <- sum(!check_if_zero(weights.matched[clabels]))
        weights.matched[clabels] <- weights.matched[clabels]*num.cs/sum(weights.matched[clabels])
    }
    
    if (any(unmatched)) {
        weights.unmatched <- rep(0, sum(unmatched))
        names(weights.unmatched) <- names(treat[unmatched])
        weights <- c(weights.matched, weights.unmatched)[names(treat)]
    }
    else {
        weights <- weights.matched
    }
    
    if (all(check_if_zero(weights))) 
        stop("No units were matched", call. = FALSE)
    else if (all(check_if_zero(weights[tlabels])))
        stop("No treated units were matched", call. = FALSE)
    else if (all(check_if_zero(weights[clabels])))
        stop("No control units were matched", call. = FALSE)
    return(weights)
}
weights.same.as.strata <- function(weights, match.strata, treat) {
    weights == match.strata2weights(match.strata, treat)
}
use.tc.fd <- function(formula = NULL, data = NULL, treat = NULL, covs = NULL, needs.treat = TRUE, needs.covs = TRUE) {
    if (is_not_null(formula) && class(formula) == "formula") {
        D <- NULL
        if (is_not_null(data)) D <- data
        if (is_not_null(covs)) if (is_not_null(D)) D <- cbind(D, covs) else D <- covs
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
        if (is_not_null(val.df)) { if (sum(is.na(val.df)) > 0) {
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
                if (is_not_null(val.df)) { if (sum(is.na(val.df)) > 0) {
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
get.X.class <- function(X) {
    if (is_not_null(X[["imp"]])) {
        if (is_not_null(X[["treat.list"]])) stop("Multiply imputed data is not yet supported with longitudinal treatments.", call. = FALSE)
        else if (is_(X[["treat"]], c("factor", "character")) && !is_binary(X[["treat"]])) stop("Multiply imputed data is not yet supported with multinomial treatments.", call. = FALSE)
        else X.class <- "imp"
    }
    else if (is_not_null(X[["treat.list"]])) X.class <- "msm"
    else if (is_binary(X[["treat"]])) X.class <- "binary"
    else if (is_(X[["treat"]], c("factor", "character"))) X.class <- "multi"
    else if (is.numeric(X[["treat"]])) X.class <- "cont"
    else probably.a.bug()
    
    return(X.class)
}
subset_X <- function(X, subset) {
    n <- length(subset)
    subset_X_internal <- function(x, subset) {
        if (is_not_null(x)) {
            if (is.factor(x) && length(x) == n) factor(x[subset])
            else if (is.atomic(x) && length(x) == n) x[subset]
            else if ((is.matrix(x) || is.data.frame(x))  && nrow(x) == n) x[subset, , drop = FALSE]
            else if (is.list(x)) lapply(x, subset_X_internal, subset = subset)
            else x
        }
        else x
    }
    lapply(X, subset_X_internal, subset)
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
        stop(paste0("The variable", s[1], " ", word.list(names(C)[covs.with.inf], quotes = TRUE), 
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
        
        if (nlevels(cluster) > 0 && qr(matrix(c(C[[i]], as.numeric(cluster)), ncol = 2))$rank == 1) {
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
                    split.points <- c(nchar(i), nchar(i)+nchar(A[["factor_sep"]]))
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
                all.possible.names <- outer(colnames(C), colnames(C), paste, sep = paste0(rep(A[["int_sep"]], nsep), collapse = ""))
                if (!any(colnames(C) %in% all.possible.names)) break
                else nsep <- nsep + 1
            }
            
            if (poly < int) poly <- int
            
            int <- TRUE
        }
        
        new <- int.poly.f(C, int = int, poly = poly, center = A[["center"]], sep = rep(A[["int_sep"]], nsep), co.names = co.names)
        C <- cbind(C, new)
        co.names <- c(co.names, attr(new, "co.names"))
    }
    
    #Add missingness indicators
    vars.w.missing <- vars.w.missing[vars.w.missing$placed.after %in% colnames(C) & vars.w.missing$has.missing, , drop = FALSE]
    if (NROW(vars.w.missing) > 0) {
        missing.ind <- apply(C[,colnames(C) %in% vars.w.missing$placed.after, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
        #colnames(missing.ind) <- paste0(rownames(vars.w.missing), ":<NA>")
        colnames(missing.ind) <- rownames(vars.w.missing)
        #missing.ind <- remove.perfect.col(missing.ind) 
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
            new <- int.poly.f(missing.ind, int = TRUE, poly = 1, sep = rep(A[["int_sep"]], nsep), co.names = miss.co.names)
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
        C <- cbind(distance, C, row.names = NULL)
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
check_if_zero_weights <- function(weights.df, treat, unique.treat = NULL, treat.type = "cat") {
    #Checks if all weights are zero in each treat group for each set of weights
    if (treat.type == "cat") {
        if (is_null(unique.treat)) unique.treat <- unique(treat)
        w.t.mat <- expand.grid(colnames(weights.df), unique.treat, stringsAsFactors = FALSE)
        if (NROW(w.t.mat) > 0) {
            problems <- vapply(seq_len(NROW(w.t.mat)), function(x) all(check_if_zero(weights.df[treat == w.t.mat[x,2], w.t.mat[x, 1]])), logical(1L))
            if (any(problems)) {
                prob.w.t.mat <- droplevels(w.t.mat[problems,])
                if (NCOL(weights.df) == 1) {
                    error <- paste0("All weights are zero when ", word.list(paste("treat =", prob.w.t.mat[, 2]), "or"), ".")
                }
                else {
                    errors <- setNames(character(nlevels(prob.w.t.mat[,1])), levels(prob.w.t.mat[,1]))
                    
                    for (i in levels(prob.w.t.mat[,1])) {
                        errors[i] <- paste0("\"", i, "\" weights are zero when ", word.list(paste("treat =", prob.w.t.mat[prob.w.t.mat[,1] == i, 2]), "or"))
                    }
                    errors <- paste(c("All", rep("all", length(errors)-1)), errors)
                    error <- paste0(word.list(errors, "and"), ".")
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
                    error <- paste0(word.list(errors, "and"), ".")
                }
                stop(error, call. = FALSE)
            }
        }
    }
    else stop("treat.type must be either \"cat\" or \"cont\".")
    
}
.col.std.diff <- function(mat, treat, weights, subclass = NULL, which.sub = NULL, x.types, continuous, binary, s.d.denom, no.weights = FALSE, s.weights = rep(1, length(treat)), pooled.sds = NULL) {
    if (no.weights) weights <- rep(1, NROW(mat))
    w <- weights*s.weights
    sw <- s.weights
    
    #Check continuous and binary
    if (missing(continuous) || is_null(continuous)) {
        continuous <- match_arg(getOption("cobalt_continuous", "std"), c("std", "raw"))
    }
    else continuous <- match_arg(continuous, c("std", "raw"))
    if (missing(binary) || is_null(binary)) {
        binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
    }
    else binary <- match_arg(binary, c("raw", "std"))
    
    no.sub <- is_null(which.sub)
    if (no.sub) ss <- sw > 0
    else {
        ss <- (!is.na(subclass) & subclass == which.sub & sw > 0)
        
        if (sum(treat==0 & ss) == 0) {
            warning(paste0("There are no control units in subclass ", which.sub, "."), call. = FALSE)
            return(rep(NA_real_, NCOL(mat)))
        }
        if (sum(treat==1 & ss) == 0) {
            warning(paste0("There are no treated units in subclass ", which.sub, "."), call. = FALSE)
            return(rep(NA_real_, NCOL(mat)))
        }
    }
    
    diffs <- col.w.m(mat[treat == 1 & ss, , drop = FALSE], w[treat == 1 & ss]) - 
        col.w.m(mat[treat == 0 & ss, , drop = FALSE], w[treat == 0 & ss])
    diffs[check_if_zero(diffs)] <- 0
    denoms <- rep(1, NCOL(mat))
    denoms.to.std <- ifelse(x.types == "Binary", binary == "std", continuous == "std")
    
    if (any(denoms.to.std)) {
        if (s.d.denom == "control") {
            denoms[denoms.to.std] <- sqrt(col.w.v(mat[treat == 0 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 0 & ss]))
        }
        else if (s.d.denom == "treated") {
            denoms[denoms.to.std] <- sqrt(col.w.v(mat[treat == 1 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 1 & ss]))
        }
        else if (s.d.denom == "pooled") {
            if (is_not_null(pooled.sds)) {
                denoms[denoms.to.std] <- pooled.sds[denoms.to.std]
            }
            else {
                denoms[denoms.to.std] <-  sqrt(.5*(col.w.v(mat[treat == 0 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 0 & ss]) +
                                                       col.w.v(mat[treat == 1 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 1 & ss])))
            }
        }
    }
    
    std.diffs <- ifelse(check_if_zero(diffs), 0, diffs/denoms)
    if (any(!is.finite(std.diffs))) {
        warning("Some standardized mean differences were not finite. This can result from no variation in one of the treatment groups.", call. = FALSE)
        std.diffs[!is.finite(std.diffs)] <- NA_real_
    }
    
    return(std.diffs)
}
std.diffs <- function(m0, s0, m1, s1, x.types, continuous, binary, s.d.denom, pooled.sds = NULL) {
    #Check continuous and binary
    if (missing(continuous) || is_null(continuous)) {
        continuous <- match_arg(getOption("cobalt_continuous", "std"), c("std", "raw"))
    }
    else continuous <- match_arg(continuous, c("std", "raw"))
    if (missing(binary) || is_null(binary)) {
        binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
    }
    else binary <- match_arg(binary, c("raw", "std"))
    
    diffs <- m1 - m0
    
    diffs[check_if_zero(diffs)] <- 0
    denoms <- rep(1, length(diffs))
    denoms.to.std <- ifelse(x.types == "Binary", binary == "std", continuous == "std")
    
    if (any(denoms.to.std)) {
        if (s.d.denom == "control") {
            denoms[denoms.to.std] <- s0[denoms.to.std]
        }
        else if (s.d.denom == "treated") {
            denoms[denoms.to.std] <- s1[denoms.to.std]
        }
        else if (s.d.denom == "pooled") {
            if (is_not_null(pooled.sds)) {
                denoms[denoms.to.std] <- pooled.sds[denoms.to.std]
            }
            else {
                denoms[denoms.to.std] <-  sqrt(.5*(s0[denoms.to.std]^2 + s1[denoms.to.std]^2))
            }
        }
    }
    
    std.diffs <- ifelse(check_if_zero(diffs), 0, diffs/denoms)
    if (any(!is.finite(std.diffs))) {
        warning("Some standardized mean differences were not finite. This can result from no variation in one of the treatment groups.", call. = FALSE)
        std.diffs[!is.finite(std.diffs)] <- NA_real_
    }
    
    return(std.diffs)
}
col.ks <- function(mat, treat, weights, x.types, no.weights = FALSE) {
    ks <- rep(NA_integer_, NCOL(mat))
    if (no.weights) weights <- rep(1, NROW(mat))
    weights[treat == 1] <- weights[treat==1]/sum(weights[treat==1])
    weights[treat == 0] <- -weights[treat==0]/sum(weights[treat==0])
    non.binary <- x.types != "Binary"
    ks[non.binary] <- apply(mat[, non.binary, drop = FALSE], 2, function(x_) {
        x <- x_[!is.na(x_)]
        ordered.index <- order(x)
        cumv <- abs(cumsum(weights[ordered.index]))[diff(x[ordered.index]) != 0]
        return(if (is_null(cumv)) 0 else max(cumv))
    })
    return(ks)
}
.col.var.ratio <- function(mat, treat, weights, x.types, no.weights = FALSE) {
    if (no.weights) weights <- rep(1, NROW(mat))
    ratios <- rep(NA_real_, NCOL(mat))
    non.binary <- x.types != "Binary"
    ratios[non.binary] <- col.w.v(mat[treat == 1, non.binary, drop = FALSE], w = weights[treat == 1]) / col.w.v(mat[treat == 0, non.binary, drop = FALSE], w = weights[treat == 0])
    return(pmax(ratios, 1/ratios))
}
var.ratios <- function(s0, s1, x.types) {
    ratios <- rep(NA_real_, length(s0))
    non.binary <- x.types != "Binary"
    ratios[non.binary] <- (s1[non.binary]^2) / (s0[non.binary]^2)
    return(ratios)
    # return(pmax(ratios, 1/ratios))
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
samplesize <- function(treat, weights = NULL, subclass = NULL, s.weights = NULL, method=c("matching", "weighting", "subclassification"), cluster = NULL, which.cluster = NULL, discarded = NULL, treat.names = c("Control", "Treated")) {
    #Computes sample size info. for unadjusted and adjusted samples.
    # method is what method the weights are to be used for. 
    # method="subclassification" is for subclass sample sizes only.
    
    if (is_not_null(cluster) && is_not_null(which.cluster)) in.cluster <- cluster == which.cluster
    else in.cluster <- rep(TRUE, length(treat))
    if (is_null(s.weights)) s.weights <- rep(1, length(treat))
    if (is_null(discarded)) discarded <- rep(0, length(treat))
    
    if (length(method) == 1 && method == "subclassification") {
        if (is_null(subclass)) stop("subclass must be a vector of subclasses.")
        qbins <- nlevels(subclass)
        
        nn <- as.data.frame(matrix(0, 3, qbins))
        
        dimnames(nn) <- list(c(treat.names[1], treat.names[2], "Total"), 
                             paste("Subclass", levels(subclass)))
        
        matched <- !is.na(subclass)
        k <- 0
        for (i in levels(subclass)) {
            qi <- subclass[matched]==i
            qt <- treat[matched][qi]
            if (sum(qt==1)<2|(sum(qt==0)<2)){
                if (sum(qt==1)<2)
                    warning("Not enough treatment units in subclass ", i, call. = FALSE)
                else if (sum(qt==0)<2)
                    warning("Not enough control units in subclass ", i, call. = FALSE)
            }
            k <- k + 1
            nn[, k] <- c(sum(qt==0), sum(qt==1), length(qt))
        }
        attr(nn, "tag") <- "Sample sizes by subclass"
    }
    else if (is_null(weights)) {
        
        t <- treat[in.cluster]
        sw <- s.weights[in.cluster]
        
        nn <- as.data.frame(matrix(0, ncol = 2, nrow = 1))
        nn[1, ] <- c(ESS(sw[t==0]), ESS(sw[t==1]))
        dimnames(nn) <- list(c("All"), 
                             c(treat.names[1], treat.names[2]))
        if (nunique.gt(s.weights, 2) || !any(s.weights==1) || !all(s.weights %in% c(0,1))) {
            attr(nn, "ss.type") <- c("ess")
            #attr(nn, "tag") <- "Effective sample sizes"
        }
        else {
            # nn <- as.data.frame(matrix(0, ncol=2, nrow=1))
            # nn[1, ] <- c(sum(in.cluster & treat==0), 
            #              sum(in.cluster & treat==1))
            # dimnames(nn) <- list(c("All"), 
            #                      c("Control", "Treated"))
            attr(nn, "ss.type") <- c("ss")
            #attr(nn, "tag") <- "Sample sizes"
        }
        
    }
    else if (NCOL(weights) == 1) {
        if (method=="matching") {
            nn <- as.data.frame(matrix(0, ncol=2, nrow=5))
            # nn[1, ] <- c(sum(in.cluster & treat==0), sum(in.cluster & treat==1))
            # nn[2, ] <- c(sum(in.cluster & treat==0 & weights>0), sum(in.cluster & treat==1 & weights>0))
            # nn[3, ] <- c(sum(in.cluster & treat==0 & weights==0 & discarded==0), sum(in.cluster & treat==1 & weights==0 & discarded==0))
            # nn[4, ] <- c(sum(in.cluster & treat==0 & weights==0 & discarded==1), sum(in.cluster & treat==1 & weights==0 & discarded==1))
            nn[1, ] <- c(sum(in.cluster & treat==0), 
                         sum(in.cluster & treat==1))
            nn[2, ] <- c(sum(weights[in.cluster & treat==0, 1]), 
                         sum(weights[in.cluster & treat==1, 1]))
            nn[3, ] <- c(sum(in.cluster & treat==0 & weights[,1]>0), 
                         sum(in.cluster & treat==1 & weights[,1]>0))
            nn[4, ] <- c(sum(in.cluster & treat==0 & weights[,1]==0 & discarded==0), 
                         sum(in.cluster & treat==1 & weights[,1]==0 & discarded==0))
            nn[5, ] <- c(sum(in.cluster & treat==0 & weights[,1]==0 & discarded==1), 
                         sum(in.cluster & treat==1 & weights[,1]==0 & discarded==1))
            dimnames(nn) <- list(c("All", "Matched", "Matched (Unweighted)", "Unmatched", "Discarded"), 
                                 c(treat.names[1], treat.names[2]))
            
            attr(nn, "ss.type") <- rep("ss", NROW(nn))
            #attr(nn, "tag") <- "Sample sizes"
        }
        else if (method == "weighting") {
            
            t <- treat[in.cluster]
            w <- weights[in.cluster, 1]
            sw <- s.weights[in.cluster]
            dc <- discarded[in.cluster]
            
            nn <- as.data.frame(matrix(0, ncol = 2, nrow = 3))
            nn[1, ] <- c(ESS(sw[t==0]), ESS(sw[t==1]))
            nn[2, ] <- c(ESS(w[t==0]*sw[t==0]), ESS(w[t==1]*sw[t==1]))
            nn[3, ] <- c(sum(t==0 & dc==1), 
                         sum(t==1 & dc==1))
            dimnames(nn) <- list(c("Unadjusted", "Adjusted", "Discarded"), 
                                 c(treat.names[1], treat.names[2]))
            attr(nn, "ss.type") <- c("ss", "ess")
            
            #attr(nn, "tag") <- "Effective sample sizes"
            
        }
    }
    else {
        t <- treat[in.cluster]
        sw <- s.weights[in.cluster]
        
        nn <- as.data.frame(matrix(0, ncol=2, nrow=1+NCOL(weights)))
        nn[1, ] <- c(ESS(sw[t==0]), ESS(sw[t==1]))
        for (i in seq_len(NCOL(weights))) {
            if (method[i] == "matching") {
                nn[1+i,] <- c(sum(weights[in.cluster & treat==0, i]), 
                              sum(weights[in.cluster & treat==1, i]))
            }
            else if (method[i] == "weighting") {
                w <- weights[in.cluster, i]
                nn[1+i,] <- c(ESS(w[t==0]*sw[t==0]), ESS(w[t==1]*sw[t==1]))
            }
            
        }
        dimnames(nn) <- list(c("All", names(weights)), 
                             c(treat.names[1], treat.names[2]))
        attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
        
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
max.imbal <- function(balance.table, col.name, thresh.col.name) {
    balance.table.clean <- balance.table[balance.table$Type != "Distance" & is.finite(balance.table[, col.name]),]
    maxed <- balance.table.clean[which.max(abs(balance.table.clean[, col.name])), match(c(col.name, thresh.col.name), names(balance.table.clean))]
    maxed <- data.frame(Variable = rownames(maxed), maxed)
    return(maxed)
    # return(balance.table[which.max(abs(balance.table[balance.table$Type != "Distance", col.name])), match(c(col.name, thresh.col.name), names(balance.table))])
}
balance.table <- function(C, weights, treat, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.sds = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, 
                          s.weights = rep(1, length(treat)), abs = FALSE, no.adj = FALSE, types = NULL, pooled.sds = NULL, disp.pop = FALSE, pop.means = NULL, pop.sds = NULL, quick = TRUE) {
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
    
    #Means for each group
    # if (!((!un || !disp.means) && quick)) {
        for (t in c(0, 1)) {
            B[[paste.("M", t, "Un")]] <- col.w.m(C[treat == t, , drop = FALSE], w = s.weights[treat==t])
        }
        # if (!(!disp.pop && quick)) {
        #     B[["M.Pop.Un"]] <- col.w.m(C, w = s.weights)
        # }
    # }
    if (!no.adj) {
    # if (!no.adj && !(!disp.means && quick)) {
        for (i in weight.names) {
            for (t in c(0, 1)) {
                B[[paste.("M", t, i)]] <- col.w.m(C[treat == t, , drop = FALSE], w = weights[[i]][treat==t]*s.weights[treat==t])
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
    
    sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else B[["Type"]] != "Binary"
    # if (!((!un || !disp.sds) && quick)) {
    for (t in c(0, 1)) {
        sds <- rep(NA_real_, NCOL(C))
        sds[sd.computable] <- sqrt(col.w.v(C[treat == t, sd.computable, drop = FALSE], w = s.weights[treat==t]))
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
            for (t in c(0, 1)) {
                sds <- rep(NA_real_, NCOL(C))
                sds[sd.computable] <- sqrt(col.w.v(C[treat == t, sd.computable, drop = FALSE], w = weights[[i]][treat==t]*s.weights[treat==t]))
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
    if (abs) a0 <- base::abs
    else a0 <- base::identity
    
    # if (!(!un && quick)) #Always compute unadjusted diffs
    # B[["Diff.Un"]] <- a0(col.std.diff(C, treat = treat, weights = NULL, x.types = B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom[1], no.weights = TRUE, s.weights = s.weights, pooled.sds = pooled.sds))
    B[["Diff.Un"]] <- a0(std.diffs(m0 = B[["M.0.Un"]], 
                                   s0 = B[["SD.0.Un"]], 
                                   m1 = B[["M.1.Un"]], 
                                   s1 = B[["SD.1.Un"]], 
                                   x.types = B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom[1], pooled.sds = pooled.sds))
    
    if (!no.adj) {
        # for (j in seq_len(NCOL(weights))) {
        #     B[[paste.("Diff", weight.names[j])]] <- a0(col.std.diff(C, treat = treat, weights = weights[[j]], x.types = B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom[j], no.weights = FALSE, s.weights = s.weights, pooled.sds = pooled.sds))
        # }
        for (i in weight.names) {
            B[[paste.("Diff", i)]] <- a0(std.diffs(m0 = B[[paste.("M.0", i)]], 
                                                   s0 = B[["SD.0.Un"]], 
                                                   m1 = B[[paste.("M.1", i)]], 
                                                   s1 = B[["SD.1.Un"]], 
                                                   x.types = B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom[i], pooled.sds = pooled.sds))
        }
    }
    
    #Variance ratios
    vabs <- function(x) pmax(x, 1/x)
    if (abs) v0 <- vabs
    else v0 <- base::identity
    
    if (!(!disp.v.ratio && quick)) {
        if (!(!un && quick)) {
            # B[["V.Ratio.Un"]] <- v0(col.var.ratio(C, treat, s.weights, B[["Type"]], no.weights = FALSE))
            B[["V.Ratio.Un"]] <- v0(var.ratios(s0 = B[["SD.0.Un"]], 
                                               s1 = B[["SD.1.Un"]], 
                                               x.types = B[["Type"]]))
        }
        if (!no.adj) {
            # for (j in seq_len(NCOL(weights))) {
            #     B[[paste.("V.Ratio", weight.names[j])]] <- v0(col.var.ratio(C, treat, weights[[j]]*s.weights, B[["Type"]], no.weights = FALSE))
            # }
            for (i in weight.names) {
                B[[paste.("V.Ratio", i)]] <- v0(var.ratios(s0 = B[[paste.("SD.0", i)]], 
                                                           s1 = B[[paste.("SD.1", i)]], 
                                                           x.types = B[["Type"]]))
            }
        }
    }
    if (!any(sapply(B[startsWith(names(B), "V.Ratio.")], is.finite))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
    
    #KS Statistics
    if (!(!disp.ks && quick)) {
        if (!(!un && quick)) {
            B[["KS.Un"]] <- col.ks(C, treat, s.weights, B[["Type"]], no.weights = FALSE)
        }
        if (!no.adj) {
            for (j in seq_len(NCOL(weights))) {
                B[[paste.("KS", weight.names[j])]] <- col.ks(C, treat, weights[[j]]*s.weights, B[["Type"]], no.weights = FALSE)
            }
        }
    }
    if (!any(sapply(B[startsWith(names(B), "KS.")], is.finite))) {disp.ks <- FALSE; ks.threshold <- NULL}
    
    
    if (is_not_null(m.threshold)) {
        if (no.adj) {
            B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Diff.Un"]]), paste0(ifelse(abs(B[["Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), round(m.threshold, 3)), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("M.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Diff", i)]]), paste0(ifelse(abs(B[[paste.("Diff", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), round(m.threshold, 3)), "")
            }
        }
        
    }
    if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
    
    if (is_not_null(v.threshold)) {
        if (no.adj) {
            B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["V.Ratio.Un"]]), paste0(ifelse(vabs(B[["V.Ratio.Un"]]) < v.threshold, "Balanced, <", "Not Balanced, >"), round(v.threshold, 3)), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("V.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("V.Ratio", i)]]), paste0(ifelse(vabs(B[[paste.("V.Ratio", i)]]) < v.threshold, "Balanced, <", "Not Balanced, >"), round(v.threshold, 3)), "")
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
        
        # if (!(!disp.means && quick)) {
            for (t in c(0,1)) {
                SB[[i]][[paste.("M", t, "Adj")]] <- colMeans(C[treat==t & in.subclass, , drop = FALSE])
            }
        # }
        # if (!(!disp.sds && quick)) {
            non.binary <- if (binary == "std") rep(TRUE, nrow(B)) else B[["Type"]] != "Binary"
            un.sds <- setNames(vector("list", 2), c("0", "1"))
            for (t in c(0, 1)) {
                sds <- rep(NA_real_, NCOL(C))
                sds[non.binary] <- apply(C[treat == t & in.subclass, non.binary, drop = FALSE], 2, sd)
                SB[[i]][[paste.("SD", t, "Adj")]] <- sds
                un.sds[[as.character(t)]] <- rep(NA_real_, NCOL(C))
                un.sds[[as.character(t)]][non.binary] <- apply(C[treat == t, non.binary, drop = FALSE], 2, sd)
            }
            
        # }
        
        
        #Mean differences
        # SB[[i]][["Diff.Adj"]] <- col.std.diff(C, treat=treat, weights=NULL, subclass=subclass, which.sub=i, x.types=B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom, no.weights = TRUE)
        SB[[i]][["Diff.Adj"]] <- std.diffs(m0 = SB[[i]][["M.0.Adj"]], 
                                       s0 = un.sds[["0"]], 
                                       m1 = SB[[i]][["M.1.Adj"]], 
                                       s1 = un.sds[["1"]], 
                                       x.types = B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom)
        
        #Variance ratios
        vabs <- function(x) pmax(x, 1/x)
        if (abs) v0 <- vabs
        else v0 <- base::identity
        
        if (!(!disp.v.ratio && quick)) {
            # SB[[i]][["V.Ratio.Adj"]] <- v0(col.var.ratio(C[in.subclass, ], treat = treat[in.subclass], weights = NULL, x.types = B[["Type"]], no.weights = TRUE))
            SB[[i]][["V.Ratio.Adj"]] <- v0(var.ratios(s0 = SB[[i]][["SD.0.Adj"]], 
                                               s1 = SB[[i]][["SD.1.Adj"]], 
                                               x.types = B[["Type"]]))
        }
        
        #KS Statistics
        if (!(!disp.ks && quick)) {
            SB[[i]][["KS.Adj"]] <- col.ks(C[in.subclass, ], treat = treat[in.subclass], weights = NULL, x.types = B[["Type"]], no.weights = TRUE)
        }
    }
    
    if (all(sapply(SB, function(x) !any(is.finite(c(x[["SD.0.Adj"]], x[["SD.1.Adj"]])))))) {
        attr(SB, "dont.disp.sds") <- TRUE
        disp.sds <- FALSE
    }
    
    if (is_not_null(m.threshold)) {
        for (i in levels(subclass)) {
            SB[[i]][["M.Threshold"]] <- ifelse(SB[[i]][["Type"]]=="Distance", "", 
                                               paste0(ifelse(is.finite(SB[[i]][["Diff.Adj"]]) & abs(SB[[i]][["Diff.Adj"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), round(m.threshold, 3)))
        }
    }
    
    if (all(sapply(SB, function(x) !any(is.finite(x[["V.Ratio.Adj"]]))))) {
        attr(SB, "dont.disp.v.ratio") <- TRUE; v.threshold <- NULL
        disp.v.ratio <- FALSE
    }
    if (is_not_null(v.threshold)) {
        for (i in levels(subclass)) {
            SB[[i]][["V.Threshold"]] <- ifelse(SB[[i]][["Type"]]!="Distance" & is.finite(SB[[i]][["V.Ratio.Adj"]]), 
                                               paste0(ifelse(v0(SB[[i]][["V.Ratio.Adj"]]) < v.threshold, "Balanced, <", "Not Balanced, >"), round(v.threshold, 3)), "")
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
                         pooled = "all", control = "control")
    }
    if (sub.by=="treat") {
        wsub <- "Treated"
    } else if (sub.by=="control") {
        wsub <- "Control"
    } else if (sub.by=="all") {
        wsub <- "Total"
    }
    
    B.A <- balance.table.subclass.list[[1]][c("M.0.Adj", "M.1.Adj", "Diff.Adj")]
    
    for(i in rownames(B.A)) {
        for(j in colnames(B.A)) {
            B.A[[i, j]] <- sum(vapply(seq_along(balance.table.subclass.list),
                                      function(s) subclass.obs[[wsub, s]]/sum(subclass.obs[wsub, ]) * (balance.table.subclass.list[[s]][[i, j]]), numeric(1)))
        }
    }
    B.A.df <- data.frame(balance.table[c("Type", "M.0.Un", "SD.0.Un", "M.1.Un", "SD.1.Un", "Diff.Un", "V.Ratio.Un", "KS.Un")], 
                         B.A, M.Threshold = NA_character_)
    if (is_not_null(m.threshold)) B.A.df[["M.Threshold"]] <- ifelse(B.A.df[["Type"]]=="Distance", "", paste0(ifelse(is.finite(B.A.df[["Diff.Adj"]]) & abs(B.A.df[["Diff.Adj"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
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
    else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(vapply(balance.table.clusters.list, function(y) y[[x, "Type"]], character(1))); return(u[!is.na(u)])}), use.names = FALSE)
    
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
                    B[[paste.(Fun, "KS", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else funs[[Fun]](sapply(balance.table.clusters.list, function(y) y[[x, paste.("KS", sample)]]), na.rm = TRUE), numeric(1))
                }            
            }
        }
    }
    
    return(B)
}

#base.bal.tab.cont
w.r <- function(x, y, w = NULL, s.weights = NULL) {
    #Computes weighted correlation but using the unweighted (s.weighted) variances
    #in the denominator.
    if (length(x) != length(y)) stop("x and y must the same length")
    
    if (is_null(w)) w <- rep(1, length(x))
    else if (length(w) != length(x)) stop("weights must be same length as x and y")
    
    if (is_null(s.weights)) s.weights <- rep(1, length(x))
    else if (length(s.weights) != length(x)) stop("s.weights must be same length as x and y")
    
    w_ <- w*s.weights
    
    r <- w.cov(x, y, w_) / (sqrt(w.v(x, s.weights) * w.v(y, s.weights)))
    #r <- w.cov(x, y, w_) / (sqrt(var(x) * var(y)))
    
    return(r)
}
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
    
    #Means
    if (!((!un || !disp.means) && quick)) {
        B[["M.Un"]] <- col.w.m(C, w = s.weights)
    }
    if (!no.adj && !(!disp.means && quick)) {
        for (i in weight.names) {
            B[[paste.("M", i)]] <- col.w.m(C, w = weights[[i]]*s.weights)
        }
    }
    
    #SDs
    non.binary <- B[["Type"]] != "Binary"
    if (!((!un || !disp.sds) && quick)) {
        sds <- rep(NA_real_, NCOL(C))
        sds[non.binary] <- sqrt(col.w.v(C[, non.binary, drop = FALSE], w = s.weights))
        B[["SD.Un"]] <- sds
    }
    if (!no.adj && !(!disp.sds && quick)) {
        for (i in weight.names) {
            sds <- rep(NA_real_, NCOL(C))
            sds[non.binary] <- sqrt(col.w.v(C[, non.binary, drop = FALSE], w = weights[[i]]*s.weights))
            B[[paste.("SD", i)]] <- sds
        }
    }
    if (!any(sapply(B[startsWith(names(B), "SD.")], is.finite))) {disp.sds <- FALSE}
    
    #Correlations
    if (abs) a0 <- base::abs
    else a0 <- base::identity
    # if (!(!un && quick)) #Always calculate unadjusted corrs
    B[["Corr.Un"]] <- a0(apply(C, 2, w.r, y = treat, s.weights = s.weights))
    if (!no.adj) {
        for (i in weight.names) {
            B[[paste.("Corr", i)]] <- a0(apply(C, 2, w.r, y = treat, w = weights[[i]], s.weights = s.weights))
        }
    }
    
    if (is_not_null(r.threshold)) {
        if (no.adj) {
            B[["R.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Corr.Un"]]), paste0(ifelse(abs(B[["Corr.Un"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), round(r.threshold, 3)), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("R.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Corr", i)]]), paste0(ifelse(abs(B[[paste.("Corr", i)]]) < r.threshold, "Balanced, <", "Not Balanced, >"), round(r.threshold, 3)), "")
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
        SB[[i]][["Corr.Adj"]] <- apply(C, 2, function(x) w.r(x[in.subclass], y = treat[in.subclass]))
        
    }
    
    if (all(sapply(SB, function(x) !any(is.finite(x[["SD.Adj"]]))))) {
        attr(SB, "dont.disp.sds") <- TRUE
        disp.sds <- FALSE
    }
    
    if (is_not_null(r.threshold)) {
        for (i in levels(subclass)) {
            SB[[i]][["R.Threshold"]] <- ifelse(SB[[i]][["Type"]]=="Distance", "", 
                                               paste0(ifelse(is.finite(SB[[i]][["Corr.Adj"]]) & abs(SB[[i]][["Corr.Adj"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), round(r.threshold, 3)))
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
    else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(sapply(bal.tab.imp.list, function(y) y[[x, "Type"]])); return(u[!is.na(u)])}), use.names = FALSE)
    
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
                    B[[paste.(Fun, "KS", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else funs[[Fun]](sapply(bal.tab.imp.list, function(y) y[[x, paste.("KS", sample)]]), na.rm = TRUE), numeric(1))
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
        cont.treat <- !is.na(charmatch("Mean.Corr.Un", Bcolnames))
        if (length(weight.names) <= 1) weight.names <- "Adj"
        #imp.functions <- c("Min", "Mean", "Median", "Max")
        imp.functions <- c("Min", "Mean", "Max")
        
        B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)))
        dimnames(B) <- list(Brownames, Bcolnames)
        
        if (is_not_null(types)) B[["Type"]] <- types
        else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(vapply(summary.tables, function(y) y[[x, "Type"]], character(1))); return(u[!is.na(u)])}), use.names = FALSE)
        
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
                        B[[paste.(Fun, "KS", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else funs[[Fun]](sapply(summary.tables, function(y) y[[x, paste.(Fun, "KS", sample)]]), na.rm = TRUE), numeric(1))
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
    else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(vapply(bal.tab.multi.list, function(y) y[[x, "Type"]], character(1))); return(u[!is.na(u)])}), use.names = FALSE)
    
    max_ <- function(x, na.rm = TRUE) {
        if (!any(is.finite(x))) NA_real_
        else max(x, na.rm = na.rm)
    }
    for (sample in c("Un", weight.names)) {
        if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
            B[[paste.("Max", "Diff", sample)]] <- vapply(Brownames, function(x) max_(sapply(bal.tab.multi.list, function(y) abs(y[[x, paste.("Diff", sample)]])), na.rm = TRUE), numeric(1))
            B[[paste.("Max", "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else max_(sapply(bal.tab.multi.list, function(y) y[[x, paste.("V.Ratio", sample)]]), na.rm = TRUE), numeric(1))
            B[[paste.("Max", "KS", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else max_(sapply(bal.tab.multi.list, function(y) y[[x, paste.("KS", sample)]]), na.rm = TRUE), numeric(1))
        }
    }
    
    if (is_not_null(m.threshold)) {
        if (no.adj) {
            B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Diff.Un"]]), paste0(ifelse(abs(B[["Max.Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("M.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.Diff", i)]]), paste0(ifelse(abs(B[[paste.("Max.Diff", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
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
samplesize.multi <- function(bal.tab.multi.list, treat.names, focal) {
    if (is_not_null(focal)) which <- c(treat.names[treat.names != focal], focal)
    else which <- treat.names
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
    else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, "Type"]] else NA_character_, character(1))); return(u[!is.na(u)])}), use.names = FALSE)
    
    B[["Times"]] <- Brownames.appear[Brownames]
    
    max_ <- function(x, na.rm = TRUE) {
        if (!any(is.finite(x))) NA_real_
        else max(x, na.rm = na.rm)
    }
    for (sample in c("Un", weight.names)) {
        if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
            if (cont.treat) {
                B[[paste.("Max", "Corr", sample)]] <- vapply(Brownames, function(x) max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) abs(y[[x, paste.("Corr", sample)]]) else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
            }
            else {
                B[[paste.("Max", "Diff", sample)]] <- vapply(Brownames, function(x) max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) abs(y[[x, paste.("Diff", sample)]]) else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
                B[[paste.("Max", "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, paste.("V.Ratio", sample)]] else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
                B[[paste.("Max", "KS", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, paste.("KS", sample)]] else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
            }
        }
    }
    
    if (is_not_null(m.threshold)) {
        if (no.adj) {
            B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Diff.Un"]]), paste0(ifelse(abs(B[["Max.Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste.("M.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.Diff", i)]]), paste0(ifelse(abs(B[[paste.("Max.Diff", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
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
samplesize.target <- function(bal.tab.target.list, treat.names, target.name) {
    which <- treat.names[treat.names != target.name]
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
assign.shapes <- function(colors, default.shape = "circle filled") {
    if (nunique(colors) < length(colors)) {
        shapes <- seq_int_cycle(21, 21 + length(colors), max = 25)
    }
    else shapes <- rep(default.shape, length(colors))
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

#To pass CRAN checks:
utils::globalVariables(c("distance", "addl", "addl.list", "distance.list",
                         "quick", "treat", "Sample", "min.stat",
                         "max.stat", "mean.stat", "count"))