#x2base
match.strata2weights <- function(match.strata, treat, covs = NULL) {
    #Process match.strata into weights (similar to weight.subclass from MatchIt)
    if (length(covs) == 0) names(treat) <- seq_along(treat)
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
    if (sum(weights.matched[clabels]) < sqrt(.Machine$double.eps)) #all control weights are 0
        weights.matched[clabels] <- rep(0, length(weights.matched[clabels]))
    else {
        ## Number of C units that were matched to at least 1 T
        num.cs <- sum(weights.matched[clabels] > sqrt(.Machine$double.eps))
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
    
    if (sum(weights) < sqrt(.Machine$double.eps)) 
        stop("No units were matched", call. = FALSE)
    else if (sum(weights[tlabels]) < sqrt(.Machine$double.eps))
        stop("No treated units were matched", call. = FALSE)
    else if (sum(weights[clabels]) < sqrt(.Machine$double.eps))
        stop("No control units were matched", call. = FALSE)
    return(weights)
}
use.tc.fd <- function(formula, data, treat, covs) {
    useWhich <- function(f, d, t, c) {
        #output: tc, fd, both
        good <- c(formula=0, data=0, covs=0, treat=0)
        if (length(f) > 0 & class(f)=="formula") good[1] <- 1
        if (length(d) > 0 & is.data.frame(d)) good[2] <- 1
        if (length(c) > 0 & is.data.frame(c)) good[3] <- 1
        if (length(t) > 0 & length(unique(t))==2) good[4] <- 1
        
        if (sum(good) <= 1) {
            stop("Either formula and data or treat and covs must be specified correctly.", call. = FALSE)
        }
        else if (sum(good) == 2) {
            if (sum(good[c("formula", "data")])==0) {
                use.which <- "tc"}
            else if (sum(good[c("formula", "data")])==1) {
                stop("Either formula and data or treat and covs must be specified correctly.", call. = FALSE)}
            else { #if (sum(good[c("formula", "data")])==2)
                use.which <- "fd"}
        }
        else if (sum(good) == 3) {
            bad <- names(good)[match(0, good)]
            if (any(c("formula", "data") == bad)) {
                warning(paste("Argument to", bad, "is missing; using treat and covs instead."), call. = FALSE, immediate.=TRUE)
                use.which <- "tc"}
            else {
                warning(paste("Argument to", bad, "is missing; using formula and data instead."), call. = FALSE, immediate.=TRUE)
                use.which <- "fd"}
        }
        else { #if (sum(good)==4)
            use.which <- "both"
        }
        return(use.which)
    }
    use.fd <- function(f, d){
        #outputs a list containing treat [1] and covs [2]
        out.list <- list(treat=NA, covs=data.frame(NA))
        tt <- terms(f)
        attr(tt, "intercept") <- 0
        mf<- tryCatch(model.frame(tt, d),
                      error = function(cond) stop(paste0(c("All right hand side variables of formula must be variables in data.\nVariables not in data: ",
                                                           paste(attr(tt, "term.labels")[which(!attr(tt, "term.labels") %in% names(d))], collapse=", "))), call. = FALSE))
        out.list$treat <- setNames(model.response(mf), rownames(d)) #treat
        out.list$covs <- d[attr(tt, "term.labels")] #covs
        attr(out.list, "which") <- "fd"
        return(out.list)
    }
    use.tc <- function(t, c) {
        if (length(t)!=nrow(c)) {
            stop("treat must be the same length as covs.", call. = FALSE)}
        out.list <- list(treat=setNames(t, rownames(c)), covs=c)
        attr(out.list, "which") <- "tc"
        return(out.list)
    }
    
    use.which <- useWhich(formula, data, treat, covs)
    
    if (use.which == "fd") {
        t.c <- use.fd(formula, data)}
    else if (use.which == "tc") {
        t.c <- use.tc(treat, covs)}
    else if (use.which == "both") {
        try.fd <- try({t.c <- use.fd(formula, data)})
        if (class(try.fd) == "try-error") {
            message("Formula, data, treat, and covs all supplied; ignoring formula and data.")
            t.c <- use.tc(treat, covs)}
        else {
            message("Formula, data, treat, and covs all supplied; ignoring treat and covs.")}
    }
    return(t.c)
}
process.val <- function(val, i, treat, covs, ...) {
    if (is.numeric(val)) {
        val.df <- setNames(data.frame(val), i)
    }
    else if (is.character(val)) {
        data.sets <- list(...)
        data.sets <- data.sets[sapply(data.sets, function(x) length(x) > 0)]
        if ((length(data.sets) > 0 && length(val) > max(sapply(data.sets, ncol))) || length(val) == nrow(covs) || length(val) == length(treat)){
            val.df <- setNames(data.frame(val), i)
        }
        else {
            if (length(data.sets) > 0) {
                val <- unique(val)
                val.df <- setNames(as.data.frame(matrix(NA, ncol = length(val), nrow = max(sapply(data.sets, nrow)))),
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
    if (length(val) > 0) {
        if (is.vector(val, mode = "list")) {
            val.list <- lapply(val, function(x) process.val(x, i, treat, covs, ...))
            val.list <- lapply(seq_along(val.list), function(x) {
                if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                val.list[[x]]})
            if (length(unique(sapply(val.list, nrow))) > 1) {
                stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
            }
            
            val.df <- setNames(do.call("cbind", val.list),
                               c(sapply(val.list, names)))
        }
        else {
            val.df <- process.val(val, i, treat, covs, ...)
        }
        if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
            stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
        }
    }
    return(val.df)
}
list.process <- function(i, List, ntimes, call.phrase, treat.list, covs.list, ...) {
    val.List <- List
    if (length(val.List) > 0) {
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
            if (length(val) > 0) {
                if (is.vector(val, mode = "list")) {
                    val.list <- lapply(val, function(x) process.val(x, strsplit(i, ".list")[[1]], treat.list[[ti]], covs.list[[ti]], ...))
                    val.list <- lapply(seq_along(val.list), function(x) {
                        if (ncol(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
                        val.list[[x]]})
                    if (length(unique(sapply(val.list, nrow))) > 1) {
                        stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
                    }
                    
                    val.df <- setNames(do.call("cbind", val.list),
                                       c(sapply(val.list, names)))
                }
                else {
                    val.df <- process.val(val, strsplit(i, ".list")[[1]], treat.list[[ti]], covs.list[[ti]], ...)
                }
                if (length(val.df) > 0) { if (sum(is.na(val.df)) > 0) {
                    stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
                }
                val.List[[ti]] <- val.df
            }
            
        }
        val.df.lengths <- sapply(val.List[lengths(val.List) > 0], nrow)
        if (max(val.df.lengths) != min(val.df.lengths)) {
            stop(paste("All columns in", i, "need to have the same number of rows."), call. = FALSE)
        }
    }
    return(val.List)
}

#get.C
#Functions to turn input covariates into usable form
#int.poly.f creates interactions and polynomials
#splitfactor splits factor variable into indicators (now in utilities)
#binarize transforms 2-value variable into binary (0,1)
#get.C controls flow and handles redunancy
#get.types gets variables types (contin./binary)

int.poly.f <- function(mat, ex=NULL, int=FALSE, poly=1, nunder=1, ncarrot=1) {
    #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
    #Only to be used in base.bal.tab; for general use see int.poly()
    #mat=matrix input
    #ex=matrix of variables to exclude in interactions and polynomials; a subset of df
    #int=whether to include interactions or not; currently only 2-way are supported
    #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
    #nunder=number of underscores between variables
    
    if (length(ex) > 0) d <- mat[, !colnames(mat) %in% colnames(ex)]
    else d <- mat
    nd <- ncol(d)
    nrd <- nrow(d)
    no.poly <- apply(d, 2, function(x) !nunique.gt(x, 2))
    npol <- nd - sum(no.poly)
    new <- matrix(ncol = (poly-1)*npol + .5*(nd)*(nd-1), nrow = nrd)
    nc <- ncol(new)
    new.names <- character(nc)
    if (poly > 1 && npol != 0) {
        for (i in 2:poly) {
            new[, (1 + npol*(i - 2)):(npol*(i - 1))] <- apply(d[, !no.poly], 2, function(x) x^i)
            new.names[(1 + npol*(i - 2)):(npol*(i - 1))] <- paste0(colnames(d)[!no.poly], paste0(replicate(ncarrot, "_"), collapse = ""), i)
        }
    }
    if (int && nd > 1) {
        new[,(nc - .5*nd*(nd-1) + 1):nc] <- matrix(t(apply(d, 1, combn, 2, prod)), nrow = nrd)
        new.names[(nc - .5*nd*(nd-1) + 1):nc] <- combn(colnames(d), 2, paste, collapse=paste0(replicate(nunder, "_"), collapse = ""))
    }
    
    single.value <- apply(new, 2, function(x) abs(max(x) - min(x)) < sqrt(.Machine$double.eps))
    colnames(new) <- new.names
    #new <- setNames(data.frame(new), new.names)[!single.value]
    return(new[, !single.value])
}
binarize <- function(variable) {
    if (length(unique(variable)) > 2) stop(paste0("Cannot binarize ", deparse(substitute(variable)), ": more than two levels."))
    variable.numeric <- as.numeric(variable)
    if (!is.na(match(0, unique(variable.numeric)))) zero <- 0
    #else if (1 %in% unique(as.numeric(variable))) zero <- unique(as.numeric(variable))[unique(as.numeric(variable)) != 1]
    else zero <- min(unique(variable.numeric))
    newvar <- setNames(ifelse(variable.numeric==zero, 0, 1), names(variable))
    return(newvar)
}
get.C <- function(covs, int = FALSE, addl = NULL, distance = NULL, cluster = NULL) {
    #gets C data.frame, which contains all variables for which balance is to be assessed. Used in balance.table.
    C <- covs
    if (!is.null(addl)) {
        if (!is.data.frame(addl)) {
            if (is.character(addl)) stop("The argument to addl must be a data.frame containing the the values of the additional variables you want to include in the balance assessment.", call. = FALSE)
            else stop("The argument to addl must be a data.frame. Wrap data.frame() around the argument if it is a matrix or vector.", call. = FALSE)
        }
        else {
            repeat.name.indices <- sapply(names(addl), function(x) !is.na(match(x, names(C))))
            if (any(repeat.name.indices)) {
                warning(paste("The following variables in addl have the same name as covariates and will be ignored:\n",
                              paste(names(addl)[repeat.name.indices], collapse = " ")), call. = FALSE)
                addl <- addl[!repeat.name.indices]
            }
            C <- cbind(C, addl)
        }
    } 
    
    for (i in names(C)) {
        if (is.character(C[[i]])) C[[i]] <- factor(C[[i]])
        else if (!nunique.gt(C[[i]], 2)) {
            if (is.logical(C[[i]])) C[[i]] <- as.numeric(C[[i]])
            else if (is.numeric(C[[i]])) C[[i]] <- binarize(C[[i]])
        }
        
        if (nlevels(cluster) > 0 && qr(matrix(c(C[[i]], as.numeric(cluster)), ncol = 2))$rank == 1) C <- C[names(C) != i] #Remove variable if it is the same (linear combo) as cluster variable
        else if (!is.numeric(C[[i]])) {
            C <- splitfactor(C, i, replace = TRUE, sep = "_", drop.first = FALSE, 
                             drop.singleton = FALSE)
        }
    }
    C <- as.matrix(C)
    single.value <- apply(C, 2, function(x) abs(max(x) - min(x)) < sqrt(.Machine$double.eps))
    C <- C[, !single.value, drop = FALSE]
    
    if (int) {
        #Prevent duplicate var names with _'s
        nunder <- ncarrot <- 1
        repeat {
            if (all(sapply(colnames(C), function(x) !x %in% do.call(paste, c(expand.grid(colnames(C), colnames(C)), list(sep = paste0(replicate(nunder, "_"), collapse = ""))))))) break
            else nunder <- nunder + 1
        }
        #Variable names don't contain carrots
        # repeat {
        #     if (all(sapply(names(C), function(x) !x %in% paste0(names(C), paste0(replicate(nunder, "_"), collapse = ""), "2")))) break
        #     else ncarrot <- ncarrot + 1
        # }
        C <- cbind(C, int.poly.f(C, int = TRUE, poly = 2, nunder = nunder, ncarrot = ncarrot))
    }
    #Remove duplicate & redundant variables
    #If many rows, select subset to test redundancy
    if (nrow(C) > 1500) {
        repeat {
            mini.C <- C[sample(seq_len(nrow(C)), 1000),,drop=FALSE]
            single.value <- apply(mini.C, 2, function(x) abs(max(x) - min(x)) < sqrt(.Machine$double.eps))
            if (all(!single.value)) break
        }
        suppressWarnings(C.cor <- cor(mini.C))
    }
    else suppressWarnings(C.cor <- cor(C))
    s <- !lower.tri(C.cor, diag=TRUE) & (1 - abs(C.cor) < sqrt(.Machine$double.eps))
    redundant.vars <- apply(s, 2, function(x) any(x))
    C <- C[, !redundant.vars, drop = FALSE]        
    #C<- as.matrix.data.frame(C)
    
    if (length(distance) > 0) {
        #distance <- data.frame(.distance = distance)
        #while (names(distance) %in% names(C)) {names(distance) <- paste0(names(distance), "_")}
        if (any(names(distance) %in% colnames(C))) stop("distance variable(s) share the same name as a covariate. Please ensure each variable name is unique.", call. = FALSE)
        C <- cbind(distance, C, row.names = NULL)
        attr(C, "distance.names") <- names(distance)
    }
    
    return(C)
    
}
get.types <- function(C) {
    sapply(colnames(C), function(x) {
        if (any(attr(C, "distance.names") == x)) "Distance"
        else if (nunique.gt(C[,x], 2))  "Contin."
        else "Binary"
    })
}

#base.bal.tab
w.m <- function(x, w = NULL, na.rm = TRUE) {
    if (length(w) == 0) w <- rep(1, length(x))
    return(sum(x*w, na.rm=na.rm)/sum(w, na.rm=na.rm))
}
col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
    if (length(w) == 0) {
        w <- 1
        w.sum <- nrow(mat)
    }
    else {
        w.sum <- sum(w, na.rm = na.rm)
    }
    return(colSums(mat*w, na.rm = na.rm)/w.sum)
}
w.cov.scale <- function(w) {
    (sum(w, na.rm = TRUE)^2 - sum(w^2, na.rm = TRUE)) / sum(w, na.rm = TRUE)
}
w.v <- function(x, w = NULL) {
    #return(sum(w*(x-w.m(x, w))^2, na.rm=TRUE)/(sum(w, na.rm=TRUE)-1))
    return(sum(w*(x-w.m(x, w))^2, na.rm=TRUE) / w.cov.scale(w))
}
col.w.v <- function(mat, w = NULL, na.rm = TRUE) {
    if (length(w) == 0) {
        w <- rep(1, nrow(mat))
    }
    return(colSums(t((t(mat) - col.w.m(mat, w, na.rm = na.rm))^2) * w, na.rm = na.rm) / w.cov.scale(w))
}
w.cov <- function(x, y , w = NULL) {
    wmx <- w.m(x, w)
    wmy <- w.m(y, w)
    #wcov <- sum(w*(x - wmx)*(y - wmy), na.rm = TRUE)/sum(w, na.rm = TRUE)
    wcov <- sum(w*(x - wmx)*(y - wmy), na.rm = TRUE) / w.cov.scale(w)
    return(wcov)
}
col.std.diff <- function(mat, treat, weights, subclass = NULL, which.sub = NULL, x.types, continuous, binary, s.d.denom, no.weights = FALSE, s.weights = rep(1, length(treat))) {
    if (no.weights) weights <- rep(1, nrow(mat))
    w <- weights*s.weights
    
    no.sub <- length(which.sub) == 0
    if (no.sub) ss <- w > 0
    else ss <- (!is.na(subclass) & subclass == which.sub & w > 0)
    
    if (sum(treat==0 & ss) == 0) {
        warning(paste0("There are no control units in subclass ", which.sub, "."), call. = FALSE)
        return(rep(NA, ncol(mat)))
    }
    if (sum(treat==1 & ss) == 0) {
        warning(paste0("There are no treated units in subclass ", which.sub, "."), call. = FALSE)
        return(rep(NA, ncol(mat)))
    }

    diffs <- col.w.m(mat[treat == 1 & ss, , drop = FALSE], w[treat == 1 & ss]) - 
        col.w.m(mat[treat == 0 & ss, , drop = FALSE], w[treat == 0 & ss])
    diffs[abs(diffs) < sqrt(.Machine$double.eps)] <- 0
    denoms <- rep(1, ncol(mat))
    denoms.to.std <- ifelse(x.types == "Binary", binary == "std", continuous == "std")
    
    if (any(denoms.to.std)) {
        if (s.d.denom == "control") {
            denoms[denoms.to.std] <- sqrt(col.w.v(mat[treat == 0 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 0 & ss]))
        }
        else if (s.d.denom == "treated") {
            denoms[denoms.to.std] <- sqrt(col.w.v(mat[treat == 1 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 1 & ss]))
        }
        else if (s.d.denom == "pooled") {
            denoms[denoms.to.std] <-  sqrt(.5*(col.w.v(mat[treat == 0 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 0 & ss]) +
                                               col.w.v(mat[treat == 1 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 1 & ss])))
        }
    }
    
    std.diffs <- diffs/denoms
    std.diffs[!is.finite(std.diffs)] <- NA
    
    return(std.diffs)
}
col.ks <- function(mat, treat, weights, x.types, no.weights = FALSE) {
    ks <- rep(NA_integer_, ncol(mat))
    if (no.weights) weights <- rep(1, nrow(mat))
    weights[treat == 1] <- weights[treat==1]/sum(weights[treat==1])
    weights[treat == 0] <- -weights[treat==0]/sum(weights[treat==0])
    non.binary <- x.types != "Binary"
    ks[non.binary] <- apply(mat[, non.binary, drop = FALSE], 2, function(x) {
        ordered.index <- order(x)
        cumv <- abs(cumsum(weights[ordered.index]))[diff(x[ordered.index]) != 0]
        return(if (length(cumv) > 0) max(cumv) else 0)
    })
    return(ks)
}
col.var.ratio <- function(mat, treat, weights, x.types, no.weights = FALSE) {
    if (no.weights) weights <- rep(1, nrow(mat))
    ratios <- rep(NA, ncol(mat))
    non.binary <- x.types != "Binary"
    ratios[non.binary] <- col.w.v(mat[treat == 1, non.binary, drop = FALSE], w = weights[treat == 1]) / col.w.v(mat[treat == 0, non.binary, drop = FALSE], w = weights[treat == 0])
    return(apply(matrix(c(ratios, 1/ratios), byrow = T, nrow = 2), 2, max))
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
    # method="subclassification is for subclass sample sizes only.
    
    if (nlevels(cluster) > 0 && length(which.cluster) > 0) in.cluster <- cluster == which.cluster
    else in.cluster <- rep(TRUE, length(treat))
    if (length(s.weights) == 0) s.weights <- rep(1, length(treat))
    if (length(discarded) == 0) discarded <- rep(0, length(treat))
    
    if (length(method) == 1 && method == "subclassification") {
        if (length(subclass) == 0) stop("subclass must be a vector of subclasses.")
        qbins <- length(levels(subclass))
        
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
    else if (length(weights) == 0) {
        sw <- s.weights[in.cluster]
        
        nn <- as.data.frame(matrix(0, ncol = 2, nrow = 1))
        nn[1, ] <- c((sum(sw[treat==0])^2)/sum(sw[treat==0]^2),
                     (sum(sw[treat==1])^2)/sum(sw[treat==1]^2))
        dimnames(nn) <- list(c("All"), 
                             c(treat.names[1], treat.names[2]))
        if (length(unique(s.weights)) > 2 || !any(s.weights==1) || !all(s.weights %in% c(0,1))) {
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
    else if (ncol(weights) == 1) {
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
            
            attr(nn, "ss.type") <- rep("ss", nrow(nn))
            #attr(nn, "tag") <- "Sample sizes"
        }
        else if (method == "weighting") {
            
            t <- treat[in.cluster]
            w <- weights[in.cluster, 1]
            sw <- s.weights[in.cluster]
            dc <- discarded[in.cluster]
            
            nn <- as.data.frame(matrix(0, ncol = 2, nrow = 3))
            nn[1, ] <- c((sum(sw[t==0])^2)/sum(sw[t==0]^2),
                         (sum(sw[t==1])^2)/sum(sw[t==1]^2))
            nn[2, ] <- c((sum(w[t==0]*sw[t==0])^2)/sum((w[t==0]*sw[t==0])^2),
                         (sum(w[t==1]*sw[t==1])^2)/sum((w[t==1]*sw[t==1])^2))
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
        
        nn <- as.data.frame(matrix(0, ncol=2, nrow=1+ncol(weights)))
        nn[1, ] <- c((sum(sw[t==0])^2)/sum(sw[t==0]^2), 
                     (sum(sw[t==1])^2)/sum(sw[t==1]^2))
        for (i in seq_len(ncol(weights))) {
            if (method[i] == "matching") {
                nn[1+i,] <- c(sum(weights[in.cluster & treat==0, i]), 
                              sum(weights[in.cluster & treat==1, i]))
            }
            else if (method[i] == "weighting") {
                w <- weights[in.cluster, i]
                nn[1+i,] <- c((sum(w[t==0]*sw[t==0])^2)/sum((w[t==0]*sw[t==0])^2),
                              (sum(w[t==1]*sw[t==1])^2)/sum((w[t==1]*sw[t==1])^2))
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
balance.table <- function(C, weights, treat, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, s.weights = rep(1, length(treat)), no.adj = FALSE, types = NULL, quick = FALSE) {
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    if (no.adj) weight.names <- "Adj"
    else weight.names <- names(weights)
    
    #B=Balance frame
    Bnames <- c("Type", 
                # "M.0.Un", 
                # "M.1.Un", 
                # "Diff.Un", 
                # "V.Ratio.Un",  
                # "KS.Un",
                apply(expand.grid(c("M.0", "M.1", "Diff", "M.Threshold", "V.Ratio", "V.Threshold", "KS", "KS.Threshold"),
                                  c("Un", weight.names)), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow = ncol(C), ncol = length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- varnames <- colnames(C)
    
    #Set var type (binary/continuous)
    if (length(types) > 0) B[,"Type"] <- types
    else B[,"Type"] <- get.types(C)
    
    if (!((!un || !disp.means) && quick)) {
        # B[,"M.0.Un"] <- apply(C[treat==0, , drop = FALSE], 2, w.m, w = s.weights[treat==0])
        # B[,"M.1.Un"] <- apply(C[treat==1, , drop = FALSE], 2, w.m, w = s.weights[treat==1])
        B[,"M.0.Un"] <- col.w.m(C[treat == 0, , drop = FALSE], w = s.weights[treat==0])
        B[,"M.1.Un"] <- col.w.m(C[treat == 1, , drop = FALSE], w = s.weights[treat==1])
    }
    if (!no.adj && !(!disp.means && quick)) {
        for (i in weight.names) {
            # B[[paste0("M.0.", i)]] <- apply(C[treat==0, , drop = FALSE], 2, w.m, w = weights[[i]][treat==0]*s.weights[treat==0])
            # B[[paste0("M.1.", i)]] <- apply(C[treat==1, , drop = FALSE], 2, w.m, w = weights[[i]][treat==1]*s.weights[treat==1])
            B[[paste0("M.0.", i)]] <- col.w.m(C[treat == 0, , drop = FALSE], w = weights[[i]][treat==0]*s.weights[treat==0])
            B[[paste0("M.1.", i)]] <- col.w.m(C[treat == 1, , drop = FALSE], w = weights[[i]][treat==1]*s.weights[treat==1])
        }
        
    }
    
    #Mean differences
    #if (!(!un && quick)) B[["Diff.Un"]] <- sapply(seq_along(varnames), function(i) diff.selector(x=C[,varnames[i]], treat=treat, weights=NULL, x.type=B[["Type"]][i], continuous=continuous, binary=binary, s.d.denom=s.d.denom[1], no.weights = TRUE, s.weights = s.weights))
    #if (!(!un && quick)) B[["Diff.Un"]] <- col.diff.selector(C, treat = treat, weights = NULL, x.types = B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom[1], no.weights = TRUE, s.weights = s.weights)
    if (!(!un && quick)) B[["Diff.Un"]] <- col.std.diff(C, treat = treat, weights = NULL, x.types = B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom[1], no.weights = TRUE, s.weights = s.weights)
    if (!no.adj) {
        for (j in seq_len(ncol(weights))) {
            #B[[paste0("Diff.", weight.names[j])]] <- sapply(seq_along(varnames), function(i) diff.selector(x=C[,varnames[i]], treat=treat, weights=weights[[j]], x.type=B[["Type"]][i], continuous=continuous, binary=binary, s.d.denom=s.d.denom[j], s.weights = s.weights))
            #B[[paste0("Diff.", weight.names[j])]] <- col.diff.selector(C, treat = treat, weights = weights[[j]], x.types = B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom[j], s.weights = s.weights)
            B[[paste0("Diff.", weight.names[j])]] <- col.std.diff(C, treat = treat, weights = weights[[j]], x.types = B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom[j], s.weights = s.weights)
        }
    }
    
    #Variance ratios
    if (!(!disp.v.ratio && quick)) {
        # if (!(!un && quick)) B[["V.Ratio.Un"]] <- sapply(seq_along(varnames), function(i) var.ratio(C[,varnames[i]], treat, s.weights, B[["Type"]][i], no.weights = FALSE))
        if (!(!un && quick)) B[["V.Ratio.Un"]] <- col.var.ratio(C, treat, s.weights, B[["Type"]], no.weights = FALSE)
        if (!no.adj) {
            for (j in seq_len(ncol(weights))) {
                # B[[paste0("V.Ratio.", weight.names[j])]] <- sapply(seq_along(varnames), function(i) var.ratio(C[,varnames[i]], treat, weights[[j]]*s.weights, B[["Type"]][i]))
                B[[paste0("V.Ratio.", weight.names[j])]] <- col.var.ratio(C, treat, weights[[j]]*s.weights, B[["Type"]], no.weights = FALSE)
            }
        }
    }
    if (!any(sapply(B[startsWith(names(B), "V.Ratio.")], is.finite))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
    
    #KS Statistics
    if (!(!disp.ks && quick)) {
        #if (!(!un && quick)) B[["KS.Un"]] <- sapply(seq_along(varnames), function(i) ks(C[,varnames[i]], treat, s.weights, B[["Type"]][i], no.weights = TRUE))
        if (!(!un && quick)) B[["KS.Un"]] <- col.ks(C, treat, s.weights, B[["Type"]], no.weights = FALSE)
        if (!no.adj) {
            for (j in seq_len(ncol(weights))) {
                #B[[paste0("KS.", weight.names[j])]] <- sapply(seq_along(varnames), function(i) ks(C[,varnames[i]], treat, weights[[j]]*s.weights, B[["Type"]][i]))
                B[[paste0("KS.", weight.names[j])]] <- col.ks(C, treat, weights[[j]]*s.weights, B[["Type"]], no.weights = FALSE)
            }
        }
    }
    if (!any(sapply(B[startsWith(names(B), "KS.")], is.finite))) {disp.ks <- FALSE; ks.threshold <- NULL}
    
    
    if (length(m.threshold) > 0) {
        if (no.adj) {
            B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Diff.Un"]]), paste0(ifelse(abs(B[["Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("M.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("Diff.", i)]]), paste0(ifelse(abs(B[[paste0("Diff.", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
            }
        }
        
    }
    if (no.adj || ncol(weights) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
    
    if (length(v.threshold) > 0) {
        if (no.adj) {
            B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["V.Ratio.Un"]]), paste0(ifelse(B[["V.Ratio.Un"]] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("V.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("V.Ratio.", i)]]), paste0(ifelse(B[[paste0("V.Ratio.", i)]] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
            }
        }
        
    }
    if (no.adj || ncol(weights) <= 1) names(B)[names(B) == "V.Threshold.Adj"] <- "V.Threshold"
    
    if (length(ks.threshold) > 0) {
        if (no.adj) {
            B[["KS.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["KS.Un"]]), paste0(ifelse(B[["KS.Un"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("KS.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("KS.", i)]]), paste0(ifelse(B[[paste0("KS.", i)]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
            }
        }
        
    }
    if (no.adj || ncol(weights) <= 1) names(B)[names(B) == "KS.Threshold.Adj"] <- "KS.Threshold"
    
    attr(B, "disp") <- c(v = disp.v.ratio,
                         ks = disp.ks)
    
    
    return(B)
}
balance.table.subclass <- function(C, weights = NULL, treat, subclass, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, s.weights = rep(1, length(treat)), types = NULL, quick = FALSE) {
    #Creates list SB of balance tables for each subclass
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    #B=Balance frame
    Bnames <- c("Type", "M.0.Adj", "M.1.Adj", "Diff.Adj", "M.Threshold", "V.Ratio.Adj", "V.Threshold", "KS.Adj", "KS.Threshold")
    B <- as.data.frame(matrix(nrow=ncol(C), ncol=length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- varnames <- colnames(C)
    #Set var type (binary/continuous)
    if (length(types) > 0) B[["Type"]] <- types
    else B[["Type"]] <- get.types(C)
    
    SB <- vector("list", length(levels(subclass)))
    names(SB) <- levels(subclass)
    
    #-------------------------------------
    for (i in levels(subclass)) {
        
        SB[[i]] <- B
        in.subclass <- !is.na(subclass) & subclass==i
        
        if (!(!disp.means && quick)) {
            SB[[i]][["M.0.Adj"]] <- colMeans(C[treat==0 & in.subclass, , drop = FALSE])
            SB[[i]][["M.1.Adj"]] <- colMeans(C[treat==1 & in.subclass, , drop = FALSE])
        }
        
        #Mean differences
        #SB[[i]][["Diff.Adj"]] <- sapply(seq_along(rownames(SB[[i]])), function(x) diff.selector(x=C[,rownames(SB[[i]])[x]], treat=treat, weights=NULL, subclass=subclass, which.sub=i, x.type=B[["Type"]][x], continuous=continuous, binary=binary, s.d.denom=s.d.denom, no.weights = TRUE))
        SB[[i]][["Diff.Adj"]] <- col.std.diff(C, treat=treat, weights=NULL, subclass=subclass, which.sub=i, x.types=B[["Type"]], continuous=continuous, binary=binary, s.d.denom=s.d.denom, no.weights = TRUE)
        
        #Variance ratios
        if (!(!disp.v.ratio && quick)) {
            #SB[[i]][["V.Ratio.Adj"]] <- sapply(seq_along(rownames(SB[[i]])), function(x) var.ratio(C[,rownames(SB[[i]])[x]][in.subclass], treat[in.subclass], weights=NULL, var.type=B[["Type"]][x], no.weights = TRUE))
            SB[[i]][["V.Ratio.Adj"]] <- col.var.ratio(C[in.subclass, ], treat = treat[in.subclass], weights = NULL, x.types = B[["Type"]], no.weights = TRUE)
        }
        
        #KS Statistics
        if (!(!disp.ks && quick)) {
            #SB[[i]][["KS.Adj"]] <- sapply(seq_along(rownames(SB[[i]])), function(x) ks(C[,rownames(SB[[i]])[x]][in.subclass], treat[in.subclass], weights=NULL, var.type=B[["Type"]][x], no.weights = TRUE))
            SB[[i]][["KS.Adj"]] <- col.ks(C[in.subclass, ], treat = treat[in.subclass], weights = NULL, x.types = B[["Type"]], no.weights = TRUE)
        }
    }
    
    if (length(m.threshold) > 0) {
        for (i in levels(subclass)) {
            SB[[i]][["M.Threshold"]] <- ifelse(SB[[i]][["Type"]]=="Distance", "", 
                                              paste0(ifelse(is.finite(SB[[i]][["Diff.Adj"]]) & abs(SB[[i]][["Diff.Adj"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
        }
    }
    
    if (all(sapply(SB, function(x) !any(is.finite(x[["V.Ratio.Adj"]]))))) {
        attr(SB, "dont.disp.v.ratio") <- TRUE; v.threshold <- NULL
    }
    if (length(v.threshold) > 0) {
        for (i in levels(subclass)) {
            SB[[i]][["V.Threshold"]] <- ifelse(SB[[i]][["Type"]]!="Distance" & is.finite(SB[[i]][["V.Ratio.Adj"]]), 
                                              paste0(ifelse(SB[[i]][["V.Ratio.Adj"]] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
        }
    }
    if (all(sapply(SB, function(x) !any(is.finite(x[["KS.Adj"]]))))) {
        attr(SB, "dont.disp.ks") <- TRUE
    }
    if (length(ks.threshold) > 0) {
        for (i in levels(subclass)) {
            SB[[i]][["KS.Threshold"]] <- ifelse(SB[[i]][["Type"]]!="Distance" & is.finite(SB[[i]][["KS.Adj"]]), 
                                               paste0(ifelse(SB[[i]][["KS.Adj"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
        }
    }
    
    attr(SB, "thresholds") <- c(m = m.threshold,
                                v = v.threshold,
                                ks = ks.threshold)
    
    return(SB)
}
balance.table.across.subclass <- function(balance.table, balance.table.subclass.list, subclass.obs, sub.by = NULL, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, s.d.denom = NULL) {
    #Variance ratio, v.threshold, and KS not yet supported
    
    if (length(s.d.denom) > 0){
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
            B.A[[i, j]] <- sum(sapply(seq_along(balance.table.subclass.list),
                                    function(s) subclass.obs[[wsub, s]]/sum(subclass.obs[wsub, ]) * (balance.table.subclass.list[[s]][[i, j]])))
        }
    }
    B.A.df <- data.frame(balance.table[c("Type", "M.0.Un", "M.1.Un", "Diff.Un")], 
                         B.A, M.Threshold = NA)
    if (length(m.threshold) > 0) B.A.df[["M.Threshold"]] <- ifelse(B.A.df[["Type"]]=="Distance", "", paste0(ifelse(is.finite(B.A.df[["Diff.Adj"]]) & abs(B.A.df[["Diff.Adj"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
    return(B.A.df)
}
balance.table.cluster.summary <- function(balance.table.clusters.list, weight.names = NULL, no.adj = FALSE, quick = quick, types = NULL) {
    
    cont.treat <- !is.na(match("Corr.Un", unique(do.call("c", lapply(balance.table.clusters.list, names)))))
    if (no.adj) weight.names <- "Adj"
    
    Brownames <- unique(do.call("c", lapply(balance.table.clusters.list, rownames)))
    cluster.functions <- c("Min", "Mean", "Median", "Max")
    stats <- if (cont.treat) "Corr" else c("Diff", "V.Ratio", "KS")
    Bcolnames <- c("Type", apply(expand.grid(cluster.functions, stats, c("Un", weight.names)), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
    names(B) <- Bcolnames
    
    if (length(types) > 0) B[["Type"]] <- types
    else B[["Type"]] <- unlist(sapply(Brownames, function(x) {u <- unique(sapply(balance.table.clusters.list, function(y) y[[x, "Type"]])); return(u[!is.na(u)])}), use.names = FALSE)
    
    funs <- structure(vector("list", length(cluster.functions)), names = cluster.functions)
    for (Fun in cluster.functions[c(!quick, TRUE, TRUE, TRUE)]) {
        funs[[Fun]] <- function(x, ...) {
            if (!any(is.finite(x))) NA
            else get(tolower(Fun))(x, ...)
        }
        for (sample in c("Un", weight.names)) {
            if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
                if (cont.treat) {
                    B[[paste(Fun, "Corr", sample, sep = ".")]] <- sapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) abs(y[[x, paste0("Corr.", sample)]])), na.rm = TRUE))
                }
                else {
                    B[[paste(Fun, "Diff", sample, sep = ".")]] <- sapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) abs(y[[x, paste0("Diff.", sample)]])), na.rm = TRUE))
                    B[[paste(Fun, "V.Ratio", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else funs[[Fun]](sapply(balance.table.clusters.list, function(y) y[[x, paste0("V.Ratio.", sample)]]), na.rm = TRUE))
                    B[[paste(Fun, "KS", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else funs[[Fun]](sapply(balance.table.clusters.list, function(y) y[[x, paste0("KS.", sample)]]), na.rm = TRUE))
                }            
            }
        }
    }
    
    return(B)
}

#base.bal.tab.cont
w.r <- function(x, y, w = NULL) {
    if (length(x) != length(y)) stop("x and y must the same length")
    if (length(w) == 0) w <- rep(1, length(x))
    else if (length(w) != length(x)) stop("weights must be same length as x and y")
    
    r <- w.cov(x, y, w) / (sqrt(w.cov(x, x, w) * w.cov(y, y, w)))
    return(r)
}
samplesize.cont <- function(treat, weights = NULL, subclass = NULL, s.weights = NULL, method=c("matching", "weighting", "subclassification"), cluster = NULL, which.cluster = NULL, discarded = NULL) {
    #Computes sample size info. for unadjusted and adjusted samples.
    # method is what method the weights are to be used for. 
    # method="subclassification" is for subclass sample sizes only.
    #method <- match.arg(method)
    if (nlevels(cluster) > 0 && length(which.cluster) > 0) in.cluster <- cluster == which.cluster
    else in.cluster <- rep(TRUE, length(treat))
    if (length(discarded) == 0) discarded <- rep(0, length(treat))
    
    if (length(method) == 1 && method == "subclassification") {
        #stop("Subclassification is not yet surpported with continuous treatments.", call. = FALSE)
        if (length(subclass) == 0) stop("subclass must be a vector of subclasses.")
        qbins <- length(levels(subclass))
        
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
    else if (length(weights) == 0) {
        nn <- as.data.frame(matrix(0, ncol = 1, nrow = 1))
        if (length(unique(s.weights)) > 2 || !any(s.weights==1) || !all(s.weights %in% c(0,1))) {
            sw <- s.weights[in.cluster]
            
            nn[1, ] <- (sum(sw)^2)/sum(sw^2)
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
            nn[1, ] <- (sum(sw)^2)/sum(sw^2)
            nn[2, ] <- (sum(w*sw)^2)/sum((w*sw)^2)
            dimnames(nn) <- list(c("Unadjusted", "Adjusted"), 
                                 c("Total"))
            attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
            #attr(nn, "tag") <- "Effective sample sizes"
        }
    }
    else {
        #t <- treat[in.cluster]
        sw <- s.weights[in.cluster]
        nn <- as.data.frame(matrix(0, ncol=1, nrow=1+ncol(weights)))
        nn[1, ] <- (sum(sw)^2)/sum(sw^2)
        for (i in seq_len(ncol(weights))) {
            if (method[i] == "matching") {
                nn[1+i,] <- c(sum(in.cluster & weights[,i] > 0))
            }
            else if (method[i] == "weighting") {
                w <- weights[in.cluster, i]
                nn[1+i,] <- (sum(w*sw)^2)/sum((w*sw)^2)
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
balance.table.cont <- function(C, weights, treat, r.threshold = NULL, un = FALSE, s.weights = rep(1, length(treat)), no.adj = FALSE, types = NULL, quick = FALSE) {
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    if (no.adj) weight.names <- "Adj"
    else weight.names <- names(weights)
    
    #B=Balance frame
    Bnames <- c("Type", 
                "Corr.Un", 
                apply(expand.grid(c("Corr", "R.Threshold"),
                                  weight.names), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow=ncol(C), ncol=length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- varnames <- colnames(C)
    
    #Set var type (binary/continuous)
    if (length(types) > 0) B[["Type"]] <- types
    else B[["Type"]] <- get.types(C)
    
    #Correlations
    if (!(!un && quick)) B[["Corr.Un"]] <- apply(C, 2, w.r, y = treat, w = s.weights)
    if (!no.adj) {
        for (i in weight.names) {
            B[[paste0("Corr.", i)]] <- apply(C, 2, w.r, y = treat, w = weights[[i]]*s.weights)
        }
    }
    
    if (length(r.threshold) > 0) {
        if (no.adj) {
            #Call Adj, but really Un. Needs to be this way.
            B[["R.Threshold.Adj"]] <- ifelse(B[["Type"]]=="Distance" | !is.finite(B[["Corr.Un"]]), "", paste0(ifelse(abs(B[["Corr.Un"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold))
        }
        else {
            for (i in weight.names) {
                B[[paste0("R.Threshold.", i)]] <- ifelse(B[["Type"]]=="Distance" | !is.finite(B[[paste0("Corr.", i)]]), "", paste0(ifelse(abs(B[[paste0("Corr.", i)]]) < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold))
            }
        }
        
    }
    
    if (no.adj || ncol(weights) <= 1) names(B)[grepl("R.Threshold", names(B), fixed = TRUE)] <- "R.Threshold"
    
    attr(B, "thresholds") <- c(r = r.threshold)
    return(B)
}
balance.table.subclass.cont <- function(C, weights = NULL, treat, subclass, r.threshold = NULL, s.weights = rep(1, length(treat)), types = NULL, quick = FALSE) {
    #Creates list SB of balance tables for each subclass
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    #B=Balance frame
    Bnames <- c("Type", "Corr.Adj", "R.Threshold")
    B <- as.data.frame(matrix(nrow=ncol(C), ncol=length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- varnames <- colnames(C)
    #Set var type (binary/continuous)
    if (length(types) > 0) B[["Type"]] <- types
    else B[["Type"]] <- get.types(C)
    
    SB <- vector("list", length(levels(subclass)))
    names(SB) <- levels(subclass)
    
    #-------------------------------------
    for (i in levels(subclass)) {
        
        SB[[i]] <- B
        in.subclass <- !is.na(subclass) & subclass==i
        
        #Correlations
        SB[[i]][["Corr.Adj"]] <- apply(C, 2, function(x) w.r(x[in.subclass], y = treat[in.subclass]))
        
    }
    
    if (length(r.threshold) > 0) {
        for (i in levels(subclass)) {
            SB[[i]][["R.Threshold"]] <- ifelse(SB[[i]][["Type"]]=="Distance", "", 
                                              paste0(ifelse(is.finite(SB[[i]][["Corr.Adj"]]) & abs(SB[[i]][["Corr.Adj"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold))
        }
    }
    
    attr(SB, "thresholds") <- c(r = r.threshold)
    
    return(SB)
}
balance.table.across.subclass.cont <- function(balance.table, balance.table.subclass.list, subclass.obs, sub.by = NULL, r.threshold = NULL) {
    #Not specified
}

#base.bal.tab.imp
balance.table.imp.summary <- function(bal.tab.imp.list, weight.names = NULL, no.adj = FALSE, quick = FALSE, types = NULL) {
    if (!is.na(match("bal.tab", unique(do.call("c", lapply(bal.tab.imp.list, class)))))) {
        bal.tab.imp.list <- lapply(bal.tab.imp.list, function(x) x[["Balance"]])}
    cont.treat <- !is.na(match("Corr.Un", unique(do.call("c", lapply(bal.tab.imp.list, names)))))
    if (length(weight.names) <= 1) weight.names <- "Adj"
    
    Brownames <- unique(do.call("c", lapply(bal.tab.imp.list, rownames)))
    imp.functions <- c("Min", "Mean", "Median", "Max")
    stats <- if (cont.treat) "Corr" else c("Diff", "V.Ratio", "KS")
    Bcolnames <- c("Type", apply(expand.grid(imp.functions, stats, c("Un", weight.names)), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
    names(B) <- Bcolnames
    
    if (length(types) > 0) B[["Type"]] <- types
    else B[["Type"]] <- unlist(sapply(Brownames, function(x) {u <- unique(sapply(bal.tab.imp.list, function(y) y[[x, "Type"]])); return(u[!is.na(u)])}), use.names = FALSE)
    
    funs <- structure(vector("list", length(imp.functions)), names = imp.functions)
    for (Fun in imp.functions[c(!quick, TRUE, TRUE, TRUE)]) {
        funs[[Fun]] <- function(x, ...) {
            if (!any(is.finite(x))) NA
            else get(tolower(Fun))(x, ...)
        }
        for (sample in c("Un", weight.names)) {
            if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
                if (cont.treat) {
                    B[[paste(Fun, "Corr", sample, sep = ".")]] <- sapply(Brownames, function(x) funs[[Fun]](sapply(bal.tab.imp.list, function(y) abs(y[x, paste0("Corr.", sample)])), na.rm = TRUE))
                }
                else {
                    B[[paste(Fun, "Diff", sample, sep = ".")]] <- sapply(Brownames, function(x) funs[[Fun]](sapply(bal.tab.imp.list, function(y) abs(y[[x, paste0("Diff.", sample)]])), na.rm = TRUE))
                    B[[paste(Fun, "V.Ratio", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else funs[[Fun]](sapply(bal.tab.imp.list, function(y) y[[x, paste0("V.Ratio.", sample)]]), na.rm = TRUE))
                    B[[paste(Fun, "KS", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else funs[[Fun]](sapply(bal.tab.imp.list, function(y) y[[x, paste0("KS.", sample)]]), na.rm = TRUE))
                }
            }
        }
    }
    return(B)
}
balance.table.clust.imp.summary <- function(summary.tables, weight.names = NULL, no.adj = FALSE, quick = FALSE, types = NULL) {
    #cont.treat <- !is.na(match("bal.tab.cont", unique(do.call("c", lapply(bal.tab.imp.list, class)))))
    #clusters <- unique(do.call("c", lapply(bal.tab.imp.list, function(x) names(x[["Cluster.Balance"]]))))
    #cluster.tables <- lapply(clusters, function(x) lapply(bal.tab.imp.list, function(y) y[["Cluster.Balance"]][[x]]))
    #cluster.balance.across.imps <- lapply(cluster.tables, balance.table.imp.summary, no.adj, quick, types)
    #names(cluster.balance.across.imps) <- clusters
    
    if (!all(sapply(summary.tables, function(x) length(x) == 0))) {
        Brownames <- unique(do.call("c", lapply(summary.tables, rownames)))
        Bcolnames <- unique(do.call("c", lapply(summary.tables, colnames)))
        cont.treat <- !is.na(charmatch("Mean.Corr.Un", Bcolnames))
        if (length(weight.names) <= 1) weight.names <- "Adj"
        imp.functions <- c("Min", "Mean", "Median", "Max")
        stats <- if (cont.treat) "Corr" else c("Diff", "V.Ratio", "KS")
        
        B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)))
        dimnames(B) <- list(Brownames, Bcolnames)
        
        if (length(types) > 0) B[["Type"]] <- types
        else B[["Type"]] <- unlist(sapply(Brownames, function(x) {u <- unique(sapply(summary.tables, function(y) y[[x, "Type"]])); return(u[!is.na(u)])}), use.names = FALSE)
        
        funs <- structure(vector("list", length(imp.functions)), names = imp.functions)
        for (Fun in imp.functions[c(!quick, TRUE, TRUE, TRUE)]) {
            funs[[Fun]] <- function(x, ...) {
                if (!any(is.finite(x))) NA
                else get(tolower(Fun))(x, ...)
            }
            for (sample in c("Un", weight.names)) {
                if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
                    if (cont.treat) {
                        B[[paste(Fun, "Corr", sample, sep = ".")]] <- sapply(Brownames, function(x) funs[[Fun]](sapply(summary.tables, function(y) abs(y[[x, paste(Fun, "Corr", sample, sep = ".")]])), na.rm = TRUE))
                    }
                    else {
                        B[[paste(Fun, "Diff", sample, sep = ".")]] <- sapply(Brownames, function(x) funs[[Fun]](sapply(summary.tables, function(y) abs(y[[x, paste(Fun, "Diff", sample, sep = ".")]])), na.rm = TRUE))
                        B[[paste(Fun, "V.Ratio", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else funs[[Fun]](sapply(summary.tables, function(y) y[[x, paste(Fun, "V.Ratio", sample, sep = ".")]]), na.rm = TRUE))
                        B[[paste(Fun, "KS", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else funs[[Fun]](sapply(summary.tables, function(y) y[[x, paste(Fun, "KS", sample, sep = ".")]]), na.rm = TRUE))
                    }
                }
            }
        }
    }
    else B <- NULL
    # out <- list(cluster.balance.across.imps, B)
    # names(out) <- c("Cluster.Balance.Across.Imputations",
    #                 "Cluster.Summary.Across.Imputations")
    # return(out)
    return(B)
}
samplesize.across.imps <- function(obs.list) {
    #obs.list <- lapply(bal.tab.imp.list, function(x) x[["Observations"]])
    
    obs <- Reduce("+", obs.list)/length(obs.list)
    attr(obs, "tag") <- paste0("Average ", tolower(attr(obs.list[[1]], "tag")), " across imputations")
    return(obs)
}

#base.bal.tab.multi
balance.table.multi.summary <- function(bal.tab.multi.list, weight.names = NULL, no.adj = FALSE, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, quick = FALSE, types = NULL) {
    if (!is.na(match("bal.tab", unique(do.call("c", lapply(bal.tab.multi.list, class)))))) {
        bal.tab.multi.list <- lapply(bal.tab.multi.list, function(x) x[["Balance"]])}
    if (length(weight.names) <= 1) weight.names <- "Adj"
    
    Brownames <- unique(do.call("c", lapply(bal.tab.multi.list, rownames)))
    stats <- c("Diff", "V.Ratio", "KS")
    Bcolnames <- c("Type", expand.grid_string(c("Max.Diff", "M.Threshold", "Max.V.Ratio", "V.Threshold", "Max.KS", "KS.Threshold"), 
                                             c("Un", weight.names), collapse = "."))
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
    names(B) <- Bcolnames
    
    if (length(types) > 0) B[["Type"]] <- types
    else B[["Type"]] <- unlist(sapply(Brownames, function(x) {u <- unique(sapply(bal.tab.multi.list, function(y) y[[x, "Type"]])); return(u[!is.na(u)])}), use.names = FALSE)
    
    max_ <- function(x, na.rm = TRUE) {
        if (!any(is.finite(x))) NA
        else max(x, na.rm = na.rm)
    }
    for (sample in c("Un", weight.names)) {
        if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
                B[[paste("Max", "Diff", sample, sep = ".")]] <- sapply(Brownames, function(x) max_(sapply(bal.tab.multi.list, function(y) abs(y[[x, paste0("Diff.", sample)]])), na.rm = TRUE))
                B[[paste("Max", "V.Ratio", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else max_(sapply(bal.tab.multi.list, function(y) y[[x, paste0("V.Ratio.", sample)]]), na.rm = TRUE))
                B[[paste("Max", "KS", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else max_(sapply(bal.tab.multi.list, function(y) y[[x, paste0("KS.", sample)]]), na.rm = TRUE))
        }
    }
   
    if (length(m.threshold) > 0) {
        if (no.adj) {
            B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Diff.Un"]]), paste0(ifelse(abs(B[["Max.Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("M.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("Max.Diff.", i)]]), paste0(ifelse(abs(B[[paste0("Max.Diff.", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
    
    if (length(v.threshold) > 0) {
        if (no.adj) {
            B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.V.Ratio.Un"]]), paste0(ifelse(B[, "Max.V.Ratio.Un"] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("V.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("Max.V.Ratio.", i)]]), paste0(ifelse(B[[paste0("Max.V.Ratio.", i)]] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "V.Threshold.Adj"] <- "V.Threshold"
    
    if (length(ks.threshold) > 0) {
        if (no.adj) {
            B[["KS.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.KS.Un"]]), paste0(ifelse(B[["Max.KS.Un"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("KS.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("Max.KS.", i)]]), paste0(ifelse(B[[paste0("Max.KS.", i)]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "KS.Threshold.Adj"] <- "KS.Threshold"

    return(B)
}
samplesize.multi <- function(bal.tab.multi.list, treat.names, focal) {
    if (length(focal) > 0) which <- c(treat.names[treat.names != focal], focal)
    else which <- treat.names
    obs <- do.call("cbind", unname(lapply(bal.tab.multi.list, function(x) x[["Observations"]])))[, which]
    attr(obs, "tag") <- attr(bal.tab.multi.list[[1]][["Observations"]], "tag")
    attr(obs, "ss.type") <- attr(bal.tab.multi.list[[1]][["Observations"]], "ss.type")
    return(obs)
}

#base.bal.tab.msm
balance.table.msm.summary <- function(bal.tab.msm.list, weight.names = NULL, no.adj = FALSE, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, quick = FALSE, types = NULL) {
    if (!is.na(match("bal.tab", unique(do.call("c", lapply(bal.tab.msm.list, class)))))) {
        bal.tab.msm.list <- lapply(bal.tab.msm.list, function(x) x[["Balance"]])}
    cont.treat <- !is.na(match("Corr.Un", unique(do.call("c", lapply(bal.tab.msm.list, names)))))
    if (length(weight.names) <= 1) weight.names <- "Adj"
    
    Brownames <- unique(do.call("c", lapply(bal.tab.msm.list, rownames)))
    Brownames.appear <- sapply(Brownames, function(x) paste(seq_along(bal.tab.msm.list)[sapply(bal.tab.msm.list, function(y) x %in% rownames(y))], collapse = ", "))
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

    if (length(types) > 0) B[["Type"]] <- types
    else B[["Type"]] <- unlist(sapply(Brownames, function(x) {u <- unique(sapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, "Type"]] else NA)); return(u[!is.na(u)])}), use.names = FALSE)
    
    B[["Times"]] <- Brownames.appear[Brownames]
    
    max_ <- function(x, na.rm = TRUE) {
        if (!any(is.finite(x))) NA
        else max(x, na.rm = na.rm)
    }
    for (sample in c("Un", weight.names)) {
        if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
            if (cont.treat) {
                B[[paste("Max", "Corr", sample, sep = ".")]] <- sapply(Brownames, function(x) max_(sapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) abs(y[[x, paste0("Corr.", sample)]]) else NA), na.rm = TRUE))
            }
            else {
                B[[paste("Max", "Diff", sample, sep = ".")]] <- sapply(Brownames, function(x) max_(sapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) abs(y[[x, paste0("Diff.", sample)]]) else NA), na.rm = TRUE))
                B[[paste("Max", "V.Ratio", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else max_(sapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, paste0("V.Ratio.", sample)]] else NA), na.rm = TRUE))
                B[[paste("Max", "KS", sample, sep = ".")]] <- sapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA else max_(sapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, paste0("KS.", sample)]] else NA), na.rm = TRUE))
            }
        }
    }
    
    if (length(m.threshold) > 0) {
        if (no.adj) {
            B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Diff.Un"]]), paste0(ifelse(abs(B[["Max.Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("M.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("Max.Diff.", i)]]), paste0(ifelse(abs(B[[paste0("Max.Diff.", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
    
    if (length(v.threshold) > 0) {
        if (no.adj) {
            B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.V.Ratio.Un"]]), paste0(ifelse(B[, "Max.V.Ratio.Un"] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("V.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("Max.V.Ratio.", i)]]), paste0(ifelse(B[[paste0("Max.V.Ratio.", i)]] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "V.Threshold.Adj"] <- "V.Threshold"
    
    if (length(ks.threshold) > 0) {
        if (no.adj) {
            B[["KS.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.KS.Un"]]), paste0(ifelse(B[["Max.KS.Un"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("KS.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("Max.KS.", i)]]), paste0(ifelse(B[[paste0("Max.KS.", i)]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "KS.Threshold.Adj"] <- "KS.Threshold"
    
    if (length(r.threshold) > 0) {
        if (no.adj) {
            B[["R.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Corr.Un"]]), paste0(ifelse(B[["Max.Corr.Un"]] < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold), "")
        }
        else {
            for (i in weight.names) {
                B[[paste0("R.Threshold.", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste0("Max.Corr.", i)]]), paste0(ifelse(B[[paste0("Max.KCorr.", i)]] < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold), "")
            }
        }
    }
    if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "R.Threshold.Adj"] <- "R.Threshold"
    
    
    return(B)
}
samplesize.msm<- function(bal.tab.msm.list) {
    obs <- do.call("cbind", lapply(bal.tab.msm.list, function(x) x[["Observations"]]))
    attr(obs, "tag") <- attr(bal.tab.msm.list[[1]][["Observations"]], "tag")
    attr(obs, "ss.type") <- attr(bal.tab.msm.list[[1]][["Observations"]], "ss.type")
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
assign.shapes <- function(colors, default.shape = 21) {
    if (length(unique(colors)) < length(colors)) {
        shapes <- seq_int_cycle(21, 21 + length(colors), max = 25)
    }
    else shapes <- rep(default.shape, length(colors))
    return(shapes)
}
shapes.ok <- function(shapes, nshapes) {
    return((length(shapes) == 1 || length(shapes) == nshapes) && is.numeric(shapes) && all(shapes %in% 1:25))
}
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#bal.plot
get.var.from.list.with.time <- function(var.name, covs.list) {
    var.name.in.covs <- sapply(covs.list, function(x) !is.na(match(var.name, names(x))))
    n.times.appeared <- sum(var.name.in.covs)
    var <- unlist(lapply(covs.list[var.name.in.covs], function(x) x[[var.name]]))
    times <- rep(var.name.in.covs, each = ncol(covs.list[[1]]))
    return(list(var = var, times = times))
}

#print.bal.tab
round_df <- function(df, digits) {
    nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
    df[, nums] <- round(df[, nums], digits = digits)
    return(df)
}
replaceNA <- function(x) {
    x[is.na(x)] <- ""
    return(x)
}
round_df_char <- function(df, digits, pad = "0") {
    nas <- is.na(df)
    if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
    rn <- rownames(df)
    cn <- colnames(df)
    df <- as.data.frame(lapply(df, function(col) {
        if (suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
            as.numeric(as.character(col))
        } else {
            col
        }
    }), stringsAsFactors = FALSE)
    nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
    df[nums] <- round(df[nums], digits = digits)
    df[nas] <- ""
    
    df <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
    
    for (i in which(nums)) {
        if (any(grepl(".", df[[i]], fixed = TRUE))) {
            s <- strsplit(df[[i]], ".", fixed = TRUE)
            lengths <- lengths(s)
            digits.r.of.. <- sapply(seq_along(s), function(x) {
                if (lengths[x] > 1) nchar(s[[x]][lengths[x]])
                else 0 })
            df[[i]] <- sapply(seq_along(df[[i]]), function(x) {
                if (df[[i]][x] == "") ""
                else if (lengths[x] <= 1) {
                    paste0(c(df[[i]][x], rep(".", pad == 0), rep(pad, max(digits.r.of..) - digits.r.of..[x] + as.numeric(pad != 0))),
                           collapse = "")
                }
                else paste0(c(df[[i]][x], rep(pad, max(digits.r.of..) - digits.r.of..[x])),
                            collapse = "")
            })
        }
    }
    
    if (length(rn) > 0) rownames(df) <- rn
    if (length(cn) > 0) names(df) <- cn
    
    return(df)
}

#To pass CRAN checks:
utils::globalVariables(c("distance", "addl", "addl.list", "distance.list",
                         "quick", "treat", "Sample", "min.stat",
                         "max.stat", "mean.stat"))