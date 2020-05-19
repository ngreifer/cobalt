#Balance functions for use in programming

col_w_mean <- function(mat, weights = NULL, s.weights = NULL, subset = NULL, na.rm = TRUE, ...) {
    
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    if (needs.splitting) {
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                      A))
    }
    
    possibly.supplied <- c("mat", "weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
    
    weights <- weights * s.weights
    
    return(col.w.m(mat[subset, , drop = FALSE], w = weights[subset], na.rm = na.rm))
}
col_w_sd <- function(mat, weights = NULL, s.weights = NULL, bin.vars, subset = NULL, na.rm = TRUE, ...) {
    
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(to.split <- vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(bin.vars, mat)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2",
                                           split.with = bin.vars),
                                      A))
        bin.vars <- attr(mat, "split.with")[[1]]
    }
    
    possibly.supplied <- c("mat", "weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
    
    weights <- weights * s.weights
    
    return(sqrt(col.w.v(mat[subset, , drop = FALSE], w = weights[subset], 
                        bin.vars = bin.vars, na.rm = na.rm)))
}

col_w_smd <- function(mat, treat, weights = NULL, std = TRUE, s.d.denom = "pooled", abs = FALSE, s.weights = NULL, bin.vars, subset = NULL, weighted.weights = weights, na.rm = TRUE, ...) {
    # allowable.s.d.denoms <- c("pooled", "all", "weighted")
    # if (length(treat_names(treat)) == 2) allowable.s.d.denoms <- c(allowable.s.d.denoms, "treated", "control")
    
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(to.split <- vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(bin.vars, mat)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2",
                                           split.with = bin.vars),
                                      A))
        bin.vars <- attr(mat, "split.with")[[1]]
    }
    
    if (missing(treat) || !(is.factor(treat) || is.atomic(treat))) stop("treat must be an atomic vector or factor.")
    if (!is.atomic(std) || anyNA(as.logical(std)) ||
        length(std) %nin% c(1L, NCOL(mat))) {
        stop("std must be a logical vector with length equal to 1 or the number of columns of mat.")
    }
    
    if (!(is.atomic(abs) && length(abs) == 1L && !anyNA(as.logical(abs)))) {
        stop("abs must be a logical of length 1.")
    }
    
    possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
    
    if (!is_binary(treat[subset])) stop("treat must be a binary variable.")
    
    weights <- weights * s.weights
    
    if (length(std) == 1L) std <- rep(std, NCOL(mat))
    
    tval1_0 <- treat[1]
    
    m1 <- col.w.m(mat[treat==tval1_0 & subset, , drop = FALSE], weights[treat==tval1_0 & subset], na.rm = na.rm)
    m0 <- col.w.m(mat[treat!=tval1_0 & subset, , drop = FALSE], weights[treat!=tval1_0 & subset], na.rm = na.rm)
    diffs <- m1 - m0
    zeros <- check_if_zero(diffs)
    
    if (any(to.sd <- std & !zeros)) {
        denoms <- compute_s.d.denom(mat, treat = treat, 
                                    s.d.denom = s.d.denom, s.weights = s.weights, 
                                    bin.vars = bin.vars, subset = subset, to.sd = to.sd,
                                    weighted.weights = weighted.weights, na.rm = na.rm)
        
        diffs[to.sd] <- diffs[to.sd]/denoms[to.sd]
    }

    if (abs) diffs <- abs(diffs)
    else {
        tval1 <- treat[subset][binarize(treat[subset])==1][1]
        if (tval1 != tval1_0) diffs <- -1*diffs
    }

    return(setNames(diffs, colnames(mat)))
    
}
col_w_vr <- function(mat, treat, weights = NULL, abs = FALSE, s.weights = NULL, bin.vars, subset = NULL, na.rm = TRUE, ...) {
    
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(to.split <- vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(bin.vars, mat)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2",
                                           split.with = bin.vars),
                                      A))
        bin.vars <- attr(mat, "split.with")[[1]]
    }
    
    if (missing(treat) || !(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
    
    if (!(is.atomic(abs) && length(abs) == 1L && !anyNA(as.logical(abs)))) {
        stop("abs must be a logical of length 1.")
    }
    
    possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
    
    weights <- weights * s.weights
    
    weights <- weights[subset]
    treat <- treat[subset]
    mat <- mat[subset, , drop = FALSE]
    
    if (abs) tval1 <- treat[1]
    else tval1 <- treat[binarize(treat)==1][1]
    
    v1 <- col.w.v(mat[treat==tval1, , drop = FALSE], weights[treat==tval1], bin.vars = bin.vars, na.rm = na.rm)
    v0 <- col.w.v(mat[treat!=tval1, , drop = FALSE], weights[treat!=tval1], bin.vars = bin.vars, na.rm = na.rm)
    
    v.ratios = v1/v0
    
    if (abs) v.ratios <- abs_(v.ratios, ratio = TRUE)
    
    return(setNames(v.ratios, colnames(mat)))
    
}
col_w_ks <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars, subset = NULL, na.rm = TRUE, ...) {
    
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(to.split <- vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(bin.vars, mat)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2",
                                           split.with = bin.vars),
                                      A))
        bin.vars <- attr(mat, "split.with")[[1]]
    }
    
    if (!(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
    
    possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
    
    weights <- weights * s.weights
    
    weights <- weights[subset]
    treat <- treat[subset]
    mat <- mat[subset, , drop = FALSE]
    
    tval1 <- treat[1]
    ks <- rep(NA_real_, NCOL(mat))
    
    if (any(!bin.vars)) {
        weights_ <- weights
        weights_[treat == tval1] <-  weights[treat == tval1]/sum(weights[treat == tval1])
        weights_[treat != tval1] <- -weights[treat != tval1]/sum(weights[treat != tval1])
        ks[!bin.vars] <- apply(mat[, !bin.vars, drop = FALSE], 2, function(x) {
            if (na.rm) x <- x[!is.na(x)]
            if (!na.rm && anyNA(x)) return(NA_real_)
            else {
                ordered.index <- order(x)
                cumv <- abs(cumsum(weights_[ordered.index]))[diff(x[ordered.index]) != 0]
                return(if (is_null(cumv)) 0 else max(cumv))
            }
        })
    }
    if (any(bin.vars)) {
        ks[bin.vars] <- abs(col.w.m(mat[treat == tval1, bin.vars, drop = FALSE], weights[treat == tval1], na.rm = na.rm) - 
                                col.w.m(mat[treat != tval1, bin.vars, drop = FALSE], weights[treat != tval1], na.rm = na.rm))
    }
    return(setNames(ks, colnames(mat)))
    
}
col_w_ovl <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars, integrate = FALSE, subset = NULL, na.rm = TRUE, ...) {
    
    A <- list(...)
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(to.split <- vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(bin.vars, mat)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        B <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2",
                                           split.with = bin.vars),
                                      B))
        bin.vars <- attr(mat, "split.with")[[1]]
    }
    
    if (!(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
    
    possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
    
    weights <- weights * s.weights
    
    weights <- weights[subset]
    treat <- treat[subset]
    mat <- mat[subset, , drop = FALSE]
    
    tval1 <- treat[1]
    
    t.sizes <- setNames(vapply(unique(treat, nmax = 2), function(x) sum(treat == x), numeric(1L)),
                        unique(treat, nmax = 2))
    smallest.t <- names(t.sizes)[which.min(t.sizes)]
    ovl <- setNames(numeric(ncol(mat)), colnames(mat))
    if (any(!bin.vars)) {
        if (is_null(A[["bw"]])) A[["bw"]] <- "nrd"
        A[names(A) %nin% names(formals(density.default))] <- NULL
        
        ovl[!bin.vars] <- apply(mat[, !bin.vars, drop = FALSE], 2, function(cov) {
            if (na.rm) cov <- cov[!is.na(cov)]
            if (!na.rm && anyNA(cov)) return(NA_real_)
            else {
                cov <- center(cov)/sd(cov)
                if (is.function(get0(paste0("bw.", A[["bw"]])))) {
                    A[["bw"]] <- get0(paste0("bw.", A[["bw"]]))(cov[treat == smallest.t])
                }
                else {
                    stop(paste(A[["bw"]], "is not an acceptable entry to bw. See ?stats::density for allowable options."), call. = FALSE)
                }
                
                f1_ <- approxfun(do.call(density.default, c(list(cov[treat==tval1], 
                                                                weights = weights[treat==tval1]/sum(weights[treat==tval1])), A)))
                f1 <- function(x) {
                    y <- f1_(x)
                    y[is.na(y)] <- 0
                    return(y)
                }
                f0_ <- approxfun(do.call(density.default, c(list(cov[treat!=tval1], 
                                                                weights = weights[treat!=tval1]/sum(weights[treat!=tval1])), A)))
                f0 <- function(x) {
                    y <- f0_(x)
                    y[is.na(y)] <- 0
                    return(y)
                }
                fn <- function(x) {
                    pmin(f1(x), f0(x))
                }
                min.c <- min(cov) - 4*A[["bw"]]
                max.c <- max(cov) + 4*A[["bw"]]
                # range <- max.c - min.c
                # min.c.ext <- min.c - .01 * range
                # max.c.ext <- max.c + .01 * range
                if (isTRUE(integrate)) {
                    s <- try(integrate(fn, lower = min.c,
                                       upper = max.c)$value,
                             silent = TRUE)
                }
                else {
                    seg <- seq(min.c, max.c, length = 1001)
                    mids <- .5 * (seg[2:length(seg)] + seg[1:(length(seg)-1)])
                    s <- sum(fn(mids))*(seg[2]-seg[1])
                }
                
                if (inherits(s, "try-error") || s > 1.2)  return(NA_real_)
                else return(1 - s) #Reverse: measure imbalance
            }
        })
    }
    if (any(bin.vars)) {
        ovl[bin.vars] <- abs(col.w.m(mat[treat == tval1, bin.vars, drop = FALSE], weights[treat == tval1]) - 
                                 col.w.m(mat[treat != tval1, bin.vars, drop = FALSE], weights[treat != tval1]))
    }
    
    return(ovl)
    
}
col_w_cov <- function(mat, treat, weights = NULL, type = "pearson", std = FALSE, s.d.denom = "all", abs = FALSE, s.weights = NULL, bin.vars, subset = NULL, weighted.weights = weights, na.rm = TRUE, ...) {
    
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(to.split <- vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(bin.vars, mat)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2",
                                           split.with = bin.vars),
                                      A))
        bin.vars <- attr(mat, "split.with")[[1]]
    }
    
    if (missing(treat) || !is.atomic(treat) || !is.numeric(treat)) stop("treat must be a numeric variable.")
    if (!is.atomic(std) || anyNA(as.logical(std)) ||
        length(std) %nin% c(1L, NCOL(mat))) {
        stop("std must be a logical vector with length equal to 1 or the number of columns of mat.")
    }
    if (!(is.atomic(abs) && length(abs) == 1L && !anyNA(as.logical(abs)))) {
        stop("abs must be a logical of length 1.")
    }
    
    possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
    
    if (length(std) == 1L) std <- rep(std, NCOL(mat))
    
    type <- match_arg(tolower(type), c("pearson", "spearman"))
    if (type == "spearman") {
        for (i in 1:ncol(mat)) if (!bin.vars[i]) mat[,i] <- rank(mat[,i], na.last = "keep")
        treat <- rank(treat, na.last = "keep")
    }
    
    weights <- weights * s.weights
    
    covars <- col.w.cov(mat[subset, , drop = FALSE], y = treat[subset], w = weights[subset], na.rm = na.rm)
    
    zeros <- check_if_zero(covars)
    
    if (any(to.sd <- std & !zeros)) {
        s.d.denom <- match_arg(s.d.denom, c("all", "weighted"))
        
        denoms <- rep(1, ncol(mat))
        
        if (s.d.denom == "all") denoms[to.sd] <- sqrt(col.w.v(mat[,to.sd, drop = FALSE], w = s.weights, bin.vars = bin.vars[to.sd], na.rm = na.rm)) * 
            sqrt(col.w.v(treat, w = s.weights, na.rm = na.rm))
        else if (s.d.denom == "weighted") denoms[to.sd] <- sqrt(col.w.v(mat[,to.sd, drop = FALSE], w = weighted.weights * s.weights, bin.vars = bin.vars[to.sd], na.rm = na.rm)) * 
            sqrt(col.w.v(treat, w = weighted.weights * s.weights, na.rm = na.rm))
        
        covars <- covars / denoms
    }
    
    if (abs) covars <- abs(covars)
    
    return(setNames(covars, colnames(mat)))
    
}
col_w_corr <- function(mat, treat, weights = NULL, type = "pearson", s.d.denom = "all", abs = FALSE, s.weights = NULL, bin.vars, subset = NULL, weighted.weights = weights, na.rm = TRUE, ...) {
    .call <- match.call(expand.dots = TRUE)
    .call[[1]] <- quote(col_w_cov)
    .call[["std"]] <- quote(TRUE)
    eval.parent(.call)
}