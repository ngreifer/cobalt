#Balance functions for use in programming

col_w_mean <- function(mat, weights = NULL, s.weights = NULL, subset = NULL, na.rm = TRUE, ...) {
    
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                A <- list(...)
                A <- A[names(A) %in% names(formals(splitfactor)) & 
                           names(A) %nin% c("data", "var.name", "drop.first",
                                            "drop.level")]
                mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                              A))
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
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
col_w_sd <- function(mat, weights = NULL, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
    
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                A <- list(...)
                A <- A[names(A) %in% names(formals(splitfactor)) & 
                           names(A) %nin% c("data", "var.name", "drop.first",
                                            "drop.level")]
                mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                              A))
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
    else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
             length(bin.vars) != NCOL(mat)) {
        stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
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

col_w_smd <- function(mat, treat, weights = NULL, std = TRUE, s.d.denom = "pooled", abs = FALSE, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
    allowable.s.d.denoms <- c("pooled", "all", "treated", "control")
    unique.treats <- as.character(unique(treat, nmax = 2))
    
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                A <- list(...)
                A <- A[names(A) %in% names(formals(splitfactor)) & 
                           names(A) %nin% c("data", "var.name", "drop.first",
                                            "drop.level")]
                mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                              A))
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    if (missing(treat) || !(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
    if (!is.atomic(std) || any(is.na(as.logical(std))) ||
        length(std) %nin% c(1L, NCOL(mat))) {
        stop("std must be a logical vector with length equal to 1 or the number of columns of mat.")
    }
    
    
    if (!(is.atomic(abs) && length(abs) == 1L && !is.na(as.logical(abs)))) {
        stop("abs must be a logical of length 1.")
    }
    
    if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
    else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
             length(bin.vars) != NCOL(mat)) {
        stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
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
    
    if (abs || !(length(s.d.denom) == 1L && as.character(s.d.denom) %in% c("treated", "control"))) tval1 <- treat[1]
    else tval1 <- treat[binarize(treat)==1][1]
    
    if (length(std) == 1L) std <- rep(std, NCOL(mat))
    
    m1 <- col.w.m(mat[treat==tval1 & subset, , drop = FALSE], weights[treat==tval1 & subset], na.rm = na.rm)
    m0 <- col.w.m(mat[treat!=tval1 & subset, , drop = FALSE], weights[treat!=tval1 & subset], na.rm = na.rm)
    diffs <- m1 - m0
    zeros <- check_if_zero(diffs)
    
    if (any(to.sd <- std & !zeros)) {
        if (is.atomic(s.d.denom) && length(s.d.denom) == 1L) {
            s.d.denom <- match_arg(as.character(s.d.denom), c(allowable.s.d.denoms, unique.treats))
            
            if (s.d.denom %in% unique.treats) s.d.denom <- sqrt(col.w.v(mat[treat == s.d.denom, to.sd, drop = FALSE], 
                                                                        w = s.weights[treat == s.d.denom],
                                                                        bin.vars = bin.vars[to.sd], na.rm = na.rm))
            else if (s.d.denom == "pooled") s.d.denom <- sqrt(.5*(col.w.v(mat[treat == tval1, to.sd, drop = FALSE], 
                                                                          w = s.weights[treat == tval1],
                                                                          bin.vars = bin.vars[to.sd], na.rm = na.rm) +
                                                                      col.w.v(mat[treat != tval1, to.sd, drop = FALSE], 
                                                                              w = s.weights[treat != tval1],
                                                                              bin.vars = bin.vars[to.sd], na.rm = na.rm)))
            else if (s.d.denom == "all") s.d.denom <- sqrt(col.w.v(mat[, to.sd, drop = FALSE], 
                                                                   w = s.weights,
                                                                   bin.vars = bin.vars[to.sd], na.rm = na.rm))
            else if (s.d.denom == "treated") s.d.denom <- sqrt(col.w.v(mat[treat == tval1, to.sd, drop = FALSE], 
                                                                       w = s.weights[treat == tval1],
                                                                       bin.vars = bin.vars[to.sd], na.rm = na.rm))
            else if (s.d.denom == "control") s.d.denom <- sqrt(col.w.v(mat[treat != tval1, to.sd, drop = FALSE], 
                                                                       w = s.weights[treat != tval1],
                                                                       bin.vars = bin.vars[to.sd], na.rm = na.rm))
        }
        else {
            if (!is.numeric(s.d.denom) || length(s.d.denom) %nin% c(sum(to.sd), length(diffs))) {
                stop(paste0("s.d.denom must be one of ", word_list(c(allowable.s.d.denoms, unique.treats), and.or = "or", quotes = TRUE), 
                            " or a numeric vector of with length equal to the number of columns of mat"))
            }
            else if (length(s.d.denom) == length(diffs)) s.d.denom <- s.d.denom[to.sd]
        }
        
        diffs[to.sd] <- diffs[to.sd]/s.d.denom
    }
    
    if (abs) diffs <- abs(diffs)
    
    return(setNames(diffs, colnames(mat)))
    
}
col_w_vr <- function(mat, treat, weights = NULL, abs = FALSE, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
    
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                A <- list(...)
                A <- A[names(A) %in% names(formals(splitfactor)) & 
                           names(A) %nin% c("data", "var.name", "drop.first",
                                            "drop.level")]
                mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                              A))
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    if (missing(treat) || !(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
    
    if (!(is.atomic(abs) && length(abs) == 1L && !is.na(as.logical(abs)))) {
        stop("abs must be a logical of length 1.")
    }
    
    if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
    else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
             length(bin.vars) != NCOL(mat)) {
        stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
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
col_w_ks <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
    
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                A <- list(...)
                A <- A[names(A) %in% names(formals(splitfactor)) & 
                           names(A) %nin% c("data", "var.name", "drop.first",
                                            "drop.level")]
                mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                              A))
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
    else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
             length(bin.vars) != NCOL(mat)) {
        stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
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
col_w_ovl <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars = NULL, integrate = FALSE, subset = NULL, na.rm = TRUE, ...) {
    
    A <- list(...)
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                A <- list(...)
                A <- A[names(A) %in% names(formals(splitfactor)) & 
                           names(A) %nin% c("data", "var.name", "drop.first",
                                            "drop.level")]
                mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                              A))
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
    else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
             length(bin.vars) != NCOL(mat)) {
        stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
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
        if (is_null(A[["bw"]])) A[["bw"]] <- "nrd" else A[["bw"]] <- A[["bw"]]
        A[names(A) %nin% names(formals(density.default))] <- NULL
        
        ovl[!bin.vars] <- apply(mat[, !bin.vars, drop = FALSE], 2, function(cov) {
            if (na.rm) cov <- cov[!is.na(cov)]
            if (!na.rm && anyNA(cov)) return(NA_real_)
            else {
                if (is.function(get0(paste0("bw.", A[["bw"]])))) {
                    A[["bw"]] <- get0(paste0("bw.", A[["bw"]]))(cov[treat == smallest.t])
                }
                else {
                    stop(paste(A[["bw"]], "is not an acceptable entry to bw. See ?stats::density for allowable options."), call. = FALSE)
                }
                
                f1 <- approxfun(do.call(density.default, c(list(cov[treat==tval1], 
                                                                weights = weights[treat==tval1]/sum(weights[treat==tval1])), A)))
                f0 <- approxfun(do.call(density.default, c(list(cov[treat!=tval1], 
                                                                weights = weights[treat!=tval1]/sum(weights[treat!=tval1])), A)))
                fn <- function(x) {
                    f1_ <- f1(x)
                    f1_[is.na(f1_)] <- 0
                    f0_ <- f0(x)
                    f0_[is.na(f0_)] <- 0
                    pmin(f1_, f0_)
                }
                min.c <- min(cov)
                max.c <- max(cov)
                range <- max.c - min.c
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
                
                if (inherits(s, "try-error"))  return(NA_real_)
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
col_w_corr <- function(mat, treat, weights = NULL, type = "pearson", abs = FALSE, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
    
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                A <- list(...)
                A <- A[names(A) %in% names(formals(splitfactor)) & 
                           names(A) %nin% c("data", "var.name", "drop.first",
                                            "drop.level")]
                mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                              A))
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    if (missing(treat) || !is.atomic(treat) || !is.numeric(treat)) stop("treat must be a numeric variable.")
    if (!(is.atomic(abs) && length(abs) == 1L && !is.na(as.logical(abs)))) {
        stop("abs must be a logical of length 1.")
    }
    if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
    else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
             length(bin.vars) != NCOL(mat)) {
        stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
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
    
    type <- match_arg(tolower(type), c("pearson", "spearman"))
    if (type == "spearman") {
        for (i in 1:ncol(mat)) mat[,i] <- rank(mat[,i], na.last = "keep")
        treat <- rank(treat, na.last = "keep")
    }
    
    weights <- weights * s.weights
    
    if (all(subset)) {
        corrs <- col.w.r(mat, treat, w = weights, s.weights = s.weights, bin.vars = bin.vars,
                         na.rm = na.rm)
    }
    else {
        cov <- col.w.cov(mat[subset, , drop = FALSE], y = treat[subset], w = weights[subset], na.rm = na.rm)
        den <- sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm)) * 
            sqrt(col.w.v(treat, w = s.weights, na.rm = na.rm))
    }
    
    if (abs) corrs <- abs(corrs)
    
    return(setNames(corrs, colnames(mat)))
    
}

#Scalar balance functions 
bal.sum <- function(mat, treat, weights = NULL, type, s.weights = NULL, check = TRUE, ...) {
    uni.type <- c("smd", "ks", "ovl")
    agg.funs <- c("max", "mean", "rms")
    sample.type <- c("mahalanobis", "gwd", "cstat", "wr2", "design.effect")
    uni.type.expanded <- expand.grid_string(uni.type, agg.funs, collapse = ".")
    shortcuts <- c("all", "rec")
    allowable.type <- c(uni.type.expanded, sample.type, shortcuts)
    if (missing(type)) stop("type must be specified.", call. = FALSE)
    else type <- match_arg(type, allowable.type, several.ok = TRUE)
    
    if ("all" %in% type) type <- unique(c(type[type != "all"], uni.type.expanded, sample.type))
    if ("rec" %in% type) type <- unique(c(type[type != "rec"], 
                                          "smd.mean", "smd.rms",
                                          "ks.mean", "ks.rms",
                                          "mahalanobis", "gwd", "wr2"))
    
    A <- list(...)
    
    if (check) {
        bad.mat <- FALSE
        if (missing(mat)) bad.mat <- TRUE
        else {
            if (is.data.frame(mat)) {
                if (any(vapply(mat, function(x) is_(x, c("character", "factor")), logical(1L)))) 
                    mat <- splitfactor(mat)
                mat <- as.matrix.data.frame(mat)
            }
            else if (is.vector(mat, "numeric")) mat <- matrix(mat, ncol = 1)
            else if (!is.matrix(mat) || !is.numeric(mat)) bad.mat <- TRUE
        }
        if (bad.mat) stop("mat must be a numeric matrix.")
        
        if (missing(treat) || !(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
        
        if (is_null(weights)) weights <- rep(1, NROW(mat))
        if (is_null(s.weights)) s.weights <- rep(1, NROW(mat))
        if (!all_the_same(c(NROW(mat), length(treat), length(weights), length(s.weights)))) {
            stop("mat, treat, weights, and s.weights (if supplied) must have the same number of units.")
        }
    }
    
    if (any(paste.("smd", agg.funs) %in% type)) {
        if (!exists("s.d.denom")) s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = A[["estimand"]], weights = data.frame(weights), treat = treat)
        smd <- col_w_smd(mat, treat, weights, std = TRUE, s.d.denom, abs = TRUE, s.weights = s.weights, check = FALSE)
        if (is_null(A[["smd.weights"]])) smd.weights <- rep(1, ncol(mat))
        else if (!is.vector(A[["smd.weights"]], "numeric")) {
            warning("smd.weights is not numeric. Ignoring smd.weights", 
                    call. = FALSE, immediate. = TRUE)
            smd.weights <- rep(1, ncol(mat))
        }
        else if (length(A[["smd.weights"]]) == ncol(mat)) {
            smd.weights <- A[["smd.weights"]]
        }
        else {
            warning("smd.weights should be of length ncol(mat). Ignoring smd.weights", 
                    call. = FALSE, immediate. = TRUE)
            smd.weights <- rep(1, ncol(mat))
        }
        smd.weights <- smd.weights/mean(smd.weights) #Make sum to ncol(mat)
        smd <- smd.weights*smd
    }
    if (any(paste.("ks", agg.funs) %in% type)) {
        ks <- col_w_ks(mat, treat, weights, bin.vars = A[["bin.vars"]], check = is_null(A[["bin.vars"]]))
    }
    if (any(paste.("ovl", agg.funs) %in% type)) {
        ovl <- do.call(col_w_ovl, c(list(mat = mat, treat = treat, weights = weights, check = is_null(A[["bin.vars"]])), A))
    }
    if ("gwd" %in% type) {
        if (!exists("s.d.denom")) s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = A[["estimand"]], weights = data.frame(weights), treat = treat)
    }
    
    bal <- setNames(vapply(type, function(m) {
        if (endsWith(m, ".mean")) {
            agg <- mean
            m <- substr(m, 1, nchar(m) - nchar(".mean"))
        }
        else if (endsWith(m, ".max")) {
            agg <- max
            m <- substr(m, 1, nchar(m) - nchar(".max"))
        }
        else if (endsWith(m, ".rms")) {
            agg <- function(x, ...) {sqrt(mean(x^2, ...))}
            m <- substr(m, 1, nchar(m) - nchar(".rms"))
        }
        
        if (m %in% uni.type) {
            return(agg(get0(m), na.rm = TRUE))
        }
        else if (m == "mahalanobis") {
            if (is_null(s.weights)) s.weights <- rep(1, nrow(mat))
            mdiff <- matrix(col_w_smd(mat, treat, weights, std = FALSE, abs = FALSE, check = FALSE), ncol = 1)
            wcov <- cov.wt(mat, s.weights)$cov
            mahal <- crossprod(mdiff, solve(wcov)) %*% mdiff
            return(mahal)
        }
        else if (m == "cstat") {
            tval1 <- treat[1]
            d <- data.frame(treat, mat)
            f <- formula(d)
            pred <- glm(f, data = d, family = quasibinomial(),
                        weights = weights)$fitted
            wi <- wilcox.test(pred ~ treat)
            cstat <- wi$statistic/(sum(treat==tval1)*sum(treat!=tval1))
            cstat <- 2*max(cstat, 1-cstat)-1
            return(cstat)
        }
        else if (m == "wr2") {
            tval1 <- treat[1]
            d <- data.frame(treat, mat)
            f <- formula(d)
            fit <- glm(f, data = d, family = quasibinomial(),
                       weights = weights)
            r2 <- 1 - (pi^2/3)/(3*var(fit$linear.predictors) + pi^2/3)
            return(r2)
        }
        else if (m == "gwd") {
            co.names <- setNames(lapply(colnames(mat), function(x) setNames(list(x, "base"), c("component", "type"))), colnames(mat))
            new <- int.poly.f(mat, int = TRUE, poly = 2, center = isTRUE(A[["center"]]), 
                              sep = getOption("cobalt_int_sep", default = " * "), co.names = co.names)
            mat_ <- cbind(mat, new)
            
            smd <- col_w_smd(mat_, treat, weights, std = TRUE, s.d.denom, abs = TRUE, check = FALSE)
            if (is_null(A[["gwd.weights"]])) gwd.weights <- c(rep(1, ncol(mat)), rep(.5, ncol(new)))
            else if (!is.vector(A[["gwd.weights"]], "numeric")) {
                warning("gwd.weights is not numeric. Ignoring gwd.weights.", 
                        call. = FALSE, immediate. = TRUE)
                gwd.weights <- c(rep(1, ncol(mat)), rep(.5, ncol(new)))
            }
            else if (length(A[["gwd.weights"]]) == 1L) {
                gwd.weights <- rep(A[["gwd.weights"]], ncol(mat_))
            }
            else if (length(A[["gwd.weights"]]) == 2L) {
                gwd.weights <- c(rep(A[["gwd.weights"]][1], ncol(mat)), rep(A[["gwd.weights"]][2], ncol(new)))
            }
            else {
                warning("gwd.weights should be of length 1 or 2. Ignoring gwd.weights.", 
                        call. = FALSE, immediate. = TRUE)
                gwd.weights <- c(rep(1, ncol(mat)), rep(.5, ncol(new)))
            }
            
            gwd.weights <- gwd.weights/sum(gwd.weights) #Make sum to 1
            return(sum(gwd.weights*smd, na.rm = TRUE))
        }
        else if (m == "design.effect") {
            tval1 <- treat[1]
            q <- sum(treat == tval1)/length(treat)
            des.eff <- function(w) length(w)*sum(w^2)/sum(w)^2
            des <- c(des.eff(weights[treat == tval1]), des.eff(weights[treat != tval1]))
            de <- des[1]*(1-q) + des[2]*q
            return(de)
        }
    }, numeric(1L)), type)
    
    return(bal)
}