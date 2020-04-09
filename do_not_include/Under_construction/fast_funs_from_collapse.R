#Fast stats from collapse package

if (check.package("collapse", TRUE)) {
    w.m <- collapse::fmean
    col.w.m <- function(mat, w = NULL, na.rm = TRUE, g = NULL) {
        collapse::fmean(x = mat, w = w, na.rm = na.rm, g = g)
    }
    col.w.v <- function(mat, w = NULL, bin.vars = NULL, na.rm = TRUE) {
        
        if (!is.matrix(mat)) {
            if (is.data.frame(mat)) {
                if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                    stop("mat must be a numeric matrix.")
                }
                else mat <- as.matrix.data.frame(mat)
            }
            else if (is.numeric(mat)) {
                mat <- matrix(mat, ncol = 1)
            }
            else stop("mat must be a numeric matrix.")
        }
        
        if (is_null(bin.vars)) bin.vars <- rep(FALSE, ncol(mat))
        else if (length(bin.vars) != ncol(mat) || anyNA(as.logical(bin.vars))) {
            stop("bin.vars must be a logical vector with length equal to the number of columns of mat.", call. = FALSE)
        }
        bin.var.present <- any(bin.vars)
        non.bin.vars.present <- any(!bin.vars)
        
        if (is_null(w)) {
            var <- collapse::fvar(mat, na.rm = na.rm)
            if (bin.var.present) {
                n <- collapse::fNobs(mat)
                var[bin.vars] <- var[bin.vars]*(n[bin.vars]-1)/n[bin.vars]
            }
        }
        else {
            var <- collapse::fvar(mat, w = w, na.rm = na.rm)
            sumw <- collapse::fsum((1-is.na(mat))*w, na.rm = na.rm)
            if (bin.var.present) {
                var[bin.vars] <- var[bin.vars]*(sumw[bin.vars]-1)/sumw[bin.vars]
            }
            if (non.bin.vars.present) {
                sumw2 <- collapse::fsum((1-is.na(mat[,!bin.vars]))*w^2, na.rm = na.rm)
                var[!bin.vars] <- var[!bin.vars]*(sumw[!bin.vars]-1)*sumw[!bin.vars]/(sumw[!bin.vars]^2-sumw2)
            }
        }
        
        return(var)
    }
    col.w.cov <- function(mat, y, w = NULL, na.rm = TRUE) {
        if (!is.matrix(mat)) {
            if (is_null(w)) return(cov(mat, y, use = if (na.rm) "pair" else "everything"))
            else mat <- matrix(mat, ncol = 1)
        }
        if (is_null(w)) {
            y <- array(y, dim = dim(mat))
            y[is.na(mat)] <- NA
            mat[is.na(y)] <- NA
            den <- collapse::fNobs(mat*y) - 1
            cov <- collapse::fsum(collapse::fwithin(mat, na.rm = na.rm)*collapse::fwithin(y, na.rm = na.rm), na.rm = na.rm)/den
        }
        else if (na.rm && anyNA(mat)) {
            y <- array(y, dim = dim(mat))
            y[is.na(mat)] <- NA
            mat[is.na(y)] <- NA
            w <- array(w, dim = dim(mat))
            w[is.na(mat)] <- NA_real_
            w_ <- collapse::fsum(w, na.rm = na.rm, TRA = "/")
            x <- w_ * center(mat, at = collapse::fsum(w_ * mat, na.rm = na.rm))
            cov <- collapse::fsum(x*y, na.rm = na.rm)/(1 - collapse::fsum(w_^2, na.rm = na.rm))
        }
        else {
            w_ <- w/sum(w)
            cov <- collapse::fmean(collapse::fwithin(mat, w = w_)*array(collapse::fwithin(y, w = w_), dim = dim(mat)), w = w_)/(1 - sum(w_^2))
        }
        return(cov)
    }
    
    nunique <- collapse::fNdistinct
    nunique.gt <- function(x, n, na.rm = TRUE) {
        nunique(x, na.rm = na.rm) > n
    }
    all_the_same <- function(x, na.rm = TRUE) {
        nunique(x, na.rm = na.rm) == 1
    }
    is_binary <- function(x, na.rm = TRUE) {
        nunique(x, na.rm = na.rm) == 2
    }
    is_binary_col <- function(dat, na.rm = TRUE) {
        if (length(dim(dat)) != 2) warning("is_binary_col shouldn't be used with objects that don't have 2 dimensions.")
        collapse::fNdistinct(dat, na.rm = na.rm) == 2
    }
    
    unique <- function(x, ...) collapse::funique(x, FALSE)
} else {
    w.m <- function(x, w = NULL, na.rm = TRUE) {
        if (is_null(w)) w <- rep(1, length(x))
        w[is.na(x)] <- NA_real_
        return(sum(x*w, na.rm=na.rm)/sum(w, na.rm=na.rm))
    }
    col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
        if (is_null(w)) w <- 1
        w.sum <- colSums(w*!is.na(mat))
        return(colSums(mat*w, na.rm = na.rm)/w.sum)
    }
    col.w.v <- function(mat, w = NULL, bin.vars = NULL, na.rm = TRUE) {
        if (!is.matrix(mat)) {
            if (is.data.frame(mat)) {
                if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                    stop("mat must be a numeric matrix.")
                }
                else mat <- as.matrix.data.frame(mat)
            }
            else if (is.numeric(mat)) {
                mat <- matrix(mat, ncol = 1)
            }
            else stop("mat must be a numeric matrix.")
        }
        
        if (is_null(bin.vars)) bin.vars <- rep(FALSE, ncol(mat))
        else if (length(bin.vars) != ncol(mat) || anyNA(as.logical(bin.vars))) {
            stop("bin.vars must be a logical vector with length equal to the number of columns of mat.", call. = FALSE)
        }
        bin.var.present <- any(bin.vars)
        non.bin.vars.present <- any(!bin.vars)
        
        var <- setNames(numeric(ncol(mat)), colnames(mat))
        if (is_null(w)) {
            if (non.bin.vars.present) {
                den <- colSums(!is.na(mat[, !bin.vars, drop = FALSE])) - 1
                var[!bin.vars] <- colSums(center(mat[, !bin.vars, drop = FALSE])^2, na.rm = na.rm)/den
            }
            if (bin.var.present) {
                means <- colMeans(mat[, bin.vars, drop = FALSE], na.rm = na.rm)
                var[bin.vars] <- means * (1 - means)
            }
        }
        else if (na.rm && anyNA(mat)) {
            # n <- nrow(mat)
            w <- array(w, dim = dim(mat))
            w[is.na(mat)] <- NA_real_
            s <- colSums(w, na.rm = na.rm)
            w <- mat_div(w, s)
            if (non.bin.vars.present) {
                x <- sqrt(w[, !bin.vars, drop = FALSE]) * center(mat[, !bin.vars, drop = FALSE],
                                                                 at = colSums(w[, !bin.vars, drop = FALSE] * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
                var[!bin.vars] <- colSums(x*x, na.rm = na.rm)/(1 - colSums(w[, !bin.vars, drop = FALSE]^2, na.rm = na.rm))
            }
            if (bin.var.present) {
                means <- colSums(w[, bin.vars, drop = FALSE] * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
                var[bin.vars] <- means * (1 - means)
            }
        }
        else {
            if (is_null(w)) w <- rep(1, nrow(mat))
            w <- w/sum(w)
            if (non.bin.vars.present) {
                x <- sqrt(w) * center(mat[, !bin.vars, drop = FALSE],
                                      at = colSums(w * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
                var[!bin.vars] <- colSums(x*x, na.rm = na.rm)/(1 - sum(w^2))
            }
            if (bin.var.present) {
                means <- colSums(w * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
                var[bin.vars] <- means * (1 - means)
            }
        }
        return(var)
    }
    col.w.cov <- function(mat, y, w = NULL, na.rm = TRUE) {
        if (!is.matrix(mat)) {
            if (is_null(w)) return(cov(mat, y, use = if (na.rm) "pair" else "everything"))
            else mat <- matrix(mat, ncol = 1)
        }
        if (is_null(w)) {
            y <- array(y, dim = dim(mat))
            y[is.na(mat)] <- NA
            mat[is.na(y)] <- NA
            den <- colSums(!is.na(mat*y)) - 1
            cov <- colSums(center(mat, na.rm = na.rm)*center(y, na.rm = na.rm), na.rm = na.rm)/den
        }
        else if (na.rm && anyNA(mat)) {
            n <- nrow(mat)
            w <- array(w, dim = dim(mat))
            w[is.na(mat)] <- NA_real_
            s <- colSums(w, na.rm = na.rm)
            w <- mat_div(w, s)
            x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
            cov <- colSums(x*y, na.rm = na.rm)/(1 - colSums(w^2, na.rm = na.rm))
        }
        else {
            n <- nrow(mat)
            w <- w/sum(w)
            x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
            cov <- colSums(x*y, na.rm = na.rm)/(1 - sum(w^2))
        }
        return(cov)
    }
    
    nunique <- function(x, nmax = NA, na.rm = TRUE) {
        if (is_null(x)) return(0)
        else {
            if (na.rm) x <- x[!is.na(x)]
            if (is.factor(x)) return(nlevels(x))
            else return(length(unique(x, nmax = nmax)))
        }
        
    }
    nunique.gt <- function(x, n, na.rm = TRUE) {
        if (missing(n)) stop("n must be supplied.")
        if (n < 0) stop("n must be non-negative.")
        if (is_null(x)) FALSE
        else {
            if (na.rm) x <- x[!is.na(x)]
            if (n == 1) !all_the_same(x)
            else if (length(x) < 2000) nunique(x) > n
            else tryCatch(nunique(x, nmax = n) > n, error = function(e) TRUE)
        }
    }
    all_the_same <- function(x, na.rm = TRUE) {
        if (na.rm && anyNA(x)) x <- x[!is.na(x)]
        if (is.double(x)) check_if_zero(abs(max_(x) - min_(x)))
        else !any(x != x[1])
    }
    is_binary <- function(x, na.rm = TRUE) {
        if (na.rm && anyNA(x)) x <- x[!is.na(x)]
        !all_the_same(x) && all_the_same(x[x != x[1]])
    }
    is_binary_col <- function(dat, na.rm = TRUE) {
        if (length(dim(dat)) != 2) stop("is_binary_col cannot be used with objects that don't have 2 dimensions.")
        apply(dat, 2, is_binary)
    }
}
