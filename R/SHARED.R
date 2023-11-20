#This document is shared across cobalt, WeightIt, and optweight

#Strings
word_list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
    #When given a vector of strings, creates a string of the form "a and b"
    #or "a, b, and c"
    #If is.are, adds "is" or "are" appropriately
    L <- length(word.list)
    word.list <- add_quotes(word.list, quotes)
    
    if (L == 0) {
        out <- ""
        attr(out, "plural") <- FALSE
        return(out)
    }
    
    word.list <- word.list[!word.list %in% c(NA_character_, "")]
    L <- length(word.list)
    if (L == 0) {
        out <- ""
        attr(out, "plural") <- FALSE
    }
    else if (L == 1) {
        out <- word.list
        if (is.are) out <- paste(out, "is")
        attr(out, "plural") <- FALSE
    }
    else {
        and.or <- match_arg(and.or)
        if (L == 2) {
            out <- paste(word.list, collapse = paste0(" ", and.or, " "))
        }
        else {
            out <- paste(paste(word.list[seq_len(L - 1)], collapse = ", "),
                         word.list[L], sep = paste0(", ", and.or, " "))
            
        }
        if (is.are) out <- paste(out, "are")
        attr(out, "plural") <- TRUE
    }
    
    out
}
add_quotes <- function(x, quotes = 2L) {
    if (isFALSE(quotes)) return(x)
    
    if (isTRUE(quotes)) quotes <- 2
    
    if (chk::vld_string(quotes)) x <- paste0(quotes, x, quotes)
    else if (chk::vld_whole_number(quotes)) {
        if (as.integer(quotes) == 0) return(x)
        
        if (as.integer(quotes) == 1) x <- paste0("\'", x, "\'")
        else if (as.integer(quotes) == 2) x <- paste0("\"", x, "\"")
        else stop("`quotes` must be boolean, 1, 2, or a string.")
    }
    else {
        stop("`quotes` must be boolean, 1, 2, or a string.")
    }
    
    x
}
firstup <- function(x) {
    #Capitalize first letter
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}
expand.grid_string <- function(..., collapse = "") {
    do.call("paste", c(expand.grid(...), sep = collapse))
}
num_to_superscript <- function(x) {
    nums <- setNames(c("\u2070",
                       "\u00B9",
                       "\u00B2",
                       "\u00B3",
                       "\u2074",
                       "\u2075",
                       "\u2076",
                       "\u2077",
                       "\u2078",
                       "\u2079"),
                     as.character(0:9))
    x <- as.character(x)
    splitx <- strsplit(x, "", fixed = TRUE)
    
    sapply(splitx, function(y) paste0(nums[y], collapse = ""))
}
ordinal <- function(x) {
    if (!is.numeric(x) || !is.vector(x) || is_null(x)) stop("'x' must be a numeric vector.")
    if (length(x) > 1) return(vapply(x, ordinal, character(1L)))
    
    x0 <- abs(x)
    out <- paste0(x0, switch(substring(x0, nchar(x0), nchar(x0)),
                             "1" = "st",
                             "2" = "nd",
                             "3" = "rd",
                             "th"))
    if (sign(x) == -1) out <- paste0("-", out)
    
    out
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
    if (NROW(df) == 0 || NCOL(df) == 0) return(df)
    if (!is.data.frame(df)) df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
    rn <- rownames(df)
    cn <- colnames(df)
    
    infs <- o.negs <- array(FALSE, dim = dim(df))
    nas <- is.na(df)
    nums <- vapply(df, is.numeric, logical(1))
    infs[,nums] <- vapply(which(nums), function(i) !nas[,i] & !is.finite(df[[i]]), logical(NROW(df)))
    
    for (i in which(!nums)) {
        if (can_str2num(df[[i]])) {
            df[[i]] <- str2num(df[[i]])
            nums[i] <- TRUE
        }
    }
    
    o.negs[,nums] <- !nas[,nums] & df[nums] < 0 & round(df[nums], digits) == 0
    df[nums] <- round(df[nums], digits = digits)
    
    for (i in which(nums)) {
        df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                          drop0trailing = !identical(as.character(pad), "0"))
        
        if (!identical(as.character(pad), "0") && any(grepl(".", df[[i]], fixed = TRUE))) {
            s <- strsplit(df[[i]], ".", fixed = TRUE)
            lengths <- lengths(s)
            digits.r.of.. <- rep(0, NROW(df))
            digits.r.of..[lengths > 1] <- nchar(vapply(s[lengths > 1], `[[`, character(1L), 2))
            max.dig <- max(digits.r.of..)
            
            dots <- ifelse(lengths > 1, "", if (as.character(pad) != "") "." else pad)
            pads <- vapply(max.dig - digits.r.of.., function(n) paste(rep(pad, n), collapse = ""), character(1L))
            
            df[[i]] <- paste0(df[[i]], dots, pads)
        }
    }
    
    df[o.negs] <- paste0("-", df[o.negs])
    
    # Insert NA placeholders
    df[nas] <- na_vals
    df[infs] <- "N/A"
    
    if (length(rn) > 0) rownames(df) <- rn
    if (length(cn) > 0) names(df) <- cn
    
    df
}
text_box_plot <- function(range.list, width = 12) {
    full.range <- range(unlist(range.list))
    ratio = diff(full.range)/(width+1)
    rescaled.range.list <- lapply(range.list, function(x) round(x/ratio))
    rescaled.full.range <- round(full.range/ratio)
    d <- make_df(c("Min", paste(rep(" ", width + 1), collapse = ""), "Max"),
                 names(range.list),
                 "character")
    d[["Min"]] <- vapply(range.list, function(x) x[1], numeric(1L))
    d[["Max"]] <- vapply(range.list, function(x) x[2], numeric(1L))
    for (i in seq_len(nrow(d))) {
        spaces1 <- rescaled.range.list[[i]][1] - rescaled.full.range[1]
        #|
        dashes <- max(c(0, diff(rescaled.range.list[[i]]) - 2))
        #|
        spaces2 <- max(c(0, diff(rescaled.full.range) - (spaces1 + 1 + dashes + 1)))
        
        d[i, 2] <- paste0(paste(rep(" ", spaces1), collapse = ""),
                          "|",
                          paste(rep("-", dashes), collapse = ""),
                          "|",
                          paste(rep(" ", spaces2), collapse = ""))
    }
    
    d
}

equivalent.factors2 <- function(f1, f2) {
    if (!chk::vld_character_or_factor(f1) ||
        !chk::vld_character_or_factor(f2)) {
        return(FALSE)
    }
    
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    ll1 <- levels(f1)
    ll2 <- levels(f2)
    f1 <- as.integer(f1)
    f2 <- as.integer(f2)
    nl1 <- length(ll1)
    nl2 <- length(ll2)
    
    dims <- c(nl1, nl2)
    dn <- list(ll1, ll2)
    
    bin <- f1 + nl1 * (f2 - 1L)
    pd <- nl1 * nl2
    
    tab_ <- array(tabulate(bin, pd), dims, dimnames = dn)
    
    all(colSums(tab_ != 0) %in% 0:1) && all(rowSums(tab_ != 0) %in% 0:1)
}

paste. <- function(..., collapse = NULL) {
    #Like paste0 but with sep = ".'
    paste(..., sep = ".", collapse = collapse)
}
wrap <- function(s, nchar, ...) {
    vapply(s, function(s_) {
        x <- strwrap(s_, width = nchar, ...)
        paste(x, collapse = "\n")
    }, character(1L))
}
strsplits <- function(x, splits, fixed = TRUE, ...) {
    #Link strsplit but takes multiple split values.
    #Only works for one string at a time (in x).
    for (split in splits) x <- unlist(strsplit(x, split, fixed = TRUE, ...))
    
    x[x != ""] # Remove empty values
}
c.factor <- function(..., recursive=TRUE) {
    #c() for factors
    unlist(list(...), recursive=recursive)
}
can_str2num <- function(x) {
    nas <- is.na(x)
    suppressWarnings(x_num <- as.numeric(as.character(x[!nas])))
    
    !anyNA(x_num)
}
str2num <- function(x) {
    nas <- is.na(x)
    suppressWarnings(x_num <- as.numeric(as.character(x)))
    x_num[nas] <- NA
    
    x_num
}
trim_string <- function(x, char = " ", symmetrical = TRUE, recursive = TRUE) {
    if (is_null(x)) return(x)
    
    sw <- startsWith(x, char)
    ew <- endsWith(x, char)
    
    if (symmetrical) {
        if (!any(sw & ew)) return(x)
        x[sw & ew] <- gsub('^.|.$', '', x[sw & ew])
    }
    else {
        asw <- any(sw)
        aew <- any(ew)
        if (!asw && !aew) return(x)
        
        if (asw) x[sw] <- gsub('^.', '', x[sw])
        if (aew) x[ew] <- gsub('.$', '', x[ew])
    }
    
    if (!recursive) return(x)
    
    trim_string(x, char, symmetrical, recursive)
}

#Numbers
check_if_zero <- function(x) {
    # this is the default tolerance used in all.equal
    tolerance <- .Machine$double.eps^0.5
    abs(x) < tolerance
}
between <- function(x, range, inclusive = TRUE, na.action = FALSE) {
    if (!all(is.numeric(x))) stop("'x' must be a numeric vector.", call. = FALSE)
    if (length(range) != 2) stop("'range' must be of length 2.", call. = FALSE)
    if (anyNA(range) || !is.numeric(range)) stop("'range' must contain numeric entries only.", call. = FALSE)
    
    if (range[2] < range[1]) range <- c(range[2], range[1])
    
    if (anyNA(x)) {
        if (length(na.action) != 1 || !is.atomic(na.action))
            stop("'na.action' must be an atomic vector of length 1.", call. = FALSE)
    }
    
    if (inclusive) out <- ifelse(is.na(x), na.action, x >= range[1] & x <= range[2])
    else out <- ifelse(is.na(x), na.action, x > range[1] & x < range[2])
    
    out
}
max_ <- function(..., na.rm = TRUE) {
    if (!any(is.finite(unlist(list(...))))) NA_real_
    else max(..., na.rm = na.rm)
}
min_ <- function(..., na.rm = TRUE) {
    if (!any(is.finite(unlist(list(...))))) NA_real_
    else min(..., na.rm = na.rm)
}
check_if_int <- function(x) {
    #Checks if integer-like
    if (is.integer(x)) rep(TRUE, length(x))
    else if (is.numeric(x)) check_if_zero(x - round(x))
    else rep(FALSE, length(x))
}

#Statistics
binarize <- function(variable, zero = NULL, one = NULL) {
    if (!is_binary(variable)) stop(paste0("Cannot binarize ", deparse1(substitute(variable)), ": more than two levels."))
    if (is.character(variable) || is.factor(variable)) {
        variable <- factor(variable, nmax = 2)
        unique.vals <- levels(variable)
    }
    else {
        unique.vals <- unique(variable, nmax = 2)
    }
    
    if (is_not_null(zero)) {
        if (zero %nin% unique.vals) {
            stop("The argument to 'zero' is not the name of a level of variable.")
        }
        return(setNames(as.integer(variable != zero), names(variable)))
    }
    
    if (is_not_null(one)) {
        if (one %nin% unique.vals) {
            stop("The argument to 'one' is not the name of a level of variable.")
        }
        return(setNames(as.integer(variable == one), names(variable)))
    }
    
    variable.numeric <- {
        if (can_str2num(unique.vals)) str2num(variable)
        else as.numeric(variable)
    }
    
    zero <- {
        if (0 %in% variable.numeric) 0
        else min(variable.numeric, na.rm = TRUE)
    }
    
    setNames(as.integer(variable.numeric != zero), names(variable))

}
ESS <- function(w) {
    sum(w)^2/sum(w^2)
}
center <- function(x, at = NULL, na.rm = TRUE) {
    if (is.data.frame(x)) {
        x <- as.matrix.data.frame(x)
        type <- "df"
    }
    
    if (!is.numeric(x)) stop("'x' must be numeric.")
    if (is.array(x) && length(dim(x)) > 2) stop("'x' must be a numeric or matrix-like (not array).")
    
    if (!is.matrix(x)) {
        x <- matrix(x, ncol = 1)
        type <- "vec"
    }
    else type <- "matrix"
    
    if (is_null(at)) at <- colMeans(x, na.rm = na.rm)
    else if (length(at) %nin% c(1, ncol(x))) stop("'at' is not the right length.")
    
    out <- x - matrix(at, byrow = TRUE, ncol = ncol(x), nrow = nrow(x))
    
    if (type == "df") out <- as.data.frame.matrix(out)
    else if (type == "vec") out <- drop(out)
    
    out
}
w.m <- function(x, w = NULL, na.rm = TRUE) {
    if (is_null(w)) w <- rep(1, length(x))
    if (anyNA(x)) is.na(w)[is.na(x)] <- TRUE
    
    sum(x * w, na.rm = na.rm)/sum(w, na.rm = na.rm)
}
col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
    if (is_null(w)) w <- rep(1, nrow(mat))
    
    w.sum <- {
        if (!na.rm || !anyNA(mat)) sum(w)
        else colSums(w*!is.na(mat))
    }
    
    colSums(mat*w, na.rm = na.rm)/w.sum
}
col.w.v <- function(mat, w = NULL, bin.vars = NULL, na.rm = TRUE) {
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                stop("'mat' must be a numeric matrix.")
            }
            mat <- data.matrix(mat)
        }
        else if (is.numeric(mat)) {
            mat <- matrix(mat, ncol = 1)
        }
        else stop("'mat' must be a numeric matrix.")
    }
    
    if (is_null(bin.vars)) bin.vars <- rep(FALSE, ncol(mat))
    else if (length(bin.vars) != ncol(mat) || anyNA(as.logical(bin.vars))) {
        stop("'bin.vars' must be a logical vector with length equal to the number of columns of 'mat'.", call. = FALSE)
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
        is.na(w)[is.na(mat)] <- TRUE
        s <- colSums(w, na.rm = na.rm)
        w <- mat_div(w, s)
        if (non.bin.vars.present) {
            x <- center(mat[, !bin.vars, drop = FALSE],
                        at = colSums(w[, !bin.vars, drop = FALSE] * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
            var[!bin.vars] <- colSums(w[, !bin.vars, drop = FALSE]*x*x, na.rm = na.rm)/(1 - colSums(w[, !bin.vars, drop = FALSE]^2, na.rm = na.rm))
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
            x <- center(mat[, !bin.vars, drop = FALSE],
                        at = colSums(w * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
            var[!bin.vars] <- colSums(w*x*x, na.rm = na.rm)/(1 - sum(w^2))
        }
        if (bin.var.present) {
            means <- colSums(w * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
            var[bin.vars] <- means * (1 - means)
        }
    }
    
    var
}
col.w.cov <- function(mat, y, w = NULL, na.rm = TRUE) {
    if (!is.matrix(mat)) {
        if (is_null(w)) {
            return(cov(mat, y, use = if (na.rm) "pair" else "everything"))
        }
        mat <- matrix(mat, ncol = 1)
    }
    
    if (is_null(w)) {
        y <- array(y, dim = dim(mat))
        if (anyNA(mat)) is.na(y)[is.na(mat)] <- TRUE
        if (anyNA(y)) is.na(mat)[is.na(y)] <- TRUE
        den <- colSums(!is.na(mat*y)) - 1
        cov <- colSums(center(mat, na.rm = na.rm)*center(y, na.rm = na.rm), na.rm = na.rm)/den
    }
    else if (na.rm && anyNA(mat)) {
        n <- nrow(mat)
        w <- array(w, dim = dim(mat))
        is.na(w)[is.na(mat)] <- TRUE
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
    
    cov
}
col.w.r <- function(mat, y, w = NULL, s.weights = NULL, bin.vars = NULL, na.rm = TRUE) {
    if (is_null(w) && is_null(s.weights)) {
        return(cor(mat, y, w, use = if (na.rm) "pair" else "everything"))
    }
    
    cov <- col.w.cov(mat, y = y, w = w, na.rm = na.rm)
    den <- sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm)) *
        sqrt(col.w.v(y, w = s.weights, na.rm = na.rm))
    
    cov/den
}
.mean_abs_dev <- function(x) {
    mean_fast(abs(x - mean_fast(x, TRUE)), TRUE)
}
rms <- function(x) {
    sqrt(mean_fast(x^2))
}
.geam_mean <- function(y) {
    exp(mean_fast(log(y[is.finite(log(y))]), TRUE))
}
mat_div <- function(mat, vec) {
    mat/vec[col(mat)]
}
abs_ <- function(x, ratio = FALSE) {
    if (ratio) pmax(x, 1/x)
    else (abs(x))
}
mean_fast <- function(x, nas.possible = FALSE) {
    #Equal to mean(x, na.rm = TRUE) but faster
    #Set no.nas = FALSE if it's possible there are NAs
    if (nas.possible && anyNA(x)) {
        s <- sum(x, na.rm = TRUE)
        n <- sum(!is.na(x))
    }
    else {
        s <- sum(x)
        n <- length(x)
    }
    
    s/n
}
bw.nrd <- function(x) {
    #R's bw.nrd doesn't always work, but bw.nrd0 does
    bw.nrd0(x)*1.06/.9
}
ave_w.m <- function(x, ..., w = NULL) {
    #ave() version of w.m() since ave() doesn't do well with multiple variables to split
    if (missing(...)) {
        x[] <- w.m(x, w)
    }
    else {
        g <- interaction(...)
        split(x, g) <- lapply(levels(g), function(i) w.m(x[g == i], w[g == i]))
    }
    
    x
}

#Formulas
subbars <- function(term) {
    if (is.name(term) || !is.language(term))
        return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }
    
    if (is.call(term) && (term[[1]] == as.name("|") || term[[1]] == as.name("||"))) {
        term[[1]] <- as.name("+")
    }
    for (j in seq_len(term)[-1]) term[[j]] <- subbars(term[[j]])
    
    term
}

#treat/covs
assign.treat.type <- function(treat, use.multi = FALSE) {
    #Returns treat with treat.type attribute
    nunique.treat <- nunique(treat)
    
    if (nunique.treat < 2) {
        stop("The treatment must have at least two unique values.", call. = FALSE)
    }
    else if (!use.multi && nunique.treat == 2) {
        treat.type <- "binary"
    }
    else if (use.multi || chk::vld_character_or_factor(treat)) {
        treat.type <- "multinomial"
        if (!inherits(treat, "processed.treat")) treat <- factor(treat)
    }
    else {
        treat.type <- "continuous"
    }
    attr(treat, "treat.type") <- treat.type
    treat
}
get.treat.type <- function(treat) {
    attr(treat, "treat.type")
}
has.treat.type <- function(treat) {
    is_not_null(get.treat.type(treat))
}
get_treated_level <- function(treat) {
    if (is.character(treat) || is.factor(treat)) {
        treat <- factor(treat, nmax = 2)
        unique.vals <- levels(treat)
    }
    else {
        unique.vals <- unique(treat, nmax = 2)
    }
    
    if (is.character(unique.vals)) {
        if (can_str2num(unique.vals)) {
            zero <- {
                if (any(is.zero <- (str2num(unique.vals) == 0))) unique.vals[is.zero]
                else min(unique.vals, na.rm = TRUE)
            }
        }
        else {
            zero <- unique.vals[1]
        }
    }
    else {
        zero <- {
            if (any(as.numeric(unique.vals) == 0)) 0
            else min(unique.vals, na.rm = TRUE)
        }
    }
  
    setdiff(unique.vals, zero)
}

#Uniqueness
nunique <- function(x, nmax = NA, na.rm = TRUE) {
    if (is_null(x)) return(0)
    
    if (na.rm && anyNA(x)) x <- na.rem(x)
    # if (is.factor(x)) return(nlevels(x))
    # else 
    length(unique(x, nmax = nmax))
}
nunique.gt <- function(x, n, na.rm = TRUE) {
    if (missing(n)) stop("'n' must be supplied.")
    
    if (n < 0) stop("'n' must be non-negative.")
    
    if (is_null(x)) return(FALSE)
    
    if (n == 1) return(!all_the_same(x, na.rm))
    
    if (length(x) < 2000) return(nunique(x, na.rm = na.rm) > n)
    
    tryCatch(nunique(x, nmax = n, na.rm = na.rm) > n, error = function(e) TRUE)
}
all_the_same <- function(x, na.rm = TRUE) {
    if (anyNA(x)) {
        x <- na.rem(x)
        if (!na.rm) return(is_null(x))
    }
    if (is.double(x)) check_if_zero(max(x) - min(x))
    else all(x == x[1])
}
is_binary <- function(x, na.rm = TRUE) {
    if (na.rm && anyNA(x)) x <- na.rem(x)
    #Only return TRUE if not a continuous number that happens to take 2 values
    is_not_null(x) && 
        (!is.numeric(x) || all(check_if_zero(round(x) - x))) && 
        !all_the_same(x) && all_the_same(x[x != x[1]])
}
is_binary_col <- function(dat, na.rm = TRUE) {
    if (length(dim(dat)) != 2) stop("`is_binary_col()` cannot be used with objects that don't have 2 dimensions.")
    apply(dat, 2, is_binary)
}
is_0_1 <- function(x) {
    is_not_null(x) && 
        all(x == 1 | x == 0)
}

#R Processing
make_list <- function(n) {
    if (length(n) == 1L && is.numeric(n)) {
        vector("list", as.integer(n))
    }
    else if (is.atomic(n)) {
        setNames(vector("list", length(n)), as.character(n))
    }
    else stop("'n' must be an integer(ish) scalar or an atomic variable.")
}
make_df <- function(ncol, nrow = 0, types = "numeric") {
    if (length(ncol) == 1L && is.numeric(ncol)) {
        col_names <- NULL
        ncol <- as.integer(ncol)
    }
    else if (is.atomic(ncol)) {
        col_names <- as.character(ncol)
        ncol <- length(ncol)
    }
    if (length(nrow) == 1L && is.numeric(nrow)) {
        row_names <- NULL
        nrow <- as.integer(nrow)
    }
    else if (is.atomic(nrow)) {
        row_names <- as.character(nrow)
        nrow <- length(nrow)
    }
    df <- as.data.frame.matrix(matrix(NA_real_, nrow = nrow, ncol = ncol))
    colnames(df) <- col_names
    rownames(df) <- row_names
    if (is_not_null(types)) {
        if (length(types) %nin% c(1, ncol)) stop("'types' must be equal to the number of columns.")
        if (any(types %nin% c("numeric", "integer", "logical", "character", NA))) {
            stop("'types' must be an acceptable type. For factors, use NA.")
        }
        if (length(types) == 1) types <- rep(types, ncol)
        for (i in seq_len(ncol)) if (!is.na(types)[i] && types[i] != "numeric") {
            df[[i]] <- get(types[i])(nrow)
        }
    }
    
    df
}
ifelse_ <- function(...) {
    dotlen <- ...length()
    if (dotlen %% 2 == 0) stop("`ifelse_()` must have an odd number of arguments: pairs of test/yes, and one no.")
    out <- ...elt(dotlen)
    
    if (dotlen <= 1) {
        if (!is.atomic(out)) stop("The first entry to `ifelse_()` must be atomic.")
    }
    
    if (!is.atomic(out)) {
        stop("The last entry to `ifelse_()` must be atomic.")
    }
    
    if (length(out) == 1) out <- rep(out, length(..1))
    n <- length(out)
    for (i in seq_len((dotlen - 1)/2)) {
        test <- ...elt(2*i - 1)
        yes <- ...elt(2*i)
        if (length(yes) == 1) yes <- rep(yes, n)
        if (length(yes) != n || length(test) != n) stop("All entries must have the same length.")
        if (!is.logical(test)) stop(paste("The", ordinal(2*i - 1), "entry to `ifelse_()` must be logical."))
        if (!is.atomic(yes)) stop(paste("The", ordinal(2*i), "entry to `ifelse_()` must be atomic."))
        pos <- which(test)
        out[pos] <- yes[pos]
    }
    
    out
}
is_ <- function(x, types, stop = FALSE, arg.to = FALSE) {
    s1 <- deparse1(substitute(x))
    if (is_not_null(x)) {
        for (i in types) {
            if (i == "list") it.is <- is.list(x) && !is.data.frame(x)
            else if (is_not_null(get0(paste0("is_", i)))) {
                it.is <- get0(paste0("is_", i))(x)
            }
            else if (is_not_null(get0(paste.("is", i)))) {
                it.is <- get0(paste.("is", i))(x)
            }
            else it.is <- inherits(x, i)
            if (it.is) break
        }
    }
    else it.is <- FALSE
    
    if (stop && !it.is) {
        s0 <- ifelse(arg.to, "The argument to ", "")
        s2 <- ifelse(any(types %in% c("factor", "character", "numeric", "logical")),
                     "vector", "")
        stop(paste0(s0, "'", s1, "' must be a ", word_list(types, and.or = "or"), " ", s2, "."), call. = FALSE)
    }
    
    it.is
}
is_mat_like <- function(x) {
    length(dim(x)) == 2
}
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
if_null_then <- function(x1 = NULL, x2 = NULL, ...) {
    if (is_not_null(x1)) return(x1)
    if (is_not_null(x2)) return(x2)
    if (...length() == 0) return(x1)
    
    for (k in seq_len(...length())) {
        elt_k <- ...elt(k)
        if (is_not_null(elt_k)) return(elt_k)
    }
    
    ..1
}
clear_null <- function(x) {
    x[vapply(x, function(i) {
        is_null(i) &&
            (is_null(dim(i)) || all(dim(i) == 0L))
    }, logical(1L))] <- NULL
    x
}
clear_attr <- function(x, all = FALSE) {
    if (all) {
        attributes(x) <- NULL
    }
    else {
        dont_clear <- c("names", "class", "dim", "dimnames", "row.names")
        attributes(x)[names(attributes(x)) %nin% dont_clear] <- NULL
    }
    x
}
probably.a.bug <- function() {
    fun <- paste(deparse1(sys.call(-1)), collapse = "\n")
    .err(paste0("An error was produced and is likely a bug. Please let the maintainer know a bug was produced by the function\n",
                fun), tidy = FALSE)
}
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
`%pin%` <- function(x, table) {
    #Partial in. TRUE if x uniquely identifies values in table.
    !is.na(pmatch(x, table))
}
`%cin%` <- function(x, table) {
    #Partial in w/ charmatch. TRUE if x at all in table.
    !is.na(charmatch(x, table))
}
is_error <- function(x) {inherits(x, "try-error")}
null_or_error <- function(x) {is_null(x) || is_error(x)}
match_arg <- function(arg, choices, several.ok = FALSE) {
    #Replaces match.arg() but gives cleaner error message and processing
    #of arg.
    if (missing(arg))
        .err("no argument was supplied to `match_arg()`", call. = FALSE)
    arg.name <- deparse1(substitute(arg))
    
    if (missing(choices)) {
        formal.args <- formals(sys.function(sysP <- sys.parent()))
        choices <- eval(formal.args[[as.character(substitute(arg))]],
                        envir = sys.frame(sysP))
    }
    
    if (is.null(arg))
        return(choices[1L])
    if (!is.character(arg))
        .err(sprintf("the argument to `%s` must be `NULL` or a character vector", arg.name))
    if (!several.ok) {
        if (identical(arg, choices))
            return(arg[1L])
        if (length(arg) > 1L)
            .err(sprintf("the argument to `%s` must have length 1", arg.name))
    }
    else if (is_null(arg))
        .err(sprintf("the argument to `%s` must have length greater than or equal to 1", arg.name))
    
    i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
        .err(sprintf("the argument to `%s` should %s %s",
                     arg.name,
                     if (length(choices) > 1) {if (several.ok) "be at least one of" else "be one of"} else "be",
                     word_list(choices, and.or = "or", quotes = 2)))
    i <- i[i > 0L]
    if (!several.ok && length(i) > 1)
        .err("there is more than one match in `match_arg()`")
    choices[i]
}
grab <- function(x, what) {
    lapply(x, function(z) z[[what]])
}
last <- function(x) {
    x[[length(x)]]
}
`last<-` <- function(x, value) {
    `[[<-`(x, length(x), value)
}
len <- function(x, recursive = TRUE) {
    if (is_null(x)) 0L
    else if (length(dim(x)) > 1) NROW(x)
    else if (is.list(x) && recursive) vapply(x, len, numeric(1L), recursive = FALSE)
    else length(x)
}
na.rem <- function(x) {
    #A faster na.omit for vectors
    x[!is.na(x)]
}
anyNA_col <- function(x) {
    colSums(is.na(x)) > 0
}
check_if_call_from_fun <- function(fun) {
    # Check if called from within function f
    if (missing(fun) || !exists(deparse1(substitute(fun)), mode = "function")) return(FALSE)
    sp <- sys.parents()
    sys.funs <- lapply(sp, sys.function)
    for (x in sys.funs) {
        if (identical(fun, x)) return(TRUE)
    }
    FALSE
}
has_method <- function(class, fun) {
    if (!is.character(fun) || length(fun) != 1) stop("'fun' must be a string of length 1.")
    if (!is.character(class)) stop("'class' must be a character vector.")
    
    vapply(class, function(cl) is_not_null(getS3method(fun, cl, optional = TRUE, envir = asNamespace(packageName()))),
           logical(1L))
}
