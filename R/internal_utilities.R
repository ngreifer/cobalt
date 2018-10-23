#Internal Utilities
word.list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
    #When given a vector of strings, creates a string of the form "a and b"
    #or "a, b, and c"
    #If is.are, adds "is" or "are" appropriately
    L <- length(word.list)
    if (quotes) word.list <- vapply(word.list, function(x) paste0("\"", x, "\""), character(1L))
    if (L == 0) {
        out <- ""
        attr(out, "plural") = FALSE
    }
    else {
        word.list <- word.list[!word.list %in% c(NA_character_, "")]
        L <- length(word.list)
        if (L == 0) {
            out <- ""
            attr(out, "plural") = FALSE
        }
        else if (L == 1) {
            out <- word.list
            if (is.are) out <- paste(out, "is")
            attr(out, "plural") = FALSE
        }
        else {
            and.or <- match.arg(and.or)
            if (L == 2) {
                out <- paste(word.list, collapse = paste0(" ", and.or," "))
            }
            else {
                out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "),
                             word.list[L], sep = paste0(", ", and.or," "))
                
            }
            if (is.are) out <- paste(out, "are")
            attr(out, "plural") = TRUE
        }
        
        
    }
    return(out)
}
expand.grid_string <- function(..., collapse = "") {
    return(apply(expand.grid(...), 1, paste, collapse = collapse))
}
nunique <- function(x, nmax = NA_real_, na.rm = TRUE) {
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
        if (n == 1 && is.numeric(x)) !check_if_zero(max(x) - min(x))
        else if (length(x) < 2000) nunique(x) > n
        else tryCatch(nunique(x, nmax = n) > n, error = function(e) TRUE)
    }
}
is_binary <- function(x) !nunique.gt(x, 2)
all_the_same <- function(x) !nunique.gt(x, 1)
is.formula <- function(f, sides = NULL) {
    res <- is.name(f[[1]])  && deparse(f[[1]]) %in% c( '~', '!') &&
        length(f) >= 2
    if (is_not_null(sides) && is.numeric(sides) && sides %in% c(1,2)) {
        res <- res && length(f) == sides + 1
    }
    return(res)
}
is.char.or.factor <- function(x) {
    return(is_not_null(x) && (is.factor(x) || is.character(x)))
}
check_if_zero <- function(x) {
    # this is the default tolerance used in all.equal
    tolerance <- .Machine$double.eps^0.5
    # If the absolute deviation between the number and zero is less than
    # the tolerance of the floating point arithmetic, then return TRUE.
    # This means, to me, that I can treat the number as 0 rather than
    # -3.20469e-16 or some such.
    abs(x - 0) < tolerance
}
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
`%+%` <- function(...) {
    a <- list(...)
    if (is.atomic(a[[1]])) do.call(crayon::`%+%`, a)
    else do.call(ggplot2::`%+%`, a)
}
strsplits <- function(x, splits, fixed = TRUE, ...) {
    #Link strsplit but takes multiple split values.
    #Only works for one string at a time (in x).
    for (split in splits) x <- unlist(strsplit(x, split, fixed = TRUE, ...))
    return(x[x != ""]) # Remove empty values
}
paste. <- function(..., collapse = NULL) {
    #Like paste0 but with sep = ".'
    paste(..., sep = ".", collapse = collapse)
}
is_ <- function(x, types, stop = FALSE, arg.to = FALSE) {
    s1 <- deparse(substitute(x))
    if (is_not_null(x)) {
        for (i in types) {
            if (is_not_null(get0(paste.("is", i)))) {
                it.is <- get0(paste.("is", i))(x)
            }
            else it.is <- inherits(x, i)
            if (it.is) break
        }
    }
    else it.is <- FALSE
    
    if (stop) {
        if (!it.is) {
            s0 <- ifelse(arg.to, "The argument to ", "")
            s2 <- ifelse(any(types %in% c("factor", "character", "numeric", "logical")),
                         "vector", "")
            stop(paste0(s0, s1, " must be a ", word.list(types, and.or = "or"), " ", s2, "."), call. = FALSE)
        }
    }
    else {
        return(it.is)
    }
}
