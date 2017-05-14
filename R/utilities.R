#Utility functions
f.build <- function(y, rhs) {
    if ((is.data.frame(rhs) || is.matrix(rhs)) && length(colnames(rhs)) > 0)
        vars <- colnames(rhs)
    else if (is.character(rhs) && length(rhs) > 0)
        vars <- rhs
    else stop("Right hand side argument to f.build() must be a vector of variable names or a data set with named variables.", call. = FALSE)
    if (!(is.character(y) && length(y) == 1)) stop ("Response argument to f.build() must be the quoted name of the response variable.", call. = FALSE)
    if (y == "") y <- NULL
    f <- reformulate(vars, y)
    return(f)
}
inxnoty <- function(x, y) {
    #Creates a list or data frame of names in x that are not in y.
    #Useful for subsetting data sets into two groups of variables.
    
    if (!(is.character(x) || is.data.frame(x) || is.matrix(x)) || !(is.character(y) || is.data.frame(y) || is.matrix(y))) {
        stop("Inputs to x and y must be either strings containing variable names or data frames or matrices with named columns.", call. = FALSE)
    }
    if (is.character(x)) X <- x else X <- colnames(x)
    if (is.character(y)) Y <- y else Y <- colnames(y)
    
    if(is.character(x)) out <- x[is.na(match(X, Y))]
    else out <- x[, is.na(match(X, Y))]
    return(out)
}
word.list <- function(word.list = NULL, and.or = "and") {
    #When given a vector of strings, creates a string of the form "a and b"
    #or "a, b, and c"
    L <- length(word.list)
    if (L == 0) {
        out <- ""
    }
    else {
        word.list <- word.list[!word.list %in% c(NA, "")]
        L <- length(word.list)
        if (L == 0) {
            out <- ""
        }
        else if (L == 1) {
            out <- word.list
        }
        else if (L == 2) {
            out <- paste(word.list, collapse = paste0(" ", and.or," "))
        }
        else if (L >= 3) {
            out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "), 
                         word.list[L], sep = paste0(", ", and.or," "))
        }
            
    }
    return(out)
}
f.recode <- function(f, ...) {
    #Simplified version of forcats::fct_recode
    f <- factor(f)
    new_levels <- unlist(list(...))
    old_levels <- levels(f)
    idx <- match(new_levels, old_levels)
    
    old_levels[idx] <- names(new_levels)
    
    levels(f) <- old_levels
    return(f)
}