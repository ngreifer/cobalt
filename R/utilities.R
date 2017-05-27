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
word.list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE) {
    #When given a vector of strings, creates a string of the form "a and b"
    #or "a, b, and c"
    #If is.are, adds "is" or "are" appropriately
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
            if (is.are) out <- paste(out, "is")
        }
        else {
            and.or <- match.arg(and.or)
            if (L == 2) {
                out <- paste(word.list, collapse = paste0(" ", and.or," "))
                if (is.are) out <- paste(out, "are")
            }
            else {
                out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "), 
                             word.list[L], sep = paste0(", ", and.or," "))
                
            }
            if (is.are) out <- paste(out, "are")
        }

            
    }
    return(out)
}
split.factor <- function(data, varname, replace = TRUE, sep = "_", drop.level = NULL, drop.first = c("if2", TRUE, FALSE), drop.singleton = FALSE) {
    #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
    #Retains all categories unless only 2 levels, in which case only the second level is retained.
    #If variable only has one level, will delete.
    #varname= the name of the variable to split when data is specified
    #data=data set to be changed
    
    if (is.data.frame(data)) {
        factor.names <- names(data)[sapply(data, is.factor)]
        if (length(factor.names) == 0) {
            stop("There are no factor variables to split in data.", call. = FALSE)
        }
        if (missing(varname)) {
            varname <- factor.names
        }
        else if (is.character(varname)) {
            if (any(varname %in% factor.names)) {
                if (any(!varname %in% factor.names)) {
                    not.in.factor.names <- varname[!varname %in% factor.names]
                    warning(paste(word.list(not.in.factor.names, "and", is.are = TRUE), 
                                  "not the name(s) of factor variable(s) in data and will not be split."), 
                            call. = FALSE)
                }
                varname <- varname[varname %in% factor.names]
            }
            else {
                stop("No names in varname are names of factor variables in data.", call. = FALSE)
            }
        }
        else {
            stop("varname must be a character vector of the name(s) of factor variable(s) in data.", call. = FALSE)
        }
    }
    else if (is.factor(data)) {
        dep <- deparse(substitute(data))
        data <- data.frame(data)
        if (missing(varname)) {
            names(data) <- dep
        }
        else if (is.vector(varname) && (is.atomic(varname) || is.factor(varname))) {
            if (length(varname) == 0) {
                names(data) <- dep
            }
            else if (length(varname) == 1) {
                names(data) <- varname
            }
            else {
                warning("Only using the first item of varname.", call. = FALSE)
                names(data) <- varname[1]
            }
        }
        else {
            stop("varname must be an atomic or factor vector of length 1 with the stem of the new variable.", call. = FALSE)
        }
        varname <- names(data)
    }
    else {
        stop("data must a be a data.frame or a factor vector.", call. = FALSE)
    }
    
    if (length(drop.level) > 0 && length(varname) > 1) {
        warning("drop.level cannot be used with multiple entries to varname. Ignoring drop.level.", call. = FALSE)
        drop.level <- NULL
    }
    
    attributes.list <- setNames(vector("list", length(varname)), varname)
    for (v in varname) {
        drop <- character(0)
        x <- factor(data[, names(data) == v])
        
        skip <- FALSE
        if (nlevels(x) > 1) {
            k <- model.matrix(as.formula(paste0("~", v, "- 1")), data = data)
        }
        else {
            if (drop.singleton) {
                data <- data[, names(data)!=v, drop = FALSE]
                skip <- TRUE
            }
            else {
                k <- matrix(1, ncol = 1, nrow = length(x))
                colnames(k) <- paste0(v, levels(x)[1])
            }
        }
        
        if (!skip) {
            colnames(k) <- paste(v, sapply(strsplit(colnames(k), v, fixed = TRUE), function(n) paste(n, collapse = "")), sep = sep)
            
            if (length(drop.level) > 0) {
                if (is.character(drop.level) && length(drop.level) == 1 && drop.level %in% levels(x)) {
                    drop <- drop.level
                }
                else {
                    stop(paste("drop must be the name of a level of", v, "which is to be dropped."), call. = FALSE)
                }
            }
            else {
                if ((ncol(k) == 2 && drop.first %in% c("if2", TRUE)) ||
                    (ncol(k) > 2 && drop.first == TRUE)) {
                    drop <- levels(x)[1]
                }
            }
            if (length(drop) > 0) k <- k[, levels(x)!=drop, drop = FALSE]
            
            if (ncol(data) == 1) {
                data <- data.frame(k)
            }
            else if (replace) {
                if (match(v, names(data)) == 1){
                    data <- data.frame(k, data[, names(data)!=v, drop = FALSE])
                }
                else if (match(v, names(data)) == ncol(data)) {
                    data <- data.frame(data[, names(data)!=v, drop = FALSE], k)
                }
                else {
                    where <- match(v, names(data))
                    data <- data.frame(data[, 1:(where-1), drop = FALSE], k, data[, (where+1):ncol(data), drop = FALSE])
                }
            }
            else {
                data <- data.frame(data, k)
            }
            #Give attr so factor can be unsplit
            attributes.list[[v]] <- list(varname = varname,
                               sep = sep,
                               levels = levels(x),
                               drop = drop)

        }

    }
    

    attr(data, "split.factor") <- attributes.list

    return(data)
}

#Under construction
un.split.factor <- function(data, varname, replace = TRUE) {
    
    return(data)
}
get.weights <- function(obj, ...) UseMethod("get.weights")
get.weights.matchit <- function(m,...) {
    
}
get.weights.ps <- function(m,...) {
    
}
get.weights.Match <- function(m,...) {
    
}
get.weights.CBPS <- function(m,...) {
    
}
get.weights.ebalance <- function(m,...) {
    
}
get.weights.optmatch <- function(m,...) {
    
}
