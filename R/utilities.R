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
split.factor0 <- function(varname, data) {
    #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
    #Retains all categories.
    #varname= the name of the variable to split when data is specified
    #data=data set to be changed
    
    x <- data[, varname] <- factor(data[, varname])
    
    if (nlevels(x) > 1) {
        k <- model.matrix(as.formula(paste0("~", varname, "- 1")), data = data)
    }
    else {
        k <- matrix(1, ncol = nlevels(x), nrow = length(x))
    }
    
    if (is.null(levels(x))) colnames(k) <- paste0(varname, 1:ncol(k))
    else colnames(k) <- paste(varname, substr(sapply(strsplit(colnames(k), varname), function(n) paste(n, collapse = "")), 1, 10), sep = "_")
    
    if (match(varname, names(data)) == 1){
        data <- data.frame(k, data[, names(data)!=varname, drop = FALSE])
    }
    else if (match(varname, names(data)) == ncol(data)) {
        data <- data.frame(data[, names(data)!=varname, drop = FALSE], k)
    }
    else {
        where <- match(varname, names(data))
        data <- data.frame(data[, 1:(where-1), drop = FALSE], k, data[, (where+1):ncol(data), drop = FALSE])
    }
    return(data)
}
split.factor <- function(data, varname, replace = TRUE, sep = "_", drop.first = c("if2", TRUE, FALSE), drop.singleton = FALSE) {
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
                                  "not names of factor variables in data and will not be split."), 
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
        data <- data.frame(data)
        if (missing(varname)) {
            names(data) <- ""
        }
        else if (is.vector(varname) && (is.atomic(varname) || is.factor(varname))) {
            if (length(varname) == 0) {
                names(data) <- ""
                sep <- ""
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
            stop("varname must be an atmomic or factor vector of length 1 with the stem of the new variable.", call. = FALSE)
        }
        varname <- names(data)
    }
    else {
        stop("data must a be a data.frame or a factor vector.", call. = FALSE)
    }

    for (v in varname) {
        x <- data[, v] <- factor(data[, v])
        
        skip <- FALSE
        if (nlevels(x) > 1) {
            k <- model.matrix(as.formula(paste0("~", varname, "- 1")), data = data)
        }
        else {
            if (drop.singleton) {
                data <- data[, names(data)!=varname]
                skip <- TRUE
            }
            else {
                k <- matrix(1, ncol = 1, nrow = length(x))
                #colnames(k) <- paste0(varname, levels(x)[1])
            }
        }
        
        if (!skip) {
            
            if (nlevels(x) == 0) colnames(k) <- paste0(varname, 1:ncol(k))
            else colnames(k) <- paste(varname, sapply(strsplit(colnames(k), varname), function(n) paste(n, collapse = "")), sep = sep)
            
            if(ncol(k) == 2){
                if (drop.first %in% c("if2", TRUE)) {
                    k <- k[,-1, drop = FALSE]
                }
            } 
            else if (ncol(k) > 2) {
                if (drop.first == TRUE) {
                    k <- k[,-1, drop = FALSE]
                }
            }
            
            if (ncol(data) == 1) {
                data <- data.frame(k)
            }
            else if (replace) {
                if (match(varname, names(data)) == 1){
                    data <- data.frame(k, data[, names(data)!=varname, drop = FALSE])
                }
                else if (match(varname, names(data)) == ncol(data)) {
                    data <- data.frame(data[, names(data)!=varname, drop = FALSE], k)
                }
                else {
                    where <- match(varname, names(data))
                    data <- data.frame(data[, 1:(where-1), drop = FALSE], k, data[, (where+1):ncol(data), drop = FALSE])
                }
            }
            else {
                data <- data.frame(data, k)
            }

        }

    }

    return(data)
}