#Utility functions
f.build <- function(y, rhs) {
    if (missing(rhs)) stop("Right hand side argument to f.build() must be a vector of variable names or a data set with named variables.", call. = FALSE)
    else if ((is.data.frame(rhs) || is.matrix(rhs)) && is_not_null(colnames(rhs))) {
        vars <- paste0("`", gsub("`", "", colnames(rhs)), "`")
    }
    else if (is.character(rhs) && is_not_null(rhs)) {
        vars <- paste0("`", gsub("`", "", rhs), "`")
    }
    else stop("Right hand side argument to f.build() must be a vector of variable names or a data set with named variables.", call. = FALSE)
    
    if (missing(y) || identical(y, "")) y <- NULL
    else if (!is.character(y) || length(y) > 1) stop ("Response argument to f.build() must be the quoted name of the response variable.", call. = FALSE)
    
    f <- reformulate(vars, y)
    return(f)
}
splitfactor <- function(data, var.name, replace = TRUE, sep = "_", drop.level = NULL, drop.first = TRUE, drop.singleton = FALSE, drop.na = TRUE, check = TRUE) {
    #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
    #Retains all categories unless only 2 levels, in which case only the second level is retained.
    #If variable only has one level, will delete.
    #var.name= the name of the variable to split when data is specified
    #data=data set to be changed
    
    if (is.data.frame(data)) {
        data <- as.data.frame(data)
        if (check) {
            factor.names <- names(data)[vapply(data, function(x) is.factor(x) || is.character(x), logical(1L))]
            if (missing(var.name)) {
                var.name <- factor.names
            }
            else if (is.character(var.name)) {
                if (any(var.name %in% factor.names)) {
                    if (any(!var.name %in% factor.names)) {
                        not.in.factor.names <- var.name[!var.name %in% factor.names]
                        warning(paste(word_list(not.in.factor.names, "and", is.are = TRUE), 
                                      "not the name(s) of factor variable(s) in data and will not be split."), 
                                call. = FALSE)
                    }
                    var.name <- var.name[var.name %in% factor.names]
                }
                else {
                    stop("No names in var.name are names of factor variables in data.", call. = FALSE)
                }
            }
            else {
                stop("var.name must be a character vector of the name(s) of factor variable(s) in data.", call. = FALSE)
            }
            if (is_null(factor.names)) {
                stop("There are no factor variables to split in data.", call. = FALSE)
            }
        }
        else {
            if (missing(var.name) || !is.character(var.name)) {
                stop("var.name must be a character vector of the names of variables in data.", call. = FALSE)
            }
            else {
                if (any(var.name %in% names(data))) {
                    if (any(var.name %nin% names(data))) {
                        not.in.data.names <- var.name[!var.name %in% names(data)]
                        warning(paste(word_list(not.in.data.names, "and", is.are = TRUE), 
                                      "not the name(s) of variable(s) in data and will not be split."), 
                                call. = FALSE)
                    }
                    var.name <- var.name[var.name %in% names(data)]
                }
                else {
                    stop("No names in var.name are names of variables in data.", call. = FALSE)
                }
            }
        }
        
    }
    else if (is.atomic(data)) {
        dep <- deparse(substitute(data))
        data <- data.frame(data)
        if (missing(var.name)) {
            names(data) <- dep
        }
        else if (is.vector(var.name) && (is.atomic(var.name) || is.factor(var.name))) {
            if (is_null(var.name)) {
                names(data) <- dep
            }
            else if (length(var.name) == 1) {
                names(data) <- var.name
            }
            else {
                warning("Only using the first item of var.name.", call. = FALSE)
                names(data) <- var.name[1]
            }
        }
        else {
            stop("var.name must be an atomic or factor vector of length 1 with the stem of the new variable.", call. = FALSE)
        }
        var.name <- names(data)
    }
    else {
        stop("data must a be a data.frame or an atomic vector.", call. = FALSE)
    }
    
    if (is_not_null(drop.level) && length(var.name) > 1) {
        warning("drop.level cannot be used with multiple entries to var.name. Ignoring drop.level.", call. = FALSE)
        drop.level <- NULL
    }
    drop.na <- setNames(rep(drop.na, length(var.name)), var.name)
    for (v in var.name) {
        drop <- character(0)
        x <- factor(data[names(data) == v][[1]], exclude = NULL)
        na.level <- is.na(levels(x))
        levels(x) <- paste0(sep, levels(x))
        data[names(data) == v][[1]] <- x
        
        skip <- FALSE
        if (nlevels(x) > 1) {
            k <- model.matrix(as.formula(paste0("~`", v, "`- 1")), data = data)
            colnames(k) <- gsub("`", colnames(k), replacement = "")
            
            if (any(na.level)) {
                if (drop.na[v]) {
                    k[k[,na.level] == 1,] <- NA_real_
                }
            }
            else drop.na[v] <- FALSE
            
        }
        else {
            if (drop.singleton) {
                data <- data[names(data)!=v]
                skip <- TRUE
            }
            else {
                k <- matrix(1, ncol = 1, nrow = length(x))
                colnames(k) <- paste0(v, levels(x)[1])
            }
        }
        
        if (!skip) {
            if (is_not_null(drop.level)) {
                if (is.character(drop.level) && length(drop.level) == 1 && drop.level %in% levels(x)) {
                    drop <- drop.level
                }
                else {
                    stop(paste("drop must be the name of a level of", v, "which is to be dropped."), call. = FALSE)
                }
            }
            else {
                if ((ncol(k) == 2 && (drop.first == "if2" || drop.first == TRUE)) ||
                    (ncol(k) > 2 && drop.first == TRUE)) {
                    drop <- levels(x)[1]
                }
            }
            
            dropl <- rep(FALSE, ncol(k))
            if (is_not_null(drop)) {
                dropl[!na.level & levels(x) %in% drop] <- TRUE
            }
            if (drop.na[v]) dropl[na.level] <- TRUE
            
            k <- k[,!dropl, drop = FALSE]
            
            if (ncol(data) == 1) {
                data <- data.frame(k, row.names = rownames(data))
            }
            else if (replace) {
                if (match(v, names(data)) == 1){
                    data <- cbind(k, data[names(data)!=v], row.names = rownames(data))
                }
                else if (match(v, names(data)) == ncol(data)) {
                    data <- cbind(data[names(data)!=v], k, row.names = rownames(data))
                }
                else {
                    where <- match(v, names(data))
                    data <- cbind(data[1:(where-1)], k, data[(where+1):ncol(data)], row.names = rownames(data))
                }
            }
            else {
                data <- data.frame(data, k, row.names = rownames(data))
            }
            
        }
        
    }
    
    return(data)
}
unsplitfactor <- function(data, var.name, replace = TRUE, sep = "_", dropped.level = NULL, dropped.na = TRUE) {
    
    if (!is.data.frame(data)) stop("data must be a data.frame containing the variables to unsplit.", call = FALSE)
    if (!is.character(var.name)) stop("var.name must be a string containing the name of the variables to unsplit.", call. = FALSE)
    if (is_not_null(dropped.level) && length(var.name) > 1) {
        warning("dropped.level cannot be used with multiple var.names and will be ignored.", call. = FALSE, immediate. = TRUE)
        dropped.level <- NULL
    }
    
    if (!is.character(var.name)) stop("var.name must be a character vector containing the name of the variable to unsplit.", call. = FALSE)
    if (length(sep) > 1 || !is.character(sep)) stop("sep must be a character vector of length 1 containing the seperating character in the names of the split variables.", call. = FALSE)
    if (length(dropped.level) > 1 && !is.atomic(dropped.level)) {
        warning("dopped.level must be an atomic vector of length 1 containing the value of the dropped category of the split variable. It will be ignored.", call. = FALSE, immediate. = TRUE)
        dropped.level <- NULL
    }
    not.the.stem <- character(0)
    
    for (v in var.name) {
        dropped.level0 <- dropped.level
        var.to.combine <- data[startsWith(names(data), paste0(v, sep))]
        if (is_null(var.to.combine)) {
            not.the.stem <- c(not.the.stem, paste0(v, sep))
            next
        }
        
        if (any(rowSums(apply(var.to.combine, 2, is.na)) %nin% c(0, ncol(var.to.combine)))) {
            stop("The variables in data selected based on var.name and sep do not seem to form a split variable based on the <NA> pattern.", call. = FALSE)
        }
        NA.column <- character(0)
        
        if (!isTRUE(dropped.na)) {
            NA.column <- paste0(v, sep, ifelse(dropped.na == FALSE, "NA", dropped.na))
            if (NA.column %in% names(var.to.combine)) {
                var.to.combine[var.to.combine[[NA.column]] == 1,] <- NA_real_
                var.to.combine <- var.to.combine[names(var.to.combine) != NA.column]
            }
            else {
                stop(paste("There is no variable called", word_list(NA.column, quotes = TRUE), "to generate the NA values."), call. = FALSE)
            }
        }
        var.sum <- rowSums(var.to.combine)
        if (isTRUE(all.equal(unique(var.sum), 1))) {
            #Already unsplit
        }
        else if (isTRUE(all.equal(sort(unique(var.sum)), c(0, 1)))) {
            #Missing category
            
            if (is_null(dropped.level)) {
                k.levels0 <- sapply(names(var.to.combine), function(x) strsplit(x, paste0(v, sep), fixed = TRUE)[[1]][2])
                
                if (suppressWarnings(all(!is.na(as.numeric(k.levels0))))) {
                    dropped.level0 <- as.character(min(as.numeric(k.levels0)) - 1)
                    dropped.name <- paste0(v, sep, dropped.level0)
                }
                else {
                    message("The dropped category will be set to NA.")
                    dropped.name <- dropped.level0 <- NA_character_
                }
                
            }
            else dropped.name <- paste0(v, sep, dropped.level)
            var.to.combine <- setNames(data.frame(1-var.sum, var.to.combine),
                                       c(dropped.name, names(var.to.combine)))
            
        }
        else {
            stop("The variables in data selected based on var.name and sep do not seem to form a split variable based on the row sums.", call. = FALSE)
        }
        
        k.levels <- vapply(names(var.to.combine), function(x) strsplit(x, paste0(v, sep), fixed = TRUE)[[1]][2], character(1L))
        
        k <- rep(NA_character_, nrow(data))
        for (i in seq_along(k.levels)) {
            k <- ifelse(var.to.combine[[i]] == 1, k.levels[i], k)
        }
        
        k <- factor(k, levels = k.levels)
        
        
        if (replace) {
            where <- which(names(data) %in% c(names(var.to.combine), NA.column))
            
            data[[min(where)]] <- k
            remove.cols <- where[where!=min(where)]
            if (is_not_null(remove.cols)) data <- data[-remove.cols]
            names(data)[min(where)] <- v
        }
        else {
            data <- cbind(data, setNames(data.frame(k), v))
        }
    }
    
    if (is_not_null(not.the.stem)) warning(paste0(word_list(not.the.stem, is.are = TRUE, quotes = TRUE), " not the stem of any variables in data and will be ignored. Ensure var.name and sep are correct."), call. = FALSE)
    
    return(data)
}

var.names <- function(b, type, file = NULL, minimal = FALSE) {
    if (is_not_null(b[["print.options"]][["co.names"]])) {
        if (minimal) vars <- unique(unlist(lapply(b[["print.options"]][["co.names"]], function(x) x[["component"]][x[["type"]] == "base"])))
        else vars <- vapply(b[["print.options"]][["co.names"]], function(x) paste(x[["component"]], collapse = ""), character(1))
    }
    else {
        vars <- NULL
        var.containers <- c(quote(b[["Balance"]]),
                            quote(b[["Cluster.Balance"]][[1]][["Balance"]]),
                            quote(b[["Subclass.Balance"]][[1]]),
                            quote(b[["Imputation.Balance"]][[1]][["Balance"]]),
                            quote(b[["Imputation.Balance"]][[1]][["Cluster.Balance"]][[1]][["Balance"]]),
                            quote(b[["Pair.Balance"]][[1]]),
                            quote(b[["Time.Balance"]][[1]][["Balance"]]))
        for (i in var.containers) {
            obj <- eval(i)
            if (is_not_null(obj)) {
                vars <- rownames(obj)
                break
            }
            else obj <- NULL
        }
        if (is_null(vars)) stop("No variable names were found in the object. It is probably not a bal.tab object.", call. = FALSE)
        if (minimal) warning("minimal is being set to FALSE because the part of the object required for it to be TRUE is missing.", call. = FALSE)
    }
    
    if (is_not_null(file)) {
        if (!endsWith(file, ".csv")) stop("The filename in file must end in \".csv\".", call. = FALSE)
    }
    
    if (missing(type)) {
        if (is_not_null(file)) type <- "df"
        else type <- "vec"
    }
    else {
        possible.types <- c("df", "vec")
        type <- possible.types[pmatch(type, possible.types)]
    }
    
    if (is.na(type)) stop("type must be \"df\" or \"vec\"")
    else if (type == "df") {
        out <- data.frame(old = vars, new = vars, stringsAsFactors = FALSE, row.names = NULL)
    }
    else {
        out <- setNames(vars, vars)
    }
    
    if (is_not_null(file)) {
        if (type == "df") {
            write.csv(out, file = file, row.names = FALSE)
            invisible(out)
        }
        else {
            warning("Only type = \"df\" is compatible with a file name.", call. = FALSE)
            out
        }
    }
    else out
}

set.cobalt.options <- function(..., default = FALSE) {
    opts <- list(...)
    if (is_not_null(opts) && (is_null(names(opts)) ||  "" %in% names(opts))) {
        stop("All arguments must be named.", call. = FALSE)
    }
    # if ("continuous" %in% names(opts)) names(opts)[names(opts) == "continuous"] <- "cont"
    # if ("binary" %in% names(opts)) names(opts)[names(opts) == "binary"] <- "bin"
    
    multiple.allowed <- c("cluster.fun", "imp.fun")
    any.string.allowed <- c("int_sep", "factor_sep")
    
    if (any(duplicates <- table(names(opts)) > 1)) {
        stop(paste0(word_list(names(duplicates)[duplicates], is.are = TRUE), " present more than once in the input to set.cobalt.options."), call. = FALSE)
    }
    
    if (any(names(opts) %nin% names(acceptable.options()))) {
        warning(paste("The following are not acceptable options and will be ignored:", word_list(unique(names(opts)[names(opts) %nin% names(acceptable.options())]))), call. = FALSE, immediate. = TRUE)
        opts <- opts[names(opts) %in% names(acceptable.options())]
    }
    
    if (default) {
        return.to.default <- setdiff(names(acceptable.options()), names(opts))
    }
    else return.to.default <- NULL

    multiple.opts <- NULL
    bad.opts <- NULL
    for (i in names(opts)) {
        if (is_null(opts[[i]])) {
            return.to.default <- c(return.to.default, i)
            opts[[i]] <- NULL
        }
        else {
            if (length(opts[[i]]) > 1 && i %nin% multiple.allowed) multiple.opts <- c(multiple.opts, i)
            if (mode(opts[[i]]) != mode(acceptable.options()[[i]]) || 
                (!(is.character(opts[[i]]) && is.character(acceptable.options()[[i]]) && (i %in% any.string.allowed || !is.na(pmatch(opts[[i]], acceptable.options()[[i]])))) &&
                 !all(opts[[i]] %in% acceptable.options()[[i]]))) bad.opts <- c(bad.opts, i)
        }
    }
    
    if (is_not_null(opts)) {
        both.opts <- intersect(multiple.opts, bad.opts)
        multiple.opts <- multiple.opts[multiple.opts %nin% both.opts]
        bad.opts <- bad.opts[bad.opts %nin% both.opts]
        problematic.opts <- setNames(vector("list", 3), c("multiple", "bad", "both"))
        problematic.opts[["multiple"]] <- setNames(lapply(multiple.opts, function(i) {
            paste(i, "must be of length 1.")
        }), multiple.opts)
        problematic.opts[["bad"]] <- setNames(lapply(bad.opts, function(i) {
            if (i %in% any.string.allowed) paste0(i, " must be a character string.")
            else paste0(i, " must be ", word_list(acceptable.options()[[i]], quotes = is.character(acceptable.options()[[i]]), and.or = "or"), ".")
        }), bad.opts)
        problematic.opts[["both"]] <- setNames(lapply(both.opts, function(i) {
            if (i %in% any.string.allowed) paste0(i, " must be a character string of length 1.")
            else paste0(i, " must be one of ", word_list(acceptable.options()[[i]], quotes = is.character(acceptable.options()[[i]]), and.or = "or"), ".")
        }), both.opts)
        
        problems <- do.call("c", unname(problematic.opts))
        problems <- problems[names(opts)[names(opts) %in% names(problems)]]
        if (is_not_null(problems)) {
            stop(do.call("paste", c(list(""), problems, list("\nNo options will be set.", sep = "\n"))), call. = FALSE)
        }
        
        names(opts) <- paste0("cobalt_", names(opts))
        options(opts)
    }
    
    if (is_not_null(return.to.default)) {
        options(setNames(replicate(length(return.to.default), NULL), paste0("cobalt_", return.to.default)))
    }
    # if ("continuous" %in% names(opts)) names(acceptable.options)[names(acceptable.options) == "continuous"] <- "cont"
    # if ("binary" %in% names(opts)) names(acceptable.options)[names(acceptable.options) == "binary"] <- "bin"
}
get.cobalt.options <- function(...) {
    opts <- list(...)
    
    opts <- clear_null(opts)
    if (is_null(opts)) opts <- names(acceptable.options())
    else {
        if (!all(vapply(opts, is.character, logical(1L)))) {
            stop("All arguments must be strings containing the name of an option to return.", call. = FALSE)
        }
        opts <- do.call("c", opts)
        if (any(not.in.accept <- opts %nin% names(acceptable.options()))) {
            plural <- sum(not.in.accept) > 1
            stop(paste0(word_list(opts[not.in.accept], is.are = TRUE, quotes = TRUE),
                        " not", ifelse(plural, "", " an"), " acceptable option", 
                        ifelse(plural, "s", ""), "."), call. = FALSE)
        }
    }
    
    out <- setNames(lapply(paste0("cobalt_", opts), getOption), opts)
    return(out)
    
}