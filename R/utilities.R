#Utility functions
f.build <- function(y, rhs) {
    if (missing(rhs)) {
        if (missing(y)) stop("The right hand side argument to f.build() must be a vector of variable names or a data set with named variables.", call. = FALSE)
        else {
            rhs <- y
            y <- ""
        }
    }
    
    tryCatch(force(y), error = function(e) stop(conditionMessage(e), call. = FALSE))
    tryCatch(force(rhs), error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    if (is_(rhs, c("matrix", "data.frame")) && is_not_null(colnames(rhs))) {
        vars <- paste0("`", gsub("`", "", colnames(rhs)), "`")
    }
    else if (is_(rhs, "character")) {
        vars <- paste0("`", gsub("`", "", rhs), "`")
    }
    else stop("Right hand side argument to f.build() must be a vector of variable names or a data set with named variables.", call. = FALSE)
    
    if (missing(y) || is_null(y) || identical(y, "")) y <- NULL
    else if (!is_(y, "atomic")) stop ("Response argument to f.build() must be a string containing the response variable.", call. = FALSE)
    
    f <- formula(paste(
        paste(as.character(y), collapse = " + ") , "~", paste(vars, collapse = " + ")
    ))
    # f <- reformulate(vars, y)
    return(f)
}
splitfactor <- function(data, var.name, drop.level = NULL, drop.first = TRUE, drop.singleton = FALSE, drop.na = TRUE, sep = "_", replace = TRUE, split.with = NULL, check = TRUE) {
    #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
    #Retains all categories unless only 2 levels, in which case only the second level is retained.
    #If variable only has one level, will delete.
    #var.name= the name of the variable to split when data is specified
    #data=data set to be changed
    
    if (is.data.frame(data)) {
        data <- as.data.frame(data)
        if (check) {
            factor.names <- names(data)[vapply(data, is_, logical(1L), c("factor", "character"))]
            if (missing(var.name)) {
                var.name <- factor.names
            }
            else if (is.character(var.name)) {
                if (any(var.name %in% factor.names)) {
                    if (any(var.name %nin% factor.names)) {
                        not.in.factor.names <- var.name[var.name %nin% factor.names]
                        warning(paste(word_list(not.in.factor.names, "and", is.are = TRUE), 
                                      "not the name(s) of factor variable(s) in 'data' and will not be split."), 
                                call. = FALSE)
                    }
                    var.name <- var.name[var.name %in% factor.names]
                }
                else {
                    stop("No names in 'var.name' are names of factor variables in 'data'.", call. = FALSE)
                }
            }
            else {
                stop("'var.name' must be a character vector of the name(s) of factor variable(s) in 'data'.", call. = FALSE)
            }
            if (is_null(factor.names)) {
                warning("There are no factor variables to split in 'data'.", call. = FALSE)
                return(data)
            }
        }
        else {
            if (missing(var.name) || !is.character(var.name)) {
                stop("'var.name' must be a character vector of the names of variables in 'data'.", call. = FALSE)
            }
            else {
                if (any(var.name %in% names(data))) {
                    if (any(var.name %nin% names(data))) {
                        not.in.data.names <- var.name[!var.name %in% names(data)]
                        warning(paste(word_list(not.in.data.names, "and", is.are = TRUE), 
                                      "not the name(s) of variable(s) in 'data' and will not be split."), 
                                call. = FALSE)
                    }
                    var.name <- var.name[var.name %in% names(data)]
                }
                else {
                    stop("No names in 'var.name' are names of variables in 'data'.", call. = FALSE)
                }
            }
        }
        
        if (is_not_null(split.with)) {
            if (is_(split.with, "list")) {
                if (any(vapply(split.with, function(x) !is_(x, "atomic"), logical(1L)))) 
                    stop("All entries in 'split.with' must must be atomic vectors or factors.", call. = FALSE)
                if (any(vapply(split.with, function(x) length(x) != ncol(data), logical(1L)))) 
                    stop("All entries in 'split.with' must have length equal to the number of columns of 'data'.", call. = FALSE)
            }
            else {
                if (!is_(split.with, "atomic"))
                    stop("'split.with' must must be an atomic vector or factor or list thereof.", call. = FALSE)
                if (length(split.with) != ncol(data))
                    stop("'split.with' must have length equal to the number of columns of 'data'.", call. = FALSE)
                split.with <- list(split.with)
            }
        }
    }
    else if (is.atomic(data)) {
        dep <- deparse1(substitute(data))
        data <- data.frame(data)
        if (missing(var.name) || is_null(var.name)) {
            names(data) <- dep
        }
        else if (is_(var.name, "atomic")) {
            if (length(var.name) == 1) {
                names(data) <- var.name
            }
            else {
                warning("Only using the first item of 'var.name'.", call. = FALSE)
                names(data) <- var.name[1]
            }
        }
        else {
            stop("'var.name' must be an atomic or factor vector of length 1 with the stem of the new variable.", call. = FALSE)
        }
        var.name <- names(data)
        
        if (is_not_null(split.with)) {
            if (is_(split.with, "list")) {
                if (any(vapply(split.with, function(x) !is_(x, "atomic"), logical(1L)))) 
                    stop("All entries in 'split.with' must must be atomic vectors or factors.", call. = FALSE)
                if (any(vapply(split.with, function(x) length(x) != ncol(data), logical(1L)))) 
                    stop("All entries in 'split.with' must have length 1.", call. = FALSE)
            }
            else {
                if (!is_(split.with, "atomic"))
                    stop("'split.with' must must be an atomic vector or factor or list thereof.", call. = FALSE)
                if (length(split.with) != ncol(data))
                    stop("'split.with' must have length 1.", call. = FALSE)
                split.with <- list(split.with)
            }
        }
    }
    else {
        stop("'data' must a be a data.frame or factor.", call. = FALSE)
    }
    
    if (is_not_null(drop.level) && length(var.name) > 1) {
        warning("'drop.level' cannot be used with multiple entries to 'var.name'. Ignoring 'drop.level'.", call. = FALSE)
        drop.level <- NULL
    }
    
    if (!rlang::is_bool(drop.na)) stop("'drop.singleton' must be TRUE or FALSE.", call. = FALSE)
    drop.na <- rlang::rep_named(var.name, drop.na)
    
    for (v in var.name) {
        x <- factor(data[v][[1]])
        na.level <- if (anyNA(x)) nlevels(x) + 1 else 0
        new.levels <- c(levels(x), if (na.level) NA_character_)
        
        if (length(new.levels) > 1) {
            k <- make_df(paste(v, new.levels, sep = sep), nrow(data), types = "integer")
            for (i in seq_len(ncol(k))) {
                if (i != na.level) k[!is.na(x) & x == new.levels[i], i] <- 1L
                else k[is.na(x), i] <- 1L
            }
            
            if (na.level > 0) {
                if (drop.na[v]) k[is.na(x), -na.level] <-  NA_integer_
            }
            else drop.na[v] <- FALSE
            
        }
        else {
            if (!rlang::is_bool(drop.singleton)) stop("'drop.singleton' must be TRUE or FALSE.", call. = FALSE)
            if (drop.singleton) {
                data[[v]] <- NULL
                next
            }
            else {
                k <- make_df(paste(v, new.levels[1], sep = sep), length(x), types = "integer")
            }
        }
        
        dropl <- rlang::rep_named(new.levels, FALSE)
        if (is_not_null(drop.level)) {
            if (is.character(drop.level) && length(drop.level) == 1 && drop.level %in% new.levels) {
                dropl[drop.level] <- TRUE
            }
            else {
                stop(paste("'drop' must be the name of a level of", v, "that is to be dropped."), call. = FALSE)
            }
        }
        else {
            if (!identical(drop.first, "if2") && !rlang::is_bool(drop.first)) stop("'drop.first' must be TRUE, FALSE, or \"if2\".", call. = FALSE)
            if ((ncol(k) == 2 && (drop.first == "if2" || isTRUE(drop.first))) ||
                (ncol(k) > 2 && isTRUE(drop.first))) {
                dropl[1] <- TRUE
            }
        }
        
        if (drop.na[v]) dropl[na.level] <- TRUE
        
        for (i in seq_len(ncol(k))) {
            attr(k[[i]], "split.var") <- v
            attr(k[[i]], "level") <- {if (i == na.level) NA_character_ else new.levels[i]}
        }
        
        k[dropl] <- NULL
        
        if (is_not_null(split.with)) {
            for (i in seq_along(split.with)) names(split.with[[i]]) <- names(data)
        }
        
        if (ncol(data) == 1) {
            data <- k
            if (is_not_null(split.with)) {
                for (i in seq_along(split.with)) {
                    split.with[[i]] <- rep(split.with[[i]][v], ncol(k))
                    names(split.with[[i]]) <- names(data)
                }
            }
        }
        else {
            if (!rlang::is_bool(replace)) stop("'replace' must be TRUE or FALSE.", call. = FALSE)
            if (replace) {
                if (match(v, names(data)) == 1){
                    data <- setNames(data.frame(k, data[names(data)!=v], row.names = rownames(data)),
                                     c(names(k), names(data)[names(data)!=v]))
                    if (is_not_null(split.with)) {
                        for (i in seq_along(split.with)) {
                            split.with[[i]] <- c(rep(split.with[[i]][v], ncol(k)), split.with[[i]][names(split.with[[i]]) != v])
                            names(split.with[[i]]) <- names(data)
                        }
                    }
                }
                else if (match(v, names(data)) == ncol(data)) {
                    data <- setNames(data.frame(data[names(data)!=v], k, row.names = rownames(data)),
                                     c(names(data)[names(data)!=v], names(k)))
                    if (is_not_null(split.with)) {
                        for (i in seq_along(split.with)) {
                            split.with[[i]] <- c(split.with[[i]][names(split.with[[i]]) != v], rep(split.with[[i]][v], ncol(k)))
                            names(split.with[[i]]) <- names(data)
                        }
                    }
                }
                else {
                    where <- match(v, names(data))
                    data <- setNames(data.frame(data[1:(where-1)], k, data[(where+1):ncol(data)], row.names = rownames(data)),
                                     c(names(data)[1:(where-1)], names(k), names(data)[(where+1):ncol(data)]))
                    if (is_not_null(split.with)) {
                        for (i in seq_along(split.with)) {
                            split.with[[i]] <- c(split.with[[i]][1:(where-1)], rep(split.with[[i]][v], ncol(k)), split.with[[i]][(where+1):length(split.with[[i]])])
                            names(split.with[[i]]) <- names(data)
                        }
                    }
                }
            }
            else {
                data <- setNames(data.frame(data, k, row.names = rownames(data)),
                                 c(names(data), names(k)))
                if (is_not_null(split.with)) {
                    for (i in seq_along(split.with)) {
                        split.with[[i]] <- c(split.with[[i]], rep(split.with[[i]][v], ncol(k)))
                        names(split.with[[i]]) <- names(data)
                    }
                }
            }
        }
        
    }
    
    attr(data, "split.with") <- split.with
    return(data)
}
unsplitfactor <- function(data, var.name, dropped.level = NULL, dropped.na = TRUE, sep = "_", replace = TRUE) {
    
    if (!is.data.frame(data)) stop("'data' must be a data.frame containing the variables to unsplit.", call. = FALSE)
    if (missing(var.name)) {
        if (any(split.dummies <- vapply(data, function(x) {
            is_not_null(attr(x, "split.var")) && is_not_null(attr(x, "level")) &&
                all(x %in% c(0L, 1L, NA_integer_))
        }, logical(1L)))) {
            var.name <- unique(vapply(data[split.dummies], function(x) attr(x, "split.var"), character(1L)))
        }
        else {
            stop("'var.name' must be a string containing the name of the variables to unsplit.", call. = FALSE)
        }
    }
    else if (!is.character(var.name)) stop("'var.name' must be a string containing the name of the variables to unsplit.", call. = FALSE)
    
    if (is_not_null(sep)) {
        if (is_(sep, "atomic")) {
            if (length(sep) == 1L) sep <- rlang::rep_named(var.name, sep)
            else if (length(sep) == length(var.name)) names(sep) <- var.name
            else {
                stop("'sep' must be a character containing the seperating character in the names of the split variables. See ?unsplitfactor for details.", call. = FALSE)
            }
        }
        else {
            stop("'sep' must be a character containing the seperating character in the names of the split variables. See ?unsplitfactor for details.", call. = FALSE)
        }
    }
    else sep <- rlang::rep_named(var.name, "")
    
    if (is_not_null(dropped.level)) {
        if (is_(dropped.level, "atomic")) {
            if (length(dropped.level) == 1L) dropped.level <- rlang::rep_named(var.name, dropped.level)
            else if (length(dropped.level) == length(var.name)) names(dropped.level) <- var.name
            else {
                stop("'dropped.level' must be an atomic vector containing the value of the dropped category of each split variable. See ?unsplitfactor for details.", call. = FALSE)
            }
        }
        else {
            stop("'dropped.level' must be an atomic vector containing the value of the dropped category of each split variable. See ?unsplitfactor for details.", call. = FALSE)
        }
    }
    else dropped.level <- NULL

    not.the.stem <- character(0)
    
    for (v in var.name) {
        dropped.level0 <- dropped.level[v]
        v.is.split <- any(split.data <- vapply(data, function(x) isTRUE(attr(x, "split.var", TRUE) == v), logical(1L)))
        
        var.to.combine <- if (v.is.split) data[split.data] else data[startsWith(names(data), paste0(v, sep[v]))]
        
        if (is_null(var.to.combine)) {
            not.the.stem <- c(not.the.stem, paste0(v, sep[v]))
            next
        }
        
        if (any(rowSums(apply(var.to.combine, 2, is.na)) %nin% c(0, ncol(var.to.combine)))) {
            stop("The variables in 'data' selected based on 'var.name' and 'sep' do not seem to form a split variable based on the <NA> pattern.", call. = FALSE)
        }
        NA.column <- character(0)
        if (v.is.split && any(na.dummy <- vapply(var.to.combine, function(x) is_not_null(s <- attr(x, "level", TRUE)) && is.na(s), logical(1L)))) {
            dropped.na <- FALSE
        }
        
        if (!rlang::is_bool(dropped.na)) stop("'dropped.na' must be TRUE or FALSE.", call. = FALSE)
        if (!dropped.na) {
            NA.column <- if (v.is.split) {
                names(var.to.combine)[na.dummy]
            } else {
                paste0(v, sep[v], {if (isFALSE(dropped.na)) "NA" else dropped.na})
            }

            if (length(NA.column) > 1) stop(paste0("There appears to be more than one NA variable for ", v, "."), call. = FALSE)
            else if (NA.column %in% names(var.to.combine)) {
                var.to.combine[var.to.combine[[NA.column]] == 1L,] <- NA_integer_
                var.to.combine[[NA.column]] <- NULL
            }
            else {
                stop(paste("There is no variable called", word_list(NA.column, quotes = 2), "to generate the NA values."), call. = FALSE)
            }
        }
        
        var.sum <- rowSums(var.to.combine)

        if (all(is.na(var.sum) | check_if_zero(var.sum - 1))) {
            #Already unsplit
        }
        else if (all(is.na(var.sum) | check_if_zero(var.sum - 1) | check_if_zero(var.sum - 0))) {
            #Missing category
            
            if (is_null(dropped.level[v])) {
                
                k.levels0 <- if (v.is.split) {
                    unlist(lapply(var.to.combine, function(x) attr(x, "level", TRUE)))
                    } else {
                        substr(names(var.to.combine), nchar(paste0(v, sep[v])), nchar(names(var.to.combine)))
                    }
                
                if (can_str2num(k.levels0)) {
                    k.levels0.num <- as.numeric(k.levels0)
                    if (any({nums <- seq(min(k.levels0.num), max(k.levels0.num))} %nin% k.levels0.num)) {
                        dropped.level0 <- nums[nums %nin% k.levels0.num][1]
                    }
                    else if (min(as.numeric(k.levels0)) == 0) dropped.level0 <- max(k.levels0.num) + 1
                    else dropped.level0 <- min(k.levels0.num) - 1
                    dropped.name <- paste0(v, sep[v], dropped.level0)
                }
                else {
                    message(paste("The dropped category for", v, "will be set to NA."))
                    dropped.name <- dropped.level0 <- NA_character_
                }
                
            }
            else dropped.name <- paste0(v, sep[v], dropped.level[v])
            var.to.combine <- setNames(data.frame(1-var.sum, var.to.combine),
                                       c(dropped.name, names(var.to.combine)))
            if (v.is.split) {
                attr(var.to.combine[[1]], "split.var") <- v
                attr(var.to.combine[[1]], "level") <- as.character(dropped.level0)
            }
            
        }
        else {
            stop("The variables in 'data' selected based on 'var.name' and 'sep' do not seem to form a split variable based on the row sums.", call. = FALSE)
        }
        
        k.levels <- if (v.is.split) {
            unlist(lapply(var.to.combine, function(x) attr(x, "level", TRUE))) 
        }
        else {
            substr(names(var.to.combine), 1 + nchar(paste0(v, sep[v])), nchar(names(var.to.combine)))
        }
        
        k <- rep(NA_character_, nrow(data))
        for (i in seq_along(k.levels)) {
            k[var.to.combine[[i]] == 1] <- k.levels[i]
        }
        
        k <- factor(k, levels = k.levels)
        
        if (!rlang::is_bool(replace)) stop("'replace' must be TRUE or FALSE.", call. = FALSE)
        if (replace) {
            where <- which(names(data) %in% c(names(var.to.combine), NA.column))
            
            data[[where[1]]] <- k
            remove.cols <- where[-1]
            if (is_not_null(remove.cols)) data[remove.cols] <- NULL
            names(data)[where[1]] <- v
        }
        else {
            data[[v]] <- k
        }
    }
    
    if (is_not_null(not.the.stem)) warning(paste0(word_list(not.the.stem, is.are = TRUE, quotes = 2), " not the stem of any variables in 'data' and will be ignored. Ensure 'var.name' and 'sep' are correct."), call. = FALSE)
    
    return(data)
}

var.names <- function(b, type, file = NULL, minimal = FALSE) {
    if (is_not_null(attr(b, "print.options")[["co.names"]])) {
        if (!rlang::is_bool(minimal)) stop("'minimal' must be TRUE or FALSE.", call. = FALSE)
        if (minimal) vars <- unique(unlist(lapply(attr(b, "print.options")[["co.names"]], function(x) x[["component"]][x[["type"]] == "base"])))
        else vars <- vapply(attr(b, "print.options")[["co.names"]], function(x) paste(x[["component"]], collapse = ""), character(1))
    }
    else {
        stop("No variable names were found in the object. It is probably not a bal.tab object.", call. = FALSE)
    }
    
    if (is_not_null(file)) {
        if (!endsWith(file, ".csv")) stop("The filename in file must end in \".csv\".", call. = FALSE)
    }
    
    if (missing(type)) {
        if (is_not_null(file)) type <- "df"
        else type <- "vec"
    }
    else {
        type <- match_arg(type, c("df", "vec"))
    }
    
    if (type == "df") {
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
    
    multiple.allowed <- c("stats", "disp", "cluster.fun", "imp.fun")
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
                (!(is.character(opts[[i]]) && is.character(acceptable.options()[[i]]) && (i %in% any.string.allowed || !anyNA(pmatch(opts[[i]], acceptable.options()[[i]])))) &&
                 !all(opts[[i]] %in% acceptable.options()[[i]]))) bad.opts <- c(bad.opts, i)
        }
    }
    
    if (is_not_null(opts)) {
        both.opts <- intersect(multiple.opts, bad.opts)
        multiple.opts <- multiple.opts[multiple.opts %nin% both.opts]
        bad.opts <- bad.opts[bad.opts %nin% both.opts]
        problematic.opts <- make_list(c("multiple", "bad", "both"))
        problematic.opts[["multiple"]] <- setNames(lapply(multiple.opts, function(i) {
            paste(i, "must be of length 1.")
        }), multiple.opts)
        problematic.opts[["bad"]] <- setNames(lapply(bad.opts, function(i) {
            if (i %in% any.string.allowed) paste0(i, " must be a character string.")
            else paste0(i, " must be ", word_list(acceptable.options()[[i]], quotes = 2*is.character(acceptable.options()[[i]]), and.or = "or"), ".")
        }), bad.opts)
        problematic.opts[["both"]] <- setNames(lapply(both.opts, function(i) {
            if (i %in% any.string.allowed) paste0(i, " must be a character string of length 1.")
            else paste0(i, " must be one of ", word_list(acceptable.options()[[i]], quotes = 2*is.character(acceptable.options()[[i]]), and.or = "or"), ".")
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
            stop(paste0(word_list(opts[not.in.accept], is.are = TRUE, quotes = 2),
                        " not", ifelse(plural, "", " an"), " acceptable option", 
                        ifelse(plural, "s", ""), "."), call. = FALSE)
        }
    }
    
    out <- setNames(lapply(paste0("cobalt_", opts), getOption), opts)
    return(out)
    
}