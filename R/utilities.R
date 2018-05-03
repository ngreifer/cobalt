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
splitfactor <- function(data, var.name, replace = TRUE, sep = "_", drop.level = NULL, drop.first = c(TRUE, FALSE, "if2"), drop.singleton = FALSE, drop.na = TRUE, check = TRUE) {
    #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
    #Retains all categories unless only 2 levels, in which case only the second level is retained.
    #If variable only has one level, will delete.
    #var.name= the name of the variable to split when data is specified
    #data=data set to be changed
    
    if (is.data.frame(data)) {
        data <- as.data.frame(data)
        if (check) {
            factor.names <- names(data)[sapply(data, function(x) is.factor(x) || is.character(x))]
            if (missing(var.name)) {
                var.name <- factor.names
            }
            else if (is.character(var.name)) {
                if (any(var.name %in% factor.names)) {
                    if (any(!var.name %in% factor.names)) {
                        not.in.factor.names <- var.name[!var.name %in% factor.names]
                        warning(paste(word.list(not.in.factor.names, "and", is.are = TRUE), 
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
            if (length(factor.names) == 0) {
                stop("There are no factor variables to split in data.", call. = FALSE)
            }
        }
        else {
            if (missing(var.name) || !is.character(var.name)) {
                stop("var.name must be a character vector of the names of variables in data.", call. = FALSE)
            }
            else {
                if (any(var.name %in% names(data))) {
                    if (any(!var.name %in% names(data))) {
                        not.in.data.names <- var.name[!var.name %in% names(data)]
                        warning(paste(word.list(not.in.data.names, "and", is.are = TRUE), 
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
            if (length(var.name) == 0) {
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
    
    if (length(drop.level) > 0 && length(var.name) > 1) {
        warning("drop.level cannot be used with multiple entries to var.name. Ignoring drop.level.", call. = FALSE)
        drop.level <- NULL
    }
    drop.na <- setNames(rep(drop.na, length(var.name)), var.name)
    for (v in var.name) {
        drop <- character(0)
        x <- factor(data[names(data) == v][[1]], exclude = NULL)
        data[names(data) == v][[1]] <- x
        
        skip <- FALSE
        if (nlevels(x) > 1) {
            k <- model.matrix(as.formula(paste0("~`", v, "`- 1")), data = data)
            colnames(k) <- gsub("`", colnames(k), replacement = "")
            
            if (any(is.na(levels(x)))) {
                
                if (drop.na[v]) {
                    k[k[,is.na(levels(x))] == 1,] <- NA
                    #k <- k[, !is.na(levels(x)), drop = FALSE]
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
            
            dropl <- rep(FALSE, ncol(k))
            if (length(drop) > 0) {
                dropl[!is.na(levels(x)) & levels(x) %in% drop] <- TRUE
            }
            if (drop.na[v]) dropl[is.na(levels(x))] <- TRUE
            k <- k[,!dropl]
            
            if (ncol(data) == 1) {
                data <- data.frame(k, row.names = NULL)
            }
            else if (replace) {
                if (match(v, names(data)) == 1){
                    data <- cbind(k, data[names(data)!=v], row.names = NULL)
                }
                else if (match(v, names(data)) == ncol(data)) {
                    data <- cbind(data[names(data)!=v], k, row.names = NULL)
                }
                else {
                    where <- match(v, names(data))
                    data <- cbind(data[1:(where-1)], k, data[(where+1):ncol(data)], row.names = NULL)
                }
            }
            else {
                data <- cbind(data, k, row.names = NULL)
            }
            
        }
        
    }
    
    return(data)
}
unsplitfactor <- function(data, var.name, replace = TRUE, sep = "_", dropped.level = NULL, dropped.na = TRUE) {
    
    if (!is.data.frame(data)) stop("data must be a data.frame containing the variables to unsplit.", call = FALSE)
    if (!is.character(var.name)) stop("var.name must be a string containing the name of the variables to unsplit.", call. = FALSE)
    if (length(dropped.level) > 0 && length(var.name) > 1) {
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
        if (length(var.to.combine) == 0) {
            not.the.stem <- c(not.the.stem, paste0(v, sep))
            next
        }
        
        if (!all(rowSums(apply(var.to.combine, 2, is.na)) %in% c(0, ncol(var.to.combine)))) {
            stop("The variables in data selected based on var.name and sep do not seem to form a split variable based on the <NA> pattern.", call. = FALSE)
        }
        NA.column <- character(0)
        
        if (!isTRUE(dropped.na)) {
            NA.column <- paste0(v, sep, ifelse(dropped.na == FALSE, "NA", dropped.na))
            if (NA.column %in% names(var.to.combine)) {
                var.to.combine[var.to.combine[[NA.column]] == 1,] <- NA
                var.to.combine <- var.to.combine[names(var.to.combine) != NA.column]
            }
            else {
                stop(paste("There is no variable called", word.list(NA.column, quotes = TRUE), "to generate the NA values."), call. = FALSE)
            }
        }
        var.sum <- rowSums(var.to.combine)
        if (isTRUE(all.equal(unique(var.sum), 1))) {
            #Already unsplit
        }
        else if (isTRUE(all.equal(sort(unique(var.sum)), c(0, 1)))) {
            #Missing category
            
            if (length(dropped.level) == 0) {
                k.levels0 <- sapply(names(var.to.combine), function(x) strsplit(x, paste0(v, sep))[[1]][2])
                
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
        
        k.levels <- sapply(names(var.to.combine), function(x) strsplit(x, paste0(v, sep))[[1]][2])
        
        k <- rep(NA_character_, nrow(data))
        for (i in seq_along(k.levels)) {
            k <- ifelse(var.to.combine[[i]] == 1, k.levels[i], k)
        }
        
        k <- factor(k, levels = k.levels)
        
        
        if (replace) {
            where <- which(names(data) %in% c(names(var.to.combine), NA.column))
            
            data[[min(where)]] <- k
            remove.cols <- where[where!=min(where)]
            if (length(remove.cols) > 0) data <- data[-remove.cols]
            names(data)[min(where)] <- v
        }
        else {
            data <- cbind(data, setNames(data.frame(k), v))
        }
    }
    
    if (length(not.the.stem) > 0) warning(paste0(word.list(not.the.stem, is.are = TRUE, quotes = TRUE), " not the stem of any variables in data and will be ignored. Ensure var.name and sep are correct."), call. = FALSE)
    
    return(data)
}

get.w <- function(...) UseMethod("get.w")
get.w.matchit <- function(m,...) {
    return(m$weights)
}
get.w.ps <- function(ps, stop.method = NULL, estimand = NULL, s.weights = FALSE, ...) {
    estimand <- tolower(estimand)
    if (length(stop.method) > 0) {
        if (any(is.character(stop.method))) {
            rule1 <- names(ps$w)[sapply(names(ps$w), function(x) any(startsWith(tolower(x), tolower(stop.method))))]
            if (length(rule1) == 0) {
                message(paste0("Warning: stop.method should be ", word.list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- names(ps$w)
            }
            # rule1 <- tryCatch(match.arg(tolower(stop.method), tolower(names(ps$w)), several.ok = TRUE),
            #                   error = function(cond) {message(paste0("Warning: stop.method should be ", word.list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."));
            #                       return(names(ps$w))})
        }
        else if (is.numeric(stop.method) && any(stop.method %in% seq_along(names(ps$w)))) {
            if (any(!stop.method %in% seq_along(names(ps$w)))) {
                message(paste0("Warning: There are ", length(names(ps$w)), " stop methods available, but you requested ", 
                               word.list(stop.method[!stop.method %in% seq_along(names(ps$w))], and.or = "and"),"."))
            }
            rule1 <- names(ps$w)[stop.method %in% seq_along(names(ps$w))]
        }
        else {
            warning("stop.method should be ", word.list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- names(ps$w)
        }
    }
    else {
        rule1 <- names(ps$w)
    }
    
    s <- names(ps$w)[match(tolower(rule1), tolower(names(ps$w)))]
    
    if (length(estimand) == 0) estimand <- setNames(substr(tolower(s), nchar(s)-2, nchar(s)), s)
    else if (!all(tolower(estimand) %in% c("att", "ate", "atc"))) {
        stop('estimand must be "ATT", "ATE", or "ATC".', call. = FALSE)
    }
    else {
        names(estimand) <- s
    }
    
    w <- setNames(lapply(seq_along(s), function(p) {
        if (estimand[p] == "att") ps$treat + (1-ps$treat)*ps$ps[,s[p]]/(1-ps$ps[,s[p]])
        else if (estimand[p] == "ate") ps$treat/ps$ps[,s[p]] + (1-ps$treat)/(1-ps$ps[,s[p]])
        else if (estimand[p] == "atc") (1-ps$treat) + ps$treat*ps$ps[,s[p]]/(1-ps$ps[,s[p]])}),
        ifelse(tolower(substr(s, nchar(s)-2, nchar(s))) == tolower(estimand), s, paste0(s, " (", toupper(estimand), ")")))
    
    if (s.weights) {
        w <- w * ps$sampw
    }
    
    if (length(w) == 1) w <- w[[1]]
    else {
        class(w) <- "data.frame"; attr(w, "row.names") <- .set_row_names(length(w[[1]]))
    }
    return(w)
}
get.w.mnps <- function(mnps, stop.method = NULL, s.weights = FALSE, ...) {
    if (length(stop.method) > 0) {
        if (any(is.character(stop.method))) {
            rule1 <- mnps$stopMethods[sapply(t(sapply(tolower(stop.method), function(x) startsWith(tolower(mnps$stopMethods), x))), any)]
            if (length(rule1) == 0) {
                message(paste0("Warning: stop.method should be ", word.list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- mnps$stopMethods
            }
        }
        else if (is.numeric(stop.method) && any(stop.method %in% seq_along(mnps$stopMethods))) {
            if (any(!stop.method %in% seq_along(mnps$stopMethods))) {
                message(paste0("Warning: There are ", length(mnps$stopMethods), " stop methods available, but you requested ", 
                               word.list(stop.method[!stop.method %in% seq_along(mnps$stopMethods)], and.or = "and"),"."))
            }
            rule1 <- mnps$stopMethods[stop.method %in% seq_along(mnps$stopMethods)]
        }
        else {
            warning("stop.method should be ", word.list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- mnps$stopMethods
        }
    }
    else {
        rule1 <- mnps$stopMethods
    }
    
    s <- paste(mnps$stopMethods[match(tolower(rule1), tolower(mnps$stopMethods))],
               mnps$estimand, sep = ".")
    
    estimand <- setNames(mnps$estimand, s)
    
    w <- setNames(as.data.frame(matrix(1, nrow = length(mnps$treatVar), ncol = length(s))),
                  s)
    
    if (estimand == "ATT") {
        for (i in mnps$levExceptTreatATT) {
            if (length(s) > 1) {
                w[mnps$treatVar == i, s] <- get.w.ps(mnps$psList[[i]])[mnps$psList[[i]]$treat == FALSE, s]
            }
            else {
                w[mnps$treatVar == i, s] <- get.w.ps(mnps$psList[[i]])[mnps$psList[[i]]$treat == FALSE]
            }
        }
    }
    else if (estimand == "ATE") {
        for (i in mnps$treatLev) {
            if (length(s) > 1) {
                w[mnps$treatVar == i, s] <- get.w.ps(mnps$psList[[i]])[mnps$psList[[i]]$treat == TRUE, s]
            }
            else {
                w[mnps$treatVar == i, s] <- get.w.ps(mnps$psList[[i]])[mnps$psList[[i]]$treat == TRUE]
            }
        }
    }
    
    if (s.weights) {
        w <- w * mnps$sampw
    }
    
    if (ncol(w) == 1) w <- w[[1]]
    
    return(w)
}
get.w.iptw <- function(iptw, stop.method = NULL, s.weights = FALSE, ...) {
    if (length(stop.method) > 0) {
        if (any(is.character(stop.method))) {
            rule1 <- names(iptw$psList[[1]]$ps)[sapply(names(iptw$psList[[1]]$ps), function(x) any(startsWith(tolower(x), tolower(stop.method))))]
            if (length(rule1) == 0) {
                message(paste0("Warning: stop.method should be ", word.list(names(iptw$psList[[1]]$ps), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- names(iptw$psList[[1]]$ps)
            }
        }
        else if (is.numeric(stop.method) && any(stop.method %in% seq_along(names(iptw$psList[[1]]$ps)))) {
            if (any(!stop.method %in% seq_along(names(iptw$psList[[1]]$ps)))) {
                message(paste0("Warning: There are ", length(names(iptw$psList[[1]]$ps)), " stop methods available, but you requested ", 
                               word.list(stop.method[!stop.method %in% seq_along(names(iptw$psList[[1]]$ps))], and.or = "and"),"."))
            }
            rule1 <- names(iptw$psList[[1]]$ps)[stop.method %in% seq_along(names(iptw$psList[[1]]$ps))]
        }
        else {
            warning("stop.method should be ", word.list(names(iptw$psList[[1]]$ps), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- names(iptw$psList[[1]]$ps)
        }
    }
    else {
        rule1 <- names(iptw$psList[[1]]$ps)
    }
    
    w <- setNames(as.data.frame(matrix(NA, nrow = nrow(iptw$psList[[1]]$ps),
                                       ncol = length(rule1))),
                  rule1)
    for (i in rule1) {
        w[i] <- Reduce("*", lapply(iptw$psList, function(x) get.w.ps(x, stop.method = i)))
    }
    
    if (s.weights) {
        w <- w * iptw$psList[[1]]$sampw
    }
    
    return(w)
}
get.w.Match <- function(M,  ...) {
    nobs <- M$orig.nobs
    weights.list <- index.list <- setNames(vector("list", 4), c("control", "treated", "unmatched", "dropped"))
    
    index.list$control <- seq_len(nobs)[seq_len(nobs) %in% M$index.control]
    index.list$treated <- seq_len(nobs)[seq_len(nobs) %in% M$index.treated]
    index.list$unmatched <- seq_len(nobs)[!seq_len(nobs) %in% c(M$index.treated, M$index.control, M$index.dropped)]
    index.list$dropped <- seq_len(nobs)[seq_len(nobs) %in% M$index.dropped]
    
    weights.list$control <- weights.list$treated <- M$weights
    weights.list$unmatched <- rep(0, sum(!seq_len(nobs) %in% c(M$index.treated, M$index.control, M$index.dropped)))
    weights.list$dropped <- rep(0, length(M$index.dropped))
    
    data.list <- lapply(1:4, function(x) cbind(data.frame(index=index.list[[x]]), data.frame(weights=weights.list[[x]])))
    o.data <- do.call(rbind, data.list)
    o.data2 <- merge(unique(o.data[is.na(match(names(o.data), "weights"))]), 
                     aggregate(weights~index, data=o.data, FUN=sum), 
                     by="index")
    return(o.data2$weights)
}
get.w.CBPS <- function(c, estimand = NULL, ...) {
    A <- list(...)
    if (length(A$use.weights) == 0) use.weights <- TRUE
    else use.weights <- A$use.weights
    
    estimand <- tolower(estimand)
    
    if ("CBPSContinuous" %in% class(c) || is.factor(c$y)) { #continuous
        return(c$weights)
    }
    else {
        if (!use.weights) {
            ps <- c$fitted.values
            t <- c$y 
            if (length(estimand) == 0) {
                if (abs(max(c$weights[t == 1], na.rm = TRUE) - 
                        min(c$weights[t == 1], na.rm = TRUE)) < 
                    sqrt(.Machine$double.eps)) {
                    estimand <- "att"
                }
                else estimand <- "ate"
            }
            
            estimand <- match.arg(tolower(estimand), c("att", "atc", "ate"))
            if (estimand == "att") {
                return(ifelse(t == 1, 1, ps/(1-ps)))
            }
            if (estimand == "atc") {
                return(ifelse(t == 1, (1-ps)/ps, 1))
            }
            else if (estimand == "ate") {
                return(ifelse(t == 1, 1/ps, 1/(1-ps)))
            }
        }
        else {
            return(c$weights)
        }
        
    }
}
get.w.npCBPS <- function(c, estimand = NULL, ...) {
    return(c$weights)
}
get.w.CBMSM <- function(c, ...) {
    return(c$weights)
}
get.w.ebalance <- function(e, treat, ...) {
    if (missing(treat)) stop("treat must be specified.", call. = FALSE)
    
    weights <- rep(1, length(treat))

    if (length(e$w) != sum(treat == 0)) {
        stop("There are more control units in treat than weights in the ebalance object.", call. = FALSE)
    }
    weights[treat == 0] <- e$w
    return(weights)
}
get.w.ebalance.trim <- get.w.ebalance
get.w.optmatch <- function(o, ...) {
    treat <- as.numeric(attr(o, "contrast.group"))
    return(match.strata2weights(o, treat = treat, covs = NULL))
}
get.w.weightit <- function(W, s.weights = FALSE, ...) {
    if (s.weights) return(W$weights * W$s.weights)
    else return(W$weights)
}

#For cobalt
word.list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
    #When given a vector of strings, creates a string of the form "a and b"
    #or "a, b, and c"
    #If is.are, adds "is" or "are" appropriately
    L <- length(word.list)
    if (quotes) word.list <- sapply(word.list, function(x) paste0("\"", x, "\""))
    if (L == 0) {
        out <- ""
        attr(out, "plural") = FALSE
    }
    else {
        word.list <- word.list[!word.list %in% c(NA, "")]
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
nunique <- function(x, nmax = NA, na.rm = TRUE) {
    if (length(x) == 0) return(0)
    else {
        if (na.rm) x <- x[!is.na(x)]
        if (is.factor(x)) return(nlevels(x))
        else return(length(unique(x, nmax = nmax)))
    }
    
}
nunique.gt <- function(x, n, na.rm = TRUE) {
    if (n < 0) stop("n must be non-negative.", call. = FALSE)
    if (length(x) == 0) FALSE
    else {
        if (na.rm) x <- x[!is.na(x)]
        if (length(x) < 2000) nunique(x) > n
        else tryCatch(nunique(x, nmax = n) > n, error = function(e) TRUE)        
    }
}
is.formula <- function(f, sides = NULL) {
    res <- is.name(f[[1]])  && deparse(f[[1]]) %in% c( '~', '!') &&
        length(f) >= 2
    if (length(sides) > 0 && is.numeric(sides) && sides %in% c(1,2)) {
        res <- res && length(f) == sides + 1
    }
    return(res)
}