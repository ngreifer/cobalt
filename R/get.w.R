get.w <- function(x, ...) {
    if (!is_(x, "cobalt.processed.obj")) {
        x <- process_obj(x)
        get.w(x, ...)
    }
    else {
        UseMethod("get.w")
    }
}
get.w.matchit <- function(x,...) {
    return(x$weights)
}
get.w.ps <- function(x, stop.method = NULL, estimand, s.weights = FALSE, ...) {

    if (!missing(estimand)) estimand <- tolower(estimand)
    else estimand <- NULL
    
    if (is_not_null(stop.method)) {
        if (any(is.character(stop.method))) {
            rule1 <- names(x$w)[pmatch(tolower(names(x$w)), tolower(stop.method), 0L)]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(names(x$w), and.or = "or", quotes = 2), ".\nUsing all available stop methods instead."))
                rule1 <- names(x$w)
            }
        }
        else if (is.numeric(stop.method) && any(stop.method %in% seq_along(names(x$w)))) {
            if (any(stop.method %nin% seq_along(names(x$w)))) {
                message(paste0("Warning: There are ", length(names(x$w)), " stop methods available, but you requested ", 
                               word_list(stop.method[stop.method %nin% seq_along(names(x$w))], and.or = "and"),"."))
            }
            rule1 <- names(x$w)[stop.method %in% seq_along(names(x$w))]
        }
        else {
            warning("stop.method should be ", word_list(names(x$w), and.or = "or", quotes = 2), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- names(x$w)
        }
    }
    else {
        rule1 <- names(x$w)
    }
    
    s <- names(x$w)[match(tolower(rule1), tolower(names(x$w)))]
    criterion <- substr(tolower(s), 1, nchar(s)-4)
    allowable.estimands <- c("ATT", "ATE", "ATC")
    
    if (is_null(estimand)) estimand <- setNames(substr(toupper(s), nchar(s)-2, nchar(s)), s)
    else if (!all(toupper(estimand) %in% allowable.estimands)) {
        stop(paste0("'estimand' must be ", word_list(allowable.estimands, "or", quotes = 1), "."), call. = FALSE)
    }
    else {
        if (length(estimand) == 1) estimand <- setNames(toupper(rep(estimand, length(s))), s)
        else if (length(estimand) >= length(s)) estimand <- setNames(toupper(estimand[seq_along(s)]), s)
        else stop("'estimand' must be the same length as the number of sets of weights requested.", call. = FALSE)
    }
    
    w <- setNames(as.data.frame(matrix(1, nrow = nrow(x$ps), ncol = length(s))), s)
    for (p in s) {
        if (estimand[p] == "ATT") w[[p]] <- x$treat + (1-x$treat)*x$ps[,p]/(1-x$ps[,p])
        else if (estimand[p] == "ATE") w[[p]] <- x$treat/x$ps[,p] + (1-x$treat)/(1-x$ps[,p])
        else if (estimand[p] == "ATC") w[[p]] <- (1-x$treat) + x$treat*x$ps[,p]/(1-x$ps[,p])
        else w[[p]] <- x$w[,p]
        if (s.weights) w[[p]] <- w[[p]] * x$sampw
    }
    
    names(w) <- ifelse(toupper(substr(s, nchar(s)-2, nchar(s))) == estimand, criterion, paste0(criterion, " (", estimand, ")"))
    if (ncol(w) == 1) w <- w[[1]]
    
    return(w)
}
get.w.mnps <- function(x, stop.method = NULL, s.weights = FALSE, ...) {

    if (is_not_null(stop.method)) {
        if (is.character(stop.method)) {
            rule1 <- x$stopMethods[pmatch(tolower(stop.method), tolower(x$stopMethods), nomatch = 0L)]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(x$stopMethods, and.or = "or", quotes = 2), ".\nUsing all available stop methods instead."))
                rule1 <- x$stopMethods
            }
        }
        else if (is.numeric(stop.method) && any(stop.method %in% seq_along(x$stopMethods))) {
            if (any(stop.method %nin% seq_along(x$stopMethods))) {
                message(paste0("Warning: There are ", length(x$stopMethods), " stop methods available, but you requested ", 
                               word_list(stop.method[stop.method %nin% seq_along(x$stopMethods)], and.or = "and"),"."))
            }
            rule1 <- x$stopMethods[stop.method %in% seq_along(x$stopMethods)]
        }
        else {
            warning("stop.method should be ", word_list(x$stopMethods, and.or = "or", quotes = 2), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- x$stopMethods
        }
    }
    else {
        rule1 <- x$stopMethods
    }
    
    s <- paste.(x$stopMethods[match(tolower(rule1), tolower(x$stopMethods))],
                x$estimand)
    
    estimand <- x$estimand
    criterion <- x$stopMethods[match(tolower(rule1), tolower(x$stopMethods))]
    
    w <- setNames(as.data.frame(matrix(1, nrow = length(x$treatVar), ncol = length(s))),
                  criterion)
    
    if (estimand == "ATT") {
        for (i in x$levExceptTreatATT) {
            if (length(s) > 1) {
                w[x$treatVar == i, criterion] <- get.w.ps(x$psList[[i]])[x$psList[[i]]$treat == FALSE, criterion]
            }
            else {
                w[x$treatVar == i, criterion] <- get.w.ps(x$psList[[i]])[x$psList[[i]]$treat == FALSE]
            }
        }
    }
    else if (estimand == "ATE") {
        for (i in x$treatLev) {
            if (length(s) > 1) {
                w[x$treatVar == i, criterion] <- get.w.ps(x$psList[[i]])[x$psList[[i]]$treat == TRUE, criterion]
            }
            else {
                w[x$treatVar == i, criterion] <- get.w.ps(x$psList[[i]])[x$psList[[i]]$treat == TRUE]
            }
        }
    }
    
    if (s.weights) {
        w <- w * x$sampw
    }
    
    names(w) <- ifelse(toupper(substr(s, nchar(s)-2, nchar(s))) == estimand, criterion, paste0(criterion, " (", estimand, ")"))
    
    if (ncol(w) == 1) w <- w[[1]]
    
    return(w)
}
get.w.ps.cont <- function(x, s.weights = FALSE, ...) {
    
    if (isTRUE(s.weights)) return(x$w * x$sampw)
    else return(x$w)
}
get.w.iptw <- function(x, stop.method = NULL, s.weights = FALSE, ...) {

    if (is_not_null(stop.method)) {
        if (any(is.character(stop.method))) {
            rule1 <- names(x$psList[[1]]$ps)[pmatch(tolower(names(x$psList[[1]]$ps)), tolower(stop.method), 0L)]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word_list(names(x$psList[[1]]$ps), and.or = "or", quotes = 2), ".\nUsing all available stop methods instead."))
                rule1 <- names(x$psList[[1]]$ps)
            }
        }
        else if (is.numeric(stop.method) && any(stop.method %in% seq_along(names(x$psList[[1]]$ps)))) {
            if (any(stop.method %nin% seq_along(names(x$psList[[1]]$ps)))) {
                message(paste0("Warning: There are ", length(names(x$psList[[1]]$ps)), " stop methods available, but you requested ", 
                               word_list(stop.method[stop.method %nin% seq_along(names(x$psList[[1]]$ps))], and.or = "and"),"."))
            }
            rule1 <- names(x$psList[[1]]$ps)[stop.method %in% seq_along(names(x$psList[[1]]$ps))]
        }
        else {
            warning("stop.method should be ", word_list(names(x$psList[[1]]$ps), and.or = "or", quotes = 2), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- names(x$psList[[1]]$ps)
        }
    }
    else {
        rule1 <- names(x$psList[[1]]$ps)
    }
    
    w <- setNames(as.data.frame(matrix(NA, nrow = nrow(x$psList[[1]]$ps),
                                       ncol = length(rule1))),
                  rule1)
    for (i in rule1) {
        w[i] <- Reduce("*", lapply(x$psList, function(x) get.w.ps(x, stop.method = i)))
    }
    
    if (s.weights) {
        w <- w * x$psList[[1]]$sampw
    }
    
    return(w)
}
get.w.Match <- function(x, ...) {
    vapply(seq_len(x$orig.nobs), function(i) {
        sum(x$weights[x$index.treated == i | x$index.control == i])
    }, numeric(1L))
}
get.w.CBPS <- function(x, estimand, ...) {
    A <- list(...)
    use.weights <- if_null_then(A$use.weights, TRUE)
    
    if (!missing(estimand)) estimand <- tolower(estimand)
    else estimand <- NULL
    
    if (is_(x, c("CBPSContinuous", "npCBPS")) || is.factor(x$y)) { #continuous, multi, or npCBPS
        return(x$weights)
    }
    else {
        if (!use.weights) {
            ps <- x$fitted.values
            t <- x$y 
            if (is_null(estimand)) {
                if (all_the_same(x$weights[t == 1])) {
                    estimand <- "att"
                }
                else estimand <- "ate"
            }
            
            estimand <- match_arg(tolower(estimand), c("att", "atc", "ate"))
            if (estimand == "att") {
                return(t + (1-t)*ps/(1-ps))
            }
            if (estimand == "atc") {
                return(t*(1-ps)/ps + (1-t))
            }
            else if (estimand == "ate") {
                return(t/ps + (1-t)/(1-ps))
            }
        }
        else {
            return(x$weights)
        }
        
    }
}
get.w.CBMSM <- function(x, ...) {
    return(x$weights[sort(unique(x$id))])
}
get.w.ebalance <- function(x, treat, ...) {
    if (missing(treat)) stop("'treat' must be specified.", call. = FALSE)
    
    if (!is_(treat, "processed.treat")) treat <- process_treat(treat)
    
    weights <- rep(1, length(treat))
    
    if (length(x$w) != sum(treat == treat_vals(treat)["Control"])) {
        stop("There are more control units in 'treat' than weights in the ebalance object.", call. = FALSE)
    }
    weights[treat == treat_vals(treat)["Control"]] <- x$w
    return(weights)
}
get.w.optmatch <- function(x, estimand, ...) {
    if (missing(estimand) || is_null(estimand)) estimand <- "ATT"
    treat <- as.numeric(attr(x, "contrast.group"))
    return(strata2weights(x, treat = treat, estimand = estimand))
}
get.w.cem.match <- function(x, estimand, ...) {
    A <- list(...)
    if (missing(estimand) || is_null(estimand)) estimand <- "ATT"
    if (isTRUE(A[["use.match.strata"]])) {
        if (is_(x, "cem.match.list")) {
            return(unlist(lapply(x[vapply(x, is_, logical(1L), "cem.match")], function(cm) strata2weights(cm[["mstrata"]], treat = cm[["groups"]], estimand = estimand)), use.names = FALSE))
        }
        else return(strata2weights(x[["mstrata"]], treat = x[["groups"]], estimand = estimand))
    }
    else {
        if (is_(x, "cem.match.list")) {
            return(unlist(lapply(x[vapply(x, is_, logical(1L), "cem.match")], `[[`, "w"), use.names = FALSE))
        }
        else return(x[["w"]])
    }
}
get.w.weightit <- function(x, s.weights = FALSE, ...) {
    if (isTRUE(s.weights)) return(x$weights * x$s.weights)
    else return(x$weights)
}
get.w.designmatch <- function(x, treat, estimand, ...) {
    if (missing(estimand) || is_null(estimand)) estimand <- "ATT"
    if (missing(treat)) stop("'treat' must be specified.", call. = FALSE)
    if (length(x[["group_id"]]) != length(x[["t_id"]]) + length(x[["c_id"]])) {
        stop("designmatch objects without 1:1 matching cannot be used.", call. = FALSE)
    }
    q <- merge(data.frame(id = seq_along(treat)), 
               data.frame(id = c(x[["t_id"]], x[["c_id"]]),
                          group = factor(x[["group_id"]])),
               all.x = TRUE, by = "id")
    q <- q[order(q$id), , drop = FALSE]
    
    return(strata2weights(q$group, treat, estimand))
}
get.w.mimids <- function(x, ...) {
    old_version <- !all(c("object", "models", "approach") %in% names(x))
    
    if (old_version) {
        weights <- unlist(lapply(x[["models"]][-1], get.w.matchit))
    }
    else {
        weights <- unlist(lapply(x[["models"]], get.w.matchit))
    }
    weights[is.na(weights)] <- 0
    return(weights)
}
get.w.wimids <- function(x, ...) {
    old_version <- !all(c("object", "models", "approach") %in% names(x))
    
    if (old_version) {
        weights <- unlist(lapply(x[["models"]][-1], get.w.weightit))
    }
    else {
        weights <- unlist(lapply(x[["models"]], get.w.weightit))
    }
    weights[is.na(weights)] <- 0
    return(weights)
}
get.w.sbwcau <- function(x, ...) {
    return(x[["dat_weights"]][[ncol(x[["dat_weights"]])]])
}