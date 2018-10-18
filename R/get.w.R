get.w <- function(...) UseMethod("get.w")
get.w.matchit <- function(m,...) {
    return(m$weights)
}
get.w.ps <- function(ps, stop.method = NULL, estimand = NULL, s.weights = FALSE, ...) {
    estimand <- tolower(estimand)
    if (is_not_null(stop.method)) {
        if (any(is.character(stop.method))) {
            rule1 <- names(ps$w)[vapply(names(ps$w), function(x) any(startsWith(tolower(x), tolower(stop.method))), logical(1L))]
            if (is_null(rule1)) {
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
    
    if (is_null(estimand)) estimand <- setNames(substr(tolower(s), nchar(s)-2, nchar(s)), s)
    else if (!all(tolower(estimand) %in% c("att", "ate", "atc"))) {
        stop('estimand must be "ATT", "ATE", or "ATC".', call. = FALSE)
    }
    else {
        names(estimand) <- s
    }
    
    w <- setNames(as.data.frame(matrix(1, nrow = nrow(ps$ps), ncol = length(s))),
                  ifelse(tolower(substr(s, nchar(s)-2, nchar(s))) == tolower(estimand), s, paste0(s, " (", toupper(estimand), ")")))
    for (p in s) {
        if (estimand[p] == "att") w[[p]] <- ps$treat + (1-ps$treat)*ps$ps[,p]/(1-ps$ps[,p])
        else if (estimand[p] == "ate") w[[p]] <- ps$treat/ps$ps[,p] + (1-ps$treat)/(1-ps$ps[,p])
        else if (estimand[p] == "atc") w[[p]] <- (1-ps$treat) + ps$treat*ps$ps[,p]/(1-ps$ps[,p])
        else w[[p]] <- ps$w[,p]
        if (s.weights) w[[p]] <- w[[p]] * ps$sampw
    }
    
    if (ncol(w) == 1) w <- w[[1]]
    
    return(w)
}
get.w.mnps <- function(mnps, stop.method = NULL, s.weights = FALSE, ...) {
    if (is_not_null(stop.method)) {
        if (any(is.character(stop.method))) {
            rule1 <- mnps$stopMethods[sapply(t(sapply(tolower(stop.method), function(x) startsWith(tolower(mnps$stopMethods), x))), any)]
            if (is_null(rule1)) {
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
    
    s <- paste.(mnps$stopMethods[match(tolower(rule1), tolower(mnps$stopMethods))],
                mnps$estimand)
    
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
get.w.ps.cont <- function(ps.cont, stop.method = NULL, s.weights = FALSE, ...) {
    if (is_not_null(stop.method)) {
        if (any(is.character(stop.method))) {
            rule1 <- names(ps.cont$w)[vapply(names(ps.cont$w), function(x) any(startsWith(tolower(x), tolower(stop.method))), logical(1L))]
            if (is_null(rule1)) {
                message(paste0("Warning: stop.method should be ", word.list(names(ps.cont$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
                rule1 <- names(ps.cont$w)
            }
            # rule1 <- tryCatch(match.arg(tolower(stop.method), tolower(names(ps$w)), several.ok = TRUE),
            #                   error = function(cond) {message(paste0("Warning: stop.method should be ", word.list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."));
            #                       return(names(ps$w))})
        }
        else if (is.numeric(stop.method) && any(stop.method %in% seq_along(names(ps.cont$w)))) {
            if (any(!stop.method %in% seq_along(names(ps.cont$w)))) {
                message(paste0("Warning: There are ", length(names(ps.cont$w)), " stop methods available, but you requested ", 
                               word.list(stop.method[!stop.method %in% seq_along(names(ps.cont$w))], and.or = "and"),"."))
            }
            rule1 <- names(ps.cont$w)[stop.method %in% seq_along(names(ps.cont$w))]
        }
        else {
            warning("stop.method should be ", word.list(names(ps.cont$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
            rule1 <- names(ps.cont$w)
        }
    }
    else {
        rule1 <- names(ps.cont$w)
    }
    
    s <- names(ps.cont$w)[match(tolower(rule1), tolower(names(ps.cont$w)))]
    
    w <- setNames(as.data.frame(matrix(1, nrow = nrow(ps.cont$w), ncol = length(s))),
                  s)
    
    for (p in s) {
        w[[p]] <- ps.cont$w[[p]]
        if (!s.weights && is_not_null(ps.cont$sampw)) w[[p]] <- w[[p]] / ps.cont$sampw
    }
    
    if (ncol(w) == 1) w <- w[[1]]
    
    return(w)
}
get.w.iptw <- function(iptw, stop.method = NULL, s.weights = FALSE, ...) {
    if (is_not_null(stop.method)) {
        if (any(is.character(stop.method))) {
            rule1 <- names(iptw$psList[[1]]$ps)[vapply(names(iptw$psList[[1]]$ps), function(x) any(startsWith(tolower(x), tolower(stop.method))), logical(1L))]
            if (is_null(rule1)) {
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
    if (is_null(A$use.weights)) use.weights <- TRUE
    else use.weights <- A$use.weights
    
    estimand <- tolower(estimand)
    
    if ("CBPSContinuous" %in% class(c) || is.factor(c$y)) { #continuous
        return(c$weights)
    }
    else {
        if (!use.weights) {
            ps <- c$fitted.values
            t <- c$y 
            if (is_null(estimand)) {
                if (nunique.gt(c$weights[t == 1], 1)) {
                    estimand <- "ate"
                }
                else estimand <- "att"
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
get.w.npCBPS <- function(c, ...) {
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
get.w.designmatch <- function(dm, treat, ...) {
    if (missing(treat)) stop("treat must be specified.", call. = FALSE)
    if (length(dm[["group_id"]]) != length(dm[["t_id"]]) + length(dm[["c_id"]])) {
        ratio <- length(dm[["c_id"]])/length(dm[["t_id"]])
        if (check_if_zero(ratio - as.integer(ratio))) {
            dm[["group_id"]] <- c(seq_along(dm[["t_id"]]),
                                  rep(seq_along(dm[["c_id"]]), each = as.integer(ratio)))
        }
        else {
            stop("There is a problem with the group_id value in the designmatch output. Matched sets cannot be determined.", call. = FALSE)
        }
    }
    q <- merge(data.frame(id = seq_along(treat)), 
               data.frame(id = c(dm[["t_id"]], dm[["c_id"]]),
                          group = dm[["group_id"]]),
               all.x = TRUE, by = "id")
    q <- q[order(q$id), , drop = FALSE]
    
    return(match.strata2weights(q$group, treat))
}