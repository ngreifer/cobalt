bal.plot <- function(obj, var.name, ..., un = FALSE, which.sub = NULL, cluster = NULL, which.cluster = NULL, imp = NULL, which.imp = NULL) {
    
    args <- list(...)
    if (missing(var.name)) {
        stop("An argument for var.name must be specified.", call. = FALSE)
    }
    
    X <- x2base(obj, ..., cluster = cluster, imp = imp, s.d.denom = "treated") #s.d.denom to avoid x2base warning
    
    if (var.name %in% names(X$covs)) var <- X$covs[, var.name]
    #else if (!is.null(args$data) && var.name %in% names(args$data)) var <- args$data[, var.name]
    else if (!is.null(X$addl) && var.name %in% names(X$addl)) var <- X$addl[, var.name]
    else if (!is.null(args$addl) && var.name %in% names(args$addl)) var <- args$addl[, var.name]
    else if (!is.null(X$distance) && var.name %in% names(X$distance)) var <- X$distance[, var.name]
    else if (!is.null(args$distance) && var.name %in% names(args$distance)) var <- args$distance[, var.name]
    else stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."))
    
    title <- paste0("Distributional Balance for \"", var.name, "\"")
    subtitle <- "Adjusted Sample"
    
    facet <- NULL
    if (length(X$subclass) > 0  && !isTRUE(un)) {
        if (length(X$cluster) > 0) stop("Subclasses are not supported with clusters.", call. = FALSE)
        if (length(X$imp) > 0) stop("Subclasses are not supported with multiple imputations.", call. = FALSE)
        if (is.null(which.sub)) { #display all subs
            X$weights <- X$weights[!is.na(X$subclass)]
            X$treat <- X$treat[!is.na(X$subclass)]
            var <- var[!is.na(X$subclass)]
            subclass <- paste("subclass", X$subclass[!is.na(X$subclass)])
            subtitle <- NULL
            #title <- paste0(title, "\nacross subclasses")
        }
        else {
            if (is.numeric(which.sub) && length(which.sub) == 1) {
                if (which.sub %in% levels(X$subclass)) {
                    X$weights <- X$weights[!is.na(X$subclass) & X$subclass==which.sub]
                    X$treat <- X$treat[!is.na(X$subclass) & X$subclass==which.sub]
                    var <- var[!is.na(X$subclass) & X$subclass==which.sub]
                    subclass <- paste("subclass", X$subclass[!is.na(X$subclass) & X$subclass==which.sub])
                    #title <- paste0(title, "\nin subclass ", which.sub)
                    subtitle <- NULL
                    if (isTRUE(un)) {
                        message("which.sub specified and un set to TRUE; setting un = FALSE.")
                        un <- FALSE
                    }
                }
                else stop(paste0("\"", which.sub, "\" does not correspond to a subclass in the object."), call. = FALSE)
                
            }
            else stop("The argument to which.sub must be a single number corresponding to the subclass for which distributions are to be displayed.", call. = FALSE)
        }
        facet <- "subclass"
    }
    else if (length(which.sub) > 0) {
        warning("which.sub was specified but no subclasses were supplied. Ignoring which.sub.", call. = FALSE)
    }
    
    #NULL: all; NA: none
    in.imp <- NULL
    if (length(X$imp) > 0) {
        if (length(which.imp) == 0) {
            in.imp <- !is.na(X$imp)
        }
        else if (all(!is.na(which.imp))) {
            if (is.numeric(which.imp)) {
                if (all(which.imp %in% seq_len(nlevels(X$imp)))) {
                    in.imp <- !is.na(X$imp) & sapply(X$imp, function(x) !is.na(match(x, levels(X$imp)[which.imp])))
                }
                else {
                    stop(paste0("The following inputs to which.imp do not correspond to given imputations:\n\t", word.list(which.imp[!which.imp %in% seq_len(nlevels(X$imp))])), call. = FALSE)
                }
            }
            else stop("The argument to which.imp must be the indices corresponding to the imputations for which distributions are to be displayed.", call. = FALSE)
        }
    }
    else if (length(which.imp) > 0) {
        warning("which.imp was specified but no imp values were supplied. Ignoring which.imp.", call. = FALSE)
    }
    
    in.cluster <- NULL
    if (length(X$cluster) > 0) {
        if (length(which.cluster) == 0) {
            in.cluster <- !is.na(X$cluster)
        }
        else if (all(!is.na(which.cluster))) {
            if (is.numeric(which.cluster)) {
                if (all(which.cluster %in% seq_len(nlevels(X$cluster)))) {
                    in.cluster <- !is.na(X$cluster) & sapply(X$cluster, function(x) !is.na(match(x, levels(X$cluster)[which.cluster])))
                }
                else {
                    stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word.list(which.cluster[!which.cluster %in% seq_len(nlevels(X$cluster))])), call. = FALSE)
                }
            }
            else if (is.character(which.cluster)) {
                if (all(!is.na(match(which.cluster, levels(X$cluster))))) {
                    in.cluster <- !is.na(X$cluster) & sapply(X$cluster, function(x) !is.na(match(x, which.cluster)))
                }
                else {
                    stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word.list(which.cluster[is.na(match(which.cluster, levels(X$cluster)))])), call. = FALSE)
                }
            }
            else stop("The argument to which.cluster must be the names or indices corresponding to the clusters for which distributions are to be displayed.", call. = FALSE)
        }
    }
    else if (length(which.cluster) > 0) {
        warning("which.cluster was specified but no cluster values were supplied. Ignoring which.cluster.", call. = FALSE)
    }
    
    if (length(unique(X$imp[in.imp])) > 1 && length(unique(X$cluster[in.cluster])) > 1) {
        stop("Only either one imputation or one cluster may be specified at a time. Ensure either which.imp or which.cluster are given one (non-NULL, non-NA) value.", call. = FALSE)
    }
    
    if (length(in.imp) > 0 && length(in.cluster) > 0) {
        X$weights <- X$weights[in.imp & in.cluster]
        X$treat <- X$treat[in.imp & in.cluster]
        var <- var[in.imp & in.cluster]
        imp <- paste("Imputation", X$imp[in.imp & in.cluster])
        cluster <- X$cluster[in.imp & in.cluster]
        facet <- c("imp", "cluster")
    }
    else if (length(in.imp) > 0) {
        X$weights <- X$weights[in.imp]
        X$treat <- X$treat[in.imp]
        var <- var[in.imp]
        imp <- paste("Imputation", X$imp[in.imp])
        facet <- "imp"
    }
    else if (length(in.cluster) > 0) {
        X$weights <- X$weights[in.cluster]
        X$treat <- X$treat[in.cluster]
        var <- var[in.cluster]
        cluster <- X$cluster[in.cluster]
        facet <- "cluster"
    }
    
    if (is.null(X$weights) || isTRUE(un)) {
        X$weights <- rep(1, length(X$treat))
        subtitle <- "Unadjusted Sample"
    }
    
    if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) { #Continuous treatments
        treat <- X$treat
        is.categorical.var <- length(unique(var)) <= 2 || is.factor(var) || is.character(var)
        if (identical(facet,"subclass")) {
            if (is.categorical.var) {
                weights <- mapply(function(w, s, v) w[var==v] / sum(X$weights[subclass==s & var==v]), w = X$weights, s = subclass, v = var)
            }
            else {
                weights <- mapply(function(w, s) w / sum(X$weights[subclass==s]), w = X$weights, s = subclass)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, subclass = subclass)
        }
        else if (identical(facet,"cluster")) {
            if (is.categorical.var) {
                weights <- mapply(function(w, c, v) w / sum(X$weights[cluster==c & var==v]), w = X$weights, c = cluster, v = var)
            }
            else {
                weights <- mapply(function(w, c) w / sum(X$weights[cluster==c]), w = X$weights, c = cluster)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, cluster = cluster)
        }
        else if ((identical(facet,"imp"))) {
            if (is.categorical.var) {
                weights <- mapply(function(w, i, v) w / sum(X$weights[imp==i & var==v]), w = X$weights, i = imp, v = var)
            }
            else {
                weights <- mapply(function(w, i) w / sum(X$weights[imp==i]), w = X$weights, i = imp)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, imp = imp)
        }
        else if (identical(facet,c("imp", "cluster"))) {
            if (is.categorical.var) {
                weights <- mapply(function(w, i, c, v) w / sum(X$weights[imp==i & cluster == c & var==v]), w = X$weights, i = imp, c = cluster, v = var)
            }
            else {
                weights <- mapply(function(w, i, c) w / sum(X$weights[imp==i & cluster == c]), w = X$weights, i = imp, c = cluster)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, imp = imp, cluster = cluster)
        }
        else {
            if (is.categorical.var) {
                weights <- mapply(function(w, v) w / sum(X$weights[var==v]), w = X$weights, v = var)
            }
            else {
                weights <- X$weights
            }
            d <- data.frame(weights = weights, treat = treat, var = var)
        }
        if (is.categorical.var) { #Categorical vars
            d$var <- factor(d$var)
            bp <- ggplot(d, mapping = aes(treat, fill = var, weight = weights)) + 
                geom_density(alpha = .4) + 
                labs(fill = var.name, y = "Density", x = "Treat", title = title, subtitle = subtitle)
        }
        else { #Continuous vars
            bp <- ggplot(d, mapping = aes(x = var, y = treat, weight = weights)) + 
                geom_point() + geom_smooth(method = "loess", se = FALSE, alpha = .1) + geom_smooth(method = "lm", se = FALSE, linetype = 2, alpha = .4) + 
                geom_hline(yintercept = w.m(treat, weights), linetype = 1, alpha = .9) + 
                labs(y = "Treat", x = var.name, title = title, subtitle = subtitle)
        }
    }
    else { #Categorical treatments (multinomial supported)
        treat <- factor(X$treat)
        
        if (length(facet) > 0) {
            if (identical(facet,"subclass")) {
                weights <- mapply(function(w, t, s) w / sum(X$weights[X$treat==t & subclass==s]), w = X$weights, t = X$treat, s = subclass)
                d <- data.frame(weights = weights, treat = treat, var = var, subclass = subclass)
            }
            else if (identical(facet,"cluster")) {
                weights <- mapply(function(w, t, c) w / sum(X$weights[X$treat==t & cluster==c]), w = X$weights, t = X$treat, c = cluster)
                d <- data.frame(weights = weights, treat = treat, var = var, cluster = cluster)
            }
            else if (identical(facet,"imp")){
                weights <- mapply(function(w, t, i) w / sum(X$weights[X$treat==t & imp==i]), w = X$weights, t = X$treat, i = imp)
                d <- data.frame(weights = weights, treat = treat, var = var, imp = imp)
            }
            else if (identical(facet, c("imp", "cluster"))) {
                weights <- mapply(function(w, t, c, i) w / sum(X$weights[X$treat==t & cluster==c & imp==i]), w = X$weights, t = X$treat, c= cluster, i = imp)
                d <- data.frame(weights = weights, treat = treat, var = var, cluster = cluster, imp = imp)
            }
        }
        else {
            weights <- mapply(function(w, t) w / sum(X$weights[X$treat==t]), w = X$weights, t = X$treat)
            d <- data.frame(weights = weights, treat = treat, var = var)
        }
        
        if (length(unique(var)) <= 2 || is.factor(var) || is.character(var)) { #Categorical vars
            d$var <- factor(d$var)
            bp <- ggplot(d, mapping = aes(var, fill = treat, weight = weights)) + 
                geom_bar(position = "dodge", alpha = .4, color = "black") + 
                labs(x = var.name, y = "Proportion", fill = "Treat", title = title, subtitle = subtitle) + 
                scale_x_discrete(drop=FALSE) + scale_fill_discrete(drop=FALSE)
        }
        else { #Continuous vars
            bp <- ggplot(d, mapping = aes(var, fill = treat, weight = weights)) + 
                geom_density(alpha = .4) + 
                labs(x = var.name, y = "Density", fill = "Treat", title = title, subtitle = subtitle)
        }
    }
    
    if (length(facet) > 0) {
        bp <- bp + facet_grid(f.build(".", facet), drop = FALSE)
    }
    
    return(bp)
}