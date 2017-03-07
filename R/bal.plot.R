bal.plot <- function(obj, var.name, ..., un = FALSE, which.sub = NULL, cluster = NULL, which.cluster = NULL) {
    
    args <- list(...)
    
    X <- x2base(obj, ..., cluster = cluster, std.ok = TRUE)
    
    #Functions
    w.m <- function(x, w = NULL) {
        if (is.null(w)) w <- rep(1, length(x)); return(sum(x*w, na.rm=TRUE)/sum(w, na.rm=TRUE))
    }
    
    if (var.name %in% names(X$covs)) var <- X$covs[, var.name]
    else if (!is.null(args$data) && var.name %in% names(args$data)) var <- args$data[, var.name]
    else if (!is.null(X$addl) && var.name %in% names(X$addl)) var <- args$addl[, var.name]
    else if (!is.null(args$addl) && var.name %in% names(args$addl)) var <- args$addl[, var.name]
    else if (var.name == ".distance" && !is.null(X$distance)) var <- X$distance
    else stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."))
    
    title <- paste0("Distributional Balance for \"", var.name, "\"")
    subtitle <- "Adjusted sample"
    
    facet <- "none"
    if (length(X$subclass) > 0  && !isTRUE(un)) {
        if (length(X$cluster) > 0) stop("Subclasses are not supported with clusters.", call. = FALSE)
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
    
    if (length(X$cluster) > 0) {
        if (length(which.cluster) == 0) {
            X$weights <- X$weights[!is.na(X$cluster)]
            X$treat <- X$treat[!is.na(X$cluster)]
            var <- var[!is.na(X$cluster)]
            cluster <- paste("cluster", X$cluster[!is.na(X$cluster)])
            facet <- "cluster"
        }
        else if (length(which.cluster) > 0) {
            if (!is.na(which.cluster)) {
                if (is.numeric(which.cluster)) {
                    if (all(which.cluster %in% seq_len(nlevels(X$cluster)))) {
                        in.cluster <- sapply(X$cluster, function(x) !is.na(match(x, levels(X$cluster)[which.cluster])))
                    }
                    else {
                        stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word.list(which.cluster[!which.cluster %in% seq_len(nlevels(X$cluster))])), call. = FALSE)
                    }
                }
                else if (is.character(which.cluster)) {
                    if (all(!is.na(match(which.cluster, levels(X$cluster))))) {
                        in.cluster <- sapply(X$cluster, function(x) !is.na(match(x, which.cluster)))
                    }
                    else {
                        stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word.list(which.cluster[is.na(match(which.cluster, levels(X$cluster)))])), call. = FALSE)
                    }
                }
                else stop("The argument to which.cluster must be the names or indices corresponding to the clusters for which distributions are to be displayed.", call. = FALSE)
                
                X$weights <- X$weights[!is.na(X$cluster) & in.cluster]
                X$treat <- X$treat[!is.na(X$cluster) & in.cluster]
                var <- var[!is.na(X$cluster) & in.cluster]
                cluster <- paste("cluster", X$cluster[!is.na(X$cluster) & in.cluster])
                facet <- "cluster"
            }
        }
    }
    else if (length(which.cluster) > 0) {
        warning("which.cluster was specified but no cluster values were supplied. Ignoring which.cluster.", call. = FALSE)
    }
    
    if (is.null(X$weights) || isTRUE(un)) {
        X$weights <- rep(1, length(X$treat))
        subtitle <- "Unadjusted Sample"
    }
    
    if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) { #Continuous treatments
        treat <- X$treat
        is.categorical.var <- length(unique(var)) <= 2 || is.factor(var) || is.character(var)
        if (facet == "subclass") {
            if (is.categorical.var) {
                weights <- mapply(function(w, s, v) w[var==v] / sum(X$weights[subclass==s & var==v]), w = X$weights, s = subclass, v = var)
            }
            else {
                weights <- mapply(function(w, s) w / sum(X$weights[subclass==s]), w = X$weights, s = subclass)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, subclass = subclass)
        }
        else if (facet == "cluster") {
            if (is.categorical.var) {
                weights <- mapply(function(w, c, v) w / sum(X$weights[cluster==c & var==v]), w = X$weights, c = cluster, v = var)
            }
            else {
                weights <- mapply(function(w, c) w / sum(X$weights[cluster==c]), w = X$weights, c = cluster)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, cluster = cluster)
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
    # else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
    #     stop("Multinomial treaments are not yet supported.", call. = FALSE)
    # }
    else { #Categorical treatments (multinomial supported)
        #treat <- factor(binarize(X$treat))
        treat <- factor(X$treat)
        if (facet == "subclass") {
            weights <- mapply(function(w, t, s) w / sum(X$weights[X$treat==t & subclass==s]), w = X$weights, t = X$treat, s = subclass)
            d <- data.frame(weights = weights, treat = treat, var = var, subclass = subclass)
        }
        else if (facet == "cluster") {
            weights <- mapply(function(w, t, c) w / sum(X$weights[X$treat==t & cluster==c]), w = X$weights, t = X$treat, c = cluster)
            d <- data.frame(weights = weights, treat = treat, var = var, cluster = cluster)
        }
        else {
            weights <- mapply(function(w, t) w / sum(X$weights[X$treat==t]), w = X$weights, t = X$treat)
            d <- data.frame(weights = weights, treat = treat, var = var)
        }
        
        if (length(unique(var)) <= 2 || is.factor(var) || is.character(var)) { #Categorical vars
            d$var <- factor(d$var)
            bp <- ggplot(d, mapping = aes(var, fill = treat, weight = weights)) + 
                geom_bar(position = "dodge", alpha = .4, color = "black") + 
                labs(x = var.name, y = "Proportion", fill = "Treat", title = title, subtitle = subtitle) 
        }
        else { #Continuous vars
            bp <- ggplot(d, mapping = aes(var, fill = treat, weight = weights)) + 
                geom_density(alpha = .4) + 
                labs(x = var.name, y = "Density", fill = "Treat", title = title, subtitle = subtitle)
        }
    }
    if (facet == "subclass") {
        bp <- bp + facet_wrap(~subclass)
    }
    else if (facet == "cluster") {
        bp <- bp + facet_wrap(~cluster)
    }
    return(bp)
}