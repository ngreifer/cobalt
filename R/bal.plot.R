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
    
    if (length(X$subclass) > 0) {
        if (length(X$cluster) > 0) stop("Subclasses are not supported with clusters.", call. = FALSE)
        if (!is.null(which.sub)) {
            if (is.numeric(which.sub) && length(which.sub) == 1) {
                if (which.sub %in% levels(X$subclass)) {
                    X$weights <- X$weights[!is.na(X$subclass) & X$subclass==which.sub]
                    X$treat <- X$treat[!is.na(X$subclass) & X$subclass==which.sub]
                    var <- var[!is.na(X$subclass) & X$subclass==which.sub]
                    title <- paste0(title, "\nin subclass ", which.sub)
                }
                else stop(paste0("\"", which.sub, "\" does not correspond to a subclass in the object."), call. = FALSE)
                
            }
            else stop("The argument to which.sub must be a single number corresponding to the subclass for which distributions are to be displayed.", call. = FALSE)
        }
        else stop("Argument contains subclasses but no which.sub value was supplied.", call. = FALSE)
    }
    else if (length(which.sub) > 0) {
        warning("which.sub was specified but no subclasses were supplied. Ignoring which.sub.", call. = FALSE)
    }
    
    if (length(X$cluster) > 0) {
        if (length(which.cluster) > 0) {
            if (!is.na(which.cluster) && length(which.cluster) == 1) {
                if (any(sapply(levels(X$cluster), function(x) isTRUE(all.equal(x, which.cluster))))) {
                    X$weights <- X$weights[!is.na(X$cluster) & sapply(X$cluster, function(x) isTRUE(all.equal(as.character(x), which.cluster)))]
                    X$treat <- X$treat[!is.na(X$cluster) & sapply(X$cluster, function(x) isTRUE(all.equal(as.character(x), which.cluster)))]
                    var <- var[!is.na(X$cluster) & sapply(X$cluster, function(x) isTRUE(all.equal(as.character(x), which.cluster)))]
                    title <- paste0(title, "\nin cluster ", which.cluster)
                }
                else if (any(sapply(unique(as.numeric(X$cluster)), function(x) isTRUE(all.equal(x, which.cluster))))) {
                    X$weights <- X$weights[!is.na(X$cluster) & sapply(X$cluster, function(x) isTRUE(all.equal(as.character(x), levels(X$cluster)[which.cluster])))]
                    X$treat <- X$treat[!is.na(X$cluster) & sapply(X$cluster, function(x) isTRUE(all.equal(as.character(x), levels(X$cluster)[which.cluster])))]
                    var <- var[!is.na(X$cluster) & sapply(X$cluster, function(x) isTRUE(all.equal(as.character(x), levels(X$cluster)[which.cluster])))]
                    title <- paste0(title, "\nin cluster ", levels(X$cluster)[which.cluster])
                }
                else stop(paste0(which.cluster, " does not correspond to a given cluster."), call. = FALSE)
                
            }
            else stop("The argument to which.cluster must be a single value corresponding to the cluster for which distributions are to be displayed.", call. = FALSE)
        }
        else stop("Clusters were specified but no which.cluster value was supplied.", call. = FALSE)
    }
    else if (length(which.cluster) > 0) {
        warning("which.cluster was specified but no cluster values were supplied. Ignoring which.cluster.", call. = FALSE)
    }
    
    un.weights <- rep(1, length(X$treat))
    if (is.null(X$weights) || isTRUE(un)) X$weights <- un.weights
    
    if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) { #Continuous treatments
        treat <- X$treat
        weights <- X$weights
        if (length(unique(var)) <= 2 || is.factor(var) || is.character(var)) { #Categorical vars
            weights <- mapply(function(w, v) w / sum(weights[var==v]), w = weights, v = var)
            var <- factor(var)
            bp <- ggplot(mapping = aes(treat, fill = var, weight = weights)) + 
                geom_density(alpha = .4) + 
                labs(fill = var.name, y = "Density", x = "Treat", title = title)
        }
        else { #Continuous vars
            bp <- ggplot(mapping = aes(x = var, y = treat, weight = weights)) + 
                geom_point() + geom_smooth(se = FALSE, alpha = .1) + geom_smooth(method = "lm", se = FALSE, linetype = 2, alpha = .4) + 
                geom_hline(yintercept = w.m(treat, weights), linetype = 1, alpha = .9) + 
                labs(y = "Treat", x = var.name, title = title)
        }
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else { #Binary treatments
        X$treat <- binarize(X$treat)
        weights <- mapply(function(w, t) w / sum(X$weights[X$treat==t]), w = X$weights, t = X$treat)
        treat <- factor(X$treat)
        
        if (length(unique(var)) <= 2 || is.factor(var) || is.character(var)) { #Categorical vars
            var <- factor(var)
            bp <- ggplot(mapping = aes(var, fill = treat, weight = weights)) + 
                geom_bar(position = "dodge", alpha = .4, color = "black") + 
                labs(x = var.name, y = "Proportion", fill = "Treat", title = title) 
        }
        else { #Continuous vars
            bp <- ggplot(mapping = aes(var, fill = treat, weight = weights)) + 
                geom_density(alpha = .4) + 
                labs(x = var.name, y = "Density", fill = "Treat", title = title)
        }
    }
    
    return(bp)
}