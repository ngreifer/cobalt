bal.plot <- function(obj, var.name, ..., which, which.sub = NULL, cluster = NULL, which.cluster = NULL, imp = NULL, which.imp = NULL, size.weight = FALSE) {
    
    args <- list(...)
    
    if (missing(var.name)) {
        stop("An argument for var.name must be specified.", call. = FALSE)
    }
    
    X <- x2base(obj, ..., cluster = cluster, imp = imp, s.d.denom = "treated") #s.d.denom to avoid x2base warning
    
    if (var.name %in% names(X$covs)) X$var <- X$covs[, var.name]
    #else if (!is.null(args$data) && var.name %in% names(args$data)) var <- args$data[, var.name]
    else if (!is.null(X$addl) && var.name %in% names(X$addl)) X$var <- X$addl[, var.name]
    else if (!is.null(args$addl) && var.name %in% names(args$addl)) X$var <- args$addl[, var.name]
    else if (!is.null(X$distance) && var.name %in% names(X$distance)) X$var <- X$distance[, var.name]
    else if (!is.null(args$distance) && var.name %in% names(args$distance)) X$var <- args$distance[, var.name]
    else stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."))
    
    if (missing(which)) {
        if (length(args$un) > 0) {
            message("Note: \'un\' is deprecated; please use \'which\' for the same and added functionality.")
            if (args$un) which <- "unadjusted"
            else which <- "adjusted"
        }
        else {
            which <- "adjusted"
        }
    }
    else {
        which <- match.arg(tolower(which), c("adjusted", "unadjusted", "both"))
    }
    
    if (length(X$weights) > 0) {
        if(ncol(X$weights) > 1) stop("Only one set of weights may be specified at a time.")
        #else X$weights <- X$weights[,1]
        # if (ncol(X$weights) == 1) {
        #     names(weights) <- "Adjusted"
        #     facet <- "facet.which"
        # }
        # else {}
    }
    
    title <- paste0("Distributional Balance for \"", var.name, "\"")
    subtitle <- "Adjusted Sample"
    
    if (length(X$weights) == 0) {
        weight.names <- "Unadjusted"
        nweights <- 1
    }
    else {
        weight.names <- names(X$weights)
        nweights <- ncol(X$weights)
    }
    
    
    facet <- NULL
    if (length(X$subclass) > 0) {
        if (which %in% c("adjusted", "both")) {
            if (length(X$cluster) > 0) stop("Subclasses are not supported with clusters.", call. = FALSE)
            if (length(X$imp) > 0) stop("Subclasses are not supported with multiple imputations.", call. = FALSE)
            if (is.null(which.sub)) { #display all subs
                weights <- X$weights[!is.na(X$subclass), 1]
                treat <- X$treat[!is.na(X$subclass)]
                var <- X$var[!is.na(X$subclass)]
                subclass <- paste("Subclass", X$subclass[!is.na(X$subclass)])
                subtitle <- NULL
                #title <- paste0(title, "\nacross subclasses")
            }
            else {
                if (is.numeric(which.sub) && length(which.sub) == 1) {
                    if (which.sub %in% levels(X$subclass)) {
                        weights <- X$weights[!is.na(X$subclass) & X$subclass==which.sub, 1]
                        treat <- X$treat[!is.na(X$subclass) & X$subclass==which.sub]
                        var <- X$var[!is.na(X$subclass) & X$subclass==which.sub]
                        subclass <- paste("Subclass", X$subclass[!is.na(X$subclass) & X$subclass==which.sub])
                        #title <- paste0(title, "\nin subclass ", which.sub)
                        subtitle <- NULL
                        if (which == "unadjusted") {
                            message("which.sub specified and unadjusted sample requested; displaying adjusted sample.")
                            which == "adjusted"
                        }
                    }
                    else stop(paste0("\"", which.sub, "\" does not correspond to a subclass in the object."), call. = FALSE)
                    
                }
                else stop("The argument to which.sub must be a single number corresponding to the subclass for which distributions are to be displayed.", call. = FALSE)
            }
            if (which == "both") {
                weights <- c(weights, X$weights[, 1])
                treat <- c(treat, X$treat)
                var <- c(var, X$var)
                subclass <- factor(c(subclass, rep("Unadjusted Sample", length(X$treat))),
                                   levels = c("Unadjusted Sample", sort(unique(subclass))))
            }
            facet <- "subclass"
        }
    }  
    else if (length(X$subclass) == 0 && length(which.sub) > 0) {
        warning("which.sub was specified but no subclasses were supplied. Ignoring which.sub.", call. = FALSE)
    }
    else if (which == "unadjusted" && length(which.sub) > 0) {
        warning("which.sub was specified but the unadjusted sample was requested. Ignoring which.sub.", call. = FALSE)
    }
    
    #NULL: all; NA: none
    in.imp <- rep(TRUE, length(X$treat))
    #in.imp <- NULL
    if (length(X$imp) > 0) {
        if (length(which.imp) == 0 || all(is.na(which.imp))) {
            in.imp <- !is.na(X$imp)
        }
        else {
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
        facet <- c(facet, "imp")
    }
    else if (length(which.imp) > 0) {
        warning("which.imp was specified but no imp values were supplied. Ignoring which.imp.", call. = FALSE)
    }
    
    in.cluster <- rep(TRUE, length(X$treat))
    #in.cluster <- NULL
    if (length(X$cluster) > 0) {
        if (length(which.cluster) == 0 || all(is.na(which.cluster))) {
            in.cluster <- !is.na(X$cluster)
        }
        else {
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
        facet <- c(facet, "cluster")
    }
    else if (length(which.cluster) > 0) {
        warning("which.cluster was specified but no cluster values were supplied. Ignoring which.cluster.", call. = FALSE)
    }
    
    if (!"subclass" %in% facet) {
        nobs <- sum(in.imp & in.cluster)
        
        imp <- cluster <- character(nobs)
        treat <- X$treat[in.imp & in.cluster]
        var <- X$var[in.imp & in.cluster]
        
        if (is.null(X$weights) || which == "unadjusted") {
            weights <- rep(1, nobs)
            facet.which <- rep("Unadjusted Sample", nobs)
            subtitle <- "Unadjusted Sample"
        }
        else {
            weights <- X$weights[in.imp & in.cluster, 1]
            facet.which <- rep("Adjusted Sample", nobs)
            subtitle <- "Adjusted Sample"
        }
        
        if ("imp" %in% facet) imp <- paste("Imputation", X$imp[in.imp & in.cluster])
        if ("cluster" %in% facet) cluster <- factor(X$cluster[in.imp & in.cluster])
        
        if (which == "both") {
            facet.which <- factor(c(facet.which,
                                    rep("Unadjusted Sample", nobs)), 
                                  levels = c("Unadjusted Sample", "Adjusted Sample"))
            weights <- c(weights, rep(1, nobs))
            treat <- rep(treat, 2)
            var <- rep(var, 2)
            if ("imp" %in% facet) imp <- rep(imp, 2)
            if ("cluster" %in% facet) cluster <- rep(cluster, 2)
            facet <- c("facet.which", facet)
            subtitle <- NULL
        }
    }
    
    if (length(unique(treat)) > 2 && is.numeric(treat)) { #Continuous treatments
        is.categorical.var <- length(unique(var)) <= 2 || is.factor(var) || is.character(var)
        if (identical(facet,"subclass")) {
            if (is.categorical.var) {
                weights <- mapply(function(w, s, v) w[var==v] / sum(weights[subclass==s & var==v]), w = weights, s = subclass, v = var)
            }
            else {
                weights <- mapply(function(w, s) w / sum(weights[subclass==s]), w = weights, s = subclass)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, subclass = subclass)
        }
        else if (identical(facet,"cluster")) {
            if (is.categorical.var) {
                weights <- mapply(function(w, c, v) w / sum(weights[cluster==c & var==v]), w = weights, c = cluster, v = var)
            }
            else {
                weights <- mapply(function(w, c) w / sum(weights[cluster==c]), w = weights, c = cluster)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, cluster = cluster)
        }
        else if ((identical(facet,"imp"))) {
            if (is.categorical.var) {
                weights <- mapply(function(w, i, v) w / sum(weights[imp==i & var==v]), w = weights, i = imp, v = var)
            }
            else {
                weights <- mapply(function(w, i) w / sum(weights[imp==i]), w = weights, i = imp)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, imp = imp)
        }
        else if (identical(facet,c("imp", "cluster"))) {
            if (is.categorical.var) {
                weights <- mapply(function(w, i, c, v) w / sum(weights[imp==i & cluster == c & var==v]), w = weights, i = imp, c = cluster, v = var)
            }
            else {
                weights <- mapply(function(w, i, c) w / sum(weights[imp==i & cluster == c]), w = weights, i = imp, c = cluster)
            }
            d <- data.frame(weights = weights, treat = treat, var = var, imp = imp, cluster = cluster)
        }
        else if (identical(facet, "facet.which")) {
            if (is.categorical.var) {
                weights <- mapply(function(w, v, f.which) w / sum(weights[var==v & facet.which==f.which]), w = weights, v = var, f.which = facet.which)
            }
            else {
                #weights <- weights
            }
            d <- data.frame(weights = weights, treat = treat, var = var, facet.which = facet.which)
        }
        else {
            if (is.categorical.var) {
                weights <- mapply(function(w, v) w / sum(weights[var==v]), w = weights, v = var)
            }
            else {
                weights <- weights
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
            bp <- ggplot(d, mapping = aes(x = var, y = treat, weight = weights))
            if (which == "unadjusted" || !isTRUE(size.weight)) bp <- bp + geom_point()
            else bp <- bp + geom_point(aes(size = weights))
            bp <- bp + geom_smooth(method = "loess", se = FALSE, alpha = .1) + geom_smooth(method = "lm", se = FALSE, linetype = 2, alpha = .4) + 
                geom_hline(yintercept = w.m(treat, weights), linetype = 1, alpha = .9) + 
                labs(y = "Treat", x = var.name, title = title, subtitle = subtitle)
        }
    }
    else { #Categorical treatments (multinomial supported)
        treat <- factor(treat)
        
        if (length(facet) > 0) {
            if ("subclass" %in% facet) {
                weights <- mapply(function(w, t, s) w / sum(weights[treat==t & subclass==s]), w = weights, t = treat, s = subclass)
                d <- data.frame(weights = weights, treat = treat, var = var, subclass = subclass)
            }
            else {
                weights <- mapply(function(w, t, c, i, f.which) w / sum(weights[treat==t & cluster==c & imp==i & facet.which==f.which]), w = weights, t = treat, c= cluster, i = imp, f.which = facet.which)
                d <- data.frame(weights = weights, treat = treat, var = var, cluster = cluster, imp = imp, facet.which = facet.which)
                
            }
            
            # if (identical(facet,"subclass")) {
            #     weights <- mapply(function(w, t, s) w / sum(weights[treat==t & subclass==s]), w = weights, t = treat, s = subclass)
            #     d <- data.frame(weights = weights, treat = treat, var = var, subclass = subclass)
            # }
            # else if (identical(facet,"cluster")) {
            #     weights <- mapply(function(w, t, c) w / sum(weights[treat==t & cluster==c]), w = weights, t = treat, c = cluster)
            #     d <- data.frame(weights = weights, treat = treat, var = var, cluster = cluster)
            # }
            # else if (identical(facet,"imp")){
            #     weights <- mapply(function(w, t, i) w / sum(weights[treat==t & imp==i]), w = weights, t = treat, i = imp)
            #     d <- data.frame(weights = weights, treat = treat, var = var, imp = imp)
            # }
            # else if (identical(facet, c("imp", "cluster"))) {
            #     weights <- mapply(function(w, t, c, i) w / sum(weights[treat==t & cluster==c & imp==i]), w = weights, t = treat, c= cluster, i = imp)
            #     d <- data.frame(weights = weights, treat = treat, var = var, cluster = cluster, imp = imp)
            # }
            # else if (identical(facet, "facet.which")) {
            #     weights <- mapply(function(w, t, f.which) w / sum(weights[treat==t & facet.which==f.which]), w = weights, t = treat, f.which = facet.which)
            #     d <- data.frame(weights = weights, treat = treat, var = var, facet.which = facet.which)
            #     
            # }
        }
        else {
            weights <- mapply(function(w, t) w / sum(weights[treat==t]), w = weights, t = treat)
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
        if (length(facet) >= 2) {
            bp <- bp + facet_grid(f.build(facet[1], facet[-1]), drop = FALSE)
        }
        else bp <- bp + facet_grid(f.build(".", facet), drop = FALSE)
    }
    
    return(bp)
}