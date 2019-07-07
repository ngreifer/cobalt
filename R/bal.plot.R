bal.plot <- function(obj, var.name, ..., which, which.sub = NULL, cluster = NULL, which.cluster = NULL, imp = NULL, which.imp = NULL, which.treat = NULL, which.time = NULL, size.weight = FALSE, mirror = FALSE, type = c("density", "histogram"), colors = NULL) {
    
    tryCatch(identity(obj), error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval(.call))
    }
    
    args <- list(...)
    
    obj <- is.designmatch(obj)
    obj <- is.time.list(obj)

    X <- x2base(obj, ..., cluster = cluster, imp = imp, s.d.denom = "treated") #s.d.denom to avoid x2base warning
    
    if (is_not_null(X$covs.list)) {
        if (missing(var.name)) {
            var.name <- names(X$covs.list[[1]])[1]
            message(paste0("No var.name was provided. Dispalying balance for ", var.name, "."))
        }
        var.list <- vector("list", length(X$covs.list))
        appears.in.time <- rep(TRUE, length(X$covs.list))
        for (i in seq_along(X$covs.list)) {
            if (var.name %in% names(X$covs.list[[i]])) var.list[[i]] <- X$covs.list[[i]][[var.name]]
            else if (!is.null(X$addl.list) && var.name %in% names(X$addl.list[[i]])) var.list[[i]] <- X$addl[[var.name]]
            else if (!is.null(X$distance.list) && var.name %in% names(X$distance.list[[i]])) var.list[[i]] <- X$distance.list[[i]][[var.name]]
            else appears.in.time[i] <- FALSE
        }
        if (all(sapply(var.list, is_null))) stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."), call. = FALSE)
        X$var <- unlist(var.list[appears.in.time])
        X$time <- rep(seq_along(X$covs.list)[appears.in.time], each = NROW(X$covs.list[[1]]))
        X$treat <- unlist(X$treat.list[appears.in.time])
        if (is_not_null(names(X$treat.list))) treat.names <- names(X$treat.list)
        else treat.names <- seq_along(X$treat.list)
        if (is_not_null(X$weights)) X$weights <- do.call("rbind", replicate(sum(appears.in.time), list(X$weights)))
    }
    else {
        if (missing(var.name)) {
            var.name <- names(X$covs)[1]
            message(paste0("No var.name was provided. Dispalying balance for ", var.name, "."))
        }
        if (var.name %in% names(X$covs)) X$var <- X$covs[[var.name]]
        else if (!is.null(X$addl) && var.name %in% names(X$addl)) X$var <- X$addl[[var.name]]
        else if (!is.null(args$addl) && var.name %in% names(args$addl)) X$var <- args$addl[[var.name]]
        else if (!is.null(X$distance) && var.name %in% names(X$distance)) X$var <- X$distance[[var.name]]
        else if (!is.null(args$distance) && var.name %in% names(args$distance)) X$var <- args$distance[[var.name]]
        else stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."), call. = FALSE)
    }
    
    #Density arguments supplied through ...
    if (is_not_null(args$bw)) bw <- args$bw else bw <- "nrd0"
    if (is_not_null(args$adjust)) adjust <- args$adjust else adjust <- 1
    if (is_not_null(args$kernel)) kernel <- args$kernel else kernel <- "gaussian"
    if (is_not_null(args$n)) n <- args$n else n <- 512
    
    if (missing(which)) {
        if (is_not_null(args$un)) {
            message("Note: \'un\' is deprecated; please use \'which\' for the same and added functionality.")
            if (args$un) which <- "unadjusted"
            else which <- "adjusted"
        }
        else {
            which <- "adjusted"
        }
    }
    else {
        if (is_null(X$weights) && is_null(X$subclass)) which <- "unadjusted"
        else {
            which <- match_arg(tolower(which), c("adjusted", "unadjusted", "both"))
        }
    }
    
    title <- paste0("Distributional Balance for \"", var.name, "\"")
    subtitle <- NULL
    
    facet <- NULL
    is.categorical.var <- is_binary(X$var) || is.factor(X$var) || is.character(X$var)
    
    if (is_not_null(X$subclass)) {
        if (which %in% c("adjusted", "both")) {
            if (is_not_null(X$cluster)) stop("Subclasses are not supported with clusters.", call. = FALSE)
            if (is_not_null(X$imp)) stop("Subclasses are not supported with multiple imputations.", call. = FALSE)
            if (is_null(which.sub)) { 
                which.sub <- levels(X$subclass)
            }
            if (is.numeric(which.sub) && any(which.sub %in% levels(X$subclass))) {
                if (any(!which.sub %in% levels(X$subclass))) {
                    w.l <- word_list(which.sub[!which.sub %in% levels(X$subclass)])
                    warning(paste(w.l, ifelse(attr(w.l, "plural"), "do", "does"), "not correspond to any subclass in the object and will be ignored."), call. = FALSE)
                    which.sub <- which.sub[which.sub %in% levels(X$subclass)]
                }
                in.sub <- !is.na(X$subclass) & X$subclass %in% which.sub
                D <- setNames(as.data.frame(matrix(0, nrow = sum(in.sub), ncol = 4)),
                              c("weights", "treat", "var", "subclass"))
                D$weights <- rep(1, NROW(D))
                D$treat <- X$treat[in.sub]
                D$var <- X$var[in.sub]
                D$subclass <- paste("Subclass", X$subclass[in.sub])
                #title <- paste0(title, "\nin subclass ", which.sub)
                subtitle <- NULL
            }
            else stop("The argument to which.sub must be a number vector corresponding to the subclass for which distributions are to be displayed.", call. = FALSE)
            facet <- "subclass"
        }
        #D$weights <- rep(1, length(treat))
        if (which == "both") {
            D2 <- setNames(as.data.frame(matrix(0, nrow = length(X$treat), ncol = 4)),
                           c("weights", "treat", "var", "subclass"))
            D2$weights <- rep(1, NROW(D2))
            D2$treat <- X$treat
            D2$var <- X$var
            D2$subclass <- rep("Unadjusted Sample", length(X$treat))
            D <- rbind(D2, D, stringsAsFactors = TRUE)
            D$subclass <- relevel(factor(D$subclass), "Unadjusted Sample")
        }
        
    }
    else if (is_null(X$subclass) && is_not_null(which.sub)) {
        warning("which.sub was specified but no subclasses were supplied. Ignoring which.sub.", call. = FALSE)
    }
    else if (which == "unadjusted" && is_not_null(which.sub)) {
        warning("which.sub was specified but the unadjusted sample was requested. Ignoring which.sub.", call. = FALSE)
    }
    
    if ("subclass" %nin% facet) {
        
        facet <- "facet.which"
        
        if (is_null(X$weights) || which == "unadjusted") {
            X$weights <- data.frame(rep(1, length(X$treat)))
            names(X$weights) <- "Unadjusted Sample"
            #nweights <- 1
        }
        else {
            if (ncol(X$weights) == 1) {
                names(X$weights) <- "Adjusted Sample"
            }
            if (which == "both") {
                X$weights <- setNames(cbind(rep(1, length(X$treat)), X$weights),
                                      c("Unadjusted Sample", names(X$weights)))
            }
        }
        
        nweights <- ncol(X$weights)
        weight.names <- names(X$weights)
        
        #NULL: all; NA: none
        in.imp <- rep(TRUE, length(X$var))
        if (is_not_null(X$imp)) {
            if (is_null(which.imp) || all(is.na(which.imp))) {
                in.imp <- !is.na(X$imp)
            }
            else {
                if (is.numeric(which.imp)) {
                    if (all(which.imp %in% seq_len(nlevels(X$imp)))) {
                        in.imp <- !is.na(X$imp) & sapply(X$imp, function(x) !is.na(match(x, levels(X$imp)[which.imp])))
                    }
                    else {
                        stop(paste0("The following inputs to which.imp do not correspond to given imputations:\n\t", word_list(which.imp[!which.imp %in% seq_len(nlevels(X$imp))])), call. = FALSE)
                    }
                }
                else stop("The argument to which.imp must be the indices corresponding to the imputations for which distributions are to be displayed.", call. = FALSE)
            }
            facet <- c("imp", facet)
        }
        else if (is_not_null(which.imp)) {
            warning("which.imp was specified but no imp values were supplied. Ignoring which.imp.", call. = FALSE)
        }
        
        in.cluster <- rep(TRUE, length(X$var))
        if (is_not_null(X$cluster)) {
            if (is_null(which.cluster)|| all(is.na(which.cluster))) {
                in.cluster <- !is.na(X$cluster)
            }
            else {
                if (is.numeric(which.cluster)) {
                    if (all(which.cluster %in% seq_len(nlevels(X$cluster)))) {
                        in.cluster <- !is.na(X$cluster) & sapply(X$cluster, function(x) !is.na(match(x, levels(X$cluster)[which.cluster])))
                    }
                    else {
                        stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word_list(which.cluster[!which.cluster %in% seq_len(nlevels(X$cluster))])), call. = FALSE)
                    }
                }
                else if (is.character(which.cluster)) {
                    if (all(!is.na(match(which.cluster, levels(X$cluster))))) {
                        in.cluster <- !is.na(X$cluster) & sapply(X$cluster, function(x) !is.na(match(x, which.cluster)))
                    }
                    else {
                        stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word_list(which.cluster[is.na(match(which.cluster, levels(X$cluster)))])), call. = FALSE)
                    }
                }
                else stop("The argument to which.cluster must be the names or indices corresponding to the clusters for which distributions are to be displayed.", call. = FALSE)
            }
            facet <- c("cluster", facet)
        }
        else if (is_not_null(which.cluster)) {
            warning("which.cluster was specified but no cluster values were supplied. Ignoring which.cluster.", call. = FALSE)
        }
        
        in.time <- rep(TRUE, length(X$var))
        if (is_not_null(X$time)) {
            if (is_null(which.time) || all(is.na(which.time))) {
                in.time <- !is.na(X$time)
            }
            else {
                if (is.numeric(which.time)) {
                    if (all(which.time %in% seq_along(X$covs.list))) {
                        if (all(which.time %in% seq_along(X$covs.list)[appears.in.time])) {
                            #nothing; which.time is good
                        }
                        else if (any(which.time %in% seq_along(X$covs.list)[appears.in.time])) {
                            warning(paste0(var.name, " does not appear in time period ", word_list(which.time[!which.time %in% seq_along(X$covs.list)[appears.in.time]], "or"), "."), call. = FALSE)
                            which.time <- which.time[which.time %in% seq_along(X$covs.list)[appears.in.time]]
                        }
                        else {
                            stop(paste0(var.name, " does not appear in time period ", word_list(which.time, "or"), "."), call. = FALSE)
                        }
                        in.time <- !is.na(X$time) & X$time %in% which.time
                    }
                    else {
                        stop(paste0("The following inputs to which.time do not correspond to given time periods:\n\t", word_list(which.time[!which.time %in% seq_along(X$covs.list)])), call. = FALSE)
                    }
                }
                else if (is.character(which.time)) {
                    if (all(which.time %in% treat.names)) {
                        if (all(which.time %in% treat.names[appears.in.time])) {
                            #nothing; which.time is good
                        }
                        else if (any(which.time %in% treat.names[appears.in.time])) {
                            time.periods <- word_list(which.time[!which.time %in% treat.names[appears.in.time]], "and")
                            warning(paste0(var.name, " does not appear in the time period", ifelse(attr(time.periods, "plural"), "s ", " "),
                                           "corresponding to treatment", ifelse(attr(time.periods, "plural"), "s ", " "),
                                           time.periods, "."), call. = FALSE)
                            which.time <- which.time[which.time %in% treat.names[appears.in.time]]
                        }
                        else {
                            time.periods <- word_list(which.time, "and")
                            stop(paste0(var.name, " does not appear in the time period", ifelse(attr(time.periods, "plural"), "s ", " "),
                                        "corresponding to treatment", ifelse(attr(time.periods, "plural"), "s ", " "),
                                        time.periods, "."), call. = FALSE)
                        }
                        in.time <- !is.na(X$time) & treat.names[X$time] %in% which.time
                        
                    }
                    else {
                        stop(paste0("The following inputs to which.time do not correspond to given time periods:\n\t", word_list(which.time[!which.time %in% treat.names])), call. = FALSE)
                    }
                }
                else stop("The argument to which.time must be the names or indices corresponding to the time periods for which distributions are to be displayed.", call. = FALSE)
            }
            facet <- c("time", facet)
        }
        else if (is_not_null(which.time)) {
            warning("which.time was specified but a point treatment was supplied. Ignoring which.time.", call. = FALSE)
        }
        
        nobs <- sum(in.imp & in.cluster & in.time)
        if (nobs == 0) stop("No observations to display.", call. = FALSE)
        
        Q <- setNames(vector("list", nweights), weight.names)
        for (i in weight.names) {
            Q[[i]] <- setNames(as.data.frame(matrix(0, ncol = 7, nrow = nobs)),
                               c("imp", "cluster", "time", "treat", "var", "weights", "facet.which"))
            Q[[i]]$imp <- Q[[i]]$cluster <- Q[[i]]$time <- character(nobs)
            Q[[i]]$treat <- X$treat[in.imp & in.cluster & in.time]
            Q[[i]]$var <- X$var[in.imp & in.cluster & in.time]
            Q[[i]]$weights <-  X$weights[in.imp & in.cluster & in.time, i]
            Q[[i]]$facet.which <- rep(i, nobs)
            
            if ("imp" %in% facet) Q[[i]]$imp <- paste("Imputation", X$imp[in.imp & in.cluster & in.time])
            if ("cluster" %in% facet) Q[[i]]$cluster <- factor(X$cluster[in.imp & in.cluster & in.time])
            if ("time" %in% facet) Q[[i]]$time <- paste("Time", X$time[in.imp & in.cluster & in.time])
        }
        D <- do.call("rbind", Q)
        D$facet.which <- factor(D$facet.which, levels = c(weight.names[weight.names == "Unadjusted Sample"],
                                                          weight.names[weight.names != "Unadjusted Sample"]))
        
    }
    
    if ("facet.which" %in% facet) {
        if (nlevels(D$facet.which) == 1) {
            subtitle <- levels(D$facet.which)[1]
            facet <- facet[!facet %in% "facet.which"]
        }
    }
    
    if (!is_binary(D$treat) && is.numeric(D$treat)) { #Continuous treatments
        if ("subclass" %in% facet) {
            if (is.categorical.var) {
                weights <- with(D, ave(weights, subclass, var, FUN = function(x) x/sum(x)))
            }
            else {
                weights <- with(D, ave(weights, subclass, FUN = function(x) x/sum(x)))
            }
            d <- data.frame(weights = weights, treat = D$treat, var = D$var, subclass = D$subclass)
        }
        else {
            if (is.categorical.var) {
                weights <- with(D, ave(weights, cluster, imp, time, facet.which, var, FUN = function(x) x/sum(x)))
            }
            else {
                weights <- with(D, ave(weights, cluster, imp, time, facet.which, FUN = function(x) x/sum(x)))
            }
            d <- data.frame(weights = weights, treat = D$treat, var = D$var, cluster = D$cluster, imp = D$imp, time = D$time, facet.which = D$facet.which)
            
        }
        
        if (is.categorical.var) { #Categorical vars
            d$var <- factor(d$var)
            cat.sizes <- tapply(rep(1, NROW(d)), d$var, sum)
            smallest.cat <- names(cat.sizes)[which.min(cat.sizes)]
            if (is.character(bw)) {
                if (is.function(get0(paste0("bw.", bw)))) {
                    bw <- get0(paste0("bw.", bw))(d$treat[d$var == smallest.cat])
                }
                else {
                    stop(paste(bw, "is not an acceptable entry to bw. See ?stats::density for allowable options."), call. = FALSE)
                }
            }
            
            #Color
            ntypes <- length(cat.sizes)
            if (is_not_null(args$colours)) colors <- args$colours
            
            if (is_null(colors)) {
                colors <- gg_color_hue(ntypes)
            }
            else {
                if (length(colors) > ntypes) {
                    colors <- colors[seq_len(ntypes)]
                    warning(paste("Only using first", ntypes, "values in colors."), call. = FALSE)
                }
                else if (length(colors) < ntypes) {
                    warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
                    colors <- gg_color_hue(ntypes)
                }
                
                if (!all(sapply(colors, isColor))) {
                    warning("The argument to colors contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
                    colors <- gg_color_hue(ntypes)
                }
                
            }
            
            bp <- ggplot(d, mapping = aes(treat, fill = var, weight = weights)) + 
                geom_density(alpha = .4, bw = bw, adjust = adjust, kernel = kernel, n = n, trim = TRUE) + 
                labs(fill = var.name, y = "Density", x = "Treat", title = title, subtitle = subtitle) +
                scale_fill_manual(values = colors) + geom_hline(yintercept = 0)
        }
        else { #Continuous vars
            bp <- ggplot(d, mapping = aes(x = var, y = treat, weight = weights))
            if (which == "unadjusted" || !isTRUE(size.weight)) bp <- bp + geom_point(alpha = .9)
            else bp <- bp + geom_point(aes(size = weights), alpha = .9)
            bp <- bp + geom_smooth(method = "loess", se = FALSE, alpha = .1) + 
                geom_smooth(method = "lm", se = FALSE, linetype = 2, alpha = .4) + 
                geom_hline(yintercept = w.m(d$treat, d$weights), linetype = 1, alpha = .9) + 
                labs(y = "Treat", x = var.name, title = title, subtitle = subtitle)
        }
    }
    else { #Categorical treatments (multinomial supported)
        D$treat <- factor(D$treat)
        
        if (is_null(which.treat)) 
            which.treat <- character(0)
        else if (is.numeric(which.treat)) {
            which.treat <- levels(D$treat)[seq_along(levels(D$treat)) %in% which.treat]
            if (is_null(which.treat)) {
                warning("No numbers in which.treat correspond to treatment values. All treatment groups will be displayed.", call. = FALSE)
                which.treat <- character(0)
            }
        }
        else if (is.character(which.treat)) {
            which.treat <- levels(D$treat)[levels(D$treat) %in% which.treat]
            if (is_null(which.treat)) {
                warning("No names in which.treat correspond to treatment values. All treatment groups will be displayed.", call. = FALSE)
                which.treat <- character(0)
            }
        }
        else if (is.na(which.treat)) {
            which.treat <- character(0)
        }
        else {
            warning("The argument to which.treat must be NA, NULL, or a vector of treatment names or indices. All treatment groups will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
        if (is_not_null(which.treat) && all(!is.na(which.treat))) D <- D[D$treat %in% which.treat,]
        
        for (i in names(D)[sapply(D, is.factor)]) D[[i]] <- factor(D[[i]])
        
        if (is_not_null(facet)) {
            if ("subclass" %in% facet) {
                weights <- with(D, ave(weights, treat, subclass, FUN = function(x) x/sum(x)))
                d <- data.frame(weights = weights, treat = D$treat, var = D$var, subclass = D$subclass)
            }
            else {
                weights <- with(D, ave(weights, treat, cluster, imp, time, facet.which, FUN = function(x) x/sum(x)))
                d <- data.frame(weights = weights, treat = D$treat, var = D$var, cluster = D$cluster, imp = D$imp, time = D$time, facet.which = D$facet.which)
            }
        }
        else {
            weights <- with(D, ave(weights, treat, FUN = function(x) x/sum(x)))
            d <- data.frame(weights = weights, treat = D$treat, var = D$var)
        }
        
        #Color
        ntypes <- nlevels.treat <- nlevels(d$treat)
        if (is_not_null(args$colours)) colors <- args$colours
        
        if (is_null(colors)) {
            colors <- gg_color_hue(ntypes)
        }
        else {
            if (length(colors) > ntypes) {
                colors <- colors[seq_len(ntypes)]
                warning(paste("Only using first", ntypes, "values in colors."), call. = FALSE)
            }
            else if (length(colors) < ntypes) {
                warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
                colors <- gg_color_hue(ntypes)
            }
            
            if (!all(sapply(colors, isColor))) {
                warning("The argument to colors contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
                colors <- gg_color_hue(ntypes)
            }
        }
        names(colors) <- levels(d$treat)
        
        if (is_binary(d$var) || is.factor(d$var) || is.character(d$var)) { #Categorical vars
            d$var <- factor(d$var)
            bp <- ggplot(d, mapping = aes(x = var, fill = treat, weight = weights)) + 
                geom_bar(position = "dodge", alpha = .8, color = "black") + 
                labs(x = var.name, y = "Proportion", fill = "Treat", title = title, subtitle = subtitle) + 
                scale_x_discrete(drop=FALSE) + scale_fill_manual(drop=FALSE, values = colors) +
                geom_hline(yintercept = 0)
        }
        else { #Continuous vars
            t.sizes <- tapply(rep(1, NROW(d)), d$treat, sum)
            smallest.t <- names(t.sizes)[which.min(t.sizes)]
            if (is.character(bw)) {
                if (is.function(get0(paste0("bw.", bw)))) {
                    bw <- get0(paste0("bw.", bw))(d$var[d$treat == smallest.t])
                }
                else {
                    stop(paste(bw, "is not an acceptable entry to bw. See ?stats::density for allowable options."), call. = FALSE)
                }
            }

            # bp <- ggplot(d, mapping = aes(var, fill = treat, weight = weights)) +
            #     geom_density(alpha = .4, bw = bw, adjust = adjust, kernel = kernel, n = n, trim = TRUE) +
            #     labs(x = var.name, y = "Density", fill = "Treat", title = title, subtitle = subtitle)
            
            #colors <- setNames(gg_color_hue(nlevels.treat), levels(d$treat))
            
            if (nlevels.treat <= 2 && mirror == TRUE) {
                posneg <- c(1, -1)
                alpha <- .8
            }
            else {
                posneg <- rep(1, nlevels.treat)
                alpha <- .4
            }
            type <- match_arg(type)
            if (type == "histogram" && nlevels.treat <= 2) {
                if (is_null(args$bins) || !is.numeric(args$bins)) args$bins <- 12
                geom_fun <- function(t) geom_histogram(data = d[d$treat == levels(d$treat)[t],],
                                                       mapping = aes(x = var, y = posneg[t]*stat(count), weight = weights,
                                                                     fill = names(colors)[t]),
                                                       alpha = alpha, bins = args$bins, color = "black")
                ylab <- "Proportion"
                legend.alpha <- alpha
            }
            else {
                geom_fun <- function(t) geom_density(data = d[d$treat == levels(d$treat)[t],],
                                                     mapping = aes(x = var, y = posneg[t]*stat(density), weight = weights,
                                                                   fill = names(colors)[t]),
                                                     alpha = alpha, bw = bw, adjust = adjust,
                                                     kernel = kernel, n = n, trim = TRUE)
                ylab <- "Density"
                legend.alpha <- alpha/nlevels.treat
            }
            
            bp <- Reduce("+", c(list(ggplot(),
                                     labs(x = var.name, y = ylab, fill = "Treat", title = title, subtitle = subtitle)),
                                lapply(seq_len(nlevels.treat), geom_fun))) +
                scale_fill_manual(values = colors, guide = guide_legend(override.aes = list(alpha = legend.alpha, size = .25))) +
                geom_hline(yintercept = 0)
        }
    }
    
    if (is_not_null(facet)) {
        if (length(facet) >= 2) {
            if ("facet.which" %in% facet) {
                facet.formula <- f.build("facet.which", facet[!facet %in% "facet.which"])
            }
            else if ("imp" %in% facet) {
                facet.formula <- f.build("imp", facet[!facet %in% "imp"])
            }
            else {
                facets <- data.frame(facet = facet, length = sapply(facet, function(x) nlevels(d[[x]])),
                                     stringsAsFactors = FALSE)
                facets <- facets[with(facets, order(length, facet, decreasing = c(FALSE, TRUE))), ]
                facet.formula <- formula(facets)
            }
        }
        else facet.formula <- f.build(".", facet)
        bp <- bp + facet_grid(facet.formula, drop = FALSE, scales = ifelse("subclass" %in% facet, "free_x", "fixed"))
    }
    return(bp)
}