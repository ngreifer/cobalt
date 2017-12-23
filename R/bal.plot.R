bal.plot <- function(obj, var.name, ..., which, which.sub = NULL, cluster = NULL, which.cluster = NULL, imp = NULL, which.imp = NULL, which.treat = NULL, which.time = NULL, size.weight = FALSE) {
    
    args <- list(...)
    
    if (missing(var.name)) {
        stop("An argument for var.name must be specified.", call. = FALSE)
    }
    
    if (all(class(obj) == "list")) {
        if (all(sapply(obj, is.formula))) class(obj) <- "formula.list"
        else if (all(sapply(obj, is.data.frame))) class(obj) <- "data.frame.list"
        else stop("If obj is a list, it must be a list of formulas specifying the treatment/covariate relationships at each time point or a list of data frames containing covariates to be assessed at each time point.", call. = FALSE)

    }
    X <- x2base(obj, ..., cluster = cluster, imp = imp, s.d.denom = "treated") #s.d.denom to avoid x2base warning
    
    if (length(X$covs.list) > 0) {
        var.list <- vector("list", length(X$covs.list))
        appears.in.time <- rep(TRUE, length(X$covs.list))
        for (i in seq_along(X$covs.list)) {
            if (var.name %in% names(X$covs.list[[i]])) var.list[[i]] <- X$covs.list[[i]][[var.name]]
            else if (!is.null(X$addl.list) && var.name %in% names(X$addl.list[[i]])) var.list[[i]] <- X$addl[[var.name]]
            else if (!is.null(X$distance.list) && var.name %in% names(X$distance.list[[i]])) var.list[[i]] <- X$distance.list[[i]][[var.name]]
            else appears.in.time[i] <- FALSE
        }
        if (all(sapply(var.list, function(x) length(x) == 0))) stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."), call. = FALSE)
        X$var <- unlist(var.list[appears.in.time])
        X$time <- rep(seq_along(X$covs.list)[appears.in.time], each = nrow(X$covs.list[[1]]))
        X$treat <- unlist(X$treat.list[appears.in.time])
        if (length(names(X$treat.list)) > 0) treat.names <- names(X$treat.list)
        else treat.names <- seq_along(X$treat.list)
        if (length(X$weights) > 0) X$weights <- do.call("rbind", replicate(sum(appears.in.time), list(X$weights)))
    }
    else {
        if (var.name %in% names(X$covs)) X$var <- X$covs[[var.name]]
        else if (!is.null(X$addl) && var.name %in% names(X$addl)) X$var <- X$addl[[var.name]]
        else if (!is.null(args$addl) && var.name %in% names(args$addl)) X$var <- args$addl[[var.name]]
        else if (!is.null(X$distance) && var.name %in% names(X$distance)) X$var <- X$distance[[var.name]]
        else if (!is.null(args$distance) && var.name %in% names(args$distance)) X$var <- args$distance[[var.name]]
        else stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."), call. = FALSE)
    }
    
    #Density arguments supplied through ...
    if (length(args$bw) > 0) bw <- args$bw else bw <- "nrd0"
    if (length(args$adjust) > 0) adjust <- args$adjust else adjust <- 1
    if (length(args$kernel) > 0) kernel <- args$kernel else kernel <- "gaussian"
    if (length(args$n) > 0) n <- args$n else n <- 512

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
        if (length(X$weights) == 0 && length(X$subclass) == 0) which <- "unadjusted"
        else {
            which <- tryCatch(match.arg(tolower(which), c("adjusted", "unadjusted", "both")),
                              error = function(cond) stop(paste("The argument to 'which' should be one of", word.list(c("adjusted", "unadjusted", "both"), "or", quotes = TRUE)), call. = FALSE))
        }
    }
    
    title <- paste0("Distributional Balance for \"", var.name, "\"")
    subtitle <- NULL
    
    facet <- NULL
    is.categorical.var <- length(unique(X$var)) <= 2 || is.factor(X$var) || is.character(X$var)
    
    if (length(X$subclass) > 0) {
        if (which %in% c("adjusted", "both")) {
            if (length(X$cluster) > 0) stop("Subclasses are not supported with clusters.", call. = FALSE)
            if (length(X$imp) > 0) stop("Subclasses are not supported with multiple imputations.", call. = FALSE)
            if (is.null(which.sub)) { 
                which.sub <- levels(X$subclass)
            }
            if (is.numeric(which.sub) && any(which.sub %in% levels(X$subclass))) {
                if (any(!which.sub %in% levels(X$subclass))) {
                    w.l <- word.list(which.sub[!which.sub %in% levels(X$subclass)])
                    warning(paste(w.l, ifelse(attr(w.l, "plural"), "do", "does"), "not correspond to any subclass in the object and will be ignored."), call. = FALSE)
                    which.sub <- which.sub[which.sub %in% levels(X$subclass)]
                }
                in.sub <- !is.na(X$subclass) & X$subclass %in% which.sub
                D <- setNames(as.data.frame(matrix(0, nrow = sum(in.sub), ncol = 4)),
                              c("weights", "treat", "var", "subclass"))
                D$weights <- rep(1, nrow(D))
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
            D2$weights <- rep(1, nrow(D2))
            D2$treat <- X$treat
            D2$var <- X$var
            D2$subclass <- rep("Unadjusted Sample", length(X$treat))
            D <- rbind(D2, D, stringsAsFactors = TRUE)
            D$subclass <- relevel(factor(D$subclass), "Unadjusted Sample")
        }
        
    }
    else if (length(X$subclass) == 0 && length(which.sub) > 0) {
        warning("which.sub was specified but no subclasses were supplied. Ignoring which.sub.", call. = FALSE)
    }
    else if (which == "unadjusted" && length(which.sub) > 0) {
        warning("which.sub was specified but the unadjusted sample was requested. Ignoring which.sub.", call. = FALSE)
    }
    
    if (!"subclass" %in% facet) {
        
        facet <- "facet.which"
        
        if (length(X$weights) == 0 || which == "unadjusted") {
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
            facet <- c("imp", facet)
        }
        else if (length(which.imp) > 0) {
            warning("which.imp was specified but no imp values were supplied. Ignoring which.imp.", call. = FALSE)
        }
        
        in.cluster <- rep(TRUE, length(X$var))
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
            facet <- c("cluster", facet)
        }
        else if (length(which.cluster) > 0) {
            warning("which.cluster was specified but no cluster values were supplied. Ignoring which.cluster.", call. = FALSE)
        }
        
        in.time <- rep(TRUE, length(X$var))
        if (length(X$time) > 0) {
            if (length(which.time) == 0 || all(is.na(which.time))) {
                in.time <- !is.na(X$time)
            }
            else {
                if (is.numeric(which.time)) {
                    if (all(which.time %in% seq_along(X$covs.list))) {
                        if (all(which.time %in% seq_along(X$covs.list)[appears.in.time])) {
                            #nothing; which.time is good
                        }
                        else if (any(which.time %in% seq_along(X$covs.list)[appears.in.time])) {
                            warning(paste0(var.name, " does not appear in time period ", word.list(which.time[!which.time %in% seq_along(X$covs.list)[appears.in.time]], "or"), "."), call. = FALSE)
                            which.time <- which.time[which.time %in% seq_along(X$covs.list)[appears.in.time]]
                        }
                        else {
                            stop(paste0(var.name, " does not appear in time period ", word.list(which.time, "or"), "."), call. = FALSE)
                        }
                        in.time <- !is.na(X$time) & X$time %in% which.time
                    }
                    else {
                        stop(paste0("The following inputs to which.time do not correspond to given time periods:\n\t", word.list(which.time[!which.time %in% seq_along(X$covs.list)])), call. = FALSE)
                    }
                }
                else if (is.character(which.time)) {
                    if (all(which.time %in% treat.names)) {
                        if (all(which.time %in% treat.names[appears.in.time])) {
                            #nothing; which.time is good
                        }
                        else if (any(which.time %in% treat.names[appears.in.time])) {
                            time.periods <- word.list(which.time[!which.time %in% treat.names[appears.in.time]], "and")
                            warning(paste0(var.name, " does not appear in the time period", ifelse(attr(time.periods, "plural"), "s ", " "),
                                           "corresponding to treatment", ifelse(attr(time.periods, "plural"), "s ", " "),
                                           time.periods, "."), call. = FALSE)
                            which.time <- which.time[which.time %in% treat.names[appears.in.time]]
                        }
                        else {
                            time.periods <- word.list(which.time, "and")
                            stop(paste0(var.name, " does not appear in the time period", ifelse(attr(time.periods, "plural"), "s ", " "),
                                        "corresponding to treatment", ifelse(attr(time.periods, "plural"), "s ", " "),
                                        time.periods, "."), call. = FALSE)
                        }
                        in.time <- !is.na(X$time) & treat.names[X$time] %in% which.time
                        
                    }
                    else {
                        stop(paste0("The following inputs to which.time do not correspond to given time periods:\n\t", word.list(which.time[!which.time %in% treat.names])), call. = FALSE)
                    }
                }
                else stop("The argument to which.time must be the names or indices corresponding to the time periods for which distributions are to be displayed.", call. = FALSE)
            }
            facet <- c("time", facet)
        }
        else if (length(which.time) > 0) {
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
    
    if (nunique(D$treat) > 2 && is.numeric(D$treat)) { #Continuous treatments
        if ("subclass" %in% facet) {
            if (is.categorical.var) {
                weights <- with(D, sapply(seq_along(weights), function(i) weights[i] / sum(weights[subclass==subclass[i] & var==var[i]])))
            }
            else {
                weights <- with(D, sapply(seq_along(weights), function(i) weights[i] / sum(weights[subclass==subclass[i]])))
            }
            
            d <- data.frame(weights = weights, treat = D$treat, var = D$var, subclass = D$subclass)
        }
        else {
            if (is.categorical.var) {
                weights <- with(D, sapply(seq_along(weights), function(i) weights[i] / sum(weights[cluster==cluster[i] & imp==imp[i] & time==time[i] & facet.which==facet.which[i] & var==var[i]])))
            }
            else {
                weights <- with(D, sapply(seq_along(weights), function(i) weights[i] / sum(weights[cluster==cluster[i] & imp==imp[i] & time==time[i] & facet.which==facet.which[i]])))
            }
            d <- data.frame(weights = weights, treat = D$treat, var = D$var, cluster = D$cluster, imp = D$imp, time = D$time, facet.which = D$facet.which)
            
        }

        if (is.categorical.var) { #Categorical vars
            d$var <- factor(d$var)
            bp <- ggplot(d, mapping = aes(treat, fill = var, weight = weights)) + 
                geom_density(alpha = .4, bw = bw, adjust = adjust, kernel = kernel, n = n) + 
                labs(fill = var.name, y = "Density", x = "Treat", title = title, subtitle = subtitle)
        }
        else { #Continuous vars
            bp <- ggplot(d, mapping = aes(x = var, y = treat, weight = weights))
            if (which == "unadjusted" || !isTRUE(size.weight)) bp <- bp + geom_point()
            else bp <- bp + geom_point(aes(size = weights))
            bp <- bp + geom_smooth(method = "loess", se = FALSE, alpha = .1) + geom_smooth(method = "lm", se = FALSE, linetype = 2, alpha = .4) + 
                geom_hline(yintercept = w.m(d$treat, d$weights), linetype = 1, alpha = .9) + 
                labs(y = "Treat", x = var.name, title = title, subtitle = subtitle)
        }
    }
    else { #Categorical treatments (multinomial supported)
        D$treat <- factor(D$treat)
        
        if (length(which.treat) == 0) 
            which.treat <- character(0)
        else if (is.numeric(which.treat)) {
            which.treat <- levels(D$treat)[seq_along(levels(D$treat)) %in% which.treat]
            if (length(which.treat) == 0) {
                warning("No numbers in which.treat correspond to treatment values. All treatment groups will be displayed.", call. = FALSE)
                which.treat <- character(0)
            }
        }
        else if (is.character(which.treat)) {
            which.treat <- levels(D$treat)[levels(D$treat) %in% which.treat]
            if (length(which.treat) == 0) {
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
        if (length(which.treat) > 0 && all(!is.na(which.treat))) D <- D[D$treat %in% which.treat,]
        
        for (i in names(D)[sapply(D, is.factor)]) D[[i]] <- factor(D[[i]])
        
        if (length(facet) > 0) {
            if ("subclass" %in% facet) {
                weights <- with(D, sapply(seq_along(weights), function(i) weights[i] / sum(weights[treat==treat[i] & subclass==subclass[i]])))
                d <- data.frame(weights = weights, treat = D$treat, var = D$var, subclass = D$subclass)
            }
            else {
                weights <- with(D, sapply(seq_along(weights), function(i) weights[i] / sum(weights[treat==treat[i] & cluster==cluster[i] & imp==imp[i] & time==time[i] & facet.which==facet.which[i]])))
                d <- data.frame(weights = weights, treat = D$treat, var = D$var, cluster = D$cluster, imp = D$imp, time = D$time, facet.which = D$facet.which)
            }
        }
        else {
            weights <- with(D, sapply(seq_along(weights), function(i) weights[i] / sum(weights[treat==treat[i]])))
            d <- data.frame(weights = weights, treat = D$treat, var = D$var)
        }
        
        if (length(unique(d$var)) <= 2 || is.factor(d$var) || is.character(d$var)) { #Categorical vars
            d$var <- factor(d$var)
            bp <- ggplot(d, mapping = aes(var, fill = treat, weight = weights)) + 
                geom_bar(position = "dodge", alpha = .4, color = "black") + 
                labs(x = var.name, y = "Proportion", fill = "Treat", title = title, subtitle = subtitle) + 
                scale_x_discrete(drop=FALSE) + scale_fill_discrete(drop=FALSE)
        }
        else { #Continuous vars
            bp <- ggplot(d, mapping = aes(var, fill = treat, weight = weights)) + 
                geom_density(alpha = .4, bw = bw, adjust = adjust, kernel = kernel, n = n) + 
                labs(x = var.name, y = "Density", fill = "Treat", title = title, subtitle = subtitle)
        }
    }
    
    if (length(facet) > 0) {
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