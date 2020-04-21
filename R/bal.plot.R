bal.plot <- function(obj, var.name, ..., which, which.sub = NULL, cluster = NULL, which.cluster = NULL, 
                     imp = NULL, which.imp = NULL, which.treat = NULL, which.time = NULL, 
                     mirror = FALSE, type = "density", colors = NULL, grid = FALSE, sample.names,
                     position = "right", facet.formula = NULL, disp.means = getOption("cobalt_disp.means", FALSE), 
                     alpha.weight = TRUE) {
    
    tryCatch(identity(obj), error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    .alls <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.all)), logical(1L))
    .nones <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.none)), logical(1L))
    if (any(c(.alls, .nones))) {
        .call[.alls] <- expression(NULL)
        .call[.nones] <- expression(NA)
        return(eval.parent(.call))
    }
    
    args <- list(...)
    
    obj <- process_designmatch(obj)
    obj <- process_time.list(obj)
    obj <- process_cem.match.list(obj)
    
    X <- x2base(obj, ..., cluster = cluster, imp = imp, s.d.denom = "all") #s.d.denom to avoid x2base warning
    
    if (is_not_null(X$covs.list)) {
        #Longitudinal
        if (missing(var.name)) {
            var.name <- names(X$covs.list[[1]])[1]
            message(paste0("No var.name was provided. Dispalying balance for ", var.name, "."))
        }
        var.list <- make_list(length(X$covs.list))
        appears.in.time <- rep.int(TRUE, length(X$covs.list))
        for (i in seq_along(X$covs.list)) {
            if (var.name %in% names(X$covs.list[[i]])) var.list[[i]] <- X$covs.list[[i]][[var.name]]
            else if (is_not_null(X$addl.list) && var.name %in% names(X$addl.list[[i]])) var.list[[i]] <- X$addl[[var.name]]
            else if (is_not_null(X$distance.list) && var.name %in% names(X$distance.list[[i]])) var.list[[i]] <- X$distance.list[[i]][[var.name]]
            else appears.in.time[i] <- FALSE
        }
        if (all(sapply(var.list, is_null))) stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."), call. = FALSE)
        X$var <- unlist(var.list[appears.in.time])
        X$time <- rep(seq_along(X$covs.list)[appears.in.time], each = NROW(X$covs.list[[1]]))
        X$treat.list <- lapply(X$treat.list, function(t) if (get.treat.type(t) != "continuous") treat_vals(t)[t] else t)
        X$treat <- unlist(X$treat.list[appears.in.time])
        if (is_not_null(names(X$treat.list))) treat.names <- names(X$treat.list)
        else treat.names <- seq_along(X$treat.list)
        if (is_not_null(X$weights)) X$weights <- X$weights[rep(seq_len(nrow(X$weights)), sum(appears.in.time)), , drop = FALSE]
        if (is_not_null(X$cluster)) X$cluster <- X$cluster[rep(seq_along(X$cluster), sum(appears.in.time))]
        if (is_not_null(X$imp)) X$imp <- X$imp[rep(seq_along(X$imp), sum(appears.in.time))]
    }
    else {
        #Point treatment
        if (missing(var.name)) {
            var.name <- names(X$covs)[1]
            message(paste0("No var.name was provided. Dispalying balance for ", var.name, "."))
        }
        if (var.name %in% names(X$covs)) X$var <- X$covs[[var.name]]
        else if (is_not_null(X$addl) && var.name %in% names(X$addl)) X$var <- X$addl[[var.name]]
        else if (is_not_null(args$addl) && var.name %in% names(args$addl)) X$var <- args$addl[[var.name]]
        else if (is_not_null(X$distance) && var.name %in% names(X$distance)) X$var <- X$distance[[var.name]]
        else if (is_not_null(args$distance) && var.name %in% names(args$distance)) X$var <- args$distance[[var.name]]
        else stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."), call. = FALSE)
        
        if (get.treat.type(X$treat) != "continuous") X$treat <- treat_vals(X$treat)[X$treat]
    }
    
    #Density arguments supplied through ...
    if (is_not_null(args$bw)) bw <- args$bw else bw <- "nrd0"
    if (is_not_null(args$adjust)) adjust <- args$adjust else adjust <- 1
    if (is_not_null(args$kernel)) kernel <- args$kernel else kernel <- "gaussian"
    if (is_not_null(args$n)) n <- args$n else n <- 512
    
    if (is_null(X$subclass)) {
        if (NCOL(X$weights) == 1L) weight.names <- "adjusted"
        else weight.names <- names(X$weights)
    }
    else weight.names <- "adjusted"
    
    if (missing(which)) {
        if (is_not_null(args$un)) {
            message("Note: \'un\' is deprecated; please use \'which\' for the same and added functionality.")
            if (args$un) which <- "unadjusted"
            else which <- weight.names
        }
        else {
            if (is_null(X$weights) && is_null(X$subclass)) which <- "unadjusted"
            else which <- weight.names
        }
    }
    else {
        if (is_null(X$weights) && is_null(X$subclass)) which <- "unadjusted"
        else {
            which <- tolower(which)
            which <- match_arg(which, unique(c("adjusted", "unadjusted", "both", weight.names)),
                               several.ok = TRUE)
        }
    }
    
    if (is_not_null(args$size.weight)) {
        message("Note: \'size.weight\' is no longer allowed; please use \'alpha.weight\' for similar functionality.")
    }
    
    title <- paste0("Distributional Balance for \"", var.name, "\"")
    subtitle <- NULL
    
    facet <- NULL
    is.categorical.var <- is.factor(X$var) || is.character(X$var) || is_binary(X$var) 
    
    if (is_null(X$subclass) || (length(which) == 1 && which == "unadjusted")) {
        if (is_not_null(which.sub) && !all(is.na(which.sub))) {
            if (is_null(X$subclass)) warning("which.sub was specified but no subclasses were supplied. Ignoring which.sub.", call. = FALSE)
            else warning("which.sub was specified but only the unadjusted sample was requested. Ignoring which.sub.", call. = FALSE)
        }
        
        facet <- "which"
        
        if ("both" %in% which) which <- c("Unadjusted Sample", weight.names)
        else {
            if ("adjusted" %in% which) which <- c(which[which != "adjusted"], weight.names)
            if ("unadjusted" %in% which) which <- c("Unadjusted Sample", which[which != "unadjusted"])
        }
        which <- unique(which)
        
        if (is_null(X$weights)) X$weights <- setNames(data.frame(rep.int(1, length(X$treat))),
                                                      "Unadjusted Sample")
        else {
            if (ncol(X$weights) == 1) {
                which[which != "Unadjusted Sample"] <- "Adjusted Sample"
                names(X$weights) <- "Adjusted Sample"
            }
            
            if ("Unadjusted Sample" %in% which) {
                X$weights <- setNames(data.frame(rep.int(max(X$weights), length(X$treat)),
                                                 X$weights),
                                      c("Unadjusted Sample", names(X$weights)))
            }
            
            
        }
        
        X$weights <- X$weights[which]
        
        #Process sample names
        ntypes <- length(which)
        nadj <- sum(which != "Unadjusted Sample")
        if (!missing(sample.names)) {
            if (!is.vector(sample.names, "character")) {
                warning("The argument to sample.names must be a character vector. Ignoring sample.names.", call. = FALSE)
                sample.names <- NULL
            }
            else if (length(sample.names) %nin% c(ntypes, nadj)) {
                warning("The argument to sample.names must contain as many names as there are sample types, or one fewer. Ignoring sample.names.", call. = FALSE)
                sample.names <- NULL
            }
        }
        else sample.names <- NULL
        
        #NULL: all; NA: none
        in.imp <- rep.int(TRUE, length(X$var))
        if (is_not_null(X$imp)) {
            if (is_null(which.imp) || all(is.na(which.imp))) {
                in.imp <- !is.na(X$imp)
            }
            else {
                if (is.numeric(which.imp)) {
                    if (all(which.imp %in% seq_len(nlevels(X$imp)))) {
                        in.imp <- !is.na(X$imp) & X$imp %in% levels(X$imp)[which.imp]
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
        
        in.cluster <- rep.int(TRUE, length(X$var))
        if (is_not_null(X$cluster)) {
            if (is_null(which.cluster)|| all(is.na(which.cluster))) {
                in.cluster <- !is.na(X$cluster)
            }
            else {
                if (is.numeric(which.cluster)) {
                    if (all(which.cluster %in% seq_len(nlevels(X$cluster)))) {
                        in.cluster <- !is.na(X$cluster) & X$cluster %in% levels(X$cluster)[which.cluster]
                    }
                    else {
                        stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word_list(which.cluster[!which.cluster %in% seq_len(nlevels(X$cluster))])), call. = FALSE)
                    }
                }
                else if (is.character(which.cluster)) {
                    if (all(which.cluster %in% levels(X$cluster))) {
                        in.cluster <- !is.na(X$cluster) & X$cluster %in% which.cluster
                    }
                    else {
                        stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word_list(which.cluster[which.cluster %nin% levels(X$cluster)])), call. = FALSE)
                    }
                }
                else stop("The argument to which.cluster must be the names or indices corresponding to the clusters for which distributions are to be displayed.", call. = FALSE)
            }
            facet <- c("cluster", facet)
        }
        else if (is_not_null(which.cluster)) {
            warning("which.cluster was specified but no cluster values were supplied. Ignoring which.cluster.", call. = FALSE)
        }
        
        in.time <- rep.int(TRUE, length(X$var))
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
        
        Q <- make_list(which)
        for (i in which) {
            Q[[i]] <- setNames(as.data.frame(matrix(0, ncol = 4, nrow = nobs)),
                               c("treat", "var", "weights", "which"))
            Q[[i]]$treat <- X$treat[in.imp & in.cluster & in.time]
            Q[[i]]$var <- X$var[in.imp & in.cluster & in.time]
            Q[[i]]$weights <-  X$weights[in.imp & in.cluster & in.time, i]
            Q[[i]]$which <- i
            
            #Add columns for additional facets
            if ("imp" %in% facet) Q[[i]]$imp <- factor(paste("Imputation", X$imp[in.imp & in.cluster & in.time]))
            if ("cluster" %in% facet) Q[[i]]$cluster <- factor(X$cluster[in.imp & in.cluster & in.time])
            if ("time" %in% facet) Q[[i]]$time <- factor(paste("Time", X$time[in.imp & in.cluster & in.time]))
        }
        D <- do.call(rbind, Q)
        
        D$which <- factor(D$which, levels = which)
        
        if (is_not_null(sample.names)) {
            if (length(sample.names) == nadj) {
                levels(D$which)[levels(D$which) != "Unadjusted Sample"] <- sample.names
            }
            else if (length(sample.names) == ntypes) {
                levels(D$which) <- sample.names
            }
        }
        
    }
    else {
        if (is_not_null(X$cluster)) stop("Subclasses are not supported with clusters.", call. = FALSE)
        if (is_not_null(X$imp)) stop("Subclasses are not supported with multiple imputations.", call. = FALSE)
        
        #Process sample names
        # ntypes <- as.numeric("both" %in% which)
        if (!missing(sample.names)) {
            if (which %nin% c("both", "unadjusted")) {
                warning("sample.names can only be used with which = \"both\" or \"unadjusted\" to rename the unadjusted sample when called with subclasses. Ignoring sample.names.", call. = FALSE)
                sample.names <- NULL
            }
            else if (!is.vector(sample.names, "character")) {
                warning("The argument to sample.names must be a character vector. Ignoring sample.names.", call. = FALSE)
                sample.names <- NULL
            }
            else if (length(sample.names) != 1) {
                warning("The argument to sample.names must be of length 1 when called with subclasses. Ignoring sample.names.", call. = FALSE)
                sample.names <- NULL
            }
        }
        else sample.names <- NULL
        
        #Get sub.names.good
        sub.names <- levels(X$subclass)
        
        if (is_null(which.sub)) {
            which.sub <- sub.names
        }
        else {
            which.sub.original <- which.sub
            if (anyNA(which.sub)) which.sub <- which.sub[!is.na(which.sub)]
            
            if (is_null(which.sub)) {
                stop(paste0("The argument to which.sub cannot be .none or NA when which = \"", which, "\"."), call. = FALSE)
            }
            else if (is.character(which.sub)) {
                which.sub <- which.sub[which.sub %in% sub.names]
            }
            else if (is.numeric(which.sub)) {
                which.sub <- sub.names[which.sub[which.sub %in% seq_along(sub.names)]]
            }
            
            if (is_null(which.sub)) {
                stop("The argument to which.sub must be .none, .all, or the valid names or indices of subclasses.", call. = FALSE)
            }
            else if (any(which.sub.original %nin% which.sub)) {
                w.l <- word_list(which.sub.original[which.sub.original %nin% which.sub])
                warning(paste(w.l, ifelse(attr(w.l, "plural"), "do", "does"), "not correspond to any subclass in the object and will be ignored."), call. = FALSE)
            }
        }
        
        in.sub <- !is.na(X$subclass) & X$subclass %in% which.sub
        D <- setNames(as.data.frame(matrix(0, nrow = sum(in.sub), ncol = 4)),
                      c("weights", "treat", "var", "subclass"))
        D$weights <- rep(1, NROW(D))
        D$treat <- X$treat[in.sub]
        D$var <- X$var[in.sub]
        D$subclass <- paste("Subclass", X$subclass[in.sub])
        #title <- paste0(title, "\nin subclass ", which.sub)
        
        if (which == "both") {
            #Make unadjusted sample
            D2 <- setNames(as.data.frame(matrix(0, nrow = length(X$treat), ncol = 4)),
                           c("weights", "treat", "var", "subclass"))
            D2$weights <- rep(1, NROW(D2))
            D2$treat <- X$treat
            D2$var <- X$var
            D2$subclass <- rep("Unadjusted Sample", length(X$treat))
            D <- rbind(D2, D, stringsAsFactors = TRUE)
            D$subclass <- relevel(factor(D$subclass), "Unadjusted Sample")
        }
        
        subtitle <- NULL
        
        facet <- "subclass"
        
        if (is_not_null(sample.names)) {
            levels(D$subclass)[levels(D$subclass) == "Unadjusted Sample"] <- sample.names
        }
        
    }
    
    # if ("which" %in% facet) {
    #     if (length(which) == 1) {
    #         subtitle <- levels(D$which)[1]
    #         facet <- facet[facet %nin% "which"]
    #     }
    # }
    
    treat.type <- get.treat.type(assign.treat.type(D$treat))
    
    D <- D[order(D$var),]
    
    if (treat.type == "continuous") { #Continuous treatments
        
        if (is.categorical.var) {
            D$weights <- ave(D[["weights"]], 
                             D[c("var", facet)], 
                             FUN = function(x) x/sum(x))
        }
        
        D$treat <- as.numeric(D$treat)
        
        if (is.categorical.var) { #Categorical vars
            D$var <- factor(D$var)
            cat.sizes <- tapply(rep(1, NROW(D)), D$var, sum)
            smallest.cat <- names(cat.sizes)[which.min(cat.sizes)]
            if (is.character(bw)) {
                if (is.function(get0(paste0("bw.", bw)))) {
                    bw <- get0(paste0("bw.", bw))(D$treat[D$var == smallest.cat])
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
                
                if (!all(vapply(colors, isColor, logical(1L)))) {
                    warning("The argument to colors contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
                    colors <- gg_color_hue(ntypes)
                }
                
            }
            
            bp <- ggplot(D, mapping = aes(x = treat, fill = var, weight = weights)) + 
                geom_density(alpha = .4, bw = bw, adjust = adjust, kernel = kernel, n = n, trim = TRUE, outline.type = "full") + 
                labs(fill = var.name, y = "Density", x = "Treat", title = title, subtitle = subtitle) +
                scale_fill_manual(values = colors) + geom_hline(yintercept = 0) +
                scale_y_continuous(expand = expansion(mult = c(0, .05)))
        }
        else { #Continuous vars
            D$var.mean <- ave(D[["var"]], D[facet], FUN = mean)
            D$treat.mean <- ave(D[["treat"]], D[facet], FUN = mean)
            
            bp <- ggplot(D, mapping = aes(x = var, y = treat, weight = weights))
            
            if (identical(which, "Unadjusted Sample") || isFALSE(alpha.weight)) bp <- bp + geom_point(alpha = .9)
            else bp <- bp + geom_point(aes(alpha = weights), show.legend = FALSE) + 
                scale_alpha(range = c(.04, 1))
            
            bp <- bp + 
                geom_smooth(method = "lm", se = FALSE, color = "firebrick2", alpha = .4, size = 1.5) + 
                geom_smooth(method = "loess", se = FALSE, color = "royalblue1", alpha = .1, size = 1.5) +
                # geom_smooth(method = "gam", se = FALSE, color = "royalblue1", alpha = .1, size = 1.5, formula = y ~ s(x, bs = "cs")) +
                geom_hline(aes(yintercept = treat.mean), linetype = 1, alpha = .9) + 
                geom_vline(aes(xintercept = var.mean), linetype = 1, alpha = .8) +
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
        else if (anyNA(which.treat)) {
            which.treat <- character(0)
        }
        else {
            warning("The argument to which.treat must be NA, NULL, or a vector of treatment names or indices. All treatment groups will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
        if (is_not_null(which.treat) && !anyNA(which.treat)) D <- D[D$treat %in% which.treat,]
        
        for (i in names(D)[sapply(D, is.factor)]) D[[i]] <- factor(D[[i]])
        
        D$weights <- ave(D[["weights"]], 
                         D[c("treat", facet)], 
                         FUN = function(x) x/sum(x))
        
        #Color
        ntypes <- nlevels.treat <- nlevels(D$treat)
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
        names(colors) <- levels(D$treat)
        
        if (is_binary(D$var) || is.factor(D$var) || is.character(D$var)) { #Categorical vars
            D$var <- factor(D$var)
            bp <- ggplot(D, mapping = aes(x = var, fill = treat, weight = weights)) + 
                geom_bar(position = "dodge", alpha = .8, color = "black") + 
                labs(x = var.name, y = "Proportion", fill = "Treatment", title = title, subtitle = subtitle) + 
                scale_x_discrete(drop=FALSE) + scale_fill_manual(drop=FALSE, values = colors) +
                geom_hline(yintercept = 0) + 
                scale_y_continuous(expand = expansion(mult = c(0, .05)))
        }
        else { #Continuous vars
            
            type <- match_arg(type, c("histogram", "density", "ecdf"))
            
            if (type %in% c("ecdf")) {
                mirror <- FALSE
                alpha <- 1
                legend.alpha <- alpha
                expandScale <- expansion(mult = .005)
            }
            else if (nlevels.treat <= 2 && mirror) {
                posneg <- c(1, -1)
                alpha <- .8
                legend.alpha <- alpha
                expandScale <- expansion(mult = .05)
            }
            else {
                mirror <- FALSE
                posneg <- rep(1, nlevels.treat)
                alpha <- .4
                legend.alpha <- alpha/nlevels.treat
                expandScale <- expansion(mult = c(0, .05))
            }
            
            if (type == "histogram") {
                if (disp.means) {
                    D$var.mean <- ave(D[c("var", "weights")], D[c("treat", facet)], 
                                      FUN = function(x) w.m(x[[1]], x[[2]]))[[1]]
                }
                
                if (is_null(args$bins) || !is.numeric(args$bins)) args$bins <- 12
                geom_fun <- function(t) {
                    out <- list(geom_histogram(data = D[D$treat == levels(D$treat)[t],],
                                               mapping = aes(x = var, y = posneg[t]*stat(count), weight = weights,
                                                             fill = names(colors)[t]),
                                               alpha = alpha, bins = args$bins, color = "black"),
                                NULL)
                    if (disp.means) out[[2]] <-
                            geom_segment(data = unique(D[D$treat == levels(D$treat)[t], c("var.mean", facet), drop = FALSE]),
                                         mapping = aes(x = var.mean, xend = var.mean, y = 0, yend = posneg[t]*Inf), 
                                         color = if (mirror) "black" else colors[[t]])
                    out
                }
                ylab <- "Proportion"
                
            }
            else if (type %in% c("ecdf")) {
                
                D$cum.pt <- ave(D[["weights"]], 
                                D[c("treat", facet)], 
                                FUN = function(x) cumsum(x)/sum(x))
                
                #Pad 0 and 1
                extra <- setNames(do.call(expand.grid, c(list(c("top", "bottom")),
                                                         lapply(c("treat", facet), function(i) levels(D[[i]])))),
                                  c("pos_", "treat", facet))
                
                merge.extra <- data.frame(pos_ = c("top", "bottom"),
                                          var = c(-Inf, Inf),
                                          cum.pt = c(0, 1))
                extra <- merge(extra, merge.extra)
                extra[["pos_"]] <- NULL
                
                D[names(D) %nin% names(extra)] <- NULL
                
                D <- rbind(extra[names(D)], D)
                
                geom_fun <- function(t) {
                    
                    geom_step(data = D[D$treat == levels(D$treat)[t],],
                              mapping = aes(x = var, y = cum.pt, color = names(colors)[t]))
                    
                }
                ylab <- "Cumulative Proportion"
            }
            else {
                t.sizes <- tapply(rep(1, NROW(D)), D$treat, sum)
                smallest.t <- names(t.sizes)[which.min(t.sizes)]
                if (is.character(bw)) {
                    if (is.function(get0(paste0("bw.", bw)))) {
                        bw <- get0(paste0("bw.", bw))(D$var[D$treat == smallest.t])
                    }
                    else {
                        stop(paste(bw, "is not an acceptable entry to bw. See ?stats::density for allowable options."), call. = FALSE)
                    }
                }
                
                if (disp.means) {
                    D$var.mean <- ave(D[c("var", "weights")], D[c("treat", facet)], 
                                      FUN = function(x) w.m(x[[1]], x[[2]]))[[1]]
                }
                
                geom_fun <- function(t) {
                    out <- list(
                        geom_density(data = D[D$treat == levels(D$treat)[t],],
                                     mapping = aes(x = var, y = posneg[t]*stat(density), weight = weights,
                                                   fill = names(colors)[t]),
                                     alpha = alpha, bw = bw, adjust = adjust,
                                     kernel = kernel, n = n, trim = TRUE,
                                     outline.type = "full"),
                        NULL)
                    if (disp.means) out[[2]] <-
                            geom_segment(data = unique(D[D$treat == levels(D$treat)[t], c("var.mean", facet), drop = FALSE]),
                                         mapping = aes(x = var.mean, xend = var.mean, y = 0, yend = posneg[t]*Inf), 
                                         color = if (mirror) "black" else colors[[t]])
                    out
                }
                ylab <- "Density"
                
            }
            
            bp <- Reduce("+", c(list(ggplot()),
                                do.call(c, lapply(seq_len(nlevels.treat), geom_fun)))) +
                scale_fill_manual(values = colors, guide = guide_legend(override.aes = list(alpha = legend.alpha, size = .25))) +
                scale_color_manual(values = colors) +
                labs(x = var.name, y = ylab, title = title, subtitle = subtitle,
                     fill = "Treatment", color = "Treatment") + 
                scale_y_continuous(expand = expandScale)
            
            if (mirror) bp <- bp + geom_hline(yintercept = 0)
        }
    }
    
    #Themes
    bp <- bp + theme(panel.background = element_rect(fill = "white", color = "black"),
                     panel.border = element_rect(fill = NA, color = "black"),
                     plot.background = element_blank(),
                     legend.position = position,
                     legend.key=element_rect(fill = "white", color = "black", size = .25))
    if (isFALSE(grid)) {
        bp <- bp + theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())
    }
    else {
        bp <- bp + theme(panel.grid.major = element_line(color = "gray87"),
                         panel.grid.minor = element_line(color = "gray90"))
    }
    
    if (is_not_null(facet)) {
        if (is_not_null(facet.formula)) {
            if (!is.formula(facet.formula)) {
                stop("facet.formula must be a formula.", call. = FALSE)
            }
            test.facet <- invisible(facet_grid(facet.formula))
            if (any(c(names(test.facet$params$rows), names(test.facet$params$cols)) %nin% facet)) {
                stop(paste("Only", word_list(facet, is.are = TRUE, quotes = TRUE), "allowed in facet.formula."), call. = FALSE)
            }
            if ("which" %nin% c(names(test.facet$params$rows), names(test.facet$params$cols))) {
                if (length(which) > 1) stop("\"which\" must be in the facet formula when the which argument refers to more than one sample.", call. = FALSE)
                else message(paste0("Displaying balance for the ", if (which %in% c("Adjusted Sample", "Unadjusted Sample")) tolower(which) else paste(which, "sample"), "."))
            }
        }
        else if (length(facet) >= 2) {
            if ("which" %in% facet) {
                facet.formula <- f.build(facet[facet %nin% "which"], "which")
            }
            else if ("imp" %in% facet) {
                facet.formula <- f.build("imp", facet[facet %nin% "imp"])
            }
            else {
                facets <- data.frame(facet = facet, length = vapply(facet, function(x) nlevels(D[[x]]), numeric(1L)),
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