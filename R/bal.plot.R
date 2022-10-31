bal.plot <- function(x, var.name, ..., which, which.sub = NULL, cluster = NULL, which.cluster = NULL, 
                     imp = NULL, which.imp = NULL, which.treat = NULL, which.time = NULL, 
                     mirror = FALSE, type = "density", colors = NULL, grid = FALSE, sample.names,
                     position = "right", facet.formula = NULL, disp.means = getOption("cobalt_disp.means", FALSE), 
                     alpha.weight = TRUE) {
    
    tryCatch(identity(x), error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    .alls <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.all)), logical(1L))
    .nones <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.none)), logical(1L))
    if (any(c(.alls, .nones))) {
        .call[.alls] <- expression(NULL)
        .call[.nones] <- expression(NA)
        return(eval.parent(.call))
    }
    
    tryCatch(force(x), error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    args <- list(...)
    
    x <- process_obj(x)
    
    X <- x2base(x, ..., cluster = cluster, imp = imp)
    
    if (is_null(X$covs.list)) {
        #Point treatment
        X$covs <- get.C2(X$covs, addl = X$addl, distance = X$distance, cluster = X$cluster, treat = X$treat,
                         drop = FALSE)
        co.names <- attr(X$covs, "co.names")
        if (missing(var.name)) {
            var.name <- NULL; k = 1
            while(is_null(var.name)) {
                x <- co.names[[k]]
                if ("isep" %nin% x[["type"]]) var.name <- x[["component"]][x[["type"]] == "base"][1]
                else {
                    if (k < length(co.names)) k <- k + 1
                    else 
                        stop("Please specify an argument to 'var.name'.", call. = FALSE)
                }
            }
            message(paste0("No 'var.name' was provided. Displaying balance for ", var.name, "."))
        }
        var.name_in_name <- vapply(co.names, function(x) var.name %in% x[["component"]][x[["type"]] == "base"] &&
                                       "isep" %nin% x[["type"]], logical(1L))
        var.name_in_name_and_factor <- vapply(seq_along(co.names), function(x) var.name_in_name[x] && "fsep" %in% co.names[[x]][["type"]], logical(1L))
        if (any(var.name_in_name_and_factor)) {
            X$var <- unsplitfactor(as.data.frame(X$covs[,var.name_in_name_and_factor, drop = FALSE]), 
                                   var.name, sep = attr(co.names, "seps")["factor"])[[1]]
        }
        else if (any(var.name_in_name)) {
            X$var <- X$covs[,var.name]
        }
        else {
            stop(paste0("\"", var.name, "\" is not the name of an available covariate."), call. = FALSE)
        }
        
        if (get.treat.type(X$treat) != "continuous") X$treat <- treat_vals(X$treat)[X$treat]
    }
    else {
        #Longitudinal
        X$covs.list <- lapply(seq_along(X$covs.list), function(i) {
            get.C2(X$covs.list[[i]], addl = X$addl.list[[i]], distance = X$distance.list[[i]], cluster = X$cluster,
                   treat = X$treat.list[[i]], drop = FALSE)
        })
        co.names.list <- lapply(X$covs.list, attr, "co.names")
        ntimes <- length(X$covs.list)
        
        if (missing(var.name)) {
            var.name <- NULL; k <- 1; time <- 1
            while (is_null(var.name)) {
                x <- co.names.list[[time]][[k]]
                if ("isep" %nin% x[["type"]]) {
                    var.name <- x[["component"]][x[["type"]] == "base"][1]
                }
                else {
                    if (time < ntimes) {
                        if (k <  length(co.names.list[[time]])) k <- k + 1
                        else {
                            k <- 1
                            time <- time + 1
                        }
                    }
                    else if (k < length(co.names.list[[time]]))  k <- k + 1
                    else stop("Please specified an argument to 'var.name'.", call. = FALSE)
                }
            }
            message(paste0("No 'var.name' was provided. Displaying balance for ", var.name, "."))
        }
        
        var.list <- make_list(length(X$covs.list))
        appears.in.time <- rep.int(TRUE, length(X$covs.list))
        for (i in seq_along(X$covs.list)) {
            var.name_in_name <- vapply(co.names.list[[i]], function(x) var.name %in% x[["component"]][x[["type"]] == "base"] &&
                                           "isep" %nin% x[["type"]], logical(1L))
            var.name_in_name_and_factor <- var.name_in_name & vapply(co.names.list[[i]], function(x) "fsep" %in% x[["type"]], logical(1L))
            if (any(var.name_in_name_and_factor)) {
                var.list[[i]] <- unsplitfactor(as.data.frame(X$covs.list[[i]][,var.name_in_name_and_factor, drop = FALSE]), 
                                               var.name, sep = attr(co.names.list[[i]], "seps")["factor"])[[1]]
            }
            else if (any(var.name_in_name)) {
                var.list[[i]] <- X$covs.list[[i]][,var.name]
            }
            else {
                appears.in.time[i] <- FALSE
            }
        }
        if (all(vapply(var.list, is_null, logical(1L)))) stop(paste0("\"", var.name, "\" is not the name of an available covariate."), call. = FALSE)
        
        X$var <- unlist(var.list[appears.in.time])
        
        X$time <- rep(which(appears.in.time), times = lengths(var.list[appears.in.time]))
        
        X$treat.list[appears.in.time] <- lapply(X$treat.list[appears.in.time], function(t) if (get.treat.type(t) != "continuous") treat_vals(t)[t] else t)
        X$treat <- unlist(X$treat.list[appears.in.time])
        if (is_not_null(names(X$treat.list)[appears.in.time])) treat.names <- names(X$treat.list)[appears.in.time]
        else treat.names <- which(appears.in.time)
        
        if (is_not_null(X$weights)) X$weights <- do.call("rbind", lapply(seq_len(sum(appears.in.time)), function(x) X$weights))
        if (is_not_null(X$cluster)) X$cluster <- rep(X$cluster, sum(appears.in.time))
        if (is_not_null(X$imp)) X$imp <- rep(X$imp, sum(appears.in.time))
    }
    
    if (is_null(X$subclass)) {
        if (NCOL(X$weights) == 1L) weight.names <- "adjusted"
        else weight.names <- names(X$weights)
    }
    else weight.names <- "adjusted"
    
    if (is_null(X$s.weights)) X$s.weights <- rep(1, length(X$treat))
    
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
            if (is_null(which.imp)) {
                in.imp <- !is.na(X$imp)
                facet <- c("imp", facet)
            }
            else if (!all(is.na(which.imp))) {
                if (is.numeric(which.imp)) {
                    if (all(which.imp %in% seq_len(nlevels(X$imp)))) {
                        in.imp <- !is.na(X$imp) & X$imp %in% levels(X$imp)[which.imp]
                    }
                    else {
                        stop(paste0("The following inputs to 'which.imp' do not correspond to given imputations:\n\t", word_list(which.imp[!which.imp %in% seq_len(nlevels(X$imp))])), call. = FALSE)
                    }
                    facet <- c("imp", facet)
                }
                else stop("The argument to 'which.imp' must be the indices corresponding to the imputations for which distributions are to be displayed.", call. = FALSE)
            }
            
        }
        else if (is_not_null(which.imp) && !all(is.na(which.imp))) {
            warning("'which.imp' was specified but no 'imp' values were supplied. Ignoring 'which.imp'.", call. = FALSE)
        }
        
        in.cluster <- rep.int(TRUE, length(X$var))
        if (is_not_null(X$cluster)) {
            if (is_null(which.cluster)) {
                in.cluster <- !is.na(X$cluster)
                facet <- c("cluster", facet)
            }
            else if (!all(is.na(which.cluster))) {
                if (is.numeric(which.cluster)) {
                    if (all(which.cluster %in% seq_len(nlevels(X$cluster)))) {
                        in.cluster <- !is.na(X$cluster) & X$cluster %in% levels(X$cluster)[which.cluster]
                    }
                    else {
                        stop(paste0("The following inputs to 'which.cluster' do not correspond to given clusters:\n\t", word_list(which.cluster[!which.cluster %in% seq_len(nlevels(X$cluster))])), call. = FALSE)
                    }
                }
                else if (is.character(which.cluster)) {
                    if (all(which.cluster %in% levels(X$cluster))) {
                        in.cluster <- !is.na(X$cluster) & X$cluster %in% which.cluster
                    }
                    else {
                        stop(paste0("The following inputs to 'which.cluster' do not correspond to given clusters:\n\t", word_list(which.cluster[which.cluster %nin% levels(X$cluster)])), call. = FALSE)
                    }
                }
                else stop("The argument to 'which.cluster' must be the names or indices corresponding to the clusters for which distributions are to be displayed.", call. = FALSE)
                facet <- c("cluster", facet)
            }
        }
        else if (is_not_null(which.cluster)) {
            warning("'which.cluster' was specified but no 'cluster' values were supplied. Ignoring 'which.cluster'.", call. = FALSE)
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
                        stop(paste0("The following inputs to 'which.time' do not correspond to given time periods:\n\t", word_list(which.time[!which.time %in% seq_along(X$covs.list)])), call. = FALSE)
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
                        stop(paste0("The following inputs to 'which.time' do not correspond to given time periods:\n\t", word_list(which.time[!which.time %in% treat.names])), call. = FALSE)
                    }
                }
                else stop("The argument to 'which.time' must be the names or indices corresponding to the time periods for which distributions are to be displayed.", call. = FALSE)
            }
            facet <- c("time", facet)
        }
        else if (is_not_null(which.time)) {
            warning("'which.time' was specified but a point treatment was supplied. Ignoring 'which.time'.", call. = FALSE)
        }
        
        nobs <- sum(in.imp & in.cluster & in.time)
        if (nobs == 0) stop("No observations to display.", call. = FALSE)
      
        D <- make_list(which)
        for (i in which) {
            D[[i]] <- make_df(c("treat", "var", "weights", "s.weights", "which"), nobs)
            D[[i]]$treat <- X$treat[in.imp & in.cluster & in.time]
            D[[i]]$var <- X$var[in.imp & in.cluster & in.time]
            D[[i]]$weights <-  X$weights[in.imp & in.cluster & in.time, i]
            D[[i]]$s.weights <- X$s.weights[in.imp & in.cluster & in.time]
            D[[i]]$which <- i
            
            #Add columns for additional facets
            if ("imp" %in% facet) D[[i]]$imp <- factor(paste("Imputation", X$imp[in.imp & in.cluster & in.time]))
            if ("cluster" %in% facet) D[[i]]$cluster <- factor(X$cluster[in.imp & in.cluster & in.time])
            if ("time" %in% facet) D[[i]]$time <- factor(paste("Time", X$time[in.imp & in.cluster & in.time]))
        }
        D <- do.call(rbind, D)
        
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
                warning("'sample.names' can only be used with which = \"both\" or \"unadjusted\" to rename the unadjusted sample when called with subclasses. Ignoring 'sample.names'.", call. = FALSE)
                sample.names <- NULL
            }
            else if (!is.vector(sample.names, "character")) {
                warning("The argument to 'sample.names' must be a character vector. Ignoring 'sample.names'.", call. = FALSE)
                sample.names <- NULL
            }
            else if (length(sample.names) != 1) {
                warning("The argument to 'sample.names' must be of length 1 when called with subclasses. Ignoring 'sample.names'.", call. = FALSE)
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
                stop(paste0("The argument to 'which.sub' cannot be .none or NA when which = \"", which, "\"."), call. = FALSE)
            }
            else if (is.character(which.sub)) {
                which.sub <- which.sub[which.sub %in% sub.names]
            }
            else if (is.numeric(which.sub)) {
                which.sub <- sub.names[which.sub[which.sub %in% seq_along(sub.names)]]
            }
            
            if (is_null(which.sub)) {
                stop("The argument to 'which.sub' must be .none, .all, or the valid names or indices of subclasses.", call. = FALSE)
            }
            else if (any(which.sub.original %nin% which.sub)) {
                w.l <- word_list(which.sub.original[which.sub.original %nin% which.sub])
                warning(paste(w.l, ifelse(attr(w.l, "plural"), "do", "does"), "not correspond to any subclass in the object and will be ignored."), call. = FALSE)
            }
        }
        
        in.sub <- !is.na(X$subclass) & X$subclass %in% which.sub
        D <- make_df(c("weights", "s.weights", "treat", "var", "subclass"), sum(in.sub))
        D$weights <- 1
        D$s.weights <- X$s.weights[in.sub]
        D$treat <- X$treat[in.sub]
        D$var <- X$var[in.sub]
        D$subclass <- paste("Subclass", X$subclass[in.sub])

        if (which == "both") {
            #Make unadjusted sample
            D2 <- make_df(c("weights", "s.weights", "treat", "var", "subclass"), length(X$treat))
            D2$weights <- 1
            D$s.weights <- X$s.weights
            D2$treat <- X$treat
            D2$var <- X$var
            D2$subclass <- rep("Unadjusted Sample", length(X$treat))
            D <- rbind(D2, D, stringsAsFactors = TRUE)
            D$subclass <- relevel(factor(D$subclass), "Unadjusted Sample")
        }
        
        facet <- "subclass"
        
        if (is_not_null(sample.names)) {
            levels(D$subclass)[levels(D$subclass) == "Unadjusted Sample"] <- sample.names
        }
        
    }
    
    treat.type <- get.treat.type(assign.treat.type(D$treat))
    
    D <- na.omit(D[order(D$var),])
    # D <- D[D$weights > 0,]
    
    if (treat.type == "continuous") { #Continuous treatments

        D$treat <- as.numeric(D$treat)
        
        if (is.categorical.var) { #Categorical vars
            #Density arguments supplied through ...
            bw <- if_null_then(args$bw, "nrd0")
            adjust <- if_null_then(args$adjust, 1)
            kernel <- if_null_then(args$kernel, "gaussian")
            n <- if_null_then(args$n, 512)
            
            D$var <- factor(D$var)
            cat.sizes <- tapply(rep(1, NROW(D)), D$var, sum)
            smallest.cat <- names(cat.sizes)[which.min(cat.sizes)]
            
            if (is.character(bw)) {
                if (is.function(get0(paste0("bw.", bw)))) {
                    bw <- get0(paste0("bw.", bw))(D$treat[D$var == smallest.cat])
                }
                else {
                    stop(paste(bw, "is not an acceptable entry to 'bw'. See ?stats::density for allowable options."), call. = FALSE)
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
                    warning(paste("Only using first", ntypes, "values in 'colors'."), call. = FALSE)
                }
                else if (length(colors) < ntypes) {
                    warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
                    colors <- gg_color_hue(ntypes)
                }
                
                if (!all(vapply(colors, isColor, logical(1L)))) {
                    warning("The argument to 'colors' contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
                    colors <- gg_color_hue(ntypes)
                }
                
            }
            
            D$weights <- ave(D[["weights"]] * D[["s.weights"]], 
                             D[c("var", facet)], 
                             FUN = function(x) x/sum(x))
            
            bp <- ggplot2::ggplot(D, mapping = aes(x = .data$treat, fill = .data$var, weight = .data$weights)) + 
                ggplot2::geom_density(alpha = .4, bw = bw, adjust = adjust, kernel = kernel, 
                                      n = n, trim = TRUE, outline.type = "full", stat = StatDensity2) + 
                ggplot2::labs(fill = var.name, y = "Density", x = "Treat", title = title, subtitle = subtitle) +
                ggplot2::scale_fill_manual(values = colors) + 
                ggplot2::geom_hline(yintercept = 0) +
                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .05)))
        }
        else { #Continuous vars
            D$var.mean <- ave(D[c("var", "s.weights")], D[facet], 
                              FUN = function(x) w.m(x[[1]], x[[2]]))[[1]]
            D$treat.mean <- ave(D[c("treat", "s.weights")], D[facet], 
                                FUN = function(x) w.m(x[[1]], x[[2]]))[[1]]
            
            bp <- ggplot2::ggplot(D, mapping = aes(x = .data$var, y = .data$treat, weight = .data$weights * .data$s.weights))
            
            if (identical(which, "Unadjusted Sample") || isFALSE(alpha.weight)) bp <- bp + ggplot2::geom_point(alpha = .9)
            else bp <- bp + ggplot2::geom_point(aes(alpha = .data$weights), show.legend = FALSE) + 
                ggplot2::scale_alpha(range = c(.04, 1))
            
            bp <- bp + 
                ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "firebrick2", alpha = .4, size = 1.5) + 
                {
                    if (nrow(D) <= 1000)
                        ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "royalblue1", alpha = .1, size = 1.5) 
                    else
                        ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "royalblue1", alpha = .1, size = 1.5)
                } +
                ggplot2::geom_hline(aes(yintercept = .data$treat.mean), linetype = 1, alpha = .9) + 
                ggplot2::geom_vline(aes(xintercept = .data$var.mean), linetype = 1, alpha = .8) +
                ggplot2::labs(y = "Treat", x = var.name, title = title, subtitle = subtitle)
        }
    }
    else { #Categorical treatments (multinomial supported)
        D$treat <- factor(D$treat)
        
        if (is_null(which.treat)) 
            which.treat <- character(0)
        else if (is.numeric(which.treat)) {
            which.treat <- levels(D$treat)[seq_along(levels(D$treat)) %in% which.treat]
            if (is_null(which.treat)) {
                warning("No numbers in 'which.treat' correspond to treatment values. All treatment groups will be displayed.", call. = FALSE)
                which.treat <- character(0)
            }
        }
        else if (is.character(which.treat)) {
            which.treat <- levels(D$treat)[levels(D$treat) %in% which.treat]
            if (is_null(which.treat)) {
                warning("No names in 'which.treat' correspond to treatment values. All treatment groups will be displayed.", call. = FALSE)
                which.treat <- character(0)
            }
        }
        else if (anyNA(which.treat)) {
            which.treat <- character(0)
        }
        else {
            warning("The argument to 'which.treat' must be NA, NULL, or a vector of treatment names or indices. All treatment groups will be displayed.", call. = FALSE)
            which.treat <- character(0)
        }
        if (is_not_null(which.treat) && !anyNA(which.treat)) D <- D[D$treat %in% which.treat,]
        
        for (i in names(D)[vapply(D, is.factor, logical(1L))]) D[[i]] <- factor(D[[i]])
        
        D$weights <- ave(D[["weights"]] * D[["s.weights"]], 
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
                warning(paste("Only using first", ntypes, "values in 'colors'."), call. = FALSE)
            }
            else if (length(colors) < ntypes) {
                warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
                colors <- gg_color_hue(ntypes)
            }
            
            if (!all(vapply(colors, isColor, logical(1L)))) {
                warning("The argument to 'colors' contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
                colors <- gg_color_hue(ntypes)
            }
        }
        names(colors) <- levels(D$treat)
        
        if (is.categorical.var) { #Categorical vars
            D$var <- factor(D$var)
            bp <- ggplot2::ggplot(D, mapping = aes(x = .data$var, fill = .data$treat, weight = .data$weights)) + 
                ggplot2::geom_bar(position = "dodge", alpha = .8, color = "black") + 
                ggplot2::labs(x = var.name, y = "Proportion", fill = "Treatment", title = title, subtitle = subtitle) + 
                ggplot2::scale_x_discrete(drop=FALSE) + 
                ggplot2::scale_fill_manual(drop=FALSE, values = colors) +
                ggplot2::geom_hline(yintercept = 0) + 
                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .05)))
        }
        else { #Continuous vars
            
            type <- match_arg(type, c("histogram", "density", "ecdf"))
            
            if (type %in% c("ecdf")) {
                mirror <- FALSE
                alpha <- 1
                legend.alpha <- alpha
                expandScale <- ggplot2::expansion(mult = .005)
            }
            else if (nlevels.treat <= 2 && isTRUE(mirror)) {
                posneg <- c(1, -1)
                alpha <- .8
                legend.alpha <- alpha
                expandScale <- ggplot2::expansion(mult = .05)
            }
            else {
                mirror <- FALSE
                posneg <- rep(1, nlevels.treat)
                alpha <- .4
                legend.alpha <- alpha/nlevels.treat
                expandScale <- ggplot2::expansion(mult = c(0, .05))
            }
            
            if (type == "histogram") {
                if (isTRUE(disp.means)) {
                    D$var.mean <- ave(D[c("var", "weights")], D[c("treat", facet)], 
                                      FUN = function(x) w.m(x[[1]], x[[2]]))[[1]]
                }
                
                if (!is_(args$bins, "numeric")) args$bins <- 12
                geom_fun <- function(t) {
                    out <- list(ggplot2::geom_histogram(data = D[D$treat == levels(D$treat)[t],],
                                                        mapping = aes(x = .data$var, y = posneg[t]*ggplot2::after_stat(!!sym("count")), 
                                                                      weight = .data$weights,
                                                                      fill = names(colors)[t]),
                                                        alpha = alpha, bins = args$bins, color = "black"),
                                NULL)
                    if (isTRUE(disp.means)) out[[2]] <-
                            ggplot2::geom_segment(data = unique(D[D$treat == levels(D$treat)[t], c("var.mean", facet), drop = FALSE]),
                                                  mapping = aes(x = .data$var.mean, xend = .data$var.mean, y = 0, yend = posneg[t]*Inf), 
                                                  color = if (isTRUE(mirror)) "black" else colors[[t]])
                    return(clear_null(out))
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
                    
                    ggplot2::geom_step(data = D[D$treat == levels(D$treat)[t],],
                                       mapping = aes(x = .data$var, y = .data$cum.pt, color = names(colors)[t]))
                    
                }
                ylab <- "Cumulative Proportion"
            }
            else {
                #Density arguments supplied through ...
                bw <- if_null_then(args$bw, "nrd0")
                adjust <- if_null_then(args$adjust, 1)
                kernel <- if_null_then(args$kernel, "gaussian")
                n <- if_null_then(args$n, 512)
                
                if (is.character(bw)) {
                    t.sizes <- tapply(rep(1, NROW(D)), D$treat, sum)
                    smallest.t <- names(t.sizes)[which.min(t.sizes)]
                    if (is.function(get0(paste0("bw.", bw)))) {
                        bw <- get0(paste0("bw.", bw))(D$var[D$treat == smallest.t])
                    }
                    else {
                        stop(paste(bw, "is not an acceptable entry to 'bw'. See ?stats::density for allowable options."), call. = FALSE)
                    }
                }
                
                if (isTRUE(disp.means)) {
                    D$var.mean <- ave(D[c("var", "weights")], D[c("treat", facet)], 
                                      FUN = function(x) w.m(x[[1]], x[[2]]))[[1]]
                }
                
                geom_fun <- function(t) {
                    out <- list(
                        ggplot2::geom_density(data = D[D$treat == levels(D$treat)[t],],
                                              mapping = aes(x = .data$var, y = posneg[t]*ggplot2::after_stat(!!sym("density")),
                                                            weight = .data$weights, fill = names(colors)[t]),
                                              alpha = alpha, bw = bw, adjust = adjust,
                                              kernel = kernel, n = n, trim = TRUE,
                                              outline.type = "full", stat = StatDensity2),
                        NULL)
                    if (isTRUE(disp.means)) out[[2]] <-
                            ggplot2::geom_segment(data = unique(D[D$treat == levels(D$treat)[t], c("var.mean", facet), drop = FALSE]),
                                                  mapping = aes(x = .data$var.mean, xend = .data$var.mean, y = 0, yend = posneg[t]*Inf), 
                                                  color = if (isTRUE(mirror)) "black" else colors[[t]])
                    return(clear_null(out))
                }
                ylab <- "Density"
                
            }
            
            bp <- Reduce("+", c(list(ggplot2::ggplot()),
                                lapply(seq_len(nlevels.treat), geom_fun))) +
                ggplot2::scale_fill_manual(values = colors, guide = ggplot2::guide_legend(override.aes = list(alpha = legend.alpha, size = .25))) +
                ggplot2::scale_color_manual(values = colors) +
                ggplot2::labs(x = var.name, y = ylab, title = title, subtitle = subtitle,
                              fill = "Treatment", color = "Treatment") + 
                ggplot2::scale_y_continuous(expand = expandScale)
            
            if (isTRUE(mirror)) bp <- bp + ggplot2::geom_hline(yintercept = 0)
        }
    }
    
    #Themes
    bp <- bp + ggplot2::theme(panel.background = element_rect(fill = "white", color = "black"),
                              panel.border = element_rect(fill = NA, color = "black"),
                              plot.background = element_blank(),
                              legend.position = position,
                              legend.key=element_rect(fill = "white", color = "black", size = .25))
    if (!isTRUE(grid)) {
        bp <- bp + ggplot2::theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())
    }
    else {
        bp <- bp + ggplot2::theme(panel.grid.major = element_line(color = "gray87"),
                                  panel.grid.minor = element_line(color = "gray90"))
    }
    
    if (is_not_null(facet)) {
        if (is_not_null(facet.formula)) {
            if (!rlang::is_formula(facet.formula)) {
                stop("'facet.formula' must be a formula.", call. = FALSE)
            }
            test.facet <- invisible(ggplot2::facet_grid(facet.formula))
            if (any(c(names(test.facet$params$rows), names(test.facet$params$cols)) %nin% facet)) {
                stop(paste("Only", word_list(facet, is.are = TRUE, quotes = 2), "allowed in 'facet.formula'."), call. = FALSE)
            }
            if ("which" %nin% c(names(test.facet$params$rows), names(test.facet$params$cols))) {
                if (length(which) > 1) stop("\"which\" must be in the facet formula when the which argument refers to more than one sample.", call. = FALSE)
                else message(paste0("Displaying balance for the ", if (which %in% c("Adjusted Sample", "Unadjusted Sample")) tolower(which) else paste(which, "sample"), "."))
            }
        }
        else if (length(facet) >= 2) {
            if ("which" %in% facet) {
                facet.formula <- reformulate("which", facet[facet %nin% "which"])
            }
            else if ("imp" %in% facet) {
                facet.formula <- reformulate(facet[facet %nin% "imp"], "imp")
            }
            else {
                facets <- data.frame(facet = facet, length = vapply(facet, function(x) nlevels(D[[x]]), numeric(1L)),
                                     stringsAsFactors = FALSE)
                facets <- facets[with(facets, order(length, facet, decreasing = c(FALSE, TRUE))), ]
                facet.formula <- formula(facets)
            }
        }
        else facet.formula <- reformulate(facet, ".")
        bp <- bp + ggplot2::facet_grid(facet.formula, drop = FALSE, scales = ifelse("subclass" %in% facet, "free_x", "fixed"))
    }
    
    return(bp)
}
