love.plot <- function(b, stat = c("mean.diffs", "variance.ratios"), threshold = NULL, abs = FALSE, var.order = NULL, no.missing = TRUE, var.names = NULL, drop.distance = TRUE, cluster.fun = c("mean", "median", "max", "range"), colors = NULL, shapes = NULL, line = FALSE, ...) {
    if (!"bal.tab" %in% class(b)) stop("The first argument must be a bal.tab object, the output of a call to bal.tab().")
    if (any(class(b) == "bal.tab.cont")) stat <- "correlation"
    else stat <- match.arg(stat)
    
    args <- list(...)
    #size
    #shape (deprecated)
    #un.color (deprecated)
    #adj.color (deprecated)
    #title
    #colors
    #shapes
    
    null.threshold <- is.null(threshold)
    if (!null.threshold) {
        if (!is.numeric(threshold) || length(threshold) > 1) stop("threshold must be a single number.", call. = FALSE)
    }
    
    which.stat <- switch(stat, mean.diffs = "Diff", variance.ratios = "V.Ratio", correlation = "Corr")
    which.stat2 <- switch(stat, mean.diffs = "Mean Difference", variance.ratios = "Variance Ratio", correlation = "Correlation")
    Cluster.Fun <- ""
    
    if (any(class(b) == "bal.tab.cluster")) {
        if (is.null(b$print.options$which.cluster) || is.na(b$print.options$which.cluster)) {
            Cluster.Fun <- switch(match.arg(cluster.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Cluster.Fun == "Range") {
                if (b$print.options$quick) stop("cluster.fun = \"range\" cannot be used when the original call to bal.tab() used quick = TRUE.", call. = FALSE)
                title <- paste0("Covariate Balance\n", which.stat2, " Range Across Clusters")
                B <- b[["Cluster.Summary"]][, c("Type",
                                                paste("Max", which.stat, "Adj", sep = "."),
                                                paste("Min", which.stat, "Adj", sep = "."),
                                                paste("Mean", which.stat, "Adj", sep = "."),
                                                paste("Max", which.stat, "Un", sep = "."),
                                                paste("Min", which.stat, "Un", sep = "."),
                                                paste("Mean", which.stat, "Un", sep = "."))]
            }
            else {
                title <- paste0("Covariate Balance\n", Cluster.Fun, " ", which.stat2, " Across Clusters")
                #which.stat <- paste(Cluster.Fun, which.stat, sep = ".")
                B <- b[["Cluster.Summary"]][, c("Type", 
                                                paste(Cluster.Fun, which.stat, "Adj", sep = "."),  
                                                paste(Cluster.Fun, which.stat, "Un", sep = "."))]
            }
            
            abs <- TRUE
        }
        else if (is.character(b$print.options$which.cluster)) {
            cluster.names.good <- sapply(names(b$Cluster.Balance), function(x) any(sapply(b$print.options$which.cluster, function(y) isTRUE(all.equal(x, y)))))
            if (any(cluster.names.good)) {
                if (sum(cluster.names.good) == 1) {
                    B <- b[["Cluster.Balance"]][cluster.names.good][[1]]
                    title <- paste0("Covariate Balance\nCluster: ", b$print.options$which.cluster)
                    
                }
                else {
                    stop("love.plot can only display balance for one cluster at a time. Make sure the argument to which.cluster in bal.tab() is the name or index of a single cluster.", call. = FALSE)
                }
            }
            else {
                stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
            }
        }
        else if (is.numeric(b$print.options$which.cluster)) {
            cluster.numbers.good <- (seq_along(b$Cluster.Balance) %in% b$print.options$which.cluster)
            if (any(cluster.numbers.good)) {
                if (sum(cluster.numbers.good) == 1) {
                    B <- b$Cluster.Balance[cluster.numbers.good][[1]]
                    title <- paste0("Covariate Balance\nCluster: ", names(b$Cluster.Balance)[b$print.options$which.cluster])
                }
                else {
                    stop("love.plot can only display balance for one cluster at a time. Make sure the argument to which.cluster in bal.tab() is the name or index of a single cluster.", call. = FALSE)
                }
            }
            else {
                stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
            }
        }
    }
    else if (any(class(b) == "bal.tab.subclass")) {
        if (which.stat=="V.Ratio") stop("Variance ratios not currently supported for subclassification.", call. = FALSE)
        B <- b[["Balance.Across.Subclass"]]
        title <- "Covariate Balance\nAcross Subclasses"
    }
    else {
        B <- b[["Balance"]]
        title <- "Covariate Balance"
    }
    
    if (drop.distance) B <- B[is.na(match(B[,"Type"], "Distance")),]
    
    if (null.threshold) {
        if (which.stat=="Diff") {
            if (!is.null(b$print.options$m.threshold)) {
                threshold <- b$print.options$m.threshold
                null.threshold <- FALSE
            }
        }
        else if (which.stat=="V.Ratio") {
            if (!is.null(b$print.options$v.threshold)) {
                threshold <- b$print.options$v.threshold
                null.threshold <- FALSE
            }
        }
        else if (which.stat=="Corr") {
            if (!is.null(b$print.options$r.threshold)) {
                threshold <- b$print.options$r.threshold
                null.threshold <- FALSE
            }
        }
    }
    
    var.labels <- row.names(B)
    if (!is.null(var.names)) {
        if (is.data.frame(var.names)) {
            if (ncol(var.names)==1) {
                if (any(sapply(var.labels, function(x) x %in% row.names(var.names)))) {
                    which. <- match(var.labels, row.names(var.names))
                    var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.,1]))
                }
                else if (nrow(var.names)==length(var.labels)) {
                    var.labels <- ifelse(is.na(var.names) | var.names=="", var.labels, as.character(var.names[,1]))
                }
                else warning("var.names has only one column, but the row names do not correspond to variable names in the bal.tab object and the length of var.names differs from the number of variables to rename.", call. = FALSE)
            }
            else if (ncol(var.names)>1) {
                if (ncol(var.names)>2) warning("Only using first 2 columns of var.names", call. = FALSE)
                if (any(sapply(var.labels, function(x) x %in% as.character(var.names[,1])))) {
                    which. <- match(var.labels, var.names[,1])
                    var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.,2]))
                }
                else warning("var.names has more than one column, but the values in the first column do not correspond to variable names in the bal.tab object.", call. = FALSE)
            } 
        }
        else if (is.character(var.names) || is.factor(var.names)) {
            if (!is.null(names(var.names))) {
                if (any(sapply(var.labels, function(x) x %in% names(var.names)))) {
                    which. <- match(var.labels, names(var.names))
                    var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.]))
                }
                else warning("var.names is a named vector but the value names do not correspond to variable names in the bal.tab object and the length of var.names differs from the number of variables to rename.", call. = FALSE)
            }
            else if (length(var.names)==length(var.labels)) {
                var.labels <- ifelse(is.na(var.names) | var.names=="", var.labels, as.character(var.names))
            }
            else warning("var.names is a vector, but its values are unnamed and its length differs from the number of variables to rename.", call. = FALSE)
        }
        else warning("Argument to var.names is not one of the accepted structures and will be ignored.\n  See help(love.plot) for details.", immediate.=TRUE, call. = FALSE)
    }
    
    Sample <- min.stat <- max.stat <- mean.stat <- NULL #To avoid CRAN checks
    if (Cluster.Fun == "Range") {
        SS <- data.frame(var = rep(var.labels, 2), 
                         min.stat = c(B[, paste("Min", which.stat, "Un", sep = ".")], B[, paste("Min", which.stat, "Adj", sep = ".")]),
                         max.stat = c(B[, paste("Max", which.stat, "Un", sep = ".")], B[, paste("Max", which.stat, "Adj", sep = ".")]), 
                         mean.stat = c(B[, paste("Mean", which.stat, "Un", sep = ".")], B[, paste("Mean", which.stat, "Adj", sep = ".")]),
                         Sample=c(rep("Unadjusted", nrow(B)), rep("Adjusted", nrow(B))))
        if (all(is.na(SS[, c("min.stat", "max.stat", "mean.stat")]))) stop("No balance statistics to display.", call. = FALSE)
        if (all(is.na(SS[SS$Sample=="Adjusted", c("min.stat", "max.stat", "mean.stat")]))) {
            gone <- "adjusted"
            SS <- SS[SS[, "Sample"]=="Unadjusted",]
        }
        else if (all(is.na(SS[SS$Sample=="Unadjusted", c("min.stat", "max.stat", "mean.stat")]))) {
            gone <- "unadjusted"
            warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
            SS <- SS[SS[, "Sample"]=="Adjusted",]
        }
        else gone <- ""            
        if (abs) dec <- FALSE
        else dec <- TRUE
        if (!is.null(var.order)) {
            ua <- c("adjusted", "unadjusted")
            var.order <- match.arg(var.order, ua)
            if (var.order == gone) {
                new.var.order <- ua[is.na(match(ua, gone))]
                warning(paste0("var.order was set to \"", var.order, "\", but no ", var.order, " ", tolower(which.stat2), "s were calculated. Using \"", new.var.order, "\" instead."), call. = FALSE, immediate. = TRUE)
                var.order <- new.var.order
            }
            SS[, "var"] <- factor(SS[, "var"], levels=SS[order(SS[tolower(SS[, "Sample"])==var.order, "mean.stat"], decreasing = dec), "var"])
        }
        else SS[, "var"] <- factor(SS[, "var"], levels = unique(SS[, "var"])[order(unique(SS[, "var"]), decreasing = TRUE)])
        SS[, "Sample"] <- factor(SS[, "Sample"], levels = c("Unadjusted", "Adjusted"))
        if (which.stat == "Diff" && any(abs(SS[, "max.stat"]) > 5)) warning("Large mean differences detected; you may not be using standardizied mean differences for continuous variables. To do so, specify continuous=\"std\" in bal.tab().", call.=FALSE, noBreaks.=TRUE)
        if (no.missing) SS <- SS[!is.na(SS[, "min.stat"]),]
        SS$stat <- SS[,"mean.stat"]
    }
    else {
        SS <- data.frame(var=rep(var.labels, 2), 
                         stat=c(B[, ifelse(Cluster.Fun == "", paste(which.stat, "Un", sep = "."),
                                           paste(Cluster.Fun, which.stat, "Un", sep = "."))], 
                                  B[, ifelse(Cluster.Fun == "", paste(which.stat, "Adj", sep = "."),
                                             paste(Cluster.Fun, which.stat, "Adj", sep = "."))]), 
                         Sample=c(rep("Unadjusted", nrow(B)), rep("Adjusted", nrow(B))))
        if (all(is.na(SS[, "stat"]))) stop("No balance statistics to display.", call. = FALSE)
        if (all(is.na(SS[SS$Sample=="Adjusted", "stat"]))) {
            gone <- "adjusted"
            SS <- SS[SS[, "Sample"]=="Unadjusted",]
        }
        else if (all(is.na(SS[SS$Sample=="Unadjusted", "stat"]))) {
            gone <- "unadjusted"
            warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
            SS <- SS[SS[, "Sample"]=="Adjusted",]
        }
        else gone <- ""
        if (abs) {
            SS[, "stat"] <- abs(SS[, "stat"])
            dec <- FALSE}
        else dec <- TRUE
        if (!is.null(var.order)) {
            ua <- c("adjusted", "unadjusted")
            var.order <- match.arg(var.order, ua)
            if (var.order == gone) {
                new.var.order <- ua[is.na(match(ua, gone))]
                warning(paste0("var.order was set to \"", var.order, "\", but no ", var.order, " ", tolower(which.stat2), "s were calculated. Using \"", new.var.order, "\" instead."), call. = FALSE, immediate. = TRUE)
                var.order <- new.var.order
            }
            SS[, "var"] <- factor(SS[, "var"], levels=SS[order(SS[tolower(SS[, "Sample"])==var.order, "stat"], decreasing = dec), "var"])
        }
        else SS[, "var"] <- factor(SS[, "var"], levels = unique(SS[, "var"])[order(unique(SS[, "var"]), decreasing = TRUE)])
        SS[, "Sample"] <- factor(SS[, "Sample"], levels = c("Unadjusted", "Adjusted"))
        if (which.stat == "Diff" && any(abs(SS[, "stat"]) > 5)) warning("Large mean differences detected; you may not be using standardizied mean differences for continuous variables. To do so, specify continuous=\"std\" in bal.tab().", call.=FALSE, noBreaks.=TRUE)
        if (no.missing) SS <- SS[!is.na(SS[, "stat"]),]
    }
    SS <- SS[order(SS[, "var"]),]
    
    #Make the plot
    #library(ggplot2)
    
    #Setting up appearance
    #Size
    if (is.null(args$size)) size <- 1
    else if (is.numeric(args$size[1])) size <- args$size[1]
    else {
        warning("The argument to size must be a number. Using 1 instead.", call. = FALSE)
        size <- 1
    }
    stroke <- .8*size
    nudge <- .15*size
    
    #Shape
    # if (is.null(args$shape)) shape <- 21 #circle
    # else if (is.numeric(args$shape[1]) && any(21:25 == args$shape[1])) shape <- args$shape[1]
    # else {
    #     warning("The argument to shape must be a number between 21 and 25. Using 21 (circle) instead.", call. = FALSE)
    #     shape <- 21
    # }
    
    #Colors
    # null.un.color <- is.null(args$un.color)
    # null.adj.color <- is.null(args$adj.color)
    # if (!null.un.color || !null.adj.color) {
    #     isColor <- function(x) {
    #         tryCatch(is.matrix(col2rgb(x)), 
    #                  error = function(e) FALSE)
    #     }
    # }
    # 
    # if (null.un.color) un.color <- "white"
    # else if (isColor(args$un.color[1])) {
    #     un.color <- args$un.color[1]
    # }
    # else {
    #     warning("The argument to un.color is not a recognized color. Using white instead.", call. = FALSE)
    #     un.color <- "white"
    # }
    # if (null.adj.color) adj.color <- "black"
    # else if (isColor(args$adj.color[1])) {
    #     adj.color <- args$adj.color[1]
    # }
    # else {
    #     warning("The argument to adj.color is not a recognized color. Using black instead.", call. = FALSE)
    #     adj.color <- "black"
    # }
    
    null.colors <- is.null(colors)
    if (!null.colors) {
        isColor <- function(x) {
            tryCatch(is.matrix(col2rgb(x)), 
                     error = function(e) FALSE)
        }
    }
    
    if (null.colors) {
        colors <- c("black", "black")
        fill <- c("white", "white")
        if (length(shapes) == 0) {
            shapes <- c(21, 24)
        }
        else if (!is.numeric(shapes) || !all(shapes %in% 21:25) || length(shapes) > 2) {
            warning("The argument to shape must contain one or two numbers between 21 and 25. \nUsing 21 (circle) and 24 (triangle) instead.", call. = FALSE)
            shapes <- c(21, 24)
        }
        else if (length(shapes) == 1) shapes <- rep(shapes, 2)
        
    }
    else {
        if (length(colors) == 1) colors <- rep(colors, 2)
        else if (length(colors) > 2) {
            colors <- colors[1:2]
            warning("Only using first two values in colors.", call. = FALSE)
        }
        if (all(sapply(colors, isColor))) {
            fill <- colors
            if (length(shapes) == 0) {
                shapes <- c(21, 21)
            }
            else if (!is.numeric(shapes) || !all(shapes %in% 21:25) || length(shapes) > 2) {
                warning("The argument to shape must contain one or two numbers between 21 and 25. \nUsing 21 (circle) for both instead.", call. = FALSE)
                shapes <- c(21, 21)
            }
            else if (length(shapes) == 1) shapes <- rep(shapes, 2)
        }
        else {
            warning("The argument to colors contains at least one value that is not a recognized color.", call. = FALSE)
            colors <- c("black", "black")
            fill <- c("white", "white")
            if (length(shapes) == 0) {
                shapes <- c(21, 24)
            }
            else if (!is.numeric(shapes) || !all(shapes %in% 21:25) || length(shapes) > 2) {
                warning("The argument to shape must contain one or two numbers between 21 and 25. \nUsing 21 (circle) and 24 (triangle) instead.", call. = FALSE)
                shapes <- c(21, 24)
            }
            else if (length(shapes) == 1) shapes <- rep(shapes, 2)
        }
    }
    
    #Title
    if (!is.null(args$title)) title <- as.character(args$title)
    
    lp <- ggplot(data = SS, aes(y = var, x = stat, group = Sample)) + theme(
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = "black")
        ) + 
        scale_shape_manual(values = shapes) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = colors) + 
        labs(title = title, y = "") 
    if (which.stat == "Corr") {
        baseline.xintercept <- 0
        if (abs) {
            lp <- lp + xlab("Absolute Treatment-Covariate Correlations")
            if (!null.threshold) threshold.xintercept <- abs(threshold)
        }
        else {
            lp <- lp + xlab("Treatment-Covariate Correlations") 
            if (!null.threshold) threshold.xintercept <- c(-threshold, threshold)
        }
    }
    else if (which.stat == "Diff") {
        baseline.xintercept <- 0
        if (abs) {
            lp <- lp + xlab("Absolute Mean Differences") 
            if (!null.threshold) threshold.xintercept <- abs(threshold)
        }
        else {
            lp <- lp + xlab("Mean Differences") 
            if (!null.threshold) threshold.xintercept <- c(-threshold, threshold)
        }
    }
    else if (which.stat == "V.Ratio") {
        baseline.xintercept <- 1
        lp <- lp + xlab("Variance Ratios")
        if (!null.threshold) threshold.xintercept <- max(threshold, 1/threshold)
    }
    
    lp <- lp + geom_vline(xintercept = baseline.xintercept, linetype = 1, color = "gray5")
    if (!null.threshold) lp <- lp + geom_vline(xintercept = threshold.xintercept, linetype=2, color = "gray8")
    
    if (Cluster.Fun == "Range") {
        position.nudge <- position_nudge(y = ifelse(SS$Sample == "Adjusted", -nudge, nudge))
        if (args$line == TRUE) {
            lp <- lp + geom_path(aes(color = Sample), size = size*.8, position = position.nudge)
        }
        lp <- lp + 
            geom_segment(aes(y = var, yend = var, x = min.stat, xend = max.stat, color = Sample), position = position.nudge, lineend = "butt", size = size) + 
            geom_point(aes(y = var, x = mean.stat, shape = Sample, color = Sample), fill = "white", size = 2*size, stroke = stroke, position = position.nudge) + 
            labs(title = title, y = "")
    }
    else {
        if (line == TRUE) {
            lp <- lp + geom_path(aes(color = Sample), size = size*.8)
        }
        lp <- lp + geom_point(aes(shape = Sample,
                                    color = Sample), 
                               size = 2*size, stroke = stroke, fill = "white") 
        
    }
    # if (Cluster.Fun == "Range") {
    #     lp <- lp + 
    #         geom_segment(aes(y = var, yend = var, x = min.stat, xend = max.stat), position = position_nudge(y = ifelse(SS$Sample == "Adjusted", -nudge, nudge)), lineend = "butt", size = size) + 
    #         geom_point(aes(y = var, x = mean.stat, fill = Sample, color = Sample), size = 2*size, stroke = stroke, shape = shape, color = "black", position = position_nudge(y = ifelse(SS$Sample == "Adjusted", -nudge, nudge))) + 
    #         scale_fill_manual(values=c(un.color, adj.color)) +
    #         scale_color_manual(values=border.colors) + 
    #         labs(title = title, y = "")
    # }
    # else {lp <- lp + geom_point(aes(y = var, x = stat, fill = Sample, color = Sample), 
    #                             size = 2*size, stroke = stroke, shape = shape) + 
    #     labs(title = title, y = "") + 
    #     scale_fill_manual(values=c(un.color, adj.color)) +
    #     scale_color_manual(values=border.colors)
    # }
    return(lp)
}