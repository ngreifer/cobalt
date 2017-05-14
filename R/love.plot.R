love.plot <- function(b, stat = c("mean.diffs", "variance.ratios"), threshold = NULL, abs = FALSE, var.order = NULL, no.missing = TRUE, var.names = NULL, drop.distance = FALSE, agg.fun = c("mean", "median", "max", "range"), 
                       colors = c("red", "blue"), shapes = NULL, line = FALSE, ...) {
    if (!"bal.tab" %in% class(b)) stop("The first argument must be a bal.tab object, the output of a call to bal.tab().")
    if (any(class(b) == "bal.tab.cont")) stat <- "correlation"
    else stat <- match.arg(stat)
    
    args <- list(...)
    #size
    #shape (deprecated)
    #un.color (deprecated)
    #adj.color (deprecated)
    #title
    #subtitle
    #limits
    #cluster.fun (deprecated)
    
    null.threshold <- is.null(threshold)
    if (!null.threshold) {
        if (!is.numeric(threshold) || length(threshold) > 1) stop("threshold must be a single number.", call. = FALSE)
    }
    
    if (length(args$cluster.fun) > 0) agg.fun <- args$cluster.fun
    which.stat <- switch(stat, mean.diffs = "Diff", variance.ratios = "V.Ratio", correlation = "Corr")
    which.stat2 <- switch(stat, mean.diffs = "Mean Difference", variance.ratios = "Variance Ratio", correlation = "Correlation")
    Agg.Fun <- NULL
    subtitle <- NULL
    title <- "Covariate Balance"
    
    facet <- NULL
    
    cluster.names.good <- NULL
    imp.numbers.good <- NULL
    null.cluster <- is.null(b$print.options$which.cluster) || is.na(b$print.options$which.cluster)
    null.imp <- (is.null(b$print.options$which.imp) || is.na(b$print.options$which.imp))
    
    if (any(class(b) == "bal.tab.imp.cluster")) {
        #Get imp.numbers.good
        imp.numbers <- seq_along(b[["Imputation.Balance"]])
        if (null.imp) {
            imp.numbers.good <- NULL
        }
        else if (is.numeric(b$print.options$which.imp)) {
            imp.numbers.good <- setNames(imp.numbers %in% b$print.options$which.imp, imp.numbers)
        }
        else stop("The argument to which.imp in bal.tab() must be NULL, NA, or the indices of imputations.", call. = FALSE)
        
        #Get cluster.names.good
        cluster.names <- names(b[["Cluster.Balance.Across.Imputations"]])
        if (null.cluster) {
            cluster.names.good <- NULL
        }
        else if (is.character(b$print.options$which.cluster)) {
            cluster.names.good <- sapply(cluster.names, function(x) any(sapply(b$print.options$which.cluster, function(y) identical(x, y))))
        }
        else if (is.numeric(b$print.options$which.cluster)) {
            cluster.names.good <- setNames(seq_along(b$Cluster.Balance) %in% b$print.options$which.cluster, cluster.names)
        }
        else stop("The argument to which.cluster in bal.tab() must be NULL, NA, or the names or indices of clusters.", call. = FALSE)
        
        #Set configuration type of B using which.imp and which.cluster
        if (null.imp) {
            if (null.cluster) {
                config <- "agg.all"
            }
            else { #1, #6
                if (any(cluster.names.good)) {
                    config <- "agg.imp"
                }
                else {
                    stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
                }
                
            }
        }
        else if (sum(imp.numbers.good) == 1) {
            if (null.cluster) {
                config <- "agg.cluster"
            }
            else { 
                if (any(cluster.names.good)) {
                    config <- "agg.none"
                }
                else {
                    stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
                }
            }
        }
        else {
            if (any(imp.numbers.good)) {
                if (null.cluster) {
                    config <- "agg.cluster"
                }
                else if (sum(cluster.names.good) == 1) {
                    config <- "agg.none"
                }
                else {
                    stop("At least one of which.cluster or which.imp must be NULL, NA, or of length 1.", call. = FALSE)
                }
            }
            else {
                stop("Make sure the arguments to which.imp in bal.tab() are valid imputation indices.", call. = FALSE)
            }
        }
        
        #Get B from b based on configuration
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(b[["Imputation.Balance"]])[imp.numbers.good],
                                         function(x) do.call("rbind", lapply(names(b[["Imputation.Balance"]][[x]][["Cluster.Balance"]])[cluster.names.good],
                                                                             function(y) cbind(b[["Imputation.Balance"]][[x]][["Cluster.Balance"]][[y]][["Balance.Table"]],
                                                                                               cluster = y,
                                                                                               imp = paste("Imputation:", x),
                                                                                               var.names = rownames(b[["Imputation.Balance"]][[x]][["Cluster.Balance"]][[y]][["Balance.Table"]]))))))
            if (sum(imp.numbers.good) == 1) {
                facet <- "cluster"
                subtitle <- paste("Imputation:", imp.numbers[imp.numbers.good])
            }
            else {
                facet <- "imp"
                subtitle <- paste("Cluster:", cluster.names[cluster.names.good])
            }
        }
        else if (config == "agg.imp") {
            Agg.Fun <- switch(match.arg(agg.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                if (b$print.options$quick) stop("\"range\" cannot be used when the original call to bal.tab() used quick = TRUE.", call. = FALSE)
                subtitle <- paste0(which.stat2, " Range Across Imputations")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Imputations", sep = " ")
            }
            B <- do.call("rbind", lapply(names(b[["Cluster.Balance.Across.Imputations"]])[cluster.names.good], 
                                         function(x) cbind(b[["Cluster.Balance.Across.Imputations"]][[x]][["Cluster.Balance"]], 
                                                           cluster = x, 
                                                           var.names = rownames(b[["Cluster.Balance.Across.Imputations"]][[x]][["Cluster.Balance"]]))))
            abs <- TRUE
            facet <- "cluster"
        }
        else if (config == "agg.cluster") {
            Agg.Fun <- switch(match.arg(agg.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                if (b$print.options$quick) stop("\"range\" cannot be used when the original call to bal.tab() used quick = TRUE.", call. = FALSE)
                subtitle <- paste0(which.stat2, " Range Across Clusters")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Clusters", sep = " ")
            }
            B <- do.call("rbind", lapply(names(b[["Imputation.Balance"]])[imp.numbers.good], 
                                         function(x) cbind(b[["Imputation.Balance"]][[x]][["Cluster.Summary"]], 
                                                           imp = paste("Imputation:", x), 
                                                           var.names = rownames(b[["Imputation.Balance"]][[x]][["Cluster.Summary"]]))))
            abs <- TRUE
            facet <- "imp"
        }
        else if (config == "agg.all") {
            #Cluster.Fun <- switch(match.arg(cluster.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            Agg.Fun <- switch(match.arg(agg.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                if (b$print.options$quick) stop("\"range\" cannot be used when the original call to bal.tab() used quick = TRUE.", call. = FALSE)
                subtitle <- paste0(which.stat2, " Range Across Clusters and Imputations")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Clusters and Imputations", sep = " ")
            }
            B <- cbind(b[["Balance.Across.Imputations"]],
                       var.names = row.names(b[["Balance.Across.Imputations"]]))
            abs <- TRUE
            facet <- NULL
        }
    }
    else if (any(class(b) == "bal.tab.imp")) {
        #Get imp.numbers.good
        imp.numbers <- seq_along(b[["Imputation.Balance"]])
        if (null.imp) {
            config <- "agg.imp"
            imp.numbers.good <- NULL
        }
        else if (is.numeric(b$print.options$which.imp)) {
            imp.numbers.good <- setNames(imp.numbers %in% b$print.options$which.imp, imp.numbers)
            if (any(imp.numbers.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.imp in bal.tab() are valid imputation indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.imp in bal.tab() must be NULL, NA, or the indices of imputations.", call. = FALSE)
        
        #Get B from b
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(b[["Imputation.Balance"]])[imp.numbers.good], 
                                         function(x) cbind(b[["Imputation.Balance"]][[x]][["Balance"]],
                                                           imp = paste("Imputation:", x),
                                                           var.names = rownames(b[["Imputation.Balance"]][[x]][["Balance"]]))))
            facet <- "imp"
        }
        else if (config == "agg.imp") {
            Agg.Fun <- switch(match.arg(agg.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                if (b$print.options$quick) stop("\"range\" cannot be used when the original call to bal.tab() used quick = TRUE.", call. = FALSE)
                subtitle <- paste0(which.stat2, " Range Across Imputations")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Imputations")
            }
            B <- cbind(b[["Balance.Across.Imputations"]],
                       var.names = rownames(b[["Balance.Across.Imputations"]]))
            abs <- TRUE
        }
    }
    else if (any(class(b) == "bal.tab.cluster")) {
        #Get cluster.names.good
        cluster.names <- names(b[["Cluster.Balance"]])
        if (null.cluster) {
            config <- "agg.cluster"
            cluster.names.good <- NULL
        }
        else if (is.character(b$print.options$which.cluster)) {
            cluster.names.good <- sapply(cluster.names, function(x) any(sapply(b$print.options$which.cluster, function(y) identical(x, y))))
            if (any(cluster.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster in bal.tab() are valid cluster names or indices.", call. = FALSE)
            }
        }
        else if (is.numeric(b$print.options$which.cluster)) {
            cluster.names.good <- setNames(seq_along(b$Cluster.Balance) %in% b$print.options$which.cluster, cluster.names)
            if (any(cluster.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster in bal.tab() are valid cluster names or indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.cluster in bal.tab() must be NULL, NA, or the names or indices of clusters.", call. = FALSE)
        
        
        #Get B from b
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(b[["Cluster.Balance"]])[cluster.names.good], 
                                         function(x) cbind(b[["Cluster.Balance"]][[x]][["Balance.Table"]],
                                                           cluster = x,
                                                           var.names = rownames(b[["Cluster.Balance"]][[x]][["Balance.Table"]]))))
            facet <- "cluster"
        }
        else if (config == "agg.cluster") {
            Agg.Fun <- switch(match.arg(agg.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                if (b$print.options$quick) stop("\"range\" cannot be used when the original call to bal.tab() used quick = TRUE.", call. = FALSE)
                subtitle <- paste0(which.stat2, " Range Across Clusters")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Clusters")
            }
            B <- cbind(b[["Cluster.Summary"]],
                       var.names = rownames(b[["Cluster.Summary"]]))
            abs <- TRUE
        }
    }
    else if (any(class(b) == "bal.tab.subclass")) {
        if (which.stat=="V.Ratio") stop("Variance ratios not currently supported for subclassification.", call. = FALSE)
        B <- cbind(b[["Balance.Across.Subclass"]], var.names = row.names(b[["Balance.Across.Subclass"]]))
        subtitle <- "Across Subclasses"
    }
    else {
        B <- cbind(b[["Balance"]], var.names = row.names(b[["Balance"]]))
    }
    
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
    
    if (!is.null(var.names)) {
        if (is.data.frame(var.names)) {
            if (ncol(var.names)==1) {
                if (length(row.names(var.names)) > 0) {
                    new.labels <- setNames(as.list(as.character(rownames(var.names))), var.names[,1])
                }
                else warning("var.names is a data.frame, but its rows are unnamed.", call. = FALSE)
                
            }
            else if (ncol(var.names)>1) {
                if (all(c("old", "new") %in% names(var.names))) {
                    new.labels <- setNames(as.list(as.character(var.names[,"old"])), var.names[,"new"])
                }
                else {
                    if (ncol(var.names)>2) warning("Only using first 2 columns of var.names", call. = FALSE)
                    new.labels <- setNames(as.list(as.character(var.names[,1])), var.names[,2])
                }
            } 
        }
        else if (is.character(var.names) || is.factor(var.names)) {
            if (length(names(var.names))>0) {
                new.labels <- setNames(as.list(names(var.names)), var.names)
            }
            else warning("var.names is a vector, but its values are unnamed.", call. = FALSE)
        }
        else if (is.list(var.names)) {
            if (all(sapply(var.names, function(x) is.character(x) || is.factor(x)))) {
                if (length(names(var.names))>0) {
                    new.labels <- setNames(names(var.names), var.names)
                }
                else warning("var.names is a list, but its values are unnamed.", call. = FALSE)
            }
            else warning("var.names is a list, but its values are not the new names of the variables.", call. = FALSE)
        }
        else warning("Argument to var.names is not one of the accepted structures and will be ignored.\n  See help(love.plot) for details.", immediate.=TRUE, call. = FALSE)
        new.labels <- new.labels[new.labels %in% B[, "var.names"]]
        B[, "var.names"] <- do.call(f.recode, c(list(B[, "var.names"]), as.list(new.labels)))
    }
    
    distance.names <- as.character(unique(B[B[,"Type"] == "Distance",  "var.names"]))
    if (drop.distance) {
        B <- B[is.na(match(B[, "var.names"], distance.names)),]
    }
    
    if (length(var.order)>0) {
        ua <- c("adjusted", "unadjusted")
        var.order <- match.arg(var.order, c(ua, "alphabetical"))
    }
    
    if (length(facet) > 0) {
        if (length(var.order) > 0 && var.order != "alphabetical" && (sum(cluster.names.good) > 1 || sum(imp.numbers.good) > 1)) {
            warning("var.order cannot be set with multiple plots (unless \"alphabetical\"). Ignoring var.order.", call. = FALSE)
            var.order <- NULL
        }
    }
    
    Sample <- min.stat <- max.stat <- mean.stat <- NULL #To avoid CRAN checks
    if (length(Agg.Fun) > 0 && Agg.Fun == "Range") {
        agg.range <- TRUE
    }
    else agg.range <- FALSE
    if (agg.range) {
        SS <- data.frame(var = rep(B[,"var.names"], 2), 
                         min.stat = c(B[, paste("Min", which.stat, "Un", sep = ".")], B[, paste("Min", which.stat, "Adj", sep = ".")]),
                         max.stat = c(B[, paste("Max", which.stat, "Un", sep = ".")], B[, paste("Max", which.stat, "Adj", sep = ".")]), 
                         mean.stat = c(B[, paste("Mean", which.stat, "Un", sep = ".")], B[, paste("Mean", which.stat, "Adj", sep = ".")]),
                         Sample=c(rep("Unadjusted", nrow(B)), rep("Adjusted", nrow(B))))
        if (length(facet) > 0 && "cluster" %in% facet) {
            SS$cluster <- rep(B[, "cluster"], 2)
        }
        if (length(facet) > 0 && "imp" %in% facet) {
            SS$imp <- rep(B[, "imp"], 2)
        }
        
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
        
        if (length(var.order)>0) {
            if (var.order %in% ua) {
                if (var.order == gone) {
                    new.var.order <- ua[is.na(match(ua, gone))]
                    warning(paste0("var.order was set to \"", var.order, "\", but no ", var.order, " ", tolower(which.stat2), "s were calculated. Using \"", new.var.order, "\" instead."), call. = FALSE, immediate. = TRUE)
                    var.order <- new.var.order
                }
                v <- as.character(SS[order(SS[tolower(SS[, "Sample"])==var.order, "mean.stat"], decreasing = dec), "var"])
                SS[, "var"] <- factor(SS[, "var"], 
                                      levels=c(v[is.na(match(v, distance.names))], 
                                               sort(distance.names, decreasing = TRUE)))
            }
            else if (var.order == "alphabetical") {
                SS[, "var"] <- factor(SS[, "var"], levels = c(sort(as.character(unique(SS[is.na(match(SS$var, distance.names)), "var"])), decreasing = TRUE), sort(distance.names, decreasing = TRUE)))
            }
        }
        else {
            SS[, "var"] <- factor(SS[, "var"], levels = c(as.character(unique(SS[is.na(match(SS$var, distance.names)), "var"])[order(unique(SS[is.na(match(SS$var, distance.names)), "var"]), decreasing = TRUE)]), sort(distance.names, decreasing = TRUE)))
        }
        SS[, "Sample"] <- factor(SS[, "Sample"], levels = c("Unadjusted", "Adjusted"))
        if (which.stat == "Diff" && any(abs(SS[, "max.stat"]) > 5)) warning("Large mean differences detected; you may not be using standardizied mean differences for continuous variables. To do so, specify continuous=\"std\" in bal.tab().", call.=FALSE, noBreaks.=TRUE)
        if (no.missing) SS <- SS[!is.na(SS[, "min.stat"]),]
        SS$stat <- SS[,"mean.stat"]
    }
    else {
        SS <- data.frame(var=rep(B[, "var.names"], 2), 
                         stat=c(B[, ifelse(length(Agg.Fun) == 0, paste(which.stat, "Un", sep = "."),
                                           paste(Agg.Fun, which.stat, "Un", sep = "."))], 
                                B[, ifelse(length(Agg.Fun) == 0, paste(which.stat, "Adj", sep = "."),
                                           paste(Agg.Fun, which.stat, "Adj", sep = "."))]), 
                         Sample=c(rep("Unadjusted", nrow(B)), rep("Adjusted", nrow(B))))
        if (length(facet) > 0 && "cluster" %in% facet) {
            SS$cluster <- rep(B[, "cluster"], 2)
        }
        if (length(facet) > 0 && "imp" %in% facet) {
            SS$imp <- rep(B[, "imp"], 2)
        }
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
        
        if (length(var.order)>0) {
            if (var.order %in% ua) {
                if (var.order == gone) {
                    new.var.order <- ua[is.na(match(ua, gone))]
                    warning(paste0("var.order was set to \"", var.order, "\", but no ", var.order, " ", tolower(which.stat2), "s were calculated. Using \"", new.var.order, "\" instead."), call. = FALSE, immediate. = TRUE)
                    var.order <- new.var.order
                }
                v <- as.character(SS[order(SS[tolower(SS[, "Sample"])==var.order, "stat"], decreasing = dec), "var"])
                SS[, "var"] <- factor(SS[, "var"], 
                                      levels=c(v[is.na(match(v, distance.names))], 
                                               sort(distance.names, decreasing = TRUE)))
            }
            else if (var.order == "alphabetical") {
                SS[, "var"] <- factor(SS[, "var"], levels = c(sort(as.character(unique(SS[is.na(match(SS$var, distance.names)), "var"])), decreasing = TRUE), sort(distance.names, decreasing = TRUE)))
            }
        }
        else {
            SS[, "var"] <- factor(SS[, "var"], levels = c(as.character(unique(SS[is.na(match(SS$var, distance.names)), "var"])[order(unique(SS[is.na(match(SS$var, distance.names)), "var"]), decreasing = TRUE)]), sort(distance.names, decreasing = TRUE)))
        }
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
    
    #Shape
    # if (is.null(args$shape)) shape <- 21 #circle
    # else if (is.numeric(args$shape[1]) && any(21:25 == args$shape[1])) shape <- args$shape[1]
    # else {
    #     warning("The argument to shape must be a number between 21 and 25. Using 21 (circle) instead.", call. = FALSE)
    #     shape <- 21
    # }
    
    if (is.null(colors)) {
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
    if (!is.null(args$subtitle)) subtitle <- as.character(args$subtitle)
    #else subtitle <- NULL
    
    if (which.stat == "Corr") {
        baseline.xintercept <- 0
        if (abs) {
            xlab <- "Absolute Treatment-Covariate Correlations"
            if (!null.threshold) threshold.xintercept <- abs(threshold)
        }
        else {
            xlab <- "Treatment-Covariate Correlations"
            if (!null.threshold) threshold.xintercept <- c(-threshold, threshold)
        }
    }
    else if (which.stat == "Diff") {
        baseline.xintercept <- 0
        if (abs) {
            xlab <- "Absolute Mean Differences"
            if (!null.threshold) threshold.xintercept <- abs(threshold)
        }
        else {
            xlab <- "Mean Differences"
            if (!null.threshold) threshold.xintercept <- c(-threshold, threshold)
        }
    }
    else if (which.stat == "V.Ratio") {
        baseline.xintercept <- 1
        xlab <- "Variance Ratios"
        if (!null.threshold) threshold.xintercept <- max(threshold, 1/threshold)
    }
    
    apply.limits <- FALSE
    if (length(args$limits) > 0) {
        limits <- args$limits
        if (length(limits) != 2 || !is.numeric(limits)) {
            warning("limits must be a numeric vector of length 2. Ignoring limits.", call. = FALSE)
        }
        else {
            if (limits[2] < limits[1]) {
                warning("The values in limits must be in ascending order. Reversing them.", call. = FALSE)
                limits <- sort(limits)
            }
            
            if (limits[1] >= 0) limits[1] <- baseline.xintercept - .05*limits[2]
            if (limits[2] <= 0) limits[2] <- baseline.xintercept - .05*limits[1]
            
            if (agg.range) {
                if (any(SS[,"min.stat"] < limits[1]) || any(SS[, "max.stat"] > limits[2])) {
                    for (i in c("min.stat", "stat", "max.stat")) {
                        SS[SS[, i] < limits[1], i] <- limits[1]
                        SS[SS[, i] > limits[2], i] <- limits[2]
                    }
                    warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
            }
            else {
                if (any(SS[,"stat"] < limits[1]) || any(SS[, "stat"] > limits[2])) {
                    warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
            }
            apply.limits <- TRUE
        }
    }
    

    lp <- ggplot(data = SS, aes(y = var, x = stat, group = Sample)) + 
        theme(panel.grid.major = element_line(color = "gray87"),
              panel.grid.minor = element_line(color = "gray90"),
              panel.background = element_rect(fill = "white", color = "black")
        ) + 
        scale_shape_manual(values = shapes) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = colors) + 
        labs(title = title, subtitle = subtitle, y = "", x = xlab) 

    lp <- lp + geom_vline(xintercept = baseline.xintercept, linetype = 1, color = "gray5")
    if (!null.threshold) lp <- lp + geom_vline(xintercept = threshold.xintercept, linetype=2, color = "gray8")
    
    if (agg.range) {
        position.dodge <- ggstance::position_dodgev(.5*(size))
        if (line == TRUE) { #Add line except to distance
            f <- function(x) {x[x$var %in% distance.names, "stat"] <- NA; x}
            lp <- lp + layer(geom = "path", data = f, position = position.dodge, stat = "identity", 
                             aes(color = Sample), params = list(size = size*.8, na.rm = TRUE))
        }
        lp <- lp + 
            ggstance::geom_linerangeh(aes(y = var, xmin = min.stat, xmax = max.stat, 
                                          color = Sample), position = position.dodge, size = size) + 
            geom_point(aes(y = var, x = mean.stat, shape = Sample, color = Sample), 
                       fill = "white", size = 2*size, stroke = stroke, na.rm = TRUE,
                       position = position.dodge) + 
            labs(title = title, y = "")
    }
    else {
        if (line == TRUE) { #Add line except to distance
            f <- function(x) {x[x$var %in% distance.names, "stat"] <- NA; x}
            lp <- lp + layer(geom = "path", data = f, position = "identity", stat = "identity", aes(color = Sample), params = list(size = size*.8, na.rm = TRUE))
        }
        lp <- lp + geom_point(data = SS, aes(shape = Sample,
                                             color = Sample), 
                              size = 2*size, stroke = stroke, fill = "white", na.rm = TRUE) 
        
    }

    if (!drop.distance && length(distance.names) > 0) {
        lp <- lp + geom_hline(linetype = 1, color = "black", 
                              yintercept = nlevels(SS[,"var"]) - length(distance.names) + .5)
    }
    if (apply.limits) {
        lp <- lp + scale_x_continuous(limits = limits, expand = c(0, 0))
    }
    if (length(facet) > 0) {
        lp <- lp + facet_grid(f.build(".", facet), drop = FALSE)
    }

    return(lp)
}