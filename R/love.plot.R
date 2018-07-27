love.plot <- function(x, stat = c("mean.diffs", "variance.ratios", "ks.statistics"), threshold = NULL, 
                      abs = TRUE, var.order = NULL, no.missing = TRUE, var.names = NULL, 
                      drop.distance = FALSE, agg.fun = c("mean", "median", "max", "range"), 
                      colors = NULL, shapes = NULL, line = FALSE, ...) {
    b <- x; rm(x)
    if ("bal.tab" %nin% class(b)) stop("The first argument must be a bal.tab object, the output of a call to bal.tab().")
    if (any(class(b) == "bal.tab.cont")) stat <- "correlation"
    else stat <- match.arg(stat)
    
    args <- list(...)
    #size
    #shape (deprecated)
    #un.color (deprecated)
    #adj.color (deprecated)
    #title
    #subtitle
    #sample.names
    #limits
    #cluster.fun (deprecated)
    
    p.ops <- c("which.cluster", "which.imp", "which.treat", "which.time", "disp.subclass")
    for (i in p.ops) {
        if (i %in% names(args)) b$print.options[[i]] <- args[[i]]
    }
    
    null.threshold <- is_null(threshold)
    if (!null.threshold) {
        if (!is.numeric(threshold) || length(threshold) > 1) stop("threshold must be a single number.", call. = FALSE)
    }
    
    if (is_not_null(args$cluster.fun)) agg.fun <- args$cluster.fun
    which.stat <- switch(stat, mean.diffs = "Diff", variance.ratios = "V.Ratio", ks.statistics = "KS", correlation = "Corr")
    which.stat2 <- switch(stat, mean.diffs = "Mean Difference", variance.ratios = "Variance Ratio", ks.statistics = "Kolmogorov-Smirnov Statistic", correlation = "Correlation")
    Agg.Fun <- NULL
    subtitle <- NULL
    title <- "Covariate Balance"
    
    facet <- NULL
    
    cluster.names.good <- NULL
    imp.numbers.good <- NULL
    treat.names.good <- NULL; disp.treat.pairs <- NULL
    time.names.good <- NULL
    subclass.names <- NULL
    
    #NA = aggregate, NULL = individual
    null.cluster <- is_null(b$print.options$which.cluster)
    na.cluster <- !null.cluster && is.na(b$print.options$which.cluster)
    null.imp <- is_null(b$print.options$which.imp)
    na.imp <- !null.imp && is.na(b$print.options$which.imp)
    null.treat <- is_null(b$print.options$which.treat)
    na.treat <- !null.treat && is.na(b$print.options$which.treat)
    null.time <- is_null(b$print.options$which.time)
    na.time <- !null.time && is.na(b$print.options$which.time)
    
    config <- "agg.none"
    
    if (any(class(b) == "bal.tab.msm")) {
        #Get time.names.good
        time.names <- names(b[["Time.Balance"]])
        if (na.time) {
            config <- "agg.time"
            time.names.good <- NULL
        }
        else if (null.time) {
            config <- "agg.none"
            time.names.good <- setNames(rep(TRUE, length(time.names)), time.names)
        }
        else if (is.character(b$print.options$which.time)) {
            time.names.good <- sapply(time.names, function(x) any(sapply(b$print.options$which.time, function(y) identical(x, y))))
            if (any(time.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.time are valid treatment names or time point indices.", call. = FALSE)
            }
        }
        else if (is.numeric(b$print.options$which.time)) {
            time.names.good <- setNames(seq_along(b$Time.Balance) %in% b$print.options$which.time, time.names)
            if (any(time.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.time are valid treatment names or time point indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.time must be NULL, NA, or the names of treatments or indices of time points.", call. = FALSE)
        
        
        #Get B from b
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(b[["Time.Balance"]])[time.names.good], 
                                         function(x) cbind(b[["Time.Balance"]][[x]][["Balance"]],
                                                           time = x,
                                                           variable.names = rownames(b[["Time.Balance"]][[x]][["Balance"]]))))
            facet <- "time"
        }
        else if (config == "agg.time") {
            if (is_null(b[["Balance.Across.Times"]])) {
                stop("Cannot aggregate across time periods without a balance summary across time periods.\nThis may be because multinomial treatments were used, multiple treatment types were used,\n or quick was set to TRUE and msm.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            #Agg.Fun <- switch(match.arg(agg.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            Agg.Fun <- "Max"
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Time Points")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Time Points")
            }
            B <- cbind(b[["Balance.Across.Times"]],
                       variable.names = rownames(b[["Balance.Across.Times"]]))
        }
    }
    else if (any(class(b) == "bal.tab.imp.cluster")) {
        #Get imp.numbers.good
        imp.numbers <- seq_along(b[["Imputation.Balance"]])
        if (na.imp) {
            imp.numbers.good <- NULL
        }
        else if (null.imp) {
            imp.numbers.good <- setNames(rep(TRUE, length(imp.numbers)), imp.numbers)
        }
        else if (is.numeric(b$print.options$which.imp)) {
            imp.numbers.good <- setNames(imp.numbers %in% b$print.options$which.imp, imp.numbers)
        }
        else stop("The argument to which.imp must be NULL, NA, or the indices of imputations.", call. = FALSE)
        
        #Get cluster.names.good
        cluster.names <- names(b[["Cluster.Balance.Across.Imputations"]])
        if (na.cluster) {
            cluster.names.good <- NULL
        }
        else if (null.cluster) {
            cluster.names.good <- setNames(rep(TRUE, length(cluster.names)), cluster.names)
        }
        else if (is.character(b$print.options$which.cluster)) {
            cluster.names.good <- sapply(cluster.names, function(x) any(sapply(b$print.options$which.cluster, function(y) identical(x, y))))
        }
        else if (is.numeric(b$print.options$which.cluster)) {
            cluster.names.good <- setNames(seq_along(b$Cluster.Balance) %in% b$print.options$which.cluster, cluster.names)
        }
        else stop("The argument to which.cluster must be NULL, NA, or the names or indices of clusters.", call. = FALSE)
        
        #Set configuration type of B using which.imp and which.cluster
        if (na.imp) { #aggregate over all imps
            if (na.cluster) {
                config <- "agg.all"
            }
            else { #1, #6
                if (any(cluster.names.good)) {
                    config <- "agg.imp"
                }
                else {
                    stop("Make sure the arguments to which.cluster are valid names or indices of clusters.", call. = FALSE)
                }
                
            }
        }
        else if (sum(imp.numbers.good) == 1) {
            if (na.cluster) {
                config <- "agg.cluster"
            }
            else { 
                if (any(cluster.names.good)) {
                    config <- "agg.none"
                }
                else {
                    stop("Make sure the arguments to which.cluster are valid names or indices of clusters.", call. = FALSE)
                }
            }
        }
        else if (sum(imp.numbers.good) > 1) {
                if (na.cluster) {
                    config <- "agg.cluster"
                }
                else if (sum(cluster.names.good) == 1) {
                    config <- "agg.none"
                }
                else {
                    stop("At least one of which.cluster or which.imp must be NA or of length 1.", call. = FALSE)
                }
            }
        else {
            stop("Make sure the arguments to which.imp are valid imputation indices.", call. = FALSE)
        }
        
        #Get B from b based on configuration
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(b[["Imputation.Balance"]])[imp.numbers.good],
                                         function(x) do.call("rbind", lapply(names(b[["Imputation.Balance"]][[x]][["Cluster.Balance"]])[cluster.names.good],
                                                                             function(y) cbind(b[["Imputation.Balance"]][[x]][["Cluster.Balance"]][[y]][["Balance"]],
                                                                                               cluster = y,
                                                                                               imp = paste("Imputation:", x),
                                                                                               variable.names = rownames(b[["Imputation.Balance"]][[x]][["Cluster.Balance"]][[y]][["Balance"]]))))))
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
            if (is_null(b[["Cluster.Balance.Across.Imputations"]])) {
                stop("Cannot aggregate across imputations without a balance summary across imputations.\nThis may be because quick was set to TRUE and imp.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            Agg.Fun <- switch(tolower(match.arg(agg.fun)), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Median") {
                warning("The median is being deprecated. Using the mean instead.", call. = FALSE)
                Agg.Fun <- "Mean"
            }
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Imputations")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Imputations", sep = " ")
            }
            B <- do.call("rbind", lapply(names(b[["Cluster.Balance.Across.Imputations"]])[cluster.names.good], 
                                         function(x) cbind(b[["Cluster.Balance.Across.Imputations"]][[x]][["Cluster.Balance"]], 
                                                           cluster = x, 
                                                           variable.names = rownames(b[["Cluster.Balance.Across.Imputations"]][[x]][["Cluster.Balance"]]))))
            facet <- "cluster"
        }
        else if (config == "agg.cluster") {
            if (is_null(b[["Imputation.Balance"]][[1]][["Cluster.Summary"]])) {
                stop("Cannot aggregate across clusters without a balance summary across clusters.\nThis may be because quick was set to TRUE and cluster.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            Agg.Fun <- switch(tolower(match.arg(agg.fun)), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Median") {
                warning("agg.fun = \"median\" is being deprecated. Using \"mean\" instead.", call. = FALSE)
                Agg.Fun <- "Mean"
            }
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Clusters")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Clusters", sep = " ")
            }
            B <- do.call("rbind", lapply(names(b[["Imputation.Balance"]])[imp.numbers.good], 
                                         function(x) cbind(b[["Imputation.Balance"]][[x]][["Cluster.Summary"]], 
                                                           imp = paste("Imputation:", x), 
                                                           variable.names = rownames(b[["Imputation.Balance"]][[x]][["Cluster.Summary"]]))))
            facet <- "imp"
        }
        else if (config == "agg.all") {
            if (is_null(b[["Balance.Across.Imputations"]])) {
                stop("Cannot aggregate across imputations without a balance summary across imputations.\nThis may be because quick was set to TRUE and cluster.summary or imp.summary were set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            #Cluster.Fun <- switch(match.arg(cluster.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            Agg.Fun <- switch(tolower(match.arg(agg.fun)), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Median") {
                warning("The median is being deprecated. Using the mean instead.", call. = FALSE)
                Agg.Fun <- "Mean"
            }
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Clusters and Imputations")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Clusters and Imputations", sep = " ")
            }
            B <- cbind(b[["Balance.Across.Imputations"]],
                       variable.names = row.names(b[["Balance.Across.Imputations"]]))
            facet <- NULL
        }
    }
    else if (any(class(b) == "bal.tab.imp")) {
        #Get imp.numbers.good
        imp.numbers <- seq_along(b[["Imputation.Balance"]])
        if (na.imp) {
            config <- "agg.imp"
            imp.numbers.good <- NULL
        }
        else if (null.imp) {
            config <- "agg.none"
            imp.numbers.good <- setNames(rep(TRUE, length(imp.numbers)), imp.numbers)
        }
        else if (is.numeric(b$print.options$which.imp)) {
            imp.numbers.good <- setNames(imp.numbers %in% b$print.options$which.imp, imp.numbers)
            if (any(imp.numbers.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.imp are valid imputation indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.imp must be NULL, NA, or the indices of imputations.", call. = FALSE)
        
        #Get B from b
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(b[["Imputation.Balance"]])[imp.numbers.good], 
                                         function(x) cbind(b[["Imputation.Balance"]][[x]][["Balance"]],
                                                           imp = paste("Imputation:", x),
                                                           variable.names = rownames(b[["Imputation.Balance"]][[x]][["Balance"]]))))
            facet <- "imp"
        }
        else if (config == "agg.imp") {
            if (is_null(b[["Balance.Across.Imputations"]])) {
                stop("Cannot aggregate across imputations without a balance summary across imputations.\nThis may be because quick was set to TRUE and imp.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            Agg.Fun <- switch(tolower(match.arg(agg.fun)), mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Median") {
                warning("The median is being deprecated. Using the mean instead.", call. = FALSE)
                Agg.Fun <- "Mean"
            }
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Imputations")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Imputations")
            }
            B <- cbind(b[["Balance.Across.Imputations"]],
                       variable.names = rownames(b[["Balance.Across.Imputations"]]))
        }
    }
    else if (any(class(b) == "bal.tab.cluster")) {
        #Get cluster.names.good
        cluster.names <- names(b[["Cluster.Balance"]])
        if (na.cluster) {
            config <- "agg.cluster"
            cluster.names.good <- NULL
        }
        else if (null.cluster) {
            config <- "agg.none"
            cluster.names.good <- setNames(rep(TRUE, length(cluster.names)), cluster.names)
        }
        else if (is.character(b$print.options$which.cluster)) {
            cluster.names.good <- sapply(cluster.names, function(x) any(sapply(b$print.options$which.cluster, function(y) identical(x, y))))
            if (any(cluster.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
            }
        }
        else if (is.numeric(b$print.options$which.cluster)) {
            cluster.names.good <- setNames(seq_along(b$Cluster.Balance) %in% b$print.options$which.cluster, cluster.names)
            if (any(cluster.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.cluster must be NULL, NA, or the names or indices of clusters.", call. = FALSE)
        
        
        #Get B from b
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(b[["Cluster.Balance"]])[cluster.names.good], 
                                         function(x) cbind(b[["Cluster.Balance"]][[x]][["Balance"]],
                                                           cluster = x,
                                                           variable.names = rownames(b[["Cluster.Balance"]][[x]][["Balance"]]))))
            facet <- "cluster"
        }
        else if (config == "agg.cluster") {
            if (is_null(b[["Cluster.Summary"]])) {
                stop("Cannot aggregate across clusters without a balance summary across clusters.\nThis may be because quick was set to TRUE and cluster.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            
            tryCatch({agg.fun <- tolower(match.arg(agg.fun))}, 
                     error = function(e) stop("agg.fun should be one of \"mean\", \"max\", or \"range\".", call. = FALSE))
            Agg.Fun <- switch(agg.fun, mean = "Mean", median = "Median", max = "Max", range = "Range")
            if (Agg.Fun == "Median") {
                warning("The median is deprecated. Using the mean instead.", call. = FALSE)
                Agg.Fun <- "Mean"
            }
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Clusters")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Clusters")
            }
            B <- cbind(b[["Cluster.Summary"]],
                       variable.names = rownames(b[["Cluster.Summary"]]))
        }
    }
    else if (any(class(b) == "bal.tab.multi")) {
        #Get cluster.names.good
        treat.names <- b$print.options$treat.names
        if (na.treat) {
            config <- "agg.pair"
            treat.names.good <- NULL
        }
        else if (null.treat) {
            config <- "agg.none"
            treat.names.good <- rep(TRUE, length(treat.names))
        }
        else if (is.character(b$print.options$which.treat)) {
            treat.names.good <- treat.names %in% b$print.options$which.treat
            if (any(treat.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.treat are valid treatment names or indices.", call. = FALSE)
            }
        }
        else if (is.numeric(b$print.options$which.treat)) {
            treat.names.good <- seq_along(treat.names) %in% b$print.options$which.treat
            if (any(treat.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.cluster must be NULL, NA, or the names or indices of clusters.", call. = FALSE)
        
        
        if (na.treat) {
            disp.treat.pairs <- character(0)
        }
        else {
            if (b$print.options$pairwise) {
                if (sum(treat.names.good) == 1) {
                    disp.treat.pairs <- names(b[["Pair.Balance"]])[sapply(names(b[["Pair.Balance"]]), function(x) any(b[["Pair.Balance"]][[x]]$print.options$treat.names == b$print.options$treat.names[treat.names.good]))]
                }
                else {
                    disp.treat.pairs <- names(b[["Pair.Balance"]])[sapply(names(b[["Pair.Balance"]]), function(x) all(b[["Pair.Balance"]][[x]]$print.options$treat.names %in% b$print.options$treat.names[treat.names.good]))]
                }
            }
            else {
                if (sum(treat.names.good) == 1) {
                    disp.treat.pairs <- names(b[["Pair.Balance"]])[sapply(names(b[["Pair.Balance"]]), function(x) {
                        treat.names <- b[["Pair.Balance"]][[x]]$print.options$treat.namestreat.names
                        any(treat.names[treat.names != "Others"] == b$print.options$treat.names[treat.names.good])})]
                }
                else {
                    disp.treat.pairs <- names(b[["Pair.Balance"]])[sapply(names(b[["Pair.Balance"]]), function(x) {
                        treat.names <- b[["Pair.Balance"]][[x]]$print.options$treat.namestreat.names
                        all(treat.names[treat.names != "Others"] %in% b$print.options$treat.names[treat.names.good])})]
                }
            }

        }
        
        #Get B from b
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(disp.treat.pairs,
                                         function(x) cbind(b[["Pair.Balance"]][[x]][["Balance"]],
                                                           treat.pair = x,
                                                           variable.names = rownames(b[["Pair.Balance"]][[x]][["Balance"]]))))
            facet <- "treat.pair"
        }
        else if (config == "agg.pair") {
            #Agg.Fun <- switch(match.arg(agg.fun), mean = "Mean", median = "Median", max = "Max", range = "Range")
            Agg.Fun <- "Max"
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Treatment", ifelse(b$print.options$pairwise, " Pairs", "s"))
            }
            else {
                subtitle <- paste0(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), " ", which.stat2, " Across Treatment", ifelse(b$print.options$pairwise, " Pairs", "s"))
            }
            B <- cbind(b[["Balance.Across.Pairs"]],
                       variable.names = rownames(b[["Balance.Across.Pairs"]]))
        }
    }
    else if (any(class(b) == "bal.tab.subclass")) {
        if (which.stat=="V.Ratio") stop("Variance ratios not currently supported for subclassification.", call. = FALSE)
        if (which.stat=="KS") stop("KS statistics not currently supported for subclassification.", call. = FALSE)
        if (any(class(b) == "bal.tab.cont")) stop("Continuous treatments not currently supported for subclassification.", call. = FALSE)
        subclass.names <- names(b[["Subclass.Balance"]])
        sub.B <- do.call("cbind", lapply(subclass.names, function(x) {
            sub <- b[["Subclass.Balance"]][[x]]
            sub.B0 <- setNames(sub[endsWith(names(sub), ".Adj")],
                               gsub(".Adj", paste0(".Subclass ", x), names(sub)[endsWith(names(sub), ".Adj")]))
            return(sub.B0) }))
        B <- cbind(b[["Balance.Across.Subclass"]], sub.B, variable.names = row.names(b[["Balance.Across.Subclass"]]))
        if (b$print.options$disp.subclass) b$print.options$weight.names <- c("Adj", paste("Subclass", subclass.names))
        else b$print.options$weight.names <- "Adj"
        subtitle <- "Across Subclasses"
    }
    else {
        B <- cbind(b[["Balance"]], variable.names = row.names(b[["Balance"]]))
    }
    
    if (config == "agg.none") {
        if (!abs && !b$print.options[["abs"]]) {
            abs <- FALSE
        }
        else {
            if (!abs && b$print.options[["abs"]])
                warning("abs is TRUE in the bal.tab object; ignoring abs in the call to love.plot().", call. = FALSE)
            abs <- TRUE
        }
    }
    else if (config %in% c("agg.time", "agg.pair")) {
        abs <- TRUE
    }
    else if (config %in% c("agg.imp", "agg.cluster", "agg.all")) {
        abs <- b$print.options[["abs"]]
    }
    
    if (null.threshold) {
        if (which.stat=="Diff") {
            if (is_not_null(b$print.options$m.threshold)) {
                threshold <- b$print.options$m.threshold
                null.threshold <- FALSE
            }
        }
        else if (which.stat=="V.Ratio") {
            if (is_not_null(b$print.options$v.threshold)) {
                threshold <- b$print.options$v.threshold
                null.threshold <- FALSE
            }
        }
        else if (which.stat=="KS") {
            if (is_not_null(b$print.options$ks.threshold)) {
                threshold <- b$print.options$ks.threshold
                null.threshold <- FALSE
            }
        }
        else if (which.stat=="Corr") {
            if (is_not_null(b$print.options$r.threshold)) {
                threshold <- b$print.options$r.threshold
                null.threshold <- FALSE
            }
        }
    }
    
    if (is_not_null(var.names)) {
        if (is.data.frame(var.names)) {
            if (ncol(var.names)==1) {
                if (is_not_null(row.names(var.names))) {
                    new.labels <- setNames(unlist(as.character(var.names[,1])), rownames(var.names))
                }
                else warning("var.names is a data.frame, but its rows are unnamed.", call. = FALSE)
            }
            else {
                if (all(c("old", "new") %in% names(var.names))) {
                    new.labels <- setNames(unlist(as.character(var.names[,"new"])), var.names[,"old"])
                }
                else {
                    if (ncol(var.names)>2) warning("Only using first 2 columns of var.names", call. = FALSE)
                    new.labels <- setNames(unlist(as.character(var.names[,2])), var.names[,1])
                }
            } 
        }
        else if (is.atomic(var.names) || is.factor(var.names)) {
            if (is_not_null(names(var.names))) {
                new.labels <- setNames(as.character(var.names), names(var.names))
            }
            else warning("var.names is a vector, but its values are unnamed.", call. = FALSE)
        }
        else if (is.list(var.names)) {
            if (all(sapply(var.names, function(x) is.character(x) || is.factor(x)))) {
                if (is_not_null(names(var.names))) {
                    new.labels <- unlist(var.names) #already a list
                }
                else warning("var.names is a list, but its values are unnamed.", call. = FALSE)
            }
            else warning("var.names is a list, but its values are not the new names of the variables.", call. = FALSE)
        }
        else warning("Argument to var.names is not one of the accepted structures and will be ignored.\n  See help(love.plot) for details.", immediate.=TRUE, call. = FALSE)

        co.names <- b[["print.options"]][["co.names"]]
        seps <- attr(co.names, "seps")
        for (i in names(co.names)) {
            comp <- co.names[[i]][["component"]]
            is.name <- co.names[[i]][["is.name"]]
            
            if (i %in% names(new.labels) && !is.na(new.labels[i])) co.names[[i]][["component"]] <- new.labels[i]
            else {
                if (seps["int"] %in% comp[!is.name]) {
                    named.vars <- character(sum(comp[!is.name] == seps["int"]) + 1)
                    sep.inds <- c(which(comp == seps["int"] & !is.name), length(comp) + 1)
                    named.vars <- lapply(seq_along(sep.inds), function(k) {
                        inds <- (if (k == 1) seq(1, sep.inds[k] - 1) 
                                 else seq(sep.inds[k-1] + 1, sep.inds[k] - 1))
                        var <- comp[inds]
                        var.is.name <- is.name[inds]
                        pasted.var <- paste(var, collapse = "")
                        if (pasted.var %in% names(new.labels)) return(new.labels[pasted.var])
                        else return(paste(ifelse(var.is.name & var %in% names(new.labels) & !is.na(new.labels[var]), new.labels[var], var), collapse = ""))
                    })
                    co.names[[i]][["component"]] <- do.call("paste", c(unname(named.vars), list(sep = seps["int"])))
                }
                else co.names[[i]][["component"]] <- ifelse(is.name & comp %in% names(new.labels) & !is.na(new.labels[comp]), new.labels[comp], comp)
            }
        }
        
        recode.labels <- setNames(names(co.names), 
                                  vapply(co.names, function(x) paste0(x[["component"]], collapse = ""), character(1L)))
        
        B[["variable.names"]] <- do.call(f.recode, c(list(B[["variable.names"]]), as.list(recode.labels)))
    }
    
    distance.names <- as.character(unique(B[["variable.names"]][B[["Type"]] == "Distance"], nmax = sum(B[["Type"]] == "Distance")))
    if (drop.distance) {
        B <- B[B[["variable.names"]] %nin% distance.names, , drop = FALSE]
    }
    
    if (is_not_null(var.order)) {
        if (b$print.options$nweights == 0) {
            ua <- c("Unadjusted", "Alphabetical")
            names(ua) <- c("unadjusted", "alphabetical")
        }
        else if (b$print.options$nweights == 1) {
            ua <- c("Adjusted", "Unadjusted", "Alphabetical")
            names(ua) <- c("adjusted", "unadjusted", "alphabetical")
        }
        else {
            ua <- c("Unadjusted", b$print.options$weight.names, "Alphabetical")
            names(ua) <- c("unadjusted", b$print.options$weight.names, "alphabetical")
        }
        var.order <- ua[match.arg(var.order, tolower(ua))]
    }

    if (is_not_null(facet)) {
        if (is_not_null(var.order) && tolower(var.order) != "alphabetical" && (sum(cluster.names.good) > 1 || sum(imp.numbers.good) > 1 || length(disp.treat.pairs) > 1 || sum(time.names.good) > 1)) {
            warning("var.order cannot be set with multiple plots (unless \"alphabetical\"). Ignoring var.order.", call. = FALSE)
            var.order <- NULL
        }
    }
    
    if (is_not_null(Agg.Fun) && Agg.Fun == "Range") {
        agg.range <- TRUE
    }
    else agg.range <- FALSE
    
    if (agg.range) {
        SS <- do.call("rbind", lapply(c("Un", b$print.options$weight.names),
                                      function(x) data.frame(var = B[["variable.names"]],
                                                             min.stat = B[[paste("Min", which.stat, x, sep = ".")]],
                                                             max.stat = B[[paste("Max", which.stat, x, sep = ".")]],
                                                             mean.stat = B[[paste("Mean", which.stat, x, sep = ".")]],
                                                             Sample = ifelse(x == "Un", "Unadjusted", 
                                                                             ifelse(x == "Adj", "Adjusted", x)))))

        if (is_not_null(facet)) {
            if ("cluster" %in% facet) {
                SS$cluster <- rep(B[["cluster"]], 1 + b$print.options$nweights)
            }
            if ("imp" %in% facet) {
                SS$imp <- rep(B[["imp"]], 1 + b$print.options$nweights)
            }
            if ("treat.pair" %in% facet) {
                SS$treat.pair <- rep(B[["treat.pair"]], 1 + b$print.options$nweights)
            }
        }
        
        if (all(sapply(SS[c("min.stat", "max.stat", "mean.stat")], is.na))) stop("No balance statistics to display.", call. = FALSE)
        gone <- character(0)
        for (i in levels(SS$Sample)) {
            if (all(sapply(SS[SS$Sample==i, c("min.stat", "max.stat", "mean.stat")], is.na))) {
                gone <- c(gone, i)
                if (i == "Unadjusted") warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
                SS <- SS[SS[["Sample"]]!=i,]
            }
        }
       
        if (abs) dec <- FALSE
        else dec <- FALSE
        
        if (is_not_null(var.order)) {
            if (var.order %in% ua) {
                if (var.order %in% gone) {
                    warning(paste0("var.order was set to \"", tolower(var.order), "\", but no ", tolower(var.order), " ", tolower(which.stat2), "s were calculated. Ignoring var.order."), call. = FALSE, immediate. = TRUE)
                    var.order <- character(0)
                }
                else {
                v <- as.character(SS[["var"]][order(SS[["mean.stat"]][SS[["Sample"]]==var.order], decreasing = dec)])
                SS[["var"]] <- factor(SS[["var"]], 
                                      levels=c(v[is.na(match(v, distance.names))], 
                                               sort(distance.names, decreasing = TRUE)))
                }
            }
            else if (var.order == "alphabetical") {
                if ("time" %in% facet) {
                    covnames0 <- vector("list", length(unique(SS[["time"]])))
                    for (i in seq_along(covnames0)) {
                        if (i == 1) {
                            covnames0[[i]] <- sort(as.character(unique(SS[["var"]][SS[["time"]] == i])))
                        }
                        else {
                            covnames0[[i]] <- sort(setdiff(as.character(unique(SS[["var"]][SS[["time"]] == i])), unlist(covnames0[seq_along(covnames0) < i])))
                        }
                    }
                    covnames <- unlist(covnames0)
                }
                else covnames <- sort(SS[["var"]])
                SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
            }
        }
        if (is_null(var.order)) {
            covnames <- as.character(unique(SS[["var"]]))
            SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
        }
        # SS[, "Sample"] <- factor(SS[, "Sample"], levels = c("Adjusted", "Unadjusted"))
        SS[["Sample"]] <- factor(SS[["Sample"]])
        if (which.stat == "Diff" && any(abs(SS[["max.stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardizied mean differences for continuous variables. To do so, specify continuous=\"std\" in bal.tab().", call.=FALSE, noBreaks.=TRUE)
        if (no.missing) SS <- SS[!is.na(SS[["min.stat"]]),]
        SS[["stat"]] <- SS[["mean.stat"]]
    }
    else {
        SS <- do.call("rbind", lapply(c("Un", b$print.options$weight.names),
                                      function(x) data.frame(var = B[["variable.names"]],
                                                             stat = B[[ifelse(is_null(Agg.Fun), paste(which.stat, x, sep = "."),
                                                                               paste(Agg.Fun, which.stat, x, sep = "."))]],
                                                             Sample = ifelse(x == "Un", "Unadjusted", 
                                                                             ifelse(x == "Adj", "Adjusted", x)))))
        
        if (is_not_null(facet)) {
            for (i in c("cluster", "imp", "treat.pair", "time")) {
                if (i %in% facet) {
                    SS[[i]] <- rep(B[[i]], 1 + b$print.options$nweights)
                }
            }
        }
        
        if (all(is.na(SS[["stat"]]))) stop("No balance statistics to display.", call. = FALSE)
        gone <- character(0)
        for (i in levels(SS$Sample)) {
            if (all(sapply(SS[["stat"]][SS$Sample==i], is.na))) {
                gone <- c(gone, i)
                if (i == "Unadjusted") warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
                SS <- SS[SS[["Sample"]]!=i,]
            }
        }

        if (abs) {
            SS[["stat"]] <- abs(SS[["stat"]])
            dec <- FALSE
        }
        else dec <- FALSE
        
        if (is_not_null(var.order)) {
            if (tolower(var.order) == "alphabetical") {
                if ("time" %in% facet) {
                    covnames0 <- vector("list", length(unique(SS[["time"]])))
                    for (i in seq_along(covnames0)) {
                        if (i == 1) {
                            covnames0[[i]] <- sort(as.character(unique(SS[["var"]][SS[["time"]] == i])))
                        }
                        else {
                            covnames0[[i]] <- sort(setdiff(as.character(unique(SS[["var"]][SS[["time"]] == i])), unlist(covnames0[seq_along(covnames0) < i])))
                        }
                    }
                    covnames <- unlist(covnames0)
                }
                else covnames <- sort(SS[["var"]])
                SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
            }
            else if (var.order %in% ua) {
                if (var.order %in% gone) {
                    warning(paste0("var.order was set to \"", tolower(var.order), "\", but no ", tolower(var.order), " ", tolower(which.stat2), "s were calculated. Ignoring var.order."), call. = FALSE, immediate. = TRUE)
                    var.order <- character(0)
                }
                else {
                    v <- as.character(SS[["var"]][order(SS[["stat"]][SS[["Sample"]]==var.order], decreasing = dec)])
                    SS[["var"]] <- factor(SS[["var"]], 
                                          levels=c(v[is.na(match(v, distance.names))], 
                                                   sort(distance.names, decreasing = TRUE)))
                }
            }
            
        }
        if (is_null(var.order)) {
            covnames <- as.character(unique(SS[["var"]]))
            SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
        }
        SS[["Sample"]] <- factor(SS[["Sample"]])
        if (which.stat == "Diff" && any(abs(SS[["stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardizied mean differences for continuous variables. To do so, specify continuous=\"std\" in bal.tab().", call.=FALSE, noBreaks.=TRUE)
        if (no.missing) SS <- SS[!is.na(SS[["stat"]]),]
    }
    SS <- SS[order(SS[["var"]]),]
    SS[["var"]] <- factor(SS[["var"]])
    
    #Make the plot
    #library(ggplot2)
    
    #Setting up appearance
    #Size
    if (is_null(args$size)) size <- 1
    else if (is.numeric(args$size[1])) size <- args$size[1]
    else {
        warning("The argument to size must be a number. Using 1 instead.", call. = FALSE)
        size <- 1
    }
    stroke <- .8*size
    
    #Color
    ntypes <- if (is_null(subclass.names)) nlevels(SS$Sample) else 2
    if (is_not_null(args$colours)) colors <- args$colours
    
    if (is_null(colors)) {
        if (shapes.ok(shapes, ntypes) && length(shapes) == ntypes) {
            colors <- rep("black", ntypes)
        }
        else colors <- gg_color_hue(ntypes)
    }
    else {
        if (length(colors) == 1) colors <- rep(colors, ntypes)
        else if (length(colors) > ntypes) {
            colors <- colors[seq_len(ntypes)]
            warning(paste("Only using first", ntypes, "value", if (ntypes>1) "s " else " ", "in colors."), call. = FALSE)
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
    fill <- colors
    
    #Shapes
    if (is_null(shapes)) {
        shapes <- assign.shapes(colors)
    }
    else {
        #check shapes
        if (!shapes.ok(shapes, ntypes)) {
            warning(paste("The argument to shape must be", ntypes, "valid shape", if (ntypes>1) "s." else ".", " See ?love.plot for more information.\nUsing default shapes instead."), call. = FALSE)
            shapes <- assign.shapes(colors)
        }
        else if (length(shapes) == 1) shapes <- rep(shapes, ntypes)
    }
    
    #Title
    if (is_not_null(args$title)) title <- as.character(args$title)
    if (is_not_null(args$subtitle)) subtitle <- as.character(args$subtitle)
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
    else if (which.stat == "KS") {
        baseline.xintercept <- 0
        xlab <- "Kolmogorov-Smirnov Statistics"
        if (!null.threshold) threshold.xintercept <- abs(threshold)
    }
    
    apply.limits <- FALSE
    if (is_not_null(args$limits)) {
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
                if (any(SS[["min.stat"]] < limits[1], na.rm = TRUE) || any(SS[["max.stat"]] > limits[2], na.rm = TRUE)) {
                    for (i in c("min.stat", "stat", "max.stat")) {
                        SS[[i]][SS[[i]] < limits[1]] <- limits[1]
                        SS[[i]][SS[[i]] > limits[2]] <- limits[2]
                    }
                    warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
            }
            else {
                if (any(SS[["stat"]] < limits[1], na.rm = TRUE) || any(SS[["stat"]] > limits[2], na.rm = TRUE)) {
                    warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
            }
            apply.limits <- TRUE
        }
    }
    if (is_not_null(args$sample.names)) {
        if (!all(is.character(args$sample.names))) {
            warning("The argument to sample.names must be a character vector. Ignoring sample.names.", call. = FALSE)
        }
        else {
            if (length(args$sample.names) == ntypes - 1) {
                levels(SS$Sample)[-1] <- args$sample.names
            }
            else if (length(args$sample.names) == ntypes) {
                levels(SS$Sample) <- args$sample.names
            }
            else {
                warning("The argument to sample.names must contain as many names as there are sample types, or one fewer. Ignoring sample.names.", call. = FALSE)
            }
        }

    }

    lp <- ggplot(data = SS, aes(y = var, x = stat, group = Sample)) + 
        theme(panel.grid.major = element_line(color = "gray87"),
              panel.grid.minor = element_line(color = "gray90"),
              panel.background = element_rect(fill = "white", color = "black"),
              axis.text.x = element_text(color = "black"),
              axis.text.y = element_text(color = "black")
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
            f <- function(x) {x[["stat"]][x$var %in% distance.names] <- NA; x}
            lp <- lp + ggplot2::layer(geom = "path", data = f, position = position.dodge, stat = "identity", 
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
        
        if (is_null(subclass.names) || !b$print.options$disp.subclass) {
            if (line == TRUE) { #Add line except to distance
                f <- function(x) {x[["stat"]][x$var %in% distance.names] <- NA; x}
                lp <- lp + ggplot2::layer(geom = "path", data = f(SS), position = "identity", stat = "identity", aes(color = Sample), params = list(size = size*.8, na.rm = TRUE))
            }
            lp <- lp + geom_point(data = SS, aes(shape = Sample,
                                                 color = Sample), 
                                  size = 2*size, stroke = stroke, fill = "white", na.rm = TRUE) 
        }
        else {
            SS.u.a <- SS[SS$Sample %in% c("Unadjusted", "Adjusted"),]
            SS.u.a$Sample <- factor(SS.u.a$Sample)
            if (line == TRUE) { #Add line except to distance
                f <- function(x) {x[["stat"]][x$var %in% distance.names] <- NA; x}
                lp <- lp + ggplot2::layer(geom = "path", data = f(SS.u.a), position = "identity", stat = "identity", aes(color = Sample), params = list(size = size*.8, na.rm = TRUE))
            }
            lp <- lp + geom_point(data = SS.u.a, 
                                  aes(shape = Sample, color = Sample), 
                                  size = 2*size, stroke = stroke, fill = "white", na.rm = TRUE) 
            lp <- lp + geom_text(data = SS[!SS$Sample %in% c("Unadjusted", "Adjusted"),], 
                                 aes(label = gsub("Subclass ", "", Sample)), 
                                 size = 2*size, na.rm = TRUE) 
        }
        
        
    }

    if (!drop.distance && is_not_null(distance.names)) {
        lp <- lp + geom_hline(linetype = 1, color = "black", 
                              yintercept = nunique(SS[["var"]]) - length(distance.names) + .5)
    }
    if (apply.limits) {
        lp <- lp + scale_x_continuous(limits = limits, expand = c(0, 0))
    }
    if (is_not_null(facet)) {
        lp <- lp + facet_grid(f.build(".", facet), drop = FALSE)
    }

    return(lp)
}
plot.bal.tab <- love.plot
