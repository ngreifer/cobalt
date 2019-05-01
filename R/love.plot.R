love.plot <- function(x, stat = c("mean.diffs", "variance.ratios", "ks.statistics"), threshold = NULL, 
                      abs = TRUE, var.order = NULL, no.missing = TRUE, var.names = NULL, 
                      drop.distance = FALSE, agg.fun = c("range", "max", "mean"), 
                      colors = NULL, shapes = NULL, line = FALSE, ...) {
    
    #Re-call bal.tab with disp.v.ratio or disp.ks if stat = "v" or "k".
    if (!exists(deparse(substitute(x)))) {
        m <- m0 <- match.call()
        if (deparse(m[["x"]][[1]]) %in% c("bal.tab", methods("bal.tab"))) {
            if (pmatch(stat[1], "variance.ratios", 0L) != 0L) {
                if (!isTRUE(eval(m[["x"]][["disp.v.ratio"]]))) {
                    m[["x"]][["un"]] <- TRUE
                    m[["x"]][["disp.v.ratio"]] <- TRUE
                }
            }
            else if (pmatch(stat[1], "ks.statistics", 0L) != 0L) {
                if (!isTRUE(eval(m[["x"]][["disp.ks"]]))) {
                    m[["x"]][["un"]] <- TRUE
                    m[["x"]][["disp.ks"]] <- TRUE
                    
                }
            }
            
            if (any(names(m[["x"]]) == "cluster")) {
                m[["x"]][["cluster.summary"]] <- TRUE
                if (any(names(m[["x"]]) == "cluster.fun")) m[["x"]][["cluster.fun"]] <- NULL
            }
            if (any(names(m[["x"]]) == "imp")) {
                m[["x"]][["imp.summary"]] <- TRUE
                if (any(names(m[["x"]]) == "imp.fun")) m[["x"]][["imp.fun"]] <- NULL
            }
            
            if (!identical(m, m0)) return(eval(m))
        }
    }
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
        .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
        return(eval(.call))
    }
    
    if (!any(class(x) == "bal.tab")) {
        #Use bal.tab on inputs first, then love.plot on that
        m <- match.call()
        m.b <- m; m.b[[1]] <- quote(bal.tab); names(m.b)[2] <- ''
        
        m.l <- m; 
        m.l[["x"]] <- m.b

        return(eval(m.l))
    }
    
    if (any(class(x) == "bal.tab.cont")) {
        if (pmatch(stat[1], "correlations", 0L) == 0L && !identical(stat, eval(formals()[["stat"]]))) {
            message("Displaying treatment-covariate correlations.")
        }
        stat <- "correlations"
    }
    else stat <- match_arg(stat)
    
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
        if (i %in% names(args)) x$print.options[[i]] <- args[[i]]
    }
    
    null.threshold <- is_null(threshold)
    if (!null.threshold) {
        if (!is.numeric(threshold) || length(threshold) > 1) stop("threshold must be a single number.", call. = FALSE)
    }
    
    if (is_not_null(args$cluster.fun) && is_null(agg.fun)) agg.fun <- args$cluster.fun
    which.stat <- switch(stat, mean.diffs = "Diff", variance.ratios = "V.Ratio", ks.statistics = "KS", correlations = "Corr")
    which.stat2 <- switch(stat, mean.diffs = "Mean Difference", variance.ratios = "Variance Ratio", ks.statistics = "Kolmogorov-Smirnov Statistic", correlations = "Correlation")
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
    null.cluster <- is_null(x$print.options$which.cluster)
    na.cluster <- !null.cluster && is.na(x$print.options$which.cluster)
    null.imp <- is_null(x$print.options$which.imp)
    na.imp <- !null.imp && is.na(x$print.options$which.imp)
    null.treat <- is_null(x$print.options$which.treat)
    na.treat <- !null.treat && is.na(x$print.options$which.treat)
    null.time <- is_null(x$print.options$which.time)
    na.time <- !null.time && is.na(x$print.options$which.time)
    
    config <- "agg.none"
    
    if (any(class(x) == "bal.tab.msm")) {
        #Get time.names.good
        time.names <- names(x[["Time.Balance"]])
        if (na.time) {
            config <- "agg.time"
            time.names.good <- NULL
        }
        else if (null.time) {
            config <- "agg.none"
            time.names.good <- setNames(rep(TRUE, length(time.names)), time.names)
        }
        else if (is.character(x$print.options$which.time)) {
            time.names.good <- sapply(time.names, function(x) any(sapply(x$print.options$which.time, function(y) identical(x, y))))
            if (any(time.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.time are valid treatment names or time point indices.", call. = FALSE)
            }
        }
        else if (is.numeric(x$print.options$which.time)) {
            time.names.good <- setNames(seq_along(x$Time.Balance) %in% x$print.options$which.time, time.names)
            if (any(time.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.time are valid treatment names or time point indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.time must be NULL, NA, or the names of treatments or indices of time points.", call. = FALSE)
        
        
        #Get B from x
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(x[["Time.Balance"]])[time.names.good], 
                                         function(t) cbind(x[["Time.Balance"]][[t]][["Balance"]],
                                                           time = t,
                                                           variable.names = rownames(x[["Time.Balance"]][[t]][["Balance"]]))))
            facet <- "time"
        }
        else if (config == "agg.time") {
            if (is_null(x[["Balance.Across.Times"]])) {
                stop("Cannot aggregate across time periods without a balance summary across time periods.\nThis may be because multinomial treatments were used, multiple treatment types were used,\n or quick was set to TRUE and msm.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            #Agg.Fun <- switch(match_arg(agg.fun), mean = "Mean", max = "Max", range = "Range")
            Agg.Fun <- "Max"
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Time Points")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Time Points")
            }
            B <- cbind(x[["Balance.Across.Times"]],
                       variable.names = rownames(x[["Balance.Across.Times"]]))
        }
    }
    else if (any(class(x) == "bal.tab.imp.cluster")) {
        #Get imp.numbers.good
        imp.numbers <- seq_along(x[["Imputation.Balance"]])
        if (na.imp) {
            imp.numbers.good <- NULL
        }
        else if (null.imp) {
            imp.numbers.good <- setNames(rep(TRUE, length(imp.numbers)), imp.numbers)
        }
        else if (is.numeric(x$print.options$which.imp)) {
            imp.numbers.good <- setNames(imp.numbers %in% x$print.options$which.imp, imp.numbers)
        }
        else stop("The argument to which.imp must be NULL, NA, or the indices of imputations.", call. = FALSE)
        
        #Get cluster.names.good
        cluster.names <- names(x[["Cluster.Balance.Across.Imputations"]])
        if (na.cluster) {
            cluster.names.good <- NULL
        }
        else if (null.cluster) {
            cluster.names.good <- setNames(rep(TRUE, length(cluster.names)), cluster.names)
        }
        else if (is.character(x$print.options$which.cluster)) {
            cluster.names.good <- sapply(cluster.names, function(n) any(sapply(x$print.options$which.cluster, function(y) identical(n, y))))
        }
        else if (is.numeric(x$print.options$which.cluster)) {
            cluster.names.good <- setNames(seq_along(x$Cluster.Balance) %in% x$print.options$which.cluster, cluster.names)
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
        
        #Get B from x based on configuration
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(x[["Imputation.Balance"]])[imp.numbers.good],
                                         function(i) do.call("rbind", lapply(names(x[["Imputation.Balance"]][[i]][["Cluster.Balance"]])[cluster.names.good],
                                                                             function(y) cbind(x[["Imputation.Balance"]][[i]][["Cluster.Balance"]][[y]][["Balance"]],
                                                                                               cluster = y,
                                                                                               imp = paste("Imputation:", i),
                                                                                               variable.names = rownames(x[["Imputation.Balance"]][[i]][["Cluster.Balance"]][[y]][["Balance"]]))))))
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
            if (is_null(x[["Cluster.Balance.Across.Imputations"]])) {
                stop("Cannot aggregate across imputations without a balance summary across imputations.\nThis may be because quick was set to TRUE and imp.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            Agg.Fun <- switch(tolower(match_arg(agg.fun)), mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Imputations")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Imputations", sep = " ")
            }
            B <- do.call("rbind", lapply(names(x[["Cluster.Balance.Across.Imputations"]])[cluster.names.good], 
                                         function(c) cbind(x[["Cluster.Balance.Across.Imputations"]][[c]][["Cluster.Balance"]], 
                                                           cluster = c, 
                                                           variable.names = rownames(x[["Cluster.Balance.Across.Imputations"]][[c]][["Cluster.Balance"]]))))
            facet <- "cluster"
        }
        else if (config == "agg.cluster") {
            if (is_null(x[["Imputation.Balance"]][[1]][["Cluster.Summary"]])) {
                stop("Cannot aggregate across clusters without a balance summary across clusters.\nThis may be because quick was set to TRUE and cluster.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            Agg.Fun <- switch(tolower(match_arg(agg.fun)), mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Clusters")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Clusters", sep = " ")
            }
            B <- do.call("rbind", lapply(names(x[["Imputation.Balance"]])[imp.numbers.good], 
                                         function(i) cbind(x[["Imputation.Balance"]][[i]][["Cluster.Summary"]], 
                                                           imp = paste("Imputation:", i), 
                                                           variable.names = rownames(x[["Imputation.Balance"]][[i]][["Cluster.Summary"]]))))
            facet <- "imp"
        }
        else if (config == "agg.all") {
            if (is_null(x[["Balance.Across.Imputations"]])) {
                stop("Cannot aggregate across imputations without a balance summary across imputations.\nThis may be because quick was set to TRUE and cluster.summary or imp.summary were set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            #Cluster.Fun <- switch(match_arg(cluster.fun), mean = "Mean", max = "Max", range = "Range")
            Agg.Fun <- switch(tolower(match_arg(agg.fun)), mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Clusters and Imputations")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Clusters and Imputations", sep = " ")
            }
            B <- cbind(x[["Balance.Across.Imputations"]],
                       variable.names = row.names(x[["Balance.Across.Imputations"]]))
            facet <- NULL
        }
    }
    else if (any(class(x) == "bal.tab.imp")) {
        #Get imp.numbers.good
        imp.numbers <- seq_along(x[["Imputation.Balance"]])
        if (na.imp) {
            config <- "agg.imp"
            imp.numbers.good <- NULL
        }
        else if (null.imp) {
            config <- "agg.none"
            imp.numbers.good <- setNames(rep(TRUE, length(imp.numbers)), imp.numbers)
        }
        else if (is.numeric(x$print.options$which.imp)) {
            imp.numbers.good <- setNames(imp.numbers %in% x$print.options$which.imp, imp.numbers)
            if (any(imp.numbers.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.imp are valid imputation indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.imp must be NULL, NA, or the indices of imputations.", call. = FALSE)
        
        #Get B from x
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(x[["Imputation.Balance"]])[imp.numbers.good], 
                                         function(i) cbind(x[["Imputation.Balance"]][[i]][["Balance"]],
                                                           imp = paste("Imputation:", i),
                                                           variable.names = rownames(x[["Imputation.Balance"]][[i]][["Balance"]]))))
            facet <- "imp"
        }
        else if (config == "agg.imp") {
            if (is_null(x[["Balance.Across.Imputations"]])) {
                stop("Cannot aggregate across imputations without a balance summary across imputations.\nThis may be because quick was set to TRUE and imp.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            Agg.Fun <- switch(tolower(match_arg(agg.fun)), mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Imputations")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Imputations")
            }
            B <- cbind(x[["Balance.Across.Imputations"]],
                       variable.names = rownames(x[["Balance.Across.Imputations"]]))
        }
    }
    else if (any(class(x) == "bal.tab.cluster")) {
        #Get cluster.names.good
        cluster.names <- names(x[["Cluster.Balance"]])
        if (na.cluster) {
            config <- "agg.cluster"
            cluster.names.good <- NULL
        }
        else if (null.cluster) {
            config <- "agg.none"
            cluster.names.good <- setNames(rep(TRUE, length(cluster.names)), cluster.names)
        }
        else if (is.character(x$print.options$which.cluster)) {
            cluster.names.good <- sapply(cluster.names, function(c) any(sapply(x$print.options$which.cluster, function(y) identical(c, y))))
            if (any(cluster.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
            }
        }
        else if (is.numeric(x$print.options$which.cluster)) {
            cluster.names.good <- setNames(seq_along(x$Cluster.Balance) %in% x$print.options$which.cluster, cluster.names)
            if (any(cluster.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.cluster must be NULL, NA, or the names or indices of clusters.", call. = FALSE)
        
        
        #Get B from x
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(names(x[["Cluster.Balance"]])[cluster.names.good], 
                                         function(c) cbind(x[["Cluster.Balance"]][[c]][["Balance"]],
                                                           cluster = c,
                                                           variable.names = rownames(x[["Cluster.Balance"]][[c]][["Balance"]]))))
            facet <- "cluster"
        }
        else if (config == "agg.cluster") {
            if (is_null(x[["Cluster.Summary"]])) {
                stop("Cannot aggregate across clusters without a balance summary across clusters.\nThis may be because quick was set to TRUE and cluster.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            
            tryCatch({agg.fun <- tolower(match_arg(agg.fun))}, 
                     error = function(e) stop("agg.fun should be one of \"mean\", \"max\", or \"range\".", call. = FALSE))
            Agg.Fun <- switch(agg.fun, mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Clusters")
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), which.stat2, "Across Clusters")
            }
            B <- cbind(x[["Cluster.Summary"]],
                       variable.names = rownames(x[["Cluster.Summary"]]))
        }
    }
    else if (any(class(x) == "bal.tab.multi")) {
        #Get cluster.names.good
        treat.names <- x$print.options$treat.names
        if (na.treat) {
            config <- "agg.pair"
            treat.names.good <- NULL
        }
        else if (null.treat) {
            config <- "agg.none"
            treat.names.good <- rep(TRUE, length(treat.names))
        }
        else if (is.character(x$print.options$which.treat)) {
            treat.names.good <- treat.names %in% x$print.options$which.treat
            if (any(treat.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.treat are valid treatment names or indices.", call. = FALSE)
            }
        }
        else if (is.numeric(x$print.options$which.treat)) {
            treat.names.good <- seq_along(treat.names) %in% x$print.options$which.treat
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
            if (x$print.options$pairwise) {
                if (sum(treat.names.good) == 1) {
                    disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) any(x[["Pair.Balance"]][[p]]$print.options$treat.names == x$print.options$treat.names[treat.names.good]))]
                }
                else {
                    disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) all(x[["Pair.Balance"]][[p]]$print.options$treat.names %in% x$print.options$treat.names[treat.names.good]))]
                }
            }
            else {
                if (sum(treat.names.good) == 1) {
                    disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) {
                        treat.names <- x[["Pair.Balance"]][[p]]$print.options$treat.namestreat.names
                        any(treat.names[treat.names != "Others"] == x$print.options$treat.names[treat.names.good])})]
                }
                else {
                    disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) {
                        treat.names <- x[["Pair.Balance"]][[p]]$print.options$treat.namestreat.names
                        all(treat.names[treat.names != "Others"] %in% x$print.options$treat.names[treat.names.good])})]
                }
            }
            
        }
        
        #Get B from x
        if (config == "agg.none") {
            B <- do.call("rbind", lapply(disp.treat.pairs,
                                         function(p) cbind(x[["Pair.Balance"]][[p]][["Balance"]],
                                                           treat.pair = p,
                                                           variable.names = rownames(x[["Pair.Balance"]][[p]][["Balance"]]))))
            facet <- "treat.pair"
        }
        else if (config == "agg.pair") {
            #Agg.Fun <- switch(match_arg(agg.fun), mean = "Mean", max = "Max", range = "Range")
            Agg.Fun <- "Max"
            if (Agg.Fun == "Range") {
                subtitle <- paste0(which.stat2, " Range Across Treatment", ifelse(x$print.options$pairwise, " Pairs", "s"))
            }
            else {
                subtitle <- paste0(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), " ", which.stat2, " Across Treatment", ifelse(x$print.options$pairwise, " Pairs", "s"))
            }
            B <- cbind(x[["Balance.Across.Pairs"]],
                       variable.names = rownames(x[["Balance.Across.Pairs"]]))
        }
    }
    else if (any(class(x) == "bal.tab.subclass")) {
        if (which.stat=="V.Ratio") stop("Variance ratios not currently supported for subclassification.", call. = FALSE)
        if (which.stat=="KS") stop("KS statistics not currently supported for subclassification.", call. = FALSE)
        if (any(class(x) == "bal.tab.cont")) stop("Continuous treatments not currently supported for subclassification.", call. = FALSE)
        subclass.names <- names(x[["Subclass.Balance"]])
        sub.B <- do.call("cbind", lapply(subclass.names, function(x) {
            sub <- x[["Subclass.Balance"]][[x]]
            sub.B0 <- setNames(sub[endsWith(names(sub), ".Adj")],
                               gsub(".Adj", paste0(".Subclass ", x), names(sub)[endsWith(names(sub), ".Adj")]))
            return(sub.B0) }))
        B <- cbind(x[["Balance.Across.Subclass"]], sub.B, variable.names = row.names(x[["Balance.Across.Subclass"]]))
        if (x$print.options$disp.subclass) x$print.options$weight.names <- c("Adj", paste("Subclass", subclass.names))
        else x$print.options$weight.names <- "Adj"
        subtitle <- "Across Subclasses"
    }
    else {
        B <- cbind(x[["Balance"]], variable.names = row.names(x[["Balance"]]))
    }
    
    if (config == "agg.none") {
        if (!abs && !x$print.options[["abs"]]) {
            abs <- FALSE
        }
        else {
            if (!abs && x$print.options[["abs"]])
                warning("abs is TRUE in the bal.tab object; ignoring abs in the call to love.plot().", call. = FALSE)
            abs <- TRUE
        }
    }
    else if (config %in% c("agg.time", "agg.pair")) {
        abs <- TRUE
    }
    else if (config %in% c("agg.imp", "agg.cluster", "agg.all")) {
        abs <- x$print.options[["abs"]]
    }
    
    if (null.threshold) {
        if (which.stat=="Diff") {
            if (is_not_null(x$print.options$m.threshold)) {
                threshold <- x$print.options$m.threshold
                null.threshold <- FALSE
            }
        }
        else if (which.stat=="V.Ratio") {
            if (is_not_null(x$print.options$v.threshold)) {
                threshold <- x$print.options$v.threshold
                null.threshold <- FALSE
            }
        }
        else if (which.stat=="KS") {
            if (is_not_null(x$print.options$ks.threshold)) {
                threshold <- x$print.options$ks.threshold
                null.threshold <- FALSE
            }
        }
        else if (which.stat=="Corr") {
            if (is_not_null(x$print.options$r.threshold)) {
                threshold <- x$print.options$r.threshold
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
        
        co.names <- x[["print.options"]][["co.names"]]
        seps <- attr(co.names, "seps")
        for (i in names(co.names)) {
            comp <- co.names[[i]][["component"]]
            type <- co.names[[i]][["type"]]
            
            if (i %in% names(new.labels) && !is.na(new.labels[i])) {
                co.names[[i]][["component"]] <- new.labels[i]
                co.names[[i]][["type"]] <- "base"
            }
            else {
                if ("isep" %in% type) {
                    named.vars <- character(sum(type == "isep") + 1)
                    sep.inds <- c(which(type == "isep"), length(comp) + 1)
                    named.vars <- lapply(seq_along(sep.inds), function(k) {
                        inds <- (if (k == 1) seq(1, sep.inds[k] - 1) 
                                 else seq(sep.inds[k-1] + 1, sep.inds[k] - 1))
                        var <- comp[inds]
                        var.is.base <- type[inds] == "base"
                        pasted.var <- paste(var, collapse = "")
                        if (pasted.var %in% names(new.labels)) return(new.labels[pasted.var])
                        else return(paste(ifelse(var.is.base & var %in% names(new.labels) & !is.na(new.labels[var]), new.labels[var], var), collapse = ""))
                    })
                    co.names[[i]][["component"]] <- do.call("paste", c(unname(named.vars), list(sep = seps["int"])))
                }
                else co.names[[i]][["component"]] <- ifelse(type == "base" & comp %in% names(new.labels) & !is.na(new.labels[comp]), new.labels[comp], comp)
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
        if (x$print.options$nweights == 0) {
            ua <- c("Unadjusted", "Alphabetical")
            names(ua) <- c("unadjusted", "alphabetical")
        }
        else if (x$print.options$nweights == 1) {
            ua <- c("Adjusted", "Unadjusted", "Alphabetical")
            names(ua) <- c("adjusted", "unadjusted", "alphabetical")
        }
        else {
            ua <- c("Unadjusted", x$print.options$weight.names, "Alphabetical")
            names(ua) <- c("unadjusted", x$print.options$weight.names, "alphabetical")
        }
        var.order <- ua[match_arg(var.order, tolower(ua))]
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
        SS <- do.call("rbind", lapply(c("Un", x$print.options$weight.names),
                                      function(w) data.frame(var = B[["variable.names"]],
                                                             min.stat = B[[paste.("Min", which.stat, w)]],
                                                             max.stat = B[[paste.("Max", which.stat, w)]],
                                                             mean.stat = B[[paste.("Mean", which.stat, w)]],
                                                             Sample = ifelse(w == "Un", "Unadjusted", 
                                                                             ifelse(w == "Adj", "Adjusted", w)))))
        
        if (is_not_null(facet)) {
            if ("cluster" %in% facet) {
                SS$cluster <- rep(B[["cluster"]], 1 + x$print.options$nweights)
            }
            if ("imp" %in% facet) {
                SS$imp <- rep(B[["imp"]], 1 + x$print.options$nweights)
            }
            if ("treat.pair" %in% facet) {
                SS$treat.pair <- rep(B[["treat.pair"]], 1 + x$print.options$nweights)
            }
        }
        
        if (all(sapply(SS[c("min.stat", "max.stat", "mean.stat")], is.na))) stop(paste("No balance statistics to display. This can occur when", switch(which.stat, V.Ratio = "disp.v.ratio", KS = "disp.ks"), "= FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
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
        if (which.stat == "Diff" && any(abs(SS[["max.stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardized mean differences for continuous variables. To do so, specify continuous=\"std\" in bal.tab().", call.=FALSE, noBreaks.=TRUE)
        if (no.missing) SS <- SS[!is.na(SS[["min.stat"]]),]
        SS[["stat"]] <- SS[["mean.stat"]]
    }
    else {
        SS <- do.call("rbind", lapply(c("Un", x$print.options$weight.names),
                                      function(w) data.frame(var = B[["variable.names"]],
                                                             stat = B[[ifelse(is_null(Agg.Fun), paste.(which.stat, w),
                                                                              paste.(Agg.Fun, which.stat, w))]],
                                                             Sample = ifelse(w == "Un", "Unadjusted", 
                                                                             ifelse(w == "Adj", "Adjusted", w)))))
        
        if (is_not_null(facet)) {
            for (i in c("cluster", "imp", "treat.pair", "time")) {
                if (i %in% facet) {
                    SS[[i]] <- rep(B[[i]], 1 + x$print.options$nweights)
                }
            }
        }
        
        if (all(is.na(SS[["stat"]]))) stop(paste("No balance statistics to display. This can occur when", switch(which.stat, V.Ratio = "disp.v.ratio", KS = "disp.ks"), "= FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
        gone <- character(0)
        for (i in levels(SS$Sample)) {
            if (all(sapply(SS[["stat"]][SS$Sample==i], is.na))) {
                gone <- c(gone, i)
                if (i == "Unadjusted") warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
                SS <- SS[SS[["Sample"]]!=i,]
            }
        }
        
        if (abs) {
            if (which.stat == "V.Ratio") SS[["stat"]] <- pmax(abs(SS[["stat"]]), 1/abs(SS[["stat"]]))
            else SS[["stat"]] <- abs(SS[["stat"]])
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
        if (which.stat == "Diff" && any(abs(SS[["stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardized mean differences for continuous variables. To do so, specify continuous=\"std\" in bal.tab().", call.=FALSE, noBreaks.=TRUE)
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
        if (abs) {
            if (!null.threshold) threshold.xintercept <- max(threshold, 1/threshold)
        }
        else {
            if (!null.threshold) threshold.xintercept <- sort(c(threshold, 1/threshold))
        }
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
            f <- function(q) {q[["stat"]][q$var %in% distance.names] <- NA; q}
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
        
        if (is_null(subclass.names) || !x$print.options$disp.subclass) {
            if (line == TRUE) { #Add line except to distance
                f <- function(q) {q[["stat"]][q$var %in% distance.names] <- NA; q}
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
                f <- function(q) {q[["stat"]][q$var %in% distance.names] <- NA; q}
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
