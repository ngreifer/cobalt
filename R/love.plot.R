love.plot <- function(x, stats = "mean.diffs", abs = TRUE, agg.fun = NULL, 
                      var.order = NULL, drop.missing = TRUE, drop.distance = FALSE, 
                      threshold = NULL, line = FALSE, stars = "none", grid = TRUE, 
                      limits = NULL, colors = NULL, shapes = NULL, alpha = 1, size = 3, 
                      wrap = 30, var.names = NULL, title, sample.names, labels = FALSE,
                      position = "right", themes = NULL, ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call(expand.dots = TRUE)
    .alls <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.all)), logical(1L))
    .nones <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.none)), logical(1L))
    if (any(c(.alls, .nones))) {
        .call[.alls] <- expression(NULL)
        .call[.nones] <- expression(NA)
        return(eval(.call))
    }
    
    #Re-call bal.tab with disp.v.ratio or disp.ks if stats = "v" or "k".
    if (!exists(deparse(substitute(x)))) { #if x is not an object (i.e., is a function call)
        mc <- match.call()
        replace.args <- function(m) {
            #m_ is bal.tab call or list (for do.call)
            if (any(sapply(stats, function(x) pmatch(x, "variance.ratios", 0L) != 0L))) {
                if (!isTRUE(eval(m[["disp.v.ratio"]]))) {
                    m[["un"]] <- TRUE
                    m[["disp.v.ratio"]] <- TRUE
                }
            }
            if (any(sapply(stats, function(x) pmatch(x, "ks.statistics", 0L) != 0L))) {
                if (!isTRUE(eval(m[["disp.ks"]]))) {
                    m[["un"]] <- TRUE
                    m[["disp.ks"]] <- TRUE
                }
            }
            
            if (any(names(m) == "cluster")) {
                m[["cluster.summary"]] <- TRUE
                if (any(names(m) == "cluster.fun")) m[["cluster.fun"]] <- NULL
            }
            if (any(names(m) == "imp")) {
                m[["imp.summary"]] <- TRUE
                if (any(names(m) == "imp.fun")) m[["imp.fun"]] <- NULL
            }
            
            m[["abs"]] <- abs
            
            return(m)
        }
        
        if (deparse(mc[["x"]][[1]]) %in% c("bal.tab", methods("bal.tab"))) { #if x i bal.tab call
            mc[["x"]] <- replace.args(mc[["x"]])
            x <- eval(mc[["x"]])
            
        }
        else if (deparse(mc[["x"]][[1]]) == "do.call") { #if x is do.call
            d <- match.call(eval(mc[["x"]][[1]]), mc[["x"]])
            if (deparse(d[["what"]]) %in% c("bal.tab", methods("bal.tab"))) {
                d[["args"]] <- replace.args(d[["args"]])
                x <- eval(d)
            }
        }
    }
    
    if ("bal.tab" %nin% class(x)) {
        #Use bal.tab on inputs first, then love.plot on that
        m <- match.call()
        m.b <- m; m.b[[1]] <- quote(bal.tab); names(m.b)[2] <- ''
        
        m.l <- m; 
        m.l[["x"]] <- m.b
        
        return(eval(m.l))
    }
    
    args <- list(...)
    
    if (any(class(x) == "bal.tab.cont")) {
        stats <- "correlations"
    }
    else stats <- match_arg(stats, c("mean.diffs", "variance.ratios", "ks.statistics"), several.ok = TRUE)
    
    which.stat <- c(mean.diffs = "Diff", variance.ratios = "V.Ratio", ks.statistics = "KS", correlations = "Corr")[stats]
    which.stat2 <- c(Diff = "Mean Difference", V.Ratio = "Variance Ratio", KS = "Kolmogorov-Smirnov Statistic", Corr = "Correlation")[which.stat]
    
    #shape (deprecated)
    #un.color (deprecated)
    #adj.color (deprecated)
    #cluster.fun (deprecated)
    #star_char
    
    p.ops <- c("which.cluster", "which.imp", "which.treat", "which.time", "disp.subclass")
    for (i in p.ops) {
        if (i %in% names(args)) attr(x, "print.options")[[i]] <- args[[i]]
    }
    
    #Using old argument names
    if (is_not_null(args$cluster.fun) && is_null(agg.fun)) agg.fun <- args$cluster.fun
    if (is_not_null(args$no.missing)) drop.missing <- args$no.missing
    
    Agg.Fun <- NULL
    subtitle <- NULL
    facet <- NULL
    
    cluster.names.good <- NULL
    imp.numbers.good <- NULL
    treat.names.good <- NULL; disp.treat.pairs <- NULL
    time.names.good <- NULL
    subclass.names <- NULL
    
    #NA = aggregate, NULL = individual
    null.cluster <- is_null(attr(x, "print.options")$which.cluster)
    na.cluster <- !null.cluster && is.na(attr(x, "print.options")$which.cluster)
    null.imp <- is_null(attr(x, "print.options")$which.imp)
    na.imp <- !null.imp && is.na(attr(x, "print.options")$which.imp)
    null.treat <- is_null(attr(x, "print.options")$which.treat)
    na.treat <- !null.treat && is.na(attr(x, "print.options")$which.treat)
    null.time <- is_null(attr(x, "print.options")$which.time)
    na.time <- !null.time && is.na(attr(x, "print.options")$which.time)
    
    config <- "agg.none"
    
    #Get B and config
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
        else if (is.character(attr(x, "print.options")$which.time)) {
            time.names.good <- sapply(time.names, function(x) any(sapply(attr(x, "print.options")$which.time, function(y) identical(x, y))))
            if (any(time.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.time are valid treatment names or time point indices.", call. = FALSE)
            }
        }
        else if (is.numeric(attr(x, "print.options")$which.time)) {
            time.names.good <- setNames(seq_along(x$Time.Balance) %in% attr(x, "print.options")$which.time, time.names)
            if (any(time.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.time are valid treatment names or time point indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.time must be .none, .all, or the names of treatments or indices of time points.", call. = FALSE)
        
        
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
                stop("Cannot aggregate across time periods without a balance summary across time periods.\nThis may be because multinomial treatments were used, multiple treatment types were used,\nor quick was set to TRUE and msm.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
            }
            #Agg.Fun <- switch(match_arg(agg.fun), mean = "Mean", max = "Max", range = "Range")
            Agg.Fun <- "Max"
            if (Agg.Fun == "Range") {
                subtitle <- "Range Across Time Points"
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Time Points")
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
        else if (is.numeric(attr(x, "print.options")$which.imp)) {
            imp.numbers.good <- setNames(imp.numbers %in% attr(x, "print.options")$which.imp, imp.numbers)
        }
        else stop("The argument to which.imp must be .none, .all, or the indices of imputations.", call. = FALSE)
        
        #Get cluster.names.good
        cluster.names <- names(x[["Cluster.Balance.Across.Imputations"]])
        if (na.cluster) {
            cluster.names.good <- NULL
        }
        else if (null.cluster) {
            cluster.names.good <- setNames(rep(TRUE, length(cluster.names)), cluster.names)
        }
        else if (is.character(attr(x, "print.options")$which.cluster)) {
            cluster.names.good <- sapply(cluster.names, function(n) any(sapply(attr(x, "print.options")$which.cluster, function(y) identical(n, y))))
        }
        else if (is.numeric(attr(x, "print.options")$which.cluster)) {
            cluster.names.good <- setNames(seq_along(x$Cluster.Balance) %in% attr(x, "print.options")$which.cluster, cluster.names)
        }
        else stop("The argument to which.cluster must be .none, .all, or the names or indices of clusters.", call. = FALSE)
        
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
                stop("At least one of which.cluster or which.imp must be .none or of length 1.", call. = FALSE)
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
            Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                              mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- "Range Across Imputations"
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Imputations", sep = " ")
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
            Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                              mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- "Range Across Clusters"
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Clusters", sep = " ")
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
            Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                              mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- "Range Across Clusters and Imputations"
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Clusters and Imputations", sep = " ")
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
        else if (is.numeric(attr(x, "print.options")$which.imp)) {
            imp.numbers.good <- setNames(imp.numbers %in% attr(x, "print.options")$which.imp, imp.numbers)
            if (any(imp.numbers.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.imp are valid imputation indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.imp must be .none, .all, or the indices of imputations.", call. = FALSE)
        
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
            Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                              mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- "Range Across Imputations"
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Imputations")
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
        else if (is.character(attr(x, "print.options")$which.cluster)) {
            cluster.names.good <- sapply(cluster.names, function(c) any(sapply(attr(x, "print.options")$which.cluster, function(y) identical(c, y))))
            if (any(cluster.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
            }
        }
        else if (is.numeric(attr(x, "print.options")$which.cluster)) {
            cluster.names.good <- setNames(seq_along(x$Cluster.Balance) %in% attr(x, "print.options")$which.cluster, cluster.names)
            if (any(cluster.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.cluster must be .none, .all, or the names or indices of clusters.", call. = FALSE)
        
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
            Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                              mean = "Mean", max = "Max", range = "Range")
            if (Agg.Fun == "Range") {
                subtitle <- "Range Across Clusters"
            }
            else {
                subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Clusters")
            }
            B <- cbind(x[["Cluster.Summary"]],
                       variable.names = rownames(x[["Cluster.Summary"]]))
        }
    }
    else if (any(class(x) == "bal.tab.multi")) {
        #Get cluster.names.good
        treat.names <- attr(x, "print.options")$treat.names
        if (na.treat) {
            config <- "agg.pair"
            treat.names.good <- NULL
        }
        else if (null.treat) {
            config <- "agg.none"
            treat.names.good <- rep(TRUE, length(treat.names))
        }
        else if (is.character(attr(x, "print.options")$which.treat)) {
            treat.names.good <- treat.names %in% attr(x, "print.options")$which.treat
            if (any(treat.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.treat are valid treatment names or indices.", call. = FALSE)
            }
        }
        else if (is.numeric(attr(x, "print.options")$which.treat)) {
            treat.names.good <- seq_along(treat.names) %in% attr(x, "print.options")$which.treat
            if (any(treat.names.good)) {
                config <- "agg.none"
            }
            else {
                stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
            }
        }
        else stop("The argument to which.cluster must be .none, .all, or the names or indices of clusters.", call. = FALSE)
        
        
        if (na.treat) {
            disp.treat.pairs <- character(0)
        }
        else {
            if (attr(x, "print.options")$pairwise) {
                if (sum(treat.names.good) == 1) {
                    disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) any(attr(x[["Pair.Balance"]][[p]], "print.options")$treat.names == attr(x, "print.options")$treat.names[treat.names.good]))]
                }
                else {
                    disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) all(attr(x[["Pair.Balance"]][[p]], "print.options")$treat.names %in% attr(x, "print.options")$treat.names[treat.names.good]))]
                }
            }
            else {
                if (sum(treat.names.good) == 1) {
                    disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) {
                        treat.names <- attr(x[["Pair.Balance"]][[p]], "print.options")$treat.namestreat.names
                        any(treat.names[treat.names != "Others"] == attr(x, "print.options")$treat.names[treat.names.good])})]
                }
                else {
                    disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) {
                        treat.names <- attr(x[["Pair.Balance"]][[p]], "print.options")$treat.namestreat.names
                        all(treat.names[treat.names != "Others"] %in% attr(x, "print.options")$treat.names[treat.names.good])})]
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
                subtitle <- paste0("Range Across Treatment", ifelse(attr(x, "print.options")$pairwise, " Pairs", "s"))
            }
            else {
                subtitle <- paste0(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), " Across Treatment", ifelse(attr(x, "print.options")$pairwise, " Pairs", "s"))
            }
            B <- cbind(x[["Balance.Across.Pairs"]],
                       variable.names = rownames(x[["Balance.Across.Pairs"]]))
        }
    }
    else if (any(class(x) == "bal.tab.subclass")) {
        if (any(stats == "variance.ratios")) stop("Variance ratios not currently supported for subclassification.", call. = FALSE)
        if (any(stats == "ks.statistics")) stop("KS statistics not currently supported for subclassification.", call. = FALSE)
        if (any(class(x) == "bal.tab.cont")) stop("Continuous treatments not currently supported for subclassification.", call. = FALSE)
        subclass.names <- names(x[["Subclass.Balance"]])
        sub.B <- do.call("cbind", lapply(subclass.names, function(s) {
            sub <- x[["Subclass.Balance"]][[s]]
            sub.B0 <- setNames(sub[endsWith(names(sub), ".Adj")],
                               gsub(".Adj", paste0(".Subclass ", s), names(sub)[endsWith(names(sub), ".Adj")]))
            return(sub.B0) }))
        B <- cbind(x[["Balance.Across.Subclass"]], sub.B, variable.names = row.names(x[["Balance.Across.Subclass"]]))
        if (attr(x, "print.options")$disp.subclass) attr(x, "print.options")$weight.names <- c("Adj", paste("Subclass", subclass.names))
        else attr(x, "print.options")$weight.names <- "Adj"
        subtitle <- "Across Subclasses"
    }
    else {
        B <- cbind(x[["Balance"]], variable.names = row.names(x[["Balance"]]))
    }
    
    if (is_not_null(facet) && length(stats) > 1) {
        stop("stats can only have a length of 1 when faceting by other dimension (e.g., cluster, treatment).", call. = FALSE)
    }
    
    if (is_not_null(agg.fun) && config == "agg.none") {
        warning("No aggregation will take place, so agg.fun will be ignored. Remember to set 'which.<ARG> = .none' to aggregate across <ARG>.", call. = FALSE)
    }
    
    #Process abs
    if (config == "agg.none") {
        if (!abs && !attr(x, "print.options")[["abs"]]) {
            abs <- FALSE
        }
        else {
            if (!abs && attr(x, "print.options")[["abs"]])
                warning("abs is TRUE in the bal.tab object; ignoring abs in the call to love.plot().", call. = FALSE)
            abs <- TRUE
        }
    }
    else if (config %in% c("agg.time", "agg.pair")) {
        abs <- TRUE
    }
    else if (config %in% c("agg.imp", "agg.cluster", "agg.all")) {
        abs <- attr(x, "print.options")[["abs"]]
    }
    
    #Process variable names
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
        
        co.names <- attr(x, "print.options")[["co.names"]]
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
    
    #Process variable order
    if (is_not_null(var.order) && "love.plot" %nin% class(var.order)) {
        if (attr(x, "print.options")$nweights == 0) {
            ua <- c("Unadjusted", "Alphabetical")
            names(ua) <- c("unadjusted", "alphabetical")
        }
        else if (attr(x, "print.options")$nweights == 1) {
            ua <- c("Adjusted", "Unadjusted", "Alphabetical")
            names(ua) <- c("adjusted", "unadjusted", "alphabetical")
        }
        else {
            ua <- c("Unadjusted", attr(x, "print.options")$weight.names, "Alphabetical")
            names(ua) <- c("unadjusted", attr(x, "print.options")$weight.names, "alphabetical")
        }
        var.order <- ua[match_arg(var.order, tolower(ua))]
    }
    
    #Process sample names
    ntypes <- if (is_null(subclass.names)) length(attr(x, "print.options")$weight.names) + 1 else 2
    if (!missing(sample.names)) {
        if (!is.vector(sample.names, "character")) {
            warning("The argument to sample.names must be a character vector. Ignoring sample.names.", call. = FALSE)
            sample.names <- NULL
        }
        else if (length(sample.names) %nin% c(ntypes, ntypes - 1)) {
            warning("The argument to sample.names must contain as many names as there are sample types, or one fewer. Ignoring sample.names.", call. = FALSE)
            sample.names <- NULL
        }
    }
    else sample.names <- NULL
    
    #Process limits
    if (is_not_null(limits)) {
        if (!is.vector(limits, "list")) {
            limits <- list(limits)
        }
        if (any(vapply(limits, 
                       function(l) !is.vector(l, "numeric") || length(l) %nin% c(0L, 2L), 
                       logical(1L)))) {
            warning("limits must be a list of numeric vectors of legnth 2. Ignoring limits.", call. = FALSE)
            limits <- NULL
        }
        
        if (is_not_null(names(limits))) {
            names(limits) <- stats[pmatch(names(limits), stats, duplicates.ok = TRUE)]
            limits <- limits[!is.na(names(limits))]
        }
        else {
            names(limits) <- stats[1:length(limits)]
        }
    }
    
    #Setting up appearance
    
    #Alpha (transparency)
    if (is.numeric(alpha[1]) && 
        !is.na(alpha[1]) && 
        between(alpha[1], c(0,1))) alpha <- alpha[1]
    else {
        warning("The argument to alpha must be a number between 0 and 1. Using 1 instead.", call. = FALSE)
        alpha <- 1
    }
    
    #Color
    if (is_not_null(args[["colours"]])) colors <- args[["colours"]]
    
    if (is_null(colors)) {
        if (shapes.ok(shapes, ntypes) && length(shapes) > 1 && length(shapes) == ntypes) {
            colors <- rep("black", ntypes)
        }
        else colors <- gg_color_hue(ntypes)
    }
    else {
        if (length(colors) == 1) colors <- rep(colors, ntypes)
        else if (length(colors) > ntypes) {
            colors <- colors[seq_len(ntypes)]
            warning(paste("Only using first", ntypes, "value", if (ntypes > 1) "s " else " ", "in colors."), call. = FALSE)
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
    # colors[] <- vapply(colors, col_plus_alpha, character(1L), alpha = alpha)
    fill <- colors
    
    #Shapes
    if (is_null(shapes)) {
        shapes <- assign.shapes(colors)
    }
    else {
        #check shapes
        if (!shapes.ok(shapes, ntypes)) {
            warning(paste("The argument to shape must be", ntypes, "valid shape", if (ntypes > 1) "s." else ".", " See ?love.plot for more information.\nUsing default shapes instead."), call. = FALSE)
            shapes <- assign.shapes(colors)
        }
        else if (length(shapes) == 1) shapes <- rep(shapes, ntypes)
    }
    
    #Size
    if (is.numeric(size[1])) size <- size[1]
    else {
        warning("The argument to size must be a number. Using 3 instead.", call. = FALSE)
        size <- 3
    }
    stroke0 <- stroke <- rep(0, ntypes)
    size0 <- size <- rep(size, ntypes)
    
    shapes.with.fill <- grepl("filled", shapes, fixed = TRUE)
    stroke[shapes.with.fill] <- size[shapes.with.fill]/3
    size[shapes.with.fill] <- size[shapes.with.fill]* .58
    
    # stroke <- .8*size
    
    if (is_not_null(facet)) {
        if (is_not_null(var.order) && "love.plot" %nin% class(var.order) && tolower(var.order) != "alphabetical" && (sum(cluster.names.good) > 1 || sum(imp.numbers.good) > 1 || length(disp.treat.pairs) > 1 || sum(time.names.good) > 1)) {
            warning("var.order cannot be set with faceted plots (unless \"alphabetical\"). Ignoring var.order.", call. = FALSE)
            var.order <- NULL
        }
    }
    
    agg.range <- is_not_null(Agg.Fun) && Agg.Fun == "Range"
    
    #Process thresholds
    stat2threshold <- c(mean.diffs = "m.threshold",
                        variance.ratios = "v.threshold",
                        ks.statistics = "ks.threshold",
                        correlations = "r.threshold")
    thresholds <- setNames(sapply(stats, function(i) {
        if (is_not_null(attr(x, "print.options")[[stat2threshold[i]]])) {
            return(attr(x, "print.options")[[stat2threshold[i]]])
        }
        else {
            return(NA_real_)
        }
    }), stats)
    
    if (is_not_null(threshold)) {
        if (!all(is.na(threshold)) && !is.numeric(threshold)) {
            stop("threshold must be numeric.", call. = FALSE)
        }
        
        if (is_not_null(names(threshold))) {
            names(threshold) <- stats[pmatch(names(threshold), stats, duplicates.ok = TRUE)]
            threshold <- threshold[!is.na(names(threshold))]
        }
        else {
            names(threshold) <- stats[1:length(threshold)]
        }
        
        thresholds[names(threshold)] <- as.numeric(threshold)
    }
    thresholds <- thresholds[!is.na(thresholds)]
    
    #Title
    if (missing(title)) title <- "Covariate Balance"
    else title <- as.character(title)
    # if (missing(subtitle)) subtitle <- as.character(subtitle)
    
    #Process themes
    if (is_not_null(themes)) {
        if (!is.vector(themes, "list")) {
            themes <- list(themes)
        }
        if (any(vapply(themes, 
                       function(t) !all(c("theme", "gg") %in% class(t)), 
                       logical(1L)))) {
            warning("themes must be a list of \"theme\" objects. Ignoring themes.", call. = FALSE)
            themes <- NULL
        }
        
        if (is_not_null(names(themes))) {
            names(themes) <- stats[pmatch(names(themes), stats, duplicates.ok = TRUE)]
            themes <- themes[!is.na(names(themes))]
        }
        else {
            names(themes) <- stats[1:length(themes)]
        }
    }
    
    plot.list <- setNames(vector("list", length(stats)), stats)
    for (s in stats) {
        variable.names <- as.character(B[["variable.names"]])
        if (s == "mean.diffs") {
            binary <- attr(x, "print.options")$binary
            continuous <- attr(x, "print.options")$continuous
            #All std, no std, some std
            if ((binary == "std" || sum(B[["Type"]] == "Binary") == 0) && 
                (continuous == "std" || sum(B[["Type"]] != "Binary") == 0)) {
                xlab.diff <- "Standardized Mean Differences"
            } 
            else if ((binary == "raw" || sum(B[["Type"]] == "Binary") == 0) && 
                     (continuous == "raw" || sum(B[["Type"]] != "Binary") == 0)) {
                xlab.diff <- "Mean Differences"
            }
            else {
                stars <- match_arg(stars, c("none", "std", "raw"))
                if (stars == "none") {
                    warning("Standardized mean differences and raw mean differences are present in the same plot. \nUse the 'stars' argument to distinguish between them and appropriately label the x-axis.", call. = FALSE)
                    xlab.diff <- "Mean Differences"
                }
                else {
                    if (length(args$star_char) == 1 && is.character(args$star_char)) star_char <- args$star_char
                    else star_char <- "*"
                    
                    vars_to_star <- setNames(rep(FALSE, nrow(B)), variable.names)
                    if (stars == "std") {
                        if (binary == "std") vars_to_star[variable.names[B[["Type"]] == "Binary"]] <- TRUE
                        if (continuous == "std") vars_to_star[variable.names[B[["Type"]] != "Binary"]] <- TRUE
                        xlab.diff <- "Mean Differences"
                    }
                    else if (stars == "raw") {
                        if (binary == "raw") vars_to_star[variable.names[B[["Type"]] == "Binary"]] <- TRUE
                        if (continuous == "raw") vars_to_star[variable.names[B[["Type"]] != "Binary"]] <- TRUE
                        xlab.diff <- "Standardized Mean Differences"
                    }
                    variable.names[vars_to_star[variable.names]] <- paste0(variable.names[vars_to_star[variable.names]], star_char)
                }
            }
        }
        
        #Get SS
        if (agg.range) {
            SS <- do.call("rbind", 
                          lapply(c("Un", attr(x, "print.options")$weight.names),
                                 function(w) data.frame(var = variable.names,
                                                        type = B[["Type"]],
                                                        min.stat = B[[paste.("Min", which.stat[s], w)]],
                                                        max.stat = B[[paste.("Max", which.stat[s], w)]],
                                                        mean.stat = B[[paste.("Mean", which.stat[s], w)]],
                                                        Sample = ifelse(w == "Un", "Unadjusted", 
                                                                        ifelse(w == "Adj", "Adjusted", w)),
                                                        row.names = NULL)))
            
            if (is_not_null(facet)) {
                if ("cluster" %in% facet) {
                    SS$cluster <- rep(B[["cluster"]], 1 + attr(x, "print.options")$nweights)
                }
                if ("imp" %in% facet) {
                    SS$imp <- rep(B[["imp"]], 1 + attr(x, "print.options")$nweights)
                }
                if ("treat.pair" %in% facet) {
                    SS$treat.pair <- rep(B[["treat.pair"]], 1 + attr(x, "print.options")$nweights)
                }
            }
            
            
            if (all(sapply(SS[c("min.stat", "max.stat", "mean.stat")], is.na))) 
                stop(paste("No balance statistics to display. This can occur when", 
                           switch(s, variance.ratios = "disp.v.ratio", ks.statistics = "disp.ks"), 
                           "= FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            
            missing.stat <- all(is.na(SS[["mean.stat"]]))
            if (missing.stat) stop(paste0(word_list(firstup(tolower(which.stat2))), 
                                          " cannot be displayed. This can occur when ", 
                                          word_list(paste.("disp", tolower(which.stat[s])), and.or = "and", is.are = TRUE), 
                                          " FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            
            gone <- character(0)
            for (i in levels(SS$Sample)) {
                if (all(sapply(SS[SS$Sample == i, c("min.stat", "max.stat", "mean.stat")], is.na))) {
                    gone <- c(gone, i)
                    if (i == "Unadjusted") warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
                    SS <- SS[SS[["Sample"]] != i,]
                }
            }
            
            dec <- FALSE
            
            if (is_not_null(var.order)) {
                if ("love.plot" %in% class(var.order)) {
                    old.vars <- levels(var.order$data$var)
                    old.vars[endsWith(old.vars, "*")] <- gsub("*", "", old.vars[endsWith(old.vars, "*")], fixed = TRUE)
                    if (any(SS[["var"]] %nin% old.vars)) {
                        warning("The love.plot object in var.order doesn't have the same variables as the current input. Ignoring var.order.", call. = FALSE)
                        var.order <- NULL
                    }
                    else {
                        SS[["var"]] <- factor(SS[["var"]], levels = old.vars[old.vars %in% SS[["var"]]])
                    }
                }
                else if (tolower(var.order) == "alphabetical") {
                    if ("time" %in% facet) {
                        covnames0 <- vector("list", length(unique(SS[["time"]])))
                        for (i in seq_along(covnames0)) {
                            if (i == 1) {
                                covnames0[[i]] <- sort(levels(SS[["var"]][SS[["time"]] == i]))
                            }
                            else {
                                covnames0[[i]] <- sort(setdiff(levels(SS[["var"]][SS[["time"]] == i]), unlist(covnames0[seq_along(covnames0) < i])))
                            }
                        }
                        covnames <- unlist(covnames0)
                    }
                    else covnames <- sort(levels(SS[["var"]]))
                    SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
                    
                }
                else if (var.order %in% ua) {
                    if (var.order %in% gone) {
                        warning(paste0("var.order was set to \"", tolower(var.order), "\", but no ", tolower(var.order), " ", tolower(which.stat2), "s were calculated. Ignoring var.order."), call. = FALSE, immediate. = TRUE)
                        var.order <- NULL
                    }
                    else {
                        v <- as.character(SS[["var"]][order(SS[["stat"]][SS[["Sample"]]==var.order], decreasing = dec, na.last = FALSE)])
                        
                        SS[["var"]] <- factor(SS[["var"]], 
                                              levels=c(v[v %nin% distance.names], 
                                                       sort(distance.names, decreasing = TRUE)))
                    }
                }
                
            }
            if (is_null(var.order)) {
                covnames <- as.character(unique(SS[["var"]]))
                SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
            }
            # SS[, "Sample"] <- factor(SS[, "Sample"], levels = c("Adjusted", "Unadjusted"))
            SS[["Sample"]] <- factor(SS[["Sample"]])
            if (s == "mean.diffs" && any(base::abs(SS[["max.stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardized mean differences for continuous variables.", call.=FALSE)
            if (length(stats) == 1 && drop.missing) SS <- SS[!is.na(SS[["min.stat"]]),]
            SS[["stat"]] <- SS[["mean.stat"]]
        }
        else {
            SS <- do.call("rbind", 
                          lapply(c("Un", attr(x, "print.options")$weight.names),
                                 function(w) data.frame(var = variable.names,
                                                        type = B[["Type"]],
                                                        stat = B[[ifelse(is_null(Agg.Fun), paste.(which.stat[s], w),
                                                                         paste.(Agg.Fun, which.stat[s], w))]],
                                                        Sample = ifelse(w == "Un", "Unadjusted", 
                                                                        ifelse(w == "Adj", "Adjusted", w)),
                                                        row.names = NULL)))
            
            if (is_not_null(facet)) {
                for (i in c("cluster", "imp", "treat.pair", "time")) {
                    if (i %in% facet) {
                        SS[[i]] <- rep(B[[i]], 1 + attr(x, "print.options")$nweights)
                    }
                }
            }
            
            missing.stat <- all(is.na(SS[["stat"]]))
            if (missing.stat) stop(paste0(word_list(firstup(tolower(which.stat2))), 
                                          " cannot be displayed. This can occur when ", 
                                          word_list(paste.("disp", tolower(which.stat[s])), and.or = "and"), 
                                          " are FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            
            gone <- character(0)
            for (i in levels(SS$Sample)) {
                if (all(is.na(SS[["stat"]][SS$Sample==i]))) {
                    gone <- c(gone, i)
                    if (i == "Unadjusted") warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
                    SS <- SS[SS[["Sample"]]!=i,]
                }
            }
            
            if (abs) {
                if (s == "variance.ratios") SS[["stat"]] <- pmax(SS[["stat"]], 1/SS[["stat"]])
                else if (s == "mean.diffs") SS[["stat"]] <- base::abs(SS[["stat"]])
            }
            dec <- FALSE
            
            if (is_not_null(plot.list[[1]])) var.order <- plot.list[[1]]
            
            if (is_not_null(var.order)) {
                if ("love.plot" %in% class(var.order)) {
                    old.vars <- levels(var.order$data$var)
                    old.vars[endsWith(old.vars, "*")] <- gsub("*", "", old.vars[endsWith(old.vars, "*")], fixed = TRUE)
                    if (any(SS[["var"]] %nin% old.vars)) {
                        warning("The love.plot object in var.order doesn't have the same variables as the current input. Ignoring var.order.", call. = FALSE)
                        var.order <- NULL
                    }
                    else {
                        SS[["var"]] <- factor(SS[["var"]], levels = old.vars[old.vars %in% SS[["var"]]])
                    }
                }
                else if (tolower(var.order) == "alphabetical") {
                    if ("time" %in% facet) {
                        covnames0 <- vector("list", length(unique(SS[["time"]])))
                        for (i in seq_along(covnames0)) {
                            if (i == 1) {
                                covnames0[[i]] <- sort(levels(SS[["var"]][SS[["time"]] == i]))
                            }
                            else {
                                covnames0[[i]] <- sort(setdiff(levels(SS[["var"]][SS[["time"]] == i]), unlist(covnames0[seq_along(covnames0) < i])))
                            }
                        }
                        covnames <- unlist(covnames0)
                    }
                    else covnames <- sort(levels(SS[["var"]]))
                    SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
                    
                }
                else if (var.order %in% ua) {
                    if (var.order %in% gone) {
                        warning(paste0("var.order was set to \"", tolower(var.order), "\", but no ", tolower(var.order), " ", tolower(which.stat2), "s were calculated. Ignoring var.order."), call. = FALSE, immediate. = TRUE)
                        var.order <- NULL
                    }
                    else {
                        v <- as.character(SS[["var"]][order(SS[["stat"]][SS[["Sample"]]==var.order], decreasing = dec, na.last = FALSE)])
                        
                        SS[["var"]] <- factor(SS[["var"]], 
                                              levels=c(v[v %nin% distance.names], 
                                                       sort(distance.names, decreasing = TRUE)))
                    }
                }
                
            }
            if (is_null(var.order)) {
                covnames <- as.character(unique(SS[["var"]])) #Don't use levels here to preserve original order
                SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[covnames %nin% distance.names]), sort(distance.names, decreasing = TRUE)))
            }
            SS[["Sample"]] <- factor(SS[["Sample"]])
            if (s == "mean.diffs" && any(base::abs(SS[["stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardized mean differences for continuous variables.", call.=FALSE)
            if (length(stats) == 1 && drop.missing) SS <- SS[!is.na(SS[["stat"]]),]
        }
        
        SS <- SS[order(SS[["var"]], na.last = FALSE),]
        SS[["var"]] <- factor(SS[["var"]])
        
        #Make the plot
        #library(ggplot2)
        
        baseline.xintercept <- switch(s, 
                                      "mean.diffs" = 0, 
                                      "variance.ratios" = 1, 
                                      "ks.statistics" = 0, 
                                      "correlations" = 0)
        threshold.xintercepts <- NULL
        if (s == "correlations") {
            if (abs) {
                xlab <- "Absolute Treatment-Covariate Correlations"
                if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = base::abs(thresholds[s]))
            }
            else {
                xlab <- "Treatment-Covariate Correlations"
                if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = -thresholds[s], upper = thresholds[s])
            }
            scale_Statistics <- scale_x_continuous
        }
        else if (s == "mean.diffs") {
            if (abs) {
                xlab <- paste("Absolute", xlab.diff)
                if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = base::abs(thresholds[s]))
            }
            else {
                xlab <- xlab.diff
                if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = -thresholds[s], upper = thresholds[s])
            }
            scale_Statistics <- scale_x_continuous
        }
        else if (s == "variance.ratios") {
            xlab <- "Variance Ratios"
            if (abs) {
                if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = max(thresholds[s], 1/thresholds[s]))
            }
            else {
                if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = min(c(thresholds[s], 1/thresholds[s])), upper = max(c(thresholds[s], 1/thresholds[s])))
            }
            scale_Statistics <- scale_x_log10
        }
        else if (s == "ks.statistics") {
            xlab <- "Kolmogorov-Smirnov Statistics"
            if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = base::abs(thresholds[s]))
            scale_Statistics <- scale_x_continuous
        }
        
        apply.limits <- FALSE
        SS[["on.border"]] <- FALSE
        if (is_not_null(limits[[s]])) {
            if (limits[[s]][2] < limits[[s]][1]) {
                limits[[s]] <- c(limits[[s]][2], limits[[s]][1])
            }
            
            if (limits[[s]][1] >= baseline.xintercept) limits[[s]][1] <- baseline.xintercept - .05*limits[[s]][2]
            if (limits[[s]][2] <= baseline.xintercept) limits[[s]][2] <- baseline.xintercept - .05*limits[[s]][1]
            
            if (identical(scale_Statistics, scale_x_log10)) limits[[s]][limits[[s]] <= 1e-2] <- 1e-2
            
            if (agg.range) {
                # for (i in c("min.stat", "mean.stat", "max.stat")) {
                #     if (any(SS[[i]] < limits[[s]][1], na.rm = TRUE)) {
                #         if (i == "mean.stat") SS[["on.border"]][SS[[i]] < limits[[s]][1]] <- TRUE
                #         if (i == "mean.stat") SS[[i]][SS[[i]] < limits[[s]][1]] <- limits[[s]][1]
                #         else SS[[i]][SS[[i]] < limits[[s]][1]] <- limits[[s]][1]
                #     }
                #     if (any(SS[[i]] > limits[[s]][2], na.rm = TRUE)) {
                #         if (i == "mean.stat") SS[["on.border"]][SS[[i]] > limits[[s]][2]] <- TRUE
                #         if (i == "mean.stat") SS[[i]][SS[[i]] > limits[[s]][2]] <- limits[[s]][2]
                #         else SS[[i]][SS[[i]] > limits[[s]][2]] <- limits[[s]][2]
                #         
                #         # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                #     }
                # }
                
                if (any(SS[["mean.stat"]] < limits[[s]][1], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["stat"]] < limits[[s]][1]] <- TRUE
                    SS[["mean.stat"]][SS[["mean.stat"]] < limits[[s]][1]] <- limits[[s]][1]
                    SS[["max.stat"]][SS[["max.stat"]] < limits[[s]][1]] <- limits[[s]][1]
                    SS[["min.stat"]][SS[["min.stat"]] < limits[[s]][1]] <- limits[[s]][1]
                }
                if (any(SS[["mean.stat"]] > limits[[s]][2], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["stat"]] > limits[[s]][2]] <- TRUE
                    SS[["mean.stat"]][SS[["mean.stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    SS[["max.stat"]][SS[["max.stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    SS[["min.stat"]][SS[["min.stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
                # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
            }
            else {
                if (any(SS[["stat"]] < limits[[s]][1], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["stat"]] < limits[[s]][1]] <- TRUE
                    SS[["stat"]][SS[["stat"]] < limits[[s]][1]] <- limits[[s]][1]
                }
                if (any(SS[["stat"]] > limits[[s]][2], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["stat"]] > limits[[s]][2]] <- TRUE
                    SS[["stat"]][SS[["stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
            }
            
            apply.limits <- TRUE
        }
        
        if (is_not_null(sample.names)) {
            if (length(sample.names) == ntypes - 1) {
                levels(SS$Sample)[-1] <- sample.names
            }
            else if (length(sample.names) == ntypes) {
                levels(SS$Sample) <- sample.names
            }
        }
        
        lp <- ggplot(aes(y = var, x = stat, group = Sample), data = SS) +
            theme(panel.background = element_rect(fill = "white"),
                  axis.text.x = element_text(color = "black"),
                  axis.text.y = element_text(color = "black"),
                  panel.border = element_rect(fill = NA, color = "black"),
                  plot.background = element_blank(),
                  legend.background = element_blank(),
                  legend.key = element_blank()
            ) +
            scale_shape_manual(values = shapes) +
            scale_size_manual(values = size) +
            scale_discrete_manual(aesthetics = "stroke", values = stroke) +
            scale_fill_manual(values = fill) +
            scale_color_manual(values = colors) +
            labs(y = NULL, x = wrap(xlab, wrap))
        
        lp <- lp + geom_vline(xintercept = baseline.xintercept,
                              linetype = 1, color = "gray5")
        
        if (is_not_null(threshold.xintercepts)) {
            lp <- lp + geom_vline(xintercept = threshold.xintercepts,
                                  linetype = 2, color = "gray8")
        }
        
        if (agg.range) {
            position.dodge <- ggstance::position_dodgev(.5*(size0[1]/3))
            if (line == TRUE) { #Add line except to distance
                f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
                lp <- lp + ggplot2::layer(geom = "path", data = f, 
                                          position = position.dodge, 
                                          stat = "identity",
                                          mapping = aes(x = mean.stat, color = Sample), 
                                          params = list(size = size0[1]*.8/3, na.rm = TRUE,
                                                        alpha = alpha))
            }
            
            lp <- lp +
                geom_errorbarh(aes(y = var, xmin = min.stat, xmax = max.stat,
                                   color = Sample), position = position.dodge,
                               size = size0[1]*.8/3,
                               alpha = alpha, 
                               height = 0,
                               na.rm = TRUE) +
                geom_point(aes(y = var, 
                               x = mean.stat, 
                               shape = Sample,
                               size = Sample,
                               stroke = Sample,
                               color = Sample),
                           fill = "white", na.rm = TRUE,
                           alpha = alpha,
                           position = position.dodge)
            
        }
        else {
            if (is_null(subclass.names) || !attr(x, "print.options")$disp.subclass) {
                if (line == TRUE) { #Add line except to distance
                    f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
                    lp <- lp + ggplot2::layer(geom = "path", data = f(SS),
                                              position = "identity", stat = "identity",
                                              mapping = aes(color = Sample),
                                              params = list(size = size0[1]*.8/3,
                                                            na.rm = TRUE,
                                                            alpha = alpha))
                }
                lp <- lp + geom_point(data = SS, aes(shape = Sample,
                                                     size = Sample,
                                                     stroke = Sample,
                                                     color = Sample),
                                      fill = "white", 
                                      na.rm = TRUE,
                                      alpha = alpha)
                
            }
            else {
                SS.u.a <- SS[SS$Sample %in% c("Unadjusted", "Adjusted"),]
                SS.u.a$Sample <- factor(SS.u.a$Sample)
                if (line == TRUE) { #Add line except to distance
                    f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
                    lp <- lp + ggplot2::layer(geom = "path", data = f(SS.u.a),
                                              position = "identity", stat = "identity",
                                              mapping = aes(color = Sample),
                                              params = list(size = size*.8,
                                                            na.rm = TRUE,
                                                            alpha = alpha))
                }
                lp <- lp + geom_point(data = SS.u.a,
                                      aes(shape = Sample,
                                          size = Sample,
                                          stroke = Sample,
                                          color = Sample),
                                      fill = "white",
                                      na.rm = TRUE)
                lp <- lp + geom_text(data = SS[SS$Sample %nin% c("Unadjusted", "Adjusted"),],
                                     mapping = aes(label = gsub("Subclass ", "", Sample)),
                                     size = 2.5*size0[1]/3, na.rm = TRUE)
            }
            
            
        }
        
        if (!drop.distance && is_not_null(distance.names)) {
            lp <- lp + geom_hline(linetype = 1, color = "black",
                                  yintercept = nunique(SS[["var"]]) - length(distance.names) + .5)
        }
        if (apply.limits) {
            lp <- lp + scale_Statistics(limits = limits[[s]], expand = c(0, 0))
        }
        else {
            lp <- lp + scale_Statistics()
        }
        
        if (isFALSE(grid)) {
            lp <- lp + theme(panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())
        }
        else {
            lp <- lp + theme(panel.grid.major = element_line(color = "gray87"),
                             panel.grid.minor = element_line(color = "gray90"))
        }
        
        if (is_not_null(facet)) {
            lp <- lp + facet_grid(f.build(".", facet), drop = FALSE) + labs(x = xlab)
        }
        
        class(lp) <- c(class(lp), "love.plot")
        plot.list[[s]] <- lp
    }
    
    if (length(stats) > 1 || isTRUE(args$use.grid)) {
        
        if (!is.character(position) || length(position) > 1) {
            position <- NA_character_
        }
        position <- match_arg(as.character(position), 
                              c("right", "left", "top", "bottom", "none"))
        
        #Process labels
        if (isTRUE(labels)) labels <- LETTERS[seq_along(plot.list)]
        else if (is_null(labels) || isFALSE(labels)) labels <- NULL
        else if (!is.atomic(labels) || length(labels) != length(plot.list)) {
            warning("labels must be TRUE or a string with the same length as stats. Ignoring labels.", call. = FALSE)
            labels <- NULL
        }
        else labels <- as.character(labels)
        
        # p <- ggpubr::ggarrange(plotlist = plot.list, common.legend = TRUE, legend = position, 
        #                align = "hv", nrow = 1)
        # if (is_not_null(subtitle)) {
        #     p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(subtitle, size = 11))
        # }
        # p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(title, size = 13.2))
        # 
        # P <- attr(P, "plots")
        
        plots.to.combine <- plot.list
        for (i in seq_along(plots.to.combine)) {
            if (i > 1) {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + 
                    theme(axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          legend.position = "none")
            }
            else {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + theme(legend.position = "none")
            }
            
            if (is_not_null(labels)) {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + labs(title = labels[i])
            }
            
            if (is_not_null(themes[[stats[i]]])) {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + themes[[stats[i]]]
            }
        }
        
        g <- ggarrange_simple(plots = plots.to.combine, nrow = 1)
        title.grob <- grid::textGrob(title, gp=grid::gpar(fontsize=13.2))
        subtitle.grob <- grid::textGrob(subtitle, gp=grid::gpar(fontsize=13.2))
        
        if (position == "none") {
            p <- gridExtra::arrangeGrob(grobs = list(g), nrow = 1)
        }
        else {
            legg <- ggplot2::ggplotGrob(plots.to.combine[[1]] + theme(legend.position = position))
            leg <- legg$grobs[[which(legg$layout$name == "guide-box")]]
            
            if (position == "left") {
                p <- gridExtra::arrangeGrob(grobs = list(leg, g), nrow = 1, 
                                            widths = grid::unit.c(sum(leg$widths), grid::unit(1, "npc") - sum(leg$widths)))
            }
            else if (position == "right") {
                p <- gridExtra::arrangeGrob(grobs = list(g, leg), nrow = 1, 
                                            widths = grid::unit.c(grid::unit(1, "npc") - sum(leg$widths), sum(leg$widths)))
            }
            else if (position == "top") {
                p <- gridExtra::arrangeGrob(grobs = list(leg, g), nrow = 2,
                                            heights = grid::unit.c(sum(leg$heights), grid::unit(1, "npc") - sum(leg$heights)))
            }
            else if (position == "bottom") {
                p <- gridExtra::arrangeGrob(grobs = list(g, leg), nrow = 2,
                                            heights = grid::unit.c(grid::unit(1, "npc") - sum(leg$heights), sum(leg$heights)))
            }
        }
        
        if (is_not_null(subtitle)) {
            p <- gridExtra::arrangeGrob(p, top = subtitle.grob)
        }
        p <- gridExtra::arrangeGrob(p, top = title.grob)
        
        grid::grid.newpage()
        grid::grid.draw(p)
        
        attr(p, "plots") <- plot.list
        class(p) <- c(class(p), "love.plot")
        
        return(invisible(p))
    }
    else {
        
        p <- plot.list[[1]] + 
            labs(title = title, subtitle = subtitle) +
            theme(plot.title = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5),
                  legend.position = position)
        
        if (is_not_null(themes[[1]])) {
            p <- p + themes[[1]]
        }
        
        return(p)
        
    }
    
}
plot.bal.tab <- love.plot
