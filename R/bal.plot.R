bal.plot <- function(obj, var.name, ...) {
    #Produces distribution plot of treated an untreated groups after adjustment. 
    #obj is the first input to bal.tab, var.name is the quoted name of the variable to be assessed.
    #
    #ps: full.stop.method (optional, will use first if blank)
    #Match: formula and data or treat and covs
    #formula: data, weights, distance (optional), method (optional, will default to matching)
    #data.frame: treat, data (optional; for naming treat and weights), weights, distance (optional), method (optional, will default to matching)
    
    args <- list(...)
    if (any(class(obj)=="matchit")) X <- matchit2base(obj)
    else if (any(class(obj)=="ps")) X <- ps2base(obj, full.stop.method = args$full.stop.method)
    else if (any(class(obj)=="Match")) X <- Match2base(obj, formula=args$formula, data=args$data, treat=args$treat, covs=args$covs)
    else if (any(class(obj)=="CBPS")) X <- CBPS2base(obj, std.ok = TRUE)
    else if (any(class(obj)=="formula")) {
        X0 <- formula2df(obj, data = args$data)
        X <- df2base(X0$covs, X0$treat, data=args$data, weights=args$weights, distance=args$distance, method=args$method)
    }
    else if (is.data.frame(obj)) X <- df2base(obj, treat=args$treat, data=args$data, weights=args$weights, distance=args$distance, method=args$method)
    
    
    if (var.name %in% names(X$covs)) var <- X$covs[, var.name]
    else if (!is.null(args$data) & var.name %in% names(args$data)) var <- args$data[, var.name]
    else if (!is.null(X$addl) & var.name %in% names(X$addl)) var <- args$data[, var.name]
    else if (!is.null(args$addl) & var.name %in% names(args$addl)) var <- args$data[, var.name]
    else stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."))
    
    if (!is.null(X$subclass)) {
        if (is.numeric(args$which.sub)) {
            X$weights <- X$weights[X$subclass==args$which.sub]
            X$treat <- X$treat[X$subclass==args$which.sub]
            var <- var[X$subclass==args$which.sub]
        }
        else stop("Argument contains subclasses but no which.sub value was supplied in \"...\".")
    }
    
    X$weights <- ifelse(X$treat==0, X$weights / sum(X$weights[X$treat==0]), X$weights / sum(X$weights[X$treat==1]))
    X$treat <- factor(X$treat)
    
    #library(ggplot2)
    if (length(unique(var)) > 2) {
        bp <- ggplot(mapping=aes(var, fill=X$treat, weight=X$weights)) + 
            geom_density(alpha=.4) + 
            labs(x=var.name, fill="Treat", title=paste0("Distributional Balance for \"", var.name, "\""))
    }
    else {
        var <- factor(var)
        bp <- ggplot(mapping = aes(var, fill = X$treat, weight = X$weights)) + 
            geom_bar(position = "dodge", alpha = .4, color = "black") + 
            labs(x = var.name, y = "Proportion", fill = "Treat", title = paste0("Distributional Balance for \"", var.name, "\"")) 
    }
    return(bp)
}