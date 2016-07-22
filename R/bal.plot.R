bal.plot <- function(obj, var.name, ..., un = FALSE, which.sub = NULL) {

    args <- list(...)

    if (any(class(obj)=="matchit")) X <- matchit2base(obj)
    else if (any(class(obj)=="ps")) X <- ps2base(obj, full.stop.method = args$full.stop.method)
    else if (any(class(obj)=="Match")) X <- Match2base(obj, formula=args$formula, data=args$data, treat=args$treat, covs=args$covs)
    else if (any(class(obj)=="CBPS")) X <- CBPS2base(obj, std.ok = TRUE)
    else if (any(class(obj)=="formula")) {
        X0 <- formula2df(obj, data = args$data)
        X <- df2base(X0$covs, X0$treat, data=args$data, weights=args$weights, distance=args$distance, subclass=args$subclass, method=args$method)
    }
    else if (is.data.frame(obj)) X <- df2base(obj, treat=args$treat, data=args$data, weights=args$weights, distance=args$distance, subclass=args$subclass, method=args$method)
    
    
    if (var.name %in% names(X$covs)) var <- X$covs[, var.name]
    else if (!is.null(args$data) && var.name %in% names(args$data)) var <- args$data[, var.name]
    else if (!is.null(X$addl) && var.name %in% names(X$addl)) var <- args$data[, var.name]
    else if (!is.null(args$addl) && var.name %in% names(args$addl)) var <- args$data[, var.name]
    else if (var.name==".distance" && !is.null(X$distance)) var <- X$distance
    else stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."))

    if (length(X$subclass)>0) {
        if (!is.null(which.sub)) {
            if (is.numeric(which.sub) && length(which.sub)==1) {
                if (which.sub %in% levels(X$subclass)) {
                    X$weights <- X$weights[!is.na(X$subclass) & X$subclass==which.sub]
                    X$treat <- X$treat[!is.na(X$subclass) & X$subclass==which.sub]
                    var <- var[!is.na(X$subclass) & X$subclass==which.sub]
                }
                else stop(paste0("\"", which.sub, "\" does not correspond to a subclass in the object."))

            }
            else stop("The argument to which.sub must be a single number corresponding to the subclass for which distributions are to be displayed.")
        }
        else stop("Argument contains subclasses but no which.sub value was supplied in \"...\".")
    }
    
    if (is.null(X$weights) || isTRUE(un)) X$weights <- rep(1, length(X$treat))
    weights <- ifelse(X$treat==0, X$weights / sum(X$weights[X$treat==0]), X$weights / sum(X$weights[X$treat==1]))
    treat <- factor(X$treat)
    
    #library(ggplot2)
    if (length(unique(var)) <= 2 || is.factor(var) || is.character(var)) {
        var <- factor(var)
        bp <- ggplot(mapping = aes(var, fill = treat, weight = weights)) + 
            geom_bar(position = "dodge", alpha = .4, color = "black") + 
            labs(x = var.name, y = "Proportion", fill = "Treat", title = paste0("Distributional Balance for \"", var.name, "\"")) 
    }
    else {
        bp <- ggplot(mapping=aes(var, fill=treat, weight=weights)) + 
            geom_density(alpha=.4) + 
            labs(x=var.name, fill="Treat", title=paste0("Distributional Balance for \"", var.name, "\""))
    }

    return(bp)
}