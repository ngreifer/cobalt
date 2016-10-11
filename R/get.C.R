#Functions to turn input covariates into usable form
#int.poly.f creates interactions and polynomials
#split.factor splits factor variable into indicators
#binarize transforms 2-value variable into binary (0,1)
#get.C controls flow and handles redunancy

int.poly.f <- function(df, ex=NULL, int=FALSE, poly=1, nunder=1, ncarrot=1) {
    #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
    #Only to be used in base.bal.tab; for general use see int.poly()
    #df=data frame input
    #ex=data frame of variables to exclude in interactions and polynomials; a subset of df
    #int=whether to include interactions or not; currently only 2-way are supported
    #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
    #nunder=number of underscores between variables
    
    if (length(ex) > 0) d <- df[, !names(df) %in% names(ex), drop = FALSE]
    else d <- df
    nd <- ncol(d)
    nrd <- nrow(d)
    no.poly <- apply(d, 2, function(x) length(unique(x)) <= 2)
    npol <- nd - sum(no.poly)
    new <- matrix(ncol = (poly-1)*npol + .5*(nd)*(nd-1), nrow = nrd)
    nc <- ncol(new)
    new.names <- vector("character", nc)
    if (poly > 1 && npol != 0) {
        for (i in 2:poly) {
            new[, (1 + npol*(i - 2)):(npol*(i - 1))] <- apply(d[, !no.poly, drop = FALSE], 2, function(x) x^i)
            new.names[(1 + npol*(i - 2)):(npol*(i - 1))] <- paste0(colnames(d)[!no.poly], strrep("_", ncarrot), i)
        }
    }
    if (int && nd > 1) {
        new[,(nc - .5*nd*(nd-1) + 1):nc] <- matrix(t(apply(d, 1, combn, 2, prod)), nrow = nrd)
        new.names[(nc - .5*nd*(nd-1) + 1):nc] <- combn(names(d), 2, paste, collapse=strrep("_", nunder))
    }
    new <- data.frame(new)
    names(new) <- new.names
    single.value <- apply(new, 2, function(x) abs(max(x) - min(x)) < sqrt(.Machine$double.eps))
    new <- new[, !single.value, drop = FALSE]
    return(new)
}
split.factor <- function(varname, data) {
    #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
    #Retains all categories unless only 2 levels, in which case only the second level is retained.
    #varname= the name of the variable to split when data is specified
    #data=data set to be changed
    
    x <- factor(data[, varname])
    if (length(unique(x)) > 1) k <- model.matrix(~ as.character(as.matrix(x)) - 1)
    else {
        k <- matrix(0, ncol = nlevels(x), nrow = length(x))
        k[, levels(x) %in% unique(x)] <- 1
    }
    if (is.null(levels(x))) colnames(k) <- paste0(varname, 1:ncol(k))
    else colnames(k) <- paste0(varname, "_", substr(levels(x), 1, 10))
    
    if (match(varname, names(data)) == 1){
        data <- data.frame(k, data[, names(data)!=varname, drop = FALSE])
    }
    else if (match(varname, names(data)) == ncol(data)) {
        data <- data.frame(data[, names(data)!=varname, drop = FALSE], k)
    }
    else {
        where <- match(varname, names(data))
        data <- data.frame(data[, 1:(where-1), drop = FALSE], k, data[, (where+1):ncol(data), drop = FALSE])
    }
    return(data)
}
binarize <- function(variable) {
    if (length(unique(variable)) > 2) stop(paste0("Cannot binarize ", deparse(substitute(variable)), ": more than two levels."))
    if (0 %in% unique(as.numeric(variable))) zero <- 0
    #else if (1 %in% unique(as.numeric(variable))) zero <- unique(as.numeric(variable))[unique(as.numeric(variable)) != 1]
    else zero <- min(unique(as.numeric(variable)))
    variable <- ifelse(as.numeric(variable)==zero, 0, 1)
    return(variable)
}
get.C <- function(covs, int = FALSE, addl = NULL, distance = NULL, cluster = NULL) {
    #gets C data.frame, which contains all variables for which balance is to be assessed. Used in balance.table.
    C <- covs
    if (!is.null(addl)) {
        if (!is.data.frame(addl)) {
            if (is.character(addl)) stop("The argument to addl must be a data.frame containing the the values of the additional variables you want to include in the balance assessment.", call. = FALSE)
            else stop("The argument to addl must be a data.frame. Wrap data.frame() around the argument if it is a matrix or vector.", call. = FALSE)
        }
        else {
            repeat.name.indices <- sapply(names(addl), function(x) x %in% names(C))
            if (any(repeat.name.indices)) {
            warning(paste("The following variables in addl have the same name as covariates and will be ignored:\n",
                          paste(names(addl)[repeat.name.indices], collapse = " ")), call. = FALSE)
            addl <- addl[, !repeat.name.indices, drop = FALSE]
            }
            C <- cbind(C, addl)
        }
    } 
    
    for (i in names(C)) {
        if (is.character(C[, i])) C[, i] <- factor(C[, i])
        if (length(unique(C[, i])) <= 2) {
            if (is.logical(C[, i])) C[, i] <- as.numeric(C[, i])
            else if (is.numeric(C[, i])) C[, i] <- binarize(C[, i])
        }
        
        if (length(cluster) > 0 && qr(matrix(c(C[, i], as.numeric(cluster)), ncol = 2))$rank == 1) C <- C[, names(C) != i] #Remove variable if it is the same (linear combo) as cluster variable
        else if (!is.numeric(C[, i])) {
            C <- split.factor(i, C)
        }
    }

    if (int) {
        #Prevent duplicate var names with _'s
        nunder <- ncarrot <- 1
        repeat {
            if (all(sapply(names(C), function(x) !x %in% do.call(paste, c(expand.grid(names(C), names(C)), list(sep = strrep("_", nunder))))))) break
            else nunder <- nunder + 1
        }
        #Variable names don't contain carrots
        # repeat {
        #     if (all(sapply(names(C), function(x) !x %in% paste0(names(C), strrep("_", nunder), "2")))) break
        #     else ncarrot <- ncarrot + 1
        # }
        C <- cbind(C, int.poly.f(C, int = TRUE, poly = 2, nunder = nunder, ncarrot = ncarrot))
    }
    #Remove duplicate & redundant variables
    suppressWarnings(C.cor <- cor(C))
    s <- (1 - abs(C.cor) < .Machine$double.eps ^ .5) & !lower.tri(C.cor, diag=TRUE)
    redundant.vars <- apply(s, 2, function(x) any(x))
    C <- C[, !redundant.vars, drop = FALSE]        
    
    if (!is.null(distance)) {
        distance <- data.frame(.distance = distance)
        while (names(distance) %in% names(C)) {names(distance) <- paste0(names(distance), "_")}
        C <- cbind(distance, C)
        attr(C, "distance.name") <- names(distance)
    }
    
    return(C)

}
get.types <- function(C) {
    types <- mapply(function(x, y) ifelse(ifelse(is.null(attr(C, "distance.name")), FALSE, x==attr(C, "distance.name")), "Distance", ifelse(length(unique(y))<=2, "Binary", "Contin.")), colnames(C), C)
    return(types)
}