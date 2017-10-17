int.poly <- function(df, int=FALSE, poly=1, ex=NULL, nunder=1) {
    #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
    #Only to be used in base.bal.tab; for general use see int.poly()
    #df=data frame input
    #ex=data frame of variables to exclude in interactions and polynomials; a subset of df
    #int=whether to include interactions or not; currently only 2-way are supported
    #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
    #nunder=number of underscores between variables
    
    if (length(ex) > 0) {
        if (is.character(ex)) d <- df[, !names(df) %in% ex, drop = FALSE]
        else if (is.data.frame(ex) || is.matrix(ex)) d <- df[, !names(df) %in% names(ex), drop = FALSE]
        else stop("The argument to ex must be a vector of variable names to exclude or a data frame or matrix with named variables to exclude.")
    }
    else d <- df
    are.factors <- sapply(d, function(x) !is.numeric(x) || !is.logical(x))
    if (any(are.factors)) stop(paste0("The following variables in df are non-numeric:\n", paste(names(df)[are.factors], collapse = " "), 
                                      "\nPlease exclude them or convert them into numeric or logical."))
    nd <- ncol(d)
    binary <- apply(d, 2, function(x) length(unique(x)) <= 2) #No polynomials for binary vars
    
    new <- matrix(ncol = (poly-1)*(nd - sum(binary)) + .5*(nd)*(nd-1), nrow = nrow(d))
    new.names <- vector("character", ncol(new))
    if (poly > 1 && sum(no.poly) < nd) {
        for (i in 2:poly) {
            new.poly <- apply(d[, !binary, drop = FALSE], 2, function(x) x^i)
            new[, (1 + ncol(new.poly)*(i - 2)):(ncol(new.poly)*(i - 1))] <- new.poly
            new.names[(1 + ncol(new.poly)*(i - 2)):(ncol(new.poly)*(i - 1))] <- paste0(colnames(d)[!binary], strrep("_", nunder), i)
        }
    }
    if (int && nd > 1) {
        non.0.1.binary <- sapply(d, function(x) (binary && any(!unique(x) %in% c(0,1))))
        if (any(non.0.1.binary)) warning(paste0("The following variables are binary but have values different from 0 or 1:",
                                                paste(names(df)[non.0.1.binary], collapse = " "),
                                                "\nBe careful interpreting the values of interactions involving these."), call. = FALSE)
        new.int <- matrix(t(apply(d, 1, combn, 2, prod)), nrow = nrow(d))
        new[,(ncol(new) - ncol(new.int) + 1):ncol(new)] <- new.int
        new.names[(ncol(new) - ncol(new.int) + 1):ncol(new)] <- combn(names(d), 2, paste, collapse=strrep("_", nunder))
    }
    new <- data.frame(new)
    names(new) <- new.names
    single.value <- apply(new, 2, function(x) abs(max(x) - min(x)) < .Machine$double.eps ^ 0.5)
    new <- new[, !single.value, drop = FALSE]
    return(new)
}
