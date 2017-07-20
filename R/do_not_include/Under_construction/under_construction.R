#base.bal.tab.mi Generate balance tables for multiply imputed datasets. 
#Choices to display balance for each one individually (somewhat similar to subclasses) or overall summary
#Overall summary might be max/min/mean/median statistics for each covariate
#Should be able to take arguments either as single objects to apply to all or multiple objects in lists, 
#or demarcated with imputation indicators (e.g., output from MICE). 

skew.diff <- function(x, group, weights, var.type) {
    #Currently unused
    #Calculates weighted skew and skew difference for groups; uses formulas at http://www.nematrian.com/WeightedMomentsAndCumulants
    #Weighted variance is different from one used for std.diff.
    skew.cols <- data.frame(Skew.C = NA, Skew.T = NA, Skew.Diff = NA)
    if (is.null(weights)) weights <- rep(1, length(x))
    if (var.type=="Contin.") {
        # skew for each group
        x0 <- x[group==0];       x1 <- x[group==1];
        w0 <- weights[group==0]; w1 <- weights[group==1]
        m0 <- w.m(x0, w0);       m1 <- w.m(x1, w1)
        wvp <- function(x, w, m) return(sum(w*(x-m)^2, na.rm=TRUE)/sum(w, na.rm=TRUE)) #weighted variance population
        Fv <- function(w) return((sum(w0, na.rm=TRUE)^2-sum(w0^2, na.rm=TRUE))/(sum(w0, na.rm=TRUE)^2)) #adjustment for sample variance
        wskp <- function(x, w, m, sp) return(sum(w*((x-m)/sp)^3, na.rm=TRUE)/sum(w, na.rm=TRUE)) #weighted skew population
        Fsk <- function(w, sp, s) return(((sum(w0, na.rm=TRUE)^3 - 3*sum(w0^2, na.rm=TRUE)*sum(w0, na.rm=TRUE) + 2*sum(w0^3, na.rm=TRUE))/(sum(w0, na.rm=TRUE)^3))*(sp/s)^3 ) #adjustment for sample skew
        #group==0
        sp0 <- sqrt(wvp(x0, w0, m0))
        s0 <- sqrt(wvp(x0, w0, m0)/Fv(w0))
        skp0 <- wskp(x0, w0, m0, sp0)
        sk0 <- skp0/Fsk(w0, sp0, s0)
        #group==1
        sp1 <- sqrt(wvp(x1, w1, m1))
        s1 <- sqrt(wvp(x1, w1, m1)/Fv(w1))
        skp1 <- wskp(x1, w1, m1, sp1)
        sk1 <- skp1/Fsk(w1, sp1, s1)
        
        skew.cols[, ] <- c(sk0, sk1, sk1-sk0)
    }
    return(skew.cols)
}

dum2factor <- function(data, varlist, newvar, other = NULL, sep = NULL) {
    #Like splitfactor; creates factor variable from series of dummies. Could be used as a utility.
    dummies <- data[, varlist]
    if (sum(rowSums(dummies)) < nrow(dummies)) dummies <- setNames(data.frame(dummies, ifelse(rowSums(dummies)==0, 1, 0)),c(names(dummies), other))
    factor <- factor(as.matrix(dummies) %*% 1:ncol(dummies), labels = colnames(dummies))
    
    if (!is.null(sep)) {
        if (sep=="") sep <- newvar
        levels(factor) <- unlist(strsplit(colnames(dummies), sep))[2 * (1:ncol(dummies))]
    }
    
    if (match(varlist[1], names(data)) == 1){
        data <- data.frame(factor, data[, !names(data) %in% varlist, drop = FALSE])
        names(data)[1] <- newvar
    }
    else if (match(varlist[length(varlist)], names(data)) == ncol(data)) {
        data <- data.frame(data[, !names(data) %in% varlist, drop = FALSE], factor)
        names(data)[ncol(data)] <- newvar
    }
    else {
        where.a <- match(varlist[1], names(data))
        where.b <- match(varlist[length(varlist)], names(data))
        data <- data.frame(data[, 1:(where.a-1), drop = FALSE], factor, data[, (where.b+1):ncol(data), drop = FALSE])
        names(data)[where.a] <- newvar
    }
    return(data)
}