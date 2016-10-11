#base.bal.tab.mi Generate balance tables for multiply imputed datasets. 
#Choices to display balance for each one individually (somewhat similar to subclasses) or overall summary
#Overall summary might be max/min/mean/median statistics for each covariate
#Should be able to take arguments either as single objects to apply to all or multiple objects in lists, 
#or demarcated with imputation indicators (e.g., output from MICE). 
base.bal.tab.mi <- function(object, weights, treat, imputation, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, un, disp.means, disp.v.ratio, disp.subclass = FALSE, method, disp.imputations = FALSE, combine = c("mean", "median", "max"), ...) {
    #object: A matchit or ps object or a data.frame of treatment and weights; for computing sample sizes
    #weights: A vector of weights for each unit in each imputation
    #treat: A vector of treatment status for each unit in each imputation
    #imputation: A vector of imputation numbers. 
    #distance: A vector of distance measures for each unit in each imputation
    #subclass: Not supported. A factor vector of subclass membership. If NULL, no subclasses are considered. 
    #covs: A data.frame of covariates for each unit in each imputation
    #call: the call of the original conditioning function
    #addl: A data.frame of additional covariates. Just 1, no imputations
    #method: the method the for conditioning: weighting, matching, or subclassifcation.
    #combine: whether aggregate table should show mean, median, or max imbalance for variable
    #*print options: un, disp.means, disp.v.ratio, disp.subclass (whether to show individual subclass balance)
    
    #Functions
    w.m <- function(x, w=NULL) {if (is.null(w)) w <- rep(1, length(x)); return(sum(x*w, na.rm=TRUE)/sum(w, na.rm=TRUE))}
    w.v <- function(x, w=NULL) {
        if (is.null(w)) w <- rep(1, length(x))
        return(sum(w*(x-w.m(x, w))^2, na.rm=TRUE)/(sum(w, na.rm=TRUE)-1))
    }
    std.diff <- function(var, group, weights, denom) {
        #Uses variance of original group, as in MatchIt.
        g0 <- var[group==0];     g1 <- var[group==1];
        w0 <- weights[group==0]; w1 <- weights[group==1]
        m0 <- w.m(g0, w0);        m1 <- w.m(g1, w1)
        v0 <- var(g0);           v1 <- var(g1)
        s.d <- switch(denom, 
                      control = sqrt(v0), 
                      treated = sqrt(v1), 
                      pooled = sqrt((v0+v1)/2))
        m.dif <- m1-m0
        if (identical(all.equal(m.dif, 0), TRUE)) s.diff <- 0
        else s.diff <- m.dif/s.d
        return(s.diff)
    }
    int.poly.f <- function(df, ex=NULL, int=FALSE, poly=1, nunder=1) {
        #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
        #Only to be used in base.bal.tab; for general use see int.poly()
        #df=data frame input
        #ex=data frame of variables to exclude in interactions and polynomials; a subset of df
        #int=whether to include interactions or not; currently only 2-way are supported
        #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
        #nunder=number of underscores between variables
        
        if (!is.null(ex)) d <- df[, !names(df) %in% names(ex)]
        else d <- df
        nd <- ncol(d)
        k <- 1
        if (poly>1) {
            for (i in 1:ncol(d)) {
                for (p in 2:poly) {
                    if (k>1) new <- data.frame(new, v = d[, i]^p)
                    else new <- data.frame(v = d[, i]^p)
                    colnames(new)[k] <- paste0(colnames(d)[i], strrep("_", nunder), p)
                    k <- k + 1
                }
            }
        }
        if (int==TRUE) {
            for (i in 1:(ncol(d)-1)) {
                for (j in (i+1):ncol(d)) {
                    if (k>1) new <- data.frame(new, v = d[, i]*d[, j])
                    else new <- data.frame(v = d[, i]*d[, j])
                    colnames(new)[k] <- paste0(colnames(d)[i], strrep("_", nunder), colnames(d)[j])
                    k <- k + 1
                }
                
            }
        }
        if (!is.null(ncol(new)) > 0) return(new)
    }
    split.factor <- function(varname, data) {
        #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
        #Retains all categories unless only 2 levels, in which case only the second level is retained.
        #varname= the name of the variable to split when data is specified
        #data=data set to be changed
        
        x <- as.factor(data[, varname])
        k <- model.matrix(~ as.character(as.matrix(x)) - 1)
        if (is.null(levels(x))) colnames(k) <- paste0(varname, 1:ncol(k))
        else colnames(k) <- paste0(varname, "_", substr(levels(x), 1, 10))
        if (ncol(k)==2) retain <- 2
        else retain <- NULL
        if (!is.null(retain)) k <- k[, retain]
        if (match(varname, names(data)) == ncol(data)){
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
        if (0 %in% unique(variable)) zero <- 0
        else zero <- min(unique(variable))
        variable <- ifelse(variable==zero, 0, 1)
        return(variable)
    }
    get.C <- function(covs, int, addl, distance) {
        #gets C data.frame, which contains all variables for which balance is to be assessed. Used in balance.table.
        C <- covs
        if (!is.null(addl)) {
            if (!is.data.frame(addl)) {
                if (is.character(addl)) stop("The argument to addl must be a data.frame containing the the values of the additional variables you want to include in the balance assessment.")
                else stop("The argument to addl must be a data.frame. Wrap data.frame() around the argument if it is a matrix or vector.")
            }
            C <- cbind(C, addl)
        } 
        for (i in names(C)) {
            if (length(unique(C[, i])) <= 2) {
                if (is.logical(C[, i])) C[, i] <- as.numeric(C[, i])
                else if (is.numeric(C[, i])) C[, i] <- binarize(C[, i])
            }
            if (!is.numeric(C[, i])) C <- split.factor(i, C)
        }
        
        if (int) {
            #Prevent duplicate var names with _'s
            nunder <- 1
            repeat {
                if (all(sapply(names(C), function(x) !x %in% paste0(names(C), strrep("_", nunder), "2"))) & all(sapply(names(C), function(x) !x %in% do.call(paste, c(expand.grid(names(C), names(C)), list(sep = strrep("_", nunder))))))) break
                else nunder <- nunder + 1
            }
            C <- cbind(C, int.poly.f(C, int = TRUE, poly = 2, nunder = nunder))
        }
        C <- C[!duplicated(as.list(as.data.frame(apply(C, 2, as.numeric))))] #Remove duplicate variables
        
        if (!is.null(distance)) {
            distance <- data.frame(.distance = distance)
            while (names(distance) %in% names(C)) {names(distance) <- paste0(names(distance), "_")}
            C <- cbind(distance, C)
            attr(C, "distance.name") <- names(distance)
        }
        return(C)
    }
    diff.selector <- function(x, group, weights = NULL, subclass = NULL, which.sub = NULL, x.type, continuous, binary, s.d.denom) {
        if (is.null(weights)) weights <- rep(1, length(x))
        if (is.null(which.sub)) ss <- rep(TRUE, length(x))
        else ss <- (subclass == which.sub & weights > 0)
        
        if (x.type=="Distance")  {
            diff <- w.m(x[group==1 & ss], w = weights[group==1 & ss]) - w.m(x[group==0 & ss], w = weights[group==0 & ss])
        }
        else if (x.type=="Binary") {
            if      (binary=="raw") diff <- w.m(x[group==1 & ss], w = weights[group==1 & ss]) - w.m(x[group==0 & ss], w = weights[group==0 & ss])
            else if (binary=="std") {
                if (is.null(which.sub)) diff <- std.diff(x, group, weights, denom = s.d.denom)
                else diff <- std.diff.subclass(x, group, weights, subclass, which.sub, denom = s.d.denom)
            }
            else diff <- NULL
        }
        else if (x.type=="Contin.") {
            if      (continuous=="raw") diff <- w.m(x[group==1 & ss], w = weights[group==1 & ss]) - w.m(x[group==0 & ss], w = weights[group==0 & ss])
            else if (continuous=="std") {
                if (is.null(which.sub)) diff <- std.diff(x, group, weights, denom = s.d.denom)
                else diff <- std.diff.subclass(x, group, weights, subclass, which.sub, denom = s.d.denom)
            }
            else diff <- NULL
        }
        return(diff)
    }  
    var.ratio <- function(x, group, weights, var.type) {
        if (is.null(weights)) weights <- rep(1, length(x))
        if (var.type=="Contin.") {
            ratio <- w.v(x[group==1], weights[group==1]) / w.v(x[group==0], weights[group==0])
            return(max(ratio, 1 / ratio))
        }
        else return(NA)
    }
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
    baltal <- function(threshold) {
        #threshold: vector of threshold values (i.e., "Balanced"/"Not Balanced")
        
        thresh.val <- substring(names(table(threshold))[2], 1+regexpr("[><]", names(table(threshold))[2]), nchar(names(table(threshold))[2]))
        b <- data.frame(count=c(sum(threshold==paste0("Balanced, <", thresh.val)), 
                                sum(threshold==paste0("Not Balanced, >", thresh.val))))
        rownames(b) <- c(paste0("Balanced, <", thresh.val), paste0("Not Balanced, >", thresh.val))
        return(b)
    }
    samplesize <- function(obj, method=c("matching", "weighting", "subclassification")) {
        #Computes sample size info. for unadjusted and adjusted samples. obj should be a matchit, ps (twang), 
        # or matching object, or a data.frame containing the variables "weights" and "treat" (and optionally 
        # "subclass") for each unit. method is what method the weights are to be used for. 
        # method="subclassification is for subclass sample sizes only.
        method <- match.arg(method)
        if (method=="matching") {
            if (any(class(obj)=="matchit")) {
                nn <- as.data.frame(matrix(0, ncol=2, nrow=4))
                nn[1, ] <- c(sum(obj$treat==0), sum(obj$treat==1))
                nn[2, ] <- c(sum(obj$treat==0 & obj$weights>0), sum(obj$treat==1 & obj$weights>0))
                nn[3, ] <- c(sum(obj$treat==0 & obj$weights==0 & obj$discarded==0), sum(obj$treat==1 & obj$weights==0 & obj$discarded==0))
                nn[4, ] <- c(sum(obj$treat==0 & obj$weights==0 & obj$discarded==1), sum(obj$treat==1 & obj$weights==0 & obj$discarded==1))
                dimnames(nn) <- list(c("All", "Matched", "Unmatched", "Discarded"), 
                                     c("Control", "Treated"))
            }
            else {
                if (all(is.na(obj$weights))) {
                    nn <- as.data.frame(matrix(0, ncol=2, nrow=1))
                    nn[1, ] <- c(sum(obj$treat==0), sum(obj$treat==1))
                    dimnames(nn) <- list(c("All"), 
                                         c("Control", "Treated"))
                }
                else {
                    nn <- as.data.frame(matrix(0, ncol=2, nrow=3))
                    nn[1, ] <- c(sum(obj$treat==0), sum(obj$treat==1))
                    nn[2, ] <- c(sum(obj$treat==0 & obj$weights>0), sum(obj$treat==1 & obj$weights>0))
                    nn[3, ] <- c(sum(obj$treat==0 & obj$weights==0), sum(obj$treat==1 & obj$weights==0))
                    dimnames(nn) <- list(c("All", "Matched", "Unmatched"), 
                                         c("Control", "Treated"))
                }
            }
            obs <- nn
            attr(obs, "tag") <- "Sample sizes:"
        }
        else if (method=="weighting") {
            t <- obj$treat
            w <- obj$weights
            obs <- as.data.frame(matrix(c(sum(t==0), 
                                          (sum(w[t==0])^2)/sum(w[t==0]^2), 
                                          sum(t==1), 
                                          (sum(w[t==1])^2)/sum(w[t==1]^2)), ncol=2, dimnames = list(c("Unadjusted", "Adjusted"), c("Control", "Treated"))))
            attr(obs, "tag") <- "Effective sample sizes:"
        }
        else if (method=="subclassification") {
            if (is.null(obj$subclass)) stop("obj must contain a vector of subclasses called \"subclass\".")
            qbins <- length(levels(obj$subclass))
            
            nn <- as.data.frame(matrix(0, 3, qbins))
            dimnames(nn) <- list(c("Control", "Treated", "Total"), 
                                 paste("Subclass", levels(obj$subclass)))
            
            matched <- obj$weights!=0
            k <- 0
            for (i in levels(obj$subclass)) {
                qi <- obj$subclass[matched]==i & (!is.na(obj$subclass[matched]))
                qt <- obj$treat[matched][qi]
                if (sum(qt==1)<2|(sum(qt==0)<2)){
                    if (sum(qt==1)<2)
                        warning("Not enough treatment units in subclass ", i, call. = FALSE)
                    else if (sum(qt==0)<2)
                        warning("Not enough control units in subclass ", i, call. = FALSE)
                }
                k <- k + 1
                nn[, k] <- c(sum(qt==0), sum(qt==1), length(qt))
            }
            obs <- nn
            attr(obs, "tag") <- "Sample sizes by subclass:"
        }
        return(obs)
    }
    max.imbal <- function(balance.table, col.name, thresh.col.name) {
        return(balance.table[which.max(abs(balance.table[, col.name])), match(c(col.name, thresh.col.name), names(balance.table))])
    }
    balance.table <- function(covs, weights, treat, distance, int, addl, continuous, binary, s.d.denom, m.threshold=NULL, v.threshold=NULL) {
        #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
        C <- get.C(covs, int, addl, distance)
        
        #B=Balance frame
        Bnames <- c("Type", 
                    "M.C.Un", 
                    "M.T.Un", 
                    "Diff.Un", 
                    "V.Ratio.Un",  
                    "M.C.Adj", 
                    "M.T.Adj", 
                    "Diff.Adj", 
                    "M.Threshold", 
                    "V.Ratio.Adj", 
                    "V.Threshold")
        B <- as.data.frame(matrix(nrow=ncol(C), ncol=length(Bnames)))
        colnames(B) <- Bnames
        rownames(B) <- names(C)
        
        #Set var type (binary/continuous)
        B[,"Type"] <- mapply(function(x, y) ifelse(ifelse(is.null(attr(C, "distance.name")), FALSE, x==attr(C, "distance.name")), "Distance", ifelse(length(unique(y))<=2, "Binary", "Contin.")), rownames(B), C)
        
        B[,"M.C.Un"] <- apply(C[treat==0, ], 2, mean)
        B[,"M.T.Un"] <- apply(C[treat==1, ], 2, mean)
        if (!is.null(weights)) {
            B[,"M.C.Adj"] <- apply(C[treat==0, ], 2, w.m, w=weights[treat==0])
            B[,"M.T.Adj"] <- apply(C[treat==1, ], 2, w.m, w=weights[treat==1])
        }
        
        #Mean differences
        if (!is.null(weights)) B[,"Diff.Adj"] <- mapply(function(x, type) diff.selector(x=C[, x], group=treat, weights=weights, x.type=type, continuous=continuous, binary=binary, s.d.denom=s.d.denom), rownames(B), B[,"Type"])
        B[,"Diff.Un"] <- mapply(function(x, type) diff.selector(x=C[, x], group=treat, weights=NULL, x.type=type, continuous=continuous, binary=binary, s.d.denom=s.d.denom), rownames(B), B[,"Type"])
        
        #Variance ratios
        if ("Contin." %in% B[,"Type"]) {
            if (!is.null(weights)) B[,"V.Ratio.Adj"] <- mapply(function(x, y) var.ratio(C[, x], treat, weights, y), rownames(B), B[,"Type"])
            B[,"V.Ratio.Un"] <- mapply(function(x, y) var.ratio(C[, x], treat, NULL, y), rownames(B), B[,"Type"])
        }
        if (prod(is.na(B[,"V.Ratio.Un"]))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
        
        if (!is.null(m.threshold)) {
            m.threshcheck <- ifelse(is.null(weights), "Diff.Un", "Diff.Adj")
            B[,"M.Threshold"] <- ifelse(B[,"Type"]=="Distance", "", paste0(ifelse(abs(B[, m.threshcheck]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
        }
        if (!is.null(v.threshold)) {
            v.threshcheck <- ifelse(is.null(weights), "V.Ratio.Un", "V.Ratio.Adj")
            B[,"V.Threshold"] <- ifelse(B[,"Type"]=="Contin.", paste0(ifelse(B[, v.threshcheck] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
        }
        return(B)
    }
    balance.table.mi <- function(bal.list, combine, m.threshold=NULL, v.threshold=NULL) {
        m.table <- array(0, dim=c(nrow(bal.list[[1]]), ncol(bal.list[[1]]), length(bal.list)),
                         dimnames = list(row.names(bal.list[[1]]), names(bal.list[[1]])))
        for (i in seq_along(bal.list)) {
            m.table[, , i] <- as.matrix(bal.list[[i]])
        }
        B <- bal.list[[1]]
        
        for (i in nrow(B)){
            for (j in names(B)) {
                if (j %in% c("M.C.Un", "M.T.Un", "Diff.Un", "V.Ratio.Un", "M.C.Adj", "M.T.Adj", "Diff.Adj", "V.Ratio.Adj")) {
                    B[i, j] <- get(combine)(m.table[i, j, ])
                }
            }
        }
        if (!is.null(m.threshold)) B[,"M.Threshold"] <- ifelse(B[,"Type"]=="Distance", "", paste0(ifelse(abs(B["Diff.Adj"]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
        if (!is.null(v.threshold)) B[,"V.Threshold"] <- ifelse(B[,"Type"]=="Contin.", paste0(ifelse(B[, "V.Ratio.Adj"] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
        return(B)
    }
    
    #Actions
    if (sum(treat !=1 & treat !=0) > 0) {
        stop("Treatment indicator must be a binary (0, 1) variable---i.e., treatment (1) or control (0)")}
    if (!is.null(m.threshold)) m.threshold <- abs(m.threshold)
    if (!is.null(v.threshold)) v.threshold <- max(v.threshold, 1/v.threshold)
    if (!is.null(v.threshold)) disp.v.ratio <- TRUE
    if (is.null(weights)) un <- TRUE
    combine <- match.arg(combine)
    
    if (method %in% "subclassification") {
        stop("Multiple imputation not supported with subclasses")
    }
    else {
        out.names <- c("Imputation.Balance", "Balance.Across.Imputations", "Observations", 
                       "call", "print.options")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        out$Imputation.Balance <- lapply(unique(imputation), function(i) balance.table(covs[imputation==i,], weights[imputation==i], treat[imputation==i], distance[imputation==i], int, addl, continuous, binary, s.d.denom, m.threshold, v.threshold))
        out$Balance.Across.Imputations <- balance.table.mi(out$Imputation.Balance, combine, m.threshold, v.threshold)
        
        #out$Balance <- balance.table(covs, weights, treat, distance, int, addl, continuous, binary, s.d.denom, m.threshold, v.threshold)
        if (!is.null(m.threshold)) {
            #m.threshcheck <- ifelse(is.null(weights), "Diff.Un", "Diff.Adj")
            #out$Balanced.Means <- baltal(out$Balance$M.Threshold)
            #out$Max.Imbalance.Means <- max.imbal(out$Balance[out$Balance[,"Type"]!="Distance", ], m.threshcheck, "M.Threshold")
        }
        #if (all(is.na(out$Balance[,"V.Ratio.Un"]))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
        if (!is.null(v.threshold)) {
            #v.threshcheck <- ifelse(is.null(weights), "V.Ratio.Un", "V.Ratio.Adj")
            #out$Balanced.Variances <- baltal(out$Balance[,"V.Threshold"])
            #out$Max.Imbalance.Variances <- max.imbal(out$Balance[out$Balance[,"Type"]!="Distance", ], v.threshcheck, "V.Threshold")
        }
        #out$Observations <- samplesize(object, method=method)
        # out$call <- call
        # out$print.options <- list(m.threshold=m.threshold, 
        #                           v.threshold=v.threshold, 
        #                           un=un, 
        #                           disp.means=disp.means, 
        #                           disp.v.ratio=disp.v.ratio, 
        #                           disp.adj=!is.null(weights),
        #                           disp.imp=disp.imputations)
        class(out) <- c("bal.tab.mi", "bal.tab")
    }
    return(out)
}

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
