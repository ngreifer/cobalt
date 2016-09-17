#Think about sample sizes for clustered data. Think about subclass summaries with clustered data.

bal.tab <- function(x, ...) UseMethod("bal.tab")

base.bal.tab <- function(object, weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, ...) {
    #object: A matchit or ps object or a data.frame of treatment and weights; for computing sample sizes. obj can have clusters within it.
    #weights: A vector of weights
    #treat: A vector of treatment status
    #distance: A vector of distance measures (i.e. propensity scores or a function thereof)
    #subclass: A factor vector of subclass membership. If NULL, no subclasses are considered. 
    #covs: A data.frame of covariates for balance checking; thosed used in the ps analysis
    #call: the call of the original conditioning function
    #method: the method the for conditioning: weighting, matching, or subclassifcation.
    #cluster: a factor vector of cluster membership; if NULL, no clusters are considered.
    #cluster.summary: whether to display across cluster summary
    #which.cluster: a vector of cluster names or numbers of which clusters to print. All clusters are calculated.
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
        if (isTRUE(all.equal(m.dif, 0))) s.diff <- 0
        else s.diff <- m.dif/s.d
        if (!is.finite(s.diff)) s.diff <- NA
        return(s.diff)
    }
    std.diff.subclass <- function(x, group, weights, subclass, which.sub, denom) {
        g0 <- x[group==0 & !is.na(subclass) & subclass==which.sub & weights==1];     g1 <- x[group==1 & !is.na(subclass) & subclass==which.sub & weights==1];
        m0 <- mean(g0);        m1 <- mean(g1)
        v0 <- var(x[group==0]) ; v1 <- var(x[group==1])
        s.d <- switch(denom, 
                      control = sqrt(v0), 
                      treated = sqrt(v1), 
                      pooled = sqrt((v0+v1)/2))
        m.dif <- m1-m0
        if (isTRUE(all.equal(m.dif, 0))) s.diff <- 0
        else s.diff <- m.dif/s.d
        if (!is.finite(s.diff)) s.diff <- NA
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
                    if (length(unique(d[, i])) > 2) {
                        if (k>1) new <- data.frame(new, v = d[, i]^p)
                        else new <- data.frame(v = d[, i]^p)
                        colnames(new)[k] <- paste0(colnames(d)[i], strrep("_", nunder), p)
                    }
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
            C <- cbind(C, addl)
        } 
        C <- as.data.frame(lapply(C, function(x) if (is.character(x)) factor(x) else x))
        for (i in names(C)) {
            if (is.character(C[, i])) C[, i] <- factor(C[, i])
            if (length(unique(C[, i])) <= 2) {
                if (is.logical(C[, i])) C[, i] <- as.numeric(C[, i])
                else if (is.numeric(C[, i])) C[, i] <- binarize(C[, i])
            }
            
            if (length(cluster) > 0 && qr(matrix(c(C[, i], as.numeric(cluster)), ncol = 2))$rank == 1) C <- C[, names(C) != i] #Remove variable if it is the same (linear combo) as cluster variable
            else if (!is.numeric(C[, i])) C <- split.factor(i, C)
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
        
        #Remove duplicate & redundant variables
        suppressWarnings(C.cor <- cor(C))
        s <- matrix(sapply(C.cor, function(x) isTRUE(all.equal(abs(x), 1))), nrow = nrow(C.cor), dimnames = dimnames(C.cor))
        remove <- c()
        for (i in seq_len(ncol(s))) {
            for (j in seq_len(nrow(s))) {
                if (i > j) {
                    if (s[j, i]) remove <- c(remove, rownames(s)[i])
                }
            }
        }
        
        C <- C[,is.na(match(names(C), remove))]
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
        else if (method == "weighting") {
            t <- obj$treat
            w <- obj$weights
            obs <- as.data.frame(matrix(c(sum(t==0), 
                                          (sum(w[t==0])^2)/sum(w[t==0]^2), 
                                          sum(t==1), 
                                          (sum(w[t==1])^2)/sum(w[t==1]^2)), ncol=2, dimnames = list(c("Unadjusted", "Adjusted"), c("Control", "Treated"))))
            attr(obs, "tag") <- "Effective sample sizes:"
        }
        else if (method == "subclassification") {
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
    balance.table <- function(covs, weights, treat, distance, int, addl, continuous, binary, s.d.denom, m.threshold=NULL, v.threshold=NULL, cluster = NULL) {
        #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
        C <- get.C(covs, int, addl, distance, cluster = cluster)
        
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
            B[,"M.Threshold"] <- ifelse(B[,"Type"]=="Distance" | is.na(B[, m.threshcheck]), "", paste0(ifelse(abs(B[, m.threshcheck]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
        }
        if (!is.null(v.threshold)) {
            v.threshcheck <- ifelse(is.null(weights), "V.Ratio.Un", "V.Ratio.Adj")
            B[,"V.Threshold"] <- ifelse(B[,"Type"]=="Contin." & !is.na(B[, v.threshcheck]), paste0(ifelse(B[, v.threshcheck] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
        }
        return(B)
    }
    balance.table.subclass <- function(covs, weights, treat, subclass, distance, int, addl, continuous, binary, s.d.denom, m.threshold=NULL, v.threshold=NULL) {
        #Creates list SB of balance tables for each subclass
        
        #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
        C <- get.C(covs, int, addl, distance)

        #B=Balance frame
        Bnames <- c("Type", "M.C.Adj", "M.T.Adj", "Diff.Adj", "M.Threshold", "V.Ratio.Adj", "V.Threshold")
        B <- as.data.frame(matrix(nrow=ncol(C), ncol=length(Bnames)))
        colnames(B) <- Bnames
        rownames(B) <- names(C)
        #Set var type (binary/continuous)
        B[,"Type"] <- mapply(function(x, y) ifelse(ifelse(is.null(attr(C, "distance.name")), FALSE, x==attr(C, "distance.name")), "Distance", ifelse(length(unique(y))<=2, "Binary", "Contin.")), rownames(B), C)
        
        SB <- vector("list", length(levels(subclass)))
        names(SB) <- levels(subclass)
        
        ############################
        for (i in levels(subclass)) {
            
            SB[[i]] <- B
            
            SB[[i]][,"M.C.Adj"] <- apply(C[treat==0 & !is.na(subclass) & subclass==i & weights>0, ], 2, mean)
            SB[[i]][,"M.T.Adj"] <- apply(C[treat==1 & !is.na(subclass) & subclass==i & weights>0, ], 2, mean)
            
            #Mean differences
            SB[[i]][,"Diff.Adj"] <- mapply(function(x, type) diff.selector(x=C[, x], group=treat, weights=as.numeric(weights>0), subclass=subclass, which.sub=i, x.type=type, continuous=continuous, binary=binary, s.d.denom=s.d.denom), x=rownames(SB[[i]]), type=B[,"Type"])
            
            #Variance ratios
            if ("Contin." %in% SB[[i]][,"Type"]) {
                SB[[i]][,"V.Ratio.Adj"] <- mapply(function(x, y) var.ratio(C[!is.na(subclass) & subclass==i & weights>0, x], treat[!is.na(subclass) & subclass==i & weights>0], weights=NULL, y), rownames(SB[[i]]), SB[[i]][,"Type"])
            }
            
            if (!is.null(m.threshold)) {
                SB[[i]][,"M.Threshold"] <- ifelse(SB[[i]][,"Type"]=="Distance", "", paste0(ifelse(abs(SB[[i]][, "Diff.Adj"]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
            }
            if (!is.null(v.threshold)) {
                SB[[i]][,"V.Threshold"] <- ifelse(SB[[i]][,"Type"]=="Contin.", paste0(ifelse(SB[[i]][, "V.Ratio.Adj"] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
            }
        }
        
        return(SB)
    }
    balance.table.across.subclass <- function(balance.table, balance.table.subclass.list, subclass.obs, sub.by=NULL, m.threshold=NULL, v.threshold=NULL) {
        #Variance ratio and v.threshold not yet supported
        q.table <- array(0, dim=c(nrow(balance.table.subclass.list[[1]]), 3, length(balance.table.subclass.list)), dimnames=list(row.names(balance.table.subclass.list[[1]]), names(balance.table.subclass.list[[1]][, 2:4])))
        for (i in seq_along(balance.table.subclass.list)) {
            q.table[, , i] <- as.matrix(balance.table.subclass.list[[i]][, 2:4])
        }
        
        if (is.null(sub.by)){
            sub.by <- "treat"
        }
        if (sub.by=="treat") {
            wsub <- subclass.obs["Treated", ]/sum(subclass.obs["Treated", ])
        } else if (sub.by=="control") {
            wsub <- subclass.obs["Control", ]/sum(subclass.obs["Control", ])
        } else if (sub.by=="all") {
            wsub <- subclass.obs["Total", ]/sum(subclass.obs["Total", ])
        }
        B.A <- array(0, dim=dim(q.table[, , 1]), dimnames = dimnames(q.table[, , 1]))
        for(i in rownames(B.A)) {
            for(j in colnames(B.A)) {
                if (j=="Diff.Adj") {
                    B.A[i, j] <- sqrt(sum((wsub^2)*(q.table[i, j, ]^2)))
                }
                # else if (j=="V.Ratio.Adj") {
                #   B.A[i, j] <- NA
                # }
                else {
                    B.A[i, j] <- sum(wsub*q.table[i, j, ])
                }
            }
        }
        B.A.df <- data.frame(balance.table[, c("Type", "M.C.Un", "M.T.Un", "Diff.Un")], 
                             B.A, M.Threshold=NA)
        if (!is.null(m.threshold)) B.A.df[,"M.Threshold"] <- ifelse(B.A.df[,"Type"]=="Distance", "", paste0(ifelse(abs(B.A.df["Diff.Adj"]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
        return(B.A.df)
    }
    balance.table.cluster.summary <- function(balance.table.clusters.list, m.threshold = NULL, v.threshold = NULL, cluster.summary = FALSE) {
        Brownames <- unique(do.call("c", lapply(balance.table.clusters.list, rownames)))
        Bcolnames <- c("Type", "Mean.Diff.Un", "Median.Diff.Un", "Max.Diff.Un", "Mean.V.Ratio.Un", "Median.V.Ratio.Un", "Max.V.Ratio.Un",
                       "Mean.Diff.Adj", "Median.Diff.Adj", "Max.Diff.Adj", "Mean.V.Ratio.Adj", "Median.V.Ratio.Adj", "Max.V.Ratio.Adj")
        B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
        names(B) <- Bcolnames
        
        B[,"Type"] <- sapply(rownames(B), function(x) {u <- unique(sapply(balance.table.clusters.list, function(y) y[x, "Type"])); return(u[!is.na(u)])})
        
        for (Fun in c("Mean", "Median", "Max")) {
            for (sample in c("Un", "Adj")) {
                B[, paste(Fun, "Diff", sample, sep = ".")] <- sapply(rownames(B), function(x) get(tolower(Fun))(sapply(balance.table.clusters.list, function(y) abs(y[x, paste0("Diff.", sample)])), na.rm = TRUE))
                B[, paste(Fun, "V.Ratio", sample, sep = ".")] <- sapply(rownames(B), function(x) if (B[x, "Type"]!="Contin.") NA else get(tolower(Fun))(sapply(balance.table.clusters.list, function(y) y[x, paste0("V.Ratio.", sample)]), na.rm = TRUE))
            }
        }
        return(B)
    }
    
    #Preparations
    if (length(unique(treat)) != 2) {
        stop("Treatment indicator must be a binary (0, 1) variable---i.e., treatment (1) or control (0)", call. = FALSE)
    }
    else {
        treat <- binarize(treat)
    }
    #if (sum(treat !=1 & treat !=0) > 0) 
    if (!is.null(m.threshold)) m.threshold <- abs(m.threshold)
    if (!is.null(v.threshold)) v.threshold <- max(v.threshold, 1/v.threshold)
    if (!is.null(v.threshold)) disp.v.ratio <- TRUE
    if (is.null(weights)) un <- TRUE
    
    #Actions
    if (length(cluster) > 0) {
        out.names <- c("Cluster.Balance", 
                       "Cluster.Balance.Across.Subclass", 
                       "Cluster.Summary", 
                       "call", "print.options")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        out[["Cluster.Balance"]] <- vector("list", length(levels(cluster)))
        names(out[["Cluster.Balance"]]) <- levels(cluster)
        
        if (method == "subclassification") {
            stop("Subclassification with clusters is not supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab")
        }
        else {
            out[["Cluster.Balance"]] <- lapply(levels(cluster), function(i) balance.table(covs = covs[cluster == i,], weights = weights[cluster == i], treat = treat[cluster == i], distance = distance[cluster == i], int = int, addl = addl[cluster == i,], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, cluster = cluster[cluster == i]))
            names(out[["Cluster.Balance"]]) <- levels(cluster)
            
            out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.table.clusters.list = out[["Cluster.Balance"]], 
                                                                 m.threshold = m.threshold, v.threshold = v.threshold, 
                                                                 cluster.summary = cluster.summary)
            if (all(sapply(out[["Cluster.Balance"]], function(x) all(is.na(x[,"V.Ratio.Un"]))))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
            out <- out[!names(out) %in% "Cluster.Balance.Across.Subclass"]
            class(out) <- c("bal.tab.cluster", "bal.tab")
        }
        
        out[["call"]] <- call
        out[["print.options"]] <- list(m.threshold=m.threshold, 
                                  v.threshold=v.threshold, 
                                  un=un, 
                                  disp.means=disp.means, 
                                  disp.v.ratio=disp.v.ratio, 
                                  disp.adj=!is.null(weights), 
                                  disp.subclass=disp.subclass,
                                  which.cluster=which.cluster,
                                  cluster.summary=cluster.summary)
    }
    else {
        if (method %in% "subclassification") {
            if (length(subclass)>0) {
                out.names <- c("Subclass.Balance", "Balance.Across.Subclass", 
                               "Balanced.Means.Subclass", "Max.Imbalance.Means.Subclass", 
                               "Balanced.Variances.Subclass", "Max.Imbalance.Variances.Subclass", 
                               "Subclass.Observations", "call", "print.options")
                out <- vector("list", length(out.names))
                names(out) <- out.names
                
                if (length(list(...)$sub.by > 0)) sub.by <- list(...)$sub.by
                else sub.by <- call$sub.by
                
                out[["Subclass.Balance"]] <- balance.table.subclass(covs=covs, weights=weights, treat=treat, subclass=subclass, distance=distance, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold)
                out[["Subclass.Observations"]] <- samplesize(object, method="subclassification")
                out[["Balance.Across.Subclass"]] <- balance.table.across.subclass(balance.table=balance.table(covs, weights, treat, distance, int, addl, continuous, binary, s.d.denom, m.threshold, v.threshold), 
                                                                             balance.table.subclass.list=out[["Subclass.Balance"]], 
                                                                             subclass.obs=out[["Subclass.Observations"]], 
                                                                             sub.by=sub.by, 
                                                                             m.threshold=m.threshold, 
                                                                             v.threshold=v.threshold)
                if (!is.null(m.threshold)) {
                    out[["Balanced.Means.Subclass"]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][,"M.Threshold"])))
                    names(out[["Balanced.Means.Subclass"]]) <- paste("Subclass", levels(subclass))
                    mims.list <- lapply(levels(subclass), function(x) {
                        mi <- max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][,"Type"]!="Distance", ], "Diff.Adj", "M.Threshold")
                        return(data.frame(Variable = row.names(mi), mi))
                    } )
                    mims <- do.call("rbind", mims.list)
                    out[["Max.Imbalance.Means.Subclass"]] <- data.frame(mims, row.names = paste("Subclass", levels(subclass)))
                }
                if (!is.null(v.threshold)) {
                    out[["Balanced.Variances.Subclass"]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][,"V.Threshold"])))
                    names(out[["Balanced.Variances.Subclass"]]) <- paste("Subclass", levels(subclass))
                    mivs.list <- lapply(c(1:max(unique(subclass), na.rm=TRUE)), function(x) {
                        mi <- max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][,"Type"]!="Distance", ], "V.Ratio.Adj", "V.Threshold")
                        return(data.frame(Variable = row.names(mi), mi))
                    } )      
                    mivs <- do.call("rbind", mivs.list)
                    out[["Max.Imbalance.Variances.Subclass"]] <- data.frame(mivs, row.names = paste("Subclass", levels(subclass)))
                }
                out[["call"]] <- call
                out[["print.options"]] <- list(m.threshold=m.threshold, 
                                          v.threshold=v.threshold, 
                                          un=un, 
                                          disp.means=disp.means, 
                                          disp.v.ratio=disp.v.ratio, 
                                          disp.adj=!is.null(weights), 
                                          disp.subclass=disp.subclass)
                class(out) <- c("bal.tab.subclass", "bal.tab")
            }
            else stop("Method specified as subclassification, but no subclasses were specified.", call. = FALSE)
        }
        else {
            out.names <- c("Balance", "Balanced.Means", 
                           "Max.Imbalance.Means", "Balanced.Variances", 
                           "Max.Imbalance.Variances", "Observations", 
                           "call", "print.options")
            out <- vector("list", length(out.names))
            names(out) <- out.names
            out[["Balance"]] <- balance.table(covs, weights, treat, distance, int, addl, continuous, binary, s.d.denom, m.threshold, v.threshold)
            if (!is.null(m.threshold)) {
                m.threshcheck <- ifelse(is.null(weights), "Diff.Un", "Diff.Adj")
                out[["Balanced.Means"]] <- baltal(out[["Balance$M.Threshold"]])
                out[["Max.Imbalance.Means"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], m.threshcheck, "M.Threshold")
            }
            if (all(is.na(out[["Balance"]][,"V.Ratio.Un"]))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
            if (!is.null(v.threshold)) {
                v.threshcheck <- ifelse(is.null(weights), "V.Ratio.Un", "V.Ratio.Adj")
                out[["Balanced.Variances"]] <- baltal(out[["Balance"]][,"V.Threshold"])
                out[["Max.Imbalance.Variances"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], v.threshcheck, "V.Threshold")
            }
            out[["Observations"]] <- samplesize(object, method = method)
            out[["call"]] <- call
            out[["print.options"]] <- list(m.threshold=m.threshold, 
                                      v.threshold=v.threshold, 
                                      un=un, 
                                      disp.means=disp.means, 
                                      disp.v.ratio=disp.v.ratio, 
                                      disp.adj=!is.null(weights))
            class(out) <- "bal.tab"
        }
    }
    
    #attr(out, "int") <- int
    return(out)
}
base.bal.tab.cont <- function(object, weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, r.threshold = NULL, un = FALSE, method = "weighting", cluster = NULL, which.cluster = NULL, cluster.summary = TRUE) {
    #object: A CBPSContinuous object
    #weights: A vector of weights
    #treat: A vector of treatment values
    #distance: A vector of distance measures (i.e. propensity scores or a function thereof)
    #subclass: NOT SUPPORTED. A factor vector of subclass membership. If NULL, no subclasses are considered. 
    #covs: A data.frame of covariates for balance checking; thosed used in the ps analysis
    #call: the call of the original conditioning function
    #r.threshold: threshold for correlation for balance
    #method: ONLY WEIGHTING FOR NOW. the method the for conditioning: weighting, matching, or subclassifcation.
    #cluster: a factor vector of cluster membership; if NULL, no clusters are considered.
    #cluster.summary: whether to display across cluster summary
    #which.cluster: a vector of cluster names or numbers of which clusters to print. All clusters are calculated.
    #*print options: un, disp.subclass (whether to show individual subclass balance)
    
    w.r <- function(x, y, w = NULL) {
        w.m <- function(x, w = NULL) {return(sum(x*w, na.rm=TRUE)/sum(w, na.rm=TRUE))}
        w.cov <- function(x, y , w = NULL) {
            wmx <- w.m(x, w)
            wmy <- w.m(y, w)
            wcov <- sum(w*(x - wmx)*(y - wmy), na.rm = TRUE)/sum(w, na.rm = TRUE)
            return(wcov)}
        
        if (length(w) == 0) w <- rep(1, length(x))
        else if (length(w) != length(x)) stop("weights must be same length as x and y")
        if (length(x) != length(y)) stop("x and y must the same length")
        
        r <- w.cov(x, y, w) / (sqrt(w.cov(x, x, w) * w.cov(y, y, w)))
        return(r)
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
        
        x <- factor(data[, varname])
        if (length(unique(x)) > 1) k <- model.matrix(~ as.character(as.matrix(x)) - 1)
        else {
            k <- matrix(0, ncol = nlevels(x), nrow = length(x))
            k[, levels(x) %in% unique(x)] <- 1
        }
        if (is.null(levels(x))) colnames(k) <- paste0(varname, 1:ncol(k))
        else colnames(k) <- paste0(varname, "_", substr(levels(x), 1, 10))
        # if (ncol(k)==2) retain <- 2
        # else retain <- NULL
        # if (!is.null(retain)) k <- k[, retain, drop = FALSE]
        
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
            C <- cbind(C, addl)
        } 
        C <- as.data.frame(lapply(C, function(x) if (is.character(x)) factor(x) else x))
        for (i in names(C)) {
            if (is.character(C[, i])) C[, i] <- factor(C[, i])
            if (length(unique(C[, i])) <= 2) {
                if (is.logical(C[, i])) C[, i] <- as.numeric(C[, i])
                else if (is.numeric(C[, i])) C[, i] <- binarize(C[, i])
            }
            
            if (length(cluster) > 0 && qr(matrix(c(C[, i], as.numeric(cluster)), ncol = 2))$rank == 1) C <- C[, names(C) != i] #Remove variable if it is the same (linear combo) as cluster variable
            else if (!is.numeric(C[, i])) C <- split.factor(i, C)
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
        
        #Remove duplicate & redundant variables
        suppressWarnings(C.cor <- cor(C))
        s <- matrix(sapply(C.cor, function(x) isTRUE(all.equal(abs(x), 1))), nrow = nrow(C.cor), dimnames = dimnames(C.cor))
        remove <- c()
        for (i in seq_len(ncol(s))) {
            for (j in seq_len(nrow(s))) {
                if (i > j) {
                    if (s[j, i]) remove <- c(remove, rownames(s)[i])
                }
            }
        }
        
        C <- C[,is.na(match(names(C), remove))]
        C <- C[!duplicated(as.list(as.data.frame(apply(C, 2, as.numeric))))] #Remove duplicate variables
        
        if (!is.null(distance)) {
            distance <- data.frame(.distance = distance)
            while (names(distance) %in% names(C)) {names(distance) <- paste0(names(distance), "_")}
            C <- cbind(distance, C)
            attr(C, "distance.name") <- names(distance)
        }
        return(C)
    }
    baltal <- function(threshold) {
        #threshold: vector of threshold values (i.e., "Balanced"/"Not Balanced")
        
        thresh.val <- substring(names(table(threshold))[2], 1+regexpr("[><]", names(table(threshold))[2]), nchar(names(table(threshold))[2]))
        b <- data.frame(count=c(sum(threshold==paste0("Balanced, <", thresh.val)), 
                                sum(threshold==paste0("Not Balanced, >", thresh.val))))
        rownames(b) <- c(paste0("Balanced, <", thresh.val), paste0("Not Balanced, >", thresh.val))
        return(b)
    }
    max.imbal <- function(balance.table, col.name, thresh.col.name) {
        return(balance.table[which.max(abs(balance.table[, col.name])), match(c(col.name, thresh.col.name), names(balance.table))])
    }
    balance.table.cont <- function(covs, weights, treat, distance, int, addl, r.threshold = NULL, cluster = NULL) {
        #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
        C <- get.C(covs, int, addl, distance, cluster = cluster)
        
        #B=Balance frame
        Bnames <- c("Type", 
                    "Corr.Un", 
                    "Corr.Adj", 
                    "R.Threshold")
        B <- as.data.frame(matrix(nrow=ncol(C), ncol=length(Bnames)))
        colnames(B) <- Bnames
        rownames(B) <- names(C)
        
        #Set var type (binary/continuous)
        B[,"Type"] <- mapply(function(x, y) ifelse(ifelse(is.null(attr(C, "distance.name")), FALSE, x==attr(C, "distance.name")), "Distance", ifelse(length(unique(y))<=2, "Binary", "Contin.")), rownames(B), C)
        
        #Correlations
        if (!is.null(weights)) B[,"Corr.Adj"] <- sapply(C, w.r, y = treat, w = weights)
        B[,"Corr.Un"] <- sapply(C, w.r, y = treat, w = NULL)
        
        if (!is.null(r.threshold)) {
            r.threshcheck <- ifelse(is.null(weights), "Corr.Un", "Corr.Adj")
            B[,"R.Threshold"] <- ifelse(B[,"Type"]=="Distance" | is.na(B[, r.threshcheck]), "", paste0(ifelse(abs(B[, r.threshcheck]) < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold))
        }
        
        return(B)
    }
    balance.table.cont.cluster.summary <- function(balance.table.clusters.list, r.threshold = NULL, cluster.summary = FALSE) {
        Brownames <- unique(do.call("c", lapply(balance.table.clusters.list, rownames)))
        Bcolnames <- c("Type", "Mean.Corr.Un", "Median.Corr.Un", "Max.Corr.Un",
                       "Mean.Corr.Adj", "Median.Corr.Adj", "Max.Corr.Adj")
        B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
        names(B) <- Bcolnames
        
        B[, "Type"] <- sapply(rownames(B), function(x) {u <- unique(sapply(balance.table.clusters.list, function(y) y[x, "Type"])); return(u[!is.na(u)])})
        
        for (Fun in c("Mean", "Median", "Max")) {
            for (sample in c("Un", "Adj")) {
                B[, paste(Fun, "Corr", sample, sep = ".")] <- sapply(rownames(B), function(x) get(tolower(Fun))(sapply(balance.table.clusters.list, function(y) abs(y[x, paste0("Corr.", sample)])), na.rm = TRUE))
            }
        }
        return(B)
    }
    
    #Preparations
    if (!is.null(r.threshold)) r.threshold <- abs(r.threshold)
    if (is.null(weights)) un <- TRUE
    
    #Actions
    if (length(cluster) > 0) {
        out.names <- c("Cluster.Balance", 
                       "Cluster.Balance.Across.Subclass", 
                       "Cluster.Summary", 
                       "call", "print.options")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        out[["Cluster.Balance"]] <- vector("list", length(levels(cluster)))
        names(out[["Cluster.Balance"]]) <- levels(cluster)
        
        
        out[["Cluster.Balance"]] <- lapply(levels(cluster), function(i) balance.table.cont(covs = covs[cluster == i,], weights = weights[cluster == i], treat = treat[cluster == i], distance = distance[cluster == i], int = int, addl = addl[cluster == i,], r.threshold = r.threshold, cluster = cluster[cluster == i]))
        names(out[["Cluster.Balance"]]) <- levels(cluster)
        
        out[["Cluster.Summary"]] <- balance.table.cont.cluster.summary(balance.table.clusters.list = out[["Cluster.Balance"]], 
                                                             r.threshold = r.threshold, 
                                                             cluster.summary = cluster.summary)
        out <- out[!names(out) %in% "Cluster.Balance.Across.Subclass"]
        class(out) <- c("bal.tab.cont.cluster", "bal.cont.tab", "bal.tab.cluster", "bal.tab")
        
        out[["call"]] <- call
        out[["print.options"]] <- list(r.threshold=r.threshold, 
                                  un=un, 
                                  disp.adj=!is.null(weights), 
                                  which.cluster=which.cluster,
                                  cluster.summary=cluster.summary)
        class(out) <- c("bal.tab.cont.cluster", "bal.tab.cluster", "bal.tab.cont", "bal.tab")
    }
    else {
        out.names <- c("Balance", "Balanced.Corr", 
                       "Max.Imbalance.Corr", "Observations", 
                       "call", "print.options")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        out[["Balance"]] <- balance.table.cont(covs, weights, treat, distance, int, addl, r.threshold, cluster)
        if (!is.null(r.threshold)) {
            r.threshcheck <- ifelse(is.null(weights), "Corr.Un", "Corr.Adj")
            out[["Balanced.Corr"]] <- baltal(out[["Balance$R.Threshold"]])
            out[["Max.Imbalance.Corr"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], r.threshcheck, "R.Threshold")
        }
        if (all(is.na(out[["Balance"]][,"Corr.Un"]))) {r.threshold <- NULL}
        #out$Observations <- samplesize(object, method=method)
        out[["call"]] <- call
        out[["print.options"]] <- list(r.threshold=r.threshold, 
                                  un=un, 
                                  disp.adj=!is.null(weights))
        class(out) <- c("bal.tab.cont", "bal.tab")
    }
    
    attr(out, "int") <- int
    return(out)
}

bal.tab.matchit <- function(x, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    if (any(sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), "")))) {
        for (arg.name in names(args)[sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), ""))]) {
            if (identical(as.character(args[arg.name]), "")) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base(x, cluster = cluster)
    
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, subclass=X$subclass, covs=X$covs, call=X$call, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.subclass=disp.subclass, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary)
    return(out)
}
bal.tab.ps <- function(x, full.stop.method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, ...) {
    #ps = ps object from twang
    #full.stop.method = stopping rule/estimand from twang, e.g. "es.mean.att"; code to use first if none or incorrect is requested
    args <- as.list(environment())[-1]
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    if (any(sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), "")))) {
        for (arg.name in names(args)[sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), ""))]) {
            if (identical(as.character(args[arg.name]), "")) {
                assign(arg.name, NULL)
            }
        }
    }
    
    
    #Initializing variables
    X <- x2base(x, 
                full.stop.method = full.stop.method, 
                s.d.denom = s.d.denom, 
                cluster = cluster)
    
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method="weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary)
    return(out)
}
bal.tab.Match <- function(x, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    if (any(sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), "")))) {
        for (arg.name in names(args)[sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), ""))]) {
            if (identical(as.character(args[arg.name]), "")) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base(x, 
                formula = formula,
                data = data, 
                treat = treat,
                covs = covs,
                s.d.denom = s.d.denom,
                cluster = cluster)
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary)
    return(out)
}
bal.tab.formula <- function(formula, data, weights = NULL, distance = NULL, subclass = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "pooled", "control"), m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    if (is.null(subclass)) disp.subclass <- FALSE
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    if (any(sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), "")))) {
        for (arg.name in names(args)[sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), ""))]) {
            if (identical(as.character(args[arg.name]), "")) {
                assign(arg.name, NULL)
            }
        }
    }
    
    X <- x2base(formula, 
                data = data,
                weights = weights,
                distance = distance,
                subclass = subclass,
                addl = addl,
                method = method,
                cluster = cluster)
    
    if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) {
        out <- base.bal.tab.cont(object=X$obj, weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = addl, r.threshold = r.threshold, un = un, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary)
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, subclass=X$subclass, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.subclass=disp.subclass, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary)
    return(out)
}
bal.tab.data.frame <- function(x, treat, data = NULL, weights = NULL, distance = NULL, subclass = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "pooled", "control"), m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    if (is.null(subclass)) disp.subclass <- FALSE
    
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    if (any(sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), "")))) {
        for (arg.name in names(args)[sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), ""))]) {
            if (identical(as.character(args[arg.name]), "")) {
                assign(arg.name, NULL)
            }
        }
    }
    
    X <- x2base(x,
                treat = treat,
                data = data,
                weights = weights,
                distance = distance,
                subclass = subclass,
                addl = addl,
                method = method,
                cluster = cluster)
    
    if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) {
        out <- base.bal.tab.cont(object=X$obj, weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = addl, r.threshold = r.threshold, un = un, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary)
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, subclass=X$subclass, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.subclass=disp.subclass, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary)
    return(out)
}
bal.tab.CBPS <- function(x, estimand, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, ...) {
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    if (any(sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), "")))) {
        for (arg.name in names(args)[sapply(formals()[-c(1, length(formals()))], function(x) identical(as.character(x), ""))]) {
            if (identical(as.character(args[arg.name]), "")) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base(x, 
                estimand = estimand, 
                s.d.denom = s.d.denom,
                cluster = cluster)
    if (any(class(x) == "CBPSContinuous")) {
        out <- base.bal.tab.cont(object=X$obj, weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = addl, r.threshold = r.threshold, un = un, method = "weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary)
    }
    else if (nlevels(as.factor(X$treat)) > 2) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method="weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary)
    return(out)
}