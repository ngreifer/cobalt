#x2base
match.strata2weights <- function(covs, treat, match.strata) {
    #Process match.strata into weights (similar to weight.subclass from MatchIt)
    names(treat) <- row.names(covs)
    matched <- !is.na(match.strata); unmatched <- !matched
    treat.matched <- treat[matched]
    match.strata.matched <- match.strata[matched]
    
    labels <- names(treat.matched)
    tlabels <- labels[treat.matched == 1]
    clabels <- labels[treat.matched == 0]
    
    weights.matched <- rep(0, length(treat.matched))
    names(weights.matched) <- labels
    weights.matched[tlabels] <- 1
    
    for (j in unique(match.strata.matched)){
        qn0 <- sum(treat.matched==0 & match.strata.matched==j)
        qn1 <- sum(treat.matched==1 & match.strata.matched==j)
        weights.matched[treat.matched==0 & match.strata.matched==j] <- qn1/qn0
    }
    if (sum(weights.matched[clabels]) < sqrt(.Machine$double.eps)) #all control weights are 0
        weights.matched[clabels] <- rep(0, length(weights.matched[clabels]))
    else {
        ## Number of C units that were matched to at least 1 T
        num.cs <- sum(weights.matched[clabels] > sqrt(.Machine$double.eps))
        weights.matched[clabels] <- weights.matched[clabels]*num.cs/sum(weights.matched[clabels])
    }
    
    
    if (any(unmatched)) {
        weights.unmatched <- rep(0, sum(unmatched))
        names(weights.unmatched) <- names(treat[unmatched])
        weights <- c(weights.matched, weights.unmatched)[names(treat)]
    }
    else {
        weights <- weights.matched
    }
    
    if (sum(weights) < sqrt(.Machine$double.eps)) 
        stop("No units were matched", call. = FALSE)
    else if (sum(weights[tlabels]) < sqrt(.Machine$double.eps))
        stop("No treated units were matched", call. = FALSE)
    else if (sum(weights[clabels]) < sqrt(.Machine$double.eps))
        stop("No control units were matched", call. = FALSE)
    return(weights)
}
use.tc.fd <- function(formula, data, treat, covs) {
    useWhich <- function(f, d, t, c) {
        #output: tc, fd, both
        good <- c(formula=0, data=0, covs=0, treat=0)
        if (!is.null(f) & class(f)=="formula") good[1] <- 1
        if (!is.null(d) & is.data.frame(d)) good[2] <- 1
        if (!is.null(c) & is.data.frame(c)) good[3] <- 1
        if (!is.null(t) & length(unique(t))==2) good[4] <- 1
        
        if (sum(good) <= 1) {
            stop("Either formula and data or treat and covs must be specified correctly.", call. = FALSE)
        }
        else if (sum(good) == 2) {
            if (sum(good[c("formula", "data")])==0) {
                use.which <- "tc"}
            else if (sum(good[c("formula", "data")])==1) {
                stop("Either formula and data or treat and covs must be specified correctly.", call. = FALSE)}
            else { #if (sum(good[c("formula", "data")])==2)
                use.which <- "fd"}
        }
        else if (sum(good) == 3) {
            bad <- names(good)[match(0, good)]
            if (bad %in% c("formula", "data")) {
                warning(paste("Argument to", bad, "is missing; using treat and covs instead."), call. = FALSE, immediate.=TRUE)
                use.which <- "tc"}
            else {
                warning(paste("Argument to", bad, "is missing; using formula and data instead."), call. = FALSE, immediate.=TRUE)
                use.which <- "fd"}
        }
        else { #if (sum(good)==4)
            use.which <- "both"
        }
        return(use.which)
    }
    use.fd <- function(f, d){
        #outputs a list containing treat [1] and covs [2]
        out.list <- list(treat=NA, covs=data.frame(NA))
        tt <- terms(f)
        attr(tt, "intercept") <- 0
        mf<- tryCatch(model.frame(tt, d),
                      error = function(cond) stop(paste0(c("All right hand side variables of formula must be variables in data.\nVariables not in data: ",
                                                           paste(attr(tt, "term.labels")[which(!attr(tt, "term.labels") %in% names(d))], collapse=", "))), call. = FALSE))
        out.list$treat <- model.response(mf) #treat
        out.list$covs <- as.data.frame(model.matrix(tt, data=mf)) #covs
        attr(out.list, "which") <- "fd"
        return(out.list)
    }
    use.tc <- function(t, c) {
        if (length(t)!=nrow(c)) {
            stop("treat must be the same length as covs.", call. = FALSE)}
        out.list <- list(treat=t, covs=c)
        attr(out.list, "which") <- "tc"
        return(out.list)
    }
    
    use.which <- useWhich(formula, data, treat, covs)
    
    if (use.which == "fd") {
        t.c <- use.fd(formula, data)}
    else if (use.which == "tc") {
        t.c <- use.tc(treat, covs)}
    else if (use.which == "both") {
        try.fd <- try({t.c <- use.fd(formula, data)})
        if (class(try.fd) == "try-error") {
            message("Formula, data, treat, and covs all supplied; ignoring formula and data.")
            t.c <- use.tc(treat, covs)}
        else {
            message("Formula, data, treat, and covs all supplied; ignoring treat and covs.")}
    }
    return(t.c)
}

#base.bal.tab
#w.m <- function(x, w) {return(sum(x*w, na.rm=TRUE)/sum(w, na.rm=TRUE))}
w.m <- function(x, w = NULL) {
    if (is.null(w)) w <- rep(1, length(x))
    return(sum(x*w, na.rm=TRUE)/sum(w, na.rm=TRUE))
}
w.v <- function(x, w) {
    return(sum(w*(x-w.m(x, w))^2, na.rm=TRUE)/(sum(w, na.rm=TRUE)-1))
}
std.diff <- function(x, group, weights, denom) {
    #Uses variance of original group, as in MatchIt.
    treated <- as.logical(group)
    g0 <- x[!treated];     g1 <- x[treated];
    w0 <- weights[!treated]; w1 <- weights[treated]
    m0 <- w.m(g0, w0);        m1 <- w.m(g1, w1)
    v0 <- var(g0);           v1 <- var(g1)
    
    s.d <- switch(denom, 
                  control = sqrt(v0), 
                  treated = sqrt(v1), 
                  pooled = sqrt((v0+v1)/2))
    m.dif <- m1-m0
    
    if (abs(m.dif) < sqrt(.Machine$double.eps)) s.diff <- 0
    else s.diff <- m.dif/s.d
    if (!is.finite(s.diff)) s.diff <- NA
    return(s.diff)
}
std.diff.subclass <- function(x, group, weights, subclass, which.sub, denom) {
    treated <- as.logical(group)
    g0 <- x[!treated & !is.na(subclass) & subclass==which.sub & weights==1]
    g1 <- x[treated & !is.na(subclass) & subclass==which.sub & weights==1]
    m0 <- sum(g0)/length(g0)
    m1 <- sum(g1)/length(g1)
    v0 <- var(x[!treated])
    v1 <- var(x[treated])
    s.d <- switch(denom, 
                  control = sqrt(v0), 
                  treated = sqrt(v1), 
                  pooled = sqrt((v0+v1)/2))
    m.dif <- m1 - m0
    if (abs(m.dif) < sqrt(.Machine$double.eps)) s.diff <- 0
    else s.diff <- m.dif/s.d
    if (!is.finite(s.diff)) s.diff <- NA
    return(s.diff)
}
diff.selector <- function(x, group, weights = NULL, subclass = NULL, which.sub = NULL, x.type, continuous, binary, s.d.denom, no.weights = FALSE) {
    if (no.weights) weights <- rep(1, length(x))
    no.sub <- is.null(which.sub)
    if (no.sub) ss <- rep(TRUE, length(x))
    else ss <- (subclass == which.sub & weights > 0)
    
    if (x.type=="Distance")  {
        diff <- w.m(x[ss & group==1], w = weights[ss & group==1]) - w.m(x[ss & group==0], w = weights[ss & group==0])
    }
    else if (x.type=="Binary") {
        if      (binary=="raw") diff <- w.m(x[ss & group==1], w = weights[ss & group==1]) - w.m(x[ss & group==0], w = weights[ss & group==0])
        else if (binary=="std") {
            if (no.sub) diff <- std.diff(x, group, weights, denom = s.d.denom)
            else diff <- std.diff.subclass(x, group, weights, subclass, which.sub, denom = s.d.denom)
        }
        else diff <- NULL
    }
    else if (x.type=="Contin.") {
        if      (continuous=="raw") diff <- w.m(x[ss & group==1], w = weights[ss & group==1]) - w.m(x[ss & group==0], w = weights[ss & group==0])
        else if (continuous=="std") {
            if (no.sub) diff <- std.diff(x, group, weights, denom = s.d.denom)
            else diff <- std.diff.subclass(x, group, weights, subclass, which.sub, denom = s.d.denom)
        }
        else diff <- NULL
    }
    return(diff)
}  
var.ratio <- function(x, group, weights, var.type, no.weights = FALSE) {
    if (no.weights) weights <- rep(1, length(x))
    if (var.type=="Contin.") {
        ratio <- w.v(x[group==1], weights[group==1]) / w.v(x[group==0], weights[group==0])
        return(max(ratio, 1 / ratio))
    }
    else return(NA)
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
            nn <- matrix(0, ncol=2, nrow=4)
            nn[1, ] <- c(sum(obj$treat==0), sum(obj$treat==1))
            nn[2, ] <- c(sum(obj$treat==0 & obj$weights>0), sum(obj$treat==1 & obj$weights>0))
            nn[3, ] <- c(sum(obj$treat==0 & obj$weights==0 & obj$discarded==0), sum(obj$treat==1 & obj$weights==0 & obj$discarded==0))
            nn[4, ] <- c(sum(obj$treat==0 & obj$weights==0 & obj$discarded==1), sum(obj$treat==1 & obj$weights==0 & obj$discarded==1))
            nn <- as.data.frame(nn)
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
balance.table <- function(C, weights, treat, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, no.adj = FALSE, types = NULL, quick = FALSE) {
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
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
    B <- as.data.frame(matrix(nrow = ncol(C), ncol = length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- varnames <- colnames(C)
    
    #Set var type (binary/continuous)
    if (length(types) > 0) B[,"Type"] <- types
    else B[,"Type"] <- get.types(C)
    
    if (!((!un || !disp.means) && quick)) {
        B[,"M.C.Un"] <- apply(C[treat==0, , drop = FALSE], 2, function(x) sum(x)/length(x))
        B[,"M.T.Un"] <- apply(C[treat==1, , drop = FALSE], 2, function(x) sum(x)/length(x))
    }
    if (!no.adj && !(!disp.means && quick)) {
        B[,"M.C.Adj"] <- apply(C[treat==0, , drop = FALSE], 2, w.m, w = weights[treat==0])
        B[,"M.T.Adj"] <- apply(C[treat==1, , drop = FALSE], 2, w.m, w = weights[treat==1])
    }
    
    #Mean differences
    if (!no.adj) B[,"Diff.Adj"] <- mapply(function(x, type) diff.selector(x=C[, x], group=treat, weights=weights, x.type=type, continuous=continuous, binary=binary, s.d.denom=s.d.denom), varnames, B[,"Type"])
    if (!(!un && quick)) B[,"Diff.Un"] <- mapply(function(x, type) diff.selector(x=C[, x], group=treat, weights=NULL, x.type=type, continuous=continuous, binary=binary, s.d.denom=s.d.denom, no.weights = TRUE), varnames, B[,"Type"])
    
    #Variance ratios
    if ("Contin." %in% B[,"Type"] && !(!disp.v.ratio && quick)) {
        if (!no.adj) B[,"V.Ratio.Adj"] <- mapply(function(x, y) var.ratio(C[, x], treat, weights, y), varnames, B[,"Type"])
        B[,"V.Ratio.Un"] <- mapply(function(x, y) var.ratio(C[, x], treat, NULL, y, no.weights = TRUE), varnames, B[,"Type"])
    }
    if (all(is.na(B[,"V.Ratio.Un"]))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
    
    if (!is.null(m.threshold)) {
        m.threshcheck <- ifelse(no.adj, "Diff.Un", "Diff.Adj")
        B[,"M.Threshold"] <- ifelse(B[,"Type"]=="Distance" | is.na(B[, m.threshcheck]), "", paste0(ifelse(abs(B[, m.threshcheck]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
    }
    if (!is.null(v.threshold)) {
        v.threshcheck <- ifelse(no.adj, "V.Ratio.Un", "V.Ratio.Adj")
        B[,"V.Threshold"] <- ifelse(B[,"Type"]=="Contin." & !is.na(B[, v.threshcheck]), paste0(ifelse(B[, v.threshcheck] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
    }
    return(B)
}
balance.table.subclass <- function(C, weights, treat, subclass, continuous, binary, s.d.denom, m.threshold=NULL, v.threshold=NULL, disp.means = FALSE, disp.v.ratio = FALSE, types = NULL, quick = FALSE) {
    #Creates list SB of balance tables for each subclass
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    #B=Balance frame
    Bnames <- c("Type", "M.C.Adj", "M.T.Adj", "Diff.Adj", "M.Threshold", "V.Ratio.Adj", "V.Threshold")
    B <- as.data.frame(matrix(nrow=ncol(C), ncol=length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- varnames <- colnames(C)
    #Set var type (binary/continuous)
    if (length(types) > 0) B[,"Type"] <- types
    else B[,"Type"] <- get.types(C)
    
    SB <- vector("list", length(levels(subclass)))
    names(SB) <- levels(subclass)
    
    ############################
    for (i in levels(subclass)) {
        
        SB[[i]] <- B
        
        if (!(!disp.means && quick)) {
            SB[[i]][,"M.C.Adj"] <- apply(C[treat==0 & !is.na(subclass) & subclass==i & weights>0, ], 2, function(x) sum(x)/length(x))
            SB[[i]][,"M.T.Adj"] <- apply(C[treat==1 & !is.na(subclass) & subclass==i & weights>0, ], 2, function(x) sum(x)/length(x))
        }
        
        #Mean differences
        SB[[i]][,"Diff.Adj"] <- mapply(function(x, type) diff.selector(x=C[, x], group=treat, weights=as.numeric(weights>0), subclass=subclass, which.sub=i, x.type=type, continuous=continuous, binary=binary, s.d.denom=s.d.denom), x=rownames(SB[[i]]), type=B[,"Type"])
        
        #Variance ratios
        if ("Contin." %in% SB[[i]][,"Type"] && !(!disp.v.ratio && quick)) {
            SB[[i]][,"V.Ratio.Adj"] <- mapply(function(x, y) var.ratio(C[!is.na(subclass) & subclass==i & weights>0, x], treat[!is.na(subclass) & subclass==i & weights>0], weights=NULL, y, no.weights = TRUE), rownames(SB[[i]]), SB[[i]][,"Type"])
        }
    }
    
    if (!is.null(m.threshold)) {
        for (i in levels(subclass)) {
            SB[[i]][,"M.Threshold"] <- ifelse(SB[[i]][,"Type"]=="Distance", "", paste0(ifelse(abs(SB[[i]][, "Diff.Adj"]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
        }
    }
    
    if (all(sapply(SB, function(x) all(is.na(x[,"V.Ratio.Adj"]))))) {
        attr(SB, "dont.disp.v.ratio") <- TRUE; v.threshold <- NULL
    }
    if (!is.null(v.threshold)) {
        for (i in levels(subclass)) {
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
                         B.A, M.Threshold = NA)
    if (!is.null(m.threshold)) B.A.df[,"M.Threshold"] <- ifelse(B.A.df[,"Type"]=="Distance", "", paste0(ifelse(abs(B.A.df["Diff.Adj"]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
    return(B.A.df)
}
balance.table.cluster.summary <- function(balance.table.clusters.list, m.threshold = NULL, v.threshold = NULL, cluster.summary = FALSE, no.adj = FALSE, types = NULL, quick = FALSE) {
    
    Brownames <- unique(do.call("c", lapply(balance.table.clusters.list, rownames)))
    cluster.functions <- c("Min", "Mean", "Median", "Max")
    Bcolnames <- c("Type", apply(expand.grid(cluster.functions, c("Diff", "V.Ratio"), c("Un", "Adj")), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
    names(B) <- Bcolnames
    
    if (length(types) > 0) B[,"Type"] <- types
    else B[,"Type"] <- unlist(sapply(Brownames, function(x) {u <- unique(sapply(balance.table.clusters.list, function(y) y[x, "Type"])); return(u[!is.na(u)])}), use.names = FALSE)
    
    funs <- structure(vector("list", length(cluster.functions)), names = cluster.functions)
    for (Fun in cluster.functions[c(!quick, TRUE, TRUE, TRUE)]) {
        funs[[Fun]] <- function(x, ...) {
            if (all(is.na(x))) NA
            else get(tolower(Fun))(x, ...)
        }
        for (sample in c("Un", "Adj")) {
            if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
                B[, paste(Fun, "Diff", sample, sep = ".")] <- sapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) abs(y[x, paste0("Diff.", sample)])), na.rm = TRUE))
                B[, paste(Fun, "V.Ratio", sample, sep = ".")] <- sapply(Brownames, function(x) if (B[x, "Type"]!="Contin.") NA else funs[[Fun]](sapply(balance.table.clusters.list, function(y) y[x, paste0("V.Ratio.", sample)]), na.rm = TRUE))
            }
        }
    }
    return(B)
}

#base.bal.tab.cont
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
samplesize.cont <- function(obj, method=c("matching", "weighting", "subclassification")) {
    #Computes sample size info. for unadjusted and adjusted samples. obj should be a cbps 
    # object, or a data.frame containing the variables "weights" and "treat" (and optionally 
    # "subclass") for each unit. method is what method the weights are to be used for. 
    # method="subclassification" is for subclass sample sizes only.
    method <- match.arg(method)
    if (method=="matching") {
        if (all(is.na(obj$weights))) {
            nn <- as.data.frame(matrix(0, ncol = 1, nrow = 1))
            nn[1, 1] <- length(obj$treat)
            dimnames(nn) <- list(c("All"), 
                                 c("Total"))
        }
        else {
            nn <- as.data.frame(matrix(0, ncol = 1, nrow = 3))
            nn[1, 1] <- length(obj$treat)
            nn[2, 1] <- sum(obj$weights > 0)
            nn[3, 1] <- sum(obj$weights == 0)
            dimnames(nn) <- list(c("All", "Matched", "Unmatched"), 
                                 c("Total"))
        }
        obs <- nn
        attr(obs, "tag") <- "Sample sizes:"
    }
    else if (method == "weighting") {
        w <- obj$weights
        obs <- as.data.frame(matrix(c(length(w), 
                                      (sum(w)^2)/sum(w^2)), 
                                    ncol=1, dimnames = list(c("Unadjusted", "Adjusted"), c("Total"))))
        attr(obs, "tag") <- "Effective sample sizes:"
    }
    else if (method == "subclassification") {
        stop("Subclassification is not yet surpported with continuous treatments.", call. = FALSE)
        # if (is.null(obj$subclass)) stop("obj must contain a vector of subclasses called \"subclass\".")
        # qbins <- length(levels(obj$subclass))
        # 
        # nn <- as.data.frame(matrix(0, 3, qbins))
        # 
        # dimnames(nn) <- list(c("Control", "Treated", "Total"), 
        #                      paste("Subclass", levels(obj$subclass)))
        # 
        # matched <- obj$weights!=0
        # k <- 0
        # for (i in levels(obj$subclass)) {
        #     qi <- obj$subclass[matched]==i & (!is.na(obj$subclass[matched]))
        #     qt <- obj$treat[matched][qi]
        #     if (sum(qt==1)<2|(sum(qt==0)<2)){
        #         if (sum(qt==1)<2)
        #             warning("Not enough treatment units in subclass ", i, call. = FALSE)
        #         else if (sum(qt==0)<2)
        #             warning("Not enough control units in subclass ", i, call. = FALSE)
        #     }
        #     k <- k + 1
        #     nn[, k] <- c(sum(qt==0), sum(qt==1), length(qt))
        # }
        # obs <- nn
        # attr(obs, "tag") <- "Sample sizes by subclass:"
    }
    return(obs)
}
balance.table.cont <- function(C, weights, treat, r.threshold = NULL, un = FALSE, no.adj = FALSE, types = NULL, quick = FALSE) {
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    
    #B=Balance frame
    Bnames <- c("Type", 
                "Corr.Un", 
                "Corr.Adj", 
                "R.Threshold")
    B <- as.data.frame(matrix(nrow=ncol(C), ncol=length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- varnames <- colnames(C)
    
    #Set var type (binary/continuous)
    if (length(types) > 0) B[,"Type"] <- types
    else B[,"Type"] <- get.types(C)
    
    #Correlations
    if (!no.adj) B[,"Corr.Adj"] <- sapply(C, w.r, y = treat, w = weights)
    if (!(!un && quick)) B[,"Corr.Un"] <- sapply(C, w.r, y = treat, w = NULL)
    
    if (!is.null(r.threshold)) {
        r.threshcheck <- ifelse(no.adj, "Corr.Un", "Corr.Adj")
        B[,"R.Threshold"] <- ifelse(B[,"Type"]=="Distance" | is.na(B[, r.threshcheck]), "", paste0(ifelse(abs(B[, r.threshcheck]) < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold))
    }
    
    return(B)
}
balance.table.cont.cluster.summary <- function(balance.table.clusters.list, r.threshold = NULL, cluster.summary = FALSE, no.adj = FALSE, types = NULL, quick = FALSE) {
    
    Brownames <- unique(do.call("c", lapply(balance.table.clusters.list, rownames)))
    cluster.functions <- c("Min", "Mean", "Median", "Max")
    Bcolnames <- c("Type", apply(expand.grid(cluster.functions, "Corr", c("Un", "Adj")), 1, paste, collapse = "."))
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
    names(B) <- Bcolnames
    
    if (length(types) > 0) B[,"Type"] <- types
    else B[,"Type"] <- unlist(sapply(Brownames, function(x) {u <- unique(sapply(balance.table.clusters.list, function(y) y[x, "Type"])); return(u[!is.na(u)])}), use.names = FALSE)
    
    funs <- structure(vector("list", length(cluster.functions)), names = cluster.functions)
    for (Fun in cluster.functions[c(!quick, TRUE, TRUE, TRUE)]) {
        funs[[Fun]] <- function(x, ...) {
            if (all(is.na(x))) NA
            else get(tolower(Fun))(x, ...)
        }
        for (sample in c("Un", "Adj")) {
            if (sample == "Un" || !no.adj) { #Only fill in Corr.Adj if no.adj = FALSE
                B[, paste(Fun, "Corr", sample, sep = ".")] <- sapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) abs(y[x, paste0("Corr.", sample)])), na.rm = TRUE))
            }
        }
    }
    return(B)
}

#love.plot
isColor <- function(x) {
    tryCatch(is.matrix(col2rgb(x)), 
             error = function(e) FALSE)
}

#print.bal.tab
round_df <- function(df, digits) {
    nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
    df[, nums] <- round(df[, nums], digits = digits)
    return(df)
}
replaceNA <- function(x) {
    x[is.na(x)] <- ""
    return(x)
}