#Think about sample sizes for clustered data. Think about subclass summaries with clustered data.

bal.tab <- function(...) UseMethod("bal.tab")
base.bal.tab <- function(object, weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    #object: A matchit, ps, or cbps object or a data.frame of treatment and weights; for computing sample sizes. obj can have clusters within it.
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
    #quick: if TRUE, values not displayed are not calculated (under construction)
    #*print options: un, disp.means, disp.v.ratio, disp.subclass (whether to show individual subclass balance)
    
    #Functions
    w.m <- function(x, w) {return(sum(x*w, na.rm=TRUE)/sum(w, na.rm=TRUE))}
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
            diff <- w.m(x[group==1 & ss], w = weights[group==1 & ss]) - w.m(x[group==0 & ss], w = weights[group==0 & ss])
        }
        else if (x.type=="Binary") {
            if      (binary=="raw") diff <- w.m(x[group==1 & ss], w = weights[group==1 & ss]) - w.m(x[group==0 & ss], w = weights[group==0 & ss])
            else if (binary=="std") {
                if (no.sub) diff <- std.diff(x, group, weights, denom = s.d.denom)
                else diff <- std.diff.subclass(x, group, weights, subclass, which.sub, denom = s.d.denom)
            }
            else diff <- NULL
        }
        else if (x.type=="Contin.") {
            if      (continuous=="raw") diff <- w.m(x[group==1 & ss], w = weights[group==1 & ss]) - w.m(x[group==0 & ss], w = weights[group==0 & ss])
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
    
    #Preparations
    if (length(unique(treat)) != 2) {
        stop("Treatment indicator must be a binary (0, 1) variable---i.e., treatment (1) or control (0)", call. = FALSE)
    }
    else {
        treat <- binarize(treat)
    }
    
    if (!is.null(m.threshold)) m.threshold <- abs(m.threshold)
    if (!is.null(v.threshold)) v.threshold <- max(v.threshold, 1/v.threshold)
    if (!is.null(v.threshold)) disp.v.ratio <- TRUE
    if (is.null(weights)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else no.adj <- FALSE
    
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
        C <- get.C(covs = covs, int = int, addl = addl, distance = distance, cluster = cluster)
        C.list <- structure(lapply(levels(cluster), function(x) C[cluster == x, ]), names = levels(cluster))
        types <- get.types(C)
        
        if (method == "subclassification") {
            stop("Subclassification with clusters is not supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab")
        }
        else {
            out[["Cluster.Balance"]] <- lapply(levels(cluster), function(i) balance.table(C = C.list[[i]], weights = weights[cluster == i], treat = treat[cluster == i], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, no.adj = no.adj, types = types, quick = quick))
            names(out[["Cluster.Balance"]]) <- levels(cluster)
            
            if (!(!cluster.summary && quick)) out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.table.clusters.list = out[["Cluster.Balance"]], 
                                                                      m.threshold = m.threshold, v.threshold = v.threshold, 
                                                                      cluster.summary = cluster.summary, no.adj = no.adj,
                                                                      types = types, quick = quick)
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
                                       disp.adj=!no.adj, 
                                       disp.subclass=disp.subclass,
                                       which.cluster=which.cluster,
                                       cluster.summary=cluster.summary,
                                       quick = quick)
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
                
                C <- get.C(covs = covs, int = int, addl = addl, distance = distance)
                
                if (length(list(...)$sub.by > 0)) sub.by <- list(...)$sub.by
                else sub.by <- call$sub.by
                
                out[["Subclass.Balance"]] <- balance.table.subclass(C, weights=weights, treat=treat, subclass=subclass, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, disp.means = disp.means, disp.v.ratio = disp.v.ratio, quick = quick)
                out[["Subclass.Observations"]] <- samplesize(object, method="subclassification")
                out[["Balance.Across.Subclass"]] <- balance.table.across.subclass(balance.table=balance.table(C, weights, treat, continuous, binary, s.d.denom, m.threshold, v.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, no.adj = T, quick = quick), 
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
                if (length(attr(out[["Subclass.Balance"]], "dont.disp.v.ratio")) > 0) {disp.v.ratio <- FALSE; v.threshold <- NULL}
                
                if (!is.null(v.threshold)) {
                    out[["Balanced.Variances.Subclass"]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][,"V.Threshold"])))
                    names(out[["Balanced.Variances.Subclass"]]) <- paste("Subclass", levels(subclass))
                    mivs.list <- lapply(levels(subclass), function(x) {
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
                                               disp.subclass=disp.subclass,
                                               quick = quick)
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
            
            C <- get.C(covs = covs, int = int, addl = addl, distance = distance, cluster = cluster)

            out[["Balance"]] <- balance.table(C, weights, treat, continuous, binary, s.d.denom, m.threshold, v.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, no.adj = no.adj, quick = quick)
            if (!is.null(m.threshold)) {
                m.threshcheck <- ifelse(no.adj, "Diff.Un", "Diff.Adj")
                out[["Balanced.Means"]] <- baltal(out[["Balance"]][,"M.Threshold"])
                out[["Max.Imbalance.Means"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], m.threshcheck, "M.Threshold")
            }
            if (all(is.na(out[["Balance"]][,"V.Ratio.Un"]))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
            if (!is.null(v.threshold)) {
                v.threshcheck <- ifelse(no.adj, "V.Ratio.Un", "V.Ratio.Adj")
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
                                           disp.adj=!no.adj,
                                           quick = quick)
            class(out) <- "bal.tab"
        }
    }
    
    #attr(out, "int") <- int
    return(out)
}
base.bal.tab.cont <- function(object, weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, r.threshold = NULL, un = FALSE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    #object: A list containing (at least) treatment status and weights for calculating sample sizes (currently unused)
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

    #Preparations
    if (!is.null(r.threshold)) r.threshold <- abs(r.threshold)
    if (is.null(weights)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else no.adj <- FALSE
    
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
        
        C <- get.C(covs = covs, int = int, addl = addl, distance = distance, cluster = cluster)
        C.list <- structure(lapply(levels(cluster), function(x) C[cluster == x, ]), names = levels(cluster))
        types <- get.types(C)
        
        if (method == "subclassification") {
            stop("Subclassification with clusters is not supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab") #add more for subclasses
        }
        else {
            out[["Cluster.Balance"]] <- lapply(levels(cluster), function(i) balance.table.cont(C = C.list[[i]], weights = weights[cluster == i], treat = treat[cluster == i], r.threshold = r.threshold, un = un, no.adj = no.adj, types = types, quick = quick))
            names(out[["Cluster.Balance"]]) <- levels(cluster)
            
            if (!(!cluster.summary && quick)) {
                out[["Cluster.Summary"]] <- balance.table.cont.cluster.summary(balance.table.clusters.list = out[["Cluster.Balance"]], 
                                                                                                             r.threshold = r.threshold, 
                                                                                                             cluster.summary = cluster.summary,
                                                                                                             no.adj = no.adj, types = types,
                                                                                                             quick = quick)
            }
            out <- out[!names(out) %in% "Cluster.Balance.Across.Subclass"]
            class(out) <- c("bal.tab.cont.cluster", "bal.cont.tab", "bal.tab.cluster", "bal.tab")
            
            out[["call"]] <- call
            out[["print.options"]] <- list(r.threshold=r.threshold, 
                                           un=un, 
                                           disp.adj=!no.adj, 
                                           which.cluster=which.cluster,
                                           cluster.summary=cluster.summary,
                                           quick = quick)
            class(out) <- c("bal.tab.cont.cluster", "bal.tab.cluster", "bal.tab.cont", "bal.tab")
        }
    }
    else {
        if (method %in% "subclassification") {
            stop("Subclassification not yet supported with continuous treatments.", call. = FALSE)
        }
        else {
            out.names <- c("Balance", "Balanced.Corr", 
                           "Max.Imbalance.Corr", 
                           "Observations", 
                           "call", "print.options")
            out <- vector("list", length(out.names))
            names(out) <- out.names
            
            C <- get.C(covs = covs, int = int, addl = addl, distance = distance, cluster = cluster)
            
            out[["Balance"]] <- balance.table.cont(C, weights, treat, r.threshold, un = un, no.adj = no.adj, quick = quick)
            if (!is.null(r.threshold)) {
                r.threshcheck <- ifelse(no.adj, "Corr.Un", "Corr.Adj")
                out[["Balanced.Corr"]] <- baltal(out[["Balance"]][, "R.Threshold"])
                out[["Max.Imbalance.Corr"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], r.threshcheck, "R.Threshold")
            }
            if (all(is.na(out[["Balance"]][,"Corr.Un"]))) {r.threshold <- NULL}
            out[["Observations"]] <- samplesize.cont(object, method = method)
            out[["call"]] <- call
            out[["print.options"]] <- list(r.threshold=r.threshold, 
                                           un=un, 
                                           disp.adj=!no.adj,
                                           quick = quick)
            class(out) <- c("bal.tab.cont", "bal.tab")
        }
    }
    
    attr(out, "int") <- int
    return(out)
}

bal.tab.matchit <- function(m, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    X <- x2base(m, cluster = cluster)
    
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, subclass=X$subclass, covs=X$covs, call=X$call, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.subclass=disp.subclass, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.ps <- function(ps, full.stop.method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
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
    X <- x2base(ps, 
                full.stop.method = full.stop.method, 
                s.d.denom = s.d.denom, 
                cluster = cluster)
    
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method="weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.Match <- function(M, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    X <- x2base(M, 
                formula = formula,
                data = data, 
                treat = treat,
                covs = covs,
                s.d.denom = s.d.denom,
                cluster = cluster)
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.formula <- function(formula, data, weights = NULL, distance = NULL, subclass = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, estimand = NULL, quick = FALSE, ...) {
    
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
                s.d.denom = s.d.denom,
                method = method,
                cluster = cluster,
                estimand = estimand)
    
    if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) {
        out <- base.bal.tab.cont(object=X$obj, weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = addl, r.threshold = r.threshold, un = un, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, subclass=X$subclass, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.subclass=disp.subclass, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.data.frame <- function(covs, treat, data = NULL, weights = NULL, distance = NULL, subclass = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, estimand = NULL, quick = FALSE, ...) {
    
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
    
    X <- x2base(covs,
                treat = treat,
                data = data,
                weights = weights,
                distance = distance,
                subclass = subclass,
                addl = addl,
                s.d.denom = s.d.denom,
                method = method,
                cluster = cluster,
                estimand = estimand)
    
    if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) {
        out <- base.bal.tab.cont(object=X$obj, weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = addl, r.threshold = r.threshold, un = un, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, subclass=X$subclass, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.subclass=disp.subclass, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.CBPS <- function(cbps, estimand, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
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
    X <- x2base(cbps, 
                estimand = estimand, 
                s.d.denom = s.d.denom,
                cluster = cluster)
    if (any(class(cbps) == "CBPSContinuous")) {
        out <- base.bal.tab.cont(object=X$obj, weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = addl, r.threshold = r.threshold, un = un, method = "weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    }
    else if (nlevels(as.factor(X$treat)) > 2) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method="weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.ebalance <- function(ebal, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    X <- x2base(ebal, 
                formula = formula,
                data = data, 
                treat = treat,
                covs = covs,
                cluster = cluster)
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
