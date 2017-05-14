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
    
    #Preparations
    if (length(unique(treat)) != 2) {
        stop("Treatment indicator must be a binary (0, 1) variable---i.e., treatment (1) or control (0)", call. = FALSE)
    }
    else {
        treat <- binarize(treat)
    }
    
    if (length(m.threshold) > 0) m.threshold <- abs(m.threshold)
    if (length(v.threshold) > 0) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (length(weights) == 0) {
        un <- TRUE
        no.adj <- TRUE
    }
    else no.adj <- FALSE
    
    #Actions
    if (nlevels(cluster) > 0) {
        out.names <- c("Cluster.Balance", 
                       "Cluster.Balance.Across.Subclass", 
                       "Cluster.Summary", "Observations",
                       "call", "print.options")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        # out[["Cluster.Balance"]] <- vector("list", length(levels(cluster)))
        # names(out[["Cluster.Balance"]]) <- levels(cluster)
        C <- get.C(covs = covs, int = int, addl = addl, distance = distance, cluster = cluster)
        C.list <- structure(lapply(levels(cluster), function(x) C[cluster == x, , drop = FALSE]), names = levels(cluster))
        types <- get.types(C)
        
        if (method == "subclassification") {
            stop("Subclassification with clusters is not yet supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab")
        }
        else {
            out[["Cluster.Balance"]] <- lapply(levels(cluster), function(c) setNames(list(balance.table(C = C.list[[c]], weights = weights[cluster == c], treat = treat[cluster == c], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, no.adj = no.adj, types = types, quick = quick),
                                                                                          samplesize(object, method = method, cluster = cluster, which.cluster = c)), 
                                                                                     c("Balance.Table", "Observations")))
            names(out[["Cluster.Balance"]]) <- levels(cluster)
            balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance.Table"]])
            observations <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Observations"]])
            
            if (!(!cluster.summary && quick)) out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.tables,
                                                                                                        no.adj = no.adj,
                                                                                                        quick = quick,
                                                                                                        types = types)
            if (all(sapply(balance.tables, function(x) !any(is.finite(x[,"V.Ratio.Un"]))))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
            out <- out[!names(out) %in% "Cluster.Balance.Across.Subclass"]
            out[["Observations"]] <- samplesize.across.clusters(observations)
            if (length(call) > 0) out[["call"]] <- call
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
            class(out) <- c("bal.tab.cluster", "bal.tab")
        }
        
    }
    else {
        if (method == "subclassification") {
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
                out[["Balance.Across.Subclass"]] <- balance.table.across.subclass(balance.table=balance.table(C, weights, treat, continuous, binary, s.d.denom, m.threshold, v.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, no.adj = TRUE, quick = quick), 
                                                                                  balance.table.subclass.list=out[["Subclass.Balance"]], 
                                                                                  subclass.obs=out[["Subclass.Observations"]], 
                                                                                  sub.by=sub.by, 
                                                                                  m.threshold=m.threshold, 
                                                                                  v.threshold=v.threshold)
                if (length(m.threshold) > 0) {
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
                
                if (length(v.threshold) > 0) {
                    out[["Balanced.Variances.Subclass"]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][,"V.Threshold"])))
                    names(out[["Balanced.Variances.Subclass"]]) <- paste("Subclass", levels(subclass))
                    mivs.list <- lapply(levels(subclass), function(x) {
                        mi <- max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][,"Type"]!="Distance", ], "V.Ratio.Adj", "V.Threshold")
                        return(data.frame(Variable = row.names(mi), mi))
                    } )      
                    mivs <- do.call("rbind", mivs.list)
                    
                    out[["Max.Imbalance.Variances.Subclass"]] <- data.frame(mivs, row.names = paste("Subclass", levels(subclass)))
                }
                if (length(call) > 0) out[["call"]] <- call
                out[["print.options"]] <- list(m.threshold=m.threshold, 
                                               v.threshold=v.threshold, 
                                               un=un, 
                                               disp.means=disp.means, 
                                               disp.v.ratio=disp.v.ratio, 
                                               disp.adj=length(weights) > 0, 
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
            if (length(m.threshold) > 0) {
                m.threshcheck <- ifelse(no.adj, "Diff.Un", "Diff.Adj")
                out[["Balanced.Means"]] <- baltal(out[["Balance"]][,"M.Threshold"])
                out[["Max.Imbalance.Means"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], m.threshcheck, "M.Threshold")
            }
            if (!any(is.finite(out[["Balance"]][,"V.Ratio.Un"]))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
            if (length(v.threshold) > 0) {
                v.threshcheck <- ifelse(no.adj, "V.Ratio.Un", "V.Ratio.Adj")
                out[["Balanced.Variances"]] <- baltal(out[["Balance"]][,"V.Threshold"])
                out[["Max.Imbalance.Variances"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], v.threshcheck, "V.Threshold")
            }
            out[["Observations"]] <- samplesize(object, method = method)
            if (length(call) > 0) out[["call"]] <- call
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
    
    #Preparations
    if (length(r.threshold) > 0) r.threshold <- abs(r.threshold)
    if (length(weights) == 0) {
        un <- TRUE
        no.adj <- TRUE
    }
    else no.adj <- FALSE
    
    #Actions
    if (nlevels(cluster) > 0) {
        out.names <- c("Cluster.Balance", 
                       "Cluster.Balance.Across.Subclass", 
                       "Cluster.Summary", "Observations",
                       "call", "print.options")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        out[["Cluster.Balance"]] <- vector("list", length(levels(cluster)))
        names(out[["Cluster.Balance"]]) <- levels(cluster)
        
        C <- get.C(covs = covs, int = int, addl = addl, distance = distance, cluster = cluster)
        C.list <- structure(lapply(levels(cluster), function(x) C[cluster == x, , drop = FALSE]), names = levels(cluster))
        types <- get.types(C)
        
        if (method == "subclassification") {
            stop("Subclassification with clusters is not yet supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab") #add more for subclasses
        }
        else {
            out[["Cluster.Balance"]] <- lapply(levels(cluster), function(c) setNames(list(balance.table.cont(C = C.list[[c]], weights = weights[cluster == c], treat = treat[cluster == c], r.threshold = r.threshold, un = un, no.adj = no.adj, types = types, quick = quick),
                                                                                          samplesize.cont(object, method = method, cluster = cluster, which.cluster = c)), 
                                                                                     c("Balance.Table", "Observations")))
            names(out[["Cluster.Balance"]]) <- levels(cluster)
            
            balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance.Table"]])
            observations <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Observations"]])
            
            if (!(!cluster.summary && quick)) out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.tables,
                                                                                                        no.adj = no.adj,
                                                                                                        quick = quick,
                                                                                                        types = types)
            out <- out[!names(out) %in% "Cluster.Balance.Across.Subclass"]
            out[["Observations"]] <- samplesize.across.clusters(observations)
            
            if (length(call) > 0) out[["call"]] <- call
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
            if (length(r.threshold) > 0) {
                r.threshcheck <- ifelse(no.adj, "Corr.Un", "Corr.Adj")
                out[["Balanced.Corr"]] <- baltal(out[["Balance"]][, "R.Threshold"])
                out[["Max.Imbalance.Corr"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], r.threshcheck, "R.Threshold")
            }
            if (!any(is.finite(out[["Balance"]][,"Corr.Un"]))) {r.threshold <- NULL}
            out[["Observations"]] <- samplesize.cont(object, method = method)
            if (length(call) > 0) out[["call"]] <- call
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
base.bal.tab.imp <- function(object, weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, quick = FALSE, ...) {
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
    #*print options: un, disp.means, disp.v.ratio, disp.subclass (whether to show individual subclass balance)
    
    #Preparations
    if (length(m.threshold) > 0) m.threshold <- abs(m.threshold)
    if (length(v.threshold) > 0) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (length(r.threshold) > 0) r.threshold <- abs(r.threshold)
    if (length(weights) == 0) {
        un <- TRUE
        no.adj <- TRUE
    }
    else no.adj <- FALSE
    
    #Setup output object
    out.names <- c("Imputation.Balance", 
                   "Cluster.Balance.Across.Imputations",
                   "Balance.Across.Imputations", 
                   "Observations", 
                   "call", "print.options")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    #Get list of bal.tabs for each imputation
    
    
    if (length(unique(treat)) > 2 && is.numeric(treat)) {#if continuous treatment
        out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) base.bal.tab.cont(object = object[imp==i, , drop = FALSE], weights = weights[imp==i], treat = treat[imp==i], distance = distance[imp==i, , drop = FALSE], subclass = subclass[imp==i], covs = covs[imp == i, , drop = FALSE], call = call, int = int, addl = addl[imp = i, , drop = FALSE], r.threshold = r.threshold, un = un, method, cluster = cluster[imp==i], which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick, ...))
    }
    else if (length(unique(treat)) > 2 && (is.factor(treat) || is.character(treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else {#if binary treatment
        out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) base.bal.tab(object = object[imp==i, , drop = FALSE], weights = weights[imp==i], treat = treat[imp==i], distance = distance[imp==i, , drop = FALSE], subclass = subclass[imp==i], covs = covs[imp==i, , drop = FALSE], call = call, int = int, addl = addl[imp==i, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.subclass = disp.subclass, method = method, cluster = cluster[imp==i], which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick, ...))
    }
    
    names(out[["Imputation.Balance"]]) <- levels(imp)
    
    #Create summary of lists
    
    if (!is.na(match("bal.tab.cluster", class(out[["Imputation.Balance"]][[1]])))) {
        if (!(!imp.summary && quick)) {
            out[["Cluster.Balance.Across.Imputations"]] <- lapply(levels(cluster), 
                                                                  function(c) setNames(list(balance.table.imp.summary(lapply(out[["Imputation.Balance"]], function(i) i[["Cluster.Balance"]][[c]][["Balance.Table"]]), 
                                                                                                                      no.adj = no.adj, quick = quick),
                                                                                            samplesize.across.imps(lapply(out[["Imputation.Balance"]], function(i) i[["Cluster.Balance"]][[c]][["Observations"]]))), 
                                                                                       c("Cluster.Balance", "Cluster.Observations")))
            names(out[["Cluster.Balance.Across.Imputations"]]) <- levels(cluster)
            balance.tables <- lapply(out[["Cluster.Balance.Across.Imputations"]], function(c) c[["Cluster.Balance"]])
            observations <- lapply(out[["Cluster.Balance.Across.Imputations"]], function(c) c[["Cluster.Observations"]])
            
            out[["Balance.Across.Imputations"]] <- balance.table.clust.imp.summary(balance.tables,
                                                                                   no.adj = no.adj,
                                                                                   quick = quick,
                                                                                   types = NULL)
            out[["Observations"]] <- samplesize.across.clusters(observations)
        }
        
        classes <- c("bal.tab.imp.cluster", "bal.tab.imp")
    }
    else {
        if (!is.na(match("bal.tab.subclass", class(out[["Imputation.Balance"]][[1]])))) {
            #Put something here
            stop("Subclassification cannot be used with multiply imputed data.", call. = FALSE)
        }
        else {
            if (!(!imp.summary && quick)) out[["Balance.Across.Imputations"]] <- balance.table.imp.summary(bal.tab.imp.list = out[["Imputation.Balance"]], 
                                                                                                           no.adj = no.adj,
                                                                                                           quick = quick,
                                                                                                           types = NULL)
            observations <- lapply(out[["Imputation.Balance"]], function(x) x[["Observations"]])
            
            out[["Observations"]] <- samplesize.across.imps(observations)
            classes <- "bal.tab.imp"
        }
    }
    
    if (length(call) > 0) out[["call"]] <- call
    out[["print.options"]] <- list(m.threshold=m.threshold,
                                   v.threshold=v.threshold,
                                   r.threshold=r.threshold,
                                   un=un, 
                                   disp.adj=!no.adj, 
                                   which.cluster=which.cluster,
                                   cluster.summary=cluster.summary,
                                   which.imp=which.imp,
                                   imp.summary=imp.summary,
                                   quick = quick,
                                   disp.means=disp.means, 
                                   disp.v.ratio=disp.v.ratio)
    class(out) <- unique(c(classes, sapply(out[["Imputation.Balance"]], class)))
    
    return(out)
}

bal.tab.matchit <- function(m, int = FALSE, distance = NULL, addl = NULL, data = NULL,  continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    blank.args <- sapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)))
    if (any(blank.args)) {
        for (arg.name in names(args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base.matchit(m, data = data, distance = distance, addl = addl, cluster = cluster)
    
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, subclass=X$subclass, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.subclass=disp.subclass, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.ps <- function(ps, full.stop.method, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    #full.stop.method = stopping rule/estimand from twang, e.g. "es.mean.att"; code to use first if none or incorrect is requested
    args <- as.list(environment())[-1]
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    blank.args <- sapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)))
    if (any(blank.args)) {
        for (arg.name in names(args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base.ps(ps, 
                   full.stop.method = full.stop.method, 
                   s.d.denom = s.d.denom, 
                   distance = distance,
                   addl = addl,
                   cluster = cluster)
    
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method="weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.Match <- function(M, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    blank.args <- sapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)))
    if (any(blank.args)) {
        for (arg.name in names(args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base.Match(M, 
                      formula = formula,
                      data = data, 
                      treat = treat,
                      covs = covs,
                      addl = addl,
                      distance = distance,
                      s.d.denom = s.d.denom,
                      cluster = cluster)
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.formula <- function(formula, data, weights = NULL, distance = NULL, subclass = NULL, match.strata = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, estimand = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    if (length(subclass) == 0) disp.subclass <- FALSE
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    blank.args <- sapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)))
    if (any(blank.args)) {
        for (arg.name in names(args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                assign(arg.name, NULL)
            }
        }
    }
    
    X <- x2base.formula(formula, 
                        data = data,
                        weights = weights,
                        distance = distance,
                        subclass = subclass,
                        match.strata = match.strata,
                        addl = addl,
                        s.d.denom = s.d.denom,
                        method = method,
                        cluster = cluster,
                        estimand = estimand,
                        imp = imp)

    if (length(X$imp) > 0) {
        out <- base.bal.tab.imp(object=X$obj, 
                                weights=X$weights, 
                                treat=X$treat, 
                                distance=X$distance, 
                                subclass=X$subclass, 
                                covs=X$covs, 
                                call=X$call, 
                                int=int, 
                                addl=X$addl, 
                                continuous=continuous, 
                                binary=binary, 
                                s.d.denom=X$s.d.denom, 
                                m.threshold=m.threshold, 
                                v.threshold=v.threshold, 
                                r.threshold=r.threshold,
                                un=un, 
                                disp.means=disp.means, 
                                disp.v.ratio=disp.v.ratio, 
                                disp.subclass=disp.subclass, 
                                method=X$method, 
                                cluster = X$cluster, 
                                which.cluster = which.cluster, 
                                cluster.summary = cluster.summary, 
                                imp = X$imp,
                                which.imp = which.imp,
                                imp.summary = imp.summary,
                                quick = quick)
    }
    else if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) {
        out <- base.bal.tab.cont(object=X$obj, weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = X$addl, r.threshold = r.threshold, un = un, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else {
        out <- base.bal.tab(object=X$obj, 
                            weights=X$weights, 
                            treat=X$treat, 
                            distance=X$distance, 
                            subclass=X$subclass, 
                            covs=X$covs, 
                            call=X$call, 
                            int=int, 
                            addl=X$addl, 
                            continuous=continuous, 
                            binary=binary, 
                            s.d.denom=X$s.d.denom, 
                            m.threshold=m.threshold, 
                            v.threshold=v.threshold, 
                            un=un, 
                            disp.means=disp.means, 
                            disp.v.ratio=disp.v.ratio, 
                            disp.subclass=disp.subclass, 
                            method=X$method, 
                            cluster = X$cluster, 
                            which.cluster = which.cluster, 
                            cluster.summary = cluster.summary, 
                            quick = quick)
    }
    return(out)
}
bal.tab.data.frame <- function(covs, treat, data = NULL, weights = NULL, distance = NULL, subclass = NULL, match.strata = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, estimand = NULL, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    if (length(subclass) == 0) disp.subclass <- FALSE
    
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    blank.args <- sapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)))
    if (any(blank.args)) {
        for (arg.name in names(args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                assign(arg.name, NULL)
            }
        }
    }
    
    X <- x2base.data.frame(covs,
                           treat = treat,
                           data = data,
                           weights = weights,
                           distance = distance,
                           subclass = subclass,
                           match.strata = match.strata,
                           addl = addl,
                           s.d.denom = s.d.denom,
                           method = method,
                           cluster = cluster,
                           estimand = estimand,
                           imp = imp)
    
    if (length(X$imp) > 0) {
        out <- base.bal.tab.imp(object=X$obj, 
                                weights=X$weights, 
                                treat=X$treat, 
                                distance=X$distance, 
                                subclass=X$subclass, 
                                covs=X$covs, 
                                call=X$call, 
                                int=int, 
                                addl=X$addl, 
                                continuous=continuous, 
                                binary=binary, 
                                s.d.denom=X$s.d.denom, 
                                m.threshold=m.threshold, 
                                v.threshold=v.threshold, 
                                r.threshold=r.threshold,
                                un=un, 
                                disp.means=disp.means, 
                                disp.v.ratio=disp.v.ratio, 
                                disp.subclass=disp.subclass, 
                                method=X$method, 
                                cluster = X$cluster, 
                                which.cluster = which.cluster, 
                                cluster.summary = cluster.summary, 
                                imp = X$imp,
                                which.imp = which.imp,
                                imp.summary = imp.summary,
                                quick = quick)
    }
    else if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) {
        out <- base.bal.tab.cont(object=X$obj, weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = X$addl, r.threshold = r.threshold, un = un, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else {
        out <- base.bal.tab(object=X$obj, 
                            weights=X$weights, 
                            treat=X$treat, 
                            distance=X$distance, 
                            subclass=X$subclass, 
                            covs=X$covs, 
                            call=X$call, 
                            int=int, 
                            addl=X$addl, 
                            continuous=continuous, 
                            binary=binary, 
                            s.d.denom=X$s.d.denom, 
                            m.threshold=m.threshold, 
                            v.threshold=v.threshold, 
                            un=un, 
                            disp.means=disp.means, 
                            disp.v.ratio=disp.v.ratio, 
                            disp.subclass=disp.subclass, 
                            method=X$method, 
                            cluster = X$cluster, 
                            which.cluster = which.cluster, 
                            cluster.summary = cluster.summary, 
                            quick = quick)
    }
    return(out)
}
bal.tab.CBPS <- function(cbps, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    blank.args <- sapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)))
    if (any(blank.args)) {
        for (arg.name in names(args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base.CBPS(cbps, 
                     s.d.denom = s.d.denom,
                     distance = distance,
                     addl = addl,
                     cluster = cluster)
    if (any(class(cbps) == "CBPSContinuous")) {
        out <- base.bal.tab.cont(object=X$obj, weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = X$addl, r.threshold = r.threshold, un = un, method = "weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    }
    else if (nlevels(as.factor(X$treat)) > 2) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method="weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.ebalance <- function(ebal, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    blank.args <- sapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)))
    if (any(blank.args)) {
        for (arg.name in names(args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base.ebalance(ebal, 
                         formula = formula,
                         data = data, 
                         treat = treat,
                         covs = covs,
                         distance = distance,
                         addl = addl,
                         cluster = cluster)
    out <- base.bal.tab(object=X$obj, weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.optmatch <- function(optmatch, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
    args <- c(as.list(environment()), list(...))[-1]
    
    #Adjustments to arguments
    args.with.choices <- names(formals()[-1])[sapply(formals()[-c(1, length(formals()))], function(x) length(x)>1)]
    for (i in seq_along(args.with.choices)) assign(args.with.choices[i], eval(parse(text=paste0("match.arg(", args.with.choices[i], ")"))))
    
    blank.args <- sapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)))
    if (any(blank.args)) {
        for (arg.name in names(args)[blank.args]) {
            if (identical(args[[arg.name]], quote(expr = ))) {
                assign(arg.name, NULL)
            }
        }
    }
    
    #Initializing variables
    X <- x2base.optmatch(optmatch, 
                         formula = formula,
                         data = data, 
                         treat = treat,
                         covs = covs,
                         distance = distance,
                         addl = addl,
                         cluster = cluster)
    out <- base.bal.tab(object=X$obj, 
                        weights=X$weights, 
                        treat=X$treat, 
                        distance=X$distance, 
                        covs=X$covs, 
                        call=X$call, 
                        int=int, 
                        addl=X$addl, 
                        continuous=continuous, 
                        binary=binary, 
                        s.d.denom=s.d.denom, 
                        m.threshold=m.threshold, 
                        v.threshold=v.threshold, 
                        un=un, 
                        disp.means=disp.means, 
                        disp.v.ratio=disp.v.ratio, 
                        method=X$method, 
                        cluster = X$cluster, 
                        which.cluster = which.cluster, 
                        cluster.summary = cluster.summary, 
                        quick = quick)
    return(out)
}