bal.tab <- function(...) UseMethod("bal.tab")
base.bal.tab <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, disp.bal.tab = TRUE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, s.weights = NULL, discarded = NULL, quick = FALSE, pooled.sds = NULL, ...) {
    
    #Preparations
    args <- list(...)
    if (length(unique(treat)) != 2) {
        stop("Treatment indicator must be a binary (0, 1) variable---i.e., treatment (1) or control (0)", call. = FALSE)
    }
    else if (is.factor(treat) || is.character(treat)) {
        if (is.factor(treat)) treat.names <- levels(treat)
        else treat.names <- unique(treat)
    }
    else {
        treat.names <- c("Control", "Treated")
    }
    treat <- binarize(treat)
    
    if (length(m.threshold) > 0) m.threshold <- abs(m.threshold)
    if (length(v.threshold) > 0) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (length(ks.threshold) == 0 && (length(args$k.threshold) > 0)) {
        ks.threshold <- args$k.threshold
    }
    if (length(ks.threshold) > 0) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (length(weights) == 0 && length(subclass) == 0) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (length(weights) > 0 && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (length(s.weights) == 0) {
        s.weights <- rep(1, length(treat))
    }
    
    #Actions
    if (nunique.gt(cluster, 1)) {
        out.names <- c("Cluster.Balance", 
                       "Cluster.Balance.Across.Subclass", 
                       "Cluster.Summary", "Observations",
                       "call", "print.options")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        C <- get.C(covs = covs, int = int, addl = addl, distance = distance, cluster = cluster)
        C.list <- setNames(lapply(levels(cluster), function(x) C[cluster == x, , drop = FALSE]), 
                           levels(cluster))
        types <- get.types(C)
        
        if (length(method) == 1 && method == "subclassification") {
            stop("Subclassification with clusters is not yet supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab")
        }
        else {
            out[["Cluster.Balance"]] <- setNames(lapply(levels(cluster), function(c) setNames(list(balance.table(C = C.list[[c]], weights = weights[cluster == c, , drop = FALSE], treat = treat[cluster == c], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, s.weights = s.weights[cluster == c], no.adj = no.adj, types = types, quick = quick),
                                                                                                   samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, cluster = cluster, which.cluster = c, discarded = discarded)), 
                                                                                              c("Balance", "Observations"))),
                                                 levels(cluster))
            balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance"]])
            observations <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Observations"]])
            
            if (!(!cluster.summary && quick)) out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.tables,
                                                                                                        weight.names = names(weights),
                                                                                                        no.adj = no.adj,
                                                                                                        quick = quick,
                                                                                                        types = types)
            if (all(sapply(balance.tables, function(x) !attr(x, "disp")["v"]))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
            if (all(sapply(balance.tables, function(x) !attr(x, "disp")["ks"]))) {disp.ks.ratio <- FALSE; ks.threshold <- NULL}
            out <- out[!names(out) %in% "Cluster.Balance.Across.Subclass"]
            out[["Observations"]] <- samplesize.across.clusters(observations)
            if (length(call) > 0) out[["call"]] <- call
            out[["print.options"]] <- list(m.threshold=m.threshold, 
                                           v.threshold=v.threshold,
                                           ks.threshold=ks.threshold,
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.means=disp.means, 
                                           disp.v.ratio=disp.v.ratio, 
                                           disp.ks=disp.ks, 
                                           disp.adj=!no.adj, 
                                           disp.subclass=disp.subclass,
                                           disp.bal.tab = disp.bal.tab, 
                                           which.cluster=which.cluster,
                                           cluster.summary=cluster.summary,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights),
                                           treat.names = treat.names)
            class(out) <- c("bal.tab.cluster", "bal.tab")
        }
        
    }
    else {
        if (length(method) == 1 && method == "subclassification") {
            if (length(subclass) > 0) {
                out.names <- c("Subclass.Balance", "Balance.Across.Subclass", 
                               "Balanced.Means.Subclass", "Max.Imbalance.Means.Subclass", 
                               "Balanced.Variances.Subclass", "Max.Imbalance.Variances.Subclass", 
                               "Balanced.KS.Subclass", "Max.Imbalance.KS.Subclass", 
                               "Subclass.Observations", "call", "print.options")
                out <- vector("list", length(out.names))
                names(out) <- out.names
                
                C <- get.C(covs = covs, int = int, addl = addl, distance = distance)
                
                if (length(list(...)$sub.by > 0)) sub.by <- list(...)$sub.by
                else sub.by <- call$sub.by
                out[["Subclass.Balance"]] <- balance.table.subclass(C, weights=weights[[1]], treat=treat, subclass=subclass, continuous=continuous, binary=binary, s.d.denom=s.d.denom[1], m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold = ks.threshold, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, quick = quick)
                out[["Subclass.Observations"]] <- samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
                out[["Balance.Across.Subclass"]] <- balance.table.across.subclass(balance.table=balance.table(C, weights[[1]], treat, continuous, binary, s.d.denom[1], m.threshold, v.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, no.adj = TRUE, quick = quick), 
                                                                                  balance.table.subclass.list=out[["Subclass.Balance"]], 
                                                                                  subclass.obs=out[["Subclass.Observations"]], 
                                                                                  sub.by=sub.by, 
                                                                                  m.threshold=m.threshold, 
                                                                                  v.threshold=v.threshold, 
                                                                                  ks.threshold=ks.threshold,
                                                                                  s.d.denom = s.d.denom[1])
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
                
                if (length(attr(out[["Subclass.Balance"]], "dont.disp.ks")) > 0) {disp.ks <- FALSE; ks.threshold <- NULL}
                if (length(ks.threshold) > 0) {
                    out[["Balanced.KS.Subclass"]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][["KS.Threshold"]])))
                    names(out[["Balanced.KS.Subclass"]]) <- paste("Subclass", levels(subclass))
                    miks.list <- lapply(levels(subclass), function(x) {
                        mi <- max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][["Type"]]!="Distance", ], "KS.Adj", "KS.Threshold")
                        return(data.frame(Variable = row.names(mi), mi))
                    } )      
                    miks <- do.call("rbind", miks.list)
                    
                    out[["Max.Imbalance.KS.Subclass"]] <- data.frame(miks, row.names = paste("Subclass", levels(subclass)))
                }
                
                if (length(call) > 0) out[["call"]] <- call
                out[["print.options"]] <- list(m.threshold=m.threshold, 
                                               v.threshold=v.threshold, 
                                               ks.threshold=ks.threshold, 
                                               imbalanced.only = imbalanced.only,
                                               un=un, 
                                               disp.means=disp.means, 
                                               disp.v.ratio=disp.v.ratio, 
                                               disp.ks=disp.ks, 
                                               disp.adj=!no.adj, 
                                               disp.subclass=disp.subclass,
                                               disp.bal.tab = disp.bal.tab, 
                                               quick = quick,
                                               treat.names = treat.names)
                class(out) <- c("bal.tab.subclass", "bal.tab")
            }
            else stop("Method specified as subclassification, but no subclasses were specified.", call. = FALSE)
        }
        else {
            out.names <- c("Balance", "Balanced.Means", 
                           "Max.Imbalance.Means", "Balanced.Variances", 
                           "Max.Imbalance.Variances", "Balanced.KS", 
                           "Max.Imbalance.KS", "Observations", 
                           "call", "print.options")
            out <- vector("list", length(out.names))
            names(out) <- out.names
            
            C <- get.C(covs = covs, int = int, addl = addl, distance = distance)
            
            out[["Balance"]] <- balance.table(C, weights, treat, continuous, binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, s.weights = s.weights, no.adj = no.adj, quick = quick, pooled.sds = pooled.sds)
            
            #Ensure comaptible with multiple weights
            if (length(m.threshold) > 0) {
                if (no.adj) {
                    out[["Balanced.Means"]] <- baltal(out[["Balance"]][,"M.Threshold.Un"])
                    out[["Max.Imbalance.Means"]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], "Diff.Un", "M.Threshold.Un")
                }
                else if (ncol(weights) == 1) {
                    out[["Balanced.Means"]] <- baltal(out[["Balance"]][,"M.Threshold"])
                    out[["Max.Imbalance.Means"]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], "Diff.Adj", "M.Threshold")
                }
                else if (ncol(weights) > 1) {
                    out[["Balanced.Means"]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][[paste0("M.Threshold.", x)]]))),
                                                        names(weights))
                    out[["Max.Imbalance.Means"]] <- cbind(Weights = names(weights),
                                                          do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], paste0("Diff.", x), paste0("M.Threshold.", x)),
                                                                                                                       c("Variable", "Diff", "M.Threshold")))),
                                                          stringsAsFactors = FALSE)
                }
            }
            if (!attr(out[["Balance"]], "disp")["v"]) {disp.v.ratio <- FALSE; v.threshold <- NULL}
            if (length(v.threshold) > 0) {
                if (no.adj) {
                    out[["Balanced.Variances"]] <- baltal(out[["Balance"]][["V.Threshold.Un"]])
                    out[["Max.Imbalance.Variances"]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], "V.Ratio.Un", "V.Threshold.Un")
                    
                }
                else if (ncol(weights) == 1) {
                    out[["Balanced.Variances"]] <- baltal(out[["Balance"]][["V.Threshold"]])
                    out[["Max.Imbalance.Variances"]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], "V.Ratio.Adj", "V.Threshold")
                }
                else {
                    out[["Balanced.Variances"]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][[paste0("V.Threshold.", x)]]))),
                                                            names(weights))
                    out[["Max.Imbalance.Variances"]] <- cbind(Weights = names(weights),
                                                              do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], paste0("V.Ratio.", x), paste0("V.Threshold.", x)),
                                                                                                                           c("Variable", "V.Ratio", "V.Threshold")))),
                                                              stringsAsFactors = FALSE)
                }
            }
            if (!attr(out[["Balance"]], "disp")["ks"]) {disp.ks <- FALSE; ks.threshold <- NULL}
            if (length(ks.threshold) > 0) {
                if (no.adj) {
                    out[["Balanced.KS"]] <- baltal(out[["Balance"]][["KS.Threshold.Un"]])
                    out[["Max.Imbalance.KS"]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], "KS.Un", "KS.Threshold.Un")
                    
                }
                else if (ncol(weights) == 1) {
                    out[["Balanced.KS"]] <- baltal(out[["Balance"]][,"KS.Threshold"])
                    out[["Max.Imbalance.KS"]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], "KS.Adj", "KS.Threshold")
                }
                else {
                    out[["Balanced.KS"]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][[paste0("KS.Threshold.", x)]]))),
                                                     names(weights))
                    out[["Max.Imbalance.KS"]] <- cbind(Weights = names(weights),
                                                       do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], paste0("KS.", x), paste0("KS.Threshold.", x)),
                                                                                                                    c("Variable", "KS", "KS.Threshold")))),
                                                       stringsAsFactors = FALSE)
                }
            }
            
            out[["Observations"]] <- samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded, treat.names = treat.names)
            if (length(call) > 0) out[["call"]] <- call
            out[["print.options"]] <- list(m.threshold=m.threshold, 
                                           v.threshold=v.threshold, 
                                           ks.threshold=ks.threshold, 
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.means=disp.means, 
                                           disp.v.ratio=disp.v.ratio, 
                                           disp.ks=disp.ks, 
                                           disp.adj=!no.adj,
                                           disp.bal.tab = disp.bal.tab, 
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights),
                                           treat.names = treat.names)
            class(out) <- "bal.tab"
        }
    }
    
    #attr(out, "int") <- int
    return(out)
}
base.bal.tab.cont <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, r.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.subclass = FALSE, disp.bal.tab = TRUE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, s.weights = NULL, discarded = NULL, quick = FALSE, ...) {
    
    #Preparations
    if (length(r.threshold) > 0) {
        r.threshold <- abs(r.threshold)
        if (r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            r.threshold <- NULL
        }
    }
    if (length(weights) == 0 && length(subclass) == 0) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (length(weights) > 0 && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (length(s.weights) == 0) {
        s.weights <- rep(1, length(treat))
    }    
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
        
        if (length(method) == 1 && method == "subclassification") {
            stop("Subclassification with clusters is not yet supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab") #add more for subclasses
        }
        else {
            out[["Cluster.Balance"]] <- lapply(levels(cluster), function(c) setNames(list(balance.table.cont(C = C.list[[c]], weights = weights[cluster == c, , drop = FALSE], treat = treat[cluster == c], r.threshold = r.threshold, un = un, s.weights = s.weights[cluster == c], no.adj = no.adj, types = types, quick = quick),
                                                                                          samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, cluster = cluster, which.cluster = c, discarded = discarded)), 
                                                                                     c("Balance", "Observations")))
            names(out[["Cluster.Balance"]]) <- levels(cluster)
            
            balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance"]])
            observations <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Observations"]])
            
            if (!(!cluster.summary && quick)) out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.tables,
                                                                                                        weight.names = names(weights),
                                                                                                        no.adj = no.adj,
                                                                                                        quick = quick,
                                                                                                        types = types)
            out <- out[!names(out) %in% "Cluster.Balance.Across.Subclass"]
            out[["Observations"]] <- samplesize.across.clusters(observations)
            
            if (length(call) > 0) out[["call"]] <- call
            out[["print.options"]] <- list(r.threshold=r.threshold, 
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.adj=!no.adj, 
                                           disp.bal.tab = disp.bal.tab,
                                           which.cluster=which.cluster,
                                           cluster.summary=cluster.summary,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights))
            class(out) <- c("bal.tab.cont.cluster", "bal.tab.cluster", "bal.tab.cont", "bal.tab")
        }
        
    }
    else {
        if (length(method) == 1 && method == "subclassification") {
            #stop("Subclassification not yet supported with continuous treatments.", call. = FALSE)
            if (length(subclass)>0) {
                out.names <- c("Subclass.Balance", 
                               "Balanced.Corr.Subclass", "Max.Imbalance.Corr.Subclass", 
                               "Subclass.Observations", "call", "print.options")
                out <- vector("list", length(out.names))
                names(out) <- out.names
                
                C <- get.C(covs = covs, int = int, addl = addl, distance = distance)
                
                # if (length(list(...)$sub.by > 0)) sub.by <- list(...)$sub.by
                # else sub.by <- call$sub.by
                
                out[["Subclass.Balance"]] <- balance.table.subclass.cont(C, weights=weights[[1]], treat=treat, subclass=subclass, r.threshold=r.threshold, quick = quick)
                out[["Subclass.Observations"]] <- samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
                #out[["Balance.Across.Subclass"]]
                if (length(r.threshold) > 0) {
                    out[["Balanced.Corr.Subclass"]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][,"R.Threshold"])))
                    names(out[["Balanced.Corr.Subclass"]]) <- paste("Subclass", levels(subclass))
                    mirs.list <- lapply(levels(subclass), function(x) {
                        mi <- max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][,"Type"]!="Distance", ], "Corr.Adj", "R.Threshold")
                        return(data.frame(Variable = row.names(mi), mi))
                    } )
                    mirs <- do.call("rbind", mirs.list)
                    out[["Max.Imbalance.Corr.Subclass"]] <- data.frame(mirs, row.names = paste("Subclass", levels(subclass)))
                }
                
                if (length(call) > 0) out[["call"]] <- call
                out[["print.options"]] <- list(r.threshold=r.threshold, 
                                               imbalanced.only = imbalanced.only,
                                               un=un,
                                               disp.adj=!no.adj, 
                                               disp.subclass=disp.subclass,
                                               disp.bal.tab = disp.bal.tab,
                                               quick = quick)
                class(out) <- c("bal.tab.subclass.cont", "bal.tab.subclass", "bal.tab.cont", "bal.tab")
            }
            else stop("Method specified as subclassification, but no subclasses were specified.", call. = FALSE)
            
        }
        else {
            out.names <- c("Balance", "Balanced.Corr", 
                           "Max.Imbalance.Corr", 
                           "Observations", 
                           "call", "print.options")
            out <- vector("list", length(out.names))
            names(out) <- out.names
            
            C <- get.C(covs = covs, int = int, addl = addl, distance = distance, cluster = cluster)
            
            out[["Balance"]] <- balance.table.cont(C, weights, treat, r.threshold, un = un, s.weights = s.weights, no.adj = no.adj, quick = quick)
            if (length(r.threshold) > 0) {
                if (no.adj) {
                    out[["Balanced.Corr"]] <- baltal(out[["Balance"]][["R.Threshold"]])
                    out[["Max.Imbalance.Corr"]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], "Corr.Un", "R.Threshold")
                }
                else if (ncol(weights) == 1) {
                    out[["Balanced.Corr"]] <- baltal(out[["Balance"]][["R.Threshold"]])
                    out[["Max.Imbalance.Corr"]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], "Corr.Adj", "R.Threshold")
                }
                else if (ncol(weights) > 1) {
                    out[["Balanced.Corr"]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][[paste0("R.Threshold.", x)]]))),
                                                       names(weights))
                    out[["Max.Imbalance.Corr"]] <- cbind(Weights = names(weights),
                                                         do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", ], paste0("Corr.", x), paste0("R.Threshold.", x)),
                                                                                                                      c("Variable", "Corr", "R.Threshold")))),
                                                         stringsAsFactors = FALSE)
                }
            }
            if (!any(is.finite(out[["Balance"]][["Corr.Un"]]))) {r.threshold <- NULL}
            out[["Observations"]] <- samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
            if (length(call) > 0) out[["call"]] <- call
            out[["print.options"]] <- list(r.threshold=r.threshold, 
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.adj=!no.adj,
                                           disp.bal.tab = disp.bal.tab,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights))
            class(out) <- c("bal.tab.cont", "bal.tab")
        }
    }
    
    attr(out, "int") <- int
    return(out)
}
base.bal.tab.imp <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, disp.bal.tab = TRUE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, s.weights = NULL, discarded = NULL, quick = FALSE, ...) {
    
    #Preparations
    if (length(m.threshold) > 0) m.threshold <- abs(m.threshold)
    if (length(v.threshold) > 0) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (length(ks.threshold) > 0) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (length(r.threshold) > 0) {
        r.threshold <- abs(r.threshold)
        if (r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            r.threshold <- NULL
        }
    }
    if (length(weights) == 0) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (length(s.weights) == 0) {
        s.weights <- rep(1, length(treat))
    }    
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
        out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) base.bal.tab.cont(weights = weights[imp==i, , drop  = FALSE], treat = treat[imp==i], distance = distance[imp==i, , drop = FALSE], subclass = subclass[imp==i], covs = covs[imp == i, , drop = FALSE], call = call, int = int, addl = addl[imp = i, , drop = FALSE], r.threshold = r.threshold, imbalanced.only = imbalanced.only, un = un, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster[imp==i], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[imp==i], discarded = discarded[imp==i], quick = quick, ...))
    }
    else if (length(unique(treat)) > 2 && (is.factor(treat) || is.character(treat))) {
        stop("Multiple categorical treaments are not yet supported with multiply imputed data.", call. = FALSE)
    }
    else {#if binary treatment
        out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) base.bal.tab(weights = weights[imp==i, , drop = FALSE], treat = treat[imp==i], distance = distance[imp==i, , drop = FALSE], subclass = subclass[imp==i], covs = covs[imp==i, , drop = FALSE], call = call, int = int, addl = addl[imp==i, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster[imp==i], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[imp==i], discarded = discarded[imp==i], quick = quick, ...))
    }
    
    names(out[["Imputation.Balance"]]) <- levels(imp)
    
    #Create summary of lists
    
    if (!is.na(match("bal.tab.cluster", class(out[["Imputation.Balance"]][[1]])))) {
        if (!(!imp.summary && quick)) {
            out[["Cluster.Balance.Across.Imputations"]] <- lapply(levels(cluster), 
                                                                  function(c) setNames(list(balance.table.imp.summary(lapply(out[["Imputation.Balance"]], function(i) i[["Cluster.Balance"]][[c]][["Balance"]]), 
                                                                                                                      weight.names = names(weights),
                                                                                                                      no.adj = no.adj, quick = quick),
                                                                                            samplesize.across.imps(lapply(out[["Imputation.Balance"]], function(i) i[["Cluster.Balance"]][[c]][["Observations"]]))), 
                                                                                       c("Cluster.Balance", "Cluster.Observations")))
            names(out[["Cluster.Balance.Across.Imputations"]]) <- levels(cluster)
            balance.tables <- lapply(out[["Cluster.Balance.Across.Imputations"]], function(c) c[["Cluster.Balance"]])
            observations <- lapply(out[["Cluster.Balance.Across.Imputations"]], function(c) c[["Cluster.Observations"]])
            
            out[["Balance.Across.Imputations"]] <- balance.table.clust.imp.summary(balance.tables,
                                                                                   weight.names = names(weights),
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
                                                                                                           weight.names = names(weights),
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
                                   ks.threshold=ks.threshold,
                                   r.threshold=r.threshold,
                                   imbalanced.only = imbalanced.only,
                                   un=un, 
                                   disp.adj=!no.adj, 
                                   which.cluster=which.cluster,
                                   cluster.summary=cluster.summary,
                                   which.imp=which.imp,
                                   imp.summary=imp.summary,
                                   quick = quick,
                                   disp.means=disp.means, 
                                   disp.v.ratio=disp.v.ratio, 
                                   disp.ks=disp.ks,
                                   disp.bal.tab = disp.bal.tab,
                                   nweights = ifelse(no.adj, 0, ncol(weights)),
                                   weight.names = names(weights))
    class(out) <- unique(c(classes, sapply(out[["Imputation.Balance"]], class)))
    
    return(out)
}
base.bal.tab.multi <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, disp.bal.tab = TRUE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = TRUE, s.weights = NULL, discarded = NULL, quick = FALSE, ...) {
    #Preparations
    
    if (length(m.threshold) > 0) m.threshold <- abs(m.threshold)
    if (length(v.threshold) > 0) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (length(ks.threshold) > 0) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (length(weights) == 0 && length(subclass) == 0) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (length(weights) > 0 && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (length(s.weights) == 0) {
        s.weights <- rep(1, length(treat))
    }
    
    #Treat is a factor variable of 3+ levels
    if (length(focal) == 0) {
        if (pairwise) treat.combinations <- combn(levels(treat), 2, list)
        else treat.combinations <- lapply(levels(treat), function(x) c(x, "Others"))
    }
    else if (length(focal) == 1) {
        if (is.numeric(focal)) {
            focal <- levels(treat)[focal]
        }
        if (is.character(focal)) {
            treat <- relevel(treat, focal)
        }
        else {
            stop("focal must be the name or index of the focal treatment group.", call. = FALSE)
        }
        treat.combinations <- lapply(levels(treat)[levels(treat) != focal], function(x) rev(c(focal, x)))
        pairwise <- TRUE
    }
    else stop("focal must be a vector of length 1 containing the name or index of the focal treatment group.", call. = FALSE)
    treat.names <- levels(treat)
    
    if (length(cluster) > 0) {
        stop("Clusters are not yet supported with multiple categorical treatments.", call. = FALSE)
    }
    else {
        #Setup output object
        out.names <- c("Pair.Balance", 
                       "Balance.Across.Pairs", 
                       "Observations", 
                       "call", "print.options")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        if (s.d.denom == "pooled") {
            C <- get.C(covs = covs, int = int, addl = addl, distance = distance)
            pooled.sds <- rowMeans(sapply(levels(treat), function(t) sqrt(col.w.v(C[treat == t, , drop = FALSE], 
                                                           s.weights[treat == t]))))
            #pooled.sds <- sqrt(col.w.v(C, s.weights)) #How twang does it
        }
        else pooled.sds <- NULL
        
        if (pairwise || length(focal) > 0) {
            balance.tables <- lapply(treat.combinations, function(t) base.bal.tab(weights = weights[treat %in% t, , drop = FALSE], treat = factor(treat[treat %in% t], t), distance = distance[treat %in% t, , drop = FALSE], subclass = subclass[treat %in% t], covs = covs[treat %in% t, , drop = FALSE], call = NULL, int = int, addl = addl[treat %in% t, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster[treat %in% t], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[treat %in% t], discarded = discarded[treat %in% t], quick = quick, pooled.sds = pooled.sds, ...))
        }
        else {
            if (any(treat.names == "Others")) stop ("\"Others\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
            balance.tables <- lapply(treat.combinations, function(t) {
                treat_ <- factor(treat, levels = c(levels(treat), "Others"))
                treat_[treat_ != t[1]] <- "Others"
                treat_ <- factor(treat_, rev(t))
                base.bal.tab(weights = weights, treat = treat_, distance = distance, subclass = subclass, covs = covs, call = NULL, int = int, addl = addl, continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights, discarded = discarded, quick = quick, pooled.sds = pooled.sds, ...)
            })
        }
        for (i in seq_along(balance.tables)) {
            names(balance.tables)[i] <- paste(rev(treat.combinations[[i]]), collapse = " vs. ")
        }
        
        out[["Pair.Balance"]] <- balance.tables
        
        out[["Observations"]] <- samplesize.multi(balance.tables, treat.names, focal)
        
        out[["Balance.Across.Pairs"]] <- balance.table.multi.summary(balance.tables, 
                                                                     weight.names = names(weights),
                                                                     m.threshold = m.threshold,
                                                                     v.threshold = v.threshold,
                                                                     ks.threshold = ks.threshold,
                                                                     no.adj = no.adj,
                                                                     quick = quick,
                                                                     types = NULL)
        
        if (length(call) > 0) out[["call"]] <- call
        
        out[["print.options"]] <- list(m.threshold=m.threshold,
                                       v.threshold=v.threshold,
                                       ks.threshold=ks.threshold,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.adj=!no.adj, 
                                       which.cluster=which.cluster,
                                       cluster.summary=cluster.summary,
                                       quick = quick,
                                       disp.means=disp.means, 
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks,
                                       disp.bal.tab = disp.bal.tab,
                                       nweights = ifelse(no.adj, 0, ncol(weights)),
                                       weight.names = names(weights),
                                       treat.names = treat.names,
                                       which.treat = which.treat,
                                       multi.summary = multi.summary,
                                       pairwise = pairwise)
        
        class(out) <- c("bal.tab.multi", "bal.tab")
    }
    return(out)
    
}
base.bal.tab.msm <- function(weights, treat.list, distance.list = NULL, subclass = NULL, covs.list, call = NULL, int = FALSE, addl.list = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.bal.tab = TRUE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, s.weights = NULL, discarded = NULL, quick = FALSE, ...) {
    #One vector of weights
    #treat.list should be a df/list of treatment vectors, one for each time period
    #cov.list should be a list of covariate data.frames, one for each time period; 
    #   should include all covs from previous time points, but no treatment statuses
    
    #Preparations
    args <- list(...)
    
    if (length(m.threshold) > 0) m.threshold <- abs(m.threshold)
    if (length(v.threshold) > 0) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (length(ks.threshold) == 0 && (length(args$k.threshold) > 0)) {
        ks.threshold <- args$k.threshold
    }
    if (length(ks.threshold) > 0) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (length(r.threshold) > 0) {
        r.threshold <- abs(r.threshold)
        if (r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            r.threshold <- NULL
        }
    }
    if (length(weights) == 0 && length(subclass) == 0) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (length(weights) > 0 && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (length(s.weights) == 0) {
        s.weights <- rep(1, length(treat.list[[1]]))
    }
    
    if (nunique.gt(cluster, 1)) {
        stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
    }
    else {
        #Setup output object
        out.names <- c("Time.Balance", 
                       "Balance.Across.Times", 
                       "Observations", 
                       "call", "print.options")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        out[["Time.Balance"]] <- vector("list", length(covs.list))
        
        treat.type <- sapply(treat.list, function(x) {
            if (nunique.gt(x, 2)) {
                 if (is.numeric(x)) "continuous"
                else if (is.factor(x) || is.character(x)) "multinomial"
            }
            else if (nunique.gt(x, 1)) {
                "binary"
            }
            else {
                stop("All treatments must have at least 2 unique values.", call. = FALSE)
            }
        })

        #Get list of bal.tabs for each time period
        out[["Time.Balance"]] <- lapply(seq_along(treat.list), function(ti) {
            if (treat.type[ti] == "continuous") {
                out_ <- base.bal.tab.cont(weights = weights, treat = treat.list[[ti]], distance = distance.list[[ti]], subclass = NULL, covs = covs.list[[ti]], call = NULL, int = int, addl = addl.list[[ti]], r.threshold = r.threshold, imbalanced.only = imbalanced.only, un = un, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights, discarded = discarded, quick = quick, ...)
            }
            else if (treat.type[ti] == "multinomial") {
                #stop("Multiple categorical treatments are not yet supported with longitudinal treatments.", call. = FALSE)
                out_ <- base.bal.tab.multi(weights = weights, treat = treat.list[[ti]], distance=distance.list[[ti]], 
                                          covs=covs.list[[ti]], call=NULL, int=int, addl=addl.list[[ti]], 
                                          continuous=continuous, binary=binary, s.d.denom=s.d.denom, 
                                          m.threshold=m.threshold, v.threshold=v.threshold, 
                                          ks.threshold=ks.threshold, 
                                          imbalanced.only = imbalanced.only,
                                          un=un, 
                                          disp.bal.tab = disp.bal.tab, 
                                          disp.means=disp.means,
                                          disp.v.ratio=disp.v.ratio, 
                                          disp.ks=disp.ks, 
                                          method=method, 
                                          cluster = cluster, which.cluster = which.cluster, 
                                          cluster.summary = cluster.summary, pairwise = pairwise, focal = NULL,
                                          which.treat = which.treat, multi.summary = multi.summary,
                                          s.weights = s.weights, quick = quick)
            }
            else if (treat.type[ti] == "binary") {
                out_ <- base.bal.tab(weights = weights, treat = treat.list[[ti]], distance = distance.list[[ti]], subclass = NULL, covs = covs.list[[ti]], call = NULL, int = int, addl = addl.list[[ti]], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = FALSE, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights, discarded = discarded, quick = quick, ...)
            }
            else stop("Each treatment must be binary, multinomial, or continuous.", call. = FALSE)
            
            return(out_)
        })

        if (length(names(treat.list)) == length(treat.list)) {
            names(out[["Time.Balance"]]) <- names(treat.list)
        }
        else names(out[["Time.Balance"]]) <- seq_along(treat.list)
        
        out[["Observations"]] <- lapply(out[["Time.Balance"]], function(x) x$Observations)
        
        if (!(quick && !msm.summary) && !nunique.gt(treat.type, 1) && !any(treat.type == "multinomial")) {
            out[["Balance.Across.Times"]] <- balance.table.msm.summary(out[["Time.Balance"]],
                                                                   weight.names = names(weights),
                                                                   no.adj = no.adj,
                                                                   m.threshold = m.threshold, 
                                                                   v.threshold = v.threshold, 
                                                                   ks.threshold = ks.threshold, 
                                                                   r.threshold = r.threshold, 
                                                                   quick = quick, 
                                                                   types = NULL)
        }
        
        out[["call"]] <- call
        
        out[["print.options"]] <- list(m.threshold=m.threshold, 
                                       v.threshold=v.threshold,
                                       ks.threshold=ks.threshold,
                                       r.threshold = r.threshold,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.means=disp.means, 
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks, 
                                       disp.adj=!no.adj, 
                                       disp.bal.tab = disp.bal.tab, 
                                       which.cluster=which.cluster,
                                       cluster.summary=cluster.summary,
                                       quick = quick,
                                       nweights = ifelse(no.adj, 0, ncol(weights)),
                                       weight.names = names(weights),
                                       which.time = which.time,
                                       msm.summary = msm.summary)
        
        class(out) <- c("bal.tab.msm", "bal.tab")
    }
    
    return(out)
}

#Point treatments
bal.tab.matchit <- function(m, int = FALSE, distance = NULL, addl = NULL, data = NULL,  continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    
    out <- base.bal.tab(weights=X$weights, 
                        treat=X$treat, 
                        distance=X$distance, 
                        subclass=X$subclass, 
                        covs=X$covs, 
                        call=X$call, 
                        int=int, 
                        addl=X$addl, 
                        continuous=continuous, 
                        binary=binary, 
                        s.d.denom=s.d.denom, 
                        m.threshold=m.threshold, 
                        v.threshold=v.threshold, 
                        ks.threshold=ks.threshold, 
                        imbalanced.only = imbalanced.only,
                        un=un, 
                        disp.means=disp.means, 
                        disp.v.ratio=disp.v.ratio, 
                        disp.ks=disp.ks, 
                        disp.subclass=disp.subclass,                                 
                        disp.bal.tab = disp.bal.tab,
                        method=X$method, 
                        cluster = X$cluster, 
                        which.cluster = which.cluster, 
                        cluster.summary = cluster.summary, 
                        discarded = X$discarded, 
                        quick = quick)
    return(out)
}
bal.tab.ps <- function(ps, stop.method, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
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
                   stop.method = stop.method, 
                   s.d.denom = s.d.denom, 
                   distance = distance,
                   addl = addl,
                   cluster = cluster,
                   ...)
    
    out <- base.bal.tab(weights=X$weights, 
                        treat=X$treat, 
                        distance=X$distance, 
                        covs=X$covs, 
                        call=X$call, 
                        int=int, 
                        addl=X$addl, 
                        continuous=continuous, 
                        binary=binary, 
                        s.d.denom=X$s.d.denom, 
                        m.threshold=m.threshold, 
                        v.threshold=v.threshold, 
                        ks.threshold=ks.threshold, 
                        imbalanced.only = imbalanced.only,
                        un=un, 
                        disp.means=disp.means, 
                        disp.v.ratio=disp.v.ratio, 
                        disp.ks=disp.ks, 
                        disp.bal.tab = disp.bal.tab,
                        method=X$method, 
                        cluster = X$cluster, 
                        which.cluster = which.cluster, 
                        cluster.summary = cluster.summary, 
                        s.weights = X$s.weights, 
                        quick = quick)
    return(out)
}
bal.tab.mnps <- function(mnps, stop.method, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = TRUE, quick = FALSE, ...) {
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
    X <- x2base.mnps(mnps,
                     stop.method = stop.method,
                     s.d.denom = s.d.denom,
                     distance = distance,
                     addl = addl,
                     cluster = cluster,
                     ...)
    
    #stop("bal.tab is not yet compatible with mnps objects.", call. = FALSE)
    out <- base.bal.tab.multi(weights=X$weights, 
                              treat=X$treat, 
                              distance=X$distance, 
                              covs=X$covs, 
                              call=X$call, 
                              int=int, 
                              addl=X$addl, 
                              continuous=continuous, 
                              binary=binary, 
                              s.d.denom=X$s.d.denom, 
                              m.threshold=m.threshold, 
                              v.threshold=v.threshold, 
                              ks.threshold=ks.threshold, 
                              imbalanced.only = imbalanced.only,
                              un=un, 
                              disp.means=disp.means, 
                              disp.v.ratio=disp.v.ratio, 
                              disp.ks=disp.ks, 
                              disp.bal.tab = disp.bal.tab,
                              method="weighting", 
                              cluster = X$cluster, 
                              which.cluster = which.cluster, 
                              cluster.summary = cluster.summary, 
                              pairwise = pairwise, focal = X$focal,
                              which.treat = which.treat,
                              multi.summary = multi.summary,
                              s.weights = X$s.weights, 
                              quick = quick)
    return(out)
}
bal.tab.Match <- function(M, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    out <- base.bal.tab(weights=X$weights, 
                        treat=X$treat, 
                        distance=X$distance, 
                        covs=X$covs, 
                        call=X$call, 
                        int=int, 
                        addl=X$addl, 
                        continuous=continuous, 
                        binary=binary, 
                        s.d.denom=X$s.d.denom, 
                        m.threshold=m.threshold, 
                        v.threshold=v.threshold, 
                        ks.threshold=ks.threshold, 
                        imbalanced.only = imbalanced.only,
                        un=un, 
                        disp.means=disp.means, 
                        disp.v.ratio=disp.v.ratio, 
                        disp.ks=disp.ks, 
                        disp.bal.tab = disp.bal.tab,
                        method=X$method, 
                        cluster = X$cluster, 
                        which.cluster = which.cluster, 
                        cluster.summary = cluster.summary, 
                        discarded = X$discarded, 
                        quick = quick)
    return(out)
}
bal.tab.formula <- function(formula, data, ...) {
    
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
    
    #Checks
    if (length(data) == 0) {
        stop("Data must be specified.", call. = FALSE)}
    if (!is.data.frame(data)) {
        stop("Data must be a data.frame.", call. = FALSE)}
    
    #Initializing variables
    tt <- terms(formula)
    attr(tt, "intercept") <- 0
    if (is.na(match(rownames(attr(tt, "factors"))[1], names(data)))) {
        stop(paste0("The given response variable, \"", rownames(attr(tt, "factors"))[1], "\", is not a variable in data."), call. = FALSE)
    }
    mf <- tryCatch(model.frame(tt, data), error = function(e) {
        stop(paste0(c("All variables of formula must be variables in data.\nVariables not in data: ",
                      paste(attr(tt, "term.labels")[is.na(match(attr(tt, "term.labels"), names(data)))], collapse=", "))), call. = FALSE)}
    )
    
    treat <- model.response(mf)
    covs <- data[!is.na(match(names(data), attr(tt, "term.labels")))]
    
    out <- do.call(bal.tab.data.frame, c(list(covs = covs, treat = treat), args))
    
    return(out)
}
bal.tab.data.frame <- function(covs, treat, data = NULL, weights = NULL, distance = NULL, subclass = NULL, match.strata = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = TRUE, s.weights = NULL, estimand = NULL, quick = FALSE, ...) {
    
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
                           imp = imp,
                           s.weights = s.weights,
                           focal = focal)
    
    if (length(X$imp) > 0) {
        out <- base.bal.tab.imp(weights=X$weights, 
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
                                v.threshold=v.threshold, ks.threshold=ks.threshold, 
                                r.threshold=r.threshold,
                                imbalanced.only = imbalanced.only,
                                un=un, 
                                disp.means=disp.means, 
                                disp.v.ratio=disp.v.ratio, disp.ks=disp.ks, 
                                disp.subclass=disp.subclass, 
                                disp.bal.tab = disp.bal.tab,
                                method=X$method, 
                                cluster = X$cluster, 
                                which.cluster = which.cluster, 
                                cluster.summary = cluster.summary, 
                                imp = X$imp,
                                which.imp = which.imp,
                                imp.summary = imp.summary,
                                s.weights = X$s.weights,
                                quick = quick,
                                ...)
    }
    else if (length(unique(X$treat)) > 2 && is.numeric(X$treat)) {
        out <- base.bal.tab.cont(weights=X$weights, 
                                 treat = X$treat, 
                                 distance = X$distance, 
                                 subclass = X$subclass,
                                 covs=X$covs, 
                                 call=X$call, 
                                 int=int, 
                                 addl = X$addl, 
                                 r.threshold = r.threshold, 
                                 imbalanced.only = imbalanced.only,
                                 un = un, 
                                 disp.subclass=disp.subclass, 
                                 disp.bal.tab = disp.bal.tab,
                                 method=X$method, 
                                 cluster = X$cluster, 
                                 which.cluster = which.cluster, 
                                 cluster.summary = cluster.summary, 
                                 s.weights = X$s.weights, 
                                 quick = quick,
                                 ...)
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        out <- base.bal.tab.multi(weights=X$weights, 
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
                                  ks.threshold=ks.threshold, 
                                  imbalanced.only = imbalanced.only,
                                  un=un, 
                                  disp.means=disp.means, 
                                  disp.v.ratio=disp.v.ratio, 
                                  disp.ks=disp.ks, 
                                  disp.subclass=disp.subclass, 
                                  disp.bal.tab = disp.bal.tab,
                                  method=X$method, 
                                  cluster = X$cluster, 
                                  which.cluster = which.cluster, 
                                  cluster.summary = cluster.summary, 
                                  pairwise = pairwise, focal = focal,
                                  which.treat = which.treat,
                                  multi.summary = multi.summary,
                                  s.weights = X$s.weights,
                                  quick = quick,
                                  ...)
    }
    else {
        out <- base.bal.tab(weights=X$weights, 
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
                            ks.threshold=ks.threshold, 
                            imbalanced.only = imbalanced.only,
                            un=un, 
                            disp.means=disp.means, 
                            disp.v.ratio=disp.v.ratio, 
                            disp.ks=disp.ks, 
                            disp.subclass=disp.subclass, 
                            disp.bal.tab = disp.bal.tab,
                            method=X$method, 
                            cluster = X$cluster, 
                            which.cluster = which.cluster, 
                            cluster.summary = cluster.summary, 
                            s.weights = X$s.weights,
                            quick = quick,
                            ...)
    }
    return(out)
}
bal.tab.CBPS <- function(cbps, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, s.weights = NULL, quick = FALSE, ...) {
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
    
    if (any(class(cbps) == "CBMSM")) {
        if (length(cluster) > 0) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
        
        #Initializing variables
        X <- x2base.CBMSM(cbps, 
                          s.d.denom = s.d.denom, 
                          distance.list = distance,
                          addl.list = addl,
                          cluster = cluster,
                          s.weights = s.weights,
                          ...)
        out <- base.bal.tab.msm(weights=X$weights, 
                                treat.list=X$treat.list, 
                                distance.list=X$distance.list, 
                                covs.list=X$covs.list, 
                                call=X$call, 
                                int=int, 
                                addl.list=X$addl.list, 
                                continuous=continuous, 
                                binary=binary, 
                                s.d.denom=X$s.d.denom, 
                                m.threshold=m.threshold, 
                                v.threshold=v.threshold, 
                                ks.threshold=ks.threshold, 
                                r.threshold = r.threshold,
                                imbalanced.only = imbalanced.only,
                                un=un, 
                                disp.means=disp.means, 
                                disp.v.ratio=disp.v.ratio, 
                                disp.ks=disp.ks, 
                                disp.bal.tab = disp.bal.tab,
                                method=X$method, 
                                pairwise = pairwise, 
                                which.treat = which.treat, 
                                multi.summary = multi.summary, 
                                s.weights = X$s.weights, 
                                which.time = which.time, 
                                msm.summary = msm.summary,
                                quick = quick, 
                                ...)
    }
    else {
        #Initializing variables
        X <- x2base.CBPS(cbps, 
                         s.d.denom = s.d.denom,
                         distance = distance,
                         addl = addl,
                         cluster = cluster, 
                         use.weights = args$use.weights,
                         s.weights = s.weights)
        if (any(class(cbps) == "CBPSContinuous")) {
            out <- base.bal.tab.cont(weights=X$weights, 
                                     treat = X$treat, 
                                     distance = X$distance, 
                                     covs=X$covs, 
                                     call=X$call, 
                                     int=int, 
                                     addl = X$addl, 
                                     r.threshold = r.threshold, 
                                     un = un, 
                                     disp.bal.tab = disp.bal.tab,
                                     method = "weighting", 
                                     cluster = X$cluster, 
                                     which.cluster = which.cluster, 
                                     cluster.summary = cluster.summary, 
                                     s.weights = X$s.weights,
                                     quick = quick)
        }
        else if (length(unique(X$treat)) > 2) {
            out <- base.bal.tab.multi(weights=X$weights, 
                                      treat=X$treat, 
                                      distance=X$distance, 
                                      covs=X$covs, 
                                      call=X$call, 
                                      int=int, 
                                      addl=X$addl, 
                                      continuous=continuous, 
                                      binary=binary, 
                                      s.d.denom=X$s.d.denom, 
                                      m.threshold=m.threshold, 
                                      v.threshold=v.threshold, 
                                      ks.threshold=ks.threshold, 
                                      imbalanced.only = imbalanced.only,
                                      un=un, 
                                      disp.means=disp.means, 
                                      disp.v.ratio=disp.v.ratio, 
                                      disp.ks=disp.ks, 
                                      disp.bal.tab = disp.bal.tab,
                                      method="weighting", 
                                      cluster = X$cluster, 
                                      which.cluster = which.cluster, 
                                      cluster.summary = cluster.summary, 
                                      pairwise = pairwise, 
                                      focal = X$focal, #NULL
                                      which.treat = which.treat,
                                      multi.summary = multi.summary,
                                      s.weights = X$s.weights,
                                      quick = quick)
        }
        else out <- base.bal.tab(weights=X$weights, 
                                 treat=X$treat, 
                                 distance=X$distance, 
                                 covs=X$covs, 
                                 call=X$call, 
                                 int=int, 
                                 addl=X$addl, 
                                 continuous=continuous, 
                                 binary=binary, 
                                 s.d.denom=X$s.d.denom, 
                                 m.threshold=m.threshold, 
                                 v.threshold=v.threshold, 
                                 ks.threshold=ks.threshold, 
                                 imbalanced.only = imbalanced.only,
                                 un=un, 
                                 disp.means=disp.means, 
                                 disp.v.ratio=disp.v.ratio, 
                                 disp.ks=disp.ks, 
                                 disp.bal.tab = disp.bal.tab,
                                 method="weighting", 
                                 cluster = X$cluster, 
                                 which.cluster = which.cluster, 
                                 cluster.summary = cluster.summary,
                                 s.weights = s.weights, 
                                 quick = quick)
    }
    return(out)
}
bal.tab.ebalance <- function(ebal, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    X <- x2base(ebal, 
                formula = formula,
                data = data, 
                treat = treat,
                covs = covs,
                distance = distance,
                addl = addl,
                cluster = cluster)
    out <- base.bal.tab(weights=X$weights, 
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
                        ks.threshold=ks.threshold, 
                        imbalanced.only = imbalanced.only,
                        un=un, 
                        disp.means=disp.means, 
                        disp.v.ratio=disp.v.ratio, 
                        disp.ks=disp.ks, 
                        disp.bal.tab = disp.bal.tab,
                        method=X$method, 
                        cluster = X$cluster, 
                        which.cluster = which.cluster, 
                        cluster.summary = cluster.summary, 
                        quick = quick)
    return(out)
}
bal.tab.ebalance.trim <- bal.tab.ebalance
bal.tab.optmatch <- function(optmatch, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    out <- base.bal.tab(weights=X$weights, 
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
                        v.threshold=v.threshold, ks.threshold=ks.threshold, 
                        imbalanced.only = imbalanced.only,
                        un=un, 
                        disp.means=disp.means, 
                        disp.v.ratio=disp.v.ratio, 
                        disp.ks=disp.ks, 
                        disp.bal.tab = disp.bal.tab,
                        method=X$method, 
                        cluster = X$cluster, 
                        which.cluster = which.cluster, 
                        cluster.summary = cluster.summary, 
                        quick = quick)
    return(out)
}
bal.tab.weightit <- function(weightit, int = FALSE, distance = NULL, addl = NULL, data = NULL,  continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, which.treat = NA, pairwise = TRUE, focal = NULL, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, quick = FALSE, ... ) {
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
    
    if (any(class(weightit) == "weightitMSM")) {
        if (length(cluster) > 0) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
        if (length(imp) > 0) stop("Multiply imputed data is not yet supported with longitudinal treatments.", call. = FALSE)
        if (length(args$addl.list) > 0) addl <- args$addl.list
        
        #Initializing variables
        X <- x2base.weightitMSM(weightit, 
                                s.d.denom = s.d.denom, 
                                distance.list = distance,
                                addl.list = addl,
                                cluster = cluster,
                                data = data,
                                ...)
        out <- base.bal.tab.msm(weights=X$weights, 
                                treat.list=X$treat.list, 
                                distance.list=X$distance.list, 
                                covs.list=X$covs.list, 
                                call=X$call, 
                                int=int, 
                                addl.list=X$addl.list, 
                                continuous=continuous, 
                                binary=binary, 
                                s.d.denom=X$s.d.denom, 
                                m.threshold=m.threshold, 
                                v.threshold=v.threshold, 
                                ks.threshold=ks.threshold, 
                                r.threshold=r.threshold,
                                imbalanced.only = imbalanced.only,
                                un=un, 
                                disp.means=disp.means, 
                                disp.v.ratio=disp.v.ratio, 
                                disp.ks=disp.ks, 
                                disp.bal.tab = disp.bal.tab,
                                method=X$method, 
                                cluster = X$cluster, 
                                which.cluster = which.cluster, 
                                cluster.summary = cluster.summary, 
                                pairwise = pairwise, 
                                focal = focal, 
                                which.treat = which.treat, 
                                multi.summary = multi.summary, 
                                s.weights = X$s.weights, 
                                which.time = which.time, 
                                msm.summary = msm.summary,
                                quick = quick, 
                                ...)
    }
    else {
        #Initializing variables
        X <- x2base.weightit(weightit, data = data, distance = distance, addl = addl, cluster = cluster, imp = imp)
        
        if (length(X$imp) > 0) {
            out <- base.bal.tab.imp(weights=X$weights, 
                                    treat=X$treat, 
                                    distance=X$distance, 
                                    subclass=NULL, 
                                    covs=X$covs, 
                                    call=X$call, 
                                    int=int, 
                                    addl=X$addl, 
                                    continuous=continuous, 
                                    binary=binary, 
                                    s.d.denom=X$s.d.denom, 
                                    m.threshold=m.threshold, 
                                    v.threshold=v.threshold, 
                                    ks.threshold=ks.threshold, 
                                    r.threshold=r.threshold,
                                    imbalanced.only = imbalanced.only,
                                    un=un, 
                                    disp.means=disp.means, 
                                    disp.v.ratio=disp.v.ratio, 
                                    disp.ks=disp.ks, 
                                    disp.subclass=FALSE, 
                                    disp.bal.tab = disp.bal.tab,
                                    method=X$method, 
                                    cluster = X$cluster, 
                                    which.cluster = which.cluster, 
                                    cluster.summary = cluster.summary, 
                                    imp = X$imp,
                                    which.imp = which.imp,
                                    imp.summary = imp.summary,
                                    s.weights = X$s.weights,
                                    discarded = X$discarded, 
                                    quick = quick,
                                    ...)
        }
        else if (weightit$treat.type == "binary") {
            out <- base.bal.tab(weights=X$weights, 
                                treat=X$treat, 
                                distance=X$distance, 
                                subclass=NULL, 
                                covs=X$covs, 
                                call=X$call, 
                                int=int, 
                                addl=X$addl, 
                                continuous=continuous, 
                                binary=binary, 
                                s.d.denom=X$s.d.denom, 
                                m.threshold=m.threshold, 
                                v.threshold=v.threshold, 
                                ks.threshold=ks.threshold, 
                                imbalanced.only = imbalanced.only,
                                un=un, 
                                disp.means=disp.means, 
                                disp.v.ratio=disp.v.ratio, 
                                disp.ks=disp.ks, 
                                disp.subclass=FALSE, 
                                disp.bal.tab = disp.bal.tab,
                                method = X$method, 
                                cluster = X$cluster, 
                                which.cluster = which.cluster, 
                                cluster.summary = cluster.summary, 
                                s.weights = X$s.weights, 
                                discarded = X$discarded, 
                                quick = quick)
        }
        else if (weightit$treat.type == "multinomial") {
            out <- base.bal.tab.multi(weights=X$weights, treat=X$treat, distance=X$distance, 
                                      covs=X$covs, call=X$call, int=int, addl=X$addl, 
                                      continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, 
                                      m.threshold=m.threshold, v.threshold=v.threshold, 
                                      ks.threshold=ks.threshold, 
                                      imbalanced.only = imbalanced.only,
                                      un=un, 
                                      disp.bal.tab = disp.bal.tab, 
                                      disp.means=disp.means,
                                      disp.v.ratio=disp.v.ratio, 
                                      disp.ks=disp.ks, 
                                      method=X$method, 
                                      cluster = X$cluster, which.cluster = which.cluster, 
                                      cluster.summary = cluster.summary, pairwise = pairwise, focal = X$focal,
                                      which.treat = which.treat, multi.summary = multi.summary,
                                      s.weights = X$s.weights, quick = quick)
        }
        else if (weightit$treat.type == "continuous") {
            out <- base.bal.tab.cont(weights=X$weights, treat = X$treat, distance = X$distance, 
                                     covs=X$covs, call=X$call, int=int, addl = X$addl, 
                                     r.threshold = r.threshold, 
                                     imbalanced.only = imbalanced.only,
                                     un=un, 
                                     disp.bal.tab = disp.bal.tab,
                                     method = X$method, 
                                     cluster = X$cluster, which.cluster = which.cluster, 
                                     cluster.summary = cluster.summary, s.weights = X$s.weights, 
                                     quick = quick)
        }
    }
    
    return(out)
}

#MSMs wth multiple time points
bal.tab.list <- function(list_, data, treat.list = NULL, weights = NULL, int = FALSE, distance.list = NULL, addl.list = NULL, method, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, pairwise = TRUE, which.treat = NA, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, s.weights = NULL, estimand = "ATE", quick = FALSE, ...) {
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
    
    if (length(args$cluster) > 0) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
    if (length(args$imp) > 0) stop("Multiply imputed data is not yet supported with longitudinal treatments.", call. = FALSE)
    
    if (all(sapply(list_, is.formula))) {
        X <- X2base.formula.list(formula.list = list_, 
                                 data = data,
                                 weights = weights,
                                 distance.list = distance.list,
                                 addl.list = addl.list,
                                 method = method,
                                 s.d.denom = s.d.denom,
                                 s.weights = s.weights,
                                 estimand = estimand)
    }
    else if (all(sapply(list_, is.data.frame))) {
        X <- X2base.data.frame.list(covs.list = list_, 
                                    treat.list = treat.list,
                                    data = data,
                                    weights = weights,
                                    distance.list = distance.list,
                                    addl.list = addl.list,
                                    method = method,
                                    s.d.denom = s.d.denom,
                                    s.weights = s.weights,
                                    estimand = estimand)
    }
    else {
        stop("If the first argument is a list, it must be a list of formulas specifying the treatment/covariate relationships at each time point or a list of data frames containing covariates to be assessed at each time point.", call. = FALSE)
    }
    
    out <- base.bal.tab.msm(weights=X$weights, 
                            treat.list=X$treat.list, 
                            distance.list=X$distance.list, 
                            covs.list=X$covs.list, 
                            call=X$call, 
                            int=int, 
                            addl.list=X$addl.list, 
                            continuous=continuous, 
                            binary=binary, 
                            s.d.denom=X$s.d.denom, 
                            m.threshold=m.threshold, 
                            v.threshold=v.threshold, 
                            ks.threshold=ks.threshold, 
                            r.threshold = r.threshold,
                            imbalanced.only = imbalanced.only,
                            un=un, 
                            disp.means=disp.means, 
                            disp.v.ratio=disp.v.ratio, 
                            disp.ks=disp.ks, 
                            disp.bal.tab = disp.bal.tab,
                            method=X$method, 
                            cluster = NULL, 
                            which.cluster = NULL, 
                            cluster.summary = FALSE, 
                            pairwise = pairwise, 
                            focal = NULL, 
                            which.treat = which.treat, 
                            multi.summary = multi.summary, 
                            s.weights = X$s.weights, 
                            which.time = which.time, 
                            msm.summary = msm.summary,
                            quick = quick, 
                            ...)
    return(out)
}
bal.tab.iptw <- function(iptw, stop.method, int = FALSE, distance.list = NULL, addl.list = NULL, data = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = FALSE, un = FALSE, disp.bal.tab = TRUE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, pairwise = TRUE, which.treat = NA, multi.summary = TRUE, which.time = NULL, msm.summary = TRUE, quick = FALSE, ...) {
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
    
    if (length(args$cluster) > 0) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
    
    #Initializing variables
    X <- x2base.iptw(iptw, 
                     stop.method = stop.method, 
                     s.d.denom = s.d.denom, 
                     distance.list = distance.list,
                     addl.list = addl.list,
                     ...)
    
    out <- base.bal.tab.msm(weights=X$weights, 
                            treat.list=X$treat.list, 
                            distance.list=X$distance.list, 
                            covs.list=X$covs.list, 
                            call=X$call, 
                            int=int, 
                            addl.list=X$addl.list, 
                            continuous=continuous, 
                            binary=binary, 
                            s.d.denom=X$s.d.denom, 
                            m.threshold=m.threshold, 
                            v.threshold=v.threshold, 
                            ks.threshold=ks.threshold, 
                            imbalanced.only = imbalanced.only,
                            un=un, 
                            disp.means=disp.means, 
                            disp.v.ratio=disp.v.ratio, 
                            disp.ks=disp.ks, 
                            disp.bal.tab = disp.bal.tab,
                            method=X$method, 
                            pairwise = pairwise, 
                            #focal = focal, 
                            which.treat = which.treat, 
                            multi.summary = multi.summary, 
                            s.weights = X$s.weights, 
                            which.time = which.time, 
                            msm.summary = msm.summary,
                            quick = quick, 
                            ...)
    return(out)
}
bal.tab.CBMSM <- bal.tab.CBPS
