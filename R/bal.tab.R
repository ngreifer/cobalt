#Make sure new scheme doesn't inerfere with subclassification
#get bal.plot to accommodate multiple weights (limit input to 1 or facet)
#Make sure imp works with multiple weights (they should)

bal.tab <- function(...) UseMethod("bal.tab")
base.bal.tab <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, s.weights = NULL, discarded = NULL, quick = FALSE, ...) {
 
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
        
        if (length(method) == 1 && method == "subclassification") {
            stop("Subclassification with clusters is not yet supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab")
        }
        else {
            out[["Cluster.Balance"]] <- setNames(lapply(levels(cluster), function(c) setNames(list(balance.table(C = C.list[[c]], weights = weights[cluster == c, , drop = FALSE], treat = treat[cluster == c], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, s.weights = s.weights[cluster == c], no.adj = no.adj, types = types, quick = quick),
                                                                                          samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, cluster = cluster, which.cluster = c, discarded = discarded)), 
                                                                                     c("Balance.Table", "Observations"))),
                                                 levels(cluster))
            balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance.Table"]])
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
                                           un=un, 
                                           disp.means=disp.means, 
                                           disp.v.ratio=disp.v.ratio, 
                                           disp.ks=disp.ks, 
                                           disp.adj=!no.adj, 
                                           disp.subclass=disp.subclass,
                                           which.cluster=which.cluster,
                                           cluster.summary=cluster.summary,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights))
            class(out) <- c("bal.tab.cluster", "bal.tab")
        }
        
    }
    else {
        if (length(method) == 1 && method == "subclassification") {
            if (length(subclass)>0) {
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
                out[["Subclass.Balance"]] <- balance.table.subclass(C, weights=weights[,1], treat=treat, subclass=subclass, continuous=continuous, binary=binary, s.d.denom=s.d.denom[1], m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold = ks.threshold, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, quick = quick)
                out[["Subclass.Observations"]] <- samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
                out[["Balance.Across.Subclass"]] <- balance.table.across.subclass(balance.table=balance.table(C, weights[,1], treat, continuous, binary, s.d.denom[1], m.threshold, v.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, no.adj = TRUE, quick = quick), 
                                                                                  balance.table.subclass.list=out[["Subclass.Balance"]], 
                                                                                  subclass.obs=out[["Subclass.Observations"]], 
                                                                                  sub.by=sub.by, 
                                                                                  m.threshold=m.threshold, 
                                                                                  v.threshold=v.threshold, 
                                                                                  ks.threshold=ks.threshold)
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
                    out[["Balanced.KS.Subclass"]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][,"KS.Threshold"])))
                    names(out[["Balanced.KS.Subclass"]]) <- paste("Subclass", levels(subclass))
                    miks.list <- lapply(levels(subclass), function(x) {
                        mi <- max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][,"Type"]!="Distance", ], "KS.Adj", "KS.Threshold")
                        return(data.frame(Variable = row.names(mi), mi))
                    } )      
                    miks <- do.call("rbind", miks.list)
                    
                    out[["Max.Imbalance.KS.Subclass"]] <- data.frame(miks, row.names = paste("Subclass", levels(subclass)))
                }
                
                if (length(call) > 0) out[["call"]] <- call
                out[["print.options"]] <- list(m.threshold=m.threshold, 
                                               v.threshold=v.threshold, 
                                               ks.threshold=ks.threshold, 
                                               un=un, 
                                               disp.means=disp.means, 
                                               disp.v.ratio=disp.v.ratio, 
                                               disp.ks=disp.ks, 
                                               disp.adj=!no.adj, 
                                               disp.subclass=disp.subclass,
                                               quick = quick)
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
            
            C <- get.C(covs = covs, int = int, addl = addl, distance = distance, cluster = cluster)
            
            out[["Balance"]] <- balance.table(C, weights, treat, continuous, binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, s.weights = s.weights, no.adj = no.adj, quick = quick)
            
            #Ensure comaptible with multiple weights
            if (length(m.threshold) > 0) {
                if (no.adj) {
                    out[["Balanced.Means"]] <- baltal(out[["Balance"]][,"M.Threshold.Un"])
                    out[["Max.Imbalance.Means"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], "Diff.Un", "M.Threshold.Un")
                }
                else if (ncol(weights) == 1) {
                    out[["Balanced.Means"]] <- baltal(out[["Balance"]][,"M.Threshold"])
                    out[["Max.Imbalance.Means"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], "Diff.Adj", "M.Threshold")
                }
                else if (ncol(weights) > 1) {
                    out[["Balanced.Means"]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][,paste0("M.Threshold.", x)]))),
                                                        names(weights))
                    out[["Max.Imbalance.Means"]] <- cbind(Weights = names(weights),
                                                          do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], paste0("Diff.", x), paste0("M.Threshold.", x)),
                                                                                                c("Diff", "M.Threshold")))),
                                                          stringsAsFactors = FALSE)
                }
            }
            if (!attr(out[["Balance"]], "disp")["v"]) {disp.v.ratio <- FALSE; v.threshold <- NULL}
            if (length(v.threshold) > 0) {
                if (no.adj) {
                    out[["Balanced.Variances"]] <- baltal(out[["Balance"]][,"V.Threshold.Un"])
                    out[["Max.Imbalance.Variances"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], "V.Ratio.Un", "V.Threshold.Un")
                    
                }
                else if (ncol(weights) == 1) {
                    out[["Balanced.Variances"]] <- baltal(out[["Balance"]][,"V.Threshold"])
                    out[["Max.Imbalance.Variances"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], "V.Ratio.Adj", "V.Threshold")
                }
                else {
                    out[["Balanced.Variances"]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][,paste0("V.Threshold.", x)]))),
                                                        names(weights))
                    out[["Max.Imbalance.Variances"]] <- cbind(Weights = names(weights),
                                                          do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], paste0("V.Ratio.", x), paste0("V.Threshold.", x)),
                                                                                                                       c("V.Ratio", "V.Threshold")))),
                                                          stringsAsFactors = FALSE)
                }
            }
            if (!attr(out[["Balance"]], "disp")["ks"]) {disp.ks <- FALSE; ks.threshold <- NULL}
            if (length(ks.threshold) > 0) {
                if (no.adj) {
                    out[["Balanced.KS"]] <- baltal(out[["Balance"]][,"KS.Threshold.Un"])
                    out[["Max.Imbalance.KS"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], "KS.Un", "KS.Threshold.Un")
                    
                }
                else if (ncol(weights) == 1) {
                    out[["Balanced.KS"]] <- baltal(out[["Balance"]][,"KS.Threshold"])
                    out[["Max.Imbalance.KS"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], "KS.Adj", "KS.Threshold")
                }
                else {
                    out[["Balanced.KS"]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][,paste0("KS.Threshold.", x)]))),
                                                            names(weights))
                    out[["Max.Imbalance.KS"]] <- cbind(Weights = names(weights),
                                                              do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], paste0("KS.", x), paste0("KS.Threshold.", x)),
                                                                                                                           c("KS", "KS.Threshold")))),
                                                              stringsAsFactors = FALSE)
                }
            }
            out[["Observations"]] <- samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
            if (length(call) > 0) out[["call"]] <- call
            out[["print.options"]] <- list(m.threshold=m.threshold, 
                                           v.threshold=v.threshold, 
                                           ks.threshold=ks.threshold, 
                                           un=un, 
                                           disp.means=disp.means, 
                                           disp.v.ratio=disp.v.ratio, 
                                           disp.ks=disp.ks, 
                                           disp.adj=!no.adj,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights))
            class(out) <- "bal.tab"
        }
    }
    
    #attr(out, "int") <- int
    return(out)
}
base.bal.tab.cont <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, r.threshold = NULL, un = FALSE, disp.subclass = FALSE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, s.weights = NULL, discarded = NULL, quick = FALSE, ...) {
    
    #Preparations
    if (length(r.threshold) > 0) r.threshold <- abs(r.threshold)
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
                                                                                     c("Balance.Table", "Observations")))
            names(out[["Cluster.Balance"]]) <- levels(cluster)
            
            balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance.Table"]])
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
                                           un=un, 
                                           disp.adj=!no.adj, 
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
                
                out[["Subclass.Balance"]] <- balance.table.subclass.cont(C, weights=weights[,1], treat=treat, subclass=subclass, r.threshold=r.threshold, quick = quick)
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
                                               un=un,
                                               disp.adj=!no.adj, 
                                               disp.subclass=disp.subclass,
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
                    out[["Balanced.Corr"]] <- baltal(out[["Balance"]][,"R.Threshold"])
                    out[["Max.Imbalance.Corr"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], "Corr.Un", "R.Threshold")
                }
                else if (ncol(weights) == 1) {
                    out[["Balanced.Corr"]] <- baltal(out[["Balance"]][,"R.Threshold"])
                    out[["Max.Imbalance.Corr"]] <- max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], "Corr.Adj", "R.Threshold")
                }
                else if (ncol(weights) > 1) {
                    out[["Balanced.Corr"]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][,paste0("R.Threshold.", x)]))),
                                                        names(weights))
                    out[["Max.Imbalance.Corr"]] <- cbind(Weights = names(weights),
                                                          do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][,"Type"]!="Distance", ], paste0("Corr.", x), paste0("R.Threshold.", x)),
                                                                                                                       c("Corr", "R.Threshold")))),
                                                          stringsAsFactors = FALSE)
                }
            }
            if (!any(is.finite(out[["Balance"]][,"Corr.Un"]))) {r.threshold <- NULL}
            out[["Observations"]] <- samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
            if (length(call) > 0) out[["call"]] <- call
            out[["print.options"]] <- list(r.threshold=r.threshold, 
                                           un=un, 
                                           disp.adj=!no.adj,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights))
            class(out) <- c("bal.tab.cont", "bal.tab")
        }
    }
    
    attr(out, "int") <- int
    return(out)
}
base.bal.tab.imp <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, method, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, s.weights = NULL, discarded = NULL, quick = FALSE, ...) {
    
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
    if (length(r.threshold) > 0) r.threshold <- abs(r.threshold)
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
        out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) base.bal.tab.cont(weights = weights[imp==i, , drop  = FALSE], treat = treat[imp==i], distance = distance[imp==i, , drop = FALSE], subclass = subclass[imp==i], covs = covs[imp == i, , drop = FALSE], call = call, int = int, addl = addl[imp = i, , drop = FALSE], r.threshold = r.threshold, un = un, method = method, cluster = cluster[imp==i], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[imp==i], discarded = discarded[imp==i], quick = quick, ...))
    }
    else if (length(unique(treat)) > 2 && (is.factor(treat) || is.character(treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
    }
    else {#if binary treatment
        out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) base.bal.tab(weights = weights[imp==i, , drop = FALSE], treat = treat[imp==i], distance = distance[imp==i, , drop = FALSE], subclass = subclass[imp==i], covs = covs[imp==i, , drop = FALSE], call = call, int = int, addl = addl[imp==i, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, un = un, disp.means = disp.means, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, method = method, cluster = cluster[imp==i], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[imp==i], discarded = discarded[imp==i], quick = quick, ...))
    }
    
    names(out[["Imputation.Balance"]]) <- levels(imp)
    
    #Create summary of lists
    
    if (!is.na(match("bal.tab.cluster", class(out[["Imputation.Balance"]][[1]])))) {
        if (!(!imp.summary && quick)) {
            out[["Cluster.Balance.Across.Imputations"]] <- lapply(levels(cluster), 
                                                                  function(c) setNames(list(balance.table.imp.summary(lapply(out[["Imputation.Balance"]], function(i) i[["Cluster.Balance"]][[c]][["Balance.Table"]]), 
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
                                   nweights = ifelse(no.adj, 0, ncol(weights)),
                                   weight.names = names(weights))
    class(out) <- unique(c(classes, sapply(out[["Imputation.Balance"]], class)))
    
    return(out)
}

bal.tab.matchit <- function(m, int = FALSE, distance = NULL, addl = NULL, data = NULL,  continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    
    out <- base.bal.tab(weights=X$weights, treat=X$treat, distance=X$distance, subclass=X$subclass, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold=ks.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.ks=disp.ks, disp.subclass=disp.subclass, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, discarded = X$discarded, quick = quick)
    return(out)
}
bal.tab.ps <- function(ps, stop.method, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
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
    
    out <- base.bal.tab(weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold=ks.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.ks=disp.ks, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = X$s.weights, quick = quick)
    return(out)
}
bal.tab.mnps <- function(mnps, stop.method, ...) {
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
                   # stop.method = stop.method, 
                   # s.d.denom = s.d.denom, 
                   # distance = distance,
                   # addl = addl,
                   # cluster = cluster,
                   ...)
    stop("bal.tab is not yet compatible with mnps objects.", call. = FALSE)
    # out <- bal.tab.multi(X, ...)
    # return(out)
}
bal.tab.Match <- function(M, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    out <- base.bal.tab(weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold=ks.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.ks=disp.ks, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, discarded = X$discarded, quick = quick)
    return(out)
}
bal.tab.formula <- function(formula, data, weights = NULL, distance = NULL, subclass = NULL, match.strata = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, s.weights = NULL, estimand = NULL, quick = FALSE, ...) {
    
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
                                un=un, 
                                disp.means=disp.means, 
                                disp.v.ratio=disp.v.ratio, disp.ks=disp.ks, 
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
        out <- base.bal.tab.cont(weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = X$addl, r.threshold = r.threshold, un = un, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
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
                            v.threshold=v.threshold, ks.threshold=ks.threshold, 
                            un=un, 
                            disp.means=disp.means, 
                            disp.v.ratio=disp.v.ratio, disp.ks=disp.ks, 
                            disp.subclass=disp.subclass, 
                            method=X$method, 
                            cluster = X$cluster, 
                            which.cluster = which.cluster, 
                            cluster.summary = cluster.summary, 
                            quick = quick)
    }
    return(out)
}
bal.tab.data.frame <- function(covs, treat, data = NULL, weights = NULL, distance = NULL, subclass = NULL, match.strata = NULL, method, int = FALSE, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, disp.subclass = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, imp = NULL, which.imp = NA, imp.summary = TRUE, s.weights = NULL, estimand = NULL, quick = FALSE, ...) {
    
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
                           s.weights = s.weights)
    
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
                                un=un, 
                                disp.means=disp.means, 
                                disp.v.ratio=disp.v.ratio, disp.ks=disp.ks, 
                                disp.subclass=disp.subclass, 
                                method=X$method, 
                                cluster = X$cluster, 
                                which.cluster = which.cluster, 
                                cluster.summary = cluster.summary, 
                                imp = X$imp,
                                which.imp = which.imp,
                                imp.summary = imp.summary,
                                s.weights = X$s.weights,
                                quick = quick)
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
                                 un = un, 
                                 disp.subclass=disp.subclass, 
                                 method=X$method, 
                                 cluster = X$cluster, 
                                 which.cluster = which.cluster, 
                                 cluster.summary = cluster.summary, 
                                 s.weights = X$s.weights, 
                                 quick = quick)
    }
    else if (length(unique(X$treat)) > 2 && (is.factor(X$treat) || is.character(X$treat))) {
        stop("Multinomial treaments are not yet supported.", call. = FALSE)
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
                            un=un, 
                            disp.means=disp.means, 
                            disp.v.ratio=disp.v.ratio, 
                            disp.ks=disp.ks, 
                            disp.subclass=disp.subclass, 
                            method=X$method, 
                            cluster = X$cluster, 
                            which.cluster = which.cluster, 
                            cluster.summary = cluster.summary, 
                            s.weights = X$s.weights,
                            quick = quick)
    }
    return(out)
}
bal.tab.CBPS <- function(cbps, int = FALSE, distance = NULL, addl = NULL, data = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
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
                     cluster = cluster, 
                     use.weights = args$use.weights)
    if (any(class(cbps) == "CBPSContinuous")) {
        out <- base.bal.tab.cont(weights=X$weights, treat = X$treat, distance = X$distance, covs=X$covs, call=X$call, int=int, addl = X$addl, r.threshold = r.threshold, un = un, method = "weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    }
    else if (nlevels(as.factor(X$treat)) > 2) {
        stop("Multinomial treatments are not yet supported.", call. = FALSE)
    }
    else out <- base.bal.tab(weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=X$s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold=ks.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.ks=disp.ks, method="weighting", cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.ebalance <- function(ebal, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
    out <- base.bal.tab(weights=X$weights, treat=X$treat, distance=X$distance, covs=X$covs, call=X$call, int=int, addl=X$addl, continuous=continuous, binary=binary, s.d.denom=s.d.denom, m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold=ks.threshold, un=un, disp.means=disp.means, disp.v.ratio=disp.v.ratio, disp.ks=disp.ks, method=X$method, cluster = X$cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, quick = quick)
    return(out)
}
bal.tab.optmatch <- function(optmatch, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, distance = NULL, addl = NULL, continuous = c("std", "raw"), binary = c("raw", "std"), s.d.denom = c("treated", "control", "pooled"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, which.cluster = NULL, cluster.summary = TRUE, quick = FALSE, ...) {
    
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
                        un=un, 
                        disp.means=disp.means, 
                        disp.v.ratio=disp.v.ratio, 
                        disp.ks=disp.ks, 
                        method=X$method, 
                        cluster = X$cluster, 
                        which.cluster = which.cluster, 
                        cluster.summary = cluster.summary, 
                        quick = quick)
    return(out)
}