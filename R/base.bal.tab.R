base.bal.tab.binary <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, addl.sds = NULL, ...) {
    #Preparations
    args <- list(...)
   
    if (nunique(treat) != 2) {
        stop("Treatment indicator must be a binary (0, 1) variable---i.e., treatment (1) or control (0)", call. = FALSE)
    }
    else if (is.factor(treat) || is.character(treat)) {
        if (is.factor(treat)) treat.names <- unique.treat <- levels(treat)
        else treat.names <- unique.treat <- unique(treat, nmax = 2)
    }
    else {
        treat.names <- c("Control", "Treated")
        unique.treat <- sort(unique(treat, nmax = 2))
    }
    names(treat.names) <- unique.treat
    
    check_if_zero_weights(weights, treat, unique.treat, treat.type = "cat")
    
    treat <- binarize(treat)
    if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
    if (is_not_null(v.threshold)) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (is_null(ks.threshold) && is_null(args$k.threshold)) {
        ks.threshold <- args$k.threshold
    }
    if (is_not_null(ks.threshold)) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (is_null(weights) && is_null(subclass)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (is_null(s.weights)) {
        s.weights <- rep(1, length(treat))
    }
    if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
    
    #Actions
    if (nunique.gt(cluster, 1)) {
        out.names <- c("Cluster.Balance", 
                       "Cluster.Balance.Across.Subclass", 
                       "Cluster.Summary", "Observations",
                       "call")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        cluster <- factor(cluster)
        
        C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, cluster = cluster, ...)
        C.list <- setNames(lapply(levels(cluster), function(x) C[cluster == x, , drop = FALSE]), 
                           levels(cluster))
        types <- get.types(C)
        co.names <- attr(C, "co.names")
        
        if (length(method) == 1 && method == "subclassification") {
            stop("Subclassification with clusters is not yet supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab")
        }
        else {
            out[["Cluster.Balance"]] <- setNames(lapply(levels(cluster), function(c) setNames(list(balance.table(C = C.list[[c]], weights = weights[cluster == c, , drop = FALSE], treat = treat[cluster == c], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, s.weights = s.weights[cluster == c], abs = abs, no.adj = no.adj, types = types, quick = quick),
                                                                                                   samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, cluster = cluster, which.cluster = c, discarded = discarded)), 
                                                                                              c("Balance", "Observations"))),
                                                 levels(cluster))
            balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance"]])
            observations <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Observations"]])
            
            if (cluster.summary || !quick) out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.tables,
                                                                                                     weight.names = names(weights),
                                                                                                     no.adj = no.adj,
                                                                                                     abs = abs,
                                                                                                     quick = quick,
                                                                                                     types = types)
            
            for (i in names(attr(balance.tables[[1]], "disp"))) {
                if (all(vapply(balance.tables, function(x) !attr(x, "disp")[i], logical(1L)))) assign(paste0("disp.", i), FALSE)
            }
            for (i in names(attr(balance.tables[[1]], "disp"))) {
                if (all(vapply(balance.tables, function(x) is_null(attr(x, "threshold")[i]), logical(1L)))) assign(paste0(i, ".threshold"), NULL)
            }
            
            # if (all(vapply(balance.tables, function(x) !attr(x, "disp")["means"], logical(1L)))) {disp.means <- FALSE}
            # if (all(vapply(balance.tables, function(x) !attr(x, "disp")["sds"], logical(1L)))) {disp.sds <- FALSE}
            # if (all(vapply(balance.tables, function(x) !attr(x, "disp")["v.ratio"], logical(1L)))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
            # if (all(vapply(balance.tables, function(x) !attr(x, "disp")["ks"], logical(1L)))) {disp.ks <- FALSE; ks.threshold <- NULL}
            
            out <- out[names(out) %nin% "Cluster.Balance.Across.Subclass"]
            out[["Observations"]] <- samplesize.across.clusters(observations)
            out[["call"]] <- call
            attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                           v.threshold=v.threshold,
                                           ks.threshold=ks.threshold,
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.means=disp.means, 
                                           disp.sds=disp.sds,
                                           disp.v.ratio=disp.v.ratio, 
                                           disp.ks=disp.ks, 
                                           disp.adj=!no.adj, 
                                           disp.subclass=disp.subclass,
                                           disp.bal.tab = disp.bal.tab, 
                                           which.cluster=which.cluster,
                                           cluster.summary=cluster.summary,
                                           cluster.fun = cluster.fun,
                                           abs = abs,
                                           continuous = continuous,
                                           binary = binary,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights),
                                           treat.names = treat.names,
                                           co.names = co.names)
            class(out) <- c("bal.tab.cluster", "bal.tab")
        }
        
    }
    else {
        if (length(method) == 1 && method == "subclassification") {
            if (is_not_null(subclass)) {
                out.names <- c("Subclass.Balance", "Balance.Across.Subclass", 
                               expand.grid_string(c("Balanced", "Max.Imbalance"),
                                                  c("Means", "Variances", "KS"),
                                                  "Subclass", collapse = "."), 
                               "Subclass.Observations", "call")
                out <- vector("list", length(out.names))
                names(out) <- out.names
                
                C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
                co.names <- attr(C, "co.names")
                
                if (is_not_null(list(...)$sub.by)) sub.by <- list(...)$sub.by
                else sub.by <- call$sub.by
                out[["Subclass.Balance"]] <- balance.table.subclass(C, weights=weights[[1]], treat=treat, subclass=subclass, continuous=continuous, binary=binary, s.d.denom=s.d.denom[1], m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold = ks.threshold, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, abs = abs, quick = quick)
                out[["Subclass.Observations"]] <- samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
                out[["Balance.Across.Subclass"]] <- balance.table.across.subclass(balance.table = balance.table(C, 
                                                                                                                weights = weights[[1]], 
                                                                                                                treat = treat, 
                                                                                                                continuous = continuous, 
                                                                                                                binary = binary, 
                                                                                                                s.d.denom = s.d.denom[1], 
                                                                                                                m.threshold = m.threshold, 
                                                                                                                v.threshold = v.threshold, 
                                                                                                                ks.threshold = ks.threshold,
                                                                                                                un = un, 
                                                                                                                disp.means = disp.means, 
                                                                                                                disp.sds = disp.sds,
                                                                                                                disp.v.ratio = disp.v.ratio, 
                                                                                                                disp.ks = disp.ks,
                                                                                                                abs = abs, 
                                                                                                                no.adj = TRUE, quick = quick), 
                                                                                  balance.table.subclass.list=out[["Subclass.Balance"]], 
                                                                                  subclass.obs=out[["Subclass.Observations"]], 
                                                                                  sub.by=sub.by, 
                                                                                  m.threshold=m.threshold, 
                                                                                  v.threshold=v.threshold, 
                                                                                  ks.threshold=ks.threshold,
                                                                                  s.d.denom = s.d.denom[1])
                #Reassign disp... and ...threshold based on balance table output
                for (i in names(attr(out[["Subclass.Balance"]], "disp"))) {
                    assign(paste.("disp", i), attr(out[["Subclass.Balance"]], "disp")[i])
                }
                for (i in names(attr(out[["Subclass.Balance"]], "threshold"))) {
                    assign(paste.(i, "threshold"), attr(out[["Subclass.Balance"]], "threshold")[i])
                }
                
                S <- list(diff = list(threshold = m.threshold,
                                      Names = "Means",
                                      Threshold = "M.Threshold",
                                      Stat = "Diff"),
                          v.ratio = list(threshold = v.threshold,
                                         Names = "Variances",
                                         Threshold = "V.Threshold",
                                         Stat = "V.Ratio"),
                          ks = list(threshold = ks.threshold,
                                    Names = "KS",
                                    Threshold = "KS.Threshold",
                                    Stat = "KS"))
                
                for (s in S) {
                    if (is_not_null(s[["threshold"]])) {
                        out[[paste.("Balanced", s[["Names"]], "Subclass")]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][[s[["Threshold"]]]])))
                        names(out[[paste.("Balanced", s[["Names"]], "Subclass")]]) <- paste("Subclass", levels(subclass))
                        max.imbal.list <- lapply(levels(subclass), function(x) {
                            return(max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][["Type"]] != "Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio"))
                        } )
                        out[[paste.("Max.Imbalance", s[["Names"]], "Subclass")]] <- data.frame(do.call("rbind", max.imbal.list), 
                                                                                              row.names = paste("Subclass", levels(subclass)))
                    }
                    else {
                        out[[paste.("Balanced", s[["Names"]], "Subclass")]] <- NULL
                        out[[paste.("Max.Imbalance", s[["Names"]], "Subclass")]] <- NULL
                    }
                }
                
                out[["call"]] <- call
                attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                               v.threshold=v.threshold, 
                                               ks.threshold=ks.threshold, 
                                               imbalanced.only = imbalanced.only,
                                               un=un, 
                                               disp.means=disp.means, 
                                               disp.sds=disp.sds,
                                               disp.v.ratio=disp.v.ratio, 
                                               disp.ks=disp.ks, 
                                               disp.adj=!no.adj, 
                                               disp.subclass=disp.subclass,
                                               disp.bal.tab = disp.bal.tab, 
                                               abs = abs,
                                               continuous = continuous,
                                               binary = binary,
                                               quick = quick,
                                               treat.names = treat.names,
                                               co.names = co.names)
                class(out) <- c("bal.tab.subclass", "bal.tab")
            }
            else stop("Method specified as subclassification, but no subclasses were specified.", call. = FALSE)
        }
        else {
            out.names <- c("Balance", 
                           expand.grid_string(c("Balanced", "Max.Imbalance"),
                                              c("Means", "Variances", "KS"),
                                              collapse = "."), 
                           "Observations", "call")
            out <- vector("list", length(out.names))
            names(out) <- out.names
            
            C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
            co.names <- attr(C, "co.names")
            
            out[["Balance"]] <- balance.table(C, weights, treat, continuous, binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, s.weights = s.weights, abs = abs, no.adj = no.adj, quick = quick, addl.sds = addl.sds)
            
            #Reassign disp... and ...threshold based on balance table output
            for (i in names(attr(out[["Balance"]], "disp"))) {
                assign(paste0("disp.", i), attr(out[["Balance"]], "disp")[i])
            }
            for (i in names(attr(out[["Balance"]], "threshold"))) {
                assign(paste0(i, ".threshold"), attr(out[["Balance"]], "threshold")[i])
            }
            
            S <- list(diff = list(threshold = m.threshold,
                                  Names = "Means",
                                  Threshold = "M.Threshold",
                                  Stat = "Diff"),
                      v.ratio = list(threshold = v.threshold,
                                     Names = "Variances",
                                     Threshold = "V.Threshold",
                                     Stat = "V.Ratio"),
                      ks = list(threshold = ks.threshold,
                                Names = "KS",
                                Threshold = "KS.Threshold",
                                Stat = "KS"))
            
            for (s in S) {
                if (is_not_null(s[["threshold"]])) {
                    if (no.adj) {
                        out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[paste.(s[["Threshold"]], "Un")]])
                        out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Un"), paste.(s[["Threshold"]], "Un"), ratio = s$Stat == "V.Ratio")
                    }
                    else if (ncol(weights) == 1) {
                        out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[s[["Threshold"]]]])
                        out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio")
                    }
                    else if (ncol(weights) > 1) {
                        out[[paste.("Balanced", s[["Names"]])]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][[paste.(s[["Threshold"]], x)]]))),
                                                                            names(weights))
                        out[[paste.("Max.Imbalance", s[["Names"]])]] <- cbind(Weights = names(weights),
                                                                              do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], x), paste.(s[["Threshold"]], x), ratio = s$Stat == "V.Ratio"),
                                                                                                                                           c("Variable", s[["Stat"]], s[["Threshold"]])))),
                                                                              stringsAsFactors = FALSE)
                    }
                }
                else {
                    out[[paste.("Balanced", s[["Names"]])]] <- NULL
                    out[[paste.("Max.Imbalance", s[["Names"]])]] <- NULL
                }
            }
            
            out[["Observations"]] <- samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded, treat.names = treat.names)
            out[["call"]] <- call
            attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                           v.threshold=v.threshold, 
                                           ks.threshold=ks.threshold, 
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.means=disp.means, 
                                           disp.sds=disp.sds,
                                           disp.v.ratio=disp.v.ratio, 
                                           disp.ks=disp.ks, 
                                           disp.adj=!no.adj,
                                           disp.bal.tab = disp.bal.tab, 
                                           abs = abs,
                                           continuous = continuous,
                                           binary = binary,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights),
                                           treat.names = treat.names,
                                           co.names = co.names)
            class(out) <- "bal.tab"
        }
    }
    
    #attr(out, "int") <- int
    return(out)
}
base.bal.tab.cont <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, r.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
    
    #Preparations
    args <- list(...)
    
    check_if_zero_weights(weights, treat.type = "cont")
    
    if (is_not_null(r.threshold)) {
        r.threshold <- abs(r.threshold)
        if (r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            r.threshold <- NULL
        }
    }
    if (is_null(weights) && is_null(subclass)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (is_null(s.weights)) {
        s.weights <- rep(1, length(treat))
    }    
    if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
    
    #Actions
    if (nlevels(cluster) > 0) {
        out.names <- c("Cluster.Balance", 
                       "Cluster.Balance.Across.Subclass", 
                       "Cluster.Summary", "Observations",
                       "call")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        cluster <- factor(cluster)
        
        out[["Cluster.Balance"]] <- vector("list", length(levels(cluster)))
        names(out[["Cluster.Balance"]]) <- levels(cluster)
        
        C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, cluster = cluster, ...)
        co.names <- attr(C, "co.names")
        C.list <- structure(lapply(levels(cluster), function(x) C[cluster == x, , drop = FALSE]), names = levels(cluster))
        types <- get.types(C)
        
        if (length(method) == 1 && method == "subclassification") {
            stop("Subclassification with clusters is not yet supported.", call. = FALSE)
            #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab") #add more for subclasses
        }
        else {
            out[["Cluster.Balance"]] <- lapply(levels(cluster), function(c) setNames(list(balance.table.cont(C = C.list[[c]], weights = weights[cluster == c, , drop = FALSE], treat = treat[cluster == c], r.threshold = r.threshold, un = un, disp.means = disp.means, disp.sds = disp.sds, s.weights = s.weights[cluster == c], abs = abs, no.adj = no.adj, types = types, quick = quick),
                                                                                          samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, cluster = cluster, which.cluster = c, discarded = discarded)), 
                                                                                     c("Balance", "Observations")))
            names(out[["Cluster.Balance"]]) <- levels(cluster)
            
            balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance"]])
            observations <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Observations"]])
            
            if (!(!cluster.summary && quick)) out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.tables,
                                                                                                        weight.names = names(weights),
                                                                                                        no.adj = no.adj,
                                                                                                        abs = abs,
                                                                                                        quick = quick,
                                                                                                        types = types)
            out <- out[names(out) %nin% "Cluster.Balance.Across.Subclass"]
            out[["Observations"]] <- samplesize.across.clusters(observations)
            
            out[["call"]] <- call
            attr(out, "print.options") <- list(r.threshold=r.threshold, 
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.means=disp.means, 
                                           disp.sds = disp.sds,
                                           disp.adj=!no.adj, 
                                           disp.bal.tab = disp.bal.tab,
                                           which.cluster=which.cluster,
                                           cluster.summary=cluster.summary,
                                           cluster.fun=cluster.fun,
                                           abs = abs,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights),
                                           co.names = co.names)
            class(out) <- c("bal.tab.cont.cluster", "bal.tab.cluster", "bal.tab.cont", "bal.tab")
        }
        
    }
    else {
        if (length(method) == 1 && method == "subclassification") {
            if (is_not_null(subclass)) {
                out.names <- c("Subclass.Balance", 
                               "Balanced.Corr.Subclass", "Max.Imbalance.Corr.Subclass", 
                               "Subclass.Observations", "call")
                out <- vector("list", length(out.names))
                names(out) <- out.names
                
                C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
                co.names <- attr(C, "co.names")
                
                # if (length(list(...)$sub.by > 0)) sub.by <- list(...)$sub.by
                # else sub.by <- call$sub.by
                
                out[["Subclass.Balance"]] <- balance.table.subclass.cont(C, weights=weights[[1]], treat=treat, subclass=subclass, r.threshold=r.threshold, disp.means = disp.means, disp.sds = disp.sds, quick = quick)
                out[["Subclass.Observations"]] <- samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
                
                #Reassign disp... and ...threshold based on balance table output
                for (i in names(attr(out[["Subclass.Balance"]], "disp"))) {
                    assign(paste.("disp", i), attr(out[["Subclass.Balance"]], "disp")[i])
                }
                for (i in names(attr(out[["Subclass.Balance"]], "threshold"))) {
                    assign(paste.(i, "threshold"), attr(out[["Subclass.Balance"]], "threshold")[i])
                }
                
                S <- list(corr = list(threshold = r.threshold,
                                      Names = "Corr",
                                      Threshold = "R.Threshold",
                                      Stat = "Corr"))
                
                for (s in S) {
                    if (is_not_null(s[["threshold"]])) {
                        out[[paste.("Balanced", s[["Stat"]], "Subclass")]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][[s[["Threshold"]]]])))
                        names(out[[paste.("Balanced", s[["Stat"]], "Subclass")]]) <- paste("Subclass", levels(subclass))
                        mi.list <- lapply(levels(subclass), function(x) {
                            return(max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][["Type"]] != "Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio"))
                        } )
                        mi <- do.call("rbind", mi.list)
                        out[[paste.("Max.Imbalance", s[["Stat"]], "Subclass")]] <- data.frame(mi, row.names = paste("Subclass", levels(subclass)))
                    }
                }
                
                out[["call"]] <- call
                attr(out, "print.options") <- list(r.threshold=r.threshold, 
                                               imbalanced.only = imbalanced.only,
                                               un=un,
                                               disp.means=disp.means, 
                                               disp.sds=disp.sds,
                                               disp.adj=!no.adj, 
                                               disp.subclass=disp.subclass,
                                               disp.bal.tab = disp.bal.tab,
                                               abs = abs,
                                               quick = quick,
                                               co.names = co.names)
                class(out) <- c("bal.tab.subclass.cont", "bal.tab.subclass", "bal.tab.cont", "bal.tab")
            }
            else stop("Method specified as subclassification, but no subclasses were specified.", call. = FALSE)
            
        }
        else {
            out.names <- c("Balance", "Balanced.Corr", 
                           "Max.Imbalance.Corr", 
                           "Observations", 
                           "call")
            out <- vector("list", length(out.names))
            names(out) <- out.names
            
            C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, cluster = cluster, ...)
            co.names <- attr(C, "co.names")
            
            out[["Balance"]] <- balance.table.cont(C, weights, treat, r.threshold, un = un, disp.means = disp.means, disp.sds = disp.sds, s.weights = s.weights, abs = abs, no.adj = no.adj, quick = quick)
            
            #Reassign disp... and ...threshold based on balance table output
            for (i in names(attr(out[["Balance"]], "disp"))) {
                assign(paste0("disp.", i), attr(out[["Balance"]], "disp")[i])
            }
            for (i in names(attr(out[["Balance"]], "threshold"))) {
                assign(paste0(i, ".threshold"), attr(out[["Balance"]], "threshold")[i])
            }
            
            S <- list(corr = list(threshold = r.threshold,
                                  Names = "Corr",
                                  Threshold = "R.Threshold",
                                  Stat = "Corr"))
            
            for (s in S) {
                if (is_not_null(s[["threshold"]])) {
                    if (no.adj) {
                        out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[paste.(s[["Threshold"]], "Un")]])
                        out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Un"), paste.(s[["Threshold"]], "Un"), ratio = s$Stat == "V.Ratio")
                    }
                    else if (ncol(weights) == 1) {
                        out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[s[["Threshold"]]]])
                        out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio")
                    }
                    else if (ncol(weights) > 1) {
                        out[[paste.("Balanced", s[["Names"]])]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][[paste.(s[["Threshold"]], x)]]))),
                                                                            names(weights))
                        out[[paste.("Max.Imbalance", s[["Names"]])]] <- cbind(Weights = names(weights),
                                                                              do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], x), paste.(s[["Threshold"]], x), ratio = s$Stat == "V.Ratio"),
                                                                                                                                           c("Variable", s[["Stat"]], s[["Threshold"]])))),
                                                                              stringsAsFactors = FALSE)
                    }
                }
            }
            
            out[["Observations"]] <- samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
            out[["call"]] <- call
            attr(out, "print.options") <- list(r.threshold=r.threshold, 
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.means=disp.means, 
                                           disp.sds = disp.sds,
                                           disp.adj=!no.adj,
                                           disp.bal.tab = disp.bal.tab,
                                           abs = abs,
                                           quick = quick,
                                           nweights = ifelse(no.adj, 0, ncol(weights)),
                                           weight.names = names(weights),
                                           co.names = co.names)
            class(out) <- c("bal.tab.cont", "bal.tab")
        }
    }
    
    attr(out, "int") <- int
    return(out)
}
base.bal.tab.imp <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), imp = NULL, which.imp = NA, imp.summary = getOption("cobalt_imp.summary", TRUE), imp.fun = getOption("cobalt_imp.fun", NULL), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
    args <- list(...)
    #Preparations
    if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
    if (is_not_null(v.threshold)) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (is_not_null(ks.threshold)) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (is_not_null(r.threshold)) {
        r.threshold <- abs(r.threshold)
        if (r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            r.threshold <- NULL
        }
    }
    if (is_null(weights)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (is_null(s.weights)) {
        s.weights <- rep(1, length(treat))
    }    
    if (is_not_null(args[["agg.fun"]])) cluster.fun <- imp.fun <- args[["agg.fun"]]
    
    #Setup output object
    out.names <- c("Imputation.Balance", 
                   "Cluster.Balance.Across.Imputations",
                   "Balance.Across.Imputations", 
                   "Observations", 
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    #Get list of bal.tabs for each imputation
    if (isTRUE(get.treat.type(treat) == "continuous") || (is.numeric(treat) && !is_binary(treat))) {#if continuous treatment
        args <- args[names(args) %nin% names(formals(base.bal.tab.cont))]
        out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) do.call("base.bal.tab.cont", 
                                                                               c(list(weights = weights[imp==i, , drop  = FALSE], 
                                                                                      treat = treat[imp==i], 
                                                                                      distance = distance[imp==i, , drop = FALSE], 
                                                                                      subclass = subclass[imp==i], 
                                                                                      covs = covs[imp == i, , drop = FALSE], 
                                                                                      call = call, 
                                                                                      int = int, 
                                                                                      poly = poly, 
                                                                                      addl = addl[imp = i, , drop = FALSE], 
                                                                                      r.threshold = r.threshold, 
                                                                                      imbalanced.only = imbalanced.only, 
                                                                                      un = un, 
                                                                                      disp.means = disp.means, 
                                                                                      disp.sds = disp.sds, 
                                                                                      disp.bal.tab = disp.bal.tab, 
                                                                                      method = method, 
                                                                                      cluster = cluster[imp==i], 
                                                                                      which.cluster = which.cluster, 
                                                                                      cluster.summary = cluster.summary, 
                                                                                      s.weights = s.weights[imp==i], 
                                                                                      discarded = discarded[imp==i], 
                                                                                      abs = abs,
                                                                                      quick = quick), 
                                                                                 args), quote = TRUE))
    }
    else if (isTRUE(get.treat.type(treat) == "multinomial") || (is_(treat, c("factor", "character")) && !is_binary(treat))) {
        args <- args[names(args) %nin% names(formals(base.bal.tab.multi))]
        stop("Multinomial treatments are not yet supported with multiply imputed data.", call. = FALSE)
    }
    else {#if binary treatment
        args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
        out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) do.call("base.bal.tab.binary", 
                                                                               c(list(weights = weights[imp==i, , drop = FALSE], 
                                                                                      treat = treat[imp==i], 
                                                                                      distance = distance[imp==i, , drop = FALSE], 
                                                                                      subclass = subclass[imp==i], 
                                                                                      covs = covs[imp==i, , drop = FALSE], 
                                                                                      call = call, 
                                                                                      int = int, 
                                                                                      poly = poly, 
                                                                                      addl = addl[imp==i, , drop = FALSE], 
                                                                                      continuous = continuous, 
                                                                                      binary = binary, 
                                                                                      s.d.denom = s.d.denom, 
                                                                                      m.threshold = m.threshold, 
                                                                                      v.threshold = v.threshold, 
                                                                                      ks.threshold = ks.threshold, 
                                                                                      imbalanced.only = imbalanced.only, 
                                                                                      un = un, 
                                                                                      disp.means = disp.means, 
                                                                                      disp.sds = disp.sds, 
                                                                                      disp.v.ratio = disp.v.ratio, 
                                                                                      disp.ks = disp.ks, 
                                                                                      disp.subclass = disp.subclass, 
                                                                                      disp.bal.tab = disp.bal.tab,
                                                                                      method = method,
                                                                                      cluster = cluster[imp==i],
                                                                                      which.cluster = which.cluster,
                                                                                      cluster.summary = cluster.summary,
                                                                                      s.weights = s.weights[imp==i],
                                                                                      discarded = discarded[imp==i],
                                                                                      abs = abs,
                                                                                      quick = quick), 
                                                                                 args), quote = TRUE))
    }
    
    names(out[["Imputation.Balance"]]) <- levels(imp)
    
    #Create summary of lists
    
    if ("bal.tab.cluster" %in% class(out[["Imputation.Balance"]][[1]])) {
        if (imp.summary || !quick) {
            out[["Cluster.Balance.Across.Imputations"]] <- lapply(levels(cluster), 
                                                                  function(c) setNames(list(balance.table.imp.summary(lapply(out[["Imputation.Balance"]], function(i) i[["Cluster.Balance"]][[c]][["Balance"]]), 
                                                                                                                      weight.names = names(weights),
                                                                                                                      no.adj = no.adj,
                                                                                                                      abs = abs, quick = quick),
                                                                                            samplesize.across.imps(lapply(out[["Imputation.Balance"]], function(i) i[["Cluster.Balance"]][[c]][["Observations"]]))), 
                                                                                       c("Cluster.Balance", "Cluster.Observations")))
            names(out[["Cluster.Balance.Across.Imputations"]]) <- levels(cluster)
            balance.tables <- lapply(out[["Cluster.Balance.Across.Imputations"]], function(c) c[["Cluster.Balance"]])
            observations <- lapply(out[["Cluster.Balance.Across.Imputations"]], function(c) c[["Cluster.Observations"]])
            
            out[["Balance.Across.Imputations"]] <- balance.table.clust.imp.summary(balance.tables,
                                                                                   weight.names = names(weights),
                                                                                   no.adj = no.adj,
                                                                                   abs = abs,
                                                                                   quick = quick,
                                                                                   types = NULL)
            out[["Observations"]] <- samplesize.across.clusters(observations)
        }
        
        classes <- c("bal.tab.imp.cluster", "bal.tab.imp")
    }
    else {
        if ("bal.tab.subclass" %in% class(out[["Imputation.Balance"]][[1]])) {
            #Put something here
            stop("Subclassification cannot be used with multiply imputed data.", call. = FALSE)
        }
        else {
            if (imp.summary || !quick) out[["Balance.Across.Imputations"]] <- balance.table.imp.summary(bal.tab.imp.list = out[["Imputation.Balance"]], 
                                                                                                        weight.names = names(weights),
                                                                                                        no.adj = no.adj,
                                                                                                        abs = abs,
                                                                                                        quick = quick,
                                                                                                        types = NULL)
            observations <- lapply(out[["Imputation.Balance"]], function(x) x[["Observations"]])
            
            out[["Observations"]] <- samplesize.across.imps(observations)
            classes <- "bal.tab.imp"
        }
    }
    
    out[["call"]] <- call
    attr(out, "print.options") <- list(m.threshold=m.threshold,
                                   v.threshold=v.threshold,
                                   ks.threshold=ks.threshold,
                                   r.threshold=r.threshold,
                                   imbalanced.only = imbalanced.only,
                                   un=un, 
                                   disp.adj=!no.adj, 
                                   which.cluster=which.cluster,
                                   cluster.summary=cluster.summary,
                                   cluster.fun = cluster.fun,
                                   which.imp=which.imp,
                                   imp.summary=imp.summary,
                                   imp.fun = imp.fun,
                                   abs = abs,
                                   continuous = continuous,
                                   binary = binary,
                                   quick = quick,
                                   disp.means=disp.means, 
                                   disp.sds = disp.sds,
                                   disp.v.ratio=disp.v.ratio, 
                                   disp.ks=disp.ks,
                                   disp.bal.tab = disp.bal.tab,
                                   nweights = ifelse(no.adj, 0, ncol(weights)),
                                   weight.names = names(weights),
                                   co.names = attr(out[["Imputation.Balance"]][[1]], "print.options")[["co.names"]])
    class(out) <- unique(c(classes, sapply(out[["Imputation.Balance"]], class)))
    
    return(out)
}
base.bal.tab.multi <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = getOption("cobalt_multi.summary", TRUE), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
    #Preparations
    args <- list(...)
    if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
    if (is_not_null(v.threshold)) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (is_not_null(ks.threshold)) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (is_null(weights) && is_null(subclass)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (is_null(s.weights)) {
        s.weights <- rep(1, length(treat))
    }
    if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
    
    #Treat is a factor variable of 3+ levels
    if (is_null(focal)) {
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
    names(treat.names) <- treat.names
    
    if (is_not_null(cluster)) {
        stop("Clusters are not yet supported with multiple categorical treatments.", call. = FALSE)
    }
    else {
        #Setup output object
        out.names <- c("Pair.Balance", 
                       "Balance.Across.Pairs", 
                       "Observations", 
                       "call")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        addl.sds <- list()
        if (any(s.d.denom %in% c("pooled", "all"))) {
            if (any(s.d.denom == "pooled")) {
                C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
                addl.sds[["pooled"]] <- sqrt(rowMeans(do.call("cbind", lapply(levels(treat), function(t) col.w.v(C[treat == t, , drop = FALSE], 
                                                                                                       s.weights[treat == t])))))
            }
            if (any(s.d.denom == "all")) {
                C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
                addl.sds[["all"]] <- sqrt(col.w.v(C, s.weights))
            }
        }
        
        if (pairwise || is_not_null(focal)) {
            args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
            balance.tables <- lapply(treat.combinations, function(t) do.call("base.bal.tab.binary", c(list(weights = weights[treat %in% t, , drop = FALSE], treat = factor(treat[treat %in% t], t), distance = distance[treat %in% t, , drop = FALSE], subclass = subclass[treat %in% t], covs = covs[treat %in% t, , drop = FALSE], call = NULL, int = int, poly = poly, addl = addl[treat %in% t, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster[treat %in% t], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[treat %in% t], discarded = discarded[treat %in% t], quick = quick, addl.sds = addl.sds), args), quote = TRUE))
        }
        else {
            if (any(treat.names == "Others")) stop ("\"Others\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
            args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
            balance.tables <- lapply(treat.combinations, function(t) {
                treat_ <- factor(treat, levels = c(levels(treat), "Others"))
                treat_[treat_ != t[1]] <- "Others"
                treat_ <- factor(treat_, rev(t))
                do.call("base.bal.tab.binary", c(list(weights = weights, treat = treat_, distance = distance, subclass = subclass, covs = covs, call = NULL, int = int, poly = poly, addl = addl, continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights, discarded = discarded, quick = quick, addl.sds = addl.sds), args), quote = TRUE)
            })
        }
        for (i in seq_along(balance.tables)) {
            names(balance.tables)[i] <- paste(rev(treat.combinations[[i]]), collapse = " vs. ")
        }
        
        out[["Pair.Balance"]] <- balance.tables
        
        out[["Observations"]] <- samplesize.multi(balance.tables, treat.names, focal)
        
        if (multi.summary || !quick) {
        out[["Balance.Across.Pairs"]] <- balance.table.multi.summary(balance.tables, 
                                                                     weight.names = names(weights),
                                                                     m.threshold = m.threshold,
                                                                     v.threshold = v.threshold,
                                                                     ks.threshold = ks.threshold,
                                                                     no.adj = no.adj,
                                                                     quick = quick,
                                                                     types = NULL)
        }
        
        out[["call"]] <- call
        
        attr(out, "print.options") <- list(m.threshold=m.threshold,
                                       v.threshold=v.threshold,
                                       ks.threshold=ks.threshold,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.adj=!no.adj, 
                                       which.cluster=which.cluster,
                                       cluster.summary=cluster.summary,
                                       cluster.fun = cluster.fun,
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       disp.means=disp.means, 
                                       disp.sds = disp.sds,
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks,
                                       disp.bal.tab = disp.bal.tab,
                                       nweights = ifelse(no.adj, 0, ncol(weights)),
                                       weight.names = names(weights),
                                       treat.names = treat.names,
                                       which.treat = which.treat,
                                       multi.summary = multi.summary,
                                       pairwise = pairwise,
                                       co.names = attr(out[["Pair.Balance"]][[1]], "print.options")[["co.names"]])
        
        class(out) <- c("bal.tab.multi", "bal.tab")
    }
    return(out)
    
}
base.bal.tab.msm <- function(weights, treat.list, distance.list = NULL, subclass = NULL, covs.list, call = NULL, int = FALSE, poly = 1, addl.list = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = getOption("cobalt_multi.summary", TRUE), which.time = NULL, msm.summary = getOption("cobalt_msm.summary", TRUE), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
    #One vector of weights
    #treat.list should be a df/list of treatment vectors, one for each time period
    #cov.list should be a list of covariate data.frames, one for each time period; 
    #   should include all covs from previous time points, but no treatment statuses
    
    #Preparations
    args <- list(...)
    
    if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
    if (is_not_null(v.threshold)) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (is_not_null(ks.threshold) && is_not_null(args$k.threshold)) {
        ks.threshold <- args$k.threshold
    }
    if (is_not_null(ks.threshold)) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (is_not_null(r.threshold)) {
        r.threshold <- abs(r.threshold)
        if (r.threshold > 1) {
            warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
            r.threshold <- NULL
        }
    }
    if (is_null(weights) && is_null(subclass)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (is_null(s.weights)) {
        s.weights <- rep(1, length(treat.list[[1]]))
    }
    if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
    
    if (nunique.gt(cluster, 1)) {
        stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
    }
    else {
        #Setup output object
        out.names <- c("Time.Balance", 
                       "Balance.Across.Times", 
                       "Observations", 
                       "call")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        out[["Time.Balance"]] <- vector("list", length(covs.list))
        
        treat.type <- vapply(treat.list, function(x) get.treat.type(assign.treat.type(x)), character(1L))
        
        #Get list of bal.tabs for each time period
        out[["Time.Balance"]] <- lapply(seq_along(treat.list), function(ti) {
            if (treat.type[ti] == "continuous") {
                args <- args[names(args) %nin% names(formals(base.bal.tab.cont))]
                out_ <- do.call("base.bal.tab.cont", c(list(weights = weights, treat = treat.list[[ti]], distance = distance.list[[ti]], subclass = NULL, covs = covs.list[[ti]], call = NULL, int = int, poly = poly, addl = addl.list[[ti]], r.threshold = r.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights, discarded = discarded, quick = quick), args),
                                quote = TRUE)
            }
            else if (treat.type[ti] == "multinomial") {
                args <- args[names(args) %nin% names(formals(base.bal.tab.multi))]
                out_ <- do.call("base.bal.tab.multi", c(list(weights = weights, treat = treat.list[[ti]], distance=distance.list[[ti]], 
                                                             covs=covs.list[[ti]], call=NULL, int=int, addl=addl.list[[ti]], 
                                                             continuous=continuous, binary=binary, s.d.denom=s.d.denom, 
                                                             m.threshold=m.threshold, v.threshold=v.threshold, 
                                                             ks.threshold=ks.threshold, 
                                                             imbalanced.only = imbalanced.only,
                                                             un=un, 
                                                             disp.bal.tab = disp.bal.tab, 
                                                             disp.means=disp.means,
                                                             disp.sds=disp.sds,
                                                             disp.v.ratio=disp.v.ratio, 
                                                             disp.ks=disp.ks, 
                                                             method=method, 
                                                             cluster = cluster, which.cluster = which.cluster, 
                                                             cluster.summary = cluster.summary, pairwise = pairwise, focal = NULL,
                                                             which.treat = which.treat, multi.summary = multi.summary,
                                                             s.weights = s.weights, quick = quick),
                                                        args),
                                quote = TRUE)
            }
            else if (treat.type[ti] == "binary") {
                args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
                out_ <- do.call("base.bal.tab.binary", c(list(weights = weights, treat = treat.list[[ti]], distance = distance.list[[ti]], subclass = NULL, covs = covs.list[[ti]], call = NULL, int = int, poly = poly, addl = addl.list[[ti]], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = FALSE, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights, discarded = discarded, quick = quick), args),
                                quote = TRUE)
            }
            else stop("Each treatment must be binary, multinomial, or continuous.", call. = FALSE)
            
            return(out_)
        })
        
        if (length(names(treat.list)) == length(treat.list)) {
            names(out[["Time.Balance"]]) <- names(treat.list)
        }
        else names(out[["Time.Balance"]]) <- seq_along(treat.list)
        
        out[["Observations"]] <- lapply(out[["Time.Balance"]], function(x) x$Observations)
        
        if (!(quick && !msm.summary) && all_the_same(treat.type) && !any(treat.type == "multinomial")) {
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
        
        attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                       v.threshold=v.threshold,
                                       ks.threshold=ks.threshold,
                                       r.threshold = r.threshold,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.means=disp.means, 
                                       disp.sds=disp.sds,
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks, 
                                       disp.adj=!no.adj, 
                                       disp.bal.tab = disp.bal.tab, 
                                       which.cluster=which.cluster,
                                       cluster.summary=cluster.summary,
                                       cluster.fun = cluster.fun,
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       nweights = ifelse(no.adj, 0, ncol(weights)),
                                       weight.names = names(weights),
                                       which.time = which.time,
                                       msm.summary = msm.summary,
                                       co.names = attr(out[["Time.Balance"]][[1]], "print.options")[["co.names"]])
        
        class(out) <- c("bal.tab.msm", "bal.tab")
    }
    
    return(out)
}
base.bal.tab.target <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), which.treat = NA, target.summary = getOption("cobalt_target.summary", TRUE), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
    #Preparations
    args <- list(...)
    
    if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
    if (is_not_null(v.threshold)) {
        v.threshold <- max(v.threshold, 1/v.threshold)
        disp.v.ratio <- TRUE
    }
    if (is_not_null(ks.threshold)) {
        if (ks.threshold > 1) {
            warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
            ks.threshold <- NULL
        }
        else disp.ks <- TRUE
    }
    if (is_null(weights) && is_null(subclass)) {
        un <- TRUE
        no.adj <- TRUE
    }
    else {
        no.adj <- FALSE
        if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
    }
    if (is_null(s.weights)) {
        s.weights <- rep(1, length(treat))
    }
    if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
    
    #Create new Target group
    target.name <- "Target"
    n <- length(treat)
    
    if (isTRUE(get.treat.type(treat) == "continuous") || (is.numeric(treat) && !is_binary(treat))) {#if continuous treatment
        covs <- data.frame(treat = treat, covs)
        treat <- factor(rep(c("All", target.name), each = n))
        target.summary <- FALSE
        which.treat <- NULL
        needs.summary <- FALSE
        treat.names <- unique.treat <- "All"
    }
    else {
        if (is.factor(treat) || is.character(treat)) {
            if (is.factor(treat)) treat.names <- unique.treat <- levels(treat)
            else treat.names <- unique.treat <- unique(treat, nmax = n - 1)
        }
        else {
            treat.names <- c("Control", "Treated")
            unique.treat <- sort(unique(treat, nmax = 2))
        }
        names(treat.names) <- unique.treat
        
        treat <- factor(c(treat.names[as.character(treat)], rep(target.name, n)))
        needs.summary <- TRUE
    }
    
    covs <- rbind(covs, covs)
    if (is_not_null(weights)) weights <- rbind(weights, as.data.frame(array(1, dim = dim(weights), 
                                                                            dimnames = dimnames(weights))))
    distance <- rbind(distance, distance)
    addl <- rbind(addl, addl)
    s.weights <- c(s.weights, s.weights)
    if (is_not_null(discarded)) discarded <- c(discarded, rep(FALSE, length(discarded)))
    s.d.denom <- "treated"
    
    treat.target.combinations <- lapply(treat.names, function(x) c(x, target.name))

    if (is_not_null(cluster)) {
        stop("Clusters are not yet supported with target balance assessment.", call. = FALSE)
    }
    else if (is_not_null(subclass)) {
        stop("Subclassification is not yet supported with target balance assessment.", call. = FALSE)
    }
    else {
        #Setup output object
        out.names <- c("Target.Balance", 
                       "Balance.Across.Treatments", 
                       "Observations", 
                       "call")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        
        if (any(treat.names == "Target")) stop ("\"Target\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
        args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
        balance.tables <- lapply(treat.target.combinations, function(t) do.call("base.bal.tab.binary", c(list(weights = weights[treat %in% t, , drop = FALSE], treat = factor(treat[treat %in% t], t), distance = distance[treat %in% t, , drop = FALSE], subclass = subclass[treat %in% t], covs = covs[treat %in% t, , drop = FALSE], call = NULL, int = int, poly = poly, addl = addl[treat %in% t, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster[treat %in% t], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[treat %in% t], discarded = discarded[treat %in% t], quick = quick), args), quote = TRUE))
        
        for (i in seq_along(balance.tables)) {
            names(balance.tables)[i] <- paste(treat.target.combinations[[i]], collapse = " vs. ")
            balance.tables[[i]][["Observations"]][[2]] <- NULL
        }
        
        out[["Target.Balance"]] <- balance.tables
        
        out[["Observations"]] <- samplesize.target(balance.tables, treat.names, target.name) 
        
        if (needs.summary && (target.summary || !quick)) {
        out[["Balance.Across.Treatments"]] <- balance.table.target.summary(balance.tables, 
                                                                     weight.names = names(weights),
                                                                     m.threshold = m.threshold,
                                                                     v.threshold = v.threshold,
                                                                     ks.threshold = ks.threshold,
                                                                     no.adj = no.adj,
                                                                     quick = quick,
                                                                     types = NULL)
        }
        
        out[["call"]] <- call
        
        attr(out, "print.options") <- list(m.threshold=m.threshold,
                                       v.threshold=v.threshold,
                                       ks.threshold=ks.threshold,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.adj=!no.adj, 
                                       which.cluster=which.cluster,
                                       cluster.summary=cluster.summary,
                                       cluster.fun = cluster.fun,
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       disp.means=disp.means, 
                                       disp.sds = disp.sds,
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks,
                                       disp.bal.tab = disp.bal.tab,
                                       nweights = ifelse(no.adj, 0, ncol(weights)),
                                       weight.names = names(weights),
                                       treat.names = treat.names,
                                       target.name = target.name,
                                       which.treat = which.treat,
                                       target.summary = target.summary,
                                       co.names = attr(out[["Target.Balance"]][[1]], "print.options")[["co.names"]])
        
        class(out) <- c("bal.tab.target", "bal.tab")
    }
    return(out)
    
}
