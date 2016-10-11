love.plot <- function(b, ...) UseMethod("love.plot")

love.plot.bal.tab <- function(b, stat = c("mean.diffs", "variance.ratios"), threshold = NULL, abs = FALSE, var.order = NULL, no.missing = TRUE, var.names = NULL, drop.distance = TRUE, cluster.fun = c("mean", "median", "max"), ...) {
  if (!"bal.tab" %in% class(b)) stop("The argument to \"b\" must be a bal.tab object, the output of a call to bal.tab().")
  stat <- match.arg(stat)
  which.stat <- switch(stat, mean.diffs = "Diff", variance.ratios = "V.Ratio")
  if (any(class(b) == "bal.tab.cluster")) {
      # if (length(b$print.options$which.cluster) > 1) stop("love.plot can only display balance for one cluster at a time. Make sure the argument to which.cluster in bal.tab() is the name or index of a single cluster.", call. = FALSE)
      # else if (!any(names(b$Cluster.Balance) == b$print.options$which.cluster[1]) && !any(seq_along(b$Cluster.Balance) == b$print.options$which.cluster[1])) stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
      # B <- b$Cluster.Balance[[b$print.options$which.cluster]]
      
      if (is.null(b$print.options$which.cluster) || is.na(b$print.options$which.cluster)) {
          #Use max, median, etc
          Cluster.Fun <- switch(match.arg(cluster.fun), mean = "Mean", median = "Median", max = "Max")
          title <- paste0("Covariate Balance\n", Cluster.Fun, " ", which.stat, " Across Clusters")
          which.stat <- paste(Cluster.Fun, which.stat, sep = ".")
          B <- b[["Cluster.Summary"]][, c("Type", 
                                     paste(which.stat, "Adj", sep = "."),  
                                     paste(which.stat, "Un", sep = "."))]
          abs <- TRUE
          
      }
      else if (is.character(b$print.options$which.cluster)) {
          cluster.names.good <- sapply(names(b$Cluster.Balance), function(x) any(sapply(b$print.options$which.cluster, function(y) isTRUE(all.equal(x, y)))))
          if (any(cluster.names.good)) {
              if (sum(cluster.names.good) == 1) {
                  B <- b[["Cluster.Balance"]][cluster.names.good][[1]]
                  title <- paste0("Covariate Balance\nCluster: ", b$print.options$which.cluster)
                  
              }
              else {
                  stop("love.plot can only display balance for one cluster at a time. Make sure the argument to which.cluster in bal.tab() is the name or index of a single cluster.", call. = FALSE)
              }
          }
          else {
              stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
          }
      }
      else if (is.numeric(b$print.options$which.cluster)) {
          cluster.numbers.good <- (seq_along(b$Cluster.Balance) %in% b$print.options$which.cluster)
          if (any(cluster.numbers.good)) {
              if (sum(cluster.numbers.good) == 1) {
                  B <- b$Cluster.Balance[cluster.numbers.good][[1]]
                  title <- paste0("Covariate Balance\nCluster: ", names(b$Cluster.Balance)[b$print.options$which.cluster])
              }
              else {
                  stop("love.plot can only display balance for one cluster at a time. Make sure the argument to which.cluster in bal.tab() is the name or index of a single cluster.", call. = FALSE)
              }
          }
          else {
              stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
          }
      }
  }
  else if (any(class(b) == "bal.tab.subclass")) {
    if (which.stat=="V.Ratio") stop("Variance ratios not currently supported for subclassification.", call. = FALSE)
    B <- b[["Balance.Across.Subclass"]]
    title <- "Covariate Balance\nAcross Subclasses"
  }
  else {
      B <- b[["Balance"]]
      title <- "Covariate Balance"
  }

  
  if (drop.distance) B <- B[is.na(match(B[,"Type"], "Distance")),]

  if (is.null(threshold)) {
    if (which.stat=="Diff") {
      if (any(names(B) == "M.Threshold") && !all(is.na(B[,"M.Threshold"]))) {
        threshold <- as.numeric(substring(names(table(B[,"M.Threshold"]))[2], 1+regexpr("[><]",names(table(B[,"M.Threshold"]))[2]), nchar(names(table(B[,"M.Threshold"]))[2])))
      }
    }
    else if (which.stat=="V.Ratio") {
      if (any(names(B) == "V.Threshold") && !all(is.na(B[,"V.Threshold"]))) {
        threshold <- as.numeric(substring(names(table(B[,"V.Threshold"]))[2], 1+regexpr("[><]",names(table(B[,"V.Threshold"]))[2]), nchar(names(table(B[,"V.Threshold"]))[2])))
      }
    }
  }
  var.labels <- row.names(B)
  if (!is.null(var.names)) {
    if (is.data.frame(var.names)) {
      if (ncol(var.names)==1) {
        if (any(sapply(var.labels, function(x) x %in% row.names(var.names)))) {
          which. <- match(var.labels, row.names(var.names))
          var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.,1]))
        }
        else if (nrow(var.names)==length(var.labels)) {
          var.labels <- ifelse(is.na(var.names) | var.names=="", var.labels, as.character(var.names[,1]))
        }
        else warning("var.names has only one column, but the row names do not correspond to variable names in the bal.tab object and the length of var.names differs from the number of variables to rename.", call. = FALSE)
      }
      else if (ncol(var.names)>1) {
        if (ncol(var.names)>2) warning("Only using first 2 columns of var.names", call. = FALSE)
        if (any(sapply(var.labels, function(x) x %in% as.character(var.names[,1])))) {
          which. <- match(var.labels, var.names[,1])
          var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.,2]))
        }
        else warning("var.names has more than one column, but the values in the first column do not correspond to variable names in the bal.tab object.", call. = FALSE)
      } 
    }
    else if (is.character(var.names) || is.factor(var.names)) {
      if (!is.null(names(var.names))) {
        if (any(sapply(var.labels, function(x) x %in% names(var.names)))) {
          which. <- match(var.labels, names(var.names))
          var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.]))
        }
        else warning("var.names is a named vector but the value names do not correspond to variable names in the bal.tab object and the length of var.names differs from the number of variables to rename.", call. = FALSE)
      }
      else if (length(var.names)==length(var.labels)) {
        var.labels <- ifelse(is.na(var.names) | var.names=="", var.labels, as.character(var.names))
      }
      else warning("var.names is a vector, but its values are unnamed and its length differs from the number of variables to rename.", call. = FALSE)
    }
    else warning("Argument to var.names is not one of the accepted structures and will be ignored.\n  See help(love.plot) for details.", immediate.=TRUE, call. = FALSE)
  }

  Sample <- NULL #To avoid CRAN checks rejecting Sample
  SS <- data.frame(var=rep(var.labels, 2), 
                   stat=c(B[, paste(which.stat, "Adj", sep = ".")], B[, paste(which.stat, "Un", sep = ".")]), 
                   Sample=c(rep("Adjusted", nrow(B)), rep("Unadjusted", nrow(B))))

  if (all(is.na(SS[, "stat"]))) stop("No balance statistics to display.", call. = FALSE)
  if (all(is.na(SS[SS$Sample=="Adjusted", "stat"]))) SS <- SS[SS[, "Sample"]=="Unadjusted",]
  if (abs) {
    SS[, "stat"] <- abs(SS[, "stat"])
    dec <- FALSE}
  else dec <- TRUE
  if (!is.null(var.order)) {
    var.order <- match.arg(var.order, c("adjusted", "unadjusted")) 
    SS[, "var"] <- factor(SS[, "var"], levels=SS[order(SS[tolower(SS[, "Sample"])==var.order, "stat"], decreasing = dec), "var"])
  }
  else SS[, "var"] <- factor(SS[, "var"], levels = unique(SS[, "var"])[order(unique(SS[, "var"]), decreasing = TRUE)])
  SS[, "Sample"] <- factor(SS[, "Sample"], levels = c("Unadjusted", "Adjusted"))
  if (stat=="mean.diffs" & any(abs(SS[, "stat"]) > 5)) warning("Large mean differences detected; you may not be using standardizied mean differences for continuous variables. To do so, specify continuous=\"std\" in bal.tab().", call.=FALSE, noBreaks.=TRUE)
  if (no.missing) SS <- SS[!is.na(SS[, "stat"]),]
  
  #Make the plot
  #library(ggplot2)
  lp <- ggplot(SS, aes(y = var, x = stat, color = Sample)) + geom_point() + labs(title = title, y = "")
  if (stat=="mean.diffs") {
    lp <- lp + geom_vline(xintercept = 0, linetype=1, alpha=.4)
    if (abs) {
      lp <- lp + xlab("Absolute Mean Differences") 
      if (!is.null(threshold)) lp <- lp + geom_vline(xintercept = abs(threshold), linetype=2)
    }
    else {
      lp <- lp + xlab("Mean Differences") 
      if (!is.null(threshold)) lp <- lp + geom_vline(xintercept = c(-threshold, threshold), linetype=2)
    }
  }
  else if (stat=="variance.ratios") {
    lp <- lp + xlab("Variance Ratios") + geom_vline(xintercept = 1, linetype=1, alpha=.4)
    if (!is.null(threshold)) lp <- lp + geom_vline(xintercept = max(threshold,1/threshold), linetype=2)
  }

  return(lp)
}
love.plot.bal.tab.cont <- function(b, threshold = NULL, abs = FALSE, var.order = NULL, no.missing = TRUE, var.names = NULL, drop.distance = TRUE, cluster.fun = c("mean", "median", "max"), ...) {
    #if (!"bal.tab" %in% class(b)) stop("The argument to \"b\" must be a bal.tab object, the output of a call to bal.tab().")
    which.stat <- "Corr"
    
    if (any(class(b) == "bal.tab.cont.cluster")) {
        # if (length(b$print.options$which.cluster) > 1) stop("love.plot can only display balance for one cluster at a time. Make sure the argument to which.cluster in bal.tab() is the name or index of a single cluster.", call. = FALSE)
        # else if (!any(names(b$Cluster.Balance) == b$print.options$which.cluster[1]) && !any(seq_along(b$Cluster.Balance) == b$print.options$which.cluster[1])) stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
        # B <- b$Cluster.Balance[[b$print.options$which.cluster]]
        
        if (is.null(b$print.options$which.cluster) || is.na(b$print.options$which.cluster)) {
            #Use max, median, etc
            Cluster.Fun <- switch(match.arg(cluster.fun), mean = "Mean", median = "Median", max = "Max")
            title <- paste0("Covariate Balance\n", Cluster.Fun, " ", which.stat, " Across Clusters")
            which.stat <- paste(Cluster.Fun, which.stat, sep = ".")
            B <- b$Cluster.Summary[, c("Type", 
                                       paste(which.stat, "Adj", sep = "."),  
                                       paste(which.stat, "Un", sep = "."))]
            abs <- TRUE
            
        }
        else if (is.character(b$print.options$which.cluster)) {
            cluster.names.good <- sapply(names(b$Cluster.Balance), function(x) any(sapply(b$print.options$which.cluster, function(y) isTRUE(all.equal(x, y)))))
            if (any(cluster.names.good)) {
                if (sum(cluster.names.good) == 1) {
                    B <- b$Cluster.Balance[cluster.names.good][[1]]
                    title <- paste0("Covariate Balance\nCluster: ", b$print.options$which.cluster)
                    
                }
                else {
                    stop("love.plot can only display balance for one cluster at a time. Make sure the argument to which.cluster in bal.tab() is the name or index of a single cluster.", call. = FALSE)
                }
            }
            else {
                stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
            }
        }
        else if (is.numeric(b$print.options$which.cluster)) {
            cluster.numbers.good <- (seq_along(b[["Cluster.Balance"]]) %in% b$print.options$which.cluster)
            if (any(cluster.numbers.good)) {
                if (sum(cluster.numbers.good) == 1) {
                    B <- b[["Cluster.Balance"]][cluster.numbers.good][[1]]
                    title <- paste0("Covariate Balance\nCluster: ", names(b[["Cluster.Balance"]])[b$print.options$which.cluster])
                }
                else {
                    stop("love.plot can only display balance for one cluster at a time. Make sure the argument to which.cluster in bal.tab() is the name or index of a single cluster.", call. = FALSE)
                }
            }
            else {
                stop("Make sure the argument to which.cluster in bal.tab() is a valid name or index of a cluster.", call. = FALSE)
            }
        }
    }
    else if (any(class(b) == "bal.tab.subclass")) {
        B <- b[["Balance.Across.Subclass"]]
        title <- "Covariate Balance\nAcross Subclasses"
    }
    else {
        B <- b[["Balance"]]
        title <- "Covariate Balance"
    }
    
    if (drop.distance) B <- B[-match("Distance", B[,"Type"]),]
    
    if (is.null(threshold)) {
        if (any(names(B) == "R.Threshold") && !all(is.na(B[,"R.Threshold"]))) {
            threshold <- as.numeric(substring(names(table(B[,"R.Threshold"]))[2], 1+regexpr("[><]",names(table(B[,"R.Threshold"]))[2]), nchar(names(table(B[,"R.Threshold"]))[2])))
        }
    }
    var.labels <- row.names(B)
    if (!is.null(var.names)) {
        if (is.data.frame(var.names)) {
            if (ncol(var.names)==1) {
                if (any(sapply(var.labels, function(x) x %in% row.names(var.names)))) {
                    which. <- match(var.labels, row.names(var.names))
                    var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.,1]))
                }
                else if (nrow(var.names)==length(var.labels)) {
                    var.labels <- ifelse(is.na(var.names) | var.names=="", var.labels, as.character(var.names[,1]))
                }
                else warning("var.names has only one column, but the row names do not correspond to variable names in the bal.tab object and the length of var.names differs from the number of variables to rename.", call. = FALSE)
            }
            else if (ncol(var.names)>1) {
                if (ncol(var.names)>2) warning("Only using first 2 columns of var.names", call. = FALSE)
                if (any(sapply(var.labels, function(x) x %in% as.character(var.names[,1])))) {
                    which. <- match(var.labels, var.names[,1])
                    var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.,2]))
                }
                else warning("var.names has more than one column, but the values in the first column do not correspond to variable names in the bal.tab object.", call. = FALSE)
            } 
        }
        else if (is.character(var.names) || is.factor(var.names)) {
            if (!is.null(names(var.names))) {
                if (any(sapply(var.labels, function(x) x %in% names(var.names)))) {
                    which. <- match(var.labels, names(var.names))
                    var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.]))
                }
                else warning("var.names is a named vector but the value names do not correspond to variable names in the bal.tab object and the length of var.names differs from the number of variables to rename.", call. = FALSE)
            }
            else if (length(var.names)==length(var.labels)) {
                var.labels <- ifelse(is.na(var.names) | var.names=="", var.labels, as.character(var.names))
            }
            else warning("var.names is a vector, but its values are unnamed and its length differs from the number of variables to rename.", call. = FALSE)
        }
        else warning("Argument to var.names is not one of the accepted structures and will be ignored.\n  See help(love.plot) for details.", immediate.=TRUE, call. = FALSE)
    }
    
    Sample <- stat <- NULL #To avoid CRAN checks rejecting Sample
    SS <- data.frame(var = rep(var.labels, 2), 
                     stat = c(B[, paste0(which.stat,".Adj")], B[, paste0(which.stat,".Un")]), 
                     Sample = c(rep("Adjusted", nrow(B)), rep("Unadjusted", nrow(B))))
    if (all(is.na(SS[, "stat"]))) stop("No balance statistics to display.", call. = FALSE)
    if (all(is.na(SS[SS[, "Sample"]=="Adjusted", "stat"]))) SS <- SS[SS[, "Sample"]=="Unadjusted",]
    if (abs) {
        SS$stat <- abs(SS[, "stat"])
        dec <- FALSE}
    else dec <- TRUE
    if (length(var.order) > 0) {
        var.order <- match.arg(var.order, c("adjusted", "unadjusted")) 
        SS[, "var"] <- factor(SS[, "var"], levels = SS[order(SS[tolower(SS[, "Sample"])==var.order, "stat"], decreasing = dec), "var"])
    }
    else SS[, "var"] <- factor(SS[, "var"], levels = unique(SS[, "var"])[order(unique(SS[, "var"]), decreasing = TRUE)])
    SS[, "Sample"] <- factor(SS[, "Sample"], levels = c("Unadjusted", "Adjusted"))
    if (no.missing) SS <- SS[!is.na(SS[, "stat"]),]
    
    #Make the plot
    #library(ggplot2)
    lp <- ggplot(SS, aes(y = var, x = stat, color = Sample)) + geom_point() + labs(title = title, y = "")
    lp <- lp + geom_vline(xintercept = 0, linetype=1, alpha=.4)
    if (abs) {
        lp <- lp + xlab("Absolute Treatment Correlations") 
        if (!is.null(threshold)) lp <- lp + geom_vline(xintercept = abs(threshold), linetype=2)
    }
    else {
        lp <- lp + xlab("Treatment Correlations") 
        if (!is.null(threshold)) lp <- lp + geom_vline(xintercept = c(-threshold, threshold), linetype=2)
    }
    
    return(lp)
}