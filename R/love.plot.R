love.plot <- function(b, stat=c("mean.diffs", "variance.ratios"), threshold=NULL, abs=FALSE, var.order=NULL, no.missing=FALSE, var.names=NULL, drop.distance=TRUE) {

  if (!"bal.tab" %in% class(b)) stop("The argument to \"b\" must be a bal.tab object, the output of a call to bal.tab().")
  stat <- match.arg(stat)
  which.stat <- switch(stat, mean.diffs="Diff", variance.ratios="V.Ratio")
  
  if("bal.tab.subclass" %in% class(b)) {
    if (which.stat=="V.Ratio") stop("Variance ratios not currently supported for subclassification.")
    B <- b$Balance.Across.Subclass
  }
  else B <- b$Balance
  
  if (drop.distance) B <- B[which(!B$Type %in% "Distance"),]
  
  if (is.null(threshold)) {
    if (which.stat=="Diff") {
      if (!all(is.na(B$M.Threshold))) {
        threshold <- as.numeric(substring(names(table(B$M.Threshold))[2], 1+regexpr("[><]",names(table(B$M.Threshold))[2]), nchar(names(table(B$M.Threshold))[2])))
      }
    }
    else if (which.stat=="V.Ratio") {
      if (!all(is.na(B$V.Threshold))) {
        threshold <- as.numeric(substring(names(table(B$V.Threshold))[2], 1+regexpr("[><]",names(table(B$V.Threshold))[2]), nchar(names(table(B$V.Threshold))[2])))
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
        else warning("var.names has only one column, but the row names do not correspond to variable names in the bal.tab object and the length of var.names differs from the number of variables to rename.")
      }
      else if (ncol(var.names)>1) {
        if (ncol(var.names)>2) warning("Only using first 2 columns of var.names")
        if (any(sapply(var.labels, function(x) x %in% as.character(var.names[,1])))) {
          which. <- match(var.labels, var.names[,1])
          var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.,2]))
        }
        else warning("var.names has more than one column, but the values in the first column do not correspond to variable names in the bal.tab object.")
      } 
    }
    else if (is.character(var.names) || is.factor(var.names)) {
      if (!is.null(names(var.names))) {
        if (any(sapply(var.labels, function(x) x %in% names(var.names)))) {
          which. <- match(var.labels, names(var.names))
          var.labels <- ifelse(is.na(which.), var.labels, as.character(var.names[which.]))
        }
        else warning("var.names is a named vector but the value names do not correspond to variable names in the bal.tab object and the length of var.names differs from the number of variables to rename.")
      }
      else if (length(var.names)==length(var.labels)) {
        var.labels <- ifelse(is.na(var.names) | var.names=="", var.labels, as.character(var.names))
      }
      else warning("var.names is a vector, but its values are unnamed and its length differs from the number of variables to rename.")
    }
    else warning("Argument to var.names is not one of the accepted structures and will be ignored.\n  See help(love.plot) for details.", immediate.=TRUE)
  }
  
  SS <- data.frame(var=rep(var.labels, 2), 
                   stat=c(B[,paste0(which.stat,".Adjusted")], B[,paste0(which.stat,".Unadjusted")]), 
                   sample=c(rep("Adjusted", nrow(B)), rep("Unadjusted", nrow(B))))
  if (all(is.na(SS$stat))) stop("No balance statistics to display.")
  if (all(is.na(SS$stat[SS$sample=="Adjusted"]))) SS <- SS[SS$sample=="Unadjusted",]
  if (abs) {
    SS$stat <- abs(SS$stat)
    dec <- FALSE}
  else dec <- TRUE
  if (!is.null(var.order)) {
    var.order <- match.arg(var.order, c("adjusted", "unadjusted")) 
    SS$var <- factor(SS$var, levels=SS$var[order(SS$stat[tolower(SS$sample)==var.order], decreasing=dec)])
  }
  else SS$var <- factor(SS$var, levels=unique(SS$var)[order(unique(SS$var), decreasing=TRUE)])
  SS$sample <- factor(SS$sample, levels=c("Unadjusted", "Adjusted"))
  if (stat=="mean.diffs" & any(abs(SS$stat)>5)) warning("Large mean differences detected; you may not be using standardizied mean differences for continuous variables. To do so, specify continuous=\"std\" in i.sum().", call.=FALSE, noBreaks.=TRUE)
  if (no.missing) SS <- SS[!is.na(SS$stat),]
  
  #Make the plot
  #library(ggplot2)
  lp <- ggplot(SS, aes(y=var, x=stat, color=sample)) + geom_point() + labs(title="Covariate Balance", y="")
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