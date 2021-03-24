#bal.tab
process_obj <- function(obj) {
  if (is_null(obj)) obj <- list()
  else {
    if (isS4(obj)) obj <- asS3(obj)
    
    #npCBPS
    if (is_(obj, "npCBPS")) {
      class(obj) <- c("CBPS", "npCBPS")
    }
    #ebalance.trim
    else if (is_(obj, "ebalance.trim")) {
      class(obj) <- "ebalance"
    }
    #Matchby
    else if (is_(obj, "Matchby")) {
      class(obj) <- "Match"
    }
    #cem.match.list
    else if (is_(obj, "cem.match.list")) {
      class(obj) <- c("cem.match", "cem.match.list")
    }
    #time.list
    else if (is_(obj, "list") && all(vapply(obj, rlang::is_formula, logical(1)))) {
      class(obj) <- c("formula.list", "time.list", class(obj))
    }
    else if (is_(obj, "list") && all(vapply(obj, is.data.frame, logical(1)))) {
      class(obj) <- c("data.frame.list", "time.list", class(obj))
    }
    #designmatch
    else if (is_(obj, "list") && length(obj) >= 6) {
      dm.b.names <- c("obj_total", "obj_dist_mat", "t_id", 
                      "c_id", "group_id", "time")
      dm.n.names <- c("obj_total", "obj_dist_mat", "id_1", 
                      "id_2", "group_id", "time")
      if (all(dm.b.names %in% names(obj)) || all(dm.n.names %in% names(obj))) {
        class(obj) <- c("designmatch")
      }
    }
  }
  class(obj) <- c(class(obj), "cobalt.processed.obj")
  return(obj)
}

#x2base
process_treat <- function(treat, datalist = list()) {
  
  if (missing(treat)) stop("'treat' must be specified.", call. = FALSE)
  
  if (inherits(treat, "unprocessed.treat")) {
    attrs <- attributes(treat)
    renamed_original <- setNames(names(treat_vals(treat)), treat_vals(treat))
    treat <- factor(renamed_original[as.character(treat)], levels = renamed_original)
    for (at in c("treat_names", "treat_vals", "treat.type", "names"))
      attr(treat, at) <- attrs[[at]]
  }
  else {
    keep_values <- has.treat.type(treat) && get.treat.type(treat) == "multinomial"
    
    treat <- vector.process(treat, name = "treat", 
                            which = "treatment statuses", 
                            datalist = datalist, missing.okay = FALSE)
    
    treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type == "binary") {
      if (!is.factor(treat)) treat <- factor(treat, levels = sort(unique(treat, nmax = 2)))
      original_values <- levels(treat)
      if (!keep_values && can_str2num(as.character(treat)) && all(original_values %in% c("0", "1"))) {
        treat_names(treat) <- setNames(c("Control", "Treated"), c("control", "treated"))
      }
      else {
        treat_names(treat) <- setNames(original_values, c("control", "treated"))
      }
      
      treat_vals(treat) <- setNames(original_values, treat_names(treat))
    }
    else if (treat.type == "multinomial") {
      treat <- factor(treat)
      treat_names(treat) <- setNames(levels(treat), levels(treat))
      treat_vals(treat) <- setNames(levels(treat), treat_names(treat))
    }
    attr(treat, "treat.type") <- treat.type
  }
  class(treat) <- c("processed.treat", class(treat))
  return(treat)
}
unprocess_treat <- function(treat) {
  if (inherits(treat, "processed.treat")) {
    attrs <- attributes(treat)
    treat <- treat_vals(treat)[as.character(treat)]
    attributes(treat) <- attrs
    class(treat) <- c("unprocessed.treat", class(treat_vals(treat)))
  }
  return(treat)
}
process_treat.list <- function(treat.list, datalist = list()) {
  if (is_null(treat.list)) stop("'treat.list' must be specified.", call. = FALSE)
  if (!is_(treat.list, "list")) {
    treat.list <- as.list(treat.list)
  }
  
  treat.list.names <- vapply(seq_along(treat.list), function(ti) {
    if (is.character(treat.list[[ti]]) && length(treat.list[[ti]])==1L && is_not_null(datalist)) {
      treat.list[[ti]]
    }
    else if (rlang::is_named(treat.list)) names(treat.list)[ti]
    else as.character(ti)
  }, character(1L))
  treat.list <- lapply(treat.list, process_treat, datalist = datalist)
  names(treat.list) <- treat.list.names
  
  return(treat.list)
}
`treat_names<-` <- function(treat, value) {
  `attr<-`(treat, "treat_names", value)
}
treat_names <- function(treat) {
  attr(treat, "treat_names")
}
`treat_vals<-` <- function(treat, value) {
  `attr<-`(treat, "treat_vals", value)
}
treat_vals <- function(treat) {
  attr(treat, "treat_vals")
}
`[.processed.treat` <- function(x, ...) {
  y <- NextMethod("[")
  treat_names(y) <- treat_names(x)
  treat_vals(y) <- treat_vals(x)
  attr(y, "treat.type") <- attr(x, "treat.type")
  class(y) <- class(x)
  y
}
`[.unprocessed.treat` <- `[.processed.treat`

initialize_X <- function() {
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "estimand",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass",
               "stats",
               "thresholds")
  X <- make_list(X.names)
  return(X)
}
initialize_X_msm <- function() {
  X.names <- c("covs.list",
               "treat.list",
               "weights",
               "distance.list",
               "addl.list",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass",
               "stats",
               "thresholds")
  X <- make_list(X.names)
  return(X)
}
weight.check <- function(w) {
  wname <- deparse1(substitute(w))
  if (!is.list(w)) w <- list(w)
  if (anyNA(w, recursive = TRUE)) stop(paste0("NAs are not allowed in the ", wname, "."), call. = FALSE)
  for (x in w) {if (!all(is.numeric(x))) stop(paste0("All ", wname, " must be numeric."), call. = FALSE)}
  for (x in w) {if (!all(is.finite(x))) stop(paste0("Infinite ", wname, " are not allowed."), call. = FALSE)}
  for (x in w) {if (any(x < 0)) {
    warning(paste0("Negative ", wname, " found. This may yield nonsensical results and errors as the square root of negative weights is encountered."), call. = FALSE)
    break
  }}
}
cluster.check <- function(cluster, treat) {
  if (!is.list(treat)) treat <- list(treat)
  
  stop_warn <- c(cw = FALSE, bw = FALSE, bs = FALSE)
  for (t in treat) {
    if (!has.treat.type(t)) t <- assign.treat.type(t)
    treat.type <- get.treat.type(t)
    
    if (treat.type == "continuous" && !stop_warn["cw"]) {
      tabu <- tabulate(cluster, nbins = nlevels(cluster))
      if (any(tabu == 1)) stop_warn["cw"] <- TRUE
    }
    else if (treat.type != "continuous" && !all(stop_warn[c("bw", "bs")])) {
      tab <- table(cluster, t)
      if (any(tab == 0)) stop_warn["bs"] <- TRUE
      else if (any(tab == 1)) stop_warn["bw"] <- TRUE
    }
  }
  if (stop_warn["cw"]) warning("Some clusters have only one unit in them, which may yield unxpected results.", call. = FALSE)
  if (stop_warn["bw"]) warning("Some clusters have only one member of a treatment group in them, which may yield unxpected results.", call. = FALSE)
  if (stop_warn["bs"]) stop("Not all treatment levels are present in all clusters.", call. = FALSE)
  
}
strata2weights <- function(strata, treat, estimand = NULL, focal = NULL) {
  #Process strata into weights (similar to weight.subclass from MatchIt)
  
  #Checks
  if (!is_(strata, "atomic") || is_not_null(dim(strata))) {
    stop("'strata' must be an atomic vector or factor.", call. = FALSE)
  }
  #Process treat
  treat <- process_treat(treat)
  
  if (get.treat.type(treat) == "continuous") {
    stop("'strata' cannot be turned into weights for continuous treatments.", call. = FALSE)
  }
  
  s.d.denom <- get.s.d.denom(NULL, estimand = estimand, subclass = strata, treat = treat, focal = focal, quietly = TRUE)
  if (s.d.denom %in% treat_vals(treat)) focal <- process_focal(s.d.denom, treat)
  else focal <- NULL
  
  NAsub <- is.na(strata)
  imat <- do.call("cbind", lapply(treat_vals(treat), function(t) treat == t & !NAsub))
  colnames(imat) <- treat_vals(treat)
  
  weights <- rep(0, length(treat))
  
  if (!is.factor(strata)) {
    strata <- factor(strata, nmax = min(colSums(imat)))
    levels(strata) <- seq_len(nlevels(strata))
  }
  
  t_by_sub <- do.call("rbind", lapply(treat_vals(treat), function(t) tabulate(strata[imat[,t]], nlevels(strata))))
  dimnames(t_by_sub) <- list(treat_vals(treat), levels(strata))
  
  total_by_sub <- colSums(t_by_sub)
  
  strata.c <- as.character(strata)
  
  if (is_not_null(focal)) {
    focal <- process_focal(focal, treat)
    for (t in treat_vals(treat)) {
      if (t == focal) weights[imat[,t]] <- 1
      else weights[imat[,t]] <- (t_by_sub[focal,]/t_by_sub[t,])[strata.c[imat[,t]]]
    }
  }
  else {
    for (t in treat_vals(treat)) {
      weights[imat[,t]] <- (total_by_sub/t_by_sub[t,])[strata.c[imat[,t]]]
    }
  }
  
  if (any(na.w <- !is.finite(weights))) {
    weights[na.w] <- 0
    warning("Some units were given weights of zero due to zeros in stratum membership.", call. = FALSE)
  }
  
  if (all(check_if_zero(weights))) 
    stop("No units were stratified", call. = FALSE)
  else {
    for (tnn in names(treat_names(treat))) {
      if (all(check_if_zero(weights[treat == treat_vals(treat)[treat_names(treat)[tnn]]])))
        stop(paste("No", tnn, "units were stratified."), call. = FALSE)
    }
  }
  
  attr(weights, "match.strata") <- strata
  return(weights)
}

use.tc.fd <- function(formula = NULL, data = NULL, treat = NULL, covs = NULL, needs.treat = TRUE, needs.covs = TRUE) {
  if (is_(formula, "formula")) {
    D <- NULL
    if (is_not_null(data)) D <- data
    if (is_not_null(covs) && is_(covs, c("data.frame", "matrix"))) {
      if (is_not_null(D)) D <- cbind(D, as.data.frame(covs)) 
      else D <- as.data.frame(covs)
    }
    treat <- get_treat_from_formula(formula, D, treat = treat)
    covs <- get_covs_from_formula(formula, D)
    t.c <- list(treat = treat, covs = covs, treat.name = attr(treat, "treat.name"))
    attr(t.c, "which") <- "fd"
  }
  else {
    if (is_not_null(covs)) {
      if (is_(covs, c("data.frame", "matrix"))) covs <- get_covs_from_formula(data = covs)
      else if (is.character(covs)) {
        if (is_not_null(data) && is_(data, c("data.frame", "matrix"))) {
          if (any(covs %nin% colnames(data))) {
            stop("All entries in 'covs' must be names of variables in 'data'.", call. = FALSE)
          }
          covs <- get_covs_from_formula(f.build(covs), data = as.data.frame(data))
        }
        else {
          stop("If 'covs' is a character vector, 'data' must be specified as a data.frame.", call. = FALSE)
        }
      }
      else stop("'covs' must be a data.frame of covariates.", call. = FALSE)
    }
    
    if (is_not_null(treat)) {
      if (!is_(treat, "atomic")) stop("'treat' must be an vector of treatment statuses.", call. = FALSE)
      if (is.character(treat) && length(treat) == 1) treat <- get_treat_from_formula(f.build(treat, "."), data = data)
      else treat <- get_treat_from_formula(treat ~ ., treat = treat)
    }
    
    t.c <- list(treat = treat, covs = covs)
    attr(t.c, "which") <- "tc"
  }
  
  if (needs.covs && is_null(t.c[["covs"]])) stop("No covariates were specified.", call. = FALSE)
  if (needs.treat && is_null(t.c[["treat"]])) stop("No treatment variable was specified.", call. = FALSE)
  
  return(t.c)
}
process.val <- function(val, i, treat = NULL, covs = NULL, addl.data = list(), ...) {
  if (is.numeric(val)) {
    val.df <- setNames(data.frame(val), i)
  }
  else if (is.character(val)) {
    addl.data <- addl.data[!vapply(addl.data, is_null, logical(1))]
    if ((is_not_null(addl.data) && length(val) > max(vapply(addl.data, ncol, numeric(1)))) || (is_null(covs) && is_null(treat)) || length(val) == NROW(covs) || length(val) == length(treat)){
      val.df <- setNames(data.frame(val), i)
    }
    else {
      if (is_not_null(addl.data)) {
        val <- unique(val)
        val.df <- make_df(val, nrow = max(vapply(addl.data, nrow, numeric(1))))
        not.found <- rlang::rep_named(val, FALSE)
        for (v in val) {
          found <- FALSE
          k <- 1
          while (found == FALSE && k <= length(addl.data)) {
            if (v %in% names(addl.data[[k]])) {
              val.df[[v]] <- addl.data[[k]][[v]]
              found <- TRUE
            }
            else k <- k + 1
          }
          if (!found) not.found[v] <- TRUE
        }
        if (any(not.found)) {
          warning(paste0("The following variable(s) named in '", i, "' are not in any available data sets and will be ignored: ",
                         paste(val[not.found])), call. = FALSE)
          val.df <- val.df[!not.found]
        }
      }
      else {
        val.df <- NULL
        warning(paste0("Names were provided to '", i, "', but no argument to 'data' was provided. Ignoring '", i,"'."), 
                call. = FALSE)
      }
    }
  }
  else if (is.data.frame(val)) {
    val.df <- val
  }
  else {
    if (i == "weights") stop("The argument supplied to 'weights' must be a named list of weights, names of variables containing weights in an available data set, or objects with a get.w() method.", call. = FALSE)
    else stop(paste0("The argument supplied to '", i, "' must be a vector, a data.frame, or the names of variables in an available data set."), call. = FALSE)
  }
  
  return(val.df)
}
data.frame.process <- function(i, df, treat = NULL, covs = NULL, addl.data = list(), ...) {
  val <- df
  val.df <- NULL
  if (is_not_null(val)) {
    if (i == "weights" && any(has_method(class(val), "get.w"))) {
      val <- list(val)
    }
    if (is_(val, "list")) {
      if (i == "weights") {
        #Use get.w() on inputs
        for (x in seq_along(val)) {
          val[[x]] <- process_obj(val[[x]])
          if (any(has_method(class(val[[x]]), "get.w"))) {
            get.w_class <- class(val[[x]])[has_method(class(val[[x]]), "get.w")][1]
            val[[x]] <- get.w(val[[x]], treat = treat, ...)
            if (rlang::names2(val)[x] == "") names(val)[x] <- get.w_class
          }
        }
        val.list <- lapply(val, function(x) process.val(x, i, treat, covs, addl.data = addl.data))
      }
      else {
        val.list <- lapply(val, function(x) process.val(x, i, treat, covs, addl.data = addl.data))
      }
      
      if (!rlang::is_named(val.list)) {
        stop(paste0("All entries in '", i, "' must have names."), call. = FALSE)
      }
      val.list <- lapply(seq_along(val.list), function(x) {
        if (NCOL(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
        return(val.list[[x]])})
      if (!all_the_same(vapply(val.list, nrow, numeric(1)))) {
        stop(paste0("Not all items in '", i, "' have the same length."), call. = FALSE)
      }
      val.df <- setNames(do.call("cbind", val.list),
                         unlist(lapply(val.list, names)))
    }
    else {
      val.df <- process.val(val, i, treat, covs, addl.data = addl.data)
    }
  }
  return(val.df)
}
list.process <- function(i, List, ntimes, call.phrase, treat.list = list(), covs.list = list(), addl.data = list(), ...) {
  val.List <- List
  if (is_not_null(val.List)) {
    if (!is_(val.List, "list")) {
      val.List <- list(val.List)
    }
    if (length(val.List) == 1) {
      val.List <- replicate(ntimes, val.List)
    }
    else if (length(val.List) != ntimes) {
      stop(paste0("The argument to '", i, "' must be a list of the same length as the number of time points in ",  call.phrase, "."), call. = FALSE)
    }
    for (ti in seq_along(val.List)) {
      val <- val.List[[ti]]
      val.df <- NULL
      if (is_not_null(val)) {
        if (is_(val, "list")) {
          val.list <- lapply(val, function(x) process.val(x, strsplit(i, ".list", fixed = TRUE)[[1]], treat.list[[ti]], covs.list[[ti]], addl.data = addl.data))
          val.list <- lapply(seq_along(val.list), function(x) {
            if (NCOL(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
            val.list[[x]]})
          if (!all_the_same(vapply(val.list, nrow, numeric(1)))) {
            stop(paste0("Not all items in '", i, "' have the same length."), call. = FALSE)
          }
          
          val.df <- setNames(do.call("cbind", val.list),
                             vapply(val.list, names, character(1)))
        }
        else {
          val.df <- process.val(val, strsplit(i, ".list", fixed = TRUE)[[1]], treat.list[[ti]], covs.list[[ti]], addl.data = addl.data)
        }
        if (is_not_null(val.df) && anyNA(val.df)) {
          stop(paste0("Missing values exist in '", i, "'."), call. = FALSE)
        }
        val.List[[ti]] <- val.df
      }
      
    }
    val.df.lengths <- vapply(val.List[lengths(val.List) > 0], nrow, numeric(1))
    if (max(val.df.lengths) != min(val.df.lengths)) {
      stop(paste0("All columns in '", i, "' need to have the same number of rows."), call. = FALSE)
    }
  }
  return(val.List)
}
vector.process <- function(vec, name = deparse1(substitute(vec)), which = name, datalist = list(), missing.okay = FALSE) {
  bad.vec <- FALSE
  if (is.character(vec) && length(vec)==1L && is_not_null(datalist)) {
    for (i in seq_along(datalist)) {
      if (is_(datalist[[i]], "matrix") && vec %in% colnames(datalist[[i]])) {
        vec <- datalist[[i]][,vec]
        break
      }
      else if (is_(datalist[[i]], "data.frame") && vec %in% names(datalist[[i]])) {
        vec <- datalist[[i]][[vec]]
        break
      }
      else if (i == length(datalist)) bad.vec <- TRUE
    }
  }
  else if (is_(vec, "atomic") && length(vec) > 1L) {
    vec <- vec
  }
  else {
    bad.vec <- TRUE
  }
  
  if (bad.vec) stop(paste0("The argument to '", name, "' must be a vector of ", which, " or the (quoted) name of a variable in 'data' that contains ", which, "."), call. = FALSE)
  
  if (!missing.okay && anyNA(vec)) stop(paste0("Missing values exist in '", name, "'."), call. = FALSE)
  
  return(vec)
} 
get.s.d.denom <- function(s.d.denom = NULL, estimand = NULL, weights = NULL, subclass = NULL, treat = NULL, focal = NULL, quietly = FALSE) {
  check.estimand <- check.weights <- check.focal <- bad.s.d.denom <- bad.estimand <- FALSE
  s.d.denom.specified <- !missing(s.d.denom) && is_not_null(s.d.denom)
  estimand.specified <- is_not_null(estimand)
  treat.is.processed <- is_(treat, "processed.treat")
  
  if (s.d.denom.specified) {
    if (!treat.is.processed) {
      treat <- process_treat(treat)
    }
    unique.treats <- as.character(treat_vals(treat))
    allowable.s.d.denoms <- c("pooled", "all", "weighted", "hedges")
    if (length(treat_names(treat)) == 2 && all(c("treated", "control") %in% names(treat_names(treat))))
      allowable.s.d.denoms <- c(allowable.s.d.denoms, "treated", "control")
    if (is_not_null(focal)) allowable.s.d.denoms <- c(allowable.s.d.denoms, "focal")
    
    if (length(s.d.denom) > 1 && length(s.d.denom) != NCOL(weights)) {
      stop(paste0("'s.d.denom' must have length 1 or equal to the number of valid sets of weights, which is ", NCOL(weights), "."), call. = FALSE)
    }
    
    s.d.denom <- match_arg(s.d.denom, unique(c(unique.treats, allowable.s.d.denoms)), 
                           several.ok = TRUE)
    
    if (any(s.d.denom %in% c("treated", "control")) && length(treat_names(treat) == 2)) {
      s.d.denom[s.d.denom %in% c("treated", "control")] <- treat_vals(treat)[treat_names(treat)[s.d.denom[s.d.denom %in% c("treated", "control")]]]
    }
    else if (any(s.d.denom %in% "focal")) check.focal <- TRUE
    
  }
  else {
    check.estimand <- TRUE
  }
  
  if (check.estimand) {
    if (estimand.specified) {
      if (!treat.is.processed) treat <- process_treat(treat)
      try.estimand <- tryCatch(match_arg(toupper(estimand), c("ATT", "ATC", "ATE", "ATO", "ATM"), several.ok = TRUE),
                               error = function(cond) NA_character_)
      if (anyNA(try.estimand) || any(try.estimand %in% c("ATC", "ATT")) && get.treat.type(treat) != "binary") {
        check.focal <- TRUE
      }
      else {
        if (length(try.estimand) > 1 && length(try.estimand) != NCOL(weights)) {
          stop(paste0("'estimand' must have length 1 or equal to the number of valid sets of weights, which is ", NCOL(weights), "."), call. = FALSE)
        }
        else s.d.denom <- vapply(try.estimand, function(x) switch(x, 
                                                                  ATT = treat_vals(treat)[treat_names(treat)["treated"]], 
                                                                  ATC = treat_vals(treat)[treat_names(treat)["control"]], 
                                                                  ATO = "weighted",
                                                                  ATM = "weighted",
                                                                  "pooled"), FUN.VALUE = character(1L))
      }
    }
    else {
      check.focal <- TRUE
    }
  }
  if (check.focal) {
    if (is_not_null(focal) && NCOL(weights) == 1) {
      s.d.denom <- focal
    }
    else check.weights <- TRUE
  }
  if (check.weights) {
    if (!treat.is.processed) treat <- process_treat(treat)
    if (is_null(weights) && is_null(subclass)) {
      s.d.denom <- "pooled"
    }
    else if (is_not_null(subclass)) {
      sub.tab <- table(treat, subclass)[treat_vals(treat), ]
      sub.tab <- rbind(sub.tab, table(subclass)[colnames(sub.tab)])
      dimnames(sub.tab) <- list(c(treat_vals(treat), "pooled"), colnames(sub.tab))
      
      ranges <- apply(sub.tab, 1, function(x) mean.abs.dev(x)/sum(x))
      s.d.denom <- rownames(sub.tab)[which.min(ranges)]
    }
    else {
      s.d.denom <- vapply(weights, function(w) {
        for (tv in treat_vals(treat)) {
          if (all_the_same(w[treat == tv]) &&
              !all_the_same(w[treat != tv])) {
            return(tv)
          }
        }
        return(NA_character_)
      }, character(1L))
      s.d.denom[is.na(s.d.denom)] <- "pooled"
    }
  }
  if (is_not_null(weights) && length(s.d.denom) == 1) s.d.denom <- rep.int(s.d.denom, NCOL(weights))
  
  if (s.d.denom.specified && bad.s.d.denom && (!estimand.specified || bad.estimand)) {
    attr(s.d.denom, "note") <- paste0("Warning: 's.d.denom' should be one of ", word_list(unique(c(unique.treats, allowable.s.d.denoms)), "or", quotes = 2), 
                                      ".\n         Using ", word_list(s.d.denom, quotes = 2), " instead.")
  }
  else if (estimand.specified && bad.estimand) {
    attr(s.d.denom, "note") <- paste0("Warning: 'estimand' should be one of ", word_list(c("ATT", "ATC", "ATE"), "or", quotes = 2), 
                                      ". Ignoring 'estimand'.")
  }
  else if ((check.focal || check.weights) && any(s.d.denom %nin% treat_vals(treat))) {
    attr(s.d.denom, "note") <- paste0("Note: 's.d.denom' not specified; assuming ", 
                                      if (all_the_same(s.d.denom)) s.d.denom[1] 
                                      else word_list(paste0("\"", vapply(s.d.denom, function(s) {
                                        if (s %in% treat_vals(treat) && all(treat_vals(treat) %in% c("0", "1"))) names(treat_names(treat))[treat_names(treat) == names(treat_vals(treat))[treat_vals(treat) == s]] else s
                                      }, character(1L)), "\" for ", names(weights))), ".")
  }
  
  if (is_not_null(weights) && length(s.d.denom) != NCOL(weights)) {
    stop(paste0("Valid inputs to 's.d.denom' or 'estimand' must have length 1 or equal to the number of valid sets of weights, which is ", NCOL(weights), "."), call. = FALSE)
  }
  
  if (!quietly && is_not_null(attr(s.d.denom, "note"))) message(attr(s.d.denom, "note"))
  
  return(s.d.denom)
}
get.s.d.denom.cont <- function(s.d.denom, weights = NULL, subclass = NULL, quietly = FALSE) {
  bad.s.d.denom <- FALSE
  s.d.denom.specified <- !missing(s.d.denom) && is_not_null(s.d.denom)
  
  if (is_not_null(subclass)) {
    s.d.denom <- "all"
  }
  else if (s.d.denom.specified) {
    allowable.s.d.denoms <- c("all", "weighted")
    
    if (length(s.d.denom) > 1 && length(s.d.denom) != NCOL(weights)) {
      stop("'s.d.denom' must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
    }
    
    s.d.denom <- match_arg(s.d.denom, unique(allowable.s.d.denoms), 
                           several.ok = TRUE)
  }
  else {
    s.d.denom <- "all"
  }
  
  if (is_not_null(weights) && NCOL(weights) > 1 && length(s.d.denom) == 1) s.d.denom <- rep.int(s.d.denom, NCOL(weights))
  
  if (!quietly) {
    if (s.d.denom.specified && bad.s.d.denom) {
      message(paste0("Warning: 's.d.denom' should be ", word_list(unique(allowable.s.d.denoms), "or", quotes = 2), 
                     ".\n         Using ", word_list(s.d.denom, quotes = 2), " instead."))
    }
  }
  
  if (is_not_null(weights) && length(s.d.denom) != NCOL(weights)) {
    stop("Valid inputs to 's.d.denom' or 'estimand' must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
  }
  
  if (is_not_null(weights)) names(s.d.denom) <- names(weights)
  
  return(s.d.denom)
}
get.estimand <- function(estimand = NULL, weights = NULL, subclass = NULL, treat = NULL, focal = NULL, quietly = TRUE) {
  check.weights <- check.focal <- FALSE
  
  if (is_not_null(estimand)) {
    try.estimand <- tryCatch(match_arg(toupper(estimand), c("ATT", "ATC", "ATE"), several.ok = TRUE),
                             error = function(cond) NA_character_)
    if (anyNA(try.estimand)) {
      check.focal <- TRUE
      # bad.estimand <- TRUE
    }
    else {
      if (length(try.estimand) > 1 && length(try.estimand) != NCOL(weights)) {
        stop("'estimand' must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
      }
      else estimand <- try.estimand
    }
  }
  else {
    check.focal <- TRUE
  }
  
  if (check.focal) {
    if (is_not_null(focal)) {
      estimand <- "ATT"
    }
    else check.weights <- TRUE
  }
  if (check.weights) {
    if (is_null(weights) && is_null(subclass)) {
      estimand <- "ATE"
    }
    else if (is_not_null(subclass)) {
      sub.tab <- table(treat, subclass)[treat_vals(treat), ]
      sub.tab <- rbind(sub.tab, table(subclass)[colnames(sub.tab)])
      dimnames(sub.tab) <- list(c(treat_vals(treat), "Total"), colnames(sub.tab))
      
      ranges <- apply(sub.tab, 1, function(x) mean.abs.dev(x)/sum(x))
      min.range <- which.min(ranges)
      if (rownames(sub.tab)[min.range] == treat_vals(treat)[treat_names(treat)["control"]]) estimand <- "ATC"
      else if (rownames(sub.tab)[min.range] == treat_vals(treat)[treat_names(treat)["treated"]]) estimand <- "ATT"
      else estimand <- "ATE"
      
    }
    else {
      estimand <- vapply(weights, function(w) {
        for (tnn in names(treat_names(treat))) {
          if (all_the_same(w[treat == treat_vals(treat)[treat_names(treat)[tnn]]]) &&
              !all_the_same(w[treat != treat_vals(treat)[treat_names(treat)[tnn]]])) {
            return(switch(tnn, "control" = "ATC", "treated" = "ATT", NA_character_))
          }
        }
        return(NA_character_)
      }, character(1L))
      estimand[is.na(estimand)] <- "ATE"
      
    }
  }
  if (is_not_null(weights) && length(estimand) == 1) estimand <- rep.int(estimand, ncol(weights))
  
  if (!quietly && (check.focal || check.weights)) {
    message("Note: 'estimand' not specified; assuming ", ifelse(all_the_same(toupper(estimand)), toupper(estimand[1]), word_list(paste("\"", toupper(estimand), "\" for ", names(weights)))), ".")
  }
  
  if (is_not_null(weights) && length(estimand) != ncol(weights)) {
    stop("Valid inputs to 'estimand' must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
  }
  return(estimand)
}
compute_s.d.denom <- function(mat, treat, s.d.denom = "pooled", s.weights = NULL, bin.vars = NULL, subset = NULL, weighted.weights = NULL, to.sd = rep(TRUE, ncol(mat)), na.rm = TRUE) {
  denoms <- setNames(rep(1, ncol(mat)), colnames(mat))
  if (is.character(s.d.denom) && length(s.d.denom) == 1L) {
    if (is_null(bin.vars)) {
      bin.vars <- rep(FALSE, ncol(mat))
      bin.vars[to.sd] <- is_binary_col(mat[subset, to.sd,drop = FALSE])
    }
    else if (!is.atomic(bin.vars) || length(bin.vars) != ncol(mat) ||
             anyNA(as.logical(bin.vars))) {
      stop("'bin.vars' must be a logical vector with length equal to the number of columns of 'mat'.")
    }
    
    possibly.supplied <- c("mat", "treat", "weighted.weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
      stop(paste(word_list(possibly.supplied[supplied], quotes = 1), "must have the same number of units."))
    }
    
    if (lengths["weighted.weights"] == 0) weighted.weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("'subset' must be a logical vector.")
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    cont.treat <- get.treat.type(treat) == "continuous"
    
    if (cont.treat) {
      unique.treats <- NULL
      s.d.denom <- get.s.d.denom.cont(as.character(s.d.denom), weights = weighted.weights[subset])
    }
    else {
      unique.treats <- if (is_(treat, "processed.treat") && all(subset)) as.character(treat_vals(treat)) else as.character(unique(treat[subset]))
      s.d.denom <- get.s.d.denom(as.character(s.d.denom), weights = weighted.weights[subset], treat = treat[subset])
      if (s.d.denom %in% c("treated", "control")) s.d.denom <- treat_vals(treat)[treat_names(treat)[s.d.denom]]
      treat <- as.character(treat)
    }
    
    if (s.d.denom %in% unique.treats)
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(col.w.v(mat[treat == s.d.denom, , drop = FALSE],
                     w = s.weights[treat == s.d.denom],
                     bin.vars = bin.vars, na.rm = na.rm))
      }
    
    else if (s.d.denom == "pooled")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(Reduce("+", lapply(unique.treats,
                                function(t) col.w.v(mat[treat == t, , drop = FALSE],
                                                    w = s.weights[treat == t],
                                                    bin.vars = bin.vars, na.rm = na.rm))) / length(unique.treats))
      }
    else if (s.d.denom == "all")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm))
      }
    else if (s.d.denom == "weighted")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(col.w.v(mat, w = weighted.weights * s.weights, bin.vars = bin.vars, na.rm = na.rm))
      }
    else if (s.d.denom == "hedges")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        (1 - 3/(4*length(treat) - 9))^-1 * sqrt(Reduce("+", lapply(unique.treats,
                                                                   function(t) (sum(treat == t) - 1) * col.w.v(mat[treat == t, , drop = FALSE],
                                                                                                               w = s.weights[treat == t],
                                                                                                               bin.vars = bin.vars, na.rm = na.rm))) / (length(treat) - 2))
      }
    else stop("s.d.denom is not an allowed value.")
    
    denoms[to.sd] <- denom.fun(mat = mat[, to.sd, drop = FALSE], treat = treat, s.weights = s.weights,
                               weighted.weights = weighted.weights, bin.vars = bin.vars[to.sd],
                               unique.treats = unique.treats, na.rm = na.rm)
    
    if (any(zero_sds <- !is.finite(denoms[to.sd]) | check_if_zero(denoms[to.sd]))) {
      denoms[to.sd][zero_sds] <- sqrt(col.w.v(mat[, to.sd, drop = FALSE][, zero_sds, drop = FALSE],
                                              w = s.weights,
                                              bin.vars = bin.vars[to.sd][zero_sds], na.rm = na.rm))
    }
    
    if (cont.treat) {
      treat.sd <- denom.fun(mat = treat, s.weights = s.weights,
                            weighted.weights = weighted.weights, bin.vars = FALSE,
                            na.rm = na.rm)
      denoms[to.sd] <- denoms[to.sd]*treat.sd
    }
  }
  else {
    if (is.numeric(s.d.denom)) {
      if (is_not_null(names(s.d.denom)) && any(colnames(mat) %in% names(s.d.denom))) {
        denoms[colnames(mat)[colnames(mat) %in% names(s.d.denom)]] <- s.d.denom[names(s.d.denom)[names(s.d.denom) %in% colnames(mat)]]
      }
      else if (length(s.d.denom) == sum(to.sd)) {
        denoms[to.sd] <- s.d.denom
      }
      else if (length(s.d.denom) == ncol(mat)) {
        denoms[] <- s.d.denom
      }
      else {
        stop("'s.d.denom' must be an allowable value or a numeric vector of with length equal to the number of columns of 'mat'. See ?col_w_smd for allowable values.")
      }
    }
    else {
      stop("'s.d.denom' must be an allowable value or a numeric vector of with length equal to the number of columns of 'mat'. See ?col_w_smd for allowable values.")
    }
  }
  return(denoms)
}
assign.X.class <- function(X) {
  X <- clear_null(X)
  
  if (is_not_null(X[["treat"]]) && !has.treat.type(X[["treat"]])) X[["treat"]] <- assign.treat.type(X[["treat"]])
  
  if (is_not_null(X[["subclass"]])) {
    if (get.treat.type(X[["treat"]]) == "binary") X.class <- "subclass.binary"
    else if (get.treat.type(X[["treat"]]) == "continuous") X.class <- "subclass.cont"
    else stop("Multi-category treatments are not currently compatible with subclasses.", call. = FALSE)
  }
  else if (is_not_null(X[["cluster"]]) && nlevels(X[["cluster"]]) > 1) X.class <- "cluster"
  else if (is_not_null(X[["covs.list"]])) X.class <- "msm"
  else if (get.treat.type(X[["treat"]]) == "multinomial") X.class <- "multi"
  else if (is_not_null(X[["imp"]]) && nlevels(X[["imp"]]) > 1) X.class <- "imp"
  else if (get.treat.type(X[["treat"]]) == "binary") X.class <- "binary"
  else if (get.treat.type(X[["treat"]]) == "continuous") X.class <- "cont"
  else probably.a.bug()
  
  class(X) <- X.class
  
  return(X)
}
get_length_X <- function(X) {
  if (is_not_null(X[["treat"]])) length(X[["treat"]])
  else if (is_not_null(X[["covs"]])) nrow(X[["covs"]])
  else if (is_not_null(X[["treat.list"]])) length(X[["treat.list"]][[1]])
  else if (is_not_null(X[["covs.list"]])) nrow(X[["covs.list"]][[1]])
  else stop("Couldn't determine length of X components.")
}
subsettable <- function() {
  c("covs",
    "treat",
    "weights",
    "distance",
    "addl",
    "cluster",
    "imp",
    "s.weights",
    "discarded",
    "subclass",
    "covs.list",
    "treat.list",
    "distance.list",
    "addl.list")
}
subset_X <- function(X, subset = NULL) {
  if (is_not_null(subset) && any(names(X) %in% subsettable())) {
    n <- get_length_X(X)
    if (is.logical(subset)) {
      if (length(subset) != n) stop("'subset' must have the same length as the other entries.")
      if (!any(subset)) stop("All 'subset' set to FALSE.", call. = FALSE)
      to_be_subset <- !all(subset)
      subset <- which(subset)
    }
    else if (is.numeric(subset)) {
      if (max(subset) > n) stop("Subset indices cannot be higher than the length of the other entries.")
      to_be_subset <- TRUE
    }
    else stop("'subset' must be logical or numeric.")
    
    if (to_be_subset) {
      subset_X_internal <- function(x, subset) {
        attrs <- attributes(x)
        attrs_to_subset <- names(attrs)[vapply(attrs, function(a) all(len(a) == n), logical(1L))]
        if (is_not_null(attrs_to_subset)) {
          subsetted_attrs <- lapply(attrs[attrs_to_subset],
                                    subset_X_internal, subset = subset)
        }
        
        if (is_null(x)) out <- x
        else if (is_(x, "processed.treat")) out <- x[subset]
        else if ((is.matrix(x) || is.data.frame(x))) out <- x[subset, , drop = FALSE]
        else if (is.factor(x)) out <- factor(x[subset], nmax = nlevels(x))
        else if (is.atomic(x)) out <- x[subset]
        else if (is.list(x)) out <- lapply(x, subset_X_internal, subset = subset)
        else out <- x
        
        if (is_not_null(attrs)) {
          if (all(c("treat_names", "treat_vals", "treat.type") %in% 
                  names(attrs))) {
            out <- process_treat(out)
          }
          else {
            for (i in names(attrs)[names(attrs) %nin% names(attributes(out))]) {
              if (i %in% attrs_to_subset) attr(out, i) <- subsetted_attrs[[i]]
              else attr(out, i) <- attrs[[i]]
            }
          }
        }
        
        out
      }
      X[names(X) %in% subsettable()] <- lapply(X[names(X) %in% subsettable()], subset_X_internal, subset)
    }
  }
  
  return(X)
}
imp.complete <- function(data) {
  if (!is_(data, "mids")) stop("'data' not of class 'mids'")
  
  single.complete <- function(data, where = NULL, imp, ell) {
    if (is_null(where)) where <- is.na(data)
    idx <- seq_len(ncol(data))[which(colSums(where) > 0)]
    for (j in idx) {
      if (is_null(imp[[j]])) data[where[, j], j] <- NA
      else data[where[, j], j] <- imp[[j]][, ell]
    }
    return(data)
  }
  
  m <- as.integer(data$m)
  idx <- seq_len(m)
  
  mylist <- lapply(idx, function(i) single.complete(data$data, data$where, data$imp, i))
  
  cmp <- data.frame(.imp = rep(idx, each = nrow(data$data)), 
                    .id = rep.int(seq_len(nrow(data$data)), length(idx)), 
                    do.call("rbind", mylist))
  
  if (is.integer(attr(data$data, "row.names"))) 
    row.names(cmp) <- seq_len(nrow(cmp))
  else row.names(cmp) <- as.character(seq_len(nrow(cmp)))
  
  return(cmp)
}
length_imp_process <- function(vectors = NULL, data.frames = NULL, lists = NULL, imp = NULL, data = NULL, original.call.to = NULL, env = sys.frame(-1)) {
  #Processes imp and assigns it to imp in env
  #Processes vectors, data.frames, and lists for length
  #If object correspond to one imputation, assigns expanded object in env 
  
  all.objects <- c(vectors, data.frames, lists)
  ensure.equal.lengths <- TRUE
  problematic <- rlang::rep_named(all.objects, FALSE)
  lengths <- setNames(c(vapply(vectors, 
                               function(x) {len(get0(x, envir = env, inherits = FALSE))
                               }, numeric(1L)), 
                        vapply(data.frames, 
                               function(x) {len(get0(x, envir = env, inherits = FALSE))
                               }, numeric(1L)),
                        vapply(lists, function(x) {
                          if (is_null(get0(x, envir = env, inherits = FALSE))) 0 
                          else max(vapply(get(x, envir = env, inherits = FALSE), len, numeric(1L)))
                        }, numeric(1L))), 
                      all.objects)
  
  #Process imp further
  if (is_not_null(imp)) {
    
    imp.lengths <- vapply(levels(imp), function(i) sum(imp == i), numeric(1L))
    
    if (all_the_same(imp.lengths)) { #all the same
      unsorted.imp <- is.unsorted(imp)
      for (i in vectors) {
        if (lengths[i] > 0 && lengths[i] != length(imp)) { 
          if (!all_the_same(imp.lengths)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
          if (lengths[i] == imp.lengths[1]) {
            i_obj <- get(i, envir = env, inherits = FALSE)
            new_i <- i_obj[rep(seq_along(i_obj), length(imp.lengths))]
            if (unsorted.imp) {for (i_ in levels(imp)) new_i[imp == i_] <- i_obj}
            assign(i, new_i, pos = env)
          }
          else {
            problematic[i] <- TRUE
          }
        }
      }
      for (i in data.frames) {
        if (lengths[i] > 0 && lengths[i] != length(imp)) {
          if (!all_the_same(imp.lengths)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
          if (lengths[i] == imp.lengths[1]) {
            i_obj <- get(i, envir = env, inherits = FALSE)
            new_i <- i_obj[rep(seq_len(nrow(i_obj)), length(imp.lengths)), , drop = FALSE]
            if (unsorted.imp) {for (i_ in levels(imp)) new_i[imp == i_,] <- i_obj}
            assign(i, new_i, pos = env)
          }
          else {
            problematic[i] <- TRUE
          }
        }
      }
      for (i in lists) {
        if (lengths[i] > 0 && lengths[i] != length(imp)) {
          if (!all_the_same(imp.lengths)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
          if (lengths[i] == imp.lengths[1]) {
            assign(i, lapply(get(i, envir = env, inherits = FALSE), function(j) {
              if (is_(j, "factor")) {
                newj <- j[rep(seq_along(j), length(imp.lengths))]
                if (unsorted.imp) {for (i_ in levels(imp)) newj[imp == i_] <- j}
                return(newj)
              }
              else if (is_(j, c("data.frame", "matrix"))) {
                newj <- j[rep(seq_len(nrow(j)), length(imp.lengths)), , drop = FALSE]
                if (unsorted.imp) {for (i_ in levels(imp)) newj[imp == i_,] <- j}
                return(newj)
              }
              else {
                stop(paste0("'", i, "' can only contain vectors or data frames."), call. = FALSE)
              }
            }), pos = env)
          }
          else {
            problematic[i] <- TRUE
          }
        }
      }
    }
    else {
      problematic <- lengths > 0 & lengths != length(imp)
    }
    if (any(problematic)) {
      stop(paste0(word_list(names(problematic)[problematic], quotes = 1), " must have the same number of observations as 'imp'."), call. = FALSE)
    }
    else ensure.equal.lengths <- FALSE
    
    assign("imp", imp, pos = env)
  }
  
  #Ensure all input lengths are the same.
  anchor <- {
    if ("treat" %in% all.objects) "treat"
    else if ("treat.list" %in% all.objects) "treat.list"
    else all.objects[which(lengths[all.objects] != 0)[1]]
  }
  if (ensure.equal.lengths) {
    problematic[lengths %nin% c(0, lengths[anchor])] <- TRUE
  }
  if (any(problematic)) {
    if (is_not_null(original.call.to)) anchor <- paste("in the original call to", original.call.to)
    
    stop(paste0(word_list(names(problematic)[problematic], quotes = 1), " must have the same number of observations as ", anchor, "."), call. = FALSE)
  }
}
process_stats <- function(stats = NULL, treat) {
  if (is_(treat, "list")) {
    stats.list <- lapply(treat, function(x) process_stats(stats, x))
    if (all_the_same(vapply(stats.list, attr, character(1L), "type")))
      type <- attr(stats.list[[1]], "type")
    else type <- NULL
    stats <- unique(unlist(stats.list))
    attr(stats, "type") <- type
  }
  else {
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type %in% c("binary", "multinomial")) {
      if (is_null(stats)) stats <- getOption("cobalt_stats", "mean.diffs")
      stats <- unique(match_arg(stats, all_STATS("bin"), several.ok = TRUE))
      attr(stats, "type") <- "bin"
    }
    else if (treat.type %in% c("continuous")) {
      if (is_null(stats)) stats <- getOption("cobalt_stats", "correlations")
      stats <- unique(match_arg(stats, all_STATS("cont"), several.ok = TRUE))
      attr(stats, "type") <- "cont"
    }
  }
  return(stats)
}
process_thresholds <- function(thresholds, stats) {
  if (is_not_null(thresholds)) {
    if (is.list(thresholds)) thresholds <- unlist(thresholds)
    if (!all(is.na(thresholds)) && !is.numeric(thresholds)) {
      stop("'thresholds' must be numeric.", call. = FALSE)
    }
    
    if (rlang::is_named(thresholds)) {
      names(thresholds) <- stats[pmatch(names(thresholds), stats, duplicates.ok = TRUE)]
      thresholds <- thresholds[!is.na(names(thresholds))]
    }
    else {
      names(thresholds) <- stats[1:min(length(stats), length(thresholds))]
    }
    
    thresholds[names(thresholds)] <- as.numeric(thresholds)
  }
  thresholds <- as.list(na.rem(thresholds))
  return(thresholds)
}
process_subset <- function(subset, n) {
  if (!is_(subset, c("logical", "numeric"))) {
    stop("The argument to 'subset' must be a logical or numeric vector.", call. = FALSE)
  }
  else if (is.numeric(subset)) {
    if (any(abs(subset) > n)) stop("Numeric values for 'subset' cannot be larger than the number of units.", call. = FALSE)
    subset <- subset[!is.na(subset) & subset != 0]
    if (any(subset < 0) && any(subset > 0)) stop("Positive and negative indices cannot be mixed with 'subset'.")
    if (any(abs(subset) > n)) stop("If 'subset' is numeric, none of its values can exceed the number of units.")
    logical.subset <- rep(any(subset < 0), n)
    logical.subset[abs(subset)] <- !logical.subset[abs(subset)]
    subset <- logical.subset
  }
  if (anyNA(subset)) {
    warning("NAs were present in 'subset'. Treating them like FALSE.", call. = FALSE)
    subset[is.na(subset)] <- FALSE
  }
  return(subset)
}
process_focal <- function(focal, treat) {
  if (is.numeric(focal)) {
    if (can_str2num(treat) && focal %in% str2num(treat)) {focal <- as.character(focal)}
    else if (focal <= length(treat_vals(treat))) focal <- treat_vals(treat)[focal]
    else 
      stop(paste0("'focal' was specified as ", focal, 
                  ", but there are only ", length(treat_vals(treat)), " treatment groups."), call. = FALSE)
  }
  else {
    if (focal %nin% treat_vals(treat)) 
      stop(paste0("The name specified to 'focal' is not the name of any treatment group."), call. = FALSE)
  }
  return(focal)
}
process_weights <- function(obj = NULL, A = NULL, treat = NULL, covs = NULL, method = character(0), addl.data = list(), ...) {
  if (is_not_null(obj)) {
    weights <- get.w(obj, treat = treat, ...)
    
    if (is_not_null(weights)) {
      if (is_(weights, c("data.frame", "matrix"))) {
        weights <- data.frame(weights)
      }
      else {
        weights <- setNames(data.frame(weights), class(obj)[has_method(class(obj), "get.w")][1])
      }
    }
    else {
      weights <- list()
    }
  }
  else weights <- list()
  
  addl.weights <- data.frame.process("weights", A[["weights"]], treat, covs, addl.data = addl.data, ...)
  if (is_not_null(addl.weights)) {
    if (is_null(A[["method"]])) addl.methods <- rep.int("weighting", ncol(addl.weights))
    else if (length(A[["method"]]) == 1) {
      addl.methods <- rep.int(match_arg(A[["method"]], c("weighting", "matching")), ncol(addl.weights))
    }
    else {
      addl.methods <- match_arg(A[["method"]], c("weighting", "matching"), several.ok = TRUE)
      if (length(addl.methods) != ncol(addl.weights)) 
        stop("Valid inputs to 'method' must have length 1 or equal to the number of valid sets of additional weights.", call. = FALSE)
    }
    
    w.names <- c(names(weights), names(addl.weights))
    unique.names <- unique(w.names)
    if (length(unique.names) != length(w.names)) {
      for (i in unique.names) {
        if (sum(w.names == i) > 1) w.names[w.names == i] <- make.unique(c(i, w.names[w.names == i]))[-1]
      } 
    }
    if (is_not_null(weights)) {
      weights <- setNames(cbind(weights, addl.weights), w.names)
      method <- setNames(c(method, addl.methods), w.names)
    }
    else {
      weights <- setNames(addl.weights, w.names)
      method <- setNames(addl.methods, w.names)
    }
  }
  weight.check(weights)
  
  attr(weights, "method") <- method
  return(weights)
}
process_disp <- function(disp = NULL, ...) {
  A <- list(...)
  if (is_not_null(disp)) {
    if (!is.character(disp)) stop("'disp' must be a character vector.")
    disp <- match_arg(disp, acceptable.options()[["disp"]], several.ok = TRUE)
  }
  else disp <- getOption("cobalt_disp")
  
  for (d in c("means", "sds")) {
    if (getOption(paste.("cobalt_disp", d), FALSE)) disp <- unique(c(disp, d))
    if (is_not_null(A[[paste.("disp", d)]])) {
      if (!rlang::is_bool(A[[paste.("disp", d)]])) stop(paste0("'disp.", d, "' must be TRUE or FALSE."), call. = FALSE)
      disp <- unique(c(disp, d[A[[paste.("disp", d)]]]))
      if (A[[paste.("disp", d)]]) disp <- unique(c(disp, d))
      else disp <- unique(disp[disp != d])
    }
  }
  return(disp)
}
process_addl <- function(addl = NULL, datalist = list()) {
  data <- do.call("data.frame", unname(clear_null(datalist)))
  if (is_(addl, "atomic") && 
      (!is_(addl, "character") || is_null(datalist) ||
       length(addl) == nrow(data))) {
    addl <- data.frame(addl = addl)
  }
  else if (is_(addl, "character")) addl <- f.build(addl)
  
  addl_t.c <- use.tc.fd(formula = addl, data = data, covs = addl, 
                        needs.treat = FALSE, needs.covs = FALSE)
  
  return(addl_t.c[["covs"]])
}
process_addl.list <- function(addl.list = NULL, datalist = list(), covs.list = list()) {
  datalist <- clear_null(c(datalist, covs.list))
  
  if (is_(addl.list, "list")) {
    if (length(addl.list) != length(covs.list)) stop("'addl' must have an entry for each time point.", call. = FALSE)
    addl.list.out <- lapply(addl.list, process_addl, datalist = datalist)
  }
  else {
    addl <- process_addl(addl.list, datalist = datalist)
    addl.list.out <- lapply(seq_along(covs.list), function(x) addl)
  }
  return(addl.list.out)
}
process_distance <- function(distance = NULL, datalist = list(), obj.distance = NULL, obj.distance.name = "distance") {
  data <- do.call("data.frame", unname(clear_null(datalist)))
  if (is_not_null(distance) && !is_(distance, c("atomic", "formula", "matrix", "data.frame"))) {
    stop("'distance' must be a formula or variable containing the distance values.", call. = FALSE)
  }
  if (is_(distance, "atomic") && 
      (!is_(distance, "character") || is_null(datalist) ||
       length(distance) == nrow(data))) {
    distance <- data.frame(distance = distance)
  }
  else if (is_(distance, "character")) distance <- f.build(distance)
  
  distance_t.c <- use.tc.fd(formula = distance, data = data, covs = distance, 
                            needs.treat = FALSE, needs.covs = FALSE)
  
  distance <- distance_t.c[["covs"]]
  
  if (is_not_null(obj.distance) && !all(is.na(obj.distance))) {
    obj.distance <- setNames(data.frame(obj.distance), obj.distance.name)
    obj.distance <- get_covs_from_formula(data = obj.distance)
    distance <- co.cbind(if_null_then(distance, NULL), obj.distance)
  }
  
  return(distance)
}
process_distance.list <- function(distance.list = NULL, datalist = list(), covs.list = list(), obj.distance = NULL, obj.distance.name = "distance") {
  datalist <- clear_null(c(datalist, covs.list))
  
  if (is_not_null(obj.distance)) {
    if (!is_(obj.distance, "list")) {
      obj.distance <- lapply(seq_along(covs.list), function(x) obj.distance)
    }
  }
  else obj.distance <- lapply(seq_along(covs.list), function(x) NULL)
  
  if (is_null(distance.list)) {
    distance.list.out <- lapply(seq_along(covs.list), function(x) process_distance(NULL, datalist = datalist, 
                                                                                   obj.distance = obj.distance[[x]], obj.distance.name = obj.distance.name))
  }
  else if (is_(distance.list, "list")) {
    if (length(distance.list) != length(covs.list)) stop("'distance' must have an entry for each time point.", call. = FALSE)
    distance.list.out <- lapply(seq_along(distance.list), function(x) process_distance(distance.list[[x]], datalist = datalist, 
                                                                                       obj.distance = obj.distance[[x]], obj.distance.name = obj.distance.name))
  }
  else {
    distance.list.out <- lapply(seq_along(covs.list), function(x) process_distance(distance.list, datalist = datalist, 
                                                                                   obj.distance = obj.distance[[x]], obj.distance.name = obj.distance.name))
  }
  
  return(distance.list.out)
}

#get.C2
get_ints_from_co.names <- function(co.names) {
  if (is_not_null(co.names)) {
    clear_null(lapply(co.names, function(co) {
      if ("isep" %in% co[["type"]]) {
        co[["type"]] <- c(co[["type"]], "isep")
        which_isep <- which(co[["type"]] == "isep")
        vapply(seq_along(which_isep), function(x) {
          if (x == 1) paste0(co[["component"]][1:(which_isep[1]-1)], collapse = "")
          else paste0(co[["component"]][(which_isep[x-1]+1):(which_isep[x]-1)], collapse = "")
        }, character(1L))
      }
      else NULL
    }))
  }
  else list()
}
get_treat_from_formula <- function(f, data = NULL, treat = NULL) {
  
  if (!rlang::is_formula(f)) stop("'f' must be a formula.", call. = FALSE)
  
  env <- rlang::f_env(f)
  
  f <- update(f, ~ 0)
  
  #Check if data exists
  if (is_not_null(data)) {
    if (is.data.frame(data)) {
      data.specified <- TRUE
    }
    else {
      warning("The argument supplied to 'data' is not a data.frame object. Ignoring 'data'.", call. = FALSE)
      data <- env
      data.specified <- FALSE
    }
  }
  else {
    data <- env
    data.specified <- FALSE
  }
  
  
  tryCatch(tt <- terms(f, data = data),
           error = function(e) {
             stop(conditionMessage(e), call. = FALSE)
           })
  
  if (rlang::is_formula(tt, lhs = TRUE)) {
    resp.vars.mentioned <- as.character(tt)[2]
    resp.vars.failed <- vapply(resp.vars.mentioned, function(v) {
      test <- tryCatch(eval(str2expression(v), data, env), error = function(e) e)
      if (inherits(test, "simpleError")) {
        if (conditionMessage(test) == paste0("object '", v, "' not found")) return(TRUE)
        else stop(conditionMessage(test), call. = FALSE)
      }
      else if (is.function(test)) stop(paste0("invalid type (function) for variable '", v, "'"), call. = FALSE)
      else return(is_null(test))
    }, logical(1L))
    
    if (any(resp.vars.failed)) {
      if (is_null(treat)) stop(paste0("The given response variable, \"", as.character(tt)[2], "\", is not a variable in ", word_list(c("'data'", "the global environment")[c(data.specified, TRUE)], "or"), "."), call. = FALSE)
      tt <- delete.response(tt)
    }
  }
  else resp.vars.failed <- TRUE
  
  if (!all(resp.vars.failed)) {
    treat.name <- resp.vars.mentioned[!resp.vars.failed][1]
    treat <- eval(str2expression(treat.name), data, env)
  }
  else {
    treat.name <- NULL
  }
  attr(treat, "treat.name") <- treat.name
  
  return(treat)
}
get_covs_from_formula <- function(f, data = NULL, factor_sep = "_", int_sep = " * ") {
  
  rebuild_f <- function(ttfactors, tics = FALSE) {
    #Set tics = TRUE if returned formula is used with tmpcovs,
    #i.e., a model.frame. If used with data, leave FALSE
    if (tics) rownames(ttfactors)[!startsWith(rownames(ttfactors), "`")] <- paste0("`", rownames(ttfactors)[!startsWith(rownames(ttfactors), "`")], "`")
    as.formula(paste("~ 0 +", paste(vapply(seq_len(ncol(ttfactors)), 
                                           function(x) paste0(rownames(ttfactors)[ttfactors[,x] > 0], collapse = ":"),
                                           character(1L)), 
                                    collapse = "+")))
  }
  append.ttfactor <- function(ttfactor, term, after) {
    addrow <- matrix(0, nrow = length(term), ncol = ncol(ttfactor), 
                     dimnames = list(term, colnames(ttfactor)))
    addcol <- matrix(0, nrow = nrow(ttfactor) + length(term), ncol = length(term),
                     dimnames = list(c(rownames(ttfactor), term), term))
    addcol[-seq_len(nrow(ttfactor)), ] <- diag(length(term))
    
    ttfactor <- rbind(ttfactor, addrow)
    if (after == 0) {
      return(cbind(addcol, ttfactor))
    }
    else if (after == ncol(ttfactor)) {
      return(cbind(ttfactor, addcol))
    }
    else {
      return(cbind(ttfactor[,seq_len(after), drop = FALSE], 
                   addcol, 
                   ttfactor[,-seq_len(after), drop = FALSE]))
    }
  }
  
  
  
  #Check if data exists
  data.specified <- FALSE
  if (is_not_null(data)) {
    if (is.matrix(data)) data <- as.data.frame.matrix(data)
    
    if (is.data.frame(data)) {
      data.specified <- TRUE
    }
    else {
      warning("The argument supplied to 'data' is not a data.frame object. Ignoring 'data'.", call. = FALSE)
    }
  }
  
  if (missing(f) && data.specified) f <- f.build(data)
  else if (!rlang::is_formula(f)) stop("'f' must be a formula.")
  
  env <- rlang::f_env(f)
  if (!data.specified) data <- env
  
  rlang::f_lhs(f) <- NULL
  
  tryCatch(tt <- terms(f, data = data),
           error = function(e) {
             if (conditionMessage(e) == "'.' in formula and no 'data' argument") {
               stop("'.' is not allowed in formulas.", call. = FALSE)
             }
             else stop(conditionMessage(e), call. = FALSE)
           })
  
  #Process RHS 
  tt.covs <- delete.response(tt)
  attr(tt.covs, "intercept") <- 0
  
  ttfactors <- attr(tt.covs, "factors")
  ttvars <- setNames(vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1], rownames(ttfactors))
  
  rhs.df.type <- setNames(vapply(ttvars, function(v) {
    if (is_(try(eval(str2expression(paste0("`", v, "`")), data, env), silent = TRUE),
            c("data.frame", "matrix", "rms"))) "lit"
    else if (is_(try(eval(str2expression(v), data, env), silent = TRUE),
                 c("data.frame", "matrix", "rms"))) "exp"
    else "not.a.df"
  }, character(1L)), ttvars)
  
  rhs.df <- setNames(rhs.df.type != "not.a.df", ttvars)
  
  while (any(rhs.df)) {
    term_is_interaction <- apply(ttfactors, 2, function(x) sum(x != 0) > 1)
    if (any(vapply(seq_along(ttvars)[rhs.df], function(x) any(ttfactors[x,] != 0 & term_is_interaction), logical(1L)))) {
      stop("Interactions with data.frames are not allowed in the input formula.", call. = FALSE)
    }
    addl.dfs <- setNames(lapply(ttvars[rhs.df], function(v) {
      if (rhs.df.type[v] == "lit") df <- eval(str2expression(paste0("`", v, "`")), data, env)
      else df <- eval(str2expression(v), data, env)
      if (is_(df, "rms")) {
        df <- setNames(as.data.frame.matrix(as.matrix(df)), colnames(df))
        return(df)
      }
      else if (can_str2num(colnames(df))) colnames(df) <- paste(v, colnames(df), sep = "_")
      return(as.data.frame(df))
    }),
    ttvars[rhs.df])
    
    for (i in colnames(ttfactors)[colnames(ttfactors) %in% names(ttvars)[rhs.df]]) {
      for (j in seq_len(ncol(addl.dfs[[ttvars[i]]]))) {
        if (names(addl.dfs[[ttvars[i]]])[j] %in% c(ttvars[!rhs.df], unlist(lapply(addl.dfs[seq_len(which(names(addl.dfs) == ttvars[i]) - 1)], names)))) {
          names(addl.dfs[[ttvars[i]]])[j] <- paste0(ttvars[i], "_", names(addl.dfs[[ttvars[i]]])[j])
        }
      }
      ind <- which(colnames(ttfactors) == i)
      ttfactors <- append.ttfactor(ttfactors,
                                   paste0("`", names(addl.dfs[[ttvars[i]]]), "`"),
                                   ind)[,-ind, drop = FALSE]
    }
    
    if (data.specified) {
      data <- do.call("cbind", unname(c(addl.dfs, list(data))))
    } else {
      data <- do.call("cbind", unname(addl.dfs))
      data.specified <- TRUE
    }
    
    new.form <- rebuild_f(ttfactors)
    tt.covs <- terms(new.form, data = data)
    
    ttfactors <- attr(tt.covs, "factors")
    ttvars <- setNames(vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1], rownames(ttfactors))
    
    rhs.df.type <- setNames(vapply(ttvars, function(v) {
      if (is_(try(eval(str2expression(paste0("`", v, "`")), data, env), silent = TRUE),
              c("data.frame", "matrix", "rms"))) "lit"
      else if (is_(try(eval(str2expression(v), data, env), silent = TRUE),
                   c("data.frame", "matrix", "rms"))) "exp"
      else "not.a.df"
    }, character(1L)), ttvars)
    
    rhs.df <- setNames(rhs.df.type != "not.a.df", ttvars)
  }
  
  #Check to make sure variables are valid
  original_ttvars <- rownames(ttfactors)
  for (i in seq_along(rownames(ttfactors))) {
    #Check if evaluable
    #If not, check if evaluable after changing to literal using ``
    #If not, stop()
    #If eventually evaluable, check if function
    
    evaled.var <- try(eval(str2expression(rownames(ttfactors)[i]), data, env), silent = TRUE)
    if (null_or_error(evaled.var)) {
      evaled.var <- try(eval(str2expression(paste0("`", rownames(ttfactors)[i], "`")), data, env), silent = TRUE)
      if (null_or_error(evaled.var)) {
        stop(conditionMessage(attr(evaled.var, "condition")), call. = FALSE)
      }
      else {
        rownames(ttfactors)[i] <- paste0("`", rownames(ttfactors)[i], "`")
      }
    }
    
    evaled.var <- tryCatch(eval(str2expression(rownames(ttfactors)[i]), data, env), error = function(e) stop(conditionMessage(e), call. = FALSE))
    if (is.function(evaled.var)) {
      stop(paste0("invalid type (function) for variable '", rownames(ttfactors)[i], "'"), call. = FALSE)
    }
  }
  
  if (!identical(original_ttvars, rownames(ttfactors))) {
    new.form <- rebuild_f(ttfactors)
    tt.covs <- terms(new.form, data = data)
    
    ttfactors <- attr(tt.covs, "factors")
    ttvars <- vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1]
  }
  
  tryCatch({tmpcovs <- stats::model.frame(tt.covs, data, na.action = "na.pass")},
           error = function(e) {stop(conditionMessage(e), call. = FALSE)})
  
  for (i in ttvars) {
    if (is_binary(tmpcovs[[i]])) tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = 2)
    else {
      if (is.character(tmpcovs[[i]])) tmpcovs[[i]] <- factor(tmpcovs[[i]])
      if (is.factor(tmpcovs[[i]])) {
        if (nlevels(tmpcovs[[i]]) == 1) tmpcovs[[i]] <- factor(tmpcovs[[i]], levels = c(paste0(".",levels(tmpcovs[[i]])), levels(tmpcovs[[i]])))
        else tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = nlevels(tmpcovs[[i]]))
      }
    }
  }
  
  #Process NAs: make NA variables
  if (anyNA(tmpcovs)) {
    has_NA <- anyNA_col(tmpcovs)
    for (i in rev(colnames(tmpcovs)[has_NA])) {
      #Find which of ttlabels i first appears, and put `i: <NA>` after it
      for (x in seq_along(colnames(ttfactors))) {
        if (i %in% c(colnames(ttfactors)[x], all.vars(str2expression(colnames(ttfactors)[x])))) {
          ind <- x
          break
        }
      }
      ttfactors <- append.ttfactor(ttfactors, 
                                   paste0("`", i, ":<NA>`"),
                                   ind)
      
      tmpcovs[[paste0(i, ":<NA>")]] <- as.numeric(is.na(tmpcovs[[i]]))
    }
    new.form <- rebuild_f(ttfactors, tics = TRUE)
    
    tt.covs <- terms(new.form, data = tmpcovs)
    ttfactors <- attr(tt.covs, "factors")
    ttvars <- vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1]
    
    na_vars <- paste0(colnames(tmpcovs)[has_NA], ":<NA>")
    
    tryCatch({tmpcovs <- stats::model.frame(tt.covs, tmpcovs, na.action = "na.pass")},
             error = function(e) {stop(conditionMessage(e), call. = FALSE)})
    
    for (i in ttvars[ttvars %nin% na_vars]) {
      if (is_binary(tmpcovs[[i]])) tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = 2)
      else {
        if (is.character(tmpcovs[[i]])) tmpcovs[[i]] <- factor(tmpcovs[[i]])
        if (is.factor(tmpcovs[[i]])) {
          if (nlevels(tmpcovs[[i]]) == 1) tmpcovs[[i]] <- factor(tmpcovs[[i]], levels = c(paste0(".",levels(tmpcovs[[i]])), levels(tmpcovs[[i]])))
          else tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = nlevels(tmpcovs[[i]]))
        }
      }
    }
  }
  else {
    na_vars <- character(0)
  }
  
  #Re-check ttfactors
  original_ttvars <- rownames(ttfactors)
  for (i in seq_along(rownames(ttfactors))) {
    #Check if evaluable in tmpcovs
    #If not, check if evaluable i tmpcovs after changing to literal using ``
    #If not, stop() (shouldn't occur)
    
    evaled.var <- try(eval(str2expression(rownames(ttfactors)[i]), tmpcovs), silent = TRUE)
    if (null_or_error(evaled.var)) {
      evaled.var <- try(eval(str2expression(paste0("`", rownames(ttfactors)[i], "`")), tmpcovs), silent = TRUE)
      if (null_or_error(evaled.var)) {
        stop(conditionMessage(attr(evaled.var, "condition")), call. = FALSE)
      }
      else {
        rownames(ttfactors)[i] <- paste0("`", rownames(ttfactors)[i], "`")
      }
    }
  }
  
  if (!identical(original_ttvars, rownames(ttfactors))) {
    new.form <- rebuild_f(ttfactors)
    tt.covs <- terms(new.form, data = data)
    
    ttfactors <- attr(tt.covs, "factors")
    ttvars <- vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1]
  }
  
  tmpcovs <- model.frame(tt.covs, data = tmpcovs, drop.unused.levels = TRUE,
                         na.action = "na.pass")
  
  #Check for infinite values
  covs.with.inf <- vapply(tmpcovs, function(x) is.numeric(x) && any(!is.na(x) & !is.finite(x)), logical(1L))
  if (any(covs.with.inf)) {
    s <- if (sum(covs.with.inf) == 1) c("", "s") else c("s", "")
    stop(paste0("The variable", s[1], " ", word_list(names(tmpcovs)[covs.with.inf], quotes = 1), 
                " contain", s[2], " non-finite values, which are not allowed."), call. = FALSE)
  }
  
  attr(tt.covs, "intercept") <- 1 #Add intercept to correctly process single-level factors
  mm <- model.matrix(tt.covs, data = tmpcovs, 
                     contrasts.arg = lapply(Filter(is.factor, tmpcovs),
                                            function(x) contrasts(x, contrasts = all_the_same(x))))
  
  rownames(ttfactors) <- trim_string(rownames(ttfactors), "`")
  
  mmassign <- attr(mm, "assign")[-1]
  mmassign2 <- setNames(factor(mmassign, levels = sort(unique(mmassign), na.last = TRUE),
                               labels = colnames(ttfactors)), colnames(mm)[-1])
  
  vars_in_each_term <- setNames(lapply(colnames(ttfactors), function(x) {
    rownames(ttfactors)[ttfactors[,x] != 0]
  }), colnames(ttfactors))
  all_factor_levels <- lapply(vars_in_each_term, function(v) {
    do.call("expand.grid", c(clear_null(setNames(lapply(v, function(fa) colnames(attr(mm, "contrasts")[[fa]])), v)),
                             list(stringsAsFactors = FALSE)))
  })
  expanded <- setNames(lapply(seq_along(mmassign2), function(x) {
    terms <- vars_in_each_term[[mmassign2[x]]]
    k <- sum(seq_along(mmassign2) <= x & mmassign2 == mmassign2[x])
    setNames(lapply(terms, function(t) {
      if (t %in% names(all_factor_levels[[mmassign2[x]]])) {
        all_factor_levels[[mmassign2[x]]][[t]][k]
      }
      else character(0)
    }), terms)
  }), names(mmassign2))
  
  #component types: base, fsep, isep, power, na, level
  co.names <- lapply(expanded, function(x) {
    Reduce(function(x1, x2 = NULL) {
      list(component = c(x1[["component"]], int_sep, x2[["component"]]),
           type = c(x1[["type"]], "isep", x2[["type"]]))
      
    }, lapply(seq_along(x), function(i) {
      base <- gsub("`", "", names(x)[i], fixed = TRUE)
      if (base %in% na_vars) {
        base <- substr(base, 1, nchar(base) - 5)
        out <- list(component = c(base, ":<NA>"),
                    type = c("base", "na"))
      }
      else {
        out <- list(component = base,
                    type = "base")
        if (is_not_null(x[[i]])) {
          out[["component"]] <- c(out[["component"]], factor_sep, x[[i]])
          out[["type"]] <- c(out[["type"]], "fsep", "level")
        }
      }
      out
    }))
  })
  
  names(co.names) <- vapply(co.names, function(x) paste0(x[["component"]], collapse = ""), character(1L))
  covs <- clear_attr(mm)[,-1, drop = FALSE] #drop the intercept
  
  attr(co.names, "seps") <- c(factor = factor_sep, int = int_sep)
  attr(covs, "co.names") <- co.names
  
  colnames(covs) <- names(co.names)
  return(covs)
}
get.C2 <- function(covs, int = FALSE, poly = 1, addl = NULL, distance = NULL, treat = NULL, cluster = NULL, drop = TRUE, ...) {
  #gets C data.frame, which contains all variables for which balance is to be assessed. Used in balance.table.
  A <- list(...)
  if (!rlang::is_bool(A[["center"]])) A[["center"]] <- getOption("cobalt_center", default = FALSE)
  
  C <- covs; rm(covs)
  
  co.names <- attr(C, "co.names")
  seps <- attr(co.names, "seps")
  
  if (is_not_null(addl)) {
    addl.co.names <- attr(addl, "co.names")
    
    same.name <- names(addl.co.names) %in% names(co.names)
    addl <- addl[,!same.name, drop = FALSE]
    addl.co.names[same.name] <- NULL
    
    #Remove variables in addl that are redundant with C
    if (drop && getOption("cobalt_remove_perfect_col", max(ncol(addl), ncol(C)) <= 900)) {
      redundant.var.indices <- find_perfect_col(addl, C)
      if (is_not_null(redundant.var.indices)) {
        addl <- addl[,-redundant.var.indices, drop = FALSE]
        addl.co.names[redundant.var.indices] <- NULL
      }
    }
    
    C <- cbind(C, addl)
    co.names <- c(co.names, addl.co.names)
  } 
  
  #Drop single_value or colinear with cluster
  if (drop) {
    test.treat <- is_not_null(treat) && get.treat.type(treat) != "continuous"
    test.cluster <- is_not_null(cluster) && !all_the_same(cluster, na.rm = FALSE)
    drop_vars <- vapply(seq_len(ncol(C)), 
                        function(i) {
                          # if (all_the_same(C[,i], na.rm = FALSE)) return(TRUE)
                          # else 
                          if (anyNA(C[,i])) return(FALSE)
                          else if (test.treat && equivalent.factors2(C[,i], treat)) return(TRUE)
                          # else if (test.cluster && equivalent.factors2(C[,i], cluster)) return(TRUE) #Note: doesn't work with multi-cat cluster vars due to splitting
                          else return(FALSE)
                        }, logical(1L))

    if (all(drop_vars)) stop("There are no variables for which to display balance.", call. = FALSE)
    C <- C[,!drop_vars, drop = FALSE]
    co.names[drop_vars] <- NULL
  }
  
  C_list <- list(C = C)
  co_list <- list(C = co.names)
  rm(C, co.names)
  
  #Process int and poly
  if (length(int) != 1L || !(rlang::is_bool(int) || is.numeric(int)) || !is.finite(int)) {
    stop("'int' must be TRUE, FALSE, or a numeric value of length 1.", call. = FALSE)
  }
  if (!rlang::is_bool(int) && (int < 0 || !check_if_int(int))) {
    stop("'int' must be TRUE, FALSE, or a numeric (integer) value greater than 1.", call. = FALSE)
  }
  int <- as.integer(round(int))
  
  if (length(poly) != 1L || !is.finite(poly) || !is.numeric(poly)) {
    stop("'poly' must be a numeric value of length 1.", call. = FALSE)
  }
  if (poly < 0 || !check_if_int(poly)) {
    stop("'poly' must be a numeric (integer) value greater than 1.", call. = FALSE)
  }
  poly <- round(poly)
  
  if (int || (poly > 1)) {
    if (int) { 
      #Prevent duplicate var names with `sep`s
      nsep <- 1
      # repeat {
      #   all.possible.names <- outer(colnames(C), colnames(C), paste, sep = paste0(rep.int(A[["int_sep"]], nsep), collapse = ""))
      #   if (!any(colnames(C) %in% all.possible.names)) break
      #   else nsep <- nsep + 1
      # }
      
      if (poly < int) poly <- int
      
      int <- TRUE
    }
    
    #Exclude NA and ints from interactions and poly
    exclude <- vapply(co_list[["C"]], function(x) any(c("na", "isep") %in% x[["type"]]), logical(1L))
    
    new <- int.poly.f2(C_list[["C"]], ex = exclude, int = int, poly = poly, center = A[["center"]], 
                       sep = rep.int(seps["int"], nsep), co.names = co_list[["C"]])
    
    C_list[["int.poly"]] <- new
    co_list[["int.poly"]] <- attr(new, "co.names")
    names(co_list[["int.poly"]]) <- vapply(co_list[["int.poly"]], 
                                           function(x) paste0(x[["component"]], collapse = ""), character(1L))
    # C <- cbind(C, new)
    # co.names <- c(co.names, attr(new, "co.names"))
    
    # names(co.names) <- vapply(co.names, function(x) paste0(x[["component"]], collapse = ""), character(1L))
    # colnames(C) <- names(co.names)
  }
  
  #Drop 0 category of 0/1 variables and rename 1 category
  if (drop) {
    drop_0_1 <- rep(NA, length(co_list[["C"]]))
    for (i in seq_along(co_list[["C"]])) {
      if (is.na(drop_0_1[i])) {
        if ("isep" %nin% co_list[["C"]][[i]][["type"]] && "fsep" %in% co_list[["C"]][[i]][["type"]]) {
          which_are_buddies <- which(vapply(co_list[["C"]], function(j) "isep" %nin% j[["type"]] && 
                                              "fsep" %in% j[["type"]] &&
                                              j[["component"]][j[["type"]] == "base"][1] == co_list[["C"]][[i]][["component"]][co_list[["C"]][[i]][["type"]] == "base"][1], 
                                            logical(1L)))
          buddies <- co_list[["C"]][which_are_buddies]
          if (length(buddies) <= 2) {
            buddy_is_0 <- vapply(buddies, function(x) x[["component"]][x[["type"]] == "level"] %in% c("0", "FALSE"), logical(1L))
            buddy_is_1 <- vapply(buddies, function(x) x[["component"]][x[["type"]] == "level"] %in% c("1", "TRUE"), logical(1L))
            if (all(buddy_is_0 | buddy_is_1)) {
              drop_0_1[which_are_buddies[buddy_is_0]] <- TRUE
              drop_0_1[which_are_buddies[buddy_is_1]] <- FALSE
              
              buddy_1 <- which_are_buddies[buddy_is_1]
              co_list[["C"]][[buddy_1]][["component"]] <- co_list[["C"]][[buddy_1]][["component"]][co_list[["C"]][[buddy_1]][["type"]] == "base"][1]
              co_list[["C"]][[buddy_1]][["type"]] <- "base"
            }
            else drop_0_1[which_are_buddies] <- c(TRUE, FALSE)
          }
          else drop_0_1[which_are_buddies] <- FALSE
        }
        else drop_0_1[i] <- FALSE
      }
    }
    
    C_list[["C"]] <- C_list[["C"]][,!drop_0_1, drop = FALSE]
    co_list[["C"]][drop_0_1] <- NULL
  }
  
  names(co_list[["C"]]) <- vapply(co_list[["C"]], function(x) paste0(x[["component"]], collapse = ""), character(1L))
  
  if (is_not_null(distance)) {
    if (anyNA(distance, recursive = TRUE)) stop("Missing values are not allowed in the distance measure.", call. = FALSE)
    
    distance.co.names <- attr(distance, "co.names")
    
    same.name <- names(distance.co.names) %in% do.call("c", lapply(co_list, names))
    distance <- distance[,!same.name, drop = FALSE]
    distance.co.names[same.name] <- NULL
    
    unique.distance.names <- unique(names(distance.co.names))
    distance <- distance[,unique.distance.names, drop = FALSE]
    distance.co.names <- distance.co.names[unique.distance.names]
    
    C_list[["distance"]] <- distance
    co_list[["distance"]] <- distance.co.names
    
  }
  
  C_list <- clear_null(C_list)
  co_list <- clear_null(co_list)
  
  #Remove duplicate & redundant variables
  if (drop) {
    for (x in setdiff(names(C_list), "distance")) {
      #Remove self-redundant variables
      if (getOption("cobalt_remove_perfect_col", ncol(C_list[[x]]) <= 900)) {
        redundant.var.indices <- find_perfect_col(C_list[[x]])
        if (is_not_null(redundant.var.indices)) {
          C_list[[x]] <- C_list[[x]][,-redundant.var.indices, drop = FALSE]
          co_list[[x]][redundant.var.indices] <- NULL
        }
      }
      if (x != "C") {
        #Remove variables in C that have same name as other variables
        if (any(dups <- names(co_list[["C"]]) %in% co_list[[x]])) {
          C_list[["C"]] <- C_list[["C"]][,!dups, drop = FALSE]
          co_list[["C"]][dups] <- NULL 
        }
        #Remove variables in C that are redundant with current piece
        if (getOption("cobalt_remove_perfect_col", max(ncol(C_list[[x]]), ncol(C_list[["C"]])) <= 900)) {
          redundant.var.indices <- find_perfect_col(C_list[["C"]], C_list[[x]])
          if (is_not_null(redundant.var.indices)) {
            C_list[["C"]] <- C_list[["C"]][,-redundant.var.indices, drop = FALSE]
            co_list[["C"]][redundant.var.indices] <- NULL
          }
        }
      }
    }
  }
  
  C <- cbind(C_list[["distance"]], C_list[["C"]], C_list[["int.poly"]])
  co.names <- c(co_list[["distance"]], co_list[["C"]], co_list[["int.poly"]])
  
  colnames(C) <- names(co.names) <- vapply(co.names, function(x) paste0(x[["component"]], collapse = ""), character(1L))
  
  attr(co.names, "seps") <- seps
  
  attr(C, "co.names") <- co.names
  
  attr(C, "missing.ind") <- colnames(C)[vapply(co.names, function(x) "na" %in% x[["type"]], logical(1L))]
  if ("distance" %in% names(C_list)) attr(C, "distance.names") <- names(co_list[["distance"]])
  
  return(C)
  
}
int.poly.f2 <- function(mat, ex = NULL, int = FALSE, poly = 1, center = FALSE, sep = " * ", co.names = NULL) {
  #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
  #Only to be used in base.bal.tab; for general use see int.poly()
  #mat=matrix input
  #ex=names of variables to exclude in interactions and polynomials; a subset of df
  #int=whether to include interactions or not; currently only 2-way are supported
  #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
  #nunder=number of underscores between variables
  
  cn <- is_not_null(co.names)
  if (is_null(ex)) ex <- rep(FALSE, ncol(mat))
  d <- mat
  
  binary.vars <- is_binary_col(d)
  interaction.vars <- if (cn) vapply(colnames(d), function(x) "isep" %in% co.names[[x]][["type"]], logical(1L)) else rep(FALSE, ncol(d))
  
  if (center) {
    d[,!binary.vars] <- center(d[, !binary.vars, drop = FALSE])
  }
  nd <- NCOL(d)
  
  if (poly > 1) {
    poly_terms <- poly_co.names <- make_list(poly-1)
    no.poly <- binary.vars | interaction.vars | ex
    npol <- nd - sum(no.poly)
    if (npol > 0) {
      for (i in 2:poly) {
        poly_terms[[i - 1]] <- apply(d[, !no.poly, drop = FALSE], 2, function(x) x^i)
        if (cn) poly_co.names[[i - 1]] <- lapply(colnames(d)[!no.poly], function(x) list(component = c(co.names[[x]][["component"]], num_to_superscript(i)), 
                                                                                         type = c(co.names[[x]][["type"]], "power")))
        else poly_co.names[[i - 1]] <- paste0(colnames(d)[!no.poly], num_to_superscript(i))
      }
    }
  }
  else poly_terms <- poly_co.names <- list()
  
  if (int && nd > 1) {
    int_terms <- int_co.names <- make_list(1)
    ints_to_make <- combn(colnames(d)[!ex], 2, simplify = FALSE)
    
    #Don't make ints out of ints that already exist
    # ints_that_already_exist <- get_ints_from_co.names(co.names[interaction.vars])
    # ints_to_make[vapply(ints_to_make, function(x) {
    #   any(vapply(ints_that_already_exist, function(y) {
    #     identical(sort(x), sort(y))
    #   }, logical(1L)))
    # }, logical(1L))] <- NULL
    
    #Don't make ints out of multiple members of the same categorical variable
    ints_to_make[vapply(ints_to_make, function(x) {
      "fsep" %in% co.names[[x[1]]][["type"]] && 
        "fsep" %in% co.names[[x[2]]][["type"]] &&
        identical(co.names[[x[1]]][["component"]][co.names[[x[1]]][["type"]] == "base"],
                  co.names[[x[2]]][["component"]][co.names[[x[2]]][["type"]] == "base"])
    }, logical(1L))] <- NULL
    
    int_terms[[1]] <- do.call("cbind", lapply(ints_to_make, function(i) d[,i[1]]*d[,i[2]]))
    
    if (cn) int_co.names[[1]] <- lapply(ints_to_make, function(x) list(component = c(co.names[[x[1]]][["component"]], sep, co.names[[x[2]]][["component"]]),
                                                                       type = c(co.names[[x[1]]][["type"]], "isep", co.names[[x[2]]][["type"]])))
    else int_co.names[[1]] <- vapply(ints_to_make, paste, character(1L), collapse = sep)
  }
  else int_terms <- int_co.names <- list()
  
  out <- do.call("cbind", c(poly_terms, int_terms))
  out_co.names <- c(do.call("c", poly_co.names), do.call("c", int_co.names))
  
  if (cn) {
    names(out_co.names) <- vapply(out_co.names, 
                                  function(x) paste0(x[["component"]], collapse = ""), character(1))
    colnames(out) <- names(out_co.names)
  }
  else {
    colnames(out) <- unlist(out_co.names)
  }
  
  #Remove single values
  single_value <- apply(out, 2, all_the_same)
  out <- out[, !single_value, drop = FALSE]
  if (cn && is_not_null(out)) attr(out, "co.names") <- out_co.names[!single_value]
  
  return(out)
}
co.cbind <- function(..., deparse.level = 1) {
  args <- clear_null(list(...))
  if (length(args) <= 1) return(args[[1]])
  
  co.names.list <- lapply(args, attr, "co.names")
  
  seps <- attr(co.names.list[[1]], "seps")
  
  out <- do.call("cbind", args)
  
  attr(out, "co.names") <- do.call("c", co.names.list)
  attr(attr(out, "co.names"), "seps") <- seps
  colnames(out) <- names(attr(out, "co.names")) <- vapply(attr(out, "co.names"), function(x) paste0(x[["component"]], collapse = ""), character(1L))
  
  out
}
co.rbind <- function(..., deparse.level = 1) {
  args <- clear_null(list(...))
  if (length(args) <= 1) return(args[[1]])
  
  co.names <- attr(args[[1]], "co.names")
  
  out <- do.call("rbind", args)
  
  attr(out, "co.names") <- co.names
  colnames(out) <- names(attr(out, "co.names")) <- vapply(attr(out, "co.names"), function(x) paste0(x[["component"]], collapse = ""), character(1L))
  
  out
}

df_clean <- function(df) {
  if (!is.list(df)) df <- as.data.frame(df)
  #If a data.frame has matrix columns, clean those into vector columns and rename
  if (any(vapply(df, function(x) is_not_null(dim(x)), logical(1L)))) {
    rn <- rownames(df)
    
    df <- do.call("cbind", lapply(seq_along(df), function(x) {
      if (is_null(dim(df[[x]]))) {
        setNames(data.frame(df[[x]]), names(df)[x])
      }
      else {
        if (is_(df[[x]], "rms")) {
          setNames(as.data.frame.matrix(as.matrix(df[[x]])), colnames(df[[x]]))
        }
        else if (can_str2num(colnames(df[[x]]))) {
          setNames(as.data.frame(df[[x]]), paste(names(df)[x], colnames(df[[x]]), sep = "_"))
        }
        else {
          setNames(as.data.frame(df[[x]]), colnames(df[[x]]))
        }
      }
    }))
    
    rownames(df) <- rn
  }
  
  df
}
get.types <- function(C) {
  vapply(colnames(C), function(x) {
    if (any(attr(C, "distance.names") == x)) "Distance"
    else if (is_binary(C[,x]))  "Binary"
    else "Contin."
  }, character(1))
}
find_perfect_col <- function(C1, C2 = NULL, fun = stats::cor) {
  #Finds indices of redundant vars in C1.
  C1.no.miss <- C1[,colnames(C1) %nin% attr(C1, "missing.ind"), drop = FALSE]
  if (is_null(C2)) {
    use <- if (anyNA(C1)) "pairwise.complete.obs" else "everything"
    suppressWarnings(C.cor <- fun(C1.no.miss, use = use))
    s <- !lower.tri(C.cor, diag=TRUE) & !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
  }
  else {
    C2.no.miss <- C2[,colnames(C2) %nin% attr(C2, "missing.ind"), drop = FALSE]
    use <- if (anyNA(C1) || anyNA(C2)) "pairwise.complete.obs" else "everything"
    suppressWarnings(C.cor <- fun(C2.no.miss, y = C1.no.miss, use = use))
    s <- !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
  }
  
  redundant.var.indices <- which(colSums(s) > 0)
  return(redundant.var.indices)
}

#base.bal.tab
check_if_zero_weights <- function(weights.df, treat = NULL) {
  #Checks if all weights are zero in each treat group for each set of weights
  if (is_not_null(treat)) {
    w.t.mat <- expand.grid(colnames(weights.df), treat_vals(treat), stringsAsFactors = TRUE)
    names(w.t.mat) <- c("weight_names", "treat_vals")
    
    if (NROW(w.t.mat) > 0) {
      problems <- vapply(seq_len(NROW(w.t.mat)), function(x) all(check_if_zero(weights.df[treat == w.t.mat[x, "treat_vals"], w.t.mat[x, "weight_names"]])), logical(1L))
      if (any(problems)) {
        prob.w.t.mat <- droplevels(w.t.mat[problems,])
        if (NCOL(weights.df) == 1) {
          error <- paste0("All weights are zero when treat is ", word_list(prob.w.t.mat[, "treat_vals"], "or"), ".")
        }
        else {
          errors <- vapply(levels(prob.w.t.mat[,"weight_names"]), function(i) {
            paste0("\"", i, "\" weights are zero when treat is ", word_list(prob.w.t.mat[prob.w.t.mat[,"weight_names"] == i, "treat_vals"], "or"))
          }, character(1L))
          errors <- paste(c("All", rep("all", length(errors)-1)), errors)
          error <- paste0(word_list(errors, "and"), ".")
        }
        stop(error, call. = FALSE)
      }
    }
  }
  else {
    if (length(colnames(weights.df)) > 0) {
      problems <- vapply(colnames(weights.df), function(wn) all(check_if_zero(weights.df[, wn])), logical(1L))
      if (any(problems)) {
        prob.wts <- colnames(weights.df)[problems]
        if (NCOL(weights.df) == 1) {
          error <- paste0("All weights are zero.")
        }
        else {
          errors <- vapply(prob.wts, function(i) {
            paste0("\"", i, "\" weights are zero")
          }, character(1L))
          errors <- paste(c("All", rep("all", length(errors)-1)), errors)
          error <- paste0(word_list(errors, "and"), ".")
        }
        stop(error, call. = FALSE)
      }
    }
  }
}
baltal <- function(threshold) {
  #threshold: vector of threshold values (i.e., "Balanced"/"Not Balanced")
  threshnames <- names(table(threshold))
  balstring <- threshnames[nchar(threshnames) > 0][1]
  thresh.val <- substring(balstring, 1 + regexpr("[><]", balstring), nchar(balstring))
  b <- data.frame(count=c(sum(threshold==paste0("Balanced, <", thresh.val)), 
                          sum(threshold==paste0("Not Balanced, >", thresh.val))))
  rownames(b) <- c(paste0("Balanced, <", thresh.val), paste0("Not Balanced, >", thresh.val))
  return(b)
}
max.imbal <- function(balance.table, col.name, thresh.col.name, abs_stat) {
  balance.table.clean <- balance.table[balance.table$Type != "Distance" & is.finite(balance.table[, col.name]),]
  maxed <- balance.table.clean[which.max(abs_stat(balance.table.clean[, col.name])), match(c(col.name, thresh.col.name), names(balance.table.clean))]
  maxed <- data.frame(Variable = rownames(maxed), maxed)
  return(maxed)
}
threshold.summary <- function(compute, thresholds, no.adj, balance.table, weight.names = NULL, agg.fun = NULL) {
  out <- do.call("c", lapply(compute, function(s) make_list(paste.(c("Balanced", "Max.Imbalance"), s))))
  
  for (s in compute) {
    
    if (is_not_null(thresholds[[s]])) {
      
      if (no.adj) {
        out[[paste.("Balanced", s)]] <- baltal(balance.table[[paste.(STATS[[s]]$Threshold, "Un")]])
        out[[paste.("Max.Imbalance", s)]] <- max.imbal(balance.table[balance.table[["Type"]]!="Distance", , drop = FALSE], 
                                                       col.name = if (is_null(agg.fun)) paste.(STATS[[s]]$bal.tab_column_prefix, "Un")
                                                       else paste.(firstup(agg.fun), STATS[[s]]$bal.tab_column_prefix, "Un"), 
                                                       thresh.col.name = paste.(STATS[[s]]$Threshold, "Un"), 
                                                       abs_stat = STATS[[s]]$abs)
      }
      else if (length(weight.names) == 1) {
        out[[paste.("Balanced", s)]] <- baltal(balance.table[[STATS[[s]]$Threshold]])
        out[[paste.("Max.Imbalance", s)]] <- max.imbal(balance.table[balance.table[["Type"]]!="Distance", , drop = FALSE], 
                                                       col.name = if (is_null(agg.fun)) paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")
                                                       else paste.(firstup(agg.fun), STATS[[s]]$bal.tab_column_prefix, "Adj"), 
                                                       thresh.col.name = STATS[[s]]$Threshold, 
                                                       abs_stat = STATS[[s]]$abs)
      }
      else if (length(weight.names) > 1) {
        out[[paste.("Balanced", s)]] <- setNames(do.call("cbind", lapply(weight.names, function(x) baltal(balance.table[[paste.(STATS[[s]]$Threshold, x)]]))),
                                                 weight.names)
        out[[paste.("Max.Imbalance", s)]] <- cbind(Weights = weight.names,
                                                   do.call("rbind", lapply(weight.names, function(x) setNames(max.imbal(balance.table[balance.table[["Type"]]!="Distance", , drop = FALSE], 
                                                                                                                        col.name = if (is_null(agg.fun)) paste.(STATS[[s]]$bal.tab_column_prefix, x)
                                                                                                                        else paste.(firstup(agg.fun), STATS[[s]]$bal.tab_column_prefix, x),  
                                                                                                                        thresh.col.name = paste.(STATS[[s]]$Threshold, x), 
                                                                                                                        abs_stat = STATS[[s]]$abs),
                                                                                                              c("Variable", 
                                                                                                                STATS[[s]]$bal.tab_column_prefix, 
                                                                                                                STATS[[s]]$Threshold)))),
                                                   stringsAsFactors = FALSE)
      }
    }
    else {
      out[[paste.("Balanced", s)]] <- NULL
      out[[paste.("Max.Imbalance", s)]] <- NULL
    }
  }
  
  out
}

balance.table <- function(C, type, weights = NULL, treat, continuous, binary, s.d.denom, 
                          thresholds = list(), un = FALSE, disp = NULL, stats = NULL, 
                          s.weights = rep(1, length(treat)), abs = FALSE, no.adj = FALSE, 
                          var_types = NULL, s.d.denom.list = NULL, quick = TRUE, ...) {
  #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
  if (no.adj) weight.names <- "Adj"
  else weight.names <- names(weights)
  
  if (is_not_null(s.d.denom.list)) names(s.d.denom.list) <- weight.names
  if (is_not_null(s.d.denom)) names(s.d.denom) <- weight.names
  
  disp <- c(disp, all_STATS(type)[all_STATS(type) %in% stats])
  compute <- if (quick) disp else c("means", "sds", all_STATS(type)[all_STATS(type) %in% stats])
  
  #B=Balance frame
  Bnames <- c("Type", 
              expand.grid_string(c(if (type == "bin") expand.grid_string(c("M"["means" %in% compute], "SD"["sds" %in% compute]), c("0", "1"), collapse = ".")
                                   else if (type == "cont") c("M"["means" %in% compute], "SD"["sds" %in% compute]), 
                                   unlist(lapply(compute[compute %in% all_STATS(type)[!get_from_STATS("adj_only")[get_from_STATS("type") == type]]], function(s) {
                                     c(STATS[[s]]$bal.tab_column_prefix,
                                       if (no.adj && is_not_null(thresholds[[s]])) STATS[[s]]$Threshold)
                                   }))
              ),
              "Un", collapse = "."),
              expand.grid_string(c(if (type == "bin") expand.grid_string(c("M"["means" %in% compute], "SD"["sds" %in% compute]), c("0", "1"), collapse = ".")
                                   else if (type == "cont") c("M"["means" %in% compute], "SD"["sds" %in% compute]), 
                                   unlist(lapply(compute[compute %in% all_STATS(type)], function(s) {
                                     c(STATS[[s]]$bal.tab_column_prefix,
                                       if (!no.adj && is_not_null(thresholds[[s]])) STATS[[s]]$Threshold)
                                   }))
              ),
              weight.names, collapse = "."))
  B <- make_df(Bnames, NCOL(C))
  rownames(B) <- colnames(C)
  
  #Set var type (binary/continuous)
  B[["Type"]] <- if_null_then(var_types, get.types(C))
  bin.vars <- B[["Type"]] == "Binary"
  
  #Means for each group
  if ("means" %in% compute) {
    if (type == "bin") {
      tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
      
      if (un || !quick) {
        for (t in c("0", "1")) {
          B[[paste.("M", t, "Un")]] <- col_w_mean(C, weights = NULL, s.weights = s.weights, subset = treat == tn01[t])
        }
      }
      
      if (!no.adj && (!quick || "means" %in% disp)) {
        for (i in weight.names) {
          for (t in c("0", "1")) {
            B[[paste.("M", t, i)]] <- col_w_mean(C, weights = weights[[i]], s.weights = s.weights, subset = treat == tn01[t])
          }
        }
      }
    }
    else if (type == "cont") {
      if (un || !quick) {
        B[["M.Un"]] <- col_w_mean(C, weights = NULL, s.weights = s.weights)
      }
      if (!no.adj && (!quick || "means" %in% disp)) {
        for (i in weight.names) {
          B[[paste.("M", i)]] <- col_w_mean(C, weights = weights[[i]], s.weights = s.weights)
        }
      }
    }
  }
  
  #SDs for each group
  binary <- match_arg(binary, c("raw", "std"))
  if ("sds" %in% compute) {
    sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else !bin.vars
    if (type == "bin") {
      tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
      
      if (un || !quick) {
        for (t in c("0", "1")) {
          sds <- rep(NA_real_, NCOL(C))
          if (any(sd.computable)) {
            sds[sd.computable] <- col_w_sd(C[, sd.computable,drop = FALSE], weights = NULL, s.weights = s.weights,
                                           bin.vars = bin.vars[sd.computable], subset = treat == tn01[t])
          }
          B[[paste.("SD", t, "Un")]] <- sds
        }
      }
      
      if (!no.adj && (!quick || "sds" %in% disp)) {
        for (i in weight.names) {
          for (t in c("0", "1")) {
            sds <- rep(NA_real_, NCOL(C))
            if (any(sd.computable)) {
              sds[sd.computable] <- col_w_sd(C[, sd.computable,drop = FALSE], weights = weights[[i]], s.weights = s.weights,
                                             bin.vars = bin.vars[sd.computable], subset = treat == tn01[t])
            }
            B[[paste.("SD", t, i)]] <- sds
          }
        }
      }
    }
    else if (type == "cont") {
      if (un || !quick) {
        sds <- rep(NA_real_, NCOL(C))
        if (any(sd.computable)) {
          sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE],
                                         weights = NULL, s.weights = s.weights,
                                         bin.vars = bin.vars[sd.computable])
        }
        B[["SD.Un"]] <- sds
      }
      if (!no.adj && (!quick || "sds" %in% disp)) {
        for (i in weight.names) {
          sds <- rep(NA_real_, NCOL(C))
          if (any(sd.computable)) {
            sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE],
                                           weights = weights[[i]], s.weights = s.weights,
                                           bin.vars = bin.vars[sd.computable])
          }
          B[[paste.("SD", i)]] <- sds
        }
      }
    }
    if (all(vapply(B[startsWith(names(B), "SD.")], function(x) all(!is.finite(x)), logical(1L)))) {
      disp <- disp[disp != "sds"]
    }
  }
  
  for (s in all_STATS(type)) {
    if (s %in% compute) {
      if (!get_from_STATS("adj_only")[s] && (!quick || un)) {
        B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]] <- STATS[[s]]$fun(C, treat = treat, weights = NULL, 
                                                                              std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                              s.d.denom = if_null_then(s.d.denom.list[[1]], s.d.denom[1]),
                                                                              abs = abs, s.weights = s.weights, bin.vars = bin.vars,
                                                                              weighted.weights = weights[[1]], ...)
      }
      
      if (!no.adj && (!quick || s %in%  disp)) {
        for (i in weight.names) {
          B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]] <- STATS[[s]]$fun(C, treat = treat, weights = weights[[i]],
                                                                             std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                             s.d.denom = if_null_then(s.d.denom.list[[i]], s.d.denom[i]),
                                                                             abs = abs, s.weights = s.weights, bin.vars = bin.vars, ...)
        }
      }
      
      
      if (all(vapply(intersect(names(B), paste.(STATS[[s]]$bal.tab_column_prefix, c("Un", weight.names))), 
                     function(x) all(!is.finite(B[[x]])), logical(1L)))) {
        disp <- disp[disp != s] 
        thresholds[[s]] <- NULL
      }
      
      if (is_not_null(thresholds[[s]])) {
        if (!get_from_STATS("adj_only")[s] && no.adj) {
          B[[paste.(STATS[[s]]$Threshold, "Un")]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]]), "", 
                                                             STATS[[s]]$abs(B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                             paste0("Not Balanced, >", round(thresholds[[s]], 3)))
        }
        else if (!no.adj) {
          for (i in weight.names) {
            B[[paste.(STATS[[s]]$Threshold, i)]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]]), "", 
                                                            STATS[[s]]$abs(B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                            paste0("Not Balanced, >", round(thresholds[[s]], 3)))
          }
        }
      }
      if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
    }
  }
  
  attr(B, "thresholds") <- thresholds
  attr(B, "disp") <- disp
  attr(B, "compute") <- compute
  
  return(B)
}

samplesize <- function(treat, type, weights = NULL, subclass = NULL, s.weights = NULL, method=c("matching", "weighting", "subclassification"), discarded = NULL) {
  #Computes sample size info. for unadjusted and adjusted samples.
  # method is what method the weights are to be used for. 
  # method="subclassification" is for subclass sample sizes only.
  
  if (is_null(s.weights)) s.weights <- rep(1, length(treat))
  if (is_null(discarded)) discarded <- rep(FALSE, length(treat))
  
  if (type == "bin") {
    if (length(method) == 1 && method == "subclassification") {
      if (is_null(subclass)) stop("'subclass' must be a vector of subclasses.")
      
      nn <- make_df(c(levels(subclass), "All"), c(treat_names(treat), "Total"))
      
      nn[["All"]] <- c(vapply(treat_vals(treat), function(tn) sum(treat==tn), numeric(1L)), length(treat))
      
      matched <- !is.na(subclass)
      for (k in levels(subclass)) {
        qt <- treat[matched & subclass == k]
        nn[[k]] <- c(vapply(treat_vals(treat), function(tn) sum(qt==tn), numeric(1L)), length(qt))
      }
      for (tnn in names(treat_names(treat))) {
        small.subclass <- nn[treat_names(treat)[tnn], levels(subclass)] <= 1
        if (any(small.subclass))
          warning(paste0("Not enough ", tnn, " units in ", ngettext(sum(small.subclass), "subclass ", "subclasses "), 
                         word_list(levels(subclass)[small.subclass]), "."), call. = FALSE)
      }
      attr(nn, "tag") <- "Sample sizes by subclass"
    }
    else {
      if (is_null(weights)) {
        
        nn <- make_df(treat_names(treat), "All")
        nn["All", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
        if (nunique.gt(s.weights, 2) || !any(s.weights==1) || any(s.weights %nin% c(0,1))) {
          attr(nn, "ss.type") <- c("ess")
        }
        else {
          attr(nn, "ss.type") <- c("ss")
        }
        
      }
      else if (NCOL(weights) == 1) {
        if (method=="matching") {
          nn <- make_df(treat_names(treat), c("All", "Matched (ESS)", "Matched (Unweighted)", "Unmatched", "Discarded"))
          nn["All", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn), numeric(1L))
          nn["Matched (ESS)", ] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, 1]), numeric(1L))
          nn["Matched (Unweighted)", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & weights[,1] > 0), numeric(1L))
          nn["Unmatched", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & weights[,1]==0 & !discarded), numeric(1L))
          nn["Discarded", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & weights[,1]==0 & discarded), numeric(1L))
          
          attr(nn, "ss.type") <- rep("ss", NROW(nn))
          
          if (!any(discarded)) {
            attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
            nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
          }
        }
        else if (method == "weighting") {
          nn <- make_df(treat_names(treat), c("Unadjusted", "Adjusted", "Discarded"))
          nn["Unadjusted", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
          nn["Adjusted", ] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, 1]*s.weights[treat==tn]), numeric(1L))
          nn["Discarded", ] <- vapply(treat_vals(treat), function(tn) sum(treat==tn & discarded), numeric(1L))
          attr(nn, "ss.type") <- c("ss", "ess", "ss")
          
          if (!any(discarded)) {
            attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
            nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
          }
        }
      }
      else {
        nn <- make_df(treat_names(treat), c("All", names(weights)))
        nn["All", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat==tn]), numeric(1L))
        for (i in seq_len(NCOL(weights))) {
          if (method[i] == "matching") {
            nn[1+i,] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, i]), numeric(1L))
          }
          else if (method[i] == "weighting") {
            nn[1+i,] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat==tn, i]*s.weights[treat==tn]), numeric(1L))
          }
          
        }
        attr(nn, "ss.type") <- c("ss", rep("ess", length(method)))
      }
      if (length(attr(nn, "ss.type")) > 1 && all(attr(nn, "ss.type")[-1] == "ess")) {
        attr(nn, "tag") <- "Effective sample sizes"
      }
      else attr(nn, "tag") <- "Sample sizes"
    }
  }
  else if (type == "cont") {
    if (length(method) == 1 && method == "subclassification") {
      if (is_null(subclass)) stop("'subclass' must be a vector of subclasses.")
      
      nn <- make_df(c(levels(subclass), "All"), c("Total"))
      
      nn[, "All"] <- length(treat)
      
      matched <- !is.na(subclass)
      for (k in levels(subclass)) {
        nn[[k]] <- sum(matched & subclass == k)
      }
      small.subclass <- nn[, levels(subclass)] <= 1
      if (any(small.subclass))
        warning(paste0("Not enough units in ", ngettext(sum(small.subclass), "subclass ", "subclasses "), 
                       word_list(levels(subclass)[small.subclass]), "."), call. = FALSE)
      attr(nn, "tag") <- "Sample sizes by subclass"
    }
    else {
      if (is_null(weights)) {
        nn <- make_df("Total", "All")
        nn["All", ] <- ESS(s.weights)
        if (nunique.gt(s.weights, 2) || !any(s.weights==1) || any(s.weights %nin% c(0,1))) {
          attr(nn, "ss.type") <- c("ess")
        }
        else {
          attr(nn, "ss.type") <- c("ss")
        }
        
      }
      else if (NCOL(weights) == 1) {
        if (method=="matching") {
          nn <- make_df("Total", c("All", "Matched (ESS)", "Matched (Unweighted)", "Unmatched", "Discarded"))
          nn["All", ] <- length(treat)
          nn["Matched (ESS)", ] <- ESS(weights[, 1])
          nn["Matched (Unweighted)", ] <- sum(weights[,1] > 0 & !discarded)
          nn["Unmatched", ] <- sum(weights[,1]==0 & !discarded)
          nn["Discarded", ] <- sum(discarded)
          
          attr(nn, "ss.type") <- rep("ss", NROW(nn))
          
          if (!any(discarded)) {
            attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
            nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
          }
        }
        else if (method == "weighting") {
          nn <- make_df("Total", c("Unadjusted", "Adjusted", "Discarded"))
          nn["Unadjusted", ] <- ESS(s.weights)
          nn["Adjusted", ] <- ESS(weights[!discarded, 1]*s.weights[!discarded])
          nn["Discarded", ] <- sum(discarded)
          attr(nn, "ss.type") <- c("ss", "ess", "ss")
          
          if (!any(discarded)) {
            attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
            nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
          }
        }
      }
      else {
        nn <- make_df("Total", c("All", names(weights)))
        nn["All", ] <- ESS(s.weights)
        for (i in seq_len(NCOL(weights))) {
          if (method[i] == "matching") {
            nn[1+i,] <- ESS(weights[!discarded, i])
          }
          else if (method[i] == "weighting") {
            nn[1+i,] <- ESS(weights[!discarded, i]*s.weights[!discarded])
          }
          
        }
        attr(nn, "ss.type") <- c("ss", rep("ess", length(method)))
        
      }
      if (length(attr(nn, "ss.type")) > 1 && all(attr(nn, "ss.type")[-1] == "ess")) {
        attr(nn, "tag") <- "Effective sample sizes"
      }
      else attr(nn, "tag") <- "Sample sizes"
    }
  }
  return(nn)
}

balance.summary <- function(bal.tab.list, agg.funs, include.times = FALSE) {
  type <- attr(bal.tab.list[[1]], "print.options")[["type"]]
  disp <- attr(bal.tab.list[[1]], "print.options")[["disp"]]
  compute <- attr(bal.tab.list[[1]], "print.options")[["compute"]]
  thresholds <- attr(bal.tab.list[[1]], "print.options")[["thresholds"]]
  quick <- attr(bal.tab.list[[1]], "print.options")[["quick"]]
  weight.names <- if_null_then(attr(bal.tab.list[[1]], "print.options")[["weight.names"]], "Adj")
  abs <- attr(bal.tab.list[[1]], "print.options")[["abs"]]
  no.adj <- attr(bal.tab.list[[1]], "print.options")[["nweights"]] == 0
  
  balance.list <- clear_null(lapply(bal.tab.list, function(x) x[["Balance"]]))
  
  Brownames <- unique(unlist(lapply(balance.list, rownames), use.names = FALSE))
  
  Agg.Funs <- firstup(if (quick) agg.funs else c("min", "mean", "max"))
  Agg.Funs.Given <- firstup(agg.funs)
  
  if (length(Agg.Funs) == 1 && Agg.Funs == "Max") abs <- TRUE
  
  Bcolnames <- c("Times", "Type", 
                 paste.(c(unlist(lapply(compute[compute %in% all_STATS(type)[!get_from_STATS("adj_only")[get_from_STATS("type") == type]]], function(s) {
                   c(paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix),
                     if (no.adj && is_not_null(thresholds[[s]]) && length(Agg.Funs.Given) == 1) STATS[[s]]$Threshold)
                 }))), "Un"),
                 unlist(lapply(weight.names, function(wn) {
                   paste.(c(unlist(lapply(compute[compute %in% all_STATS(type)], function(s) {
                     c(paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix),
                       if (!no.adj && is_not_null(thresholds[[s]]) && length(Agg.Funs.Given) == 1) STATS[[s]]$Threshold)
                   }))), wn)
                 })))
  
  B <- make_df(Bcolnames, Brownames)
  
  B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(vapply(balance.list, function(y) if (x %in% rownames(y)) y[[x, "Type"]] else NA_character_, character(1))))), use.names = FALSE)
  
  if (include.times) B[["Times"]] <- vapply(Brownames, function(x) paste(seq_along(balance.list)[vapply(balance.list, function(y) x %in% rownames(y), logical(1L))], collapse = ", "), character(1))[Brownames]
  else B[["Times"]] <- NULL
  
  for (Agg.Fun in Agg.Funs) {
    for (s in compute[compute %in% all_STATS(type)]) {
      abs0 <- function(x) {if (is_null(x)) NA_real_ else if (abs) STATS[[s]]$abs(x) else (x)}
      agg <- function(x, ...) {
        if (!any(is.finite(x))) NA_real_
        else if (s == "variance.ratios" && tolower(Agg.Fun) == "mean") geom.mean(x)
        else if (tolower(Agg.Fun) == "rms") sqrt(mean_fast(STATS[[s]]$abs(x)^2, TRUE))
        else get(tolower(Agg.Fun))(x, ...)
      }
      for (sample in c("Un", weight.names)) {
        if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
          if (sample != "Un" || !get_from_STATS("adj_only")[s]) {
            B[[paste.(Agg.Fun, STATS[[s]]$bal.tab_column_prefix, sample)]] <- vapply(Brownames, function(x) agg(unlist(lapply(balance.list, function(y) if (x %in% rownames(y)) abs0(y[[x, paste.(STATS[[s]]$bal.tab_column_prefix, sample)]]))), na.rm = TRUE), numeric(1))
          }
        }
      }
    }
  }
  
  if (length(Agg.Funs.Given) == 1) {
    #Assign X.Threshold values
    for (s in compute[compute %in% all_STATS(type)]) {
      if (is_not_null(thresholds[[s]])) {
        if (!get_from_STATS("adj_only")[s] && no.adj) {
          B[[paste.(STATS[[s]]$Threshold, "Un")]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, "Un")]]), "", 
                                                             STATS[[s]]$abs(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, "Un")]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                             paste0("Not Balanced, >", round(thresholds[[s]], 3)))
        }
        else {
          for (i in weight.names) {
            B[[paste.(STATS[[s]]$Threshold, i)]] <- ifelse_(B[["Type"]]=="Distance" | !is.finite(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, i)]]), "", 
                                                            STATS[[s]]$abs(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, i)]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3)),
                                                            paste0("Not Balanced, >", round(thresholds[[s]], 3)))
          }
        }
      }
      if (no.adj || length(weight.names) <= 1) names(B)[names(B) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
    }
  }
  
  return(B)
}


#base.bal.tab.imp
samplesize.across.imps <- function(obs.list) {
  obs.list <- clear_null(obs.list)
  
  obs <- Reduce("+", obs.list)/length(obs.list)
  attr(obs, "tag") <- paste0("Average ", tolower(attr(obs.list[[1]], "tag")), " across imputations")
  return(obs)
}

#base.bal.tab.multi
samplesize.multi <- function(bal.tab.multi.list, treat_names, focal = NULL) {
  if (is_not_null(focal)) which <- c(treat_names[treat_names != focal], focal)
  else which <- treat_names
  bal.tab.multi.list <- clear_null(bal.tab.multi.list)
  obs <- do.call("cbind", unname(lapply(bal.tab.multi.list, function(x) x[["Observations"]])))[, which]
  attr(obs, "tag") <- attr(bal.tab.multi.list[[1]][["Observations"]], "tag")
  attr(obs, "ss.type") <- attr(bal.tab.multi.list[[1]][["Observations"]], "ss.type")
  return(obs)
}

#base.bal.tab.msm
samplesize.msm <- function(bal.tab.msm.list) {
  obs <- do.call("cbind", lapply(bal.tab.msm.list, function(x) x[["Observations"]]))
  attr(obs, "tag") <- attr(bal.tab.msm.list[[1]][["Observations"]], "tag")
  attr(obs, "ss.type") <- attr(bal.tab.msm.list[[1]][["Observations"]], "ss.type")
  return(obs)
}

#base.bal.tab.cluster
samplesize.across.clusters <- function(obs.list) {
  obs.list <- clear_null(obs.list)
  obs <- Reduce("+", obs.list)
  attr(obs, "tag") <- paste0("Total ", tolower(attr(obs.list[[1]], "tag")), " across clusters")
  return(obs)
}

#base.bal.tab.subclass
balance.table.subclass <- function(C, type, weights = NULL, treat, subclass,
                                   continuous, binary, s.d.denom, 
                                   thresholds = list(), un = FALSE, disp = NULL, stats = NULL, 
                                   s.weights = rep(1, length(treat)), abs = FALSE, 
                                   var_types = NULL, quick = TRUE, ...) {
  #Creates list SB of balance tables for each subclass
  #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
  
  disp <- unique(c(disp, all_STATS(type)[all_STATS(type) %in% stats]))
  compute <- if (quick) disp else c("means", "sds", all_STATS(type))
  
  #B=Balance frame
  Bnames <- c("Type", 
              expand.grid_string(c(if (type == "bin") expand.grid_string(c("M"["means" %in% compute], "SD"["sds" %in% compute]), c("0", "1"), collapse = ".")
                                   else if (type == "cont") c("M"["means" %in% compute], "SD"["sds" %in% compute]), 
                                   unlist(lapply(compute[compute %in% all_STATS(type)], function(s) {
                                     c(STATS[[s]]$bal.tab_column_prefix,
                                       if (is_not_null(thresholds[[s]])) STATS[[s]]$Threshold)
                                   }))
              ),
              c("Adj"), collapse = "."))
  B <- make_df(Bnames, colnames(C))
  
  #Set var type (binary/continuous)
  B[["Type"]] <- if_null_then(var_types, get.types(C))
  bin.vars <- B[["Type"]] == "Binary"
  
  SB <- setNames(lapply(levels(subclass), function(i) B), levels(subclass))
  
  binary <- match_arg(binary, c("raw", "std"))
  sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else !bin.vars
  
  if (type == "bin") 
    subclass_w_empty <- vapply(levels(subclass), function(i) {
      any(vapply(treat_vals(treat), function(t) !any(treat == t & subclass == i), logical(1L)))
    }, logical(1L))
  else {
    subclass_w_empty <- vapply(levels(subclass), function(i) !any(subclass == i), logical(1L))
  }
  
  for (i in levels(subclass)) {
    
    in.subclass <- !is.na(subclass) & subclass==i
    
    #Means for each group
    if ("means" %in% compute) {
      if (type == "bin") {
        tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
        for (t in c("0", "1")) {
          SB[[i]][[paste.("M", t, "Adj")]] <- col_w_mean(C, subset = treat==tn01[t] & in.subclass, s.weights = s.weights)
        }
      }
      else if (type == "cont") {
        SB[[i]][["M.Adj"]] <- col_w_mean(C, subset = in.subclass, s.weights = s.weights)
      }
    }
    
    #SDs for each group
    if ("sds" %in% compute) {
      if (type == "bin") {
        for (t in c("0", "1")) {
          sds <- rep(NA_real_, NCOL(C))
          sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == tn01[t] & in.subclass, s.weights = s.weights)
          SB[[i]][[paste.("SD", t, "Adj")]] <- sds
        }
      }
      else if (type == "cont") {
        sds <- rep(NA_real_, NCOL(C))
        sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == in.subclass, s.weights = s.weights)
        SB[[i]][["SD.Adj"]] <- sds
      }
    }
    
    for (s in all_STATS(type)) {
      if (s %in% compute && !subclass_w_empty[i]) {
        SB[[i]][[paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")]] <- STATS[[s]]$fun(C, treat = treat, weights = NULL, 
                                                                                     std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                                     s.d.denom = s.d.denom,
                                                                                     abs = FALSE, s.weights = s.weights, 
                                                                                     bin.vars = bin.vars, subset = in.subclass)
        
        if (all(vapply(SB[[i]][paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")], 
                       function(x) all(!is.finite(x)), logical(1L)))) {
          disp <- disp[disp != s] 
          thresholds[[s]] <- NULL
        }
        
        
        if (is_not_null(thresholds[[s]])) {
          
          SB[[i]][[paste.(STATS[[s]]$Threshold, "Adj")]] <- ifelse(SB[[i]][["Type"]]!="Distance" & is.finite(SB[[i]][[paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")]]), 
                                                                   paste0(ifelse(abs_(SB[[i]][[paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")]]) < thresholds[[s]], "Balanced, <", "Not Balanced, >"), round(thresholds[[s]], 3)), "")
        }
        names(SB[[i]])[names(SB[[i]]) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
      }
      
    }
    
  }
  
  if (all(vapply(SB, function(sb) all(vapply(sb[startsWith(names(sb), "SD.")], function(x) all(!is.finite(x)), logical(1L))), logical(1L)))) {
    disp <- disp[disp != "sds"]
  }
  
  attr(SB, "thresholds") <- thresholds
  attr(SB, "disp") <- disp
  attr(SB, "compute") <- compute
  
  return(SB)
}

# !!! NEEDS TO BE UPDATED !!!
balance.table.across.subclass.cont <- function(balance.table, balance.table.subclass.list, subclass.obs, r.threshold = NULL) {
  #NEEDS TO BE UPDATED
  
  B.A <- balance.table.subclass.list[[1]][names(balance.table.subclass.list[[1]]) %in% c("M.Adj", "SD.Adj", "Corr.Adj")]
  
  for(i in rownames(B.A)) {
    for(j in colnames(B.A)) {
      if (startsWith(j, "SD.")) {
        B.A[[i, j]] <- sqrt(sum(vapply(seq_along(balance.table.subclass.list),
                                       function(s) subclass.obs[[s]]/sum(subclass.obs) * (balance.table.subclass.list[[s]][[i, j]]^2), numeric(1))))
      }
      else {
        B.A[[i, j]] <- sum(vapply(seq_along(balance.table.subclass.list),
                                  function(s) subclass.obs[[s]]/sum(subclass.obs) * (balance.table.subclass.list[[s]][[i, j]]), numeric(1)))
        
      }
    }
  }
  B.A.df <- cbind(balance.table[c("Type", "M.Un", "SD.Un", "Corr.Un", "R.Threshold.Un")], 
                       B.A, R.Threshold = NA_character_)
  if (is_not_null(r.threshold)) B.A.df[["R.Threshold"]] <- ifelse(B.A.df[["Type"]]=="Distance", "", paste0(ifelse(is.finite(B.A.df[["Corr.Adj"]]) & abs_(B.A.df[["Corr.Adj"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold))
  return(B.A.df)
}

#love.plot
isColor <- function(x) {
  tryCatch(is.matrix(grDevices::col2rgb(x)), 
           error = function(e) FALSE)
}
f.recode <- function(f, ...) {
  #Simplified version of forcats::fct_recode
  f <- factor(f)
  new_levels <- unlist(list(...), use.names = TRUE)
  old_levels <- levels(f)
  idx <- match(new_levels, old_levels)
  
  old_levels[idx] <- names(new_levels)
  
  levels(f) <- old_levels
  return(f)
}
seq_int_cycle <- function(begin, end, max) {
  seq(begin, end, by = 1) - max*(seq(begin-1, end-1, by = 1) %/% max)
}
assign.shapes <- function(colors, default.shape = "circle") {
  if (nunique(colors) < length(colors)) {
    shapes <- seq_int_cycle(19, 19 + length(colors) - 1, max = 25)
  }
  else shapes <- rep.int(default.shape, length(colors))
  return(shapes)
}
shapes.ok <- function(shapes, nshapes) {
  shape_names <- c(
    "circle", paste("circle", c("open", "filled", "cross", "plus", "small")), "bullet",
    "square", paste("square", c("open", "filled", "cross", "plus", "triangle")),
    "diamond", paste("diamond", c("open", "filled", "plus")),
    "triangle", paste("triangle", c("open", "filled", "square")),
    paste("triangle down", c("open", "filled")),
    "plus", "cross", "asterisk"
  )
  shape_nums <- 1:25
  return((length(shapes) == 1 || length(shapes) == nshapes) && ((is.numeric(shapes) && all(shapes %in% shape_nums)) || (is.character(shapes) && all(shapes %in% shape_names))))
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
ggarrange_simple <- function (plots, nrow = NULL, ncol = NULL) {
  #A thin version of egg:ggarrange
  
  gtable_frame <- function (g, width = grid::unit(1, "null"), height = grid::unit(1, "null")) {
    panels <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]),]
    pargins <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]),]
    ll <- unique(panels$l)
    margins <- if (length(ll) == 1) grid::unit(0, "pt") else g$widths[ll[-length(ll)] + 2]
    tt <- unique(panels$t)
    fixed_ar <- g$respect
    if (fixed_ar) {
      ar <- as.numeric(g$heights[tt[1]])/as.numeric(g$widths[ll[1]])
      height <- width * (ar/length(ll))
      g$respect <- FALSE
    }
    core <- g[seq(min(tt), max(tt)), seq(min(ll), max(ll))]
    top <- g[seq(1, min(tt) - 1), seq(min(ll), max(ll))]
    bottom <- g[seq(max(tt) + 1, nrow(g)), seq(min(ll), max(ll))]
    left <- g[seq(min(tt), max(tt)), seq(1, min(ll) - 1)]
    right <- g[seq(min(tt), max(tt)), seq(max(ll) + 1, ncol(g))]
    fg <- grid::nullGrob()
    if (length(left)) {
      lg <- gtable::gtable_add_cols(left, grid::unit(1, "null"), 0)
      lg <- gtable::gtable_add_grob(lg, fg, 1, l = 1)
    }
    else {
      lg <- fg
    }
    if (length(right)) {
      rg <- gtable::gtable_add_cols(right, grid::unit(1, "null"))
      rg <- gtable::gtable_add_grob(rg, fg, 1, l = ncol(rg))
    }
    else {
      rg <- fg
    }
    if (length(top)) {
      tg <- gtable::gtable_add_rows(top, grid::unit(1, "null"), 0)
      tg <- gtable::gtable_add_grob(tg, fg, t = 1, l = 1)
    }
    else {
      tg <- fg
    }
    if (length(bottom)) {
      bg <- gtable::gtable_add_rows(bottom, grid::unit(1, "null"), 
                                    -1)
      bg <- gtable::gtable_add_grob(bg, fg, t = nrow(bg), l = 1)
    }
    else {
      bg <- fg
    }
    grobs <- list(fg, tg, fg, lg, core, rg, fg, bg, fg)
    widths <- grid::unit.c(sum(left$widths), width, sum(right$widths))
    heights <- grid::unit.c(sum(top$heights), height, sum(bottom$heights))
    all <- gtable::gtable_matrix("all", grobs = matrix(grobs, ncol = 3, nrow = 3, byrow = TRUE), 
                                 widths = widths, heights = heights)
    
    all[["layout"]][5, "name"] <- "panel"
    if (fixed_ar) 
      all$respect <- TRUE
    all
  }
  
  n <- length(plots)
  
  grobs <- lapply(plots, ggplot2::ggplotGrob)
  
  if (is_null(nrow) && is_null(ncol)) {
    nm <- grDevices::n2mfrow(n)
    nrow <- nm[1]
    ncol <- nm[2]
  }
  
  hw <- lapply(rep(1, n), grid::unit, "null")
  
  fg <- lapply(seq_along(plots), function(i) gtable_frame(g = grobs[[i]], 
                                                          width = hw[[i]], height = hw[[i]]))
  
  spl <- split(fg, rep(1, n))
  
  rows <- lapply(spl, function(r) do.call(gridExtra::gtable_cbind, r))
  
  gt <- do.call(gridExtra::gtable_rbind, rows)
  
  invisible(gt)
}
bal.tab_class_sequence <- function(b) {
  if (is_(b, c("bal.tab.bin", "bal.tab.cont"))) return(NULL)
  else {
    b_ <- b[[which(endsWith(names(b), ".Balance"))]][[1]]
    return(c(class(b)[1], bal.tab_class_sequence(b_)))
  }
}
unpack_bal.tab <- function(b) {
  unpack_bal.tab_internal <- function(b) {
    if (is_(b, c("bal.tab.bin", "bal.tab.cont"))) return(b[["Balance"]])
    else {
      b_ <- b[[which(endsWith(names(b), ".Balance"))]]
      
      b_list <- lapply(b_, function(i) {
        if (is_(i, c("bal.tab.bin", "bal.tab.cont"))) return(i[["Balance"]])
        else return(unpack_bal.tab_internal(i))
      })
      return(b_list)
    }
  }
  LinearizeNestedList <- function(NList, NameSep) {
    # LinearizeNestedList:
    #
    # https://sites.google.com/site/akhilsbehl/geekspace/
    #         articles/r/linearize_nested_lists_in_r
    #
    # Akhil S Bhel
    # 
    
    if (is.data.frame(NList)) return(NList)
    
    A <- 1
    B <- length(NList)
    
    while (A <= B) {
      Element <- NList[[A]]
      EName <- names(NList)[A]
      if (is.list(Element)) {
        
        if (A == 1) {
          Before <- NULL
        } else {
          Before <- NList[1:(A - 1)]
        }
        if (A == B) {
          After <- NULL
        } else {
          After <- NList[(A + 1):B]
        }
        
        if (is.data.frame(Element)) {
          Jump <- 1
        } else {
          NList[[A]] <- NULL
          
          Element <- LinearizeNestedList(Element, NameSep)
          names(Element) <- paste(EName, names(Element), sep=NameSep)
          Jump <- length(Element)
          NList <- c(Before, Element, After)
        }
      } else {
        Jump <- 1
      }
      A <- A + Jump
      B <- length(NList)
    }
    return(NList)
  }
  
  namesep <- paste(c("|", do.call(c, lapply(1:20, function(i) sample(LETTERS, 1))), "|"), collapse = "")
  
  out_ <- unpack_bal.tab_internal(b)
  out <- LinearizeNestedList(out_, NameSep = namesep)
  
  attr(out, "namesep") <- namesep
  attr(out, "class_sequence") <- bal.tab_class_sequence(b)
  
  return(out)
}

#bal.plot

#print.bal.tab
print.data.frame_ <- function(x, ...) {
  if (is_not_null(x) && NROW(x) > 0 && NCOL(x) > 0) {
    print.data.frame(x, ...)
  }
}

#set.cobalt.options
acceptable.options <- function() {
  TF <- c(TRUE, FALSE)
  return(list(stats = c("mean.diffs"),
              un = TF,
              continuous = c("raw", "std"),
              binary = c("raw", "std"),
              imbalanced.only = TF,
              disp = c("means", "sds"),
              disp.means = TF,
              disp.sds = TF,
              disp.v.ratio = TF,
              disp.ks = TF,
              disp.subclass = TF,
              disp.bal.tab = TF,
              cluster.summary = TF,
              cluster.fun = c("min", "mean", "max"),
              imp.summary = TF,
              imp.fun = c("min", "mean", "max"),
              multi.summary = TF,
              msm.summary = TF,
              target.summary = TF,
              subclass.summary = TF,
              int_sep = " * ",
              factor_sep = "_",
              center = TF,
              remove_perfect_col = TF,
              disp.call = TF))
}

#Misc
`%+%` <- function(...) {
  if (is_(..1, "atomic") && is_(..2, "atomic")) crayon::`%+%`(as.character(..1), as.character(..2))
  else ggplot2::`%+%`(...)
}

#On attach
.onAttach <- function(...) {
  cobaltLib <- dirname(system.file(package = "cobalt"))
  version <- packageDescription("cobalt", lib.loc = cobaltLib)$Version
  BuildDate <- packageDescription("cobalt", lib.loc = cobaltLib)$Date
  
  foo <- paste0(" cobalt (Version ", version, ", Build Date: ", BuildDate, ")")
  packageStartupMessage(foo)
}

.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}
