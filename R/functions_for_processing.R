#bal.tab
process_obj <- function(obj) {
  if (is_null(obj)) {
    return(set_class(list(), "cobalt.processed.obj"))
  }
  
  if (isS4(obj)) {
    obj <- asS3(obj)
  }
  
  new_class <- character()
  
  #npCBPS
  if (inherits(obj, "npCBPS")) {
    new_class <- c("CBPS", "npCBPS")
  }
  #ebalance.trim
  else if (inherits(obj, "ebalance.trim")) {
    new_class <- "ebalance"
  }
  #Matchby
  else if (inherits(obj, "Matchby")) {
    new_class <- "Match"
  }
  #cem.match.list
  else if (inherits(obj, "cem.match.list")) {
    new_class <- c("cem.match", "cem.match.list")
  }
  #time.list
  else if (is.list(obj) && !is.data.frame(obj) &&
           all_apply(obj, rlang::is_formula)) {
    new_class <- c("formula.list", "time.list", class(obj))
  }
  else if (is.list(obj) && !is.data.frame(obj) &&
           all_apply(obj, is.data.frame)) {
    new_class <- c("data.frame.list", "time.list", class(obj))
  }
  #designmatch
  else if (is.list(obj) && length(obj) >= 6L) {
    dm.b.names <- c("obj_total", "obj_dist_mat", "t_id", 
                    "c_id", "group_id", "time")
    dm.n.names <- c("obj_total", "obj_dist_mat", "id_1", 
                    "id_2", "group_id", "time")
    if (all(dm.b.names %in% names(obj)) || all(dm.n.names %in% names(obj))) {
      new_class <- "designmatch"
    }
  }
  
  if (is_null(new_class)) {
    set_class(obj, "cobalt.processed.obj", .replace = FALSE)
  }
  else {
    set_class(obj, c(new_class, "cobalt.processed.obj"))
  }
}

#x2base
process_treat <- function(treat, ..., keep_values = FALSE) {
  
  .chk_not_missing(treat, "`treat`")
  
  if (inherits(treat, "unprocessed.treat")) {
    attrs <- attributes(treat)
    renamed_original <- setNames(names(treat_vals(treat)), treat_vals(treat))
    treat <- factor(renamed_original[as.character(treat)], levels = renamed_original)
    for (at in c("treat_names", "treat_vals", "keep_values", "treat.type", "names")) {
      attr(treat, at) <- attrs[[at]]
    }
  }
  else {
    # keep_values <- isTRUE(attr(treat, "keep_values")) || 
    #     (has.treat.type(treat) && get.treat.type(treat) == "multinomial")
    
    treat <- .process_vector(treat, name = "treat", 
                             which = "treatment statuses", 
                             datalist = list(...), missing.okay = FALSE) |>
      assign.treat.type()
    
    treat.type <- get.treat.type(treat)
    
    if (treat.type == "binary") {
      if (!is.factor(treat)) {
        treat <- factor(treat, levels = sort(unique(treat, nmax = 2L)))
      }
      
      original_values <- intersect(levels(treat), unique(treat, nmax = 2L))
      
      treat_names(treat) <- {
        if (!keep_values && can_str2num(as.character(treat)) &&
            all(original_values %in% c("0", "1"))) {
          setNames(c("Control", "Treated"), c("control", "treated"))
        }
        else {
          setNames(original_values, c("control", "treated"))
        }
      }
      
      treat_vals(treat) <- setNames(original_values, treat_names(treat))
    }
    else if (treat.type == "multinomial") {
      treat <- factor(treat, ordered = FALSE)
      treat_names(treat) <- setNames(levels(treat), levels(treat))
      treat_vals(treat) <- setNames(levels(treat), treat_names(treat))
    }
    attr(treat, "treat.type") <- treat.type
    # attr(treat, "keep_values") <- keep_values
  }
  
  set_class(treat, "processed.treat", .replace = FALSE)
}
unprocess_treat <- function(treat) {
  if (!inherits(treat, "processed.treat")) {
    return(treat)
  }
  
  attrs <- attributes(treat)
  treat <- treat_vals(treat)[as.character(treat)]
  attributes(treat) <- attrs[setdiff(names(attrs), "class")]
  
  set_class(treat, c("unprocessed.treat", class(treat_vals(treat))))
}
process_treat.list <- function(treat.list, ...) {
  .chk_not_missing(treat.list, "`treat.list`")
  
  if (!is.list(treat.list)) {
    treat.list <- as.list(treat.list)
  }
  
  hasdots <- ...length() > 0L
  
  treat.list.names <- vapply(seq_along(treat.list), function(ti) {
    if (hasdots && is.character(treat.list[[ti]]) && length(treat.list[[ti]]) == 1L) {
      treat.list[[ti]]
    }
    else if (rlang::is_named(treat.list)) names(treat.list)[ti]
    else as.character(ti)
  }, character(1L))
  
  lapply(treat.list, process_treat, ...) |>
    setNames(treat.list.names)
}
`treat_names<-` <- function(treat, value) {
  `attr<-`(treat, "treat_names", value)
}
treat_names <- function(treat) {
  attr(treat, "treat_names", TRUE)
}
`treat_vals<-` <- function(treat, value) {
  `attr<-`(treat, "treat_vals", value)
}
treat_vals <- function(treat) {
  attr(treat, "treat_vals", TRUE)
}
subset_processed.treat <- function(x, index) {
  y <- x[index]
  treat_names(y) <- treat_names(x)[treat_vals(x) %in% unique(y)]
  treat_vals(y) <- treat_vals(x)[treat_vals(x) %in% unique(y)]
  
  assign.treat.type(y) |>
    set_class(class(x))
}

initialize_X <- function() {
  c("covs",
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
    "thresholds") |>
    make_list()
}
initialize_X_msm <- function() {
  c("covs.list",
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
    "thresholds") |>
    make_list()
}
.weight_check <- function(w) {
  wname <- deparse1(substitute(w))
  if (!is.list(w)) w <- list(w)
  
  if (anyNA(w, recursive = TRUE)) {
    .err(sprintf("`NA`s are not allowed in the %s", wname))
  }
  
  for (x in w) {
    if (!all(is.numeric(x))) {
      .err(sprintf("all %s must be numeric", wname))
    }
    
    if (!all(is.finite(x))) {
      .err(sprintf("infinite %s are not allowed", wname))
    }
  }
}
.cluster_check <- function(cluster = NULL, treat) {
  if (is_null(cluster)) {
    return(invisible())
  }
  
  if (!is.list(treat)) {
    treat <- list(treat)
  }
  
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
  
  if (stop_warn["cw"]) {
    .wrn("some clusters have only one unit in them, which may yield unxpected results")
  }
  
  if (stop_warn["bw"]) {
    .wrn("some clusters have only one member of a treatment group in them, which may yield unxpected results")
  }
  
  if (stop_warn["bs"]) {
    .err("not all treatment levels are present in all clusters")
  }
  
}
strata2weights <- function(strata, treat, estimand = NULL, focal = NULL) {
  #Process strata into weights (similar to weight.subclass from MatchIt)
  
  #Checks
  if (!is.atomic(strata) || is_not_null(dim(strata))) {
    .err("`strata` must be an atomic vector or factor")
  }
  
  #Process treat
  treat <- process_treat(treat)
  
  if (get.treat.type(treat) == "continuous") {
    .err("`strata` cannot be turned into weights for continuous treatments")
  }
  
  s.d.denom <- .get_s.d.denom(NULL, estimand = estimand, subclass = strata, treat = treat,
                              focal = focal, quietly = TRUE)
  focal <- {
    if (s.d.denom %in% treat_vals(treat)) process_focal(s.d.denom, treat)
    else NULL
  }
  
  NAsub <- is.na(strata)
  imat <- do.call("cbind", lapply(treat_vals(treat), function(t) treat == t & !NAsub))
  colnames(imat) <- treat_vals(treat)
  
  weights <- rep.int(0, length(treat))
  
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
      weights[imat[,t]] <- {
        if (t == focal)  1
        else (t_by_sub[focal,] / t_by_sub[t,])[strata.c[imat[,t]]]
      }
    }
  }
  else {
    for (t in treat_vals(treat)) {
      weights[imat[,t]] <- (total_by_sub / t_by_sub[t,])[strata.c[imat[,t]]]
    }
  }
  
  na.w <- !is.finite(weights)
  if (any(na.w)) {
    weights[na.w] <- 0
    .wrn("some units were given weights of zero due to zeros in stratum membership")
  }
  
  if (all(check_if_zero(weights))) {
    .err("no units were stratified")
  }
  
  for (tnn in names(treat_names(treat))) {
    if (all(check_if_zero(weights[treat == treat_vals(treat)[treat_names(treat)[tnn]]]))) {
      .err(sprintf("no %s units were stratified", tnn))
    }
  }
  
  attr(weights, "match.strata") <- strata
  weights
}

.use_tc_fd <- function(formula = NULL, data = NULL, treat = NULL, covs = NULL, needs.treat = TRUE, needs.covs = TRUE) {
  
  treat_f <- treat_c <- covs_f <- covs_c <- NULL
  
  if (is_not_null(formula)) {
    formula <- try(as.formula(formula), silent = TRUE)
    if (rlang::is_formula(formula)) {
      D <- NULL
      
      if (is_not_null(data)) {
        D <- data
      }
      
      if (is_mat_like(covs)) {
        D <- {
          if (is_not_null(D)) cbind(D, as.data.frame(covs)) 
          else as.data.frame(covs)
        }
      }
      
      treat_f <- try(get_treat_from_formula(formula, D, treat = treat), silent = TRUE)
      covs_f <- try(get_covs_from_formula(formula, D), silent = TRUE)
    }
  }
  
  if (is_not_null(covs)) {
    if (is_mat_like(covs)) {
      covs_c <- try(get_covs_from_formula(data = covs), silent = TRUE)
    }
    else if (is.character(covs)) {
      if (!is_mat_like(data)) {
        .err("if `covs` is a character vector, `data` must be specified as a data.frame")
      }
      
      if (!all(covs %in% colnames(data))) {
        .err("all entries in `covs` must be names of variables in `data`")
      }
      
      covs_c <- try(get_covs_from_formula(f.build(covs), data = as.data.frame(data)), silent = TRUE)
    }
    else {
      .err("`covs` must be a data.frame of covariates")
    }
  }
  
  if (is_not_null(treat)) {
    treat_c <- try({
      if (!is.atomic(treat)) .err("`treat` must be an vector of treatment statuses")
      if (is.character(treat) && length(treat) == 1L) get_treat_from_formula(reformulate(".", treat), data = data)
      else get_treat_from_formula(treat ~ ., treat = treat)
    }, silent = TRUE)
  }
  
  covs_to_use <- treat_to_use <- "c"
  if (is_error(covs_c)) {
    if (is_error(covs_f)) {
      .err(attr(covs_c, "condition")$message, tidy = FALSE)
    }
    covs_to_use <- "f"
  }
  else if (is_null(covs_c)) {
    if (is_error(covs_f)) {
      .err(attr(covs_f, "condition")$message, tidy = FALSE)
    }
    if (is_null(covs_f) && needs.covs) {
      .err("no covariates were specified")
    }
    covs_to_use <- "f"
  }
  
  if (is_error(treat_c)) {
    if (is_error(treat_f)) {
      .err(attr(treat_c, "condition")$message, tidy = FALSE)
    }
    treat_to_use <- "f"
  }
  else if (is_null(treat_c)) {
    if (is_error(treat_f)) {
      .err(attr(treat_f, "condition")$message, tidy = FALSE)
    }
    if (is_null(treat_f) && needs.treat) {
      .err("no treatment variable was specified")
    }
    treat_to_use <- "f"
  }
  
  t.c <- list(treat = switch(treat_to_use, "f" = treat_f, "c" = treat_c), 
              covs = switch(covs_to_use, "f" = covs_f, "c" = covs_c))
  t.c$treat.name <- attr(t.c$treat, "treat.name")
  
  t.c
}
.process_val <- function(val, i, treat = NULL, covs = NULL, addl.data = list(), ...) {
  
  if (is.data.frame(val)) {
    return(val)
  }
  
  if (is.numeric(val)) {
    return(setNames(data.frame(val), i))
  }
  
  if (!is.character(val)) {
    if (i == "weights") {
      .err("the argument supplied to `weights` must be a named list of weights, names of variables containing weights in an available data set, or objects with a `get.w()` method")
    }
    else {
      .err(sprintf("the argument supplied to %s must be a vector, a data.frame, or the names of variables in an available data set",
                   add_quotes(i, "`")))
    }
  }
  
  addl.data <- clear_null(addl.data)
  if ((is_null(covs) && is_null(treat)) ||
      length(val) == NROW(covs) ||
      length(val) == length(treat) ||
      (is_not_null(addl.data) && length(val) > max(vapply(addl.data, ncol, numeric(1L))))) {
    return(setNames(data.frame(val), i))
  }
  
  if (is_null(addl.data)) {
    .wrn(sprintf("names were provided to %s, but no argument to `data` was provided. Ignoring %s",
                 add_quotes(i, "`"), add_quotes(i, "`")))
    return(NULL)
  }
  
  val <- unique(val)
  val.df <- make_df(val, nrow = max(vapply(addl.data, nrow, numeric(1L))))
  not.found <- rlang::rep_named(val, TRUE)
  for (v in val) {
    for (k in seq_along(addl.data)) {
      if (utils::hasName(addl.data[[k]], v)) {
        val.df[[v]] <- addl.data[[k]][[v]]
        not.found[v] <- FALSE
        break
      }
    }
  }
  
  if (!any(not.found)) {
    return(val.df)
  }
  
  .wrn(sprintf("the following variable(s) named in %s are not in any available data sets and will be ignored: %s",
               add_quotes(i, "`"), paste(val[not.found])))
  
  val.df[!not.found]
}
.process_data_frame <- function(i, df, treat = NULL, covs = NULL, addl.data = list(), ...) {
  if (is_null(df)) {
    return(NULL)
  }
  
  val <- df
  val.df <- NULL
  
  if (i == "weights" && any(has_method(class(val), "get.w"))) {
    val <- list(val)
  }
  
  if (!is.list(val) || is.data.frame(val)) {
    return(.process_val(val, i, treat, covs, addl.data = addl.data))
  }
  
  if (i == "weights") {
    #Use get.w() on inputs
    for (x in seq_along(val)) {
      val[[x]] <- process_obj(val[[x]])
      has_get.w_method <- has_method(class(val[[x]]), "get.w")
      
      if (!any(has_get.w_method)) {
        next
      }
      
      get.w_class <- class(val[[x]])[has_get.w_method][1L]
      
      val[[x]] <- get.w(val[[x]], treat = treat, ...)
      
      if (!nzchar(rlang::names2(val)[x])) {
        names(val)[x] <- get.w_class
      }
    }
  }
  
  val.list <- lapply(val, .process_val, i, treat, covs, addl.data = addl.data)
  
  if (!rlang::is_named(val.list)) {
    .err(sprintf("all entries in `%s` must have names", i))
  }
  
  if (!all_the_same(vapply(val.list, NROW, numeric(1L)))) {
    .err(sprintf("not all items in `%s` have the same length", i))
  }
  
  for (j in which(vapply(val.list, NCOL, numeric(1L)) == 1L)) {
    names(val.list[[j]]) <- names(val.list)[j]
  }
  
  do.call("cbind", val.list) |>
    setNames(unlist(lapply(val.list, names)))
}
.process_list <- function(i, List, ntimes, call.phrase, treat.list = list(), covs.list = list(), addl.data = list(), ...) {
  if (is_null(List)) {
    return(NULL)
  }
  
  val.List <- List
  
  if (!is.list(val.List)) {
    val.List <- list(val.List)
  }
  
  if (length(val.List) == 1L) {
    val.List <- replicate(ntimes, val.List)
  }
  else if (length(val.List) != ntimes) {
    .err(sprintf("the argument to %s must be a list of the same length as the number of time points in %s",
                 add_quotes(i, "`"), call.phrase))
  }
  
  for (ti in which(lengths(val.List) > 0L)) {
    val <- val.List[[ti]]
    val.df <- NULL
    
    if (is.list(val) && !is.data.frame(val)) {
      val.list <- lapply(val, .process_val, strsplit(i, ".list", fixed = TRUE)[[1L]],
                         treat.list[[ti]], covs.list[[ti]], addl.data = addl.data)
      
      if (!all_the_same(vapply(val.list, NROW, numeric(1L)))) {
        .err(sprintf("not all items in `%s` have the same length", i))
      }
      
      for (j in which(vapply(val.list, NCOL, numeric(1L)) == 1L)) {
        names(val.list[[j]]) <- names(val.list)[j]
      }
      
      val.df <- setNames(do.call("cbind", val.list),
                         vapply(val.list, names, character(1L)))
    }
    else {
      val.df <- .process_val(val, strsplit(i, ".list", fixed = TRUE)[[1L]],
                             treat.list[[ti]], covs.list[[ti]], addl.data = addl.data)
    }
    
    if (is_not_null(val.df) && anyNA(val.df)) {
      .err(sprintf("missing values exist in %s", add_quotes(i, "`")))
    }
    
    val.List[[ti]] <- val.df
  }
  
  val.df.lengths <- val.List |> clear_null() |>
    vapply(nrow, numeric(1L))
  
  if (!all_the_same(val.df.lengths)) {
    .err(sprintf("all columns in `%s` need to have the same number of rows", i))
  }
  
  val.List
}
.process_vector <- function(vec, name = deparse1(substitute(vec)), which = name, datalist = list(), missing.okay = FALSE) {
  bad.vec <- FALSE
  if (is.character(vec) && length(vec) == 1L && is_not_null(datalist)) {
    for (i in seq_along(datalist)) {
      if (is.matrix(datalist[[i]]) && vec %in% colnames(datalist[[i]])) {
        vec <- datalist[[i]][,vec]
        break
      }
      
      if (is.data.frame(datalist[[i]]) && utils::hasName(datalist[[i]], vec)) {
        vec <- datalist[[i]][[vec]]
        break
      }
      
      if (i == length(datalist)) {
        bad.vec <- TRUE
      }
    }
  }
  else if (!is.atomic(vec) || length(vec) <= 1L) {
    bad.vec <- TRUE
  }
  
  if (bad.vec) {
    .err(sprintf("the argument to %s must be a vector of %s or the (quoted) name of a variable in `data` that contains %s",
                 add_quotes(name, "`"), which, which))
  }
  
  if (!missing.okay && anyNA(vec)) {
    .err(sprintf("missing values exist in %s", add_quotes(name, "`")))
  }
  
  vec
} 
.get_s.d.denom <- function(s.d.denom = NULL, estimand = NULL, weights = NULL,
                           subclass = NULL, treat = NULL, focal = NULL, quietly = FALSE) {
  check.estimand <- check.weights <- check.focal <- bad.s.d.denom <- bad.estimand <- FALSE
  s.d.denom.specified <- !missing(s.d.denom) && is_not_null(s.d.denom)
  estimand.specified <- is_not_null(estimand)
  
  if (s.d.denom.specified) {
    if (isTRUE(attr(s.d.denom, "checked"))) {
      return(s.d.denom)
    }
    
    if (length(s.d.denom) > 1L && length(s.d.denom) != NCOL(weights)) {
      .err(sprintf("`s.d.denom` must have length 1 or equal to the number of valid sets of weights, which is %s",
                   NCOL(weights)))
    }
    
    allowable.s.d.denoms <- c("pooled", "all", "weighted", "hedges")
    
    if (!all(s.d.denom %in% allowable.s.d.denoms)) {
      
      if (!inherits(treat, "processed.treat")) {
        treat <- process_treat(treat)
      }
      
      unique.treats <- as.character(treat_vals(treat))
      
      if (length(treat_names(treat)) == 2L && all(c("treated", "control") %in% names(treat_names(treat)))) {
        allowable.s.d.denoms <- c(allowable.s.d.denoms, "treated", "control")
      }
      
      if (is_not_null(focal)) {
        allowable.s.d.denoms <- c(allowable.s.d.denoms, "focal")
      }
      
      s.d.denom <- match_arg(s.d.denom, unique(c(unique.treats, allowable.s.d.denoms)), 
                             several.ok = length(weights) > 1L)
      s.d.t.c <- s.d.denom %in% c("treated", "control")
      
      if (any(s.d.t.c) && length(treat_names(treat)) == 2L) {
        s.d.denom[s.d.t.c] <- treat_vals(treat)[treat_names(treat)[s.d.denom[s.d.t.c]]]
      }
      else if (any(s.d.denom == "focal")) {
        check.focal <- TRUE
      }
    }
  }
  else {
    check.estimand <- TRUE
  }
  
  if (check.estimand) {
    if (estimand.specified) {
      if (!inherits(treat, "processed.treat")) {
        treat <- process_treat(treat)
      }
      
      try.estimand <- tryCatch(match_arg(toupper(estimand), c("ATT", "ATC", "ATE", "ATO", "ATM"), several.ok = TRUE),
                               error = function(cond) NA_character_)
      if (anyNA(try.estimand) || any(try.estimand %in% c("ATC", "ATT")) && get.treat.type(treat) != "binary") {
        check.focal <- TRUE
      }
      else {
        if (length(try.estimand) > 1L && length(try.estimand) != NCOL(weights)) {
          .err(sprintf("`estimand` must have length 1 or equal to the number of valid sets of weights, which is %s",
                       NCOL(weights)))
        }
        
        s.d.denom <- vapply(try.estimand, switch, character(1L),
                            ATT = treat_vals(treat)[treat_names(treat)["treated"]], 
                            ATC = treat_vals(treat)[treat_names(treat)["control"]], 
                            ATO = "weighted",
                            ATM = "weighted",
                            "pooled")
      }
    }
    else {
      check.focal <- TRUE
    }
  }
  
  if (check.focal) {
    if (is_not_null(focal) && NCOL(weights) == 1L) {
      s.d.denom <- focal
    }
    else {
      check.weights <- TRUE
    }
  }
  
  if (check.weights) {
    if (!inherits(treat, "processed.treat")) {
      treat <- process_treat(treat)
    }
    
    if (is_null(weights) && is_null(subclass)) {
      s.d.denom <- "pooled"
    }
    else if (is_not_null(subclass)) {
      sub.tab <- table(treat, subclass)[treat_vals(treat), ]
      sub.tab <- rbind(sub.tab, table(subclass)[colnames(sub.tab)])
      dimnames(sub.tab) <- list(c(treat_vals(treat), "pooled"), colnames(sub.tab))
      
      ranges <- apply(sub.tab, 1L, function(x) .mean_abs_dev(x) / sum(x))
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
        "pooled"
      }, character(1L))
    }
  }
  
  if (is_not_null(weights)) {
    if (length(s.d.denom) == 1L && NCOL(weights) > 1L) {
      s.d.denom <- rep.int(s.d.denom, NCOL(weights))
    }
    
    if (length(s.d.denom) != NCOL(weights)) {
      .err(sprintf("valid inputs to `s.d.denom` or `estimand` must have length 1 or equal to the number of valid sets of weights, which is %s",
                   NCOL(weights)))
    }
    
    names(s.d.denom) <- names(weights)
  }
  
  if (s.d.denom.specified && is_null(weights) && "weighted" %in% s.d.denom) {
    attr(s.d.denom, "note") <- 'note: `s.d.denom` specified as "weighted", but no weights supplied; setting to "all"'
  }
  else if (s.d.denom.specified && bad.s.d.denom && (!estimand.specified || bad.estimand)) {
    attr(s.d.denom, "note") <- sprintf("warning: `s.d.denom` should be one of %s. Using %s instead",
                                       word_list(unique(c(unique.treats, allowable.s.d.denoms)), "or", quotes = 2L),
                                       word_list(s.d.denom, quotes = 2L))
  }
  else if (estimand.specified && bad.estimand) {
    attr(s.d.denom, "note") <- sprintf("warning: `estimand` should be one of %s. Ignoring `estimand`",
                                       word_list(c("ATT", "ATC", "ATE"), "or", quotes = 2L))
  }
  else if ((check.focal || check.weights) && !all(s.d.denom %in% treat_vals(treat))) {
    attr(s.d.denom, "note") <- sprintf("note: `s.d.denom` not specified; assuming %s", 
                                       if (all_the_same(s.d.denom)) add_quotes(s.d.denom[1L])
                                       else word_list(paste0(add_quotes(vapply(s.d.denom, function(s) {
                                         if (s %in% treat_vals(treat) && all(treat_vals(treat) %in% c("0", "1"))) {
                                           names(treat_names(treat))[treat_names(treat) == names(treat_vals(treat))[treat_vals(treat) == s]]
                                         }
                                         else s
                                       }, character(1L))), " for ", names(weights))))
  }
  
  if (!quietly && is_not_null(attr(s.d.denom, "note"))) {
    .msg(attr(s.d.denom, "note"))
  }
  
  attr(s.d.denom, "checked") <- TRUE
  
  s.d.denom
}
.get_s.d.denom.cont <- function(s.d.denom, weights = NULL, subclass = NULL, quietly = FALSE) {
  bad.s.d.denom <- FALSE
  s.d.denom.specified <- !missing(s.d.denom) && is_not_null(s.d.denom)
  
  if (s.d.denom.specified) {
    if (isTRUE(attr(s.d.denom, "checked"))) {
      return(s.d.denom)
    }
    
    if (length(s.d.denom) > 1L && length(s.d.denom) != NCOL(weights)) {
      .err("`s.d.denom` must have length 1 or equal to the number of valid sets of weights")
    }
    
    allowable.s.d.denoms <- {
      if (is_not_null(subclass)) "all"
      else c("all", "weighted")
    }
    
    if (!all(s.d.denom %in% allowable.s.d.denoms)) {
      s.d.denom <- match_arg(s.d.denom, unique(allowable.s.d.denoms), 
                             several.ok = length(weights) > 1L)
    }
  }
  else {
    s.d.denom <- "all"
  }
  
  if (is_not_null(weights)) {
    if (length(s.d.denom) == 1L && NCOL(weights) > 1L) {
      s.d.denom <- rep.int(s.d.denom, NCOL(weights))
    }
    
    if (length(s.d.denom) != NCOL(weights)) {
      .err("valid inputs to `s.d.denom` or `estimand` must have length 1 or equal to the number of valid sets of weights")
    }
    
    names(s.d.denom) <- names(weights)
  }
  
  if (s.d.denom.specified && is_null(weights) && any(s.d.denom == "weighted")) {
    attr(s.d.denom, "note") <- 'note: `s.d.denom` specified as "weighted", but no weights supplied; setting to "all"'
  }
  else if (s.d.denom.specified && bad.s.d.denom) {
    attr(s.d.denom, "note") <- sprintf("warning: `s.d.denom` should be %s.\n         Using %s instead",
                                       word_list(unique(allowable.s.d.denoms), "or", quotes = 2L),
                                       word_list(s.d.denom, quotes = 2L))
  }
  
  if (!quietly && is_not_null(attr(s.d.denom, "note"))) {
    .msg(attr(s.d.denom, "note"))
  }
  
  attr(s.d.denom, "checked") <- TRUE
  
  s.d.denom
}
.compute_s.d.denom <- function(mat, treat = NULL, s.d.denom = "pooled", s.weights = NULL,
                               bin.vars = NULL, subset = NULL, weighted.weights = NULL,
                               to.sd = rep.int(TRUE, ncol(mat)), na.rm = TRUE) {
  denoms <- setNames(rep.int(1, ncol(mat)), colnames(mat))
  
  if (is.character(s.d.denom) && length(s.d.denom) == 1L) {
    if (is_null(bin.vars)) {
      bin.vars <- rep.int(FALSE, ncol(mat))
      bin.vars[to.sd] <- is_binary_col(mat[subset, to.sd,drop = FALSE])
    }
    else if (!is.atomic(bin.vars) || length(bin.vars) != ncol(mat) ||
             anyNA(as.logical(bin.vars))) {
      .err("`bin.vars` must be a logical vector with length equal to the number of columns of `mat`")
    }
    
    possibly.supplied <- c("mat", "treat", "weighted.weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, numeric(1L)),
                        possibly.supplied)
    supplied <- lengths > 0L
    if (!all_the_same(lengths[supplied])) {
      .err(sprintf("%s must have the same number of units",
                   word_list(possibly.supplied[supplied], quotes = "`")))
    }
    
    if (lengths["weighted.weights"] == 0L) weighted.weights <- rep.int(1, NROW(mat))
    if (lengths["s.weights"] == 0L) s.weights <- rep.int(1, NROW(mat))
    
    if (lengths["subset"] == 0L) {
      subset <- rep.int(TRUE, NROW(mat))
    }
    else if (anyNA(as.logical(subset))) {
      .err("`subset` must be a logical vector")
    }
    
    if (is_null(treat)) {
      cont.treat <- TRUE
    }
    else {
      if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
      cont.treat <- get.treat.type(treat) == "continuous"
    }
    
    if (cont.treat) {
      unique.treats <- NULL
    }
    else {
      unique.treats <- {
        if (inherits(treat, "processed.treat")) as.character(treat_vals(treat))
        else as.character(unique(treat))
      }
      
      if (s.d.denom %in% c("treated", "control") && s.d.denom %nin% unique.treats) {
        s.d.denom <- treat_vals(treat)[treat_names(treat)[s.d.denom]]
      }
      
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
        df <- length(treat) - length(unique.treats)
        (1 - 3 / (4 * df - 1))^-1 *
          sqrt(Reduce("+", lapply(unique.treats,
                                  function(t) (sum(treat == t) - 1) * col.w.v(mat[treat == t, , drop = FALSE],
                                                                              w = s.weights[treat == t],
                                                                              bin.vars = bin.vars, na.rm = na.rm))) / df)
      }
    else {
      .err("`s.d.denom` is not an allowed value")
    }
    
    denoms[to.sd] <- denom.fun(mat = mat[, to.sd, drop = FALSE], treat = treat, s.weights = s.weights,
                               weighted.weights = weighted.weights, bin.vars = bin.vars[to.sd],
                               unique.treats = unique.treats, na.rm = na.rm)
    
    zero_sds <- !is.finite(denoms[to.sd]) | check_if_zero(denoms[to.sd])
    if (any(zero_sds)) {
      denoms[to.sd][zero_sds] <- sqrt(col.w.v(mat[, to.sd, drop = FALSE][, zero_sds, drop = FALSE],
                                              w = s.weights,
                                              bin.vars = bin.vars[to.sd][zero_sds], na.rm = na.rm))
    }
    
    if (cont.treat && is_not_null(treat)) {
      treat.sd <- denom.fun(mat = treat, treat = NULL, s.weights = s.weights,
                            weighted.weights = weighted.weights, bin.vars = FALSE,
                            unique.treats = NULL, na.rm = na.rm)
      denoms[to.sd] <- denoms[to.sd] * treat.sd
    }
  }
  else if (is.numeric(s.d.denom)) {
    if (is_not_null(names(s.d.denom)) && any(colnames(mat) %in% names(s.d.denom))) {
      nm <- colnames(mat)[colnames(mat) %in% names(s.d.denom)]
      denoms[nm] <- s.d.denom[nm]
    }
    else if (length(s.d.denom) == sum(to.sd)) {
      denoms[to.sd] <- s.d.denom
    }
    else if (length(s.d.denom) == ncol(mat)) {
      denoms[] <- s.d.denom
    }
    else {
      .err("`s.d.denom` must be an allowable value or a numeric vector of with length equal to the number of columns of `mat`. See ?col_w_smd for allowable values")
    }
  }
  else {
    .err("`s.d.denom` must be an allowable value or a numeric vector of with length equal to the number of columns of `mat`. See ?col_w_smd for allowable values")
  }
  
  denoms
}
.assign_X_class <- function(X) {
  X <- clear_null(X)
  
  if (is_not_null(X[["treat"]]) && !has.treat.type(X[["treat"]])) {
    X[["treat"]] <- assign.treat.type(X[["treat"]])
  }
  
  if (is_not_null(X[["subclass"]])) {
    if (get.treat.type(X[["treat"]]) == "binary") X.class <- "subclass.binary"
    else if (get.treat.type(X[["treat"]]) == "continuous") X.class <- "subclass.cont"
    else .err("multi-category treatments are not currently compatible with subclasses")
  }
  else if (is_not_null(X[["cluster"]]) && nlevels(X[["cluster"]]) > 1L) X.class <- "cluster"
  else if (is_not_null(X[["covs.list"]])) X.class <- "msm"
  else if (get.treat.type(X[["treat"]]) == "multinomial") X.class <- "multi"
  else if (is_not_null(X[["imp"]]) && nlevels(X[["imp"]]) > 1L) X.class <- "imp"
  else if (get.treat.type(X[["treat"]]) == "binary") X.class <- "binary"
  else if (get.treat.type(X[["treat"]]) == "continuous") X.class <- "cont"
  else probably.a.bug()
  
  attr(X, "X.class") <- X.class
  
  X
}
.get_length_X <- function(X) {
  if (is_not_null(X[["treat"]])) length(X[["treat"]])
  else if (is_not_null(X[["covs"]])) nrow(X[["covs"]])
  else if (is_not_null(X[["treat.list"]])) length(X[["treat.list"]][[1L]])
  else if (is_not_null(X[["covs.list"]])) nrow(X[["covs.list"]][[1L]])
  else .err("couldn't determine length of `X` components")
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
  if (is_null(subset) || !any(names(X) %in% subsettable())) {
    return(X)
  }
  
  n <- .get_length_X(X)
  
  if (is.logical(subset)) {
    if (length(subset) != n) {
      .err("`subset` must have the same length as the other entries")
    }
    
    if (!any(subset)) {
      .err("all `subset` set to `FALSE`")
    }
    
    if (all(subset)) {
      return(X)
    }
    
    subset <- which(subset)
  }
  else if (!is.numeric(subset)) {
    .err("`subset` must be logical or numeric")
  }
  else if (max(subset) > n) {
    .err("subset indices cannot be higher than the length of the other entries")
  }
  
  subset_X_internal <- function(x, subset) {
    attrs <- attributes(x)
    attrs_to_subset <- names(attrs)[vapply(attrs, function(a) all(len(a) == n), logical(1L))]
    
    out <- {
      if (is_null(x)) x
      else if (inherits(x, "processed.treat")) subset_processed.treat(x, subset)
      else if ((is.matrix(x) || is.data.frame(x))) x[subset, , drop = FALSE]
      else if (is.factor(x)) factor(x[subset], nmax = nlevels(x))
      else if (is.atomic(x)) x[subset]
      else if (is.list(x)) lapply(x, subset_X_internal, subset = subset)
      else x
    }
    
    if (is_null(attrs)) {
      return(out)
    }
    
    if (inherits(x, "processed.treat")) {
      return(process_treat(out, keep_values = TRUE))
    }
    
    if (is_null(attrs_to_subset)) {
      for (i in setdiff(names(attrs), names(attributes(out)))) {
        attr(out, i) <- attrs[[i]]
      }
      
      return(out)
    }
    
    subsetted_attrs <- lapply(attrs[attrs_to_subset],
                              subset_X_internal, subset = subset)
    for (i in setdiff(names(attrs), names(attributes(out)))) {
      attr(out, i) <- {
        if (i %in% attrs_to_subset) subsetted_attrs[[i]]
        else attrs[[i]]
      }
    }
    
    out
  }
  
  for (i in intersect(names(X), subsettable())) {
    X[[i]] <- subset_X_internal(X[[i]], subset)
  }
  
  X
}
.mids_complete <- function(data) {
  if (!inherits(data, "mids")) {
    .err("`data` not of class `mids`")
  }
  
  single.complete <- function(ell, data, where = NULL, imp) {
    if (is_null(where)) where <- is.na(data)
    idx <- seq_col(data)[which(colSums(where) > 0L)]
    
    for (j in idx) {
      if (is_null(imp[[j]])) is.na(data[where[, j], j]) <- TRUE
      else data[where[, j], j] <- imp[[j]][, ell]
    }
    
    data
  }
  
  m <- as.integer(data$m)
  idx <- seq_len(m)
  
  mylist <- lapply(idx, single.complete, data = data$data,
                   where = data$where, imp = data$imp)
  
  cmp <- cbind(.imp = rep(idx, each = nrow(data$data)), 
               .id = rep.int(seq_row(data$data), length(idx)), 
               as.data.frame(do.call("rbind", mylist)))
  
  row.names(cmp) <- {
    if (is.integer(attr(data$data, "row.names"))) seq_row(cmp)
    else as.character(seq_row(cmp))
  }
  
  cmp
}

length_imp_process <- function(vectors = NULL, data.frames = NULL, lists = NULL,
                               imp = NULL, data = NULL, original.call.to = NULL,
                               env = sys.frame(-1L)) {
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
      
      for (i in all.objects[lengths > 0 & lengths != length(imp)]) {
        if (lengths[i] != imp.lengths[1L]) {
          problematic[i] <- TRUE
          next
        }
        
        i_obj <- get(i, envir = env, inherits = FALSE)
        
        if (i %in% vectors) {
          new_i <- i_obj[rep(seq_along(i_obj), length(imp.lengths))]
          
          if (unsorted.imp) {
            for (i_ in levels(imp)) {
              new_i[imp == i_] <- i_obj
            }
          }
        }
        else if (i %in% data.frames) {
          new_i <- i_obj[rep(seq_row(i_obj), length(imp.lengths)), , drop = FALSE]
          
          if (unsorted.imp) {
            for (i_ in levels(imp)) {
              new_i[imp == i_,] <- i_obj
            }
          }
        }
        else if (i %in% lists) {
          new_i <- lapply(i_obj, function(j) {
            if (!is.factor(j) && !is_mat_like(j)) {
              .err(sprintf("%s can only contain vectors or data frames",
                           add_quotes(i, "`")))
            }
            
            if (is.factor(j)) {
              newj <- j[rep(seq_along(j), length(imp.lengths))]
              if (unsorted.imp) {
                for (i_ in levels(imp)) {
                  newj[imp == i_] <- j
                }
              }
              
              return(newj)
            }
            
            newj <- j[rep(seq_row(j), length(imp.lengths)), , drop = FALSE]
            if (unsorted.imp) {
              for (i_ in levels(imp)) {
                newj[imp == i_,] <- j
              }
            }
            
            newj
          })
        }
        
        assign(i, new_i, pos = env)
        
      }
    }
    else {
      problematic <- lengths > 0L & lengths != length(imp)
    }
    
    if (any(problematic)) {
      .err(sprintf("%s must have the same number of observations as `imp`",
                   word_list(names(problematic)[problematic], quotes = "`")))
    }
    
    ensure.equal.lengths <- FALSE
    
    assign("imp", imp, pos = env)
  }
  
  #Ensure all input lengths are the same.
  anchor <- {
    if ("treat" %in% all.objects) "treat"
    else if ("treat.list" %in% all.objects) "treat.list"
    else all.objects[which(lengths[all.objects] != 0L)[1L]]
  }
  
  if (ensure.equal.lengths) {
    problematic[lengths %nin% c(0L, lengths[anchor])] <- TRUE
  }
  
  if (any(problematic)) {
    if (is_not_null(original.call.to)) {
      anchor <- sprintf("in the original call to %s",
                        add_quotes(original.call.to, "`"))
    }
    
    .err("%s must have the same number of observations as %s",
         word_list(names(problematic)[problematic], quotes = "`"),
         anchor)
  }
}
process_imp <- function(imp = NULL, ...) {
  if (is_null(imp)) {
    return(NULL)
  }
  
  .process_vector(imp, 
                  datalist = list(...),
                  name = "imp", 
                  which = "imputation identifiers",
                  missing.okay = FALSE) |>
    factor()
}
process_stats <- function(stats = NULL, treat) {
  if (is.list(treat) && !is.data.frame(treat)) {
    stats.list <- lapply(treat, function(x) process_stats(stats, x))
    
    type <- {
      if (all_the_same(vapply(stats.list, attr, character(1L), "type")))
        attr(stats.list[[1L]], "type")
      else NULL
    }
    
    stats <- unique(unlist(stats.list))
    attr(stats, "type") <- type
    
    return(stats)
  }
  
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
  
  stats
}
process_thresholds <- function(thresholds, stats) {
  if (is_null(thresholds)) {
    return(list())
  }
  
  if (is.list(thresholds)) {
    thresholds <- unlist(thresholds)
  }
  
  if (!all(is.na(thresholds)) && !is.numeric(thresholds)) {
    .err("`thresholds` must be numeric")
  }
  
  if (rlang::is_named(thresholds)) {
    names(thresholds) <- stats[pmatch(names(thresholds), stats, duplicates.ok = TRUE)]
    thresholds <- thresholds[!is.na(names(thresholds))]
  }
  else {
    names(thresholds) <- stats[seq_len(min(length(stats), length(thresholds)))]
  }
  
  thresholds[names(thresholds)] <- as.numeric(thresholds)
  
  as.list(na.rem(thresholds))
}
process_subset <- function(subset = NULL, n) {
  if (is_null(subset)) {
    return(NULL)
  }
  
  if (!is.logical(subset) && !is.numeric(subset)) {
    .err("the argument to `subset` must be a logical or numeric vector")
  }
  
  if (is.numeric(subset)) {
    if (any(abs(subset) > n)) {
      .err("numeric values for `subset` cannot be larger than the number of units")
    }
    subset <- subset[!is.na(subset) & subset != 0]
    
    if (any(subset < 0) && any(subset > 0)) {
      .err("positive and negative indices cannot be mixed with `subset`")
    }
    
    if (any(abs(subset) > n)) {
      .err("if `subset` is numeric, none of its values can exceed the number of units")
    }
    
    logical.subset <- rep.int(any(subset < 0), n)
    logical.subset[abs(subset)] <- !logical.subset[abs(subset)]
    subset <- logical.subset
  }
  
  if (anyNA(subset)) {
    .wrn("NAs were present in `subset`. Treating them like FALSE")
    subset[is.na(subset)] <- FALSE
  }
  
  subset
}
process_cluster <- function(cluster = NULL, ...) {
  if (is_null(cluster)) {
    return(NULL)
  }
  
  .process_vector(cluster, 
                  datalist = list(...),
                  name = "cluster", 
                  which = "cluster membership",
                  missing.okay = FALSE) |>
    factor()
}
process_s.weights <- function(s.weights = NULL, ...) {
  if (is_null(s.weights)) {
    return(NULL)
  }
  
  s.weights <- .process_vector(s.weights, 
                               datalist = list(...),
                               name = "s.weights", 
                               which = "sampling weights",
                               missing.okay = FALSE)
  .weight_check(s.weights)
  
  s.weights
}
process_focal <- function(focal = NULL, treat) {
  if (is_null(focal) || get.treat.type(treat) == "continuous") {
    return(NULL)
  }
  
  if (!is.numeric(focal)) {
    if (focal %in% treat_vals(treat)) {
      return(focal)
    }
    
    .err("the name specified to `focal` is not the name of any treatment group")
  }
  
  if (can_str2num(treat) && focal %in% str2num(treat)) {
    return(as.character(focal))
  }
  
  if (focal > length(treat_vals(treat))) {
    .err(sprintf("`focal` was specified as %s, but there are only %s treatment groups",
                 focal, length(treat_vals(treat))))
  }
  
  treat_vals(treat)[focal]
}
process_weights <- function(obj = NULL, A = NULL, treat = NULL, covs = NULL,
                            method = character(), addl.data = list(), ...) {
  A[["x"]] <- NULL
  A[["treat"]] <- NULL
  
  weights <- list()
  if (is_not_null(obj)) {
    if (any(has_method(class(obj), "get.w"))) {
      weights <- do.call("get.w", c(list(obj, treat = treat), A))
      
      weights <- {
        if (is_null(weights)) list()
        else if (is_mat_like(weights)) as.data.frame(weights)
        else setNames(data.frame(weights), class(obj)[has_method(class(obj), "get.w")][1L])
      }
    }
    else {
      weights <- obj[["weights"]]
      weights <- {
        if (is_null(weights)) list()
        else if (is_mat_like(weights)) as.data.frame(weights)
        else setNames(data.frame(weights), class(obj)[1L])
      }
    }
  }
  
  addl.weights <- .process_data_frame("weights", A[["weights"]], treat, covs,
                                      addl.data = addl.data, ...)
  if (is_not_null(addl.weights)) {
    if (is_null(A[["method"]])) {
      addl.methods <- rep.int("weighting", ncol(addl.weights))
    }
    else if (length(A[["method"]]) == 1L) {
      addl.methods <- rep.int(match_arg(A[["method"]], c("weighting", "matching")), ncol(addl.weights))
    }
    else {
      addl.methods <- match_arg(A[["method"]], c("weighting", "matching"), several.ok = TRUE)
      
      if (length(addl.methods) != ncol(addl.weights)) {
        .err("Valid inputs to `method` must have length 1 or equal to the number of valid sets of additional weights")
      }
    }
    
    w.names <- c(names(weights), names(addl.weights))
    unique.names <- unique(w.names)
    if (length(unique.names) != length(w.names)) {
      for (i in unique.names) {
        if (sum(w.names == i) > 1L) {
          w.names[w.names == i] <- make.unique(c(i, w.names[w.names == i]))[-1L]
        }
      } 
    }
    
    if (is_null(weights)) {
      weights <- setNames(addl.weights, w.names)
      method <- setNames(addl.methods, w.names)
    }
    else {
      weights <- setNames(cbind(weights, addl.weights), w.names)
      method <- setNames(c(method, addl.methods), w.names)
    }
  }
  
  .weight_check(weights)
  
  attr(weights, "method") <- method
  
  weights
}
process_disp <- function(disp = NULL, ...) {
  .chk_null_or(disp, .chk_character)
  disp <- {
    if (is_null(disp)) getOption("cobalt_disp")
    else match_arg(disp, acceptable.options()[["disp"]], several.ok = TRUE)
  }
  
  for (d in c("means", "sds")) {
    if (getOption(paste.("cobalt_disp", d), FALSE)) {
      disp <- unique(c(disp, d))
    }
    
    disp.d <- ...get(paste.("disp", d))
    
    if (is_not_null(disp.d)) {
      if (!chk::vld_flag(disp.d)) {
        .err(sprintf("`disp.%s` must be `TRUE` or `FALSE`", d))
      }
      
      disp <- unique(c(disp, d[disp.d]))
      
      if (disp.d) disp <- unique(c(disp, d))
      else disp <- unique(disp[disp != d])
    }
  }
  
  disp
}
process_addl <- function(addl = NULL, datalist = list()) {
  if (is_not_null(addl) && !is.atomic(addl) && !rlang::is_formula(addl) &&
      !is.matrix(addl) && !is.data.frame(addl)) {
    .err("`addl` must be a formula or variable containing the distance values")
  }
  
  data <- do.call("data.frame", unname(clear_null(datalist)))
  
  if (is.atomic(addl) && 
      (!is.character(addl) || is_null(datalist) ||
       length(addl) == nrow(data))) {
    addl <- data.frame(addl = addl)
  }
  else if (is.character(addl)) {
    addl <- reformulate(addl)
  }
  
  addl_t.c <- {
    if (is.data.frame(addl)) {
      .use_tc_fd(data = data, covs = addl, 
                 needs.treat = FALSE, needs.covs = FALSE)
    }
    else {
      .use_tc_fd(formula = addl, data = data, 
                 needs.treat = FALSE, needs.covs = FALSE)
    }
  }
  
  addl_t.c[["covs"]]
}
process_addl.list <- function(addl.list = NULL, datalist = list(), covs.list = list()) {
  datalist <- clear_null(c(datalist, covs.list))
  
  if (is.list(addl.list) && !is.data.frame(addl.list)) {
    if (length(addl.list) != length(covs.list)) {
      .err("`addl` must have an entry for each time point")
    }
    
    return(lapply(addl.list, process_addl, datalist = datalist))
  }
  
  addl <- process_addl(addl.list, datalist = datalist)
  
  lapply(seq_along(covs.list), function(x) addl)
}
process_distance <- function(distance = NULL, datalist = list(), obj.distance = NULL,
                             obj.distance.name = "distance") {
  
  if (is_not_null(distance) && !is.atomic(distance) && !rlang::is_formula(distance) &&
      !is.matrix(distance) && !is.data.frame(distance)) {
    .err("`distance` must be a formula or variable containing the distance values")
  }
  
  data <- do.call("data.frame", unname(clear_null(datalist)))
  
  if (is.atomic(distance) && 
      (!is.character(distance) || is_null(datalist) ||
       length(distance) == nrow(data))) {
    distance <- data.frame(distance = distance)
  }
  else if (is.character(distance)) {
    distance <- reformulate(distance)
  }
  
  distance_t.c <- {
    if (is.data.frame(distance)) {
      .use_tc_fd(data = data, covs = distance, 
                 needs.treat = FALSE, needs.covs = FALSE)
    }
    else {
      .use_tc_fd(formula = distance, data = data, 
                 needs.treat = FALSE, needs.covs = FALSE)
    }
  }
  
  if (!(is_not_null(obj.distance) && !all(is.na(obj.distance)))) {
    return(distance_t.c[["covs"]])
  }
  
  obj.distance <- setNames(data.frame(obj.distance), obj.distance.name)
  obj.distance <- get_covs_from_formula(data = obj.distance)
  
  if (is_null(distance_t.c[["covs"]])) {
    return(obj.distance)
  }
  
  co.cbind(distance_t.c[["covs"]], obj.distance)
}
process_distance.list <- function(distance.list = NULL, datalist = list(),
                                  covs.list = list(), obj.distance = NULL, obj.distance.name = "distance") {
  datalist <- clear_null(c(datalist, covs.list))
  
  if (is_null(obj.distance)) {
    obj.distance <- lapply(seq_along(covs.list), function(x) NULL)
  }
  else if (!is.list(obj.distance) || is.data.frame(obj.distance)) {
    obj.distance <- lapply(seq_along(covs.list), function(x) obj.distance)
  }
  
  if (is_null(distance.list)) {
    distance.list.out <- lapply(seq_along(covs.list), function(x) {
      process_distance(NULL, datalist = datalist, 
                       obj.distance = obj.distance[[x]], obj.distance.name = obj.distance.name)})
  }
  else if (is.list(distance.list) && !is.data.frame(distance.list)) {
    if (length(distance.list) != length(covs.list)) {
      .err("`distance` must have an entry for each time point")
    }
    
    distance.list.out <- lapply(seq_along(distance.list), function(x) {
      process_distance(distance.list[[x]], datalist = datalist, 
                       obj.distance = obj.distance[[x]],
                       obj.distance.name = obj.distance.name)})
  }
  else {
    distance.list.out <- lapply(seq_along(covs.list), function(x) {
      process_distance(distance.list, datalist = datalist, 
                       obj.distance = obj.distance[[x]],
                       obj.distance.name = obj.distance.name)})
  }
  
  distance.list.out
}
process_focal_and_estimand <- function(focal, estimand, treat, treated = NULL) {
  reported.estimand <- estimand
  
  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)
  
  unique.treat <- unique(treat, nmax = switch(treat.type, "binary" = 2L, "multinomial" = length(treat) / 4))
  
  #Check focal
  if (is_not_null(focal) && (length(focal) > 1L || focal %nin% unique.treat)) {
    .err("the argument supplied to `focal` must be the name of a level of treatment")
  }
  
  if (treat.type == "multinomial") {
    
    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      .wrn(sprintf('%s is not compatible with `focal`. Setting `estimand` to "ATT"',
                   add_quotes(estimand)))
      reported.estimand <- estimand <- "ATT"
    }
    
    if (estimand == "ATT") {
      if (is_null(focal)) {
        if (is_null(treated) || treated %nin% unique.treat) {
          .err('when estimand = "ATT" for multinomial treatments, an argument must be supplied to `focal`')
        }
        focal <- treated
      }
    }
    else if (estimand == "ATC" && is_null(focal)) {
      .err('when estimand = "ATC" for multinomial treatments, an argument must be supplied to `focal`')
    }
  }
  else if (treat.type == "binary") {
    unique.treat.bin <- unique(binarize(treat), nmax = 2L)
    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      .wrn(sprintf("%s is not compatible with `focal`. Setting `estimand` to \"ATT\"",
                   add_quotes(estimand)))
      reported.estimand <- estimand <- "ATT"
    }
    
    if (is_not_null(treated) && treated %in% unique.treat) {
      if (is_null(focal)) {
        if (estimand == "ATT")
          focal <- treated
        else if (estimand == "ATC")
          focal <- setdiff(unique.treat, treated)
      }
      
      if (estimand == "ATC") {
        estimand <- "ATT"
      }
    }
    else if (is_null(focal)) {
      if (all(as.character(unique.treat.bin) == as.character(unique.treat))) {
        treated <- unique.treat[unique.treat.bin == 1]
      }
      else {
        if (is.factor(treat)) treated <- levels(treat)[2L]
        else treated <- unique.treat[unique.treat.bin == 1]
        
        if (estimand == "ATT") {
          .msg(sprintf("assuming %s the treated level. If not, supply an argument to `focal`",
                       word_list(treated, quotes = !is.numeric(treat), is.are = TRUE)))
          
        }
        else if (estimand == "ATC") {
          .msg(sprintf("assuming %s the control level. If not, supply an argument to `focal`",
                       word_list(setdiff(unique.treat, treated),
                                 quotes = !is.numeric(treat), is.are = TRUE)))
        }
      }
      
      if (estimand == "ATT")
        focal <- treated
      else if (estimand == "ATC")
        focal <- setdiff(unique.treat, treated)
    }
    else if (estimand == "ATT") {
      treated <- focal
    }
    else if (estimand == "ATC") {
      treated <- setdiff(unique.treat, focal)
    }
    
    if (estimand == "ATC") {
      estimand <- "ATT"
    }
  }
  
  list(focal = as.character(focal),
       estimand = estimand,
       reported.estimand = reported.estimand,
       treated = if (is.factor(treated)) as.character(treated) else treated)
}

#.get_C2
get_ints_from_co.names <- function(co.names) {
  if (is_null(co.names)) {
    return(list())
  }
  
  clear_null(lapply(co.names, function(co) {
    if ("isep" %in% co[["type"]]) {
      co[["type"]] <- c(co[["type"]], "isep")
      which_isep <- which(co[["type"]] == "isep")
      vapply(seq_along(which_isep), function(x) {
        if (x == 1L) paste(co[["component"]][1:(which_isep[1L] - 1L)], collapse = "")
        else paste(co[["component"]][(which_isep[x - 1L] + 1L):(which_isep[x] - 1L)], collapse = "")
      }, character(1L))
    }
    else NULL
  }))
}
get_treat_from_formula <- function(f, data = NULL, treat = NULL) {
  
  if (is.character(f)) {
    f <- try(as.formula(f), silent = TRUE)
  }
  
  .chk_formula(f)
  
  env <- rlang::f_env(f)
  
  f <- update(f, ~ 0)
  
  #Check if data exists
  if (is_null(data)) {
    data <- env
    data.specified <- FALSE
  }
  else if (is.data.frame(data)) {
    data.specified <- TRUE
  }
  else {
    .wrn("the argument supplied to `data` is not a data.frame object. Ignoring `data`")
    data <- env
    data.specified <- FALSE
  }
  
  tt <- try_chk(terms(f, data = data))
  
  if (rlang::is_formula(tt, lhs = TRUE)) {
    resp.vars.mentioned <- as.character(rlang::f_lhs(tt))
    resp.vars.failed <- vapply(resp.vars.mentioned, function(v) {
      test <- tryCatch(eval(str2expression(v), data, env), error = function(e) e)
      
      if (inherits(test, "simpleError")) {
        if (!identical(conditionMessage(test), sprintf("object '%s' not found", v))) {
          .err(conditionMessage(test), tidy = FALSE)
        }
        
        return(TRUE)
      }
      
      if (is.function(test)) {
        .err(sprintf("invalid type (function) for variable '%s'", v))
      }
      
      is_null(test)
    }, logical(1L))
    
    if (any(resp.vars.failed)) {
      if (is_null(treat)) {
        .err(sprintf("the given response variable, %s, is not a variable in %s",
                     add_quotes(resp.vars.mentioned),
                     word_list(c("`data`", "the global environment")[c(data.specified, TRUE)], "or")))
      }
      tt <- delete.response(tt)
    }
  }
  else {
    resp.vars.failed <- TRUE
  }
  
  if (all(resp.vars.failed)) {
    treat.name <- NULL
  }
  else {
    treat.name <- resp.vars.mentioned[!resp.vars.failed][1L]
    treat <- eval(str2expression(treat.name), data, env)
  }
  
  attr(treat, "treat.name") <- treat.name
  
  treat
}
get_covs_from_formula <- function(f, data = NULL, factor_sep = "_", int_sep = " * ") {
  
  rebuild_f <- function(ttfactors, tics = FALSE) {
    #Set tics = TRUE if returned formula is used with tmpcovs,
    #i.e., a model.frame. If used with data, leave FALSE
    if (tics) {
      no_tics <- !startsWith(rownames(ttfactors), "`")
      rownames(ttfactors)[no_tics] <- add_quotes(rownames(ttfactors)[no_tics], "`")
    }
    
    as.formula(paste("~ 0 +", paste(vapply(seq_col(ttfactors), 
                                           function(x) paste(rownames(ttfactors)[ttfactors[, x] > 0], collapse = ":"),
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
    
    if (after == ncol(ttfactor)) {
      return(cbind(ttfactor, addcol))
    }
    
    cbind(ttfactor[,seq_len(after), drop = FALSE], 
          addcol, 
          ttfactor[,-seq_len(after), drop = FALSE])
    
  }
  
  #Check if data exists
  data.specified <- FALSE
  if (is_not_null(data)) {
    if (is.matrix(data)) {
      data <- as.data.frame.matrix(data)
    }
    
    if (!is.data.frame(data)) {
      .err("the argument supplied to `data` must be a data.frame object")
    }
    
    data.specified <- TRUE
  }
  
  if (missing(f)) {
    f <- f.build(names(data))
  }
  else {
    if (is.character(f)) {
      f <- try(as.formula(f), silent = TRUE)
    }
    .chk_formula(f)
  }
  
  env <- rlang::f_env(f)
  if (!data.specified) data <- env
  
  # rlang::f_lhs(f) <- NULL
  
  tt <- tryCatch(terms(f, data = data),
                 error = function(e) {
                   if (conditionMessage(e) == "'.' in formula and no 'data' argument") {
                     .err("'.' is not allowed in formulas")
                   }
                   .err(conditionMessage(e))
                 })
  
  #Process RHS 
  tt.covs <- delete.response(tt)
  attr(tt.covs, "intercept") <- 0
  
  ttfactors <- attr(tt.covs, "factors")
  ttvars <- setNames(vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1L], rownames(ttfactors))
  
  rhs.df.type <- setNames(vapply(ttvars, function(v) {
    if (length(dim(try(eval(str2expression(add_quotes(v, "`")), data, env), silent = TRUE))) == 2L) "lit"
    else if (length(dim(try(eval(str2expression(v), data, env), silent = TRUE))) == 2L) "exp"
    else "not.a.df"
  }, character(1L)), ttvars)
  
  rhs.df <- setNames(rhs.df.type != "not.a.df", ttvars)
  
  if (any(rhs.df)) {
    term_is_interaction <- colSums(ttfactors != 0) > 1L
    
    if (any_apply(which(rhs.df), function(x) any(ttfactors[x,] != 0 & term_is_interaction))) {
      .err("interactions with data.frames are not allowed in the input formula")
    }
    
    addl.dfs <- setNames(lapply(ttvars[rhs.df], function(v) {
      df <- switch(rhs.df.type[v],
                   "lit" = eval(str2expression(add_quotes(v, "`")), data, env),
                   eval(str2expression(v), data, env))
      
      if (inherits(df, "rms")) {
        return(setNames(as.data.frame.matrix(as.matrix(df)), colnames(df)))
      }
      
      if (is.data.frame(df)) {
        #Deal with the fact that data.frames may contain matrices and data.frames, which
        #may contain data/frames, and so on
        non.vec.col <- which(vapply(df, function(x) is_not_null(dim(x)), logical(1L)))
        while (is_not_null(non.vec.col)) {
          for (i in non.vec.col) {
            if (NCOL(df[[i]]) == 1L && is_null(colnames(df[[i]]))) colnames(df[[i]]) <- names(df)[i]
            else if (can_str2num(colnames(df[[i]]))) colnames(df[[i]]) <- paste(names(df)[i], colnames(df[[i]]), sep = "_")
          }
          names(df)[non.vec.col] <- ""
          
          df <- as.data.frame(do.call("cbind", df))
          non.vec.col <- which(vapply(df, function(x) is_not_null(dim(x)), logical(1L)))
        }
        
        if (ncol(df) == 1L && is_null(colnames(df))) colnames(df) <- v
        else if (can_str2num(colnames(df))) colnames(df) <- paste(v, colnames(df), sep = "_")
        
        return(df)
      }
      
      if (ncol(df) == 1L && is_null(colnames(df))) colnames(df) <- v
      else if (can_str2num(colnames(df))) colnames(df) <- paste(v, colnames(df), sep = "_")
      
      as.data.frame(df)
    }),
    ttvars[rhs.df])
    
    for (i in colnames(ttfactors)[colnames(ttfactors) %in% names(ttvars)[rhs.df]]) {
      for (j in seq_col(addl.dfs[[ttvars[i]]])) {
        if (names(addl.dfs[[ttvars[i]]])[j] %in% c(ttvars[!rhs.df], unlist(lapply(addl.dfs[seq_len(which(names(addl.dfs) == ttvars[i]) - 1L)], names)))) {
          names(addl.dfs[[ttvars[i]]])[j] <- paste0(ttvars[i], "_", names(addl.dfs[[ttvars[i]]])[j])
        }
      }
      ind <- which(colnames(ttfactors) == i)
      ttfactors <- append.ttfactor(ttfactors,
                                   add_quotes(names(addl.dfs[[ttvars[i]]]), "`"),
                                   ind)[, -ind, drop = FALSE]
    }
    
    if (data.specified) {
      data <- do.call("cbind", unname(c(addl.dfs, list(data))))
    }
    else {
      data <- do.call("cbind", unname(addl.dfs))
      data.specified <- TRUE
    }
    
    tt.covs <- rebuild_f(ttfactors) |>
      terms(data = data)
    
    ttfactors <- attr(tt.covs, "factors")
    ttvars <- setNames(vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1L], rownames(ttfactors))
  }
  
  #Check to make sure variables are valid
  original_ttvars <- rownames(ttfactors)
  for (i in seq_along(rownames(ttfactors))) {
    #Check if evaluable
    #If not, check if evaluable after changing to literal using ``
    #If not, stop()
    #If eventually evaluable, check if function
    
    evaled.var <- try(eval(str2expression(rownames(ttfactors)[i]), data, env),
                      silent = TRUE)
    
    if (null_or_error(evaled.var)) {
      evaled.var <- try(eval(str2expression(add_quotes(rownames(ttfactors)[i], "`")), data, env),
                        silent = TRUE)
      
      if (null_or_error(evaled.var)) {
        ee <- conditionMessage(attr(evaled.var, "condition"))
        
        if (startsWith(ee, "object '") && endsWith(ee, "' not found")) {
          v <- sub("object '([^']+)' not found", "\\1", ee)
          .err(sprintf("the variable %s cannot be found. Be sure it is entered correctly or supply a dataset that contains this varialble to `data`",
                       add_quotes(v, 2L)))
        }
        
        .err(ee)
      }
      
      rownames(ttfactors)[i] <- add_quotes(rownames(ttfactors)[i], "`")
    }
    
    if (is.function(evaled.var)) {
      .err(sprintf("invalid type (function) for variable '%s'",
                   rownames(ttfactors)[i]))
    }
  }
  
  if (!identical(original_ttvars, rownames(ttfactors))) {
    tt.covs <- rebuild_f(ttfactors) |>
      terms(data = data)
    
    ttfactors <- attr(tt.covs, "factors")
    ttvars <- vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1L]
  }
  
  tmpcovs <- try_chk(model.frame2(tt.covs, data))
  
  for (i in ttvars) {
    if (is_binary(tmpcovs[[i]])) {
      tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = 2L)
      next
    }
    
    if (is.character(tmpcovs[[i]])) {
      tmpcovs[[i]] <- factor(tmpcovs[[i]])
    }
    
    if (is.factor(tmpcovs[[i]])) {
      tmpcovs[[i]] <- {
        if (nlevels(tmpcovs[[i]]) == 1L) factor(tmpcovs[[i]], levels = c(paste0(".", levels(tmpcovs[[i]])), levels(tmpcovs[[i]])))
        else factor(tmpcovs[[i]], nmax = nlevels(tmpcovs[[i]]))
      }
    }
  }
  
  #Process NAs: make NA variables
  if (anyNA(tmpcovs)) {
    vars_with_NA <- colnames(tmpcovs)[anyNA_col(tmpcovs)]
    
    for (i in rev(vars_with_NA)) {
      #Find which of ttlabels i first appears, and put `i: <NA>` after it
      ind <- 1L
      for (x in seq_along(colnames(ttfactors))) {
        if (i %in% c(colnames(ttfactors)[x], all.vars(str2expression(colnames(ttfactors)[x])))) {
          ind <- x
          break
        }
      }
      ttfactors <- append.ttfactor(ttfactors, 
                                   sprintf("`%s:<NA>`", i),
                                   ind)
      
      tmpcovs[[sprintf("%s:<NA>", i)]] <- as.numeric(is.na(tmpcovs[[i]]))
    }
    
    tt.covs <- rebuild_f(ttfactors, tics = TRUE) |>
      terms(data = tmpcovs)
    
    ttfactors <- attr(tt.covs, "factors")
    
    ttvars <- vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1L]
    
    na_vars <- paste0(vars_with_NA, ":<NA>")
    
    tmpcovs <- try_chk(model.frame2(tt.covs, tmpcovs))
    
    for (i in setdiff(ttvars, na_vars)) {
      if (is_binary(tmpcovs[[i]])) {
        tmpcovs[[i]] <- factor(tmpcovs[[i]], nmax = 2L)
        next
      }
      
      if (is.character(tmpcovs[[i]])) {
        tmpcovs[[i]] <- factor(tmpcovs[[i]])
      }
      
      if (is.factor(tmpcovs[[i]])) {
        tmpcovs[[i]] <- {
          if (nlevels(tmpcovs[[i]]) == 1L) factor(tmpcovs[[i]], levels = c(paste0(".", levels(tmpcovs[[i]])), levels(tmpcovs[[i]])))
          else factor(tmpcovs[[i]], nmax = nlevels(tmpcovs[[i]]))
        }
      }
    }
  }
  else {
    na_vars <- character()
  }
  
  #Re-check ttfactors
  original_ttvars <- rownames(ttfactors)
  for (i in seq_along(rownames(ttfactors))) {
    #Check if evaluable in tmpcovs
    #If not, check if evaluable i tmpcovs after changing to literal using ``
    #If not, stop() (shouldn't occur)
    
    evaled.var <- try(eval(str2expression(rownames(ttfactors)[i]), tmpcovs), silent = TRUE)
    if (null_or_error(evaled.var)) {
      evaled.var <- try(eval(str2expression(add_quotes(rownames(ttfactors)[i], "`")), tmpcovs), silent = TRUE)
      if (null_or_error(evaled.var)) {
        .err(conditionMessage(attr(evaled.var, "condition")))
      }
      rownames(ttfactors)[i] <- add_quotes(rownames(ttfactors)[i], "`")
    }
  }
  
  if (!identical(original_ttvars, rownames(ttfactors))) {
    tt.covs <- rebuild_f(ttfactors) |>
      terms(data = data)
    
    ttfactors <- attr(tt.covs, "factors")
    ttvars <- vapply(attr(tt.covs, "variables"), deparse1, character(1L))[-1L]
  }
  
  tmpcovs <- model.frame2(tt.covs, data = tmpcovs, drop.unused.levels = TRUE)
  
  #Check for infinite values
  covs.with.inf <- vapply(tmpcovs, function(x) is.numeric(x) && any(!is.na(x) & !is.finite(x)), logical(1L))
  if (any(covs.with.inf)) {
    s <- if (sum(covs.with.inf) == 1L) c("", "s") else c("s", "")
    .err(sprintf("the variable%s %s contain%s non-finite values, which are not allowed",
                 s[1L],
                 word_list(names(tmpcovs)[covs.with.inf], quotes = 1L),
                 s[2L]))
  }
  
  attr(tt.covs, "intercept") <- 1 #Add intercept to correctly process single-level factors
  mm <- model.matrix(tt.covs, data = tmpcovs, 
                     contrasts.arg = lapply(Filter(is.factor, tmpcovs),
                                            function(x) contrasts(x, contrasts = all_the_same(x))))
  
  rownames(ttfactors) <- trim_string(rownames(ttfactors), "`")
  
  mmassign <- attr(mm, "assign")[-1L]
  mmassign2 <- setNames(factor(mmassign, levels = sort(unique(mmassign), na.last = TRUE),
                               labels = colnames(ttfactors)), colnames(mm)[-1L])
  
  vars_in_each_term <- setNames(lapply(colnames(ttfactors), function(x) {
    rownames(ttfactors)[ttfactors[,x] != 0]
  }), colnames(ttfactors))
  
  all_factor_levels <- lapply(vars_in_each_term, function(v) {
    do.call("expand.grid", c(clear_null(setNames(lapply(v, function(fa) colnames(attr(mm, "contrasts")[[fa]])), v)),
                             list(stringsAsFactors = FALSE)))
  })
  
  expanded <- setNames(lapply(seq_along(mmassign2), function(x) {
    .terms <- vars_in_each_term[[mmassign2[x]]]
    k <- sum(seq_along(mmassign2) <= x & mmassign2 == mmassign2[x])
    setNames(lapply(.terms, function(t) {
      if (t %in% names(all_factor_levels[[mmassign2[x]]])) {
        all_factor_levels[[mmassign2[x]]][[t]][k]
      }
      else character()
    }), .terms)
  }), names(mmassign2))
  
  #component types: base, fsep, isep, power, na, level
  co.names <- lapply(expanded, function(x) {
    Reduce(function(x1, x2 = NULL) {
      list(component = c(x1[["component"]], int_sep, x2[["component"]]),
           type = c(x1[["type"]], "isep", x2[["type"]]))
    },
    lapply(seq_along(x), function(i) {
      base <- gsub("`", "", names(x)[i], fixed = TRUE)
      if (base %in% na_vars) {
        base <- substr(base, 1L, nchar(base) - 5L)
        list(component = c(base, ":<NA>"),
             type = c("base", "na"))
      }
      else if (is_null(x[[i]])) {
        list(component = base,
             type = "base")
      }
      else {
        list(component = c(base, factor_sep, x[[i]]),
             type = c("base", "fsep", "level"))
      }
    }))
  })
  
  names(co.names) <- vapply(co.names, function(x) paste(x[["component"]], collapse = ""), character(1L))
  covs <- clear_attr(mm)[, -1L, drop = FALSE] #drop the intercept
  
  attr(co.names, "seps") <- c(factor = factor_sep, int = int_sep)
  attr(covs, "co.names") <- co.names
  
  colnames(covs) <- names(co.names)
  
  covs
}
.get_C2 <- function(covs = NULL, int = FALSE, poly = 1, addl = NULL, distance = NULL,
                    treat = NULL, cluster = NULL, drop = TRUE, factor_sep = "_",
                    int_sep = " * ", ...) {
  #gets C data.frame, which contains all variables for which balance is to be assessed. Used in balance.table.
  if (inherits(covs, "processed_C")) {
    return(covs)
  }
  
  if (is_null(covs)) {
    drop <- FALSE
  }
  
  .chk_string(factor_sep)
  .chk_string(int_sep)
  
  #Process int and poly
  .chk_whole_number(poly)
  .chk_gte(poly, 1)
  poly <- round(poly)
  
  if (is.numeric(int)) {
    if (!chk::vld_whole_number(int) || !chk::vld_gt(int, 1)) {
      .err("`int` must be TRUE, FALSE, or a numeric (integer) value greater than 1")
    }
    
    if (int > poly) {
      poly <- int
    }
    
    int <- TRUE
  }
  else {
    .chk_flag(int)
  }
  
  center <- ...get("center", getOption("cobalt_center", default = FALSE))
  .chk_flag(center)
  orth <- ...get("orth", getOption("cobalt_orth", default = FALSE))
  .chk_flag(orth)
  
  co.names <- attr(covs, "co.names")
  seps <- attr(co.names, "seps")
  
  if (is_not_null(addl)) {
    addl.co.names <- attr(addl, "co.names")
    
    same.name <- names(addl.co.names) %in% names(co.names)
    addl <- addl[, !same.name, drop = FALSE]
    addl.co.names[same.name] <- NULL
    
    #Remove variables in addl that are redundant with covs
    if (drop && getOption("cobalt_remove_perfect_col", max(ncol(addl), ncol(covs)) <= 900)) {
      redundant.var.indices <- find_perfect_col(addl, covs)
      if (is_not_null(redundant.var.indices)) {
        addl <- addl[,-redundant.var.indices, drop = FALSE]
        addl.co.names[redundant.var.indices] <- NULL
      }
    }
    
    covs <- cbind(covs, addl)
    co.names <- c(co.names, addl.co.names)
  } 
  
  #Drop colinear with treat
  if (drop) {
    test.treat <- is_not_null(treat) && get.treat.type(treat) != "continuous"
    
    # test.cluster <- is_not_null(cluster) && !all_the_same(cluster, na.rm = FALSE)
    
    drop_vars <- vapply(seq_col(covs), 
                        function(i) {
                          # if (all_the_same(covs[,i], na.rm = FALSE)) return(TRUE)
                          test.treat && !anyNA(covs[,i]) && equivalent.factors2(covs[,i], treat)
                        }, logical(1L))
    
    if (any(drop_vars)) {
      covs <- covs[, !drop_vars, drop = FALSE]
      co.names[drop_vars] <- NULL
    }
  }
  
  C_list <- list(C = covs)
  co_list <- list(C = co.names)
  rm(co.names)
  
  if (int || (poly > 1)) {
    nsep <- 1L
    
    #Exclude NA and ints from interactions and poly
    exclude <- vapply(co_list[["C"]],
                      function(x) any(c("na", "isep") %in% x[["type"]]),
                      logical(1L))
    
    new <- .int_poly_f2(C_list[["C"]], ex = exclude, int = int, poly = poly, center = center, 
                        orth = orth, sep = rep.int(seps["int"], nsep), co.names = co_list[["C"]])
    
    C_list[["int.poly"]] <- new
    co_list[["int.poly"]] <- attr(new, "co.names")
    names(co_list[["int.poly"]]) <- vapply(co_list[["int.poly"]], 
                                           function(x) paste(x[["component"]], collapse = ""),
                                           character(1L))
  }
  
  #Drop 0 category of 0/1 variables and rename 1 category
  if (drop) {
    drop_0_1 <- rep.int(NA, length(co_list[["C"]]))
    
    for (i in seq_along(co_list[["C"]])) {
      if (!is.na(drop_0_1[i])) {
        next
      }
      
      if ("isep" %in% co_list[["C"]][[i]][["type"]] || "fsep" %nin% co_list[["C"]][[i]][["type"]]) {
        drop_0_1[i] <- FALSE
        next
      }
      
      which_are_buddies <- which(vapply(co_list[["C"]], function(j) "isep" %nin% j[["type"]] && 
                                          "fsep" %in% j[["type"]] &&
                                          j[["component"]][j[["type"]] == "base"][1L] == co_list[["C"]][[i]][["component"]][co_list[["C"]][[i]][["type"]] == "base"][1L], 
                                        logical(1L)))
      
      buddies <- co_list[["C"]][which_are_buddies]
      
      if (is_not_null(cluster)) {
        #Remove variables perfectly redundant with cluster
        unsplit_var <- unsplitfactor(as.data.frame(C_list[["C"]][, which_are_buddies, drop = FALSE]),
                                     buddies[[1L]][["component"]][buddies[[1L]][["type"]] == "base"],
                                     sep = attr(co_list[["C"]], "seps")["factor"])[[1L]]
        tab <- table(cluster, unsplit_var)
        tab <- tab[rowSums(tab == 0) < ncol(tab),,
                   drop = FALSE]
        
        if (all(rowSums(tab > 0) == 1)) {
          drop_0_1[which_are_buddies] <- TRUE
          next
        }
      }
      
      if (length(buddies) > 2L) {
        drop_0_1[which_are_buddies] <- FALSE
        next
      }
      
      buddy_is_0 <- vapply(buddies,
                           function(x) x[["component"]][x[["type"]] == "level"] %in% c("0", "FALSE"),
                           logical(1L))
      
      buddy_is_1 <- vapply(buddies,
                           function(x) x[["component"]][x[["type"]] == "level"] %in% c("1", "TRUE"),
                           logical(1L))
      
      if (!all(buddy_is_0 | buddy_is_1)) {
        drop_0_1[which_are_buddies] <- c(TRUE, FALSE)
        next
      }
      
      drop_0_1[which_are_buddies[buddy_is_0]] <- TRUE
      drop_0_1[which_are_buddies[buddy_is_1]] <- FALSE
      
      buddy_1 <- which_are_buddies[buddy_is_1]
      co_list[["C"]][[buddy_1]][["component"]] <- co_list[["C"]][[buddy_1]][["component"]][co_list[["C"]][[buddy_1]][["type"]] == "base"][1L]
      co_list[["C"]][[buddy_1]][["type"]] <- "base"
    }
    
    if (any(drop_0_1)) {
      C_list[["C"]] <- C_list[["C"]][, !drop_0_1, drop = FALSE]
      co_list[["C"]][drop_0_1] <- NULL
    }
  }
  
  if (is_not_null(co_list[["C"]])) {
    names(co_list[["C"]]) <- vapply(co_list[["C"]], function(x) paste(x[["component"]], collapse = ""), character(1L))
  }
  
  if (is_not_null(distance)) {
    if (anyNA(distance, recursive = TRUE)) {
      .err("missing values are not allowed in the distance measure")
    }
    
    distance.co.names <- attr(distance, "co.names")
    
    same.name <- names(distance.co.names) %in% unlist(lapply(co_list, names))
    if (any(same.name)) {
      distance <- distance[,!same.name, drop = FALSE]
      distance.co.names[same.name] <- NULL
    }
    
    unique.distance.names <- unique(names(distance.co.names))
    distance <- distance[,unique.distance.names, drop = FALSE]
    distance.co.names <- distance.co.names[unique.distance.names]
    
    C_list[["distance"]] <- distance
    co_list[["distance"]] <- distance.co.names
  }
  
  # C_list <- clear_null(C_list)
  # co_list <- clear_null(co_list)
  
  #Remove duplicate & redundant variables
  if (drop) {
    for (x in setdiff(names(C_list), "distance")) {
      
      #Remove self-redundant variables
      # if (getOption("cobalt_remove_perfect_col", ncol(C_list[[x]]) <= 900)) {
      #     
      #     redundant.var.indices <- find_perfect_col(C_list[[x]])
      #     if (is_not_null(redundant.var.indices)) {
      #         C_list[[x]] <- C_list[[x]][,-redundant.var.indices, drop = FALSE]
      #         co_list[[x]][redundant.var.indices] <- NULL
      #     }
      # }
      if (x != "C") {
        #Remove variables in C that have same name as other variables
        dups <- names(co_list[["C"]]) %in% co_list[[x]]
        if (any(dups)) {
          C_list[["C"]] <- C_list[["C"]][, !dups, drop = FALSE]
          co_list[["C"]][dups] <- NULL 
        }
        #Remove variables in C that are redundant with current piece
        # if (getOption("cobalt_remove_perfect_col", max(ncol(C_list[[x]]), ncol(C_list[["C"]])) <= 900)) {
        #     redundant.var.indices <- find_perfect_col(C_list[["C"]], C_list[[x]])
        #     if (is_not_null(redundant.var.indices)) {
        #         C_list[["C"]] <- C_list[["C"]][,-redundant.var.indices, drop = FALSE]
        #         co_list[["C"]][redundant.var.indices] <- NULL
        #     }
        # }
      }
    }
  }
  
  C <- do.call("cbind", clear_null(C_list[c("distance", "C", "int.poly")]))
  
  
  if (is_null(C)) {
    C <- matrix(0, nrow = length(treat), ncol = 0L,
                dimnames = list(rownames(covs), NULL))
  }
  else {
    co.names <- do.call("c", co_list[c("distance", "C", "int.poly")])
    
    for (i in seq_along(co.names)) {
      co.names[[i]]$component[co.names[[i]]$type == "fsep"] <- factor_sep
      co.names[[i]]$component[co.names[[i]]$type == "isep"] <- int_sep
    }
    
    seps["factor"] <- factor_sep
    seps["int"] <- int_sep
    
    colnames(C) <- names(co.names) <- vapply(co.names, function(x) paste(x[["component"]], collapse = ""), character(1L))
    
    
    attr(co.names, "seps") <- seps
    
    attr(C, "co.names") <- co.names
    
    attr(C, "missing.ind") <- colnames(C)[vapply(co.names, function(x) "na" %in% x[["type"]], logical(1L))]
    
    if ("distance" %in% names(C_list)) {
      attr(C, "distance.names") <- names(co_list[["distance"]])
    }
    
    attr(C, "var_types") <- .get_types(C)
  }
  
  set_class(C, "processed_C", .replace = FALSE)
}
.int_poly_f2 <- function(mat, ex = NULL, int = FALSE, poly = 1, center = FALSE,
                         orth = FALSE, sep = " * ", co.names = NULL) {
  #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
  #Only to be used in base.bal.tab; for general use see int.poly()
  #mat=matrix input
  #ex=names of variables to exclude in interactions and polynomials; a subset of df
  #int=whether to include interactions or not; currently only 2-way are supported
  #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
  #orth=use orthogonal polynomials
  #nunder=number of underscores between variables
  
  cn <- is_not_null(co.names)
  if (is_null(ex)) ex <- rep.int(FALSE, ncol(mat))
  d <- mat
  
  binary.vars <- is_binary_col(d)
  interaction.vars <- {
    if (cn) vapply(colnames(d), function(x) "isep" %in% co.names[[x]][["type"]], logical(1L))
    else rep(FALSE, ncol(d))
  }
  
  if (center && !all(binary.vars)) {
    d[, !binary.vars] <- center(d[, !binary.vars, drop = FALSE])
  }
  nd <- NCOL(d)
  
  if (poly > 1) {
    poly_terms <- poly_co.names <- make_list(poly - 1L)
    no.poly <- binary.vars | interaction.vars | ex
    npol <- nd - sum(no.poly)
    
    if (npol > 0) {
      for (i in seq_len(poly)[-1L]) {
        poly_terms[[i - 1L]] <- apply(d[, !no.poly, drop = FALSE], 2L, function(x) {
          if (orth) poly(x, degree = poly)[,i] else x^i
        })
        poly_co.names[[i - 1L]] <- {
          if (cn) {
            lapply(colnames(d)[!no.poly], function(x) {
              list(component = c(co.names[[x]][["component"]], num_to_superscript(i)), 
                   type = c(co.names[[x]][["type"]], "power"))
            })
          }
          else {
            paste0(colnames(d)[!no.poly], num_to_superscript(i))
          }
        }
      }
    }
  }
  else {
    poly_terms <- poly_co.names <- list()
  }
  
  if (int && nd > 1) {
    int_terms <- int_co.names <- make_list(1L)
    ints_to_make <- utils::combn(colnames(d)[!ex], 2L, simplify = FALSE)
    
    #Don't make ints out of ints that already exist
    # ints_that_already_exist <- get_ints_from_co.names(co.names[interaction.vars])
    # ints_to_make[vapply(ints_to_make, function(x) {
    #   any(vapply(ints_that_already_exist, function(y) {
    #     identical(sort(x), sort(y))
    #   }, logical(1L)))
    # }, logical(1L))] <- NULL
    
    #Don't make ints out of multiple members of the same categorical variable
    ints_to_make[vapply(ints_to_make, function(x) {
      "fsep" %in% co.names[[x[1L]]][["type"]] && 
        "fsep" %in% co.names[[x[2L]]][["type"]] &&
        identical(co.names[[x[1L]]][["component"]][co.names[[x[1L]]][["type"]] == "base"],
                  co.names[[x[2L]]][["component"]][co.names[[x[2L]]][["type"]] == "base"])
    }, logical(1L))] <- NULL
    
    int_terms[[1L]] <- do.call("cbind", lapply(ints_to_make, function(i) d[,i[1L]] * d[,i[2L]]))
    
    int_co.names[[1L]] <- {
      if (cn) lapply(ints_to_make, function(x) list(component = c(co.names[[x[1L]]][["component"]], sep, co.names[[x[2L]]][["component"]]),
                                                    type = c(co.names[[x[1L]]][["type"]], "isep", co.names[[x[2L]]][["type"]])))
      else vapply(ints_to_make, paste, character(1L), collapse = sep)
    }
  }
  else {
    int_terms <- int_co.names <- list()
  }
  
  out <- do.call("cbind", c(poly_terms, int_terms))
  out_co.names <- c(do.call("c", poly_co.names), do.call("c", int_co.names))
  
  if (cn) {
    names(out_co.names) <- vapply(out_co.names, 
                                  function(x) paste(x[["component"]], collapse = ""), character(1L))
    colnames(out) <- names(out_co.names)
  }
  else {
    colnames(out) <- unlist(out_co.names)
  }
  
  #Remove single values
  single_value <- apply(out, 2L, all_the_same)
  out <- out[, !single_value, drop = FALSE]
  if (cn && is_not_null(out)) {
    attr(out, "co.names") <- out_co.names[!single_value]
  }
  
  out
}
co.cbind <- function(..., deparse.level = 1) {
  if (...length() == 0L) {
    return(NULL)
  }
  
  if (...length() == 1L) {
    return(...elt(1L))
  }
  
  args <- clear_null(list(...))
  if (length(args) <= 1L) {
    return(args[[1L]])
  }
  
  co.names.list <- lapply(args, attr, "co.names")
  
  seps <- attr(co.names.list[[1L]], "seps")
  
  out <- do.call("cbind", args)
  
  attr(out, "co.names") <- do.call("c", co.names.list)
  attr(attr(out, "co.names"), "seps") <- seps
  colnames(out) <- names(attr(out, "co.names")) <- vapply(attr(out, "co.names"), function(x) paste(x[["component"]], collapse = ""), character(1L))
  
  out
}
co.rbind <- function(..., deparse.level = 1) {
  if (...length() == 0L) {
    return(NULL)
  }
  
  if (...length() == 1L) {
    return(...elt(1L))
  }
  
  args <- clear_null(list(...))
  if (length(args) <= 1L) {
    return(args[[1L]])
  }
  
  co.names <- attr(args[[1L]], "co.names")
  
  out <- do.call("rbind", args)
  
  attr(out, "co.names") <- co.names
  colnames(out) <- names(attr(out, "co.names")) <- vapply(attr(out, "co.names"), function(x) paste(x[["component"]], collapse = ""), character(1L))
  
  out
}

.get_types <- function(C) {
  vapply(colnames(C), function(x) {
    if (any(attr(C, "distance.names") == x)) "Distance"
    else if (all_the_same(C[,x]) || is_binary(C[,x]))  "Binary"
    else "Contin."
  }, character(1L))
}
find_perfect_col <- function(C1, C2 = NULL, fun = stats::cor) {
  
  #Finds indices of redundant vars in C1.
  C1.no.miss <- C1[,colnames(C1) %nin% attr(C1, "missing.ind"), drop = FALSE]
  if (is_null(C2)) {
    .use <- if (anyNA(C1)) "pairwise.complete.obs" else "everything"
    C.cor <- suppressWarnings(fun(C1.no.miss, use = .use))
    s <- !lower.tri(C.cor, diag = TRUE) & !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
  }
  else {
    C2.no.miss <- C2[,colnames(C2) %nin% attr(C2, "missing.ind"), drop = FALSE]
    .use <- if (anyNA(C1) || anyNA(C2)) "pairwise.complete.obs" else "everything"
    C.cor <- suppressWarnings(fun(C2.no.miss, y = C1.no.miss, use = .use))
    s <- !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
  }
  
  which(colSums(s) > 0)
}

model.frame2 <- function(formula, data = NULL, na.action = "na.pass", ...) {
  data <- withCallingHandlers(force(data),
                              error = function(e) .err(conditionMessage(e)),
                              warning = function(w) .wrn(conditionMessage(w)))
  
  tryCatch({
    stats::model.frame(formula, data = data, na.action = na.action, ...)
  },
  error = function(e) {
    ee <- conditionMessage(e)
    if (startsWith(ee, "object '") && endsWith(ee, "' not found")) {
      v <- sub("object '([^']+)' not found", "\\1", ee)
      .err(sprintf("the variable %s cannot be found. Be sure it is entered correctly or supply a dataset that contains this varialble to `data`",
                   add_quotes(v, 2L)))
    }
    
    .err(ee)
  })
}

#base.bal.tab
check_if_zero_weights <- function(weights.df, treat = NULL) {
  #Checks if all weights are zero in each treat group for each set of weights
  if (is_not_null(treat)) {
    w.t.mat <- expand.grid(weight_names = colnames(weights.df), 
                           treat_vals = treat_vals(treat), 
                           stringsAsFactors = FALSE)
    if (NROW(w.t.mat) > 0L) {
      problems <- vapply(seq_row(w.t.mat),
                         function(x) all(check_if_zero(weights.df[treat == w.t.mat[x, "treat_vals"], w.t.mat[x, "weight_names"]])),
                         logical(1L))
      
      if (any(problems)) {
        prob.w.t.mat <- w.t.mat[problems,]
        if (NCOL(weights.df) == 1L) {
          error <- sprintf("All weights are zero when treat is %s.",
                           word_list(prob.w.t.mat[, "treat_vals"], "or", quotes = 2L))
        }
        else {
          errors <- vapply(unique(prob.w.t.mat[,"weight_names"]), function(i) {
            sprintf("%s weights are zero when treat is %s",
                    add_quotes(i),
                    word_list(prob.w.t.mat[prob.w.t.mat[,"weight_names"] == i, "treat_vals"], "or",
                              quotes = 2L))
          }, character(1L))
          errors <- paste(c("All", rep.int("all", length(errors) - 1L)), errors)
          error <- paste0(word_list(errors, "and"), ".")
        }
        .err(error, tidy = FALSE)
      }
    }
  }
  else if (is_not_null(colnames(weights.df))) {
    problems <- vapply(colnames(weights.df),
                       function(wn) all(check_if_zero(weights.df[, wn])),
                       logical(1L))
    
    if (any(problems)) {
      prob.wts <- colnames(weights.df)[problems]
      if (NCOL(weights.df) == 1L) {
        error <- "All weights are zero."
      }
      else {
        errors <- vapply(prob.wts, function(i) {
          sprintf("%s weights are zero", add_quotes(i))
        }, character(1L))
        errors <- paste(c("All", rep.int("all", length(errors) - 1L)), errors)
        error <- paste0(word_list(errors, "and"), ".")
      }
      .err(error, tidy = FALSE)
    }
  }
}
.baltal <- function(threshold) {
  #threshold: vector of threshold values (i.e., "Balanced"/"Not Balanced")
  threshnames <- names(table(threshold))
  balstring <- threshnames[nzchar(threshnames)][1L]
  thresh.val <- substring(balstring, 1L + regexpr("[><]", balstring), nchar(balstring))
  
  b <- data.frame(count = c(sum(threshold == sprintf("Balanced, <%s", thresh.val)), 
                            sum(threshold == sprintf("Not Balanced, >%s", thresh.val))))
  
  rownames(b) <- c(sprintf("Balanced, <%s", thresh.val),
                   sprintf("Not Balanced, >%s", thresh.val))
  
  b
}
.max_imbal <- function(balance.table, col.name, thresh.col.name, abs_stat) {
  balance.table.clean <- balance.table[balance.table$Type != "Distance" & is.finite(balance.table[, col.name]),]
  maxed <- balance.table.clean[which.max(abs_stat(balance.table.clean[, col.name])), match(c(col.name, thresh.col.name), names(balance.table.clean))]
  
  cbind(Variable = rownames(maxed), as.data.frame(maxed))
}
threshold_summary <- function(compute, thresholds, no.adj, balance.table, weight.names = NULL, agg.fun = NULL) {
  out <- do.call("c", lapply(compute, function(s) make_list(paste.(c("Balanced", "Max.Imbalance"), s))))
  
  for (s in compute) {
    if (is_null(thresholds[[s]])) {
      out[[paste.("Balanced", s)]] <- NULL
      out[[paste.("Max.Imbalance", s)]] <- NULL
    }
    else if (no.adj) {
      out[[paste.("Balanced", s)]] <- .baltal(balance.table[[paste.(STATS[[s]]$Threshold, "Un")]])
      out[[paste.("Max.Imbalance", s)]] <- .max_imbal(balance.table[balance.table[["Type"]] != "Distance", , drop = FALSE], 
                                                      col.name = {
                                                        if (is_null(agg.fun)) paste.(STATS[[s]]$bal.tab_column_prefix, "Un")
                                                        else paste.(firstup(agg.fun), STATS[[s]]$bal.tab_column_prefix, "Un")
                                                      }, 
                                                      thresh.col.name = paste.(STATS[[s]]$Threshold, "Un"), 
                                                      abs_stat = STATS[[s]]$abs)
    }
    else if (length(weight.names) == 1L) {
      out[[paste.("Balanced", s)]] <- .baltal(balance.table[[STATS[[s]]$Threshold]])
      out[[paste.("Max.Imbalance", s)]] <- .max_imbal(balance.table[balance.table[["Type"]] != "Distance", , drop = FALSE], 
                                                      col.name = {
                                                        if (is_null(agg.fun)) paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")
                                                        else paste.(firstup(agg.fun), STATS[[s]]$bal.tab_column_prefix, "Adj")
                                                      }, 
                                                      thresh.col.name = STATS[[s]]$Threshold, 
                                                      abs_stat = STATS[[s]]$abs)
    }
    else if (length(weight.names) > 1L) {
      out[[paste.("Balanced", s)]] <- setNames(do.call("cbind", lapply(weight.names, function(x) .baltal(balance.table[[paste.(STATS[[s]]$Threshold, x)]]))),
                                               weight.names)
      out[[paste.("Max.Imbalance", s)]] <- cbind(Weights = weight.names,
                                                 do.call("rbind", lapply(weight.names, function(x) setNames(.max_imbal(balance.table[balance.table[["Type"]] != "Distance", , drop = FALSE], 
                                                                                                                       col.name = {
                                                                                                                         if (is_null(agg.fun)) paste.(STATS[[s]]$bal.tab_column_prefix, x)
                                                                                                                         else paste.(firstup(agg.fun), STATS[[s]]$bal.tab_column_prefix, x)
                                                                                                                       },  
                                                                                                                       thresh.col.name = paste.(STATS[[s]]$Threshold, x), 
                                                                                                                       abs_stat = STATS[[s]]$abs),
                                                                                                            c("Variable", 
                                                                                                              STATS[[s]]$bal.tab_column_prefix, 
                                                                                                              STATS[[s]]$Threshold)))),
                                                 stringsAsFactors = FALSE)
    }
  }
  
  out
}

balance_table <- function(C, type, weights = NULL, treat, continuous, binary, s.d.denom, 
                          thresholds = list(), un = FALSE, disp = NULL, stats = NULL, 
                          s.weights = rep.int(1, length(treat)), abs = FALSE, no.adj = FALSE, 
                          var_types = NULL, s.d.denom.list = NULL, quick = TRUE, ...) {
  #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
  weight.names <- if (no.adj) "Adj" else names(weights)
  
  if (is_not_null(s.d.denom.list)) names(s.d.denom.list) <- weight.names
  if (is_not_null(s.d.denom)) names(s.d.denom) <- weight.names
  
  disp <- c(disp, all_STATS(type)[all_STATS(type) %in% stats])
  compute <- if (quick) disp else c("means", "sds", all_STATS(type)[all_STATS(type) %in% stats])
  
  #B=Balance frame
  Bnames <- c("Type", 
              expand_grid_string(c(switch(type,
                                          "bin" = expand_grid_string(c("M"["means" %in% compute],
                                                                       "SD"["sds" %in% compute]),
                                                                     c("0", "1"), collapse = "."),
                                          "cont" = c("M"["means" %in% compute], "SD"["sds" %in% compute])), 
                                   unlist(lapply(intersect(compute, all_STATS(type)[!get_from_STATS("adj_only")[get_from_STATS("type") == type]]), function(s) {
                                     c(STATS[[s]]$bal.tab_column_prefix,
                                       if (no.adj && is_not_null(thresholds[[s]])) STATS[[s]]$Threshold)
                                   }))
              ),
              "Un", collapse = "."),
              expand_grid_string(c(switch(type,
                                          "bin" = expand_grid_string(c("M"["means" %in% compute], "SD"["sds" %in% compute]),
                                                                     c("0", "1"), collapse = "."),
                                          "cont" = c("M"["means" %in% compute], "SD"["sds" %in% compute])), 
                                   unlist(lapply(intersect(compute, all_STATS(type)), function(s) {
                                     c(STATS[[s]]$bal.tab_column_prefix,
                                       if (!no.adj && is_not_null(thresholds[[s]])) STATS[[s]]$Threshold)
                                   }))
              ),
              weight.names, collapse = "."))
  B <- make_df(Bnames, NCOL(C))
  rownames(B) <- colnames(C)
  
  #Set var type (binary/continuous)
  B[["Type"]] <- if_null_then(var_types, .get_types(C))
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
    sd.computable <- if (binary == "std") rep.int(TRUE, nrow(B)) else !bin.vars
    if (type == "bin") {
      tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
      
      if (un || !quick) {
        for (t in c("0", "1")) {
          sds <- rep.int(NA_real_, NCOL(C))
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
            sds <- rep.int(NA_real_, NCOL(C))
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
        sds <- rep.int(NA_real_, NCOL(C))
        if (any(sd.computable)) {
          sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE],
                                         weights = NULL, s.weights = s.weights,
                                         bin.vars = bin.vars[sd.computable])
        }
        B[["SD.Un"]] <- sds
      }
      if (!no.adj && (!quick || "sds" %in% disp)) {
        for (i in weight.names) {
          sds <- rep.int(NA_real_, NCOL(C))
          if (any(sd.computable)) {
            sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE],
                                           weights = weights[[i]], s.weights = s.weights,
                                           bin.vars = bin.vars[sd.computable])
          }
          B[[paste.("SD", i)]] <- sds
        }
      }
    }
    
    if (all_apply(B[startsWith(names(B), "SD.")], function(x) !any(is.finite(x)))) {
      disp <- disp[disp != "sds"]
    }
  }
  
  for (s in all_STATS(type)) {
    if (s %in% compute) {
      if (!get_from_STATS("adj_only")[s] && (!quick || un)) {
        B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]] <- STATS[[s]]$fun(C, treat = treat, weights = NULL, 
                                                                              std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                              s.d.denom = if_null_then(s.d.denom.list[[1L]], s.d.denom[1L]),
                                                                              abs = abs, s.weights = s.weights, bin.vars = bin.vars,
                                                                              weighted.weights = weights[[1L]], ...)
      }
      
      if (!no.adj && (!quick || s %in% disp)) {
        for (i in weight.names) {
          B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]] <- STATS[[s]]$fun(C, treat = treat, weights = weights[[i]],
                                                                             std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                                                             s.d.denom = if_null_then(s.d.denom.list[[i]], s.d.denom[i]),
                                                                             abs = abs, s.weights = s.weights, bin.vars = bin.vars, ...)
        }
      }
      
      
      if (all_apply(intersect(names(B), paste.(STATS[[s]]$bal.tab_column_prefix, c("Un", weight.names))), 
                    function(x) !any(is.finite(B[[x]])))) {
        disp <- disp[disp != s] 
        thresholds[[s]] <- NULL
      }
      
      if (is_not_null(thresholds[[s]])) {
        if (!get_from_STATS("adj_only")[s] && no.adj) {
          B[[paste.(STATS[[s]]$Threshold, "Un")]] <- ifelse_(B[["Type"]] == "Distance" | !is.finite(B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]]), "", 
                                                             STATS[[s]]$abs(B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3L)),
                                                             paste0("Not Balanced, >", round(thresholds[[s]], 3L)))
        }
        else if (!no.adj) {
          for (i in weight.names) {
            B[[paste.(STATS[[s]]$Threshold, i)]] <- ifelse_(B[["Type"]] == "Distance" | !is.finite(B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]]), "", 
                                                            STATS[[s]]$abs(B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3L)),
                                                            paste0("Not Balanced, >", round(thresholds[[s]], 3L)))
          }
        }
      }
      
      if (no.adj || NCOL(weights) <= 1L) {
        names(B)[names(B) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
      }
    }
  }
  
  attr(B, "thresholds") <- thresholds
  attr(B, "disp") <- disp
  attr(B, "compute") <- compute
  
  B
}

samplesize <- function(treat, type, weights = NULL, subclass = NULL, s.weights = NULL,
                       method = c("matching", "weighting", "subclassification"), discarded = NULL) {
  #Computes sample size info. for unadjusted and adjusted samples.
  # method is what method the weights are to be used for. 
  # method="subclassification" is for subclass sample sizes only.
  
  if (is_null(s.weights)) s.weights <- rep.int(1, length(treat))
  if (is_null(discarded)) discarded <- rep.int(FALSE, length(treat))
  
  if (type == "bin") {
    if (length(method) == 1L && method == "subclassification") {
      if (is_null(subclass)) {
        .err("`subclass` must be a vector of subclasses")
      }
      
      nn <- make_df(c(levels(subclass), "Discarded", "All"), c(treat_names(treat), "Total"))
      
      nn[["All"]] <- c(vapply(treat_vals(treat), function(tn) sum(treat == tn), numeric(1L)), length(treat))
      nn[["Discarded"]] <- {
        if (any(discarded)) c(vapply(treat_vals(treat), function(tn) sum(treat[discarded] == tn), numeric(1L)), length(treat))
        else NULL
      }
      
      matched <- !is.na(subclass)
      for (k in levels(subclass)) {
        tk <- treat[matched & subclass == k]
        nn[[k]] <- c(vapply(treat_vals(treat), function(tn) sum(tk == tn), numeric(1L)),
                     length(tk))
      }
      
      for (tnn in names(treat_names(treat))) {
        small.subclass <- nn[treat_names(treat)[tnn], levels(subclass)] <= 1L
        if (any(small.subclass))
          .wrn(sprintf("Not enough %s units in subclass%s %s",
                       tnn,
                       ngettext(sum(small.subclass), "", "es"), 
                       word_list(levels(subclass)[small.subclass])))
      }
      attr(nn, "tag") <- "Sample sizes by subclass"
    }
    else {
      if (is_null(weights)) {
        nn <- make_df(treat_names(treat), "All")
        nn["All", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat == tn]), numeric(1L))
        if (nunique.gt(s.weights, 2L) || !any(s.weights == 1) || !all(s.weights %in% c(0, 1))) {
          attr(nn, "ss.type") <- c("ess")
        }
        else {
          attr(nn, "ss.type") <- c("ss")
        }
      }
      else if (NCOL(weights) == 1L) {
        if (method == "matching") {
          nn <- make_df(treat_names(treat), c("All (ESS)", "All (Unweighted)", "Matched (ESS)",
                                              "Matched (Unweighted)", "Unmatched", "Discarded"))
          nn["All (ESS)", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat == tn]), numeric(1L))
          nn["All (Unweighted)", ] <- vapply(treat_vals(treat), function(tn) sum(treat == tn & s.weights > 0), numeric(1L))
          nn["Matched (ESS)", ] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat == tn, 1L] * s.weights[treat == tn]), numeric(1L))
          nn["Matched (Unweighted)", ] <- vapply(treat_vals(treat), function(tn) sum(treat == tn & weights[,1L] > 0 & s.weights > 0), numeric(1L))
          nn["Unmatched", ] <- vapply(treat_vals(treat), function(tn) sum(treat == tn & weights[,1L] == 0 & !discarded), numeric(1L))
          nn["Discarded", ] <- vapply(treat_vals(treat), function(tn) sum(treat == tn & discarded), numeric(1L))
          
          attr(nn, "ss.type") <- rep.int("ss", NROW(nn))
          
          if (!any(discarded)) {
            attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
            nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
          }
        }
        else if (method == "weighting") {
          nn <- make_df(treat_names(treat), c("Unadjusted", "Adjusted", "Discarded"))
          nn["Unadjusted", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat == tn]), numeric(1L))
          nn["Adjusted", ] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat == tn, 1L] * s.weights[treat == tn]), numeric(1L))
          nn["Discarded", ] <- vapply(treat_vals(treat), function(tn) sum(treat == tn & discarded), numeric(1L))
          attr(nn, "ss.type") <- c("ss", "ess", "ss")
          
          if (!any(discarded)) {
            attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
            nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
          }
        }
      }
      else {
        nn <- make_df(treat_names(treat), c("All", names(weights)))
        nn["All", ] <- vapply(treat_vals(treat), function(tn) ESS(s.weights[treat == tn]), numeric(1L))
        for (i in seq_col(weights)) {
          nn[1L + i,] <- vapply(treat_vals(treat), function(tn) ESS(weights[treat == tn, i] * s.weights[treat == tn]), numeric(1L))
        }
        attr(nn, "ss.type") <- c("ss", rep.int("ess", length(method)))
      }
      
      attr(nn, "tag") <- {
        if (length(attr(nn, "ss.type")) > 1L && all(attr(nn, "ss.type")[-1L] == "ess")) {
          "Effective sample sizes"
        }
        else "Sample sizes"
      }
    }
  }
  else if (type == "cont") {
    if (length(method) == 1L && method == "subclassification") {
      if (is_null(subclass)) {
        .err("`subclass` must be a vector of subclasses")
      }
      
      nn <- make_df(c(levels(subclass), "All"), c("Total"))
      
      nn[, "All"] <- length(treat)
      
      matched <- !is.na(subclass)
      for (k in levels(subclass)) {
        nn[[k]] <- sum(matched & subclass == k)
      }
      small.subclass <- nn[, levels(subclass)] <= 1L
      if (any(small.subclass))
        .wrn(sprintf("not enough units in subclass%s %s",
                     ngettext(sum(small.subclass), "", "es"), 
                     word_list(levels(subclass)[small.subclass])))
      attr(nn, "tag") <- "Sample sizes by subclass"
    }
    else {
      if (is_null(weights)) {
        nn <- make_df("Total", "All")
        nn["All", ] <- ESS(s.weights)
        if (nunique.gt(s.weights, 2L) || !any(s.weights == 1) || !all(s.weights %in% c(0, 1))) {
          attr(nn, "ss.type") <- c("ess")
        }
        else {
          attr(nn, "ss.type") <- c("ss")
        }
        
      }
      else if (NCOL(weights) == 1L) {
        if (method == "matching") {
          nn <- make_df("Total", c("All (ESS)", "All (Unweighted)", "Matched (ESS)", "Matched (Unweighted)", "Unmatched", "Discarded"))
          nn["All (ESS)", ] <- ESS(s.weights)
          nn["All (Unweighted)", ] <- sum(s.weights > 0)
          nn["Matched (ESS)", ] <- ESS(weights[, 1L] * s.weights)
          nn["Matched (Unweighted)", ] <- sum(weights[, 1L] > 0 & s.weights > 0 & !discarded)
          nn["Unmatched", ] <- sum(weights[, 1L] == 0 & !discarded)
          nn["Discarded", ] <- sum(discarded)
          
          attr(nn, "ss.type") <- rep.int("ss", NROW(nn))
          
          if (!any(discarded)) {
            attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
            nn <- nn[rownames(nn) != "Discarded",, drop = FALSE]
          }
        }
        else if (method == "weighting") {
          nn <- make_df("Total", c("Unadjusted", "Adjusted", "Discarded"))
          nn["Unadjusted", ] <- ESS(s.weights)
          nn["Adjusted", ] <- ESS(weights[!discarded, 1L] * s.weights[!discarded])
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
        
        for (i in seq_col(weights)) {
          nn[1L + i,] <- switch(method[i],
                                "matching" = ESS(weights[!discarded, i]),
                                "weighting" = ESS(weights[!discarded, i] * s.weights[!discarded]))
          
        }
        
        attr(nn, "ss.type") <- c("ss", rep.int("ess", length(method)))
      }
      
      attr(nn, "tag") <- {
        if (length(attr(nn, "ss.type")) > 1L && all(attr(nn, "ss.type")[-1L] == "ess")) {
          "Effective sample sizes"
        }
        else "Sample sizes"
      }
    }
  }
  
  nn
}

balance_summary <- function(bal.tab.list, agg.funs, include.times = FALSE) {
  type <- attr(bal.tab.list[[1L]], "print.options")[["type"]]
  disp <- attr(bal.tab.list[[1L]], "print.options")[["disp"]]
  compute <- attr(bal.tab.list[[1L]], "print.options")[["compute"]]
  thresholds <- attr(bal.tab.list[[1L]], "print.options")[["thresholds"]]
  quick <- attr(bal.tab.list[[1L]], "print.options")[["quick"]]
  weight.names <- if_null_then(attr(bal.tab.list[[1L]], "print.options")[["weight.names"]], "Adj")
  abs <- attr(bal.tab.list[[1L]], "print.options")[["abs"]]
  no.adj <- attr(bal.tab.list[[1L]], "print.options")[["nweights"]] == 0
  
  balance.list <- clear_null(grab(bal.tab.list, "Balance"))
  
  Brownames <- unique(unlist(lapply(balance.list, rownames), use.names = FALSE))
  
  Agg.Funs <- firstup(if (quick) agg.funs else c("min", "mean", "max"))
  Agg.Funs.Given <- firstup(agg.funs)
  
  if (length(Agg.Funs) == 1L && Agg.Funs == "Max") {
    abs <- TRUE
  }
  
  Bcolnames <- c("Times", "Type", 
                 paste.(c(unlist(lapply(compute[compute %in% all_STATS(type)[!get_from_STATS("adj_only")[get_from_STATS("type") == type]]], function(s) {
                   c(paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix),
                     if (no.adj && is_not_null(thresholds[[s]]) && length(Agg.Funs.Given) == 1L) STATS[[s]]$Threshold)
                 }))), "Un"),
                 unlist(lapply(weight.names, function(wn) {
                   paste.(c(unlist(lapply(compute[compute %in% all_STATS(type)], function(s) {
                     c(paste.(Agg.Funs, STATS[[s]]$bal.tab_column_prefix),
                       if (!no.adj && is_not_null(thresholds[[s]]) && length(Agg.Funs.Given) == 1L) STATS[[s]]$Threshold)
                   }))), wn)
                 })))
  
  B <- make_df(Bcolnames, Brownames)
  
  B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(vapply(balance.list, function(y) if (x %in% rownames(y)) y[[x, "Type"]] else NA_character_, character(1L))))), use.names = FALSE)
  
  B[["Times"]] <- {
    if (include.times) vapply(Brownames, function(x) toString(seq_along(balance.list)[vapply(balance.list, function(y) x %in% rownames(y), logical(1L))]), character(1L))[Brownames]
    else NULL
  }
  
  for (Agg.Fun in Agg.Funs) {
    for (s in compute[compute %in% all_STATS(type)]) {
      abs0 <- function(x) {
        if (is_null(x)) NA_real_
        else if (abs) STATS[[s]]$abs(x)
        else x
      }
      
      agg <- function(x, ...) {
        if (!any(is.finite(x))) NA_real_
        else if (s == "variance.ratios" && tolower(Agg.Fun) == "mean") .geom_mean(x)
        else if (tolower(Agg.Fun) == "rms") sqrt(mean_fast(STATS[[s]]$abs(x)^2, TRUE))
        else get(tolower(Agg.Fun))(x, ...)
      }
      
      for (sample in c("Un", weight.names)) {
        if ((sample == "Un" || !no.adj) && (sample != "Un" || !get_from_STATS("adj_only")[s])) {
          B[[paste.(Agg.Fun, STATS[[s]]$bal.tab_column_prefix, sample)]] <- vapply(Brownames, function(x) agg(unlist(lapply(balance.list, function(y) if (x %in% rownames(y)) abs0(y[[x, paste.(STATS[[s]]$bal.tab_column_prefix, sample)]]))), na.rm = TRUE), numeric(1L))
        }
      }
    }
  }
  
  if (length(Agg.Funs.Given) == 1L) {
    #Assign X.Threshold values
    for (s in compute[compute %in% all_STATS(type)]) {
      if (is_not_null(thresholds[[s]])) {
        if (!get_from_STATS("adj_only")[s] && no.adj) {
          B[[paste.(STATS[[s]]$Threshold, "Un")]] <- ifelse_(B[["Type"]] == "Distance" | !is.finite(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, "Un")]]), "", 
                                                             STATS[[s]]$abs(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, "Un")]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3L)),
                                                             paste0("Not Balanced, >", round(thresholds[[s]], 3)))
        }
        else {
          for (i in weight.names) {
            B[[paste.(STATS[[s]]$Threshold, i)]] <- ifelse_(B[["Type"]] == "Distance" | !is.finite(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, i)]]), "", 
                                                            STATS[[s]]$abs(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix, i)]]) < thresholds[[s]], paste0("Balanced, <", round(thresholds[[s]], 3L)),
                                                            paste0("Not Balanced, >", round(thresholds[[s]], 3)))
          }
        }
      }
      
      if (no.adj || length(weight.names) <= 1L) {
        names(B)[names(B) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
      }
    }
  }
  
  B
}


#base.bal.tab.imp
samplesize_across_imps <- function(obs.list) {
  obs.list <- clear_null(obs.list)
  
  obs <- Reduce("+", obs.list) / length(obs.list)
  attr(obs, "tag") <- sprintf("Average %s across imputations",
                              tolower(attr(obs.list[[1L]], "tag")))
  obs
}

#base.bal.tab.multi
samplesize_multi <- function(bal.tab.multi.list, treat_names, focal = NULL) {
  which <- {
    if (is_null(focal)) treat_names
    else c(setdiff(treat_names, focal), focal)
  }
  
  bal.tab.multi.list <- clear_null(bal.tab.multi.list)
  obs <- do.call("cbind", unname(grab(bal.tab.multi.list, "Observations")))[, which]
  attr(obs, "tag") <- attr(bal.tab.multi.list[[1L]][["Observations"]], "tag")
  attr(obs, "ss.type") <- attr(bal.tab.multi.list[[1L]][["Observations"]], "ss.type")
  
  obs
}

#base.bal.tab.cluster
samplesize_across_clusters <- function(obs.list) {
  obs.list <- clear_null(obs.list)
  obs <- Reduce("+", obs.list)
  attr(obs, "tag") <- sprintf("Total %s across clusters",
                              tolower(attr(obs.list[[1L]], "tag")))
  
  obs
}

#base.bal.tab.subclass
balance_table_subclass <- function(C, type, weights = NULL, treat, subclass,
                                   continuous, binary, s.d.denom, 
                                   thresholds = list(), un = FALSE, disp = NULL, stats = NULL, 
                                   s.weights = rep.int(1, length(treat)), abs = FALSE, 
                                   var_types = NULL, quick = TRUE, ...) {
  #Creates list SB of balance tables for each subclass
  #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
  
  disp <- unique(c(disp, intersect(all_STATS(type), stats)))
  compute <- if (quick) disp else c("means", "sds", all_STATS(type))
  
  #B=Balance frame
  Bnames <- c("Type", 
              expand_grid_string(c(switch(type,
                                          "bin" = expand_grid_string(c("M"["means" %in% compute],
                                                                       "SD"["sds" %in% compute]),
                                                                     c("0", "1"), collapse = "."),
                                          "cont" = c("M"["means" %in% compute], "SD"["sds" %in% compute])), 
                                   unlist(lapply(intersect(compute, all_STATS(type)), function(s) {
                                     c(STATS[[s]]$bal.tab_column_prefix,
                                       if (is_not_null(thresholds[[s]])) STATS[[s]]$Threshold)
                                   }))
              ),
              c("Adj"), collapse = "."))
  B <- make_df(Bnames, colnames(C))
  
  #Set var type (binary/continuous)
  B[["Type"]] <- if_null_then(var_types, .get_types(C))
  bin.vars <- B[["Type"]] == "Binary"
  
  SB <- setNames(rep(list(B), nlevels(subclass)), levels(subclass))
  
  binary <- match_arg(binary, c("raw", "std"))
  sd.computable <- if (binary == "std") rep.int(TRUE, nrow(B)) else !bin.vars
  
  subclass_w_empty <- {
    if (type == "bin") 
      vapply(levels(subclass), function(i) {
        any_apply(treat_vals(treat), function(t) !any(treat == t & subclass == i))
      }, logical(1L))
    else
      vapply(levels(subclass), function(i) !any(subclass == i), logical(1L))
  }
  
  for (i in levels(subclass)) {
    
    in.subclass <- !is.na(subclass) & subclass == i
    
    #Means for each group
    if ("means" %in% compute) {
      if (type == "bin") {
        tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
        for (t in c("0", "1")) {
          SB[[i]][[paste.("M", t, "Adj")]] <- col_w_mean(C, subset = treat == tn01[t] & in.subclass, s.weights = s.weights)
        }
      }
      else if (type == "cont") {
        SB[[i]][["M.Adj"]] <- col_w_mean(C, subset = in.subclass, s.weights = s.weights)
      }
    }
    
    #SDs for each group
    if ("sds" %in% compute) {
      if (type == "bin") {
        tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
        for (t in c("0", "1")) {
          sds <- rep.int(NA_real_, NCOL(C))
          sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == tn01[t] & in.subclass,
                                         s.weights = s.weights)
          SB[[i]][[paste.("SD", t, "Adj")]] <- sds
        }
      }
      else if (type == "cont") {
        sds <- rep.int(NA_real_, NCOL(C))
        sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == in.subclass,
                                       s.weights = s.weights)
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
        
        if (all_apply(SB[[i]][paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")], 
                      function(x) !any(is.finite(x)))) {
          disp <- disp[disp != s] 
          thresholds[[s]] <- NULL
        }
        
        
        if (is_not_null(thresholds[[s]])) {
          
          SB[[i]][[paste.(STATS[[s]]$Threshold, "Adj")]] <- ifelse(SB[[i]][["Type"]] != "Distance" & is.finite(SB[[i]][[paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")]]), 
                                                                   paste0(ifelse(abs_(SB[[i]][[paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")]]) < thresholds[[s]], "Balanced, <", "Not Balanced, >"), round(thresholds[[s]], 3L)), "")
        }
        names(SB[[i]])[names(SB[[i]]) == paste.(STATS[[s]]$Threshold, "Adj")] <- STATS[[s]]$Threshold
      }
      
    }
    
  }
  
  if (all_apply(SB, function(sb) all_apply(sb[startsWith(names(sb), "SD.")], function(x) !any(is.finite(x))))) {
    disp <- disp[disp != "sds"]
  }
  
  attr(SB, "thresholds") <- thresholds
  attr(SB, "disp") <- disp
  attr(SB, "compute") <- compute
  
  SB
}

# !!! NEEDS TO BE UPDATED !!!
balance_table_across_subclass_cont <- function(balance.table, balance.table.subclass.list, subclass.obs, r.threshold = NULL) {
  #NEEDS TO BE UPDATED
  
  B.A <- balance.table.subclass.list[[1L]][names(balance.table.subclass.list[[1L]]) %in% c("M.Adj", "SD.Adj", "Corr.Adj")]
  
  for (i in rownames(B.A)) {
    for (j in colnames(B.A)) {
      B.A[[i, j]] <- {
        if (startsWith(j, "SD."))
          sqrt(sum(vapply(seq_along(balance.table.subclass.list),
                          function(s) subclass.obs[[s]] / sum(subclass.obs) * (balance.table.subclass.list[[s]][[i, j]]^2),
                          numeric(1L))))
        else
          sum(vapply(seq_along(balance.table.subclass.list),
                     function(s) subclass.obs[[s]] / sum(subclass.obs) * (balance.table.subclass.list[[s]][[i, j]]),
                     numeric(1L)))
      }
    }
  }
  
  B.A.df <- cbind(balance.table[c("Type", "M.Un", "SD.Un", "Corr.Un", "R.Threshold.Un")], 
                  B.A, R.Threshold = NA_character_)
  
  if (is_not_null(r.threshold)) {
    B.A.df[["R.Threshold"]] <- ifelse(B.A.df[["Type"]] == "Distance", "", paste0(ifelse(is.finite(B.A.df[["Corr.Adj"]]) & abs_(B.A.df[["Corr.Adj"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold))
  }
  
  B.A.df
}

#Misc
`%+%` <- function(...) {
  if (is.atomic(..1) && is.atomic(..2)) crayon::`%+%`(as.character(..1), as.character(..2))
  else ggplot2::`%+%`(...)
}

check_arg_lengths <- function(...) {
  dots_names <- vapply(match.call(expand.dots = FALSE)$..., deparse1,
                       character(1L))
  
  lens <- setNames(integer(...length()), dots_names)
  for (i in seq_along(lens)) {
    lens[i] <- len(...elt(i))
  }
  
  supplied <- lens > 0L
  if (!all_the_same(lens[supplied])) {
    .err(sprintf("%s must have the same number of units",
                 word_list(dots_names[supplied], quotes = "`")))
  }
}

intapprox <- function(f, from, to, steps, method = "midpoint") {
  method <- match_arg(method, c("midpoint", "trapezoidal", "simpsons"))
  
  seg <- seq(from, to, length = steps)
  delta <- seg[2L] - seg[1L]
  
  if (method == "midpoint") {
    mids <- (seg[-1L] + seg[-steps]) / 2
    s <- sum(f(mids)) * delta
  }
  else if (method == "trapezoidal") {
    s <- (f(from) + 2 * sum(f(seg[-c(1L, steps)])) + f(to)) * delta / 2
  }
  else if (method == "simpsons") {
    mids <- (seg[-1L] + seg[-steps]) / 2
    sm <- sum(f(mids)) * delta
    st <- (f(from) + 2 * sum(f(seg[-c(1L, steps)])) + f(to)) * delta / 2
    s <- (2 * sm + st) / 3
  }
  
  s
}