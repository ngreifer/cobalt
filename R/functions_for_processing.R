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
  
  arg::arg_supplied(treat)
  
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
process_treat.list <- function(treat.list, ...) {
  arg::arg_supplied(treat.list)
  
  if (!is.list(treat.list)) {
    treat.list <- as.list(treat.list)
  }
  
  hasdots <- ...length() > 0L
  
  treat.list.names <- vapply(seq_along(treat.list), function(ti) {
    if (hasdots && rlang::is_string(treat.list[[ti]])) treat.list[[ti]]
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
  .attr(treat, "treat_names")
}
`treat_vals<-` <- function(treat, value) {
  `attr<-`(treat, "treat_vals", value)
}
treat_vals <- function(treat) {
  .attr(treat, "treat_vals")
}
subset_processed.treat <- function(x, index) {
  y <- x[index]
  treat_names(y) <- treat_names(x)[treat_vals(x) %in% unique(y)]
  treat_vals(y) <- treat_vals(x)[treat_vals(x) %in% unique(y)]
  
  assign.treat.type(y) |>
    set_class(class(x))
}

#Every slot of `X`, in the order `bal.tab()` returns them. `.finish_X()` fills each
#from the like-named local in the calling `x2base()` method, so a slot missing here is
#a slot silently dropped.
X_SLOTS <- c("covs", "treat", "weights", "distance", "addl", "s.d.denom", "estimand",
             "call", "cluster", "imp", "s.weights", "focal", "discarded", "method",
             "subclass", "stats", "thresholds")

X_MSM_SLOTS <- c("covs.list", "treat.list", "weights", "distance.list", "addl.list",
                 "s.d.denom", "call", "cluster", "imp", "s.weights", "focal",
                 "discarded", "method", "subclass", "stats", "thresholds")

initialize_X <- function() {
  make_list(X_SLOTS)
}
initialize_X_msm <- function() {
  make_list(X_MSM_SLOTS)
}
.weight_check <- function(w) {
  wname <- deparse1(substitute(w))
  
  if (!is.list(w)) {
    w <- list(w)
  }
  
  if (anyNA(w, recursive = TRUE)) {
    arg::err("{.val {NA}}s are not allowed in the {wname}")
  }
  
  for (x in w) {
    if (!all(is.numeric(x))) {
      arg::err("all {wname} must be numeric")
    }
    
    if (!all(is.finite(x))) {
      arg::err("infinite {wname} are not allowed")
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
    arg::wrn("some clusters have only one unit in them, which may yield unexpected results")
  }
  
  if (stop_warn["bw"]) {
    arg::wrn("some clusters have only one member of a treatment group in them, which may yield unexpected results")
  }
  
  if (stop_warn["bs"]) {
    arg::err("not all treatment levels are present in all clusters")
  }
  
}
strata2weights <- function(strata, treat, estimand = NULL, focal = NULL) {
  #Process strata into weights (similar to weight.subclass from MatchIt)
  
  #Checks
  arg::arg_vector(strata)
  
  #Process treat
  treat <- process_treat(treat)
  
  if (get.treat.type(treat) == "continuous") {
    arg::err("{.arg strata} cannot be turned into weights for continuous treatments")
  }
  
  s.d.denom <- .get_s.d.denom(NULL, estimand = estimand, subclass = strata, treat = treat,
                              focal = focal, quietly = TRUE)
  
  focal <- if (s.d.denom %in% treat_vals(treat)) process_focal(s.d.denom, treat)
  
  NAsub <- is.na(strata)
  imat <- do.call("cbind", lapply(treat_vals(treat), function(t) treat == t & !NAsub))
  colnames(imat) <- treat_vals(treat)
  
  weights <- rep_with(0, treat)
  
  if (!is.factor(strata)) {
    strata <- factor(strata, nmax = min(colSums(imat)))
    levels(strata) <- seq_len(nlevels(strata))
  }
  
  t_by_sub <- do.call("rbind", lapply(treat_vals(treat), function(t) tabulate(strata[imat[, t]], nlevels(strata))))
  dimnames(t_by_sub) <- list(treat_vals(treat), levels(strata))
  
  total_by_sub <- colSums(t_by_sub)
  
  strata.c <- as.character(strata)
  
  if (is_not_null(focal)) {
    focal <- process_focal(focal, treat)
    for (t in treat_vals(treat)) {
      weights[imat[, t]] <- {
        if (t == focal)  1
        else (t_by_sub[focal, ] / t_by_sub[t, ])[strata.c[imat[, t]]]
      }
    }
  }
  else {
    for (t in treat_vals(treat)) {
      weights[imat[, t]] <- (total_by_sub / t_by_sub[t, ])[strata.c[imat[, t]]]
    }
  }
  
  na.w <- !is.finite(weights)
  if (any(na.w)) {
    weights[na.w] <- 0
    arg::wrn("some units were given weights of zero due to zeros in stratum membership")
  }
  
  if (all(check_if_zero(weights))) {
    arg::err("no units were stratified")
  }
  
  for (tnn in names(treat_names(treat))) {
    if (all(check_if_zero(weights[treat == treat_vals(treat)[treat_names(treat)[tnn]]]))) {
      arg::err("no {tnn} units were stratified")
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
          if (is_null(D)) as.data.frame(covs)
          else cbind(D, as.data.frame(covs)) 
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
        arg::err("if {.arg covs} is a character vector, {.arg data} must be specified as a data frame")
      }
      
      if (!all(covs %in% colnames(data))) {
        arg::err("all entries in {.arg covs} must be names of variables in {.arg data}")
      }
      
      covs_c <- try(get_covs_from_formula(f.build(covs), data = as.data.frame(data)), silent = TRUE)
    }
    else {
      arg::err("{.arg covs} must be a data.frame of covariates")
    }
  }
  
  if (is_not_null(treat)) {
    treat_c <- try({
      if (!is.atomic(treat)) arg::err("{.arg treat} must be a vector of treatment statuses")
      if (rlang::is_string(treat)) get_treat_from_formula(reformulate(".", treat), data = data)
      else get_treat_from_formula(treat ~ ., treat = treat)
    }, silent = TRUE)
  }
  
  covs_to_use <- treat_to_use <- "c"
  if (is_error(covs_c)) {
    if (is_error(covs_f)) {
      arg::err("{(.attr(covs_c, 'condition')$message)}")
    }
    covs_to_use <- "f"
  }
  else if (is_null(covs_c)) {
    if (is_error(covs_f)) {
      arg::err("{(.attr(covs_f, 'condition')$message)}")
    }
    if (is_null(covs_f) && needs.covs) {
      arg::err("no covariates were specified")
    }
    covs_to_use <- "f"
  }
  
  if (is_error(treat_c)) {
    if (is_error(treat_f)) {
      arg::err("{(.attr(treat_c, 'condition')$message)}")
    }
    treat_to_use <- "f"
  }
  else if (is_null(treat_c)) {
    if (is_error(treat_f)) {
      arg::err("{(.attr(treat_f, 'condition')$message)}")
    }
    if (is_null(treat_f) && needs.treat) {
      arg::err("no treatment variable was specified")
    }
    treat_to_use <- "f"
  }
  
  t.c <- list(treat = switch(treat_to_use, "f" = treat_f, "c" = treat_c), 
              covs = switch(covs_to_use, "f" = covs_f, "c" = covs_c))
  
  t.c$treat.name <- .attr(t.c$treat, "treat.name")
  
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
      arg::err("the argument supplied to {.arg {i}} must be a named list of weights, names of variables containing weights in an available data set, or objects with a {.fun get.w} method")
    }
    else {
      arg::err("the argument supplied to {.arg {i}} must be a vector, a data frame, or the names of variables in an available data set")
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
    arg::wrn("names were provided to {.arg {i}}, but no argument to {.arg data} was provided. Ignoring {.arg {i}}")
    
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
  
  arg::wrn("the following {cli::qty(sum(not.found))} variable{?s} named in {.arg {i}} {?is/are} not in any available data sets and will be ignored: {.var {val[not.found]}}")
  
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
    arg::err("all entries in {.arg {i}} must have names")
  }
  
  if (!all_the_same(vapply(val.list, NROW, numeric(1L)))) {
    arg::err("all items in {.arg {i}} must have the same length")
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
    arg::err("the argument to {.arg {i}} must be a list of the same length as the number of time points in {call.phrase}")
  }
  
  for (ti in which(lengths(val.List) > 0L)) {
    val <- val.List[[ti]]
    val.df <- NULL
    
    if (is.list(val) && !is.data.frame(val)) {
      val.list <- lapply(val, .process_val, strsplit(i, ".list", fixed = TRUE)[[1L]],
                         treat.list[[ti]], covs.list[[ti]], addl.data = addl.data)
      
      if (!all_the_same(vapply(val.list, NROW, numeric(1L)))) {
        arg::err("all items in {.arg {i}} must have the same length")
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
      arg::err("missing values must not exist in {.arg {i}}")
    }
    
    val.List[[ti]] <- val.df
  }
  
  val.df.lengths <- val.List |>
    clear_null() |>
    vapply(nrow, numeric(1L))
  
  if (!all_the_same(val.df.lengths)) {
    arg::err("all columns in {.arg {i}} must have the same number of rows")
  }
  
  val.List
}
.process_vector <- function(vec, name = deparse1(substitute(vec)), which = name, datalist = list(), missing.okay = FALSE) {
  bad.vec <- FALSE
  if (rlang::is_string(vec) && is_not_null(datalist)) {
    for (i in seq_along(datalist)) {
      if (is.matrix(datalist[[i]]) && vec %in% colnames(datalist[[i]])) {
        vec <- datalist[[i]][, vec]
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
    arg::err("the argument to {.arg {name}} must be a vector of {which} or the (quoted) name of a variable in {.arg data} that contains {which}")
  }
  
  if (!missing.okay && anyNA(vec)) {
    arg::err("missing values must not exist in {.arg {name}}")
  }
  
  vec
} 
#The standardization factor for a binary or multi-category treatment.
#
#Four sources are consulted in order -- an explicit `s.d.denom`, the `estimand`, the
#focal group, and the weights themselves -- each falling through to the next when it
#has nothing to say. This used to be expressed with three boolean flags threaded
#between blocks; it is now the chain it always was, so each step can be read on its
#own and the note at the end knows whether anything was inferred.
.get_s.d.denom <- function(s.d.denom = NULL, estimand = NULL, weights = NULL,
                           subclass = NULL, treat = NULL, focal = NULL,
                           quietly = FALSE) {
  #A denominator already resolved by a wrapper is not re-derived by its children.
  if (isTRUE(.attr(s.d.denom, "checked"))) {
    return(s.d.denom)
  }

  n.weights <- NCOL(weights)
  specified <- is_not_null(s.d.denom)

  .treat <- function() {
    if (inherits(treat, "processed.treat")) treat
    else process_treat(treat)
  }

  out <- {
    if (specified) .s.d.denom_from_arg(s.d.denom, n.weights, .treat(), focal)
    else if (is_not_null(estimand)) {
      .s.d.denom_from_estimand(estimand, n.weights, .treat(), quietly)
    }
    else NULL
  }

  #An explicit `s.d.denom` or `estimand` naming the focal group defers to it, as does
  #having neither.
  inferred <- is_null(out)

  if (inferred) {
    treat <- .treat()

    out <- {
      if (is_not_null(focal) && n.weights <= 1L) focal
      else .s.d.denom_from_weights(weights, subclass, treat)
    }
  }

  if (is_not_null(weights)) {
    if (length(out) == 1L) {
      out <- rep.int(out, n.weights)
    }

    if (length(out) != n.weights) {
      arg::err("valid inputs to {.arg s.d.denom} or {.arg estimand} must have length 1 or equal to the number of valid sets of weights, which is {n.weights}")
    }

    names(out) <- names(weights)
  }

  if (!quietly) {
    .s.d.denom_note(out, specified, inferred, weights, .treat())
  }

  attr(out, "checked") <- TRUE

  out
}

#An explicit `s.d.denom`. Returns NULL when it names the focal group, which the caller
#resolves.
.s.d.denom_from_arg <- function(s.d.denom, n.weights, treat, focal) {
  if (length(s.d.denom) > 1L && length(s.d.denom) != n.weights) {
    arg::err("{.arg s.d.denom} must have length 1 or equal to the number of valid sets of weights, which is {n.weights}")
  }

  #These four need nothing from the treatment, so a value made only of them is taken
  #as given.
  allowable <- c("pooled", "all", "weighted", "hedges")

  if (all(s.d.denom %in% allowable)) {
    return(s.d.denom)
  }

  binary <- length(treat_names(treat)) == 2L &&
    all(c("treated", "control") %in% names(treat_names(treat)))

  allowable <- c(allowable, c("treated", "control")[binary], "focal"[is_not_null(focal)])

  s.d.denom <- arg::match_arg(s.d.denom,
                              unique(c(as.character(treat_vals(treat)), allowable)),
                              several.ok = n.weights > 1L)

  t.c <- s.d.denom %in% c("treated", "control")

  if (any(t.c)) {
    s.d.denom[t.c] <- treat_vals(treat)[treat_names(treat)[s.d.denom[t.c]]]
    return(s.d.denom)
  }

  if (any(s.d.denom == "focal")) {
    return(NULL)
  }

  s.d.denom
}

#The `estimand`. Returns NULL when the estimand names no group of its own -- an
#unrecognized value, or an ATT or ATC on a treatment with no "treated" or "control".
.s.d.denom_from_estimand <- function(estimand, n.weights, treat, quietly = FALSE) {
  allowable <- c("ATT", "ATC", "ATE", "ATO", "ATM")

  estimand <- tryCatch(arg::match_arg(estimand, allowable, several.ok = TRUE),
                       error = function(cond) NA_character_)

  if (anyNA(estimand)) {
    #Ignoring an unrecognized estimand silently would substitute a different
    #standardization factor without saying so: a typo such as "ATTT" would read as ATE.
    if (!quietly) {
      arg::wrn("{.arg estimand} should be {.or {.str {allowable}}}; ignoring it")
    }

    return(NULL)
  }

  #An ATT or ATC with a multi-category treatment is legitimate; the focal group
  #determines the denominator instead.
  if (any(estimand %in% c("ATT", "ATC")) && get.treat.type(treat) != "binary") {
    return(NULL)
  }

  if (length(estimand) > 1L && length(estimand) != n.weights) {
    arg::err("{.arg estimand} must have length 1 or equal to the number of valid sets of weights, which is {n.weights}")
  }

  vapply(estimand, switch, character(1L),
         ATT = treat_vals(treat)[treat_names(treat)["treated"]],
         ATC = treat_vals(treat)[treat_names(treat)["control"]],
         ATO = "weighted",
         ATM = "weighted",
         "pooled")
}

#What the weights or subclasses imply, when nothing was specified: the group left
#unweighted is the one the others were weighted toward.
.s.d.denom_from_weights <- function(weights, subclass, treat) {
  if (is_not_null(subclass)) {
    #The group whose share is most even across subclasses is the one subclassification
    #held fixed.
    sub.tab <- table(treat, subclass)[treat_vals(treat), ]
    sub.tab <- rbind(sub.tab, table(subclass)[colnames(sub.tab)])
    dimnames(sub.tab) <- list(c(treat_vals(treat), "pooled"), colnames(sub.tab))

    evenness <- apply(sub.tab, 1L, function(x) .mean_abs_dev(x) / sum(x))

    return(rownames(sub.tab)[which.min(evenness)])
  }

  if (is_null(weights)) {
    return("pooled")
  }

  vapply(weights, function(w) {
    for (tv in treat_vals(treat)) {
      if (all_the_same(w[treat == tv]) && !all_the_same(w[treat != tv])) {
        return(tv)
      }
    }

    "pooled"
  }, character(1L))
}

#Says what was assumed, when anything was.
.s.d.denom_note <- function(s.d.denom, specified, inferred, weights, treat) {
  #`weighted` needs weights to mean anything. Without them `.compute_s.d.denom()`
  #reaches the same factor as `all`, which is what this reports.
  if (specified && is_null(weights) && any(s.d.denom == "weighted")) {
    arg::msg("note: {.arg s.d.denom} specified as {.str weighted}, but no weights supplied; setting to {.str all}")
    return(invisible())
  }

  #Reading a treatment group off the weights is unambiguous, so it goes unremarked;
  #falling back to `pooled` or similar is a guess and is announced.
  if (!inferred || all(s.d.denom %in% treat_vals(treat))) {
    return(invisible())
  }

  if (all_the_same(s.d.denom)) {
    arg::msg("note: {.arg s.d.denom} not specified; assuming {.str {s.d.denom[1L]}}")
    return(invisible())
  }

  #With several sets of weights, name each one -- reporting a treatment value by the
  #role it plays where the treatment is 0/1.
  as.role <- function(s) {
    if (s %nin% treat_vals(treat) || !all(treat_vals(treat) %in% c("0", "1"))) {
      return(s)
    }

    names(treat_names(treat))[treat_names(treat) ==
                                names(treat_vals(treat))[treat_vals(treat) == s]]
  }

  wt_strs <- sprintf("{.str %s} for {.var %s}",
                     vapply(s.d.denom, as.role, character(1L)),
                     names(weights)) |>
    vapply(cli::format_inline, character(1L))

  arg::msg("note: {.arg s.d.denom} not specified; assuming {wt_strs}")
}
.get_s.d.denom.cont <- function(s.d.denom, weights = NULL, subclass = NULL, quietly = FALSE) {
  s.d.denom.specified <- !missing(s.d.denom) && is_not_null(s.d.denom)
  
  if (s.d.denom.specified) {
    if (isTRUE(.attr(s.d.denom, "checked"))) {
      return(s.d.denom)
    }
    
    if (length(s.d.denom) > 1L && length(s.d.denom) != NCOL(weights)) {
      arg::err("{.arg s.d.denom} must have length 1 or equal to the number of valid sets of weights")
    }
    
    allowable.s.d.denoms <- {
      if (is_not_null(subclass)) "all"
      else c("all", "weighted")
    }
    
    if (!all(s.d.denom %in% allowable.s.d.denoms)) {
      s.d.denom <- arg::match_arg(s.d.denom, unique(allowable.s.d.denoms), 
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
      arg::err("valid inputs to {.arg s.d.denom} or {.arg estimand} must have length 1 or equal to the number of valid sets of weights")
    }
    
    names(s.d.denom) <- names(weights)
  }
  
  if (!quietly) {
    if (s.d.denom.specified && is_null(weights) && any(s.d.denom == "weighted")) {
      arg::msg("note: {.arg s.d.denom} specified as {.str weighted}, but no weights supplied; setting to {.str all}")
    }
  }

  attr(s.d.denom, "checked") <- TRUE
  
  s.d.denom
}
#The `continuous`/`binary` defaults. `type == "cont"` and
#`get.treat.type(treat) == "continuous"` are the same condition -- `process_stats()`
#sets the type attribute from the treatment type -- so one expression serves both
#the leaf and the wrapper call sites, which spelled it four different ways.
.get_std_defaults <- function(treat, continuous = NULL, binary = NULL) {
  list(continuous = continuous %or% getOption("cobalt_continuous", "std"),
       binary = binary %or% getOption("cobalt_binary",
                                      switch(get.treat.type(treat),
                                             "continuous" = "std", "raw")))
}

#Resolve the standardization factor, or `NULL` when nothing will be standardized.
#
#Returns the value rather than assigning it, because the leaf and the wrappers want
#different things when it comes back `NULL`: the leaf clears `X$s.d.denom`, while a
#wrapper must keep whatever the user supplied so that each per-stratum child does
#not re-resolve it from scratch -- and possibly differently, since a covariate that
#is continuous overall can be constant within one stratum.
.resolve_s.d.denom <- function(X, var_types, continuous, binary) {
  #A multi-category treatment precomputes one denominator per weight set across the
  #whole sample and passes it down as `s.d.denom.list`; that supersedes this.
  if (is_not_null(X[["s.d.denom.list"]])) {
    return(NULL)
  }

  any_std <- (binary == "std" && any(var_types == "Binary")) ||
    (continuous == "std" && !all(var_types == "Binary"))

  if (!any_std) {
    return(NULL)
  }

  if (!any(vapply(X[["stats"]], function(s) STATS[[s]]$needs_s.d.denom, logical(1L)))) {
    return(NULL)
  }

  if (get.treat.type(X[["treat"]]) == "continuous") {
    return(.get_s.d.denom.cont(X[["s.d.denom"]], weights = X[["weights"]], subclass = X[["subclass"]]))
  }

  .get_s.d.denom(X[["s.d.denom"]], estimand = X[["estimand"]], weights = X[["weights"]],
                 subclass = X[["subclass"]], treat = X[["treat"]], focal = X[["focal"]])
}
#The two verdicts a threshold produces, rounded to three places independently of
#`print()`'s `digits`. `.baltal()` counts against these same strings, so the tally's
#row names cannot drift from the labels in the table.
.threshold_verdicts <- function(threshold) {
  c(paste0("Balanced, <", round(threshold, 3L)),
    paste0("Not Balanced, >", round(threshold, 3L)))
}

#The verdict for one statistic column. Distance rows and non-finite values get an empty
#label rather than a verdict.
.threshold_label <- function(x, var_types, threshold, abs_stat) {
  verdicts <- .threshold_verdicts(threshold)

  ifelse_(var_types == "Distance" | !is.finite(x), "",
          abs_stat(x) < threshold, verdicts[1L], verdicts[2L])
}
.compute_s.d.denom <- function(mat, treat = NULL, s.d.denom = "pooled", s.weights = NULL,
                               bin.vars = NULL, subset = NULL, weighted.weights = NULL,
                               to.sd = rep.int(TRUE, ncol(mat)), na.rm = TRUE) {
  denoms <- setNames(rep.int(1, ncol(mat)), colnames(mat))
  
  if (rlang::is_string(s.d.denom)) {
    if (is_null(bin.vars)) {
      bin.vars <- rep.int(FALSE, ncol(mat))
      bin.vars[to.sd] <- is_binary_col(mat[subset, to.sd, drop = FALSE])
    }
    else if (!is.atomic(bin.vars) || length(bin.vars) != ncol(mat) ||
             anyNA(as.logical(bin.vars))) {
      arg::err("{.arg bin.vars} must be a logical vector with length equal to the number of columns of {.arg mat}")
    }
    
    possibly.supplied <- c("mat", "treat", "weighted.weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, numeric(1L)),
                        possibly.supplied)
    supplied <- lengths > 0L
    if (!all_the_same(lengths[supplied])) {
      arg::err("{.arg {possibly.supplied[supplied]}} must have the same number of units")
    }
    
    if (lengths["weighted.weights"] == 0L) {
      weighted.weights <- rep.int(1, NROW(mat))
    }
    
    if (lengths["s.weights"] == 0L) {
      s.weights <- rep.int(1, NROW(mat))
    }
    
    if (lengths["subset"] == 0L) {
      subset <- rep.int(TRUE, NROW(mat))
    }
    else if (anyNA(as.logical(subset))) {
      arg::err("{.arg subset} must be a logical vector")
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
      arg::err("{.arg s.d.denom} is not an allowed value")
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
      arg::err("{.arg s.d.denom} must be an allowable value or a numeric vector of with length equal to the number of columns of {.arg mat}. See {.help [?col_w_smd](cobalt::col_w_smd)} for allowable values")
    }
  }
  else {
    arg::err("{.arg s.d.denom} must be an allowable value or a numeric vector of with length equal to the number of columns of {.arg mat}. See {.help [?col_w_smd](cobalt::col_w_smd)} for allowable values")
  }
  
  denoms
}
.assign_X_class <- function(X) {
  X <- clear_null(X)
  
  if (is_not_null(X[["treat"]]) && !has.treat.type(X[["treat"]])) {
    X[["treat"]] <- assign.treat.type(X[["treat"]])
  }
  
  #`cluster` is checked before `subclass` so that clustered subclassification is
  #handled by nesting the subclass tables within each cluster, as is done for
  #multi-category treatments and imputations.
  if (is_not_null(X[["cluster"]]) && nlevels(X[["cluster"]]) > 1L) X.class <- "cluster"
  else if (is_not_null(X[["subclass"]])) {
    X.class <- switch(get.treat.type(X[["treat"]]),
                      binary = "subclass.binary",
                      continuous = "subclass.cont",
                      arg::err("multi-category treatments are not currently compatible with subclasses"))
  }
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
  else arg::err("couldn't determine length of {.arg X} components")
}
#The slots `subset_X()` knows how to subset: everything with one entry per unit. The
#others -- `estimand`, `focal`, `method`, `stats`, and so on -- describe the analysis
#rather than the units.
SUBSETTABLE_SLOTS <- c("covs", "treat", "weights", "distance", "addl", "cluster",
                       "imp", "s.weights", "discarded", "subclass", "covs.list",
                       "treat.list", "distance.list", "addl.list")
subset_X <- function(X, subset = NULL) {
  if (is_null(subset) || !any(names(X) %in% SUBSETTABLE_SLOTS)) {
    return(X)
  }
  
  n <- .get_length_X(X)
  
  if (is.logical(subset)) {
    if (length(subset) != n) {
      arg::err("{.arg subset} must have the same length as the other entries")
    }
    
    if (!any(subset)) {
      arg::err("all {.arg subset} set to {.val {FALSE}}")
    }
    
    if (all(subset)) {
      return(X)
    }
    
    subset <- which(subset)
  }
  else if (!is.numeric(subset)) {
    arg::err("{.arg subset} must be logical or numeric")
  }
  else if (max(subset) > n) {
    arg::err("subset indices cannot be higher than the length of the other entries")
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
  
  for (i in intersect(names(X), SUBSETTABLE_SLOTS)) {
    X[[i]] <- subset_X_internal(X[[i]], subset)
  }
  
  X
}
.mids_complete <- function(data) {
  arg::arg_is(data, "mids")
  
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
    if (is.integer(.attr(data$data, "row.names"))) seq_row(cmp)
    else as.character(seq_row(cmp))
  }
  
  cmp
}

#Expands anything that was supplied for a single imputation to the full stacked data,
#and checks that everything has a compatible number of observations. `objects` is a
#named list of the values to check, named after the `vectors`/`data.frames`/`lists`
#arguments that say how each should be expanded; the (possibly expanded) list comes
#back. `.finish_X()` is the only caller, and it reads and writes the method's locals
#itself, so this does not touch its caller's frame.
length_imp_process <- function(objects, vectors = NULL, data.frames = NULL,
                               lists = NULL, imp = NULL, original.call.to = NULL) {
  all.objects <- c(vectors, data.frames, lists)
  ensure.equal.lengths <- TRUE
  problematic <- rlang::rep_named(all.objects, FALSE)

  lengths <- vapply(all.objects, function(x) {
    if (x %in% lists) {
      if (is_null(objects[[x]])) 0
      else max(vapply(objects[[x]], len, numeric(1L)))
    }
    else len(objects[[x]])
  }, numeric(1L))

  #Process imp further
  if (is_not_null(imp)) {
    imp.lengths <- vapply(levels(imp), function(i) sum(imp == i), numeric(1L))

    if (all_the_same(imp.lengths)) {
      unsorted.imp <- is.unsorted(imp)

      #Repeats one imputation's worth of values across all of them, keeping each
      #imputation's block in the order `imp` gives.
      .stack <- function(x, i_) {
        new_x <- x[rep(i_, length(imp.lengths))]

        if (unsorted.imp) {
          for (i in levels(imp)) {
            new_x[imp == i] <- x
          }
        }

        new_x
      }

      .stack_rows <- function(x) {
        new_x <- x[rep(seq_row(x), length(imp.lengths)), , drop = FALSE]

        if (unsorted.imp) {
          for (i in levels(imp)) {
            new_x[imp == i, ] <- x
          }
        }

        new_x
      }

      for (i in all.objects[lengths > 0 & lengths != length(imp)]) {
        if (lengths[i] != imp.lengths[1L]) {
          problematic[i] <- TRUE
          next
        }

        objects[[i]] <- {
          if (i %in% vectors) .stack(objects[[i]], seq_along(objects[[i]]))
          else if (i %in% data.frames) .stack_rows(objects[[i]])
          else lapply(objects[[i]], function(j) {
            if (!is.factor(j) && !is_mat_like(j)) {
              arg::err("{.arg {i}} can only contain vectors or data frames")
            }

            if (is.factor(j)) .stack(j, seq_along(j))
            else .stack_rows(j)
          })
        }
      }
    }
    else {
      problematic <- lengths > 0L & lengths != length(imp)
    }

    if (any(problematic)) {
      arg::err("{.arg {names(problematic)[problematic]}} must have the same number of observations as {.arg imp}")
    }

    ensure.equal.lengths <- FALSE
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
      arg::err(sprintf("{.arg {names(problematic)[problematic]}} must have the same number of observations as in the original call to %s", original.call.to))
    }
    else {
      arg::err("{.arg {names(problematic)[problematic]}} must have the same number of observations as {.arg {anchor}}")
    }
  }

  objects
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
    if (is_null(stats)) {
      stats <- getOption("cobalt_stats", "mean.diffs")
    }
    
    stats <- unique(arg::match_arg(stats, all_STATS("bin"), several.ok = TRUE))
    attr(stats, "type") <- "bin"
  }
  else if (treat.type == "continuous") {
    if (is_null(stats)) {
      stats <- getOption("cobalt_stats", "correlations")
    }
    
    stats <- unique(arg::match_arg(stats, all_STATS("cont"), several.ok = TRUE))
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
    arg::err("{.arg thresholds} must be numeric")
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
#Resolve `stats`, `thresholds`, and `s.d.denom` from the arguments `bal.tab()`
#forwards to `x2base()`. Each `x2base` method calls this once, passing its own
#`treat` (or `treat.list`) and `...`.
#
#`...` is forwarded rather than `list(...)` so that only the elements actually
#looked up are forced, matching the inline code this replaces. The legacy
#per-statistic arguments -- `disp.<stat>` and `<stat>.threshold`, named in the
#`STATS` registry -- are among them, which is why they are read through `...get()`
#rather than declared.
#
#`bal.plot()` also goes through `x2base()` but needs no balance statistics, so it
#gets an empty list back and all three slots stay `NULL`.
#
#The first formal is `.treat`, not `treat`: `...` carries the user's own `treat`
#argument, and a named element of `...` matches a formal of the same name even when
#the caller supplied that formal positionally -- so a formal named `treat` would be
#handed the unprocessed argument (a column name, say) instead of the processed
#treatment vector.
process_stats_and_thresholds <- function(.treat, ...) {
  if (check_if_call_from_fun(bal.plot)) {
    return(list())
  }

  stats <- process_stats(...get("stats"), treat = .treat)
  type <- .attr(stats, "type")

  thresholds <- ...get("thresholds", list())

  if (is_not_null(thresholds)) {
    thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
    stats <- unique(c(stats, names(thresholds)))
  }

  for (s in all_STATS(type)) {
    #A `disp.<stat>` flag adds the statistic or removes it, falling back to the
    #corresponding global option as `process_disp()` does for `disp.means` and
    #`disp.sds`.
    disp.stat <- ...get(STATS[[s]][["disp_stat"]],
                        getOption(paste0("cobalt_", STATS[[s]][["disp_stat"]])))

    if (isTRUE(disp.stat)) {
      stats <- unique(c(stats, s))
    }
    else if (isFALSE(disp.stat)) {
      stats <- setdiff(stats, s)
    }

    #A `<stat>.threshold` overrides `thresholds` and implies the statistic.
    if (is_not_null(...get(STATS[[s]][["threshold"]]))) {
      thresholds[[s]] <- ...get(STATS[[s]][["threshold"]])
    }

    if (is_not_null(thresholds[[s]])) {
      thresholds[[s]] <- STATS[[s]][["abs"]](thresholds[[s]])

      if (between(thresholds[[s]], STATS[[s]][["threshold_range"]])) {
        stats <- unique(c(stats, s))
      }
      else {
        thresholds[[s]] <- NULL
        arg::wrn('{.arg {STATS[[s]][["threshold"]]}} must be between {STATS[[s]][["threshold_range"]]}; ignoring {.arg {STATS[[s]][["threshold"]]}}')
      }
    }
  }

  list(stats = process_stats(stats, treat = .treat),
       thresholds = thresholds,
       s.d.denom = ...get("s.d.denom"))
}
process_subset <- function(subset = NULL, n) {
  if (is_null(subset)) {
    return(NULL)
  }
  
  if (!is.logical(subset) && !is.numeric(subset)) {
    arg::err("the argument to {.arg subset} must be a logical or numeric vector")
  }
  
  if (is.numeric(subset)) {
    if (any(abs(subset) > n)) {
      arg::err("numeric values for {.arg subset} cannot be larger than the number of units")
    }
    subset <- subset[!is.na(subset) & subset != 0]

    if (any(subset < 0) && any(subset > 0)) {
      arg::err("positive and negative indices cannot be mixed with {.arg subset}")
    }

    logical.subset <- rep.int(any(subset < 0), n)
    logical.subset[abs(subset)] <- !logical.subset[abs(subset)]
    subset <- logical.subset
  }
  
  if (anyNA(subset)) {
    arg::wrn("{.val {NA}}s were present in {.arg subset}. Treating them like {.val {FALSE}}")
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
    
    arg::err("the name specified to {.arg focal} is not the name of any treatment group")
  }
  
  if (can_str2num(treat) && focal %in% str2num(treat)) {
    return(as.character(focal))
  }
  
  if (focal > length(treat_vals(treat))) {
    arg::err("{.arg focal} was specified as {focal}, but there are only {length(treat_vals(treat))} treatment groups")
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
      addl.methods <- rep.int(arg::match_arg(A[["method"]], c("weighting", "matching")), ncol(addl.weights))
    }
    else {
      addl.methods <- arg::match_arg(A[["method"]], c("weighting", "matching"), several.ok = TRUE)
      
      if (length(addl.methods) != ncol(addl.weights)) {
        arg::err("valid inputs to {.arg method} must have length 1 or equal to the number of valid sets of additional weights")
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
  arg::when_not_null(disp, arg::arg_character)
  disp <- {
    if (is_null(disp)) getOption("cobalt_disp")
    else arg::match_arg(disp, acceptable.options()[["disp"]], several.ok = TRUE)
  }
  
  for (d in c("means", "sds")) {
    if (getOption(paste.("cobalt_disp", d), FALSE)) {
      disp <- unique(c(disp, d))
    }
    
    disp.d <- ...get(paste.("disp", d))
    
    if (is_not_null(disp.d)) {
      arg::arg_flag(disp.d, .arg = sprintf("disp.%s", d))

      disp <- unique(c(disp, d[disp.d]))
      
      disp <-  {
        if (disp.d) unique(c(disp, d))
        else unique(disp[disp != d])
      }
    }
  }
  
  disp
}
process_addl <- function(addl = NULL, datalist = list()) {
  if (is_not_null(addl) && !is.atomic(addl) && !rlang::is_formula(addl) &&
      !is.matrix(addl) && !is.data.frame(addl)) {
    arg::err("{.arg addl} must be a formula or variable containing the additional covariates")
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
      arg::err("{.arg addl.list} must have an entry for each time point")
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
    arg::err("{.arg distance} must be a formula or variable containing the distance values")
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
  
  if (is_null(obj.distance) || all(is.na(obj.distance))) {
    return(distance_t.c[["covs"]])
  }

  #Callers may pass `obj.distance.name = NULL` (e.g., when the name attribute was
  #dropped while coercing the component to a data frame). Falling through to
  #`setNames(., NULL)` would leave the column unnamed and silently discard the
  #distance variable, so recover the name from the object or use the default.
  if (is_null(obj.distance.name)) {
    obj.distance.name <- colnames(as.data.frame(obj.distance)) %or% "distance"
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
      arg::err("{.arg distance} must have an entry for each time point")
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
    arg::err("the argument supplied to {.arg focal} must be the name of a level of treatment")
  }
  
  if (treat.type == "multinomial") {
    
    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      arg::wrn('{.val {estimand}} is not compatible with {.arg focal}. Setting {.arg estimand} to {.val ATT}')
      reported.estimand <- estimand <- "ATT"
    }
    
    #With more than two groups there is no single control group for an ATC to name, so
    #an ATC is an ATT against whichever group `focal` identifies. Either way that group
    #has to be named.
    if (estimand %in% c("ATT", "ATC")) {
      if (is_null(focal)) {
        if (is_null(treated) || treated %nin% unique.treat) {
          arg::err('when {.code estimand = "{reported.estimand}"} for multinomial treatments, an argument must be supplied to {.arg focal}')
        }

        focal <- treated
      }

      estimand <- "ATT"
    }
  }
  else if (treat.type == "binary") {
    unique.treat.bin <- unique(binarize(treat), nmax = 2L)
    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      arg::wrn("{.val {estimand}} is not compatible with {.arg focal}. Setting {.arg estimand} to {.val ATT}")
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
        treated <- {
          if (is.factor(treat)) levels(treat)[2L]
          else unique.treat[unique.treat.bin == 1]
        }
        
        if (estimand == "ATT") {
          arg::msg("assuming {.val {treated}} {?is/are} the treated level{?s}. If not, supply an argument to {.arg focal}")
          
        }
        else if (estimand == "ATC") {
          arg::msg("assuming {.val {setdiff(unique.treat, treated)}} {?is/are} the control level{?s}. If not, supply an argument to {.arg focal}")
        }
      }
      
      if (estimand == "ATT") {
        focal <- treated
      }
      else if (estimand == "ATC") {
        focal <- setdiff(unique.treat, treated)
      }
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
.process_stop_method <- function(sm, sm_avail) {
  if (is_null(sm)) {
    rule1 <- sm_avail
  }
  else if (any(is.character(sm))) {
    rule1 <- sm_avail[vapply(tolower(sm_avail), function(x) any(startsWith(x, tolower(sm))), logical(1L))]
    if (is_null(rule1)) {
      arg::wrn("{.arg stop.method} should be {.or {.val {sm_avail}}}. Using all available stop methods instead")
      rule1 <- sm_avail
    }
  }
  else if (is.numeric(sm) && any(sm %in% seq_along(sm_avail))) {
    if (!all(sm %in% seq_along(sm_avail))) {
      arg::wrn("there {?is/are} {length(sm_avail)} stop method{?s} available, but you requested {setdiff(sm, seq_along(sm_avail))}")
    }
    rule1 <- sm_avail[sm %in% seq_along(sm_avail)]
  }
  else {
    arg::wrn("{.arg stop.method} should be {.or {.val {sm_avail}}}. Using all available stop methods instead")
    rule1 <- sm_avail
  }
  
  rule1
}

#.get_C2
get_treat_from_formula <- function(f, data = NULL, treat = NULL) {
  
  if (is.character(f)) {
    f <- try(as.formula(f), silent = TRUE)
  }
  
  arg::arg_formula(f)
  
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
    arg::wrn("the argument supplied to {.arg data} is not a data frame. Ignoring {.arg data}")
    data <- env
    data.specified <- FALSE
  }
  
  tt <- try_arg(terms(f, data = data))
  
  if (rlang::is_formula(tt, lhs = TRUE)) {
    resp.vars.mentioned <- as.character(rlang::f_lhs(tt))
    resp.vars.failed <- vapply(resp.vars.mentioned, function(v) {
      test <- tryCatch(eval(str2expression(v), data, env), error = function(e) e)
      
      if (inherits(test, "simpleError")) {
        if (!identical(conditionMessage(test), sprintf("object '%s' not found", v))) {
          arg::err("{conditionMessage(test)}")
        }
        
        return(TRUE)
      }
      
      if (is.function(test)) {
        arg::err("invalid type (function) for variable {.var {v}}")
      }
      
      is_null(test)
    }, logical(1L))
    
    if (any(resp.vars.failed)) {
      if (is_not_null(treat)) {
        tt <- delete.response(tt)
      }
      else if (data.specified) {
        arg::err("the given response variable, {.var {resp.vars.mentioned}}, is not a variable in {.arg data} or the global environment")
      }
      else {
        arg::err("the given response variable, {.var {resp.vars.mentioned}}, is not a variable in the global environment")
      }
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
#Each row name of a terms object's `factors` matrix is a variable as it was written in
#the formula. One that will not evaluate as written may still evaluate quoted in
#backticks -- a column called `my var` does, `I(x^2)` does not -- so that is tried
#before giving up. Returns the names, backtick-quoted where that was needed.
#
#`...` is passed to `eval()`, so the caller says where to look: the data and the
#formula's environment on the first pass, the assembled model frame on the second.
.backtick_unevaluable <- function(vars, ..., reject.functions = FALSE) {
  for (i in seq_along(vars)) {
    evaled <- try(eval(str2expression(vars[i]), ...), silent = TRUE)

    if (null_or_error(evaled)) {
      quoted <- add_quotes(vars[i], "`")
      evaled <- try(eval(str2expression(quoted), ...), silent = TRUE)

      if (null_or_error(evaled)) {
        .err_unevaluable(evaled, name.missing = reject.functions)
      }

      vars[i] <- quoted
    }

    #A function used where a variable belongs would otherwise be evaluated silently.
    if (reject.functions && is.function(evaled)) {
      arg::err("invalid type (function) for variable {.var {vars[i]}}")
    }
  }

  vars
}

#`name.missing` is set only on the first pass, where a variable genuinely may not
#exist; by the second it is in the model frame and anything failing is a bug.
.err_unevaluable <- function(evaled, name.missing) {
  ee <- conditionMessage(.attr(evaled, "condition"))

  if (name.missing && startsWith(ee, "object '") && endsWith(ee, "' not found")) {
    v <- sub("object '([^']+)' not found", "\\1", ee)
    arg::err("the variable {.val {v}} cannot be found. Be sure it is entered correctly or supply a dataset that contains this variable to {.arg data}")
  }

  arg::err("{ee}")
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
    
    cbind(ttfactor[, seq_len(after), drop = FALSE], 
          addcol, 
          ttfactor[, -seq_len(after), drop = FALSE])
    
  }
  
  #Check if data exists
  data.specified <- FALSE
  if (is_not_null(data)) {
    if (is.matrix(data)) {
      data <- as.data.frame.matrix(data)
    }
    
    if (!is.data.frame(data)) {
      arg::err("the argument supplied to {.arg data} must be a data frame")
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
    arg::arg_formula(f)
  }
  
  env <- rlang::f_env(f)
  if (!data.specified) data <- env
  
  # rlang::f_lhs(f) <- NULL
  
  tt <- tryCatch(terms(f, data = data),
                 error = function(e) {
                   if (conditionMessage(e) == "'.' in formula and no 'data' argument") {
                     arg::err("'.' is not allowed in formulas")
                   }
                   arg::err("{conditionMessage(e)}")
                 })
  
  #Process RHS 
  tt.covs <- delete.response(tt)
  attr(tt.covs, "intercept") <- 0
  
  ttfactors <- .attr(tt.covs, "factors")
  ttvars <- setNames(vapply(.attr(tt.covs, "variables"), deparse1, character(1L))[-1L], rownames(ttfactors))
  
  rhs.df.type <- setNames(vapply(ttvars, function(v) {
    if (length(dim(try(eval(str2expression(add_quotes(v, "`")), data, env), silent = TRUE))) == 2L) "lit"
    else if (length(dim(try(eval(str2expression(v), data, env), silent = TRUE))) == 2L) "exp"
    else "not.a.df"
  }, character(1L)), ttvars)
  
  rhs.df <- setNames(rhs.df.type != "not.a.df", ttvars)
  
  if (any(rhs.df)) {
    term_is_interaction <- colSums(ttfactors != 0) > 1L
    
    if (any_apply(which(rhs.df), function(x) any(ttfactors[x, ] != 0 & term_is_interaction))) {
      arg::err("interactions with data frames are not allowed in the input formula")
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
        
        if (ncol(df) == 1L && is_null(colnames(df))) {
          colnames(df) <- v
        }
        else if (can_str2num(colnames(df))) {
          colnames(df) <- paste(v, colnames(df), sep = "_")
        }
        
        return(df)
      }
      
      if (ncol(df) == 1L && is_null(colnames(df))) {
        colnames(df) <- v
      }
      else if (can_str2num(colnames(df))) {
        colnames(df) <- paste(v, colnames(df), sep = "_")
      }
      
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
    
    ttfactors <- .attr(tt.covs, "factors")
    ttvars <- vapply(.attr(tt.covs, "variables"), deparse1, character(1L))[-1L] |>
      setNames(rownames(ttfactors))
  }
  
  #Check to make sure variables are valid
  #Check to make sure variables are valid
  vars <- .backtick_unevaluable(rownames(ttfactors), data, env, reject.functions = TRUE)
  
  if (!identical(vars, rownames(ttfactors))) {
    rownames(ttfactors) <- vars
    
    tt.covs <- rebuild_f(ttfactors) |>
      terms(data = data)
    
    ttfactors <- .attr(tt.covs, "factors")
    ttvars <- vapply(.attr(tt.covs, "variables"), deparse1, character(1L))[-1L]
  }
  
  tmpcovs <- try_arg(model.frame2(tt.covs, data))
  
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
    
    ttfactors <- .attr(tt.covs, "factors")
    
    ttvars <- vapply(.attr(tt.covs, "variables"), deparse1, character(1L))[-1L]
    
    na_vars <- paste0(vars_with_NA, ":<NA>")
    
    tmpcovs <- try_arg(model.frame2(tt.covs, tmpcovs))
    
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
  
  #The missingness indicators added above are new variables, so re-check.
  vars <- .backtick_unevaluable(rownames(ttfactors), tmpcovs)
  
  if (!identical(vars, rownames(ttfactors))) {
    rownames(ttfactors) <- vars

    tt.covs <- rebuild_f(ttfactors) |>
      terms(data = data)
    
    ttfactors <- .attr(tt.covs, "factors")
    ttvars <- vapply(.attr(tt.covs, "variables"), deparse1, character(1L))[-1L]
  }
  
  tmpcovs <- model.frame2(tt.covs, data = tmpcovs, drop.unused.levels = TRUE)
  
  #Check for infinite values
  covs.with.inf <- vapply(tmpcovs, function(x) is.numeric(x) && any(!is.na(x) & !is.finite(x)), logical(1L))
  if (any(covs.with.inf)) {
    s <- if (sum(covs.with.inf) == 1L) c("", "s") else c("s", "")
    arg::err("the variable{?s} {.var {names(tmpcovs)[covs.with.inf]}} contain{?s/} non-finite values, which are not allowed")
  }
  
  attr(tt.covs, "intercept") <- 1 #Add intercept to correctly process single-level factors
  mm <- model.matrix(tt.covs, data = tmpcovs, 
                     contrasts.arg = lapply(Filter(is.factor, tmpcovs),
                                            function(x) contrasts(x, contrasts = all_the_same(x))))
  
  rownames(ttfactors) <- trim_string(rownames(ttfactors), "`")
  
  mmassign <- .attr(mm, "assign")[-1L]
  mmassign2 <- setNames(factor(mmassign, levels = sort(unique(mmassign), na.last = TRUE),
                               labels = colnames(ttfactors)),
                        colnames(mm)[-1L])
  
  vars_in_each_term <- setNames(lapply(colnames(ttfactors), function(x) {
    rownames(ttfactors)[ttfactors[, x] != 0]
  }), colnames(ttfactors))
  
  all_factor_levels <- lapply(vars_in_each_term, function(v) {
    do.call("expand.grid", c(clear_null(setNames(lapply(v, function(fa) colnames(.attr(mm, "contrasts")[[fa]])), v)),
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
#A covariate matrix carries the parsed form of its column names in a `co.names`
#attribute. `[` drops attributes, which is why `.get_C2()` used to carry the matrices
#and their names in two parallel lists and subset both at every step -- eight places
#where the two could fall out of step. These keep them together instead.

#The columns `keep` selects, with their names.
.C_keep <- function(C, keep) {
  co.names <- .attr(C, "co.names")[keep]
  C <- C[, keep, drop = FALSE]

  .C_rename(C, co.names)
}

#Attaches `co.names` and derives the column names from it.
.C_rename <- function(C, co.names) {
  if (is_not_null(co.names)) {
    names(co.names) <- vapply(co.names, function(x) paste(x[["component"]], collapse = ""),
                              character(1L))

    if (NCOL(C) > 0L) {
      colnames(C) <- names(co.names)
    }
  }

  attr(C, "co.names") <- co.names

  C
}

#Merges `addl` into `covs`, dropping anything already there under the same name or
#perfectly collinear with a covariate.
.C_add_addl <- function(C, addl, drop) {
  same.name <- names(.attr(addl, "co.names")) %in% names(.attr(C, "co.names"))
  addl <- .C_keep(addl, !same.name)

  if (drop && getOption("cobalt_remove_perfect_col", max(ncol(addl), ncol(C)) <= 900)) {
    redundant <- find_perfect_col(addl, C)

    if (is_not_null(redundant)) {
      addl <- .C_keep(addl, -redundant)
    }
  }

  co.cbind(C, addl)
}

#For a 0/1 or FALSE/TRUE variable that was split into two dummies, keeps the `1` column
#alone and renames it after the variable. Also drops both when they are perfectly
#redundant with `cluster`.
.C_drop_0_1 <- function(C, cluster) {
  co.names <- .attr(C, "co.names")
  drop_0_1 <- rep.int(NA, length(co.names))

  base_of <- function(x) x[["component"]][x[["type"]] == "base"][1L]

  for (i in seq_along(co.names)) {
    if (!is.na(drop_0_1[i])) {
      next
    }

    #Only the dummies of a split factor are candidates, and never an interaction.
    if ("isep" %in% co.names[[i]][["type"]] || "fsep" %nin% co.names[[i]][["type"]]) {
      drop_0_1[i] <- FALSE
      next
    }

    #The other dummies split from the same variable.
    buddy.i <- which(vapply(co.names, function(j) {
      "isep" %nin% j[["type"]] && "fsep" %in% j[["type"]] &&
        base_of(j) == base_of(co.names[[i]])
    }, logical(1L)))

    buddies <- co.names[buddy.i]

    if (is_not_null(cluster)) {
      unsplit_var <- unsplitfactor(as.data.frame(C[, buddy.i, drop = FALSE]),
                                  base_of(buddies[[1L]]),
                                  sep = .attr(co.names, "seps")["factor"])[[1L]]

      tab <- table(cluster, unsplit_var)
      tab <- tab[rowSums(tab == 0) < ncol(tab), , drop = FALSE]

      #Each cluster takes a single level, so the variable adds nothing to it.
      if (all(rowSums(tab > 0) == 1)) {
        drop_0_1[buddy.i] <- TRUE
        next
      }
    }

    #More than two levels is not a 0/1 variable.
    if (length(buddies) > 2L) {
      drop_0_1[buddy.i] <- FALSE
      next
    }

    is_0 <- vapply(buddies, function(x) {
      x[["component"]][x[["type"]] == "level"] %in% c("0", "FALSE")
    }, logical(1L))

    is_1 <- vapply(buddies, function(x) {
      x[["component"]][x[["type"]] == "level"] %in% c("1", "TRUE")
    }, logical(1L))

    #Two levels that are not 0/1: keep the first, drop the second.
    if (!all(is_0 | is_1)) {
      drop_0_1[buddy.i] <- c(TRUE, FALSE)
      next
    }

    drop_0_1[buddy.i[is_0]] <- TRUE
    drop_0_1[buddy.i[is_1]] <- FALSE

    #The surviving dummy is named after the variable rather than the level.
    kept <- buddy.i[is_1]
    co.names[[kept]][["component"]] <- base_of(co.names[[kept]])
    co.names[[kept]][["type"]] <- "base"
  }

  .C_rename(C, co.names) |>
    .C_keep(!drop_0_1)
}

.get_C2 <- function(covs = NULL, int = FALSE, poly = 1, addl = NULL, distance = NULL,
                    treat = NULL, cluster = NULL, drop = TRUE, factor_sep = "_",
                    int_sep = " * ", ...) {
  #Gets the C matrix, holding every variable balance is assessed on. Used in
  #`balance_table()`. Each piece carries its own `co.names`, so `.C_keep()` is the only
  #thing that removes a column and its name cannot be left behind.
  if (inherits(covs, "processed_C")) {
    return(covs)
  }

  if (is_null(covs)) {
    drop <- FALSE
  }

  arg::arg_string(factor_sep)
  arg::arg_string(int_sep)

  #Process int and poly
  arg::arg_whole_number(poly)
  arg::arg_gte(poly, 1)
  poly <- round(poly)

  arg::arg_or(int,
              arg::arg_flag,
              arg::arg_and(
                arg::arg_whole_number,
                arg::arg_gt(1)
              ))

  if (is.numeric(int)) {
    if (int > poly) {
      poly <- int
    }

    int <- TRUE
  }

  center <- ...get("center", getOption("cobalt_center", default = FALSE))
  arg::arg_flag(center)
  orth <- ...get("orth", getOption("cobalt_orth", default = FALSE))
  arg::arg_flag(orth)

  seps <- .attr(.attr(covs, "co.names"), "seps")

  C <- covs

  if (is_not_null(addl)) {
    C <- .C_add_addl(C, addl, drop)
  }

  #Drop anything that is just the treatment under another name.
  if (drop && is_not_null(treat) && get.treat.type(treat) != "continuous") {
    collinear <- vapply(seq_col(C), function(i) {
      !anyNA(C[, i]) && equivalent.factors2(C[, i], treat)
    }, logical(1L))

    C <- .C_keep(C, !collinear)
  }

  int.poly <- NULL

  if (int || poly > 1) {
    #Interactions and polynomials are not built out of missingness indicators or out of
    #existing interactions.
    exclude <- vapply(.attr(C, "co.names"), function(x) {
      any(c("na", "isep") %in% x[["type"]])
    }, logical(1L))

    #`sep` must be unnamed: it ends up inside each term's `component` vector, and a
    #stray "int" name there shows up in the object `bal.tab()` returns.
    int.poly <- .int_poly_f2(C, ex = exclude, int = int, poly = poly, center = center,
                             orth = orth, sep = unname(seps["int"]),
                             co.names = .attr(C, "co.names"))

    #NULL when there is nothing to add.
    if (is_not_null(int.poly)) {
      int.poly <- .C_rename(int.poly, .attr(int.poly, "co.names"))
    }
  }

  if (drop) {
    C <- .C_drop_0_1(C, cluster)
  }
  else {
    C <- .C_rename(C, .attr(C, "co.names"))
  }

  if (is_not_null(distance)) {
    if (anyNA(distance, recursive = TRUE)) {
      arg::err("missing values are not allowed in the distance measure")
    }

    taken <- c(names(.attr(C, "co.names")), names(.attr(int.poly, "co.names")))
    distance.co.names <- .attr(distance, "co.names")

    distance <- .C_keep(distance, names(distance.co.names) %nin% taken)

    #A distance may be supplied twice under the same name.
    distance <- .C_keep(distance, unique(names(.attr(distance, "co.names"))))
  }

  #Drop covariates that an interaction or polynomial term already carries by name.
  if (drop && is_not_null(int.poly)) {
    C <- .C_keep(C, names(.attr(C, "co.names")) %nin% .attr(int.poly, "co.names"))
  }

  C <- co.cbind(distance, C, int.poly)

  if (is_null(C)) {
    return(matrix(0, nrow = length(treat), ncol = 0L,
                  dimnames = list(rownames(covs), NULL)) |>
             set_class("processed_C", .replace = FALSE))
  }

  #The separators are only fixed now, so that everything above compares names built
  #with the separators the pieces arrived with.
  co.names <- .attr(C, "co.names")

  for (i in seq_along(co.names)) {
    co.names[[i]][["component"]][co.names[[i]][["type"]] == "fsep"] <- factor_sep
    co.names[[i]][["component"]][co.names[[i]][["type"]] == "isep"] <- int_sep
  }

  seps["factor"] <- factor_sep
  seps["int"] <- int_sep
  attr(co.names, "seps") <- seps

  C <- .C_rename(C, co.names)

  attr(C, "missing.ind") <- colnames(C)[vapply(co.names, function(x) {
    "na" %in% x[["type"]]
  }, logical(1L))]

  if (is_not_null(distance)) {
    attr(C, "distance.names") <- names(.attr(distance, "co.names"))
  }

  attr(C, "var_types") <- .get_types(C)

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
          if (orth) poly(x, degree = poly)[, i] else x^i
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
    
    #Don't make ints out of multiple members of the same categorical variable
    ints_to_make[vapply(ints_to_make, function(x) {
      "fsep" %in% co.names[[x[1L]]][["type"]] && 
        "fsep" %in% co.names[[x[2L]]][["type"]] &&
        identical(co.names[[x[1L]]][["component"]][co.names[[x[1L]]][["type"]] == "base"],
                  co.names[[x[2L]]][["component"]][co.names[[x[2L]]][["type"]] == "base"])
    }, logical(1L))] <- NULL
    
    int_terms[[1L]] <- do.call("cbind", lapply(ints_to_make, function(i) d[, i[1L]] * d[, i[2L]]))
    
    int_co.names[[1L]] <- {
      if (cn)
        lapply(ints_to_make, function(x) list(component = c(co.names[[x[1L]]][["component"]], sep, co.names[[x[2L]]][["component"]]),
                                              type = c(co.names[[x[1L]]][["type"]], "isep", co.names[[x[2L]]][["type"]])))
      else
        vapply(ints_to_make, paste, character(1L), collapse = sep)
    }
  }
  else {
    int_terms <- int_co.names <- list()
  }
  
  out <- do.call("cbind", c(poly_terms, int_terms))
  out_co.names <- c(do.call("c", poly_co.names), do.call("c", int_co.names))

  #No terms to add; e.g., `poly = 1` and `int = FALSE`, `int = TRUE` with only one
  #covariate, or `poly > 1` when every covariate is binary or excluded.
  if (is_null(out)) {
    return(NULL)
  }

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
co.cbind <- function(...) {
  args <- clear_null(list(...))
  
  if (length(args) == 0L) {
    return(NULL)
  }
  
  if (length(args) == 1L) {
    return(args[[1L]])
  }
  
  co.names.list <- lapply(args, attr, "co.names")
  
  seps <- .attr(co.names.list[[1L]], "seps")
  
  out <- do.call("cbind", args)
  
  attr(out, "co.names") <- do.call("c", co.names.list)
  attr(attr(out, "co.names"), "seps") <- seps
  colnames(out) <- names(attr(out, "co.names")) <- vapply(.attr(out, "co.names"),
                                                          function(x) paste(x[["component"]], collapse = ""),
                                                          character(1L))
  
  out
}

.get_types <- function(C) {
  distance.names <- .attr(C, "distance.names")

  vapply(seq_col(C), function(i) {
    if (colnames(C)[i] %in% distance.names) "Distance"
    else if (all_the_same(C[, i]) || is_binary(C[, i])) "Binary"
    else "Contin."
  }, character(1L)) |>
    setNames(colnames(C))
}
find_perfect_col <- function(C1, C2 = NULL, fun = stats::cor) {
  
  #Finds indices of redundant vars in C1.
  C1.no.miss <- C1[, colnames(C1) %nin% .attr(C1, "missing.ind"), drop = FALSE]
  if (is_null(C2)) {
    .use <- if (anyNA(C1)) "pairwise.complete.obs" else "everything"
    C.cor <- suppressWarnings(fun(C1.no.miss, use = .use))
    s <- !lower.tri(C.cor, diag = TRUE) & !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
  }
  else {
    C2.no.miss <- C2[, colnames(C2) %nin% .attr(C2, "missing.ind"), drop = FALSE]
    .use <- if (anyNA(C1) || anyNA(C2)) "pairwise.complete.obs" else "everything"
    C.cor <- suppressWarnings(fun(C2.no.miss, y = C1.no.miss, use = .use))
    s <- !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
  }
  
  which(colSums(s) > 0)
}

model.frame2 <- function(formula, data = NULL, na.action = "na.pass", ...) {
  data <- rlang::try_fetch(force(data),
                           error = function(e) arg::err("{conditionMessage(e)}"),
                           warning = function(w) arg::wrn("{conditionMessage(w)}"))
  
  rlang::try_fetch({
    stats::model.frame(formula, data = data, na.action = na.action, ...)
  },
  error = function(e) {
    ee <- conditionMessage(e)
    if (startsWith(ee, "object '") && endsWith(ee, "' not found")) {
      v <- sub("object '([^']+)' not found", "\\1", ee)
      arg::err("the variable {.val {v}} cannot be found. Be sure it is entered correctly or supply a dataset that contains this variable to {.arg data}")
    }
    
    arg::err("{ee}")
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
        prob.w.t.mat <- w.t.mat[problems, ]
        if (NCOL(weights.df) == 1L) {
          arg::err('all weights are zero when the treatment is {.or {.val {prob.w.t.mat[, "treat_vals"]}}}')
        }
        
        errors <- vapply(unique(prob.w.t.mat[, "weight_names"]), function(i) {
          cli::format_inline('all {.var {i}} weights are zero when the treatment is {.or {.val {prob.w.t.mat[prob.w.t.mat[, "weight_names"] == i, "treat_vals"]}}}')
        }, character(1L))
        
        arg::err("{errors}")
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
        arg::err("all weights are zero")
      }
      
      errors <- vapply(prob.wts, function(i) {
        cli::format_inline('all {.var {i}} weights are zero')
      }, character(1L))
      
      arg::err("{errors}")
    }
  }
}
#How many covariates met a threshold and how many did not. Takes the threshold rather
#than recovering it by regex from the labels the package generated a moment earlier.
.baltal <- function(labels, threshold) {
  verdicts <- .threshold_verdicts(threshold)

  data.frame(count = c(sum(labels == verdicts[1L]), sum(labels == verdicts[2L])),
             row.names = verdicts)
}
.max_imbal <- function(balance.table, col.name, thresh.col.name, abs_stat) {
  clean <- balance.table[balance.table[["Type"]] != "Distance" &
                           is.finite(balance.table[[col.name]]), , drop = FALSE]
  
  maxed <- clean[which.max(abs_stat(clean[[col.name]])),
                 intersect(c(col.name, thresh.col.name), names(clean)), drop = FALSE]
  
  cbind(Variable = rownames(maxed), maxed)
}
#Resolve the covariates and the standardization defaults once in the wrapper, so that
#every child inherits them instead of deriving its own. `s.d.denom` is left to the
#caller: a leaf clears it when nothing is standardized, a wrapper keeps the user's
#value, and the multi-category wrapper precomputes one denominator per weight set.
.bal.tab_prepare <- function(X, A) {
  X[["covs"]] <- do.call(".get_C2", c(X, A[setdiff(names(A), names(X))]), quote = TRUE)

  std.defaults <- .get_std_defaults(X[["treat"]], A[["continuous"]], A[["binary"]])
  A[["continuous"]] <- std.defaults$continuous
  A[["binary"]] <- std.defaults$binary

  list(X = X, A = A, var_types = .attr(X[["covs"]], "var_types"))
}

#The tail the four recursing wrappers share: summarize the children's balance tables,
#tally thresholds against that summary, and combine the children's sample sizes. They
#differ only in what the summary is called, which aggregating functions apply, how
#sample sizes combine, and whether time points are listed -- so those are arguments,
#and the four heads stay as four honest functions.
.bal.tab_summarize <- function(child.list, summary.name, agg.funs, obs.fun,
                               include.times = FALSE) {
  p.ops <- .attr(child.list[[1L]], "print.options")
  balance <- child.list[[1L]][["Balance"]]

  out <- setNames(list(balance_summary(child.list, agg.funs = agg.funs,
                                      include.times = include.times)),
                  summary.name)

  #A threshold labels one aggregate, so the tally is only meaningful when a single
  #aggregating function was asked for.
  if (length(agg.funs) == 1L) {
    out <- c(out,
             threshold_summary(compute = .attr(balance, "compute"),
                               thresholds = .attr(balance, "thresholds"),
                               no.adj = !p.ops[["disp.adj"]],
                               balance.table = out[[summary.name]],
                               weight.names = p.ops[["weight.names"]],
                               agg.fun = agg.funs))
  }

  out[["Observations"]] <- obs.fun(child.list)

  out
}

#The balance tally and greatest-imbalance tables for each statistic with a threshold.
#
#These used to be written out three times -- once for an unadjusted-only table, once
#for a single weight set, once for several -- differing only in which sample's columns
#to read and whether the result is one column or one per weight set. The sample list
#says which, so the body is written once.
threshold_summary <- function(compute, thresholds, no.adj, balance.table,
                              weight.names = NULL, agg.fun = NULL) {
  samples <- if (no.adj) "Un" else weight.names

  #With a single weight set the threshold column carries no suffix, matching how
  #`.bal_tab_col_spec()` named it.
  bare.threshold <- !no.adj && length(samples) == 1L

  out <- do.call("c", lapply(compute, function(s) make_list(paste.(c("Balanced", "Max.Imbalance"), s))))

  no.distance <- balance.table[balance.table[["Type"]] != "Distance", , drop = FALSE]

  for (s in compute) {
    if (is_null(thresholds[[s]])) {
      out[[paste.("Balanced", s)]] <- NULL
      out[[paste.("Max.Imbalance", s)]] <- NULL
      next
    }

    prefix <- STATS[[s]][["bal.tab_column_prefix"]]
    Threshold <- STATS[[s]][["Threshold"]]

    thresh.col <- if (bare.threshold) rep_with(Threshold, samples) else paste.(Threshold, samples)
    stat.col <- .paste_col(if (is_null(agg.fun)) NULL else firstup(agg.fun), prefix) |>
      paste.(samples)

    tallies <- lapply(thresh.col, function(tc) {
      .baltal(balance.table[[tc]], thresholds[[s]])
    })

    imbalances <- lapply(seq_along(samples), function(i) {
      .max_imbal(no.distance, stat.col[i], thresh.col[i], STATS[[s]][["abs"]])
    })

    out[[paste.("Balanced", s)]] <- {
      if (length(samples) == 1L) tallies[[1L]]
      else setNames(do.call("cbind", tallies), samples)
    }

    out[[paste.("Max.Imbalance", s)]] <- {
      if (length(samples) == 1L) imbalances[[1L]]
      else cbind(Weights = samples,
                 do.call("rbind", lapply(imbalances, setNames,
                                         c("Variable", prefix, Threshold))),
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
  spec <- .bal_tab_col_spec(type, compute, thresholds,
                            samples = c("Un", weight.names),
                            threshold.samples = if (no.adj) "Un" else weight.names)

  B <- make_df(spec[["name"]], NCOL(C))
  rownames(B) <- colnames(C)
  
  #Set var type (binary/continuous)
  B[["Type"]] <- var_types %or% .get_types(C)
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
  binary <- arg::match_arg(binary, c("raw", "std"))
  if ("sds" %in% compute) {
    sd.computable <- if (binary == "std") rep.int(TRUE, nrow(B)) else !bin.vars
    if (type == "bin") {
      tn01 <- setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]], 0:1)
      
      if (un || !quick) {
        for (t in c("0", "1")) {
          sds <- rep.int(NA_real_, NCOL(C))
          if (any(sd.computable)) {
            sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], weights = NULL, s.weights = s.weights,
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
              sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], weights = weights[[i]], s.weights = s.weights,
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
        B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Un")]] <- STATS[[s]]$fun(
          C, treat = treat, weights = NULL, 
          std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
          s.d.denom = s.d.denom.list[[1L]] %or% s.d.denom[1L],
          abs = abs, s.weights = s.weights, bin.vars = bin.vars,
          weighted.weights = weights[[1L]], ...
        )
      }
      
      if (!no.adj && (!quick || s %in% disp)) {
        for (i in weight.names) {
          B[[paste.(STATS[[s]]$bal.tab_column_prefix, i)]] <- STATS[[s]]$fun(
            C, treat = treat, weights = weights[[i]],
            std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
            s.d.denom = s.d.denom.list[[i]] %or% s.d.denom[i],
            abs = abs, s.weights = s.weights, bin.vars = bin.vars, ...
          )
        }
      }
      
      
      if (all_apply(intersect(names(B), paste.(STATS[[s]]$bal.tab_column_prefix, c("Un", weight.names))), 
                    function(x) !any(is.finite(B[[x]])))) {
        disp <- disp[disp != s] 
        thresholds[[s]] <- NULL
      }
      
      if (is_not_null(thresholds[[s]])) {
        thr <- .threshold_cols(spec, s)

        for (k in seq_len(nrow(thr))) {
          B[[thr[["name"]][k]]] <-
            .threshold_label(B[[paste.(STATS[[s]]$bal.tab_column_prefix, thr[["sample"]][k])]],
                             B[["Type"]], thresholds[[s]], STATS[[s]]$abs)
        }
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
  
  if (is_null(s.weights)) {
    s.weights <- rep_with(1, treat)
  }
  
  if (is_null(discarded)) {
    discarded <- rep_with(FALSE, treat)
  }
  
  #The subclassification tables are transposed relative to the others -- one column
  #per subclass -- so they are built on their own.
  if (length(method) == 1L && method == "subclassification") {
    if (is_null(subclass)) {
      arg::err("{.arg subclass} must be a vector of subclasses")
    }
    
    matched <- !is.na(subclass)
    
    if (type == "bin") {
      nn <- make_df(c(levels(subclass), "Discarded", "All"), c(treat_names(treat), "Total"))
      
      nn[["All"]] <- c(vapply(treat_vals(treat), function(tn) sum(treat == tn), numeric(1L)), length(treat))
      nn[["Discarded"]] <- {
        if (any(discarded)) c(vapply(treat_vals(treat), function(tn) sum(treat[discarded] == tn), numeric(1L)), sum(discarded))
        else NULL
      }
      
      for (k in levels(subclass)) {
        tk <- treat[matched & subclass == k]
        nn[[k]] <- c(vapply(treat_vals(treat), function(tn) sum(tk == tn), numeric(1L)),
                     length(tk))
      }
      
      for (tnn in names(treat_names(treat))) {
        small.subclass <- nn[treat_names(treat)[tnn], levels(subclass)] <= 1L
        if (any(small.subclass)) {
          arg::wrn("not enough {tnn} units in {cli::qty(sum(small.subclass))} subclass{?es} {levels(subclass)[small.subclass]}")
        }
      }
    }
    else if (type == "cont") {
      nn <- make_df(c(levels(subclass), "All"), "Total")
      
      nn[, "All"] <- length(treat)
      
      for (k in levels(subclass)) {
        nn[[k]] <- sum(matched & subclass == k)
      }
      
      small.subclass <- nn[, levels(subclass)] <= 1L
      if (any(small.subclass)) {
        arg::wrn("not enough units in subclass{?es} {levels(subclass)[small.subclass]}")
      }
    }
    
    attr(nn, "tag") <- "Sample sizes by subclass"
    
    return(nn)
  }
  
  #A binary treatment reports one column per treatment group and a continuous
  #treatment a single `Total` column; the counting is otherwise the same, so it is
  #expressed once over a list of per-column unit sets.
  cols <- if (type == "bin") treat_names(treat) else "Total"
  in.col <- {
    if (type == "bin") unname(lapply(treat_vals(treat), function(tn) treat == tn))
    else list(rep_with(TRUE, treat))
  }
  
  #Fills one row of `nn` by counting within each column's units.
  by_col <- function(f) vapply(in.col, f, numeric(1L), USE.NAMES = FALSE)
  
  if (is_null(weights)) {
    nn <- make_df(cols, "All")
    nn["All", ] <- by_col(function(i) ESS(s.weights[i]))
    
    attr(nn, "ss.type") <- {
      if (nunique.gt(s.weights, 2L) || !any(s.weights == 1) || !all(s.weights %in% c(0, 1))) "ess"
      else "ss"
    }
  }
  else if (NCOL(weights) == 1L) {
    if (method == "matching") {
      nn <- make_df(cols, c("All (ESS)", "All (Unweighted)", "Matched (ESS)",
                            "Matched (Unweighted)", "Unmatched", "Discarded"))
      
      nn["All (ESS)", ] <- by_col(function(i) ESS(s.weights[i]))
      nn["All (Unweighted)", ] <- by_col(function(i) sum(i & s.weights > 0))
      nn["Matched (ESS)", ] <- by_col(function(i) ESS((weights[, 1L] * s.weights)[i & !discarded]))
      nn["Matched (Unweighted)", ] <- by_col(function(i) sum(i & weights[, 1L] > 0 & s.weights > 0 & !discarded))
      nn["Unmatched", ] <- by_col(function(i) sum(i & weights[, 1L] == 0 & !discarded))
      nn["Discarded", ] <- by_col(function(i) sum(i & discarded))
      
      attr(nn, "ss.type") <- rep.int("ss", NROW(nn))
    }
    else if (method == "weighting") {
      nn <- make_df(cols, c("Unadjusted", "Adjusted", "Discarded"))
      
      nn["Unadjusted", ] <- by_col(function(i) ESS(s.weights[i]))
      nn["Adjusted", ] <- by_col(function(i) ESS((weights[, 1L] * s.weights)[i & !discarded]))
      nn["Discarded", ] <- by_col(function(i) sum(i & discarded))
      
      attr(nn, "ss.type") <- c("ss", "ess", "ss")
    }
  }
  else {
    nn <- make_df(cols, c("All", names(weights)))
    
    nn["All", ] <- by_col(function(i) ESS(s.weights[i]))
    for (j in seq_col(weights)) {
      nn[1L + j, ] <- by_col(function(i) ESS((weights[, j] * s.weights)[i & !discarded]))
    }
    
    attr(nn, "ss.type") <- c("ss", rep_with("ess", method))
  }
  
  #A `Discarded` row that no unit reached is not shown at all.
  if ("Discarded" %in% rownames(nn) && !any(discarded)) {
    attr(nn, "ss.type") <- .attr(nn, "ss.type")[rownames(nn) != "Discarded"]
    nn <- nn[rownames(nn) != "Discarded", , drop = FALSE]
  }
  
  attr(nn, "tag") <- {
    if (length(.attr(nn, "ss.type")) > 1L && all(.attr(nn, "ss.type")[-1L] == "ess")) {
      "Effective sample sizes"
    }
    else "Sample sizes"
  }
  
  nn
}

balance_summary <- function(bal.tab.list, agg.funs, include.times = FALSE) {
  type <- .attr(bal.tab.list[[1L]], "print.options")[["type"]]
  disp <- .attr(bal.tab.list[[1L]], "print.options")[["disp"]]
  compute <- .attr(bal.tab.list[[1L]], "print.options")[["compute"]]
  thresholds <- .attr(bal.tab.list[[1L]], "print.options")[["thresholds"]]
  quick <- .attr(bal.tab.list[[1L]], "print.options")[["quick"]]
  weight.names <- .attr(bal.tab.list[[1L]], "print.options")[["weight.names"]] %or% "Adj"
  abs <- .attr(bal.tab.list[[1L]], "print.options")[["abs"]]
  #Read from `disp.adj`, as `threshold_summary()`'s callers do: a subclassified child
  #carries no `nweights`, and `NULL == 0` is `logical(0)`, not `TRUE`.
  no.adj <- !isTRUE(.attr(bal.tab.list[[1L]], "print.options")[["disp.adj"]])
  
  balance.list <- clear_null(grab(bal.tab.list, "Balance"))
  
  Brownames <- unique(unlist(lapply(balance.list, rownames), use.names = FALSE))
  
  Agg.Funs <- firstup(if (quick) agg.funs else c("min", "mean", "max"))
  Agg.Funs.Given <- firstup(agg.funs)
  
  if (length(Agg.Funs) == 1L && Agg.Funs == "Max") {
    abs <- TRUE
  }
  
  #A threshold labels one aggregate, so it is only meaningful when a single
  #aggregating function was asked for.
  spec <- .bal_tab_col_spec(type, compute, thresholds,
                            samples = c("Un", weight.names),
                            quantities = NULL,
                            agg.funs = Agg.Funs,
                            threshold.samples = {
                              if (length(Agg.Funs.Given) != 1L) character(0L)
                              else if (no.adj) "Un"
                              else weight.names
                            },
                            threshold.agg.fun = Agg.Funs.Given,
                            include.times = include.times)

  B <- make_df(spec[["name"]], Brownames)
  
  B[["Type"]] <- unlist(lapply(Brownames, function(x) na.rem(unique(vapply(balance.list, function(y) if (x %in% rownames(y)) y[[x, "Type"]] else NA_character_, character(1L))))), use.names = FALSE)
  
  if (include.times) {
    B[["Times"]] <- vapply(Brownames, function(x) toString(seq_along(balance.list)[vapply(balance.list, function(y) x %in% rownames(y), logical(1L))]), character(1L))[Brownames]
  }
  
  for (Agg.Fun in Agg.Funs) {
    for (s in compute[compute %in% all_STATS(type)]) {
      abs0 <- function(x) {
        if (is_null(x)) NA_real_
        else if (abs) STATS[[s]]$abs(x)
        else x
      }
      
      #A statistic may name its own aggregator for a given function -- variance
      #ratios are averaged geometrically, not arithmetically. Indexing by a name the
      #vector does not carry yields NA, and indexing NULL yields NULL, so both mean
      #"no override".
      own <- unname(STATS[[s]]$agg_fun[tolower(Agg.Fun)])

      if (anyNA(own)) {
        own <- NULL
      }

      agg <- function(x, ...) {
        if (!any(is.finite(x))) NA_real_
        else if (is_not_null(own)) get(own)(x)
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
  
  #Assign X.Threshold values
  for (s in compute[compute %in% all_STATS(type)]) {
    if (is_null(thresholds[[s]])) next

    thr <- .threshold_cols(spec, s)

    for (k in seq_len(nrow(thr))) {
      B[[thr[["name"]][k]]] <-
        .threshold_label(B[[paste.(Agg.Funs.Given, STATS[[s]]$bal.tab_column_prefix,
                                   thr[["sample"]][k])]],
                         B[["Type"]], thresholds[[s]], STATS[[s]]$abs)
    }
  }

  B
}


#base.bal.tab.imp
samplesize_across_imps <- function(obs.list) {
  obs.list <- clear_null(obs.list)
  
  obs <- Reduce("+", obs.list) / length(obs.list)
  attr(obs, "tag") <- sprintf("Average %s across imputations",
                              tolower(.attr(obs.list[[1L]], "tag")))
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
  attr(obs, "tag") <- .attr(bal.tab.multi.list[[1L]][["Observations"]], "tag")
  attr(obs, "ss.type") <- .attr(bal.tab.multi.list[[1L]][["Observations"]], "ss.type")
  
  obs
}

#base.bal.tab.cluster
samplesize_across_clusters <- function(obs.list) {
  obs.list <- clear_null(obs.list)
  obs <- Reduce("+", obs.list)
  attr(obs, "tag") <- sprintf("Total %s across clusters",
                              tolower(.attr(obs.list[[1L]], "tag")))
  
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
  spec <- .bal_tab_col_spec(type, compute, thresholds, samples = "Adj")

  B <- make_df(spec[["name"]], colnames(C))
  
  #Set var type (binary/continuous)
  B[["Type"]] <- var_types %or% .get_types(C)
  bin.vars <- B[["Type"]] == "Binary"
  
  SB <- setNames(rep(list(B), nlevels(subclass)), levels(subclass))
  
  binary <- arg::match_arg(binary, c("raw", "std"))
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
        sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = in.subclass,
                                       s.weights = s.weights)
        SB[[i]][["SD.Adj"]] <- sds
      }
    }
    
    for (s in all_STATS(type)) {
      if (s %in% compute && !subclass_w_empty[i]) {
        stat.col <- paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")

        SB[[i]][[stat.col]] <- STATS[[s]]$fun(C, treat = treat, weights = NULL,
                                              std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                              s.d.denom = s.d.denom,
                                              abs = FALSE, s.weights = s.weights,
                                              bin.vars = bin.vars, subset = in.subclass)

        if (all_apply(SB[[i]][stat.col], function(x) !any(is.finite(x)))) {
          disp <- disp[disp != s]
          thresholds[[s]] <- NULL
        }

        if (is_not_null(thresholds[[s]])) {
          SB[[i]][[.threshold_cols(spec, s)[["name"]]]] <-
            .threshold_label(SB[[i]][[stat.col]], SB[[i]][["Type"]],
                             thresholds[[s]], STATS[[s]]$abs)
        }
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

#The balance summary across subclasses for a continuous treatment.
#
#A binary treatment expresses subclassification as weights and hands the whole job to
#`balance_table()`. That is unavailable here: no set of unit weights makes a continuous
#treatment independent of the covariates within a subclass. So the subclass-specific
#statistics are combined instead, weighting each subclass by its share of the
#subclassified units -- the share `strata2weights()` gives it in the binary case.
#Standard deviations are combined in quadrature so that the result is itself a standard
#deviation rather than an average of them.
#
#`balance.table` is an unadjusted `balance_table()`, which supplies the layout, the
#variable types, and the unadjusted half of the values. The only columns of the layout
#it leaves out are the adjusted thresholds, since it had nothing to adjust.
balance_table_across_subclass <- function(balance.table, subclass.balance, subclass,
                                          type, s.weights = NULL, thresholds = list(),
                                          abs = FALSE) {
  #`threshold.samples` is what a `balance_table()` with something to adjust would use,
  #and has to be: `print()` decides which columns to show from the print options alone,
  #so a column the layout does not predict shifts every column after it.
  spec <- .bal_tab_col_spec(type, .attr(balance.table, "compute"), thresholds,
                            samples = c("Un", "Adj"), threshold.samples = "Adj")

  B <- make_df(spec[["name"]], rownames(balance.table))

  shared <- intersect(names(B), names(balance.table))
  B[shared] <- balance.table[shared]

  if (is_null(s.weights)) {
    s.weights <- rep_with(1, subclass)
  }

  #Sampling weights make a subclass's share its share of the population rather than of
  #the sample, which is what keeps the aggregated means equal to the unadjusted ones.
  in.subclass <- !is.na(subclass)
  p <- vapply(levels(subclass), function(k) {
    sum(s.weights[in.subclass & subclass == k])
  }, numeric(1L)) |>
    prop.table()

  across <- function(col, quadrature = FALSE) {
    x <- do.call("cbind", lapply(subclass.balance, `[[`, col))

    if (quadrature) sqrt(drop(x^2 %*% p)) else drop(x %*% p)
  }

  for (i in which(spec[["sample"]] == "Adj")) {
    nm <- spec[["name"]][i]
    s <- spec[["stat"]][i]

    #A statistic is folded the way that statistic defines, and a threshold is read off
    #the statistic column filled just above it -- the spec orders each statistic before
    #its own threshold.
    B[[nm]] <- switch(spec[["quantity"]][i],
                      "sds" = across(nm, quadrature = TRUE),
                      "stat" = if (abs) STATS[[s]]$abs(across(nm)) else across(nm),
                      "threshold" = .threshold_label(B[[paste.(STATS[[s]]$bal.tab_column_prefix, "Adj")]],
                                                     B[["Type"]], thresholds[[s]], STATS[[s]]$abs),
                      across(nm))
  }

  B
}

#Misc

check_arg_lengths <- function(...) {
  dots_names <- vapply(match.call(expand.dots = FALSE)$..., deparse1,
                       character(1L))
  
  lens <- setNames(vapply(seq_len(...length()), function(i) len(...elt(i)), numeric(1L)),
                   dots_names)
  
  supplied <- lens > 0L
  if (!all_the_same(lens[supplied])) {
    arg::err("{.arg {dots_names[supplied]}} must have the same number of units")
  }
}

#The midpoint rule, used as `col_w_ovl()`'s fallback when `integrate()` fails. The
#trapezoidal and Simpson's rules this also offered were never asked for.
intapprox <- function(f, from, to, steps) {
  seg <- seq(from, to, length.out = steps)
  mids <- (seg[-1L] + seg[-steps]) / 2

  sum(f(mids)) * (seg[2L] - seg[1L])
}
