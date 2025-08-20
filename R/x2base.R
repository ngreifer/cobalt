#Functions to convert object to base.bal.tab input

x2base <- function(x, ...) {
  UseMethod("x2base")
}

#' @exportS3Method NULL
x2base.matchit <- function(x, ...) {
  #Process matchit
  
  #Process data and get imp
  m.data <- if (NROW(x[["model"]][["data"]]) == length(x[["treat"]])) x[["model"]][["data"]] else NULL 
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, m.data)
  
  #Process treat
  treat <- process_treat(x[["treat"]], data, m.data)
  
  #Process covs
  if (is.data.frame(x[["X"]])) {
    covs <- get_covs_from_formula(data = x[["X"]])
  }
  else if (is_not_null(x[["model"]][["model"]])) {
    if (nrow(x[["model"]][["model"]]) == length(treat)) {
      covs <- get_covs_from_formula(x[["formula"]], data = x[["model"]][["model"]])
    }
    else {
      #Recreating covs from model object and x[["X"]]. Have to do this because when 
      #discard != NULL and reestimate = TRUE, cases are lost. This recovers them.
      
      # if (is_not_null(data)) {
      #     covs <- get_covs_from_formula(x[["formula"]], data = x[["model"]][["model"]])
      # }
      # else {
      order <- setNames(attr(x[["model"]][["terms"]], "order"),
                        attr(x[["model"]][["terms"]], "term.labels"))
      assign <- setNames(attr(x[["X"]], "assign"), colnames(x[["X"]]))
      assign1 <- assign[assign %in% which(order == 1)] #Just main effects
      
      dataClasses <- attr(x[["model"]][["terms"]], "dataClasses")
      factors.to.unsplit <- names(dataClasses)[dataClasses %in% c("factor", "character", "logical")]
      f0 <- setNames(lapply(factors.to.unsplit, 
                            function(z) {
                              if (dataClasses[z] == "factor")
                                list(levels = levels(x[["model"]][["model"]][[z]]),
                                     faclev = paste0(z, levels(x[["model"]][["model"]][[z]])))
                              else 
                                list(levels = unique(x[["model"]][["model"]][[z]]),
                                     faclev = paste0(z, unique(x[["model"]][["model"]][[z]])))
                            }),
                     factors.to.unsplit)
      covs <- as.data.frame(x[["X"]][, names(assign1)])
      for (i in factors.to.unsplit) {
        covs <- unsplitfactor(covs, i, sep = "",
                              dropped.level = f0[[i]][["levels"]][f0[[i]][["faclev"]] %nin% colnames(x[["X"]])])
        if (dataClasses[i] == "logical") {
          covs[[i]] <- as.logical(covs[[i]])
        }
      }
      covs <- get_covs_from_formula(x[["formula"]], data = covs)
      # }
    }
  }
  else if (is_not_null(x[["formula"]]) && is_not_null(data)) {
    covs <- get_covs_from_formula(x[["formula"]], data = data)
  }
  else {
    covs <- get_covs_from_formula(data = x[["X"]])
  }
  
  #Get estimand
  estimand <- ...get("estimand", x[["estimand"]])
  
  #Get method
  if (inherits(x, "matchit.subclass")) {
    method <- ...get("method")
    
    method <- {
      if (is_null(method) || !chk::vld_string(method)) "subclassification"
      else match_arg(method, c("weighting", "subclassification"))
    }
  }
  else {
    method <- "matching"
  }
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, m.data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, m.data, covs),
                               obj.distance = x[["distance"]], 
                               obj.distance.name = "distance")
  
  #Process focal
  focal <- x[["focal"]]
  if (get.treat.type(treat) == "binary") {
    if (is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with matchit objects")
  }
  subclass <- switch(method,
                     "subclassification" = as.factor(x[["subclass"]]),
                     NULL)
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with matchit objects")
  }
  
  #Process weights
  if (is_not_null(x[["weights"]]) && !all_the_same(x[["weights"]])) {
    weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, m.data))
    method <- attr(weights, "method")
  }
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["s.weights"]]),
                                 data, m.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, m.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  discarded <- x[["discarded"]]
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "matchit()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) stats <- unique(c(stats, names(thresholds)))
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- x[["call"]]
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names)
}

#' @exportS3Method NULL
x2base.ps <- function(x, ...) {
  #Process ps
  stop.method <- ...get("stop.method")
  if (is_null(stop.method) && ...length() > 0L && !nzchar(...names()[1L])) {
    stop.method <- ...elt(1L)
  }
  
  if (is_null(stop.method) && is_not_null(...get("full.stop.method"))) {
    stop.method <- ...get("full.stop.method")
  }
  
  available.stop.methods <- names(x[["w"]])
  
  if (is_null(stop.method)) {
    rule1 <- available.stop.methods
  }
  else if (is.character(stop.method)) {
    rule1 <- available.stop.methods[vapply(tolower(available.stop.methods),
                                           function(z) any(startsWith(z, tolower(stop.method))), logical(1L))]
    if (is_null(rule1)) {
      .wrn(sprintf("`stop.method` should be %s. Using all available stop methods instead",
                   word_list(available.stop.methods, and.or = "or", quotes = 2L)))
      rule1 <- available.stop.methods
    }
  }
  else if (is.numeric(stop.method) && any(stop.method %in% seq_along(available.stop.methods))) {
    if (!all(stop.method %in% seq_along(available.stop.methods))) {
      .wrn(sprintf("there are %s stop methods available, but you requested %s",
                   length(available.stop.methods), 
                   word_list(setdiff(stop.method, seq_along(available.stop.methods)), and.or = "and")))
    }
    rule1 <- available.stop.methods[stop.method %in% seq_along(available.stop.methods)]
  }
  else {
    .wrn("`stop.method` should be %s. Using all available stop methods instead",
         word_list(available.stop.methods, and.or = "or", quotes = 2L))
    rule1 <- available.stop.methods
  }
  
  s <- available.stop.methods[match(tolower(rule1), tolower(available.stop.methods))]
  
  #Process data and get imp
  ps.data <- x[["data"]]
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, ps.data)
  
  #Process treat
  treat <- process_treat(x[["treat"]], data, ps.data)
  
  #Process covs
  if (is_not_null(x[["gbm.obj"]][["var.names"]])) {
    covs <- reformulate(x[["gbm.obj"]][["var.names"]]) |>
      get_covs_from_formula(data = ps.data)
  }
  else if (is_not_null(...get("formula")) || is_not_null(...get("covs"))) {
    t.c <- .use_tc_fd(formula = ...get("formula"),
                      data = if_null_then(data, ps.data),
                      covs = ...get("covs"),
                      needs.treat = FALSE)
    
    #Process covs
    covs <- t.c[["covs"]]
  }
  else if (all(x[["gbm.obj"]][["feature_names"]] %in% names(ps.data))) {
    covs <- reformulate(x[["gbm.obj"]][["feature_names"]]) |>
      get_covs_from_formula(data = ps.data)
  }
  else {
    .err('when `version = "xgboost"` in the call to `ps()` and any variables are categorical, `formula` or `covs` must be supplied')
  }
  
  #Get estimand
  estimand <- x[["estimand"]]
  
  #Get method
  method <- rep_with("weighting", s)
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, ps.data))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, ps.data, covs),
                               obj.distance = x[["ps"]][s], 
                               obj.distance.name = {
                                 if (length(s) > 1L) paste.("prop.score", substr(s, 1L, nchar(s) - 4L))
                                 else "prop.score"
                               })
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with ps objects")
  }
  
  if (get.treat.type(treat) == "binary") {
    if (is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with ps objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with ps objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, ps.data), 
                             stop.method = s, estimand = estimand)
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["sampw"]]),
                                 data, ps.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, ps.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "ps()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) stats <- unique(c(stats, names(thresholds)))
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  # call <- ps[["parameters"]]
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("binary")
}

#' @exportS3Method NULL
x2base.mnps <- function(x, ...) {
  #Process mnps
  stop.method <- ...get("stop.method")
  if (is_null(stop.method) && ...length() > 0L && !nzchar(...names()[1L])) {
    stop.method <- ...elt(1L)
  }
  
  if (is_null(stop.method) && is_not_null(...get("full.stop.method"))) {
    stop.method <- ...get("full.stop.method")
  }
  
  available.stop.methods <- x[["stopMethods"]]
  
  if (is_null(stop.method)) {
    rule1 <- available.stop.methods
  }
  else if (is.character(stop.method)) {
    rule1 <- available.stop.methods[vapply(tolower(available.stop.methods),
                                           function(z) any(startsWith(z, tolower(stop.method))), logical(1L))]
    if (is_null(rule1)) {
      .wrn(sprintf("`stop.method` should be %s. Using all available stop methods instead",
                   word_list(available.stop.methods, and.or = "or", quotes = 2L)))
      rule1 <- available.stop.methods
    }
  }
  else if (is.numeric(stop.method) && any(stop.method %in% seq_along(available.stop.methods))) {
    if (!all(stop.method %in% seq_along(available.stop.methods))) {
      .wrn(sprintf("there are %s stop methods available, but you requested %s",
                   length(available.stop.methods), 
                   word_list(setdiff(stop.method, seq_along(available.stop.methods)), and.or = "and")))
    }
    rule1 <- available.stop.methods[stop.method %in% seq_along(available.stop.methods)]
  }
  else {
    .wrn("`stop.method` should be %s. Using all available stop methods instead",
         word_list(available.stop.methods, and.or = "or", quotes = 2))
    rule1 <- available.stop.methods
  }
  
  s <- available.stop.methods[match(tolower(rule1), tolower(available.stop.methods))]
  
  #Process data and get imp
  mnps.data <- x[["data"]]
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, mnps.data)
  
  #Process treat
  treat <- process_treat(x[["treatVar"]], data, mnps.data)
  
  #Process covs
  .v <- if_null_then(x[["balanceVars"]],
                     x[["psList"]][[1L]][["gbm.obj"]][["var.names"]])
  
  covs <- get_covs_from_formula(reformulate(.v), mnps.data)
  
  #Get estimand
  estimand <- x[["estimand"]]
  
  #Get method
  method <- rep_with("weighting", s)
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, mnps.data))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, mnps.data))
  
  #Process focal
  focal <- x[["treatATT"]]
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with mnps objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with mnps objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, mnps.data), 
                             stop.method = s)
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["sampw"]]),
                                 data, mnps.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, mnps.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "mnps()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) stats <- unique(c(stats, names(thresholds)))
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- NULL
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("multi")
}

#' @exportS3Method NULL
x2base.ps.cont <- function(x, ...) {
  #Process data and get imp
  ps.data <- x[["data"]]
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, ps.data)
  
  #Process treat
  treat <- process_treat(x[["treat"]], data, ps.data)
  
  #Process covs
  covs <- reformulate(x[["gbm.obj"]][["var.names"]]) |>
    get_covs_from_formula(ps.data)
  
  #Get estimand
  
  #Get method
  method <- "weighting"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, ps.data))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, ps.data))
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with ps.cont objects")
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with ps.cont objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with ps.cont objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, ps.data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["sampw"]]),
                                 data, ps.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, ps.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "ps.cont()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) stats <- unique(c(stats, names(thresholds)))
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  # call <- ps.cont[["parameters"]]
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("cont")
}

#' @exportS3Method NULL
x2base.Match <- function(x, ...) {
  #Process Match
  if (is_not_null(x) && !is.list(x)) {
    .err("the supplied Match object contains no valid matches")
  }
  
  #Process data and get imp
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data)
  
  #Process treat
  t.c <- .use_tc_fd(...get("formula"), data, ...get("treat"), ...get("covs"))
  treat <- process_treat(t.c[["treat"]], data)
  
  #Process covs
  covs <- t.c[["covs"]]
  
  #Get estimand
  estimand <- x[["estimand"]]
  
  #Get method
  method <- "matching"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, covs))
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with Match objects")
  }
  
  if (get.treat.type(treat) == "binary") {
    if (is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with Match objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with Match objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  discarded <- rep_with(FALSE, treat)
  if (is_not_null(x[["index.dropped"]])) {
    discarded[x[["index.dropped"]]] <- TRUE
  }
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "Match()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) stats <- unique(c(stats, names(thresholds)))
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- NULL
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("binary")
}

#' @exportS3Method NULL
x2base.formula <- function(x, ...) {
  x2base.data.frame(x = x, ...)
}

#' @exportS3Method NULL
x2base.data.frame <- function(x, ...) {
  #Process data.frame
  
  #Process data and get imp
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data)
  
  #Process treat
  treat <- {
    if (rlang::is_formula(x))
      get_treat_from_formula(x, data, treat = ...get("treat"))
    else
      ...get("treat")
  }
  treat <- process_treat(treat, data)
  
  #Process covs
  covs <- x
  if (is.null(covs)) {
    .err("`covs` data.frame must be specified")
  }
  
  if (rlang::is_formula(covs)) {
    covs <- get_covs_from_formula(covs, data = data)
    # if (is_null(covs)) {
    #     .err("The right hand side of the formula must contain covariates for which balance is to be assessed")
    # }
  }
  
  if (is_null(attr(covs, "co.names"))) {
    if (is.matrix(covs)) covs <- as.data.frame.matrix(covs)
    covs <- get_covs_from_formula(data = covs)
  }
  
  #Get estimand
  estimand <- ...get("estimand")
  
  #Get method
  specified <- setNames(rep.int(FALSE, 3L), c("match.strata", "subclass", "weights"))
  
  for (i in names(specified)) {
    specified[i] <- is_not_null(...get(i))
  }
  
  .specified_method <- character()
  .specified_args <- character()
  .assuming <- character()
  .ignoring <- character()
  .not_present <- character()
  .using <- character()
  
  method <- ...get("method")
  if (is_null(method)) {
    if (any(specified)) {
      .using <- names(specified)[specified][1L]
      method <- switch(.using,
                       match.strata = "matching",
                       subclass = "subclassification",
                       weights = "weighting")
      
      if (sum(specified) > 1L) {
        .specified_args <- names(specified)[specified]
        .assuming <- method
        .ignoring <- setdiff(names(specified)[specified], .using)
      }
    }
    else {
      method <- "weighting"
    }
  }
  else {
    .specified_method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    
    if (length(method) == 1L) {
      if (.specified_method == "weighting") {
        if (specified["weights"]) {
          method <- "weighting"
          .using <- "weights"
          
          if (sum(specified) > 1L) {
            .specified_args <- names(specified)[specified]
            .ignoring <- setdiff(names(specified)[specified], .using)
          }
        }
        else if (any(specified)) {
          .using <- names(specified)[specified][1L]
          method <- switch(.using,
                           match.strata = "matching",
                           subclass = "subclassification",
                           weights = "weighting")
          .assuming <- method
          .not_present <- "weights"
        }
        else {
          method <- "matching"
        }
      }
      else if (.specified_method == "matching") {
        if (specified["match.strata"]) {
          method <- "matching"
          .using <- "match.strata"
          
          if (sum(specified) > 1L) {
            .specified_args <- names(specified)[specified]
            .ignoring <- setdiff(names(specified)[specified], .using)
          }
        }
        else if (specified["weights"]) {
          method <- "matching"
          .using <- "weights"
          
          if (sum(specified) > 1L) {
            .specified_args <- names(specified)[specified]
            .ignoring <- setdiff(names(specified)[specified], .using)
          }
        }
        else if (specified["subclass"]) {
          method <- "subclassification"
          .using <- "subclass"
          .not_present <- c("weights", "match.strata")
          .assuming <- method
        }
        else {
          method <- "matching"
        }
      }
      else if (.specified_method == "subclassification") {
        if (specified["subclass"]) {
          method <- "subclassification"
          .using <- "subclass"
          
          if (sum(specified) > 1L) {
            .specified_args <- names(specified)[specified]
            .ignoring <- setdiff(names(specified)[specified], .using)
          }
        }
        else if (any(specified)) {
          .using <- names(specified)[specified][1L]
          method <- switch(.using,
                           match.strata = "matching",
                           subclass = "subclassification",
                           weights = "weighting")
          .assuming <- method
          .not_present <- "subclass"
        }
        else {
          method <- "matching"
        }
      }
    }
    else {
      if (specified["subclass"] || any(.specified_method == "subclassification")) {
        .err("subclassification cannot be specified along with other methods")
      }
      
      if (specified["match.strata"]) {
        .err("only weights can be specified with multiple methods")
      }
      
      if (specified["weights"]) {
        method <- .specified_method
        .using <- "weights"
      }
      else {
        .wrn("multiple methods were specified, but no weights were provided. Providing unadjusted data only")
        method <- "matching"
      }
    }
  }
  
  .m <- NULL
  if (is_not_null(.using)) {
    if (is_not_null(.specified_args) && is_not_null(.ignoring)) {
      if (is_not_null(.assuming)) {
        .m <- sprintf("%s specified. Assuming %s and using %s and ignoring %s",
                      word_list(.specified_args, and.or = "and", is.are = TRUE, quotes = "`"),
                      add_quotes(.assuming),
                      word_list(.using, and.or = "and", quotes = "`"),
                      word_list(.ignoring, and.or = "and", quotes = "`"))
      }
      else {
        .m <- sprintf("%s specified. Using %s and ignoring %s",
                      word_list(.specified_args, and.or = "and", is.are = TRUE, quotes = "`"),
                      word_list(.using, and.or = "and", quotes = "`"),
                      word_list(.ignoring, and.or = "and", quotes = "`"))
      }
    }
    else if (is_not_null(.specified_method) && is_not_null(.not_present) && is_not_null(.assuming)) {
      .m <- sprintf("`method = %s` is specified, but %s was not supplied. Assuming %s and using %s instead",
                    add_quotes(.specified_method, 2L),
                    add_quotes(.not_present, "`"),
                    add_quotes(.assuming),
                    word_list(.using, and.or = "and", quotes = "`"))
    }
  }
  
  if (is_not_null(.m)) {
    .msg(.m)
  }
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, covs))
  
  #Process focal
  focal <- process_focal(...get("focal"), treat)
  
  if (get.treat.type(treat) == "binary") {
    if (is_null(focal) && is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process subclass
  if ("subclass" %in% .using) {
    subclass <- .process_vector(...get("subclass"), 
                                datalist = list(data),
                                name = "subclass", 
                                which = "subclass membership",
                                missing.okay = TRUE) |>
      factor()
    weights <- NULL
  }
  
  #Process match.strata
  if ("match.strata" %in% .using) {
    match.strata <- .process_vector(...get("match.strata"), 
                                    datalist = list(data),
                                    name = "match.strata", 
                                    which = "matching strata membership",
                                    missing.okay = TRUE)
    
    weights <- data.frame(weights = strata2weights(match.strata,
                                                   treat = treat,
                                                   estimand = estimand,
                                                   focal = focal))
  }
  
  #Process weights
  if ("weights" %in% .using) {
    weights <- process_weights(NULL, list(...), treat, covs, method, addl.data = list(data))
    method <- attr(weights, "method")
  }
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  discarded <- ...get("discarded")
  
  #Process length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp)
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names)
}

#' @exportS3Method NULL
x2base.CBPS <- function(x, ...) {
  #Process CBPS
  
  #Process data and get imp
  c.data <- x[["data"]]
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, c.data)
  
  #Process treat
  treat <- get_treat_from_formula(x[["formula"]], c.data) |>
    process_treat(data, c.data)
  
  #Process covs
  covs <- get_covs_from_formula(x[["formula"]], c.data)
  
  #Get estimand
  estimand <- ...get("estimand")
  if (is_not_null(estimand)) {
    .err("`estimand` is not allowed with CBPS objects")
  }
  
  #Get method
  method <- "weighting"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, c.data))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, c.data),
                               obj.distance = if (get.treat.type(treat) == "binary") x[["fitted.values"]], 
                               obj.distance.name = "prop.score")
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with CBPS objects")
  }
  
  if (get.treat.type(treat) == "binary") {
    if (is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with CBPS objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with CBPS objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, c.data), 
                             use.weights = ...get("use.weights"))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"),
                                 data, c.data)
  if (is_not_null(s.weights)) {
    weights <- weights / s.weights #Because CBPS weights contain s.weights in them
  }
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, c.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "CBPS()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- x[["call"]]
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names)
}

#' @exportS3Method NULL
x2base.ebalance <- function(x, ...) {
  #Process ebalance
  
  #Process data and get imp
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data)
  
  #Process treat
  t.c <- .use_tc_fd(...get("formula"), data, ...get("treat"), ...get("covs"))
  treat <- process_treat(t.c[["treat"]], data)
  
  #Process covs
  covs <- t.c[["covs"]]
  
  #Get estimand
  estimand <- "ATT"
  
  #Get method
  method <- "weighting"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, covs))
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with ebalance objects")
  }
  
  focal <- switch(toupper(estimand), 
                  "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                  "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                  NULL)
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with ebalance objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with ebalance objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "ebalance()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) stats <- unique(c(stats, names(thresholds)))
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("binary")
}

#' @exportS3Method NULL
x2base.optmatch <- function(x, ...) {
  #Process optmatch
  if (all(is.na(x))) {
    .err("the supplied optmatch abject contains no valid matches")
  }
  
  #Process data and get imp
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data)
  
  #Process treat
  t.c <- .use_tc_fd(...get("formula"), data = data, covs = ...get("covs"),
                    treat = ...get("treat", as.numeric(attr(x, "contrast.group"))))
  treat <- process_treat(t.c[["treat"]], data)
  
  #Process covs
  covs <- t.c[["covs"]]
  
  #Get estimand
  estimand <- ...get("estimand", "ATT")
  
  #Get method
  method <- "matching"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, covs))
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with optmatch objects")
  }
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with optmatch objects")
  }
  
  if (get.treat.type(treat) == "binary") {
    if (is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with optmatch objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = sprintf("%s()", deparse1(attr(x, "call")[[1L]])))
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- attr(x, "call")
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("binary")
}

#' @exportS3Method NULL
x2base.cem.match <- function(x, ...) {
  #Process cem.match
  if (inherits(x, "cem.match.list")) {
    x[["vars"]] <- x[[1L]][["vars"]]
    x[["baseline.group"]] <- x[[1L]][["baseline.group"]]
    x[["groups"]] <- unlist(grab(x[vapply(x, inherits, logical(1L), "cem.match")], "groups"))
    x[["w"]] <- get.w.cem.match(x)
  }
  
  if (all(check_if_zero(x[["w"]]))) {
    .err("the supplied cem.match object contains no valid matches")
  }
  
  #Process data and get imp
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  if (is_null(data)) {
    .err("an argument to `data` must be specified with cem.match objects")
  }
  
  #Process imp
  imp <- process_imp(imp, data)
  
  if (is_null(imp) && inherits(x, "cem.match.list") &&
      sum(vapply(x, is_, logical(1L), "cem.match")) != 1L) {
    .err("an argument to `imp` must be specified or the argument to data must be a mids object")
  }
  
  #Process treat
  t.c <- .use_tc_fd(data = data, treat = x[["groups"]], 
                    covs = x[["vars"]])
  treat <- process_treat(t.c[["treat"]], data)
  
  #Process covs
  covs <- t.c[["covs"]]
  
  #Get estimand
  estimand <- ...get("estimand")
  
  #Get method
  method <- "matching"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, covs))
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with cem.match objects")
  }
  
  #Process focal
  focal <- x[["baseline.group"]]
  
  #Process pairwise
  if (get.treat.type(treat) == "binary" && is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
    attr(treat, "treat.type") <- "multinomial"
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with cem.match objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "cem()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("binary")
}

#' @exportS3Method NULL
x2base.weightit <- function(x, ...) {
  #Process weightit
  
  #Process data and get imp
  d.e.in.w <- vapply(c("covs", "exact", "by", "moderator"), function(z) is_not_null(x[[z]]), logical(1L))
  if (any(d.e.in.w)) weightit.data <- do.call("data.frame", unname(x[c("covs", "exact", "by", "moderator")[d.e.in.w]]))
  else weightit.data <- NULL
  
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, weightit.data)
  
  #Process treat
  treat <- process_treat(x[["treat"]], data, weightit.data)
  
  #Process covs
  covs <- x[["covs"]]
  if (is_not_null(covs)) {
    covs <- get_covs_from_formula(data = covs)
  }
  
  #Get estimand
  estimand <- x[["estimand"]]
  
  #Get method
  method <- "weighting"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, weightit.data))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, weightit.data),
                               obj.distance = if (get.treat.type(treat) == "binary") x[["ps"]], 
                               obj.distance.name = "prop.score")
  
  #Process focal
  focal <- x[["focal"]]
  
  #Process pairwise
  if (get.treat.type(treat) == "binary" && is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
    attr(treat, "treat.type") <- "multinomial"
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with weightit objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with weightit objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, weightit.data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["s.weights"]]),
                                 data, weightit.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, weightit.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  discarded <- x[["discarded"]]
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "weightit()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- x[["call"]]
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names)
}

#' @exportS3Method NULL
x2base.designmatch <- function(x, ...) {
  #Process designmatch
  if (all(c("id_1", "id_2") %in% names(x))) {
    .err("Balance cannot currently be checked on a nonbipartite match")
  }
  
  #Process data and get imp
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data)
  
  #Process treat
  t.c <- .use_tc_fd(...get("formula"), data, ...get("treat"), ...get("covs"))
  treat <- process_treat(t.c[["treat"]], data)
  if (is.unsorted(rev(treat))) {
    .wrn("designmatch requires the input data to be sorted by treatment; the data supplied to `bal.tab()` was not, indicating a possible coding error")
  }
  
  #Process covs
  covs <- t.c[["covs"]]
  
  #Get estimand
  estimand <- ...get("estimand")
  
  #Get method
  method <- "matching"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, covs))
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with designmatch objects")
  }
  
  if (get.treat.type(treat) == "binary") {
    if (is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with designmatch objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with designmatch objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "the matching function in designmatch")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- NULL
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("binary")
}

#' @exportS3Method NULL
x2base.mimids <- function(x, ...) {
  #Process mimids
  old_version <- !all(c("object", "models", "approach") %in% names(x))
  models <- if (old_version) x[["models"]][-1L] else x[["models"]]
  
  #Process data and get imp
  m.data <- {
    if (!old_version) .mids_complete(x[["object"]])
    else if (inherits(x[["original.datasets"]], "mids")) .mids_complete(x[["original.datasets"]])
    else .mids_complete(x[["others"]][["source"]])
  }
  
  imp <- m.data[[".imp"]]
  data <- ...get("data")
  
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      # .wrn("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt")
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, m.data)
  
  #Process treat
  treat <- process_treat(unlist(grab(models, "treat")))
  
  #Process covs
  covs <- do.call("rbind", grab(models, "X"))
  covs <- get_covs_from_formula(data = covs)
  
  #Get estimand
  estimand <- models[[1L]][["estimand"]]
  
  #Get method
  method <- "matching"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, m.data))
  
  #Process distance
  m.distance <- unlist(grab(models, "distance"))
  
  if (all(is.na(m.distance))) m.distance <- NULL
  
  distance <- process_distance(...get("distance"), datalist = list(data, m.data),
                               obj.distance = m.distance, 
                               obj.distance.name = "distance")
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with mimids objects")
  }
  
  if (get.treat.type(treat) == "binary") {
    if (is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with mimids objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with mimids objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, m.data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", unlist(grab(models, "s.weights"))),
                                 data, m.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, m.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(table(imp)))
  
  #Process discarded
  discarded <- unlist(grab(models, "discarded"))
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "matchthem()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- NULL
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("imp")
}

#' @exportS3Method NULL
x2base.wimids <- function(x, ...) {
  #Process wimids
  old_version <- !all(c("object", "models", "approach") %in% names(x))
  models <- if (old_version) x[["models"]][-1L] else x[["models"]]
  
  #Process data and get imp
  w.data <- {
    if (!old_version) .mids_complete(x[["object"]])
    else if (inherits(x[["original.datasets"]], "mids")) .mids_complete(x[["original.datasets"]])
    else .mids_complete(x[["others"]][["source"]])
  }
  
  imp <- w.data[[".imp"]]
  data <- ...get("data")
  
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      # .wrn("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt")
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, w.data)
  
  #Process treat
  treat <- process_treat(unlist(grab(models, "treat")))
  
  #Process covs
  covs <- do.call("rbind", grab(models, "covs"))
  covs <- get_covs_from_formula(data = covs)
  
  #Get estimand
  estimand <- unique(unlist(grab(models, "estimand")))
  
  #Get method
  method <- "weighting"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, w.data))
  
  #Process distance
  w.distance <- unlist(grab(models, "ps"))
  if (all(is.na(w.distance))) w.distance <- NULL
  
  distance <- process_distance(...get("distance"), datalist = list(data, w.data),
                               obj.distance = if (get.treat.type(treat) == "binary") w.distance, 
                               obj.distance.name = "prop.score")
  
  #Process focal
  focal <- unique(unlist(grab(models, "focal")))
  
  #Process pairwise
  if (get.treat.type(treat) == "binary" && is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
    attr(treat, "treat.type") <- "multinomial"
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with wimids objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with wimids objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, w.data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", unlist(grab(models, "s.weights"))),
                                 data, w.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, w.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(table(imp)))
  
  #Process discarded
  discarded <- unlist(grab(models, "discarded"))
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "weightthem()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- NULL
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("imp")
}

#' @exportS3Method NULL
x2base.sbwcau <- function(x, ...) {
  #Process sbwcau
  
  #Process data and get imp
  sbw.data <- x[["dat_weights"]][names(x[["dat_weights"]]) != "weights"]
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, sbw.data)
  
  #Process treat
  treat <- process_treat(x[["ind"]], data, sbw.data)
  
  #Process covs
  covs <- reformulate(x[["bal"]][["bal_cov"]]) |>
    get_covs_from_formula(data = sbw.data)
  
  #Get estimand
  estimand <- x[["par"]][["par_est"]]
  
  #Get method
  method <- "weighting"
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, sbw.data))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, sbw.data))
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with sbwcau objects")
  }
  
  if (get.treat.type(treat) == "binary") {
    if (is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with sbwcau objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with sbwcau objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, sbw.data))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data, sbw.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, sbw.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process imp and length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp,
                     original.call.to = "sbw()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("binary")
}

#MSMs wth multiple time points

#' @exportS3Method NULL
x2base.iptw <- function(x, ...) {
  #Process iptw
  stop.method <- ...get("stop.method")
  if (is_null(stop.method) && ...length() > 0L && !nzchar(...names()[1L])) {
    stop.method <- ...elt(1L)
  }
  
  if (is_null(stop.method) && is_not_null(...get("full.stop.method"))) {
    stop.method <- ...get("full.stop.method")
  }
  
  available.stop.methods <- names(x[["psList"]][[1L]][["ps"]])
  
  if (is_null(stop.method)) {
    rule1 <- available.stop.methods
  }
  else if (is.character(stop.method)) {
    rule1 <- available.stop.methods[vapply(tolower(available.stop.methods),
                                           function(z) any(startsWith(z, tolower(stop.method))), logical(1L))]
    if (is_null(rule1)) {
      .wrn(sprintf("`stop.method` should be %s. Using all available stop methods instead",
                   word_list(available.stop.methods, and.or = "or", quotes = 2L)))
      rule1 <- available.stop.methods
    }
  }
  else if (is.numeric(stop.method) && any(stop.method %in% seq_along(available.stop.methods))) {
    if (!all(stop.method %in% seq_along(available.stop.methods))) {
      .wrn(sprintf("there are %s stop methods available, but you requested %s",
                   length(available.stop.methods), 
                   word_list(setdiff(stop.method, seq_along(available.stop.methods)), and.or = "and")))
    }
    rule1 <- available.stop.methods[stop.method %in% seq_along(available.stop.methods)]
  }
  else {
    .wrn("`stop.method` should be %s. Using all available stop methods instead",
         word_list(available.stop.methods, and.or = "or", quotes = 2L))
    rule1 <- available.stop.methods
  }
  
  s <- available.stop.methods[match(tolower(rule1), tolower(available.stop.methods))]
  
  #Process data and get imp
  ps.data <- x[["psList"]][[1L]][["data"]]
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, ps.data)
  
  #Process treat.list
  treat.list <- process_treat.list(grab(x[["psList"]], "treat"), data, ps.data)
  
  #Process covs.list
  if (all_apply(x[["psList"]], function(z) is_not_null(z[["gbm.obj"]][["var.names"]]))) {
    covs.list <- lapply(x[["psList"]], function(z) {
      reformulate(z[["gbm.obj"]][["var.names"]]) |>
        get_covs_from_formula(data = z[["data"]])
    })
  }
  else if (is_not_null(...get("formula.list")) || is_not_null(...get("covs.list"))) {
    covs.list <- ...get("covs.list")
    if (is_not_null(covs.list)) {
      if (!is.list(covs.list) || is.data.frame(covs.list)) {
        .err("`covs.list` must be a list of covariates for which balance is to be assessed at each time point")
      }
      
      if (!all_apply(covs.list, is_mat_like)) {
        .err("each item in `covs.list` must be a data frame")
      }
      
      if (length(covs.list) != length(x[["psList"]])) {
        .err("`covs.list` must have as many entries as time points in the call to `iptw()`")
      }
    }
    
    formula.list <- ...get("formula.list")
    if (is_not_null(formula.list)) {
      if (!is.list(formula.list) || !all_apply(formula.list, rlang::is_formula)) {
        .err("`formula.list` must be a list of formulas identifying the covariates for which balance is to be assessed at each time point")
      }
      
      if (length(formula.list) != length(x[["psList"]])) {
        .err("`formula.list` must have as many entries as time points in the call to `iptw()`")
      }
    }
    
    covs.list <- lapply(seq_along(x[["psList"]]), function(i) {
      .use_tc_fd(formula = formula.list[[i]],
                 data = if_null_then(data, x[["psList"]][[i]][["data"]]),
                 covs = covs.list[[i]],
                 needs.treat = FALSE)[["covs"]]
    })
  }
  else if (all_apply(x[["psList"]], function(z) all(z[["gbm.obj"]][["feature_names"]] %in% names(z[["data"]])))) {
    covs.list <- lapply(x[["psList"]], function(z) {
      reformulate(z[["gbm.obj"]][["feature_names"]]) |>
        get_covs_from_formula(data = z[["data"]])
    })
  }
  else {
    .err('when `version = "xgboost"` in the call to `iptw()` and any variables are categorical, `formula.list` or `covs.list` must be supplied')
  }
  
  #Get estimand
  estimand <- substr(toupper(s), nchar(s) - 2L, nchar(s))
  
  #Get method
  method <- rep_with("weighting", s)
  
  #Process addl.list 
  addl.list <- process_addl.list(...get("addl.list", ...get("addl")),
                                 datalist = list(data, ps.data),
                                 covs.list = covs.list)
  
  #Process distance
  # ntimes <- iptw[["nFits"]]
  # distance.list <- .process_list("distance.list", ...get("distance.list"), ntimes, 
  #                               "the original call to iptw()",
  #                               treat.list,
  #                               covs.list,
  #                               list(data, ps.data))
  # if (is_not_null(distance.list)) {
  #     for (ti in seq_along(distance.list)) {
  #         if (length(s) == 1) {
  #             distance.list[[ti]] <- data.frame(distance[[ti]], prop.score = iptw[["psList"]][[ti]][["ps"]][[s]])
  #         }
  #         else {
  #             distance.list[[ti]] <- data.frame(distance[[ti]], prop.score = iptw[["psList"]][[ti]][["ps"]][s])
  #         }
  #     }
  #     
  # }
  # else {
  #     distance.list <- make_list(ntimes)
  #     for (ti in seq_along(distance.list)) {
  #         if (length(s) == 1) {
  #             distance.list[[ti]] <- data.frame(prop.score = iptw[["psList"]][[ti]][["ps"]][[s]])
  #         }
  #         else {
  #             distance.list[[ti]] <- data.frame(prop.score = iptw[["psList"]][[ti]][["ps"]][s])
  #         }
  #     }
  # }
  # if (is_not_null(distance.list)) distance.list <- lapply(distance.list, function(z) get_covs_from_formula(~z))
  # 
  distance.list <- process_distance.list(...get("distance.list", ...get("distance")),
                                         datalist = list(data, ps.data),
                                         covs.list = covs.list,
                                         obj.distance = lapply(x[["psList"]], function(z) z[["ps"]][, s, drop = FALSE]),
                                         obj.distance.name = {
                                           if (length(s) > 1L) paste.("prop.score", substr(s, 1L, nchar(s) - 4L))
                                           else "prop.score"
                                         })
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with iptw objects")
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with iptw objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with iptw objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat.list[[1L]], covs.list[[1L]],
                             method, addl.data = list(data, ps.data), 
                             stop.method = s)
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["psList"]][[1L]][["sampw"]]),
                                 data, ps.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, ps.data)
  .cluster_check(cluster, treat.list)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(lengths(treat.list)))
  
  #Process discarded
  
  #Process length
  length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("weights"),
                     lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                     imp = imp,
                     original.call.to = "iptw()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat.list)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat.list)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- NULL
  
  #Process output
  X <- initialize_X_msm()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names)
}

#' @exportS3Method NULL
x2base.data.frame.list <- function(x, ...) {
  #Process data and get imp
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data)
  
  #Process treat.list
  treat.list <- process_treat.list(...get("treat.list"), data)
  
  #Process covs.list
  covs.list <- x
  if (is_null(covs.list)) {
    .err("`covs.list` must be specified")
  }
  
  if (!is.list(covs.list) || is.data.frame(covs.list)) {
    .err("`covs.list` must be a list of covariates for which balance is to be assessed at each time point")
  }
  
  if (!all_apply(covs.list, is_mat_like)) {
    .err("each item in `covs.list` must be a data frame")
  }
  
  if (any_apply(covs.list, function(z) is_null(attr(z, "co.names")))) {
    covs.list <- lapply(covs.list, function(z) get_covs_from_formula(data = z))
  }
  
  if (length(treat.list) != length(covs.list)) {
    .err("`treat.list` must be a list of treatment statuses at each time point")
  }
  
  #Get estimand
  estimand <- "ATE"
  
  #Get method
  specified <- setNames(rep.int(FALSE, 1L), "weights")
  if (is_not_null(...get("weights"))) {
    if (!is_(...get("weights"), c("character", "numeric", "data.frame", "list"))) {
      .err("the argument to `weights` must be a vector, list, or data frame of weights or the (quoted) names of variables in `data` that contain weights")
    }
    specified["weights"] <- TRUE
  }
  
  method <- ...get("method")
  if (is_null(method)) {
    method <- if (specified["weights"]) "weighting" else "matching"
  }
  else if (length(method) == 1L) {
    specified.method <- match_arg(method, c("weighting", "matching", "subclassification"))
    if (specified.method == "weighting") {
      method <- if (specified["weights"]) "weighting" else "matching"
    }
    else if (specified["weights"]) {
      .wrn("only weighting is allowed with multiple treatment time points. Assuming weighting instead")
      method <- "weighting"
    }
    else {
      method <- "matching"
    }
  }
  else {
    specified.method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    if (any(specified.method == "subclassification") || specified["subclass"] || specified["match.strata"]) {
      .wrn("only weighting is allowed with multiple treatment time points. Assuming weighting instead")
      method <- "weighting"
    }
    else if (specified["weights"]) {
      .wrn("only weighting is allowed with multiple treatment time points. Assuming weighting instead")
      method <- "weighting"
    }
    # else if (!specified["weights"]) {
    #   #Should never happen
    #   .wrn("multiple methods were specified, but no weights were provided. Providing unadjusted data only")
    #   method <- "matching"
    # }
    else {
      method <- "matching"
    }
  }
  
  #Process addl.list 
  addl.list <- process_addl.list(...get("addl.list", ...get("addl")),
                                 datalist = list(data),
                                 covs.list = covs.list)
  
  #Process distance
  distance.list <- process_distance.list(...get("distance.list", ...get("distance")),
                                         datalist = list(data),
                                         covs.list = covs.list)
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with longitudinal treatments")
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with longitudinal treatments")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with longitudinal treatments")
  }
  
  #Process weights
  weights <- ...get("weights")
  if (is_not_null(weights)) {
    weights <- process_weights(NULL, list(...), treat.list[[1L]], covs.list[[1L]],
                               method, addl.data = list(data))
    method <- attr(weights, "method")
  }
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  .cluster_check(cluster, treat.list)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(lengths(treat.list)))
  
  #Process discarded
  
  #Process length
  length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("weights"),
                     lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                     imp = imp)
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat.list)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat.list)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- NULL
  
  #Process output
  X <- initialize_X_msm()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("msm")
}

#' @exportS3Method NULL
x2base.formula.list <- function(x, ...) {
  
  treat.list <- covs.list <- make_list(length(x))
  
  for (i in seq_along(x)) {
    treat.list[[i]] <- get_treat_from_formula(x[[i]], data = ...get("data"))
    covs.list[[i]] <- get_covs_from_formula(x[[i]], data = ...get("data"))
    names(treat.list)[i] <- attr(treat.list[[i]], "treat.name")
  }
  
  if ("treat.list" %in% ...names()) {
    A <- list(...)
    
    A[["x"]] <- covs.list
    A[["treat.list"]] <- treat.list
    
    return(do.call("x2base.data.frame.list", A))
  }
  
  x2base.data.frame.list(covs.list, treat.list = treat.list, ...)
}

#' @exportS3Method NULL
x2base.CBMSM <- function(x, ...) {
  #Process CBMSM
  ID <- sort(unique(x[["id"]]))
  times <- sort(unique(x[["time"]]))
  x[["data"]] <- x[["data"]][order(x[["id"]], x[["time"]]), , drop = FALSE]
  
  #Process data and get imp
  cbmsm.data <- x[["data"]][x[["time"]] == 1, , drop = FALSE]
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data)
  
  #Process treat.list
  treat.list <- process_treat.list(lapply(times, function(z) x[["treat.hist"]][ID, z]), 
                                   data, cbmsm.data)
  
  #Process covs.list
  covs.list <- make_list(times)
  for (i in seq_along(times)) {
    ti <- times[i]
    cov_i <- get_covs_from_formula(x[["formula"]], data = x[["data"]][x[["time"]] == ti, , drop = FALSE])
    for (co in seq_along(attr(cov_i, "co.names"))) {
      attr(cov_i, "co.names")[[co]][["component"]][attr(cov_i, "co.names")[[co]][["type"]] == "base"] <-
        paste0(attr(cov_i, "co.names")[[co]][["component"]][attr(cov_i, "co.names")[[co]][["type"]] == "base"], "_T", ti)
    }
    names(attr(cov_i, "co.names")) <- vapply(attr(cov_i, "co.names"), function(z) paste(z[["component"]], collapse = ""), character(1L))
    colnames(cov_i) <- names(attr(cov_i, "co.names"))
    covs.list[[i]] <- {
      if (i == 1L) cov_i
      else co.cbind(covs.list[[i - 1L]], cov_i)
    }
  }
  
  #Get estimand
  estimand <- "ATE"
  
  #Get method
  method <- "weighting"
  
  #Process addl.list 
  addl.list <- process_addl.list(...get("addl.list", ...get("addl")),
                                 datalist = list(data, cbmsm.data),
                                 covs.list = covs.list)
  
  #Process distance
  distance.list <- process_distance.list(...get("distance.list", ...get("distance")),
                                         datalist = list(data, cbmsm.data),
                                         covs.list = covs.list, obj.distance = x[["fitted.values"]],
                                         obj.distance.name = "prop.score")
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with CBMSM objects")
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with CBMSM objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with CBMSM objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat.list[[1L]], covs.list[[1L]], method,
                             addl.data = list(data, cbmsm.data))
  method <- attr(weights, "method")
  
  #Process s.weights
  if (is_not_null(...get("s.weights"))) {
    .err("sampling weights are not allowed with CBMSM objects")
  }
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, cbmsm.data)
  .cluster_check(cluster, treat.list)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(lengths(treat.list)))
  
  #Process discarded
  
  #Process length
  length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("weights"),
                     lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                     imp = imp,
                     original.call.to = "CBMSM()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat.list)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat.list)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- x[["call"]]
  
  #Process output
  X <- initialize_X_msm()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names)
}

#' @exportS3Method NULL
x2base.weightitMSM <- function(x, ...) {
  #Process weightitMSM
  
  #Process data and get imp
  weightitMSM.data <- x[["data"]]
  
  d.e.in.w <- vapply(c("exact", "by", "moderator"), function(z) is_not_null(x[[z]]), logical(1L))
  weightitMSM.data2 <- {
    if (any(d.e.in.w)) do.call("data.frame", unname(x[c("exact", "by", "moderator")[d.e.in.w]]))
    else NULL
  }
  
  imp <- ...get("imp")
  data <- ...get("data")
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, weightitMSM.data, weightitMSM.data2)
  
  #Process treat.list
  treat.list <- process_treat.list(x[["treat.list"]],
                                   data, weightitMSM.data, weightitMSM.data2)
  #Process covs.list
  covs.list <- lapply(x[["covs.list"]], function(z) get_covs_from_formula(data = z))
  
  #Get estimand
  estimand <- x[["estimand"]]
  
  #Get method
  method <- "weighting"
  
  #Process addl.list 
  addl.list <- process_addl.list(...get("addl.list", ...get("addl")), 
                                 datalist = list(data, weightitMSM.data,
                                                 weightitMSM.data2),
                                 covs.list = covs.list)
  
  #Process distance
  # ntimes <- length(covs.list)
  # distance.list <- .process_list("distance.list", ...get("distance.list"), ntimes, 
  #                               "the original call to weightitMSM()",
  #                               treat.list,
  #                               covs.list,
  #                               list(data, weightitMSM.data,
  #                                    weightitMSM.data2))
  # if (is_not_null(distance.list)) distance.list <- lapply(seq_along(distance.list), function(z) data.frame(distance.list[[z]], prop.score = weightitMSM[["ps.list"]][[z]]))
  # else if (is_not_null(weightitMSM[["ps.list"]])) distance.list <- lapply(seq_along(weightitMSM[["ps.list"]]), function(z) data.frame(prop.score = weightitMSM[["ps.list"]][[z]]))
  # else distance.list <- NULL
  # if (is_not_null(distance.list)) distance.list <- lapply(distance.list, function(z) get_covs_from_formula(~z))
  distance.list <- process_distance.list(...get("distance.list", ...get("distance")),
                                         datalist = list(data, weightitMSM.data, weightitMSM.data2),
                                         covs.list = covs.list, obj.distance = x[["ps.list"]],
                                         obj.distance.name = "prop.score")
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with weightitMSM objects")
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with weightitMSM objects")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with weightitMSM objects")
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat.list[[1L]], covs.list[[1L]], method, 
                             addl.data = list(data, weightitMSM.data, weightitMSM.data2))
  method <- attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["s.weights"]]),
                                 data, weightitMSM.data, weightitMSM.data2)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, weightitMSM.data, weightitMSM.data2)
  .cluster_check(cluster, treat.list)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(lengths(treat.list)))
  
  #Process discarded
  
  #Process length
  length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("weights"),
                     lists = c("treat.list", "covs.list", "addl.list", "distance.list"),
                     imp = imp,
                     original.call.to = "weightitMSM()")
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat.list)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat.list)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- x[["call"]]
  
  #Process output
  X <- initialize_X_msm()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names)
}

#' @exportS3Method NULL
x2base.default <- function(x, ...) {
  
  if (!is.list(x)) {
    .err("the input object must be an appropriate list, data.frame, formula, or the output of one of the supported packages")
  }
  
  if (...length() > 0L && (is_null(...names()) || !all(nzchar(...names())))) {
    .err("all arguments to `...` must be named")
  }
  
  Q <- list(treat = list(name = c("treat", "tr"), 
                         type = c("numeric", "character", "factor", "logical")),
            treat.list = list(name = c("treat.list", "treat", "tr"),
                              type = c("list", "data.frame")),
            covs = list(name = c("covs", "covariates", "x"), 
                        type = c("data.frame")),
            covs.list = list(name = c("covs.list", "covs", "covariates"),
                             type = c("list")),
            formula = list(name = c("formula", "form"), 
                           type = c("formula")),
            formula.list = list(name = c("formula.list", "formula", "form"),
                                type = c("list")),
            data = list(name = c("data"),
                        type = c("data.frame", "mids")),
            weights = list(name = c("weights", "w", "wts"),
                           type = c("data.frame", "matrix", "numeric")),
            distance = list(name = c("distance", "distance.list", "ps", "pscore", "p.score", "propensity.score"),
                            type = c("data.frame", "matrix", "numeric", "list")),
            subclass = list(name = c("subclass", "strata"),
                            type = c("factor", "character", "numeric")),
            match.strata = list(name = c("match.strata"),
                                type = c("factor", "character", "numeric")),
            estimand = list(name = c("estimand", "target", "att", "ate"),
                            type = c("character", "numeric", "logical")),
            s.weights = list(name = c("s.weights", "sw", "sweights", "sampw"),
                             type = c("numeric")),
            focal = list(name = c("focal", "treatATT"), 
                         type = c("character", "numeric")),
            call = list(name = c("call"),
                        type = c("call")))
  
  P <- make_list(names(Q))
  names(x) <- tolower(names(x))
  
  #Make a new list, P, containing the extracted components of obj; P acts as
  #new object for future steps.
  for (i in setdiff(names(Q), ...names())) {
    for (j in Q[[i]][["name"]]) {
      if (is_null(P[[i]]) && is_not_null(x[[j]])) {
        which.type <- vapply(Q[[i]][["type"]], function(z) is_(x[[j]], z), logical(1L))
        if (any(which.type)) {
          P[[i]] <- x[[j]]
          attr(P[[i]], "name") <- j
          attr(P[[i]], "type") <- Q[[i]][["type"]][which.type]
        }
      }
    }
  }
  
  #treat OK
  
  #treat.list
  if (is_not_null(P[["treat.list"]]) &&
      !all_apply(P[["treat.list"]], function(z) any_apply(Q[["treat"]][["type"]],
                                                          function(c) is_(z, c)))) {
    P[["treat.list"]] <- NULL
  }
  
  #covs 
  if (is_not_null(P[["covs"]])) {
    P[["covs"]] <- as.data.frame(P[["covs"]])
  }
  
  #covs.list
  if (is_not_null(P[["covs.list"]]) &&
      !all_apply(P[["covs.list"]], function(z) any_apply(Q[["covs"]][["type"]],
                                                         function(c) is_(z, c)))) {
    P[["covs.list"]] <- NULL
  }
  
  #formula
  
  #formula.list
  if (is_not_null(P[["formula.list"]]) &&
      !all_apply(P[["formula.list"]], function(z) any_apply(Q[["formula"]][["type"]],
                                                            function(c) is_(z, c)))) {
    P[["formula.list"]] <- NULL
  }
  
  #data
  #model (only to extract data)
  if (is_null(P[["data"]]) && is_not_null(x[["model"]]) && utils::hasName(x[["model"]], "data")) {
    P[["data"]] <- x[["model"]][["data"]]
  }
  
  #weights
  
  #distance
  if (is_not_null(P[["distance"]])) {
    if (is.list(P[["distance"]]) && !is.data.frame(P[["distance"]])) {
      if (!all_apply(P[["distance"]], function(z) any_apply(Q[["distance"]][["type"]],
                                                            function(c) is_(z, c)))) {
        P[["distance"]] <- NULL
      }
    }
    else if (is.numeric(P[["distance"]])) {
      P[["distance"]] <- setNames(data.frame(P[["distance"]]),
                                  if_null_then(attr(P[["distance"]], "name", TRUE), "distance"))
    }
    else {
      P[["distance"]] <- as.data.frame(P[["distance"]])
    }
  }
  
  #distance.list
  if (is_not_null(P[["distance.list"]]) &&
      !all_apply(P[["distance.list"]], function(z) any_apply(Q[["distance"]][["type"]],
                                                             function(c) is_(z, c)))) {
    P[["distance.list"]] <- NULL
  }
  
  #subclass
  if (is_not_null(P[["subclass"]])) {
    P[["subclass"]] <- factor(P[["subclass"]])
  }
  
  #match.strata
  if (is_not_null(P[["match.strata"]])) {
    P[["match.strata"]] <- factor(P[["match.strata"]])
  }
  
  #estimand
  if (is_not_null(P[["estimand"]])) {
    estimand.name <- attr(P[["estimand"]], "name", TRUE)
    
    P[["estimand"]] <- {
      if (is_not_null(estimand.name) && toupper(estimand.name) == "ATT") {
        if (as.numeric(P[["estimand"]]) == 0) "ATE" else "ATT"
      }
      else if (is_not_null(estimand.name) && toupper(estimand.name) == "ATE") {
        if (as.numeric(P[["estimand"]]) == 0) "ATT" else "ATE"
      }
      else if (tolower(P[["estimand"]]) %in% c("att", "treat", "treated", "tr", "t", "atet")) {
        "ATT"
      }
      else if (tolower(P[["estimand"]]) %in% c("ate", "all")) {
        "ATE"
      }
      else if (tolower(P[["estimand"]]) %in% c("atc", "control", "untreated", "u", "c", "ctrl", "atu", "atec", "ateu")) {
        "ATC"
      }
      else {
        NULL
      }
    }
  }
  
  #s.weights
  
  #focal
  
  #call
  
  msm <- is_not_null(P[["treat.list"]]) || is_not_null(P[["formula.list"]])
  
  if (msm) {
    .x2base_default_msm(P, ...)
  }
  else {
    .x2base_default_point(P, ...)
  }
}

.x2base_default_point <- function(obj, ...) {
  #Process data and get imp
  imp <- ...get("imp")
  data <- ...get("data")
  
  o.data <- obj[["data"]]
  if (is_null(o.data) && is_not_null(obj[["model"]]) && utils::hasName(obj[["model"]], "data")) {
    o.data <- obj[["model"]][["data"]]
  }
  if (inherits(o.data, "mids")) {
    o.data <- .mids_complete(o.data)
  }
  
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, o.data)
  
  #Process treat
  treat <- ...get("treat")
  if (is_null(treat)) {
    formula <- ...get("formula")
    if (is_not_null(formula) && rlang::is_formula(formula, lhs = TRUE)) {
      treat <- get_treat_from_formula(formula, data = if_null_then(data, o.data))
    }
    else if (is_not_null(obj[["treat"]])) {
      treat <- obj[["treat"]]
    }
    else {
      formula.obj <- obj[["formula"]]
      
      if (is_not_null(formula.obj) && rlang::is_formula(formula.obj, lhs = TRUE)) {
        treat <- get_treat_from_formula(formula.obj, data = if_null_then(data, o.data))
      }
    }
  }
  treat <- process_treat(treat, data, o.data)
  
  #Process covs
  covs <- ...get("covs")
  if (is_null(covs)) {
    formula <- ...get("formula")
    if (is_not_null(formula) && rlang::is_formula(formula)) {
      covs <- get_covs_from_formula(formula, data = if_null_then(data, o.data))
    }
    else if (is_not_null(obj[["covs"]])) {
      covs <- obj[["covs"]]
    }
    else {
      formula.obj <- obj[["formula"]]
      
      if (is_not_null(formula.obj) && rlang::is_formula(formula.obj)) {
        covs <- get_covs_from_formula(formula.obj, data = if_null_then(data, o.data))
      }
    }
  }
  
  if (is_null(covs)) {
    .err("`covs` data.frame must be specified")
  }
  
  if (is_null(attr(covs, "co.names", TRUE))) {
    if (is.matrix(covs)) covs <- as.data.frame.matrix(covs)
    covs <- get_covs_from_formula(data = covs)
  }
  
  #Get estimand
  estimand <- ...get("estimand", obj[["estimand"]])
  
  #Get method
  specified <- setNames(rep.int(FALSE, 3L), c("match.strata", "subclass", "weights"))
  
  for (i in names(specified)) {
    specified[i] <- is_not_null(...get(i, obj[[i]]))
  }
  
  .specified_method <- character()
  .specified_args <- character()
  .assuming <- character()
  .ignoring <- character()
  .not_present <- character()
  .using <- character()
  
  method <- ...get("method")
  if (is_null(method)) {
    if (any(specified)) {
      .using <- names(specified)[specified][1L]
      method <- switch(.using,
                       match.strata = "matching",
                       subclass = "subclassification",
                       weights = "weighting")
      
      if (sum(specified) > 1L) {
        .specified_args <- names(specified)[specified]
        .assuming <- method
        .ignoring <- setdiff(names(specified)[specified], .using)
      }
    }
    else {
      method <- "matching"
    }
  }
  else {
    .specified_method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    
    if (length(method) == 1L) {
      if (.specified_method == "weighting") {
        if (specified["weights"]) {
          method <- "weighting"
          .using <- "weights"
          
          if (sum(specified) > 1L) {
            .specified_args <- names(specified)[specified]
            .ignoring <- setdiff(names(specified)[specified], .using)
          }
        }
        else if (any(specified)) {
          .using <- names(specified)[specified][1L]
          method <- switch(.using,
                           match.strata = "matching",
                           subclass = "subclassification",
                           weights = "weighting")
          .assuming <- method
          .not_present <- "weights"
        }
        else {
          method <- "matching"
        }
      }
      else if (.specified_method == "matching") {
        if (specified["match.strata"]) {
          method <- "matching"
          .using <- "match.strata"
          
          if (sum(specified) > 1L) {
            .specified_args <- names(specified)[specified]
            .ignoring <- setdiff(names(specified)[specified], .using)
          }
        }
        else if (specified["weights"]) {
          method <- "matching"
          .using <- "weights"
          
          if (sum(specified) > 1L) {
            .specified_args <- names(specified)[specified]
            .ignoring <- setdiff(names(specified)[specified], .using)
          }
        }
        else if (specified["subclass"]) {
          method <- "subclassification"
          .using <- "subclass"
          .not_present <- c("weights", "match.strata")
          .assuming <- method
        }
        else {
          method <- "matching"
        }
      }
      else if (.specified_method == "subclassification") {
        if (specified["subclass"]) {
          method <- "subclassification"
          .using <- "subclass"
          
          if (sum(specified) > 1L) {
            .specified_args <- names(specified)[specified]
            .ignoring <- setdiff(names(specified)[specified], .using)
          }
        }
        else if (any(specified)) {
          .using <- names(specified)[specified][1L]
          method <- switch(.using,
                           match.strata = "matching",
                           subclass = "subclassification",
                           weights = "weighting")
          .assuming <- method
          .not_present <- "subclass"
        }
        else {
          method <- "matching"
        }
      }
    }
    else {
      if (specified["subclass"] || any(.specified_method == "subclassification")) {
        .err("subclassification cannot be specified along with other methods")
      }
      
      if (specified["match.strata"]) {
        .err("only weights can be specified with multiple methods")
      }
      
      if (specified["weights"]) {
        method <- .specified_method
        .using <- "weights"
      }
      else {
        .wrn("multiple methods were specified, but no weights were provided. Providing unadjusted data only")
        method <- "matching"
      }
    }
  }
  
  .m <- NULL
  if (is_not_null(.using)) {
    if (is_not_null(.specified_args) && is_not_null(.ignoring)) {
      if (is_not_null(.assuming)) {
        .m <- sprintf("%s specified. Assuming %s and using %s and ignoring %s",
                      word_list(.specified_args, and.or = "and", is.are = TRUE, quotes = "`"),
                      add_quotes(.assuming),
                      word_list(.using, and.or = "and", quotes = "`"),
                      word_list(.ignoring, and.or = "and", quotes = "`"))
      }
      else {
        .m <- sprintf("%s specified. Using %s and ignoring %s",
                      word_list(.specified_args, and.or = "and", is.are = TRUE, quotes = "`"),
                      word_list(.using, and.or = "and", quotes = "`"),
                      word_list(.ignoring, and.or = "and", quotes = "`"))
      }
    }
    else if (is_not_null(.specified_method) && is_not_null(.not_present) && is_not_null(.assuming)) {
      .m <- sprintf("`method = %s` is specified, but %s was not supplied. Assuming %s and using %s instead",
                    add_quotes(.specified_method, 2L),
                    add_quotes(.not_present, "`"),
                    add_quotes(.assuming),
                    word_list(.using, and.or = "and", quotes = "`"))
    }
  }
  
  if (is_not_null(.m)) {
    .msg(.m)
  }
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, o.data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, o.data, covs),
                               obj.distance = obj[["distance"]], 
                               obj.distance.name = attr(obj[["distance"]], "name", TRUE))
  
  #Process focal
  focal <- process_focal(...get("focal", obj[["focal"]]), treat)
  
  if (get.treat.type(treat) == "binary") {
    if (is_null(focal) && is_not_null(estimand)) {
      focal <- switch(toupper(estimand), 
                      "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                      "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                      NULL)
    }
    
    #Process pairwise
    if (is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
      attr(treat, "treat.type") <- "multinomial"
    }
  }
  
  #Process subclass
  if ("subclass" %in% .using) {
    subclass <- .process_vector(...get("subclass", obj[["subclass"]]), 
                                datalist = list(data, o.data),
                                name = "subclass", 
                                which = "subclass membership",
                                missing.okay = TRUE) |>
      factor()
    weights <- NULL
  }
  
  #Process match.strata
  if ("match.strata" %in% .using) {
    match.strata <- .process_vector(...get("match.strata", obj[["match.strata"]]), 
                                    datalist = list(data, o.data),
                                    name = "match.strata", 
                                    which = "matching strata membership",
                                    missing.okay = TRUE)
    
    weights <- data.frame(weights = strata2weights(match.strata,
                                                   treat = treat,
                                                   estimand = estimand,
                                                   focal = focal))
  }
  
  #Process weights
  if ("weights" %in% .using) {
    weights <- process_weights(obj, list(...), treat, covs, method, addl.data = list(data, o.data))
    method <- attr(weights, "method")
  }
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", obj[["s.weights"]]), data, o.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, o.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  discarded <- ...get("discarded", obj[["discarded"]])
  
  #Process length
  length_imp_process(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("covs", "weights", "distance", "addl"),
                     imp = imp)
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs) || anyNA(addl)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  
  #Process output
  X <- initialize_X()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names)
}

.x2base_default_msm <- function(obj, ...) {
  #Process data and get imp
  imp <- ...get("imp")
  data <- ...get("data")
  
  o.data <- obj[["data"]]
  if (is_null(o.data) && is_not_null(obj[["model"]]) && utils::hasName(obj[["model"]], "data")) {
    o.data <- obj[["model"]][["data"]]
  }
  if (inherits(o.data, "mids")) {
    o.data <- .mids_complete(o.data)
  }
  
  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      if (is_null(imp)) imp <- data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }
  
  #Process imp
  imp <- process_imp(imp, data, o.data)
  
  #Process treat.list
  treat.list <- process_treat.list(...get("treat.list"), data)
  
  treat.list <- ...get("treat.list")
  if (is_null(treat.list)) {
    formula.list <- ...get("formula.list")
    if (is_not_null(formula.list) && all_apply(formula.list, rlang::is_formula, lhs = TRUE)) {
      treat.list <- make_list(length(formula.list))
      
      for (i in seq_along(formula.list)) {
        treat.list[[i]] <- get_treat_from_formula(formula.list[[i]], data = if_null_then(data, o.data))
        names(treat.list)[i] <- attr(treat.list[[i]], "treat.name")
      }
    }
    else if (is_not_null(obj[["treat.list"]])) {
      treat.list <- obj[["treat.list"]]
    }
    else {
      formula.list.obj <- obj[["formula.list.obj"]]
      
      if (is_not_null(formula.list.obj) && all_apply(formula.list.obj, rlang::is_formula, lhs = TRUE)) {
        treat.list <- make_list(length(formula.list.obj))
        
        for (i in seq_along(formula.list.obj)) {
          treat.list[[i]] <- get_treat_from_formula(formula.list.obj[[i]], data = if_null_then(data, o.data))
          names(treat.list)[i] <- attr(treat.list[[i]], "treat.name")
        }
      }
    }
  }
  treat.list <- process_treat.list(treat.list, data, o.data)
  
  #Process covs.list
  covs.list <- ...get("covs.list")
  if (is_null(covs.list)) {
    formula.list <- ...get("formula.list")
    if (is_not_null(formula.list) && all_apply(formula.list, rlang::is_formula)) {
      covs.list <- make_list(length(formula.list))
      
      for (i in seq_along(formula.list)) {
        covs.list[[i]] <- get_covs_from_formula(formula.list[[i]], data = if_null_then(data, o.data))
      }
    }
    else if (is_not_null(obj[["covs.list"]])) {
      covs.list <- obj[["covs.list"]]
    }
    else {
      formula.list.obj <- obj[["formula.list.obj"]]
      
      if (is_not_null(formula.list.obj) && all_apply(formula.list.obj, rlang::is_formula)) {
        covs.list <- make_list(length(formula.list.obj))
        
        for (i in seq_along(formula.list.obj)) {
          covs.list[[i]] <- get_covs_from_formula(formula.list.obj[[i]], data = if_null_then(data, o.data))
        }
      }
    }
  }
  
  if (is_null(covs.list)) {
    .err("`covs.list` must be specified")
  }
  
  if (!is.list(covs.list) || is.data.frame(covs.list)) {
    .err("`covs.list` must be a list of covariates for which balance is to be assessed at each time point")
  }
  
  if (!all_apply(covs.list, is_mat_like)) {
    .err("each item in `covs.list` must be a data frame")
  }
  
  if (any_apply(covs.list, function(z) is_null(attr(z, "co.names")))) {
    covs.list <- lapply(covs.list, function(z) get_covs_from_formula(data = z))
  }
  
  if (length(treat.list) != length(covs.list)) {
    .err("`treat.list` must be a list of treatment statuses at each time point")
  }
  
  #Get estimand
  estimand <- "ATE"
  
  #Get method
  specified <- setNames(rep.int(FALSE, 1L), "weights")
  
  for (i in names(specified)) {
    specified[i] <- is_not_null(...get(i, obj[[i]]))
  }
  
  .using <- character()
  
  method <- ...get("method")
  if (is_null(method)) {
    method <- if (specified["weights"]) "weighting" else "matching"
    if (any(specified)) {
      .using <- names(specified)[specified][1L]
    }
  }
  else {
    specified.method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    
    method <- if (specified["weights"]) "weighting" else "matching"
    
    if (any(specified)) {
      .using <- names(specified)[specified][1L]
    }
    
    if (any(specified.method == "subclassification")) {
      if (specified["weights"]) {
        .wrn("only weighting is allowed with multiple treatment time points. Assuming weighting instead")
      }
      else {
        .wrn("only weighting is allowed with multiple treatment time points. Providing unadjusted data only")
      }
    }
  }
  
  #Process addl.list 
  addl.list <- process_addl.list(...get("addl.list", ...get("addl")),
                                 datalist = list(data, o.data),
                                 covs.list = covs.list)
  
  #Process distance
  distance.list <- process_distance.list(...get("distance.list", ...get("distance", obj[["distance.list"]])),
                                         datalist = list(data, o.data),
                                         covs.list = covs.list)
  
  #Process focal
  focal <- ...get("focal")
  if (is_not_null(focal)) {
    .err("`focal` is not allowed with longitudinal treatments")
  }
  
  #Process subclass
  if (is_not_null(...get("subclass"))) {
    .err("subclasses are not allowed with longitudinal treatments")
  }
  
  #Process match.strata
  if (is_not_null(...get("match.strata"))) {
    .err("matching strata are not allowed with longitudinal treatments")
  }
  
  #Process weights
  if ("weights" %in% .using) {
    weights <- process_weights(obj, list(...), treat.list[[1L]], covs.list[[1L]],
                               method, addl.data = list(data, o.data))
    method <- attr(weights, "method")
  }
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", obj[["s.weights"]]),
                                 data, o.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, o.data)
  .cluster_check(cluster, treat.list)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(lengths(treat.list)))
  
  #Process discarded
  
  #Process length
  length_imp_process(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                     data.frames = c("weights"),
                     lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                     imp = imp)
  
  #Process stats and thresholds
  if (!check_if_call_from_fun(bal.plot)) {
    stats <- process_stats(...get("stats"), treat = treat.list)
    type <- attr(stats, "type")
    thresholds <- ...get("thresholds", list())
    
    if (is_not_null(thresholds)) {
      thresholds <- process_thresholds(thresholds, c(stats, setdiff(all_STATS(type), stats)))
      if (!all(names(thresholds) %in% stats)) {
        stats <- unique(c(stats, names(thresholds)))
      }
    }
    
    for (s in all_STATS(type)) {
      #If disp.stat is TRUE, add stat to stats
      if (isTRUE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- unique(c(stats, s))
      }
      else if (isFALSE(...get(STATS[[s]][["disp_stat"]]))) {
        stats <- setdiff(stats, s)
      }
      
      #Process and check thresholds
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
          .wrn(sprintf("%s must be between %s; ignoring %s",
                       add_quotes(STATS[[s]][["threshold"]], "`"),
                       word_list(STATS[[s]][["threshold_range"]]),
                       add_quotes(STATS[[s]][["threshold"]], "`")))
        }
      }
    }
    
    stats <- process_stats(stats, treat = treat.list)
    
    #Get s.d.denom
    s.d.denom <- ...get("s.d.denom")
  }
  
  #Missing values warning
  if (anyNA(covs.list, recursive = TRUE) || anyNA(addl.list, recursive = TRUE)) {
    .wrn("missing values exist in the covariates. Displayed values omit these observations")
  }
  
  #Get call
  call <- NULL
  
  #Process output
  X <- initialize_X_msm()
  X.names <- names(X)
  
  for (i in X.names) {
    X[[i]] <- get0(i, inherits = FALSE)
  }
  
  X <- subset_X(X, subset)
  
  setNames(X[X.names], X.names) |>
    set_class("msm")
}
