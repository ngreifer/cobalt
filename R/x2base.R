#Functions to convert object to base.bal.tab input

x2base <- function(x, ...) {
  UseMethod("x2base")
}

#Reject arguments that do not apply to the input class. `.what` completes the
#sentence "... not allowed with <.what>" and carries its own cli markup, so it is
#pasted into the template rather than interpolated -- cli does not re-render markup
#found inside an interpolated value.
#
#The subject is per-argument so the message names the concept the user supplied
#rather than the internal slot. `...` is forwarded rather than materialized so only
#the arguments actually checked are forced.
#Resolve `data` and `imp` from the arguments `bal.tab()` forwards. A `mids` object
#given as `data` is completed to a long data frame, and its `.imp` column supplies
#`imp` unless the caller named one; anything else that is not a data frame is
#discarded.
#
#`.datalist` holds the data frames the input object carries, against which a
#character `imp` is resolved. `.imp` is for the classes whose imputation index comes
#from the object itself rather than from the user; supplying it replaces the `imp`
#argument entirely.
#
#As elsewhere, `...` is forwarded rather than materialized and the named formals are
#dot-prefixed, since `data` and `imp` are themselves elements of `...`.
.x2base_data <- function(..., .datalist = list(), .imp) {
  imp <- if (missing(.imp)) ...get("imp") else .imp
  data <- ...get("data")

  if (is_not_null(data)) {
    if (inherits(data, "mids")) {
      data <- .mids_complete(data)
      imp <- imp %or% data[[".imp"]]
    }
    else if (!is.data.frame(data)) {
      data <- NULL
    }
  }

  list(data = data,
       imp = do.call("process_imp", c(list(imp, data), .datalist)))
}

.reject_args <- function(.args, .what, ...) {
  for (a in .args) {
    if (is_null(...get(a))) next

    subject <- switch(a,
                      "subclass" = "subclasses are",
                      "match.strata" = "matching strata are",
                      "s.weights" = "sampling weights are",
                      sprintf("{.arg %s} is", a))

    arg::err(paste(subject, "not allowed with", .what))
  }
}

#`estimand` and `focal` name a treatment group to target. A censoring model's target is
#the full at-risk sample, so neither has anything to say; they are ignored rather than
#rejected, because a `weightit` object may carry an `estimand` slot it never used.
.reject_cens_args <- function(...) {
  supplied <- c("estimand", "focal")[c(is_not_null(...get("estimand")),
                                       is_not_null(...get("focal")))]

  if (is_not_null(supplied)) {
    arg::wrn("{.arg {supplied}} {?does/do} not apply to a censoring indicator, whose target is the full at-risk sample; ignoring {?it/them}")
  }

  #Subclassification is not a way of estimating censoring weights.
  .reject_args(c("subclass", "match.strata"), "a censoring indicator", ...)
}

#The tail every method shares: expand anything supplied for a single imputation,
#assemble `X` from the locals the method has built, resolve the requested statistics,
#warn about missing covariate values, and apply `subset`.
#
#`.length` holds `length_imp_process()`'s arguments, naming which of the method's
#locals must be length-compatible and how each expands. All formals are dot-prefixed
#because `...` carries the user's own arguments, and a named element of `...` matches a
#formal of the same name even when the caller supplied that formal positionally. `...`
#must reach `process_stats_and_thresholds()` unforced: `bal.plot()` passes its own lazy
#dots through, and materializing them changes when an erroring unused argument fails.
.finish_X <- function(.length, ..., .call = NULL, .msm = FALSE,
                      .env = parent.frame()) {
  .get <- function(nm) get0(nm, envir = .env, inherits = FALSE)

  checked <- unlist(.length[c("vectors", "data.frames", "lists")], use.names = FALSE)

  expanded <- do.call("length_imp_process",
                      c(list(setNames(lapply(checked, .get), checked)), .length,
                        list(imp = .get("imp"))))

  X <- if (.msm) initialize_X_msm() else initialize_X()

  #`X`'s slots are named after the locals that fill them, taking the expanded value
  #where there is one. A slot the method never assigned stays NULL, which is what
  #every consumer expects of an absent one.
  for (i in names(X)) {
    X[[i]] <- if (i %in% checked) expanded[[i]] else .get(i)
  }

  .stats <- process_stats_and_thresholds(X[[if (.msm) "treat.list" else "treat"]], ...)

  X[["stats"]] <- .stats[["stats"]]
  X[["thresholds"]] <- .stats[["thresholds"]]
  X[["s.d.denom"]] <- .stats[["s.d.denom"]]
  X[["call"]] <- .call

  #A longitudinal method holds lists of data frames, which `anyNA()` only descends
  #into when asked.
  missings <- {
    if (.msm) anyNA(X[["covs.list"]], recursive = TRUE) ||
      anyNA(X[["addl.list"]], recursive = TRUE)
    else anyNA(X[["covs"]]) || anyNA(X[["addl"]])
  }

  if (missings) {
    arg::wrn("missing values exist in the covariates. Displayed values omit these observations")
  }

  X <- subset_X(X, .get("subset"))

  setNames(X[names(X)], names(X))
}

#' @exportS3Method NULL
x2base.matchit <- function(x, ...) {
  #Process matchit
  
  #Process data and get imp
  m.data <- if (NROW(x[["model"]][["data"]]) == length(x[["treat"]])) x[["model"]][["data"]]
  .d <- .x2base_data(..., .datalist = list(m.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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

      order <- setNames(.attr(x[["model"]][["terms"]], "order"),
                        .attr(x[["model"]][["terms"]], "term.labels"))
      assign <- setNames(.attr(x[["X"]], "assign"), colnames(x[["X"]]))
      assign1 <- assign[assign %in% which(order == 1)] #Just main effects
      
      dataClasses <- .attr(x[["model"]][["terms"]], "dataClasses")
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
    method <- ...get("method", "subclassification")

    method <- arg::match_arg(method, c("weighting", "subclassification"))
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
  
  #Reject arguments that do not apply
  .reject_args("subclass", "{.cls matchit} objects", ...)
  subclass <- switch(method,
                     "subclassification" = as.factor(x[["subclass"]]),
                     NULL)
  
  #Reject arguments that do not apply
  .reject_args("match.strata", "{.cls matchit} objects", ...)
  
  #Process weights
  if (is_not_null(x[["weights"]]) && !all_the_same(x[["weights"]])) {
    weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, m.data))
    method <- .attr(weights, "method")
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
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun matchit}"),
            ...,
            .call = x[["call"]])
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
  
  rule1 <- .process_stop_method(stop.method, available.stop.methods)
  
  s <- available.stop.methods[match(tolower(rule1), tolower(available.stop.methods))]
  
  #Process data and get imp
  ps.data <- x[["data"]]
  .d <- .x2base_data(..., .datalist = list(ps.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]

  #Process treat
  treat <- process_treat(x[["treat"]], data, ps.data)
  
  #Process covs
  if (is_not_null(x[["gbm.obj"]][["var.names"]])) {
    covs <- reformulate(x[["gbm.obj"]][["var.names"]]) |>
      get_covs_from_formula(data = ps.data)
  }
  else if (is_not_null(...get("formula")) || is_not_null(...get("covs"))) {
    t.c <- .use_tc_fd(formula = ...get("formula"),
                      data = data %or% ps.data,
                      covs = ...get("covs"),
                      needs.treat = FALSE)
    
    #Process covs
    covs <- t.c[["covs"]]
  }
  #`is_not_null()` is required as well as the `%in%` test: `all()` of an empty set
  #is TRUE, so a model object with no `feature_names` (e.g., an xgboost 3.x
  #booster, which exposes them elsewhere) would otherwise enter this branch and
  #fail inside `reformulate()` instead of reaching the informative error below.
  else if (is_not_null(x[["gbm.obj"]][["feature_names"]]) &&
           all(x[["gbm.obj"]][["feature_names"]] %in% names(ps.data))) {
    covs <- reformulate(x[["gbm.obj"]][["feature_names"]]) |>
      get_covs_from_formula(data = ps.data)
  }
  else {
    arg::err('when {.code version = "xgboost"} in the call to {.fn ps} and any variables are categorical, {.arg formula} or {.arg covs} must be supplied')
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
  .reject_args("focal", "{.cls ps} objects", ...)
  
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
  
  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.cls ps} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, ps.data), 
                             stop.method = s, estimand = estimand)
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["sampw"]]),
                                 data, ps.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, ps.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun ps}"),
            ...)
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
  
  rule1 <- .process_stop_method(stop.method, available.stop.methods)
  
  s <- available.stop.methods[match(tolower(rule1), tolower(available.stop.methods))]
  
  #Process data and get imp
  mnps.data <- x[["data"]]
  .d <- .x2base_data(..., .datalist = list(mnps.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
  #Process treat
  treat <- process_treat(x[["treatVar"]], data, mnps.data)
  
  #Process covs
  .v <- x[["balanceVars"]] %or% x[["psList"]][[1L]][["gbm.obj"]][["var.names"]]
  
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
  
  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.cls mnps} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, mnps.data), 
                             stop.method = s)
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["sampw"]]),
                                 data, mnps.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, mnps.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun mnps}"),
            ...,
            .call = NULL)
}

#' @exportS3Method NULL
x2base.ps.cont <- function(x, ...) {
  #Process data and get imp
  ps.data <- x[["data"]]
  .d <- .x2base_data(..., .datalist = list(ps.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
  
  #Reject arguments that do not apply
  .reject_args(c("focal", "subclass", "match.strata"), "{.cls ps.cont} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, ps.data))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["sampw"]]),
                                 data, ps.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, ps.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun ps.cont}"),
            ...)
}

#' @exportS3Method NULL
x2base.Match <- function(x, ...) {
  #Process Match
  if (is_not_null(x) && !is.list(x)) {
    arg::err("the supplied {.cls Match} object contains no valid matches")
  }
  
  #Process data and get imp
  .d <- .x2base_data(...)
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
  .reject_args("focal", "{.cls Match} objects", ...)
  
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
  
  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.cls Match} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- .attr(weights, "method")
  
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
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun Match}"),
            ...,
            .call = NULL)
}

#' @exportS3Method NULL
x2base.formula <- function(x, ...) {
  x2base.data.frame(x = x, ...)
}

#' @exportS3Method NULL
x2base.data.frame <- function(x, ...) {
  #Process data.frame
  
  #Process data and get imp
  .d <- .x2base_data(...)
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
    arg::err("{.arg covs} data frame must be specified")
  }
  
  if (rlang::is_formula(covs)) {
    covs <- get_covs_from_formula(covs, data = data)
    # if (is_null(covs)) {
    # }
  }
  
  if (is_null(.attr(covs, "co.names"))) {
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
    .specified_method <- arg::match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    
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
        arg::err("subclassification cannot be specified along with other methods")
      }
      
      if (specified["match.strata"]) {
        arg::err("only weights can be specified with multiple methods")
      }
      
      if (specified["weights"]) {
        method <- .specified_method
        .using <- "weights"
      }
      else {
        arg::wrn("multiple methods were specified, but no weights were provided. Providing unadjusted data only")
        method <- "matching"
      }
    }
  }
  
  if (is_not_null(.using)) {
    if (is_not_null(.specified_args) && is_not_null(.ignoring)) {
      if (is_not_null(.assuming)) {
        arg::msg("{.arg {(.specified_args)}} {?is/are} specified. Assuming {.val {(.assuming)}} and using {.arg {(.using)}} and ignoring {.arg {(.ignoring)}}")
      }
      else {
        arg::msg("{.arg {(.specified_args)}} {?is/are} specified. Using {.arg {(.using)}} and ignoring {.arg {(.ignoring)}}")
      }
    }
    else if (is_not_null(.specified_method) && is_not_null(.not_present) && is_not_null(.assuming)) {
      arg::msg('{.code method = "{(.specified_method)}"} is specified, but {.arg {.not_present}} {?was/were} not supplied. Assuming {.val {(.assuming)}} and using {.arg {(.using)}} instead')
    }
  }
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, covs))
  
  #Process focal
  #A censoring indicator has no treatment group to be focal, and no estimand.
  if (get.treat.type(treat) == "censoring") {
    .reject_cens_args(...)
    estimand <- focal <- NULL
  }
  else {
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
    method <- .attr(weights, "method")
  }
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  discarded <- ...get("discarded")
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl")),
            ...)
}

#' @exportS3Method NULL
x2base.CBPS <- function(x, ...) {
  #Process CBPS
  
  #Process data and get imp
  c.data <- x[["data"]]
  .d <- .x2base_data(..., .datalist = list(c.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
  #Process treat
  treat <- get_treat_from_formula(x[["formula"]], c.data) |>
    process_treat(data, c.data)
  
  #Process covs
  covs <- get_covs_from_formula(x[["formula"]], c.data)
  
  #Process estimand
  estimand <- ...get("estimand")
  .reject_args("estimand", "{.cls CBPS} objects", ...)

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
  .reject_args("focal", "{.cls CBPS} objects", ...)
  
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
  
  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.cls CBPS} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, c.data), 
                             use.weights = ...get("use.weights"))
  method <- .attr(weights, "method")
  
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
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun CBPS}"),
            ...,
            .call = x[["call"]])
}

#' @exportS3Method NULL
x2base.ebalance <- function(x, ...) {
  #Process ebalance
  
  #Process data and get imp
  .d <- .x2base_data(...)
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
  .reject_args("focal", "{.cls ebalance} objects", ...)
  
  focal <- switch(toupper(estimand), 
                  "ATT" = treat_vals(treat)[treat_names(treat)["treated"]], 
                  "ATC" = treat_vals(treat)[treat_names(treat)["control"]], 
                  NULL)
  
  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.cls ebalance} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun ebalance}"),
            ...)
}

#' @exportS3Method NULL
x2base.optmatch <- function(x, ...) {
  #Process optmatch
  if (all(is.na(x))) {
    arg::err("the supplied {.cls optmatch} object contains no valid matches")
  }
  
  #Process data and get imp
  .d <- .x2base_data(...)
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
  #Process treat
  t.c <- .use_tc_fd(...get("formula"), data = data, covs = ...get("covs"),
                    treat = ...get("treat", as.numeric(.attr(x, "contrast.group"))))
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
  
  #Process focal
  focal <- ...get("focal")
  .reject_args(c("subclass", "focal"), "{.cls optmatch} objects", ...)
  
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
  
  #Reject arguments that do not apply
  .reject_args("match.strata", "{.cls optmatch} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = sprintf("{.fun %s}", deparse1(.attr(x, "call")[[1L]]))),
            ...,
            .call = .attr(x, "call"))
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
    arg::err("the supplied {.cls cem.match} object contains no valid matches")
  }
  
  #Process data and get imp
  .d <- .x2base_data(...)
  data <- .d[["data"]]
  imp <- .d[["imp"]]

  if (is_null(data)) {
    arg::err("an argument to {.arg data} must be specified with {.cls cem.match} objects")
  }
  
  if (is_null(imp) && inherits(x, "cem.match.list") &&
      sum(vapply(x, is_, logical(1L), "cem.match")) != 1L) {
    arg::err("an argument to {.arg imp} must be specified or the argument to {.arg data} must be a {.cls mids} object")
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
  
  #Reject arguments that do not apply
  .reject_args("subclass", "{.cls cem.match} objects", ...)
  
  #Process focal
  focal <- x[["baseline.group"]]
  
  #Process pairwise
  if (get.treat.type(treat) == "binary" && is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
    attr(treat, "treat.type") <- "multinomial"
  }
  
  #Reject arguments that do not apply
  .reject_args("match.strata", "{.cls cem.match} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun cem}"),
            ...)
}

#' @exportS3Method NULL
x2base.weightit <- function(x, ...) {
  #Process weightit
  
  #Process data and get imp
  d.e.in.w <- vapply(c("covs", "exact", "by", "moderator"), function(z) is_not_null(x[[z]]), logical(1L))
  weightit.data <- if (any(d.e.in.w)) do.call("data.frame", unname(x[c("covs", "exact", "by", "moderator")[d.e.in.w]]))
  
  .d <- .x2base_data(..., .datalist = list(weightit.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
  #A censoring model's propensity score is P(C = 1 | X); balance on it is as
  #informative here as it is for a binary treatment.
  distance <- process_distance(...get("distance"), datalist = list(data, weightit.data),
                               obj.distance = if (get.treat.type(treat) %in% c("binary", "censoring")) x[["ps"]],
                               obj.distance.name = "prop.score")

  #Process focal
  focal <- x[["focal"]]

  #Process pairwise
  if (get.treat.type(treat) == "binary" && is_null(focal) && isFALSE(...get("pairwise", TRUE))) {
    attr(treat, "treat.type") <- "multinomial"
  }

  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.cls weightit} objects", ...)

  if (get.treat.type(treat) == "censoring") {
    .reject_cens_args(...)
    focal <- NULL
    estimand <- NULL
  }
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, weightit.data))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["s.weights"]]),
                                 data, weightit.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, weightit.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  discarded <- x[["discarded"]]
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun weightit}"),
            ...,
            .call = x[["call"]])
}

#' @exportS3Method NULL
x2base.designmatch <- function(x, ...) {
  #Process designmatch
  if (all(c("id_1", "id_2") %in% names(x))) {
    arg::err("balance cannot currently be checked on a nonbipartite match")
  }
  
  #Process data and get imp
  .d <- .x2base_data(...)
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
  #Process treat
  t.c <- .use_tc_fd(...get("formula"), data, ...get("treat"), ...get("covs"))
  treat <- process_treat(t.c[["treat"]], data)
  if (is.unsorted(rev(treat))) {
    arg::wrn("{.pkg designmatch} requires the input data to be sorted by treatment; the data supplied to {.fn bal.tab} was not, indicating a possible coding error")
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
  .reject_args("focal", "{.pkg designmatch} objects", ...)
  
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
  
  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.pkg designmatch} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "the matching function in {.pkg designmatch}"),
            ...,
            .call = NULL)
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
  
  .d <- .x2base_data(..., .datalist = list(m.data), .imp = m.data[[".imp"]])
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
  .reject_args("focal", "{.cls mimids} objects", ...)
  
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
  
  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.cls mimids} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, m.data))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", unlist(grab(models, "s.weights"))),
                                 data, m.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, m.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(table(imp)))
  
  #Process discarded
  discarded <- unlist(grab(models, "discarded"))
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun matchthem}"),
            ...,
            .call = NULL)
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
  
  .d <- .x2base_data(..., .datalist = list(w.data), .imp = w.data[[".imp"]])
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
  
  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.cls wimids} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, w.data))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", unlist(grab(models, "s.weights"))),
                                 data, w.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, w.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(table(imp)))
  
  #Process discarded
  discarded <- unlist(grab(models, "discarded"))
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun weightthem}"),
            ...,
            .call = NULL)
}

#' @exportS3Method NULL
x2base.sbwcau <- function(x, ...) {
  #Process sbwcau
  
  #Process data and get imp
  sbw.data <- x[["dat_weights"]][names(x[["dat_weights"]]) != "weights"]
  .d <- .x2base_data(..., .datalist = list(sbw.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
  .reject_args("focal", "{.cls sbwcau} objects", ...)
  
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
  
  #Reject arguments that do not apply
  .reject_args(c("subclass", "match.strata"), "{.cls sbwcau} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat, covs, method, addl.data = list(data, sbw.data))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data, sbw.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, sbw.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl"),
                 original.call.to = "{.fun sbw}"),
            ...)
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
  
  rule1 <- .process_stop_method(stop.method, available.stop.methods)
  
  s <- available.stop.methods[match(tolower(rule1), tolower(available.stop.methods))]
  
  #Process data and get imp
  ps.data <- x[["psList"]][[1L]][["data"]]
  .d <- .x2base_data(..., .datalist = list(ps.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
        arg::err("{.arg covs.list} must be a list of covariates for which balance is to be assessed at each time point")
      }
      
      if (!all_apply(covs.list, is_mat_like)) {
        arg::err("each item in {.arg covs.list} must be a data frame")
      }
      
      if (length(covs.list) != length(x[["psList"]])) {
        arg::err("{.arg covs.list} must have as many entries as time points in the call to {.fn iptw}")
      }
    }
    
    formula.list <- ...get("formula.list")
    if (is_not_null(formula.list)) {
      if (!is.list(formula.list) || !all_apply(formula.list, rlang::is_formula)) {
        arg::err("{.arg formula.list} must be a list of formulas identifying the covariates for which balance is to be assessed at each time point")
      }
      
      if (length(formula.list) != length(x[["psList"]])) {
        arg::err("{.arg formula.list} must have as many entries as time points in the call to {.fn iptw}")
      }
    }
    
    covs.list <- lapply(seq_along(x[["psList"]]), function(i) {
      .use_tc_fd(formula = formula.list[[i]],
                 data = data %or% x[["psList"]][[i]][["data"]],
                 covs = covs.list[[i]],
                 needs.treat = FALSE)[["covs"]]
    })
  }
  #See the note in `x2base.ps()`: `feature_names` must be non-empty for this
  #branch, or `all()` of an empty set sends a model object that does not record
  #them into `reformulate(NULL)`.
  else if (all_apply(x[["psList"]], function(z) {
    is_not_null(z[["gbm.obj"]][["feature_names"]]) &&
      all(z[["gbm.obj"]][["feature_names"]] %in% names(z[["data"]]))
  })) {
    covs.list <- lapply(x[["psList"]], function(z) {
      reformulate(z[["gbm.obj"]][["feature_names"]]) |>
        get_covs_from_formula(data = z[["data"]])
    })
  }
  else {
    arg::err('when {.code version = "xgboost"} in the call to {.fn iptw} and any variables are categorical, {.arg formula.list} or {.arg covs.list} must be supplied')
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
  distance.list <- process_distance.list(...get("distance.list", ...get("distance")),
                                         datalist = list(data, ps.data),
                                         covs.list = covs.list,
                                         obj.distance = lapply(x[["psList"]], function(z) z[["ps"]][, s, drop = FALSE]),
                                         obj.distance.name = {
                                           if (length(s) > 1L) paste.("prop.score", substr(s, 1L, nchar(s) - 4L))
                                           else "prop.score"
                                         })
  
  #Reject arguments that do not apply
  .reject_args(c("focal", "subclass", "match.strata"), "{.cls iptw} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat.list[[1L]], covs.list[[1L]],
                             method, addl.data = list(data, ps.data), 
                             stop.method = s)
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["psList"]][[1L]][["sampw"]]),
                                 data, ps.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, ps.data)
  .cluster_check(cluster, treat.list)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(lengths(treat.list)))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("weights"),
                 lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                 original.call.to = "{.fun iptw}"),
            ...,
            .call = NULL,
            .msm = TRUE)
}

#' @exportS3Method NULL
x2base.data.frame.list <- function(x, ...) {
  #Process data and get imp
  .d <- .x2base_data(...)
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
  #Process treat.list
  treat.list <- process_treat.list(...get("treat.list"), data)
  
  #Process covs.list
  covs.list <- x
  if (is_null(covs.list)) {
    arg::err("{.arg covs.list} must be specified")
  }
  
  if (!is.list(covs.list) || is.data.frame(covs.list)) {
    arg::err("{.arg covs.list} must be a list of covariates for which balance is to be assessed at each time point")
  }
  
  if (!all_apply(covs.list, is_mat_like)) {
    arg::err("each item in {.arg covs.list} must be a data frame")
  }
  
  if (any_apply(covs.list, function(z) is_null(.attr(z, "co.names")))) {
    covs.list <- lapply(covs.list, function(z) get_covs_from_formula(data = z))
  }
  
  if (length(treat.list) != length(covs.list)) {
    arg::err("{.arg treat.list} must be a list of treatment statuses at each time point")
  }
  
  #Get estimand
  estimand <- "ATE"
  
  #Get method
  specified <- setNames(rep.int(FALSE, 1L), "weights")
  if (is_not_null(...get("weights"))) {
    if (!is_(...get("weights"), c("character", "numeric", "data.frame", "list"))) {
      arg::err("the argument to {.arg weights} must be a vector, list, or data frame of weights or the (quoted) names of variables in {.arg data} that contain weights")
    }
    specified["weights"] <- TRUE
  }
  
  method <- ...get("method")
  if (is_null(method)) {
    method <- if (specified["weights"]) "weighting" else "matching"
  }
  else if (length(method) == 1L) {
    specified.method <- arg::match_arg(method, c("weighting", "matching", "subclassification"))
    if (specified.method == "weighting") {
      method <- if (specified["weights"]) "weighting" else "matching"
    }
    else if (specified["weights"]) {
      arg::wrn("only weighting is allowed with multiple treatment time points. Assuming weighting instead")
      method <- "weighting"
    }
    else {
      method <- "matching"
    }
  }
  else {
    specified.method <- arg::match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    if (any(specified.method == "subclassification") || specified["subclass"] || specified["match.strata"]) {
      arg::wrn("only weighting is allowed with multiple treatment time points. Assuming weighting instead")
      method <- "weighting"
    }
    else if (specified["weights"]) {
      arg::wrn("only weighting is allowed with multiple treatment time points. Assuming weighting instead")
      method <- "weighting"
    }
    # else if (!specified["weights"]) {
    #   #Should never happen
    #   arg::wrn("multiple methods were specified, but no weights were provided. Providing unadjusted data only")
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
  
  #Reject arguments that do not apply
  .reject_args(c("focal", "subclass", "match.strata"), "longitudinal treatments", ...)
  
  #Process weights
  weights <- ...get("weights")
  if (is_not_null(weights)) {
    weights <- process_weights(NULL, list(...), treat.list[[1L]], covs.list[[1L]],
                               method, addl.data = list(data))
    method <- .attr(weights, "method")
  }
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights"), data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data)
  .cluster_check(cluster, treat.list)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(lengths(treat.list)))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("weights"),
                 lists = c("covs.list", "treat.list", "addl.list", "distance.list")),
            ...,
            .call = NULL,
            .msm = TRUE)
}

#' @exportS3Method NULL
x2base.formula.list <- function(x, ...) {
  
  treat.list <- covs.list <- make_list(length(x))
  
  for (i in seq_along(x)) {
    treat.list[[i]] <- get_treat_from_formula(x[[i]], data = ...get("data"))
    covs.list[[i]] <- get_covs_from_formula(x[[i]], data = ...get("data"))
    names(treat.list)[i] <- .attr(treat.list[[i]], "treat.name")
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
  .d <- .x2base_data(...)
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
  #Process treat.list
  treat.list <- process_treat.list(lapply(times, function(z) x[["treat.hist"]][ID, z]), 
                                   data, cbmsm.data)
  
  #Process covs.list
  covs.list <- make_list(times)
  for (i in seq_along(times)) {
    ti <- times[i]
    cov_i <- get_covs_from_formula(x[["formula"]], data = x[["data"]][x[["time"]] == ti, , drop = FALSE])
    for (co in seq_along(.attr(cov_i, "co.names"))) {
      attr(cov_i, "co.names")[[co]][["component"]][.attr(cov_i, "co.names")[[co]][["type"]] == "base"] <-
        paste0(.attr(cov_i, "co.names")[[co]][["component"]][.attr(cov_i, "co.names")[[co]][["type"]] == "base"], "_T", ti)
    }
    names(attr(cov_i, "co.names")) <- vapply(.attr(cov_i, "co.names"), function(z) paste(z[["component"]], collapse = ""), character(1L))
    colnames(cov_i) <- names(.attr(cov_i, "co.names"))
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
  
  #Reject arguments that do not apply
  .reject_args(c("focal", "subclass", "match.strata"), "{.cls CBMSM} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat.list[[1L]], covs.list[[1L]], method,
                             addl.data = list(data, cbmsm.data))
  method <- .attr(weights, "method")
  
  #Reject arguments that do not apply
  .reject_args("s.weights", "{.cls CBMSM} objects", ...)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, cbmsm.data)
  .cluster_check(cluster, treat.list)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(lengths(treat.list)))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("weights"),
                 lists = c("covs.list", "treat.list", "addl.list", "distance.list"),
                 original.call.to = "{.fun CBMSM}"),
            ...,
            .call = x[["call"]],
            .msm = TRUE)
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
  
  .d <- .x2base_data(..., .datalist = list(weightitMSM.data, weightitMSM.data2))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
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
  distance.list <- process_distance.list(...get("distance.list", ...get("distance")),
                                         datalist = list(data, weightitMSM.data, weightitMSM.data2),
                                         covs.list = covs.list, obj.distance = x[["ps.list"]],
                                         obj.distance.name = "prop.score")
  
  #Reject arguments that do not apply
  .reject_args(c("focal", "subclass", "match.strata"), "{.cls weightitMSM} objects", ...)
  
  #Process weights
  weights <- process_weights(x, list(...), treat.list[[1L]], covs.list[[1L]], method, 
                             addl.data = list(data, weightitMSM.data, weightitMSM.data2))
  method <- .attr(weights, "method")
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", x[["s.weights"]]),
                                 data, weightitMSM.data, weightitMSM.data2)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, weightitMSM.data, weightitMSM.data2)
  .cluster_check(cluster, treat.list)
  
  #Process subset
  subset <- process_subset(...get("subset"), min(lengths(treat.list)))
  
  #Process discarded
  
  #Process output
  .finish_X(list(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("weights"),
                 lists = c("treat.list", "covs.list", "addl.list", "distance.list"),
                 original.call.to = "{.fun weightitMSM}"),
            ...,
            .call = x[["call"]],
            .msm = TRUE)
}

#' @exportS3Method NULL
x2base.default <- function(x, ...) {
  
  if (!is.list(x)) {
    arg::err("the input object must be an appropriate list, data frame, formula, or the output of one of the supported packages")
  }
  
  if (...length() > 0L && (is_null(...names()) || !all(nzchar(...names())))) {
    arg::err("all arguments to {.arg ...} must be named")
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
            #`names(x)` is lowercased below before these are matched, so the
            #aliases must be lowercase too.
            focal = list(name = c("focal", "treatatt"),
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
                                  .attr(P[["distance"]], "name") %or% "distance")
    }
    else {
      P[["distance"]] <- as.data.frame(P[["distance"]])
    }
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
    estimand.name <- .attr(P[["estimand"]], "name")
    
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

.x2base_default_msm <- function(obj, ...) {
  #Process data and get imp
  o.data <- obj[["data"]]
  if (is_null(o.data) && is_not_null(obj[["model"]]) && utils::hasName(obj[["model"]], "data")) {
    o.data <- obj[["model"]][["data"]]
  }
  if (inherits(o.data, "mids")) {
    o.data <- .mids_complete(o.data)
  }

  .d <- .x2base_data(..., .datalist = list(o.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
  #Process treat.list
  treat.list <- process_treat.list(...get("treat.list"), data)
  
  treat.list <- ...get("treat.list")
  if (is_null(treat.list)) {
    formula.list <- ...get("formula.list")
    if (is_not_null(formula.list) && all_apply(formula.list, rlang::is_formula, lhs = TRUE)) {
      treat.list <- make_list(length(formula.list))
      
      for (i in seq_along(formula.list)) {
        treat.list[[i]] <- get_treat_from_formula(formula.list[[i]], data = data %or% o.data)
        names(treat.list)[i] <- .attr(treat.list[[i]], "treat.name")
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
          treat.list[[i]] <- get_treat_from_formula(formula.list.obj[[i]], data = data %or% o.data)
          names(treat.list)[i] <- .attr(treat.list[[i]], "treat.name")
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
        covs.list[[i]] <- get_covs_from_formula(formula.list[[i]], data = data %or% o.data)
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
          covs.list[[i]] <- get_covs_from_formula(formula.list.obj[[i]], data = data %or% o.data)
        }
      }
    }
  }
  
  if (is_null(covs.list)) {
    arg::err("{.arg covs.list} must be specified")
  }
  
  if (!is.list(covs.list) || is.data.frame(covs.list)) {
    arg::err("{.arg covs.list} must be a list of covariates for which balance is to be assessed at each time point")
  }
  
  if (!all_apply(covs.list, is_mat_like)) {
    arg::err("each item in {.arg covs.list} must be a data frame")
  }
  
  if (any_apply(covs.list, function(z) is_null(.attr(z, "co.names")))) {
    covs.list <- lapply(covs.list, function(z) get_covs_from_formula(data = z))
  }
  
  if (length(treat.list) != length(covs.list)) {
    arg::err("{.arg treat.list} must be a list of treatment statuses at each time point")
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
    specified.method <- arg::match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    
    method <- if (specified["weights"]) "weighting" else "matching"
    
    if (any(specified)) {
      .using <- names(specified)[specified][1L]
    }
    
    if (any(specified.method == "subclassification")) {
      if (specified["weights"]) {
        arg::wrn("only weighting is allowed with multiple treatment time points. Assuming weighting instead")
      }
      else {
        arg::wrn("only weighting is allowed with multiple treatment time points. Providing unadjusted data only")
      }
    }
  }
  
  #Process addl.list 
  addl.list <- process_addl.list(...get("addl.list", ...get("addl")),
                                 datalist = list(data, o.data),
                                 covs.list = covs.list)
  
  #Process distance
  #A distance supplied in the input object is stored in `obj[["distance"]]`,
  #whatever name it was given (`distance.list` is among the accepted aliases), so
  #that is the component to fall back on.
  distance.list <- process_distance.list(...get("distance.list", ...get("distance", obj[["distance"]])),
                                         datalist = list(data, o.data),
                                         covs.list = covs.list)
  
  #Reject arguments that do not apply
  .reject_args(c("focal", "subclass", "match.strata"), "longitudinal treatments", ...)
  
  #Process weights
  if ("weights" %in% .using) {
    weights <- process_weights(obj, list(...), treat.list[[1L]], covs.list[[1L]],
                               method, addl.data = list(data, o.data))
    method <- .attr(weights, "method")
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
  
  #Process output
  .finish_X(list(vectors = c("subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("weights"),
                 lists = c("covs.list", "treat.list", "addl.list", "distance.list")),
            ...,
            .msm = TRUE)
}

.x2base_default_point <- function(obj, ...) {
  #Process data and get imp
  o.data <- obj[["data"]]
  if (is_null(o.data) && is_not_null(obj[["model"]]) && utils::hasName(obj[["model"]], "data")) {
    o.data <- obj[["model"]][["data"]]
  }
  if (inherits(o.data, "mids")) {
    o.data <- .mids_complete(o.data)
  }

  .d <- .x2base_data(..., .datalist = list(o.data))
  data <- .d[["data"]]
  imp <- .d[["imp"]]
  
  #Process treat
  treat <- ...get("treat")
  if (is_null(treat)) {
    formula <- ...get("formula")
    if (is_not_null(formula) && rlang::is_formula(formula, lhs = TRUE)) {
      treat <- get_treat_from_formula(formula, data = data %or% o.data)
    }
    else if (is_not_null(obj[["treat"]])) {
      treat <- obj[["treat"]]
    }
    else {
      formula.obj <- obj[["formula"]]
      
      if (is_not_null(formula.obj) && rlang::is_formula(formula.obj, lhs = TRUE)) {
        treat <- get_treat_from_formula(formula.obj, data = data %or% o.data)
      }
    }
  }
  treat <- process_treat(treat, data, o.data)
  
  #Process covs
  covs <- ...get("covs")
  if (is_null(covs)) {
    formula <- ...get("formula")
    if (is_not_null(formula) && rlang::is_formula(formula)) {
      covs <- get_covs_from_formula(formula, data = data %or% o.data)
    }
    else if (is_not_null(obj[["covs"]])) {
      covs <- obj[["covs"]]
    }
    else {
      formula.obj <- obj[["formula"]]
      
      if (is_not_null(formula.obj) && rlang::is_formula(formula.obj)) {
        covs <- get_covs_from_formula(formula.obj, data = data %or% o.data)
      }
    }
  }
  
  if (is_null(covs)) {
    arg::err("{.arg covs} data frame must be specified")
  }
  
  if (is_null(.attr(covs, "co.names"))) {
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
    .specified_method <- arg::match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    
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
        arg::err("subclassification cannot be specified along with other methods")
      }
      
      if (specified["match.strata"]) {
        arg::err("only weights can be specified with multiple methods")
      }
      
      if (specified["weights"]) {
        method <- .specified_method
        .using <- "weights"
      }
      else {
        arg::wrn("multiple methods were specified, but no weights were provided. Providing unadjusted data only")
        method <- "matching"
      }
    }
  }
  
  if (is_not_null(.using)) {
    if (is_not_null(.specified_args) && is_not_null(.ignoring)) {
      if (is_not_null(.assuming)) {
        arg::msg("{.arg {(.specified_args)}} {?is/are} specified. Assuming {.val {(.assuming)}} and using {.arg {(.using)}} and ignoring {.arg {(.ignoring)}}")
      }
      else {
        arg::msg("{.arg {(.specified_args)}} {?is/are} specified. Using {.arg {(.using)}} and ignoring {.arg {(.ignoring)}}")
      }
    }
    else if (is_not_null(.specified_method) && is_not_null(.not_present) && is_not_null(.assuming)) {
      arg::msg('{.code method = "{(.specified_method)}"} is specified, but {.arg {.not_present}} {?was/were} not supplied. Assuming {.val {(.assuming)}} and using {.arg {(.using)}} instead')
    }
  }
  
  #Process addl 
  addl <- process_addl(...get("addl"), datalist = list(data, o.data, covs))
  
  #Process distance
  distance <- process_distance(...get("distance"), datalist = list(data, o.data, covs),
                               obj.distance = obj[["distance"]], 
                               obj.distance.name = .attr(obj[["distance"]], "name"))
  
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
    method <- .attr(weights, "method")
  }
  
  #Process s.weights
  s.weights <- process_s.weights(...get("s.weights", obj[["s.weights"]]), data, o.data)
  
  #Process cluster
  cluster <- process_cluster(...get("cluster"), data, o.data)
  
  #Process subset
  subset <- process_subset(...get("subset"), length(treat))
  
  #Process discarded
  discarded <- ...get("discarded", obj[["discarded"]])
  
  #Process output
  .finish_X(list(vectors = c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset", "discarded"),
                 data.frames = c("covs", "weights", "distance", "addl")),
            ...)
}
