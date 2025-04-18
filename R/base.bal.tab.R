base.bal.tab <- function(X, ...) {
  fun <- switch(attr(X, "X.class"),
                "binary" = base.bal.tab.binary,
                "cont" = base.bal.tab.cont,
                "subclass.binary" = base.bal.tab.subclass.binary,
                "subclass.cont" = base.bal.tab.subclass.cont,
                "cluster" = base.bal.tab.cluster,
                "msm" = base.bal.tab.msm,
                "multi" = base.bal.tab.multi,
                "imp" = base.bal.tab.imp)
  
  fun(X, ...)
}

base.bal.tab.binary <- function(X, ...) {
  base.bal.tab.base(X, type = "bin", ...)
}
base.bal.tab.cont <- function(X, ...) {
  base.bal.tab.base(X, type = "cont", ...)
}

base.bal.tab.base <- function(X,
                              type,
                              int = FALSE,
                              poly = 1,
                              continuous = NULL,
                              binary = NULL,
                              imbalanced.only = getOption("cobalt_imbalanced.only", FALSE),
                              un = getOption("cobalt_un", FALSE),
                              disp = NULL,
                              disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE),
                              disp.call = getOption("cobalt_disp.call", FALSE),
                              abs = FALSE,
                              quick = TRUE,
                              ...) {
  #Preparations
  A <- clear_null(list(...))
  A$subset <- NULL
  
  if (type == "bin" && get.treat.type(X$treat) != "binary") {
    .err("the treatment must be a binary variable")
  }
  
  if (is_null(continuous)) continuous <- getOption("cobalt_continuous", "std")
  if (is_null(binary)) binary <- getOption("cobalt_binary", switch(type, bin = "raw", "std"))
  
  if (is_null(X$weights)) {
    un <- TRUE
    no.adj <- TRUE
  }
  else {
    no.adj <- FALSE
    
    if (type == "bin") check_if_zero_weights(X$weights, X$treat)
    else if (type == "cont") check_if_zero_weights(X$weights)
    
    if (ncol(X$weights) == 1L) names(X$weights) <- "Adj"
  }
  
  if (is_null(X$s.weights)) {
    X$s.weights <- rep.int(1, length(X$treat))
  }
  
  disp <- process_disp(disp, ...)
  
  #Actions
  out <- list()
  
  C <- do.call(".get_C2", c(X, A[setdiff(names(A), names(X))], list(int = int, poly = poly)), quote = TRUE)
  
  co.names <- attr(C, "co.names")
  
  var_types <- attr(C, "var_types")
  
  if (is_not_null(X$s.d.denom.list)) {
    X$s.d.denom <- NULL
  }
  else if (get.treat.type(X$treat) != "continuous" &&
           "mean.diffs" %in% X$stats &&
           ((binary == "std" && any(var_types == "Binary")) ||
            (continuous == "std" && !all(var_types == "Binary")))) {
    X$s.d.denom <- .get_s.d.denom(X$s.d.denom,
                                  estimand = X$estimand,
                                  weights = X$weights, 
                                  subclass = X$subclass,
                                  treat = X$treat,
                                  focal = X$focal)
  }
  else if (get.treat.type(X$treat) == "continuous" &&
           any(c("correlations", "spearman.correlations", "distance.correlations") %in% X$stats) &&
           ((binary == "std" && any(var_types == "Binary")) ||
            (continuous == "std" && !all(var_types == "Binary")))) {
    X$s.d.denom <- .get_s.d.denom.cont(X$s.d.denom,
                                       weights = X$weights,
                                       subclass = X$subclass)
  }
  else {
    X$s.d.denom <- NULL
  }
  
  out[["Balance"]] <- do.call("balance_table",
                              c(list(C, type = type, weights = X$weights, treat = X$treat, 
                                     s.d.denom = X$s.d.denom, s.weights = X$s.weights, 
                                     continuous = continuous, binary = binary, 
                                     thresholds = X$thresholds,
                                     un = un, disp = disp, 
                                     stats = X$stats, abs = abs, 
                                     no.adj = no.adj, quick = quick, 
                                     var_types = var_types,
                                     s.d.denom.list = X$s.d.denom.list),
                                A),
                              quote = TRUE)
  
  #Reassign disp... and ...threshold based on balance table output
  compute <- attr(out[["Balance"]], "compute")
  thresholds <- attr(out[["Balance"]], "thresholds")
  disp <- attr(out[["Balance"]], "disp")
  
  out <- c(out,
           threshold_summary(compute = compute,
                             thresholds = thresholds,
                             no.adj = no.adj,
                             balance.table = out[["Balance"]],
                             weight.names = names(X$weights)))
  
  out[["Observations"]] <- samplesize(treat = X$treat, type = type, weights = X$weights,
                                      s.weights = X$s.weights, method = X$method,
                                      discarded = X$discarded)
  
  out[["call"]] <- X$call
  attr(out, "print.options") <- list(thresholds = thresholds,
                                     imbalanced.only = imbalanced.only,
                                     un = un, 
                                     compute = compute, 
                                     disp = disp,
                                     disp.adj = !no.adj,
                                     disp.bal.tab = disp.bal.tab,
                                     disp.call = disp.call, 
                                     abs = abs,
                                     continuous = continuous,
                                     binary = binary,
                                     quick = quick,
                                     nweights = if (no.adj) 0L else ncol(X$weights),
                                     weight.names = names(X$weights),
                                     treat_names = treat_names(X$treat),
                                     type = type,
                                     co.names = co.names)
  
  set_class(out, c(paste.("bal.tab", type), "bal.tab"))
}
