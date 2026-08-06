#' Using `bal.tab()` with Longitudinal Treatments
#' @name class-bal.tab.msm
#' 
#' @description
#'     When using [bal.tab()] with longitudinal treatments, the output will be different from the case with point treatments, and there are some options that are common across all `bal.tab()` methods for dealing with longitudinal data. This page outlines the outputs and options in this case.
#'     
#'     There are two main components of the output of `bal.tab()` with longitudinal treatments: the time-point-specific balance summary and across-time-points balance summary. The time-point-specific balance summaries are standard point treatment balance summaries at each time point.
#'
#'     The across-time-points balance summary is, for each variable, the greatest imbalance across all time-point-specific balance summaries. If the greatest observed imbalance is tolerable, then all other imbalances for that variable will be tolerable too, so focusing on reducing the greatest imbalance is sufficient for reducing imbalance overall. The balance summary will not be computed if multi-category treatments or multiply imputed data are used, or if the time points are not all of the same kind (see below).
#'
#' @section Mixing treatments and censoring:
#' One entry of the list may be a censoring indicator marked with [.cens()] rather than a treatment, as in `list(A1 ~ x, .cens(C1) ~ x, A2 ~ x)`, which is how a joint treatment-and-censoring model is written for `WeightIt::weightitMSM()`. `bal.tab()` produces one table per entry, of whichever kind that entry is: an ordinary balance table for a treatment, and the censoring balance table described at [`class-bal.tab.cens`] for an indicator.
#'
#' Each entry is assessed among the units still under observation entering it. A unit is under observation until a censoring indicator earlier in the list marks it censored; the risk set is accumulated from the indicators themselves rather than read off the treatments, so it makes no difference whether the data records a treatment for a unit that has already dropped out or leaves it missing. A missing treatment for a unit that is still under observation is an error naming the time point it appeared in. A censoring entry's own comparison uses the risk set as it stands entering it, so it includes the units it is about to remove; those are the full sample it compares against.
#'
#' A censoring balance table and a treatment balance table say different things about different samples, so a list that mixes them gets no balance summary across time points, exactly as a list mixing continuous and binary treatments does not. The sample sizes are still reported for every time point. A list in which every entry is a censoring indicator is not a mixture and is summarized as usual.
#'
#' Each time point's default `s.d.denom` is the one its own kind of model implies -- `"pooled"` for a binary or multi-category treatment, `"all"` for a continuous one, and `"full"` for a censoring indicator -- so that a model gives the same numbers in a list as it would on its own. A value supplied to `s.d.denom` is shared by every time point, which is one reason it is recommended not to set it for longitudinal treatments.
#'
#' Note that `weightitMSM()` returns a single set of weights, the product across all the models, so balance at an early time point is assessed with weights that include later censoring.
#'
#' @section Allowable arguments:
#' 
#' There are two additional arguments for each `bal.tab()` method that can handle longitudinal treatments: `which.time` and `msm.summary`.
#' 
#' \describe{
#'     \item{`which.time`}{This is a display option that does not affect computation. If `.all` (the default), all time points will be displayed. If `.none`, no time points will be displayed. Otherwise, can be a vector of treatment names or indices for which to display balance.}
#'     \item{`msm.summary`}{This is a display option that does not affect computation. If `TRUE`, the balance summary across time points will be displayed. The default is `TRUE`, and if `which.time` is `.none`, it will automatically be set to `TRUE`.}
#' }
#' 
#' @section Output: 
#' The output is a `bal.tab.msm` object, which inherits from `bal.tab`. It has the following elements:
#'         
#' * `Time.Balance`: For each time point, a regular `bal.tab` object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified.
#' * `Balance.Across.Times`: The balance summary across time points. This will include the maximum balance statistic(s) for each covariate across all time points. Absent when the time points are not all of the same kind.
#' * `Observations`: A list with a table of sample sizes or effective sample sizes for each time point, before and after adjustment. Always present, even when there is no balance summary.
#'
#' As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
#'
#' @note The balance tables presented here are not the same as those recommended by Jackson (2016) and computed in his R package, \CRANpkg{confoundr}, as these do not take into account treatment history. The balance statistics presented here should be used with caution and may not reflect balance in an accurate way.
#' 
#' @references 
#' Jackson, J. W. (2016). Diagnostics for Confounding of Time-varying and Other Joint Exposures: *Epidemiology*, 27(6), 859–869. \doi{10.1097/EDE.0000000000000547}
#' 
#' @seealso
#' * [bal.tab()]
#' * [bal.tab.time.list()]
#' * [print.bal.tab()]
#' * [.cens()] and [`class-bal.tab.cens`] for a censoring indicator among the time points
#' * `vignette("longitudinal-treat")` for examples
#'
NULL

base.bal.tab.msm <- function(X,
                             which.time,
                             msm.summary = getOption("cobalt_msm.summary"),
                             ...) {
  A <- list(...)
  
  #Preparations
  
  if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
  
  treat.types <- vapply(X[["treat.list"]], get.treat.type, character(1L))
  
  if (missing(which.time)) {
    which.time <- {
      if (all_the_same(treat.types) && "multinomial" %nin% treat.types && is_null(X[["imp"]])) NA
      else NULL
    }
  }
  
  if (is_null(msm.summary)) {
    msm.summary <- is_not_null(which.time) &&
      (anyNA(which.time) ||
         !(is.character(which.time) || is.numeric(which.time)) ||
         (is.numeric(which.time) && !any(which.time %in% seq_along(X[["treat.list"]]))) ||
         (is.character(which.time) && !any(which.time %in% names(X[["treat.list"]]))))
  }
  
  #Setup output object
  out <- list()

  time.names <- {
    if (length(names(X[["treat.list"]])) == length(X[["treat.list"]])) names(X[["treat.list"]])
    else as.character(seq_along(X[["treat.list"]]))
  }

  #The risk set entering each entry of the list. All `TRUE` unless the list contains a
  #censoring indicator, in which case `subset_X()` returns each time point untouched and
  #nothing below this changes.
  at.risk <- .msm_at_risk(X[["treat.list"]])

  #Get list of bal.tabs for each time period
  out[["Time.Balance"]] <- lapply(seq_along(X[["covs.list"]]), function(ti) {
    #The whole body is inside `tryCatch()` so that an error raised for one time point --
    #by the restriction below, or by a leaf that rejects a shared argument -- says which
    #one it came from, as `base.bal.tab.cluster()` does for a cluster.
    tryCatch({
      X_ti <- X

      X_ti[["covs"]] <- X_ti[["covs.list"]][[ti]]
      X_ti[["treat"]] <- X_ti[["treat.list"]][[ti]]
      X_ti[["addl"]] <- X_ti[["addl.list"]][[ti]]
      X_ti[["distance"]] <- X_ti[["distance.list"]][[ti]]

      X_ti[c("covs.list", "treat.list", "addl.list", "distance.list", "call")] <- NULL

      #A unit censored earlier has no treatment here and nothing left to balance, so
      #this time point is about the units still under observation. A missing treatment
      #among *those* is a different thing -- an unknown rather than an absence -- and
      #cannot be put in either group.
      if (!.is_cens(X_ti[["treat"]]) && anyNA(X_ti[["treat"]][at.risk[[ti]]])) {
        arg::err("missing values must not exist in {.arg treat} for units still under observation")
      }

      X_ti <- subset_X(X_ti, at.risk[[ti]]) |>
        .assign_X_class()

      #A longitudinal treatment targets the ATE, so each time point's default
      #denominator is the ATE's rather than one inferred from that time point's own
      #weights; a censoring model's is instead its target, the full at-risk sample.
      #This reads the treatment type rather than the class of `X_ti` because a time
      #point wrapped in clusters or imputations has the class of the wrapper, and the
      #vocabulary `s.d.denom` is written in belongs to the treatment. A value the user
      #supplied still takes precedence.
      X_ti[["s.d.denom"]] <- X_ti[["s.d.denom"]] %or%
        switch(get.treat.type(X_ti[["treat"]]),
               continuous = "all",
               censoring = "full",
               "pooled")

      do.call("base.bal.tab", c(list(X_ti), A[setdiff(names(A), names(X_ti))]), quote = TRUE)
    },
    error = function(e) {
      arg::err("in time point {.str {time.names[ti]}}: {conditionMessage(e)}")
    })
  })

  names(out[["Time.Balance"]]) <- time.names

  if ((!A[["quick"]] || msm.summary) && is_null(X[["imp"]]) && all_the_same(treat.types) &&
      !any(treat.types == "multinomial")) {
    summ <- .bal.tab_summarize(out[["Time.Balance"]], "Balance.Across.Times",
                               agg.funs = "max",
                               obs.fun = function(cl) grab(cl, "Observations"),
                               include.times = TRUE)

    out[names(summ)] <- summ
  }

  #Sample sizes describe each time point on its own, so they are produced whether or not
  #the time points can be summarized together -- a list mixing kinds of model has no
  #summary but still has sample sizes. `.bal.tab_summarize()` produces the same table
  #when there is a summary; assigning over it keeps the element where it was.
  out[["Observations"]] <- grab(out[["Time.Balance"]], "Observations")

  out[["call"]] <- X[["call"]]
  
  attr(out, "print.options") <- c(.attr(out[["Time.Balance"]][[1L]], "print.options"),
                                  list(which.time = which.time,
                                       msm.summary = msm.summary))
  
  set_class(out, c("bal.tab.msm", "bal.tab"))
}

#Which units are still under observation entering each entry of a longitudinal treatment
#list. Everyone starts at risk, and each censoring indicator removes the units censored
#there from every entry after it -- so a censoring entry's own risk set still contains
#the units it is about to remove, which are half of what that entry compares.
#
#This is `WeightIt::weightitMSM()`'s own rule, and it reads only the indicators. Inferring
#it instead from wherever a later treatment happens to be `NA` would depend on a data
#convention: a user who records post-censoring treatments rather than blanking them would
#silently get censored units counted in later time points.
.msm_at_risk <- function(treat.list) {
  n <- len(treat.list[[1L]])

  #One risk set indexes every time point, so they have to describe the same units.
  if (!all(lengths(treat.list) == n)) {
    arg::err("{.arg treat.list} must have the same number of units at each time point")
  }

  at.risk <- make_list(length(treat.list))

  in.risk.set <- rep.int(TRUE, n)

  for (ti in seq_along(treat.list)) {
    at.risk[[ti]] <- in.risk.set

    if (.is_cens(treat.list[[ti]])) {
      C <- as.numeric(treat.list[[ti]])

      #A missing indicator marks a unit that was not at risk to begin with, so it is not
      #at risk afterwards either.
      in.risk.set <- in.risk.set & !is.na(C) & C == 0
    }
  }

  at.risk
}
