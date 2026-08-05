#' Using `bal.tab()` with a Censoring Indicator
#' @name class-bal.tab.cens
#'
#' @description
#' When [bal.tab()] is given a censoring indicator marked with [.cens()] rather than a treatment, it asks a different question. A censoring model has no second treatment group to compare against; its target is the full at-risk sample, and what matters is whether the units still under observation resemble that sample once weighted. This page outlines the output in this case.
#'
#' Balance is therefore assessed between two samples built from the same units:
#'
#' * **Uncensored**: the units with `C == 0`, carrying the censoring weights.
#' * **Full**: every at-risk unit, i.e. every unit with a non-missing `C`, carrying a weight of 1.
#'
#' Setting `un = TRUE` adds the same comparison with the uncensored units unweighted, which is the balance the weights were estimated to remove.
#'
#' The two samples are compared using the same statistics available for binary treatments (`mean.diffs`, `variance.ratios`, `ks.statistics`, and `ovl.coefficients`), because internally they are stacked into a binary comparison: group `0` is the uncensored sample and group `1` the full sample. So in the balance table, `M.0` and `SD.0` describe the uncensored sample, `M.1` and `SD.1` the full sample, and `Diff` is the difference between the full sample and the uncensored one. This is the same arrangement described under "Assessing balance" in `WeightIt::.cens()`, so the values agree with the manual computation documented there.
#'
#' @section Allowable arguments:
#'
#' Every argument to `bal.tab()` applies as it does with a binary treatment, with the following exceptions.
#'
#' \describe{
#'     \item{`s.d.denom`}{allowable values are `"full"` (the default), which uses the standard deviation of the covariate in the full at-risk sample; `"uncensored"`, which uses that of the uncensored sample; and `"pooled"`, `"all"`, `"weighted"`, and `"hedges"`, which behave as they do for a binary treatment. Because a censoring model's target is fixed by its design rather than inferred, no note is printed about the value used.}
#'     \item{`estimand` and `focal`}{do not apply and are ignored with a warning. The target is the full at-risk sample.}
#'     \item{`subclass` and `match.strata`}{do not apply, as subclassification is not a way of estimating censoring weights.}
#' }
#'
#' `cluster` and `imp` apply as usual, and produce a `bal.tab.cluster` or `bal.tab.imp` object whose per-cluster or per-imputation components are the censoring balance tables described here.
#'
#' @section Output:
#' The output is a `bal.tab.cens` object, which inherits from `bal.tab.bin` and `bal.tab`, and has the same elements as an ordinary binary-treatment `bal.tab` object. Its `Observations` component differs: it has a single `Total` column, because the target sample is the same in every row, and the rows
#'
#' * `Full`: the number (or effective number) of at-risk units,
#' * `Uncensored`: the number of units still under observation,
#' * `Adjusted` (or one row per set of weights): their effective number once weighted,
#' * `Censored`: the number of units censored, omitted when there are none.
#'
#' `Uncensored` and `Censored` sum to `Full`.
#'
#' @seealso
#' * [bal.tab()]
#' * [.cens()] for marking an indicator as censoring
#' * [`class-bal.tab.cluster`] and [`class-bal.tab.imp`] for the segmented cases
#'
NULL

base.bal.tab.cens <- function(X, ...) {
  C <- X[["treat"]]

  #A unit with a missing indicator is not at risk -- in a longitudinal model it was
  #censored at an earlier time point -- so it is in neither sample.
  at.risk <- which(!is.na(C))
  uncensored <- which(!is.na(C) & C == 0)

  if (is_null(uncensored)) {
    arg::err("every unit is censored, so there is no sample left to compare against the full one")
  }

  n.uncensored <- length(uncensored)

  obs <- samplesize_cens(C, weights = X[["weights"]], s.weights = X[["s.weights"]])

  #The two samples become one stacked data set: the uncensored units, then every
  #at-risk unit. `subset_X()` does the stacking, since it takes indices and nothing
  #stops an index from repeating; that way every per-unit slot is stacked the same way,
  #including ones added later.
  X <- subset_X(X, c(uncensored, at.risk))

  X[["treat"]] <- .cens_pseudo_treat(n.uncensored, length(at.risk))

  #The full sample is the target, so it is never reweighted.
  if (is_not_null(X[["weights"]])) {
    X[["weights"]][-seq_len(n.uncensored), ] <- 1
  }

  #The target of a censoring model is fixed by its design, so the denominator is settled
  #here rather than inferred from the weights further down. Handing it on as if the user
  #had supplied it is also what keeps `.get_s.d.denom()` from printing a note about a
  #choice nobody made.
  s.d.denom <- .cens_s.d.denom(X[["s.d.denom"]])

  X[["s.d.denom"]] <- switch(s.d.denom,
                             "full" = "treated",
                             "uncensored" = "control",
                             s.d.denom)

  X[c("estimand", "focal")] <- NULL

  out <- base.bal.tab.base(X, type = "bin", .obs = obs, ...)

  set_class(out, c("bal.tab.cens", class(out)))
}

#The stacked pseudo-treatment. Built as a processed binary treatment so that every
#binary statistic, and the whole binary column layout, apply without change. `0` is the
#uncensored sample and `1` the full one, matching the stacking documented in
#`WeightIt::.cens()` and cobalt's own ATT convention, in which the target is the treated
#group and the reweighted group is the control.
.cens_pseudo_treat <- function(n.uncensored, n.full) {
  treat <- process_treat(rep(c(0, 1), times = c(n.uncensored, n.full)))

  treat_names(treat) <- setNames(c("Uncensored", "Full"), c("control", "treated"))
  treat_vals(treat) <- setNames(c("0", "1"), treat_names(treat))

  treat
}

#`s.d.denom` for a censoring model, in the vocabulary of the two samples rather than of
#treated and control, defaulting to the target. Stays in that vocabulary so that it can
#be resolved once by a wrapper and handed down: a value already resolved must survive
#being resolved again.
.cens_s.d.denom <- function(s.d.denom = NULL) {
  if (is_null(s.d.denom)) {
    return("full")
  }

  arg::match_arg(s.d.denom,
                 c("full", "uncensored", "pooled", "all", "weighted", "hedges"))
}

#Sample sizes for a censoring model. One column, because the target sample is the same
#in every row, and the rows describe the samples rather than the adjustment: the
#at-risk sample, the units still under observation before and after weighting, and how
#many were censored.
samplesize_cens <- function(C, weights = NULL, s.weights = NULL) {
  if (is_null(s.weights)) {
    s.weights <- rep_with(1, C)
  }

  at.risk <- !is.na(C)
  uncensored <- at.risk & C == 0

  adj.rows <- {
    if (is_null(weights)) character()
    else if (NCOL(weights) == 1L) "Adjusted"
    else names(weights)
  }

  nn <- make_df("Total", c("Full", "Uncensored", adj.rows, "Censored"))

  nn["Full", ] <- ESS(s.weights[at.risk])
  nn["Uncensored", ] <- ESS(s.weights[uncensored])

  for (j in seq_along(adj.rows)) {
    nn[adj.rows[j], ] <- ESS((weights[[j]] * s.weights)[uncensored])
  }

  nn["Censored", ] <- sum(at.risk & C == 1)

  attr(nn, "ss.type") <- c("ss", "ss", rep_with("ess", adj.rows), "ss")

  #A `Censored` row that no unit reached says nothing, as with `Discarded`.
  if (nn["Censored", ] == 0) {
    attr(nn, "ss.type") <- .attr(nn, "ss.type")[rownames(nn) != "Censored"]
    nn <- nn[rownames(nn) != "Censored", , drop = FALSE]
  }

  attr(nn, "tag") <- {
    if (is_null(adj.rows)) "Sample sizes"
    else "Effective sample sizes"
  }

  nn
}
