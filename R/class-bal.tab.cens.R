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
#' The two samples are compared using the same statistics available for binary treatments (`mean.diffs`, `variance.ratios`, `ks.statistics`, and `ovl.coefficients`), because internally they are stacked into a binary comparison. The two samples name their own columns rather than taking a binary treatment's positional `0` and `1`, so the balance table has `M.Uncensored` and `SD.Uncensored` for the uncensored sample and `M.Full` and `SD.Full` for the full one, and `Diff` is the difference between the full sample and the uncensored one. That direction is the arrangement described under "Assessing balance" in `WeightIt::.cens()`, so the values agree with the manual computation documented there.
#'
#' [bal.plot()] works too, showing the weighted uncensored sample against the unweighted full sample.
#'
#' @section Allowable arguments:
#'
#' Every argument to `bal.tab()` applies as it does with a binary treatment, with the following exceptions.
#'
#' \describe{
#'     \item{`s.d.denom`}{allowable values are `"full"` (the default), which uses the standard deviation of the covariate in the full at-risk sample; `"uncensored"`, which uses that of the uncensored sample; and `"pooled"`, `"all"`, `"weighted"`, and `"hedges"`, which behave as they do for a binary treatment. Because a censoring model's target is fixed by its design rather than inferred, no note is printed about the value used.}
#'     \item{`estimand` and `focal`}{do not apply and are ignored with a warning. The target is the full at-risk sample.}
#'     \item{`subclass`}{applies: subclassifying is itself a way of solving a censoring problem, in that within each subclass the uncensored units should resemble every at-risk unit in it. See "With subclasses" below.}
#'     \item{`match.strata`}{does not apply. It turns strata into weights before the two samples exist; the same strata say the same thing supplied to `subclass`.}
#' }
#'
#' `cluster` and `imp` apply as usual, and produce a `bal.tab.cluster` or `bal.tab.imp` object whose per-cluster or per-imputation components are the censoring balance tables described here.
#'
#' @section Among longitudinal treatments:
#' A censoring indicator can appear among a list of longitudinal treatments, as in `list(A1 ~ x, .cens(C1) ~ x, A2 ~ x)`, which is how a joint treatment-and-censoring model is written for `WeightIt::weightitMSM()`; `bal.tab()` also accepts such a `weightitMSM` object directly. Each entry of the list gets a table of its own kind, and each is assessed among the units still under observation entering it, so the full sample a censoring indicator is compared against is the risk set at that time point rather than the original cohort, and a treatment after it is assessed only among the units it did not censor. A list that mixes censoring with treatment gets no balance summary across time points. See [`class-bal.tab.msm`].
#'
#' @section With subclasses:
#' Subclassification is an alternative to weighting for solving a censoring problem: within each subclass, the units still under observation should resemble every at-risk unit in that subclass. Supplying `subclass` therefore produces a `bal.tab.cens` object that also inherits from `bal.tab.subclass`, with a balance table for each subclass and a summary across them, as described at [`class-bal.tab.subclass`].
#'
#' The summary across subclasses is subclassification expressed as censoring weights: a unit still under observation in subclass \eqn{k} receives \eqn{n_k / n_{k1}}, where \eqn{n_k} is the number of at-risk units in the subclass and \eqn{n_{k1}} the number of them still under observation, and the full sample is left unweighted. This is the same summary one would get from supplying those weights to a censoring model directly.
#'
#' The sample sizes have one column per subclass, with `Full`, `Uncensored`, and `Censored` rows.
#'
#' @section Output:
#' The output is a `bal.tab.cens` object, which inherits from `bal.tab.bin` and `bal.tab` (or from `bal.tab.subclass` when `subclass` is supplied), and has the same elements as an ordinary binary-treatment `bal.tab` object. Its `Observations` component differs: it has a single `Total` column, because the target sample is the same in every row, and the rows
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
#' * [`class-bal.tab.cluster`], [`class-bal.tab.imp`], [`class-bal.tab.subclass`], and [`class-bal.tab.msm`] for the segmented and longitudinal cases
#' * [`treat-class`] for the attributes that decide how a treatment is compared and what its groups are called
#'
NULL

base.bal.tab.cens <- function(X, ...) {
  s <- .stack_cens_X(X)

  out <- base.bal.tab.base(s[["X"]], type = "bin", .obs = s[["obs"]], ...)

  set_class(out, c("bal.tab.cens", class(out)))
}

base.bal.tab.subclass.cens <- function(X, ...) {
  s <- .stack_cens_X(X)

  X <- s[["X"]]

  #Subclassification is itself the censoring solution: within each subclass the
  #uncensored units should look like every at-risk unit in it. Expressed as weights that
  #is `n_k / n_{k,uncensored}` for the uncensored units and 1 for the full sample, which
  #is exactly what `strata2weights()` produces for an ATT whose focal group is the full
  #sample. So the summary across subclasses needs nothing new either; it just needs to
  #be told which group is the target.
  X[["estimand"]] <- "ATT"

  out <- base.bal.tab.subclass(X, type = "bin", .obs = s[["obs"]], ...)

  set_class(out, c("bal.tab.cens", class(out)))
}

#The two samples of a censoring model, stacked into one `X` with a binary
#pseudo-treatment, together with the sample sizes of the units they were built from.
#Both the plain and the subclassified leaf start here.
.stack_cens_X <- function(X, .count = TRUE) {
  C <- X[["treat"]]

  s <- .cens_stack_index(C)

  #`bal.plot()` wants the stacking but not the counting, which would also warn about
  #subclasses too small to compute a balance statistic in -- irrelevant to a plot.
  obs <- {
    if (.count) samplesize_cens(C, weights = X[["weights"]], s.weights = X[["s.weights"]],
                                subclass = X[["subclass"]])
    else NULL
  }

  #`subset_X()` does the stacking, since it takes indices and nothing stops an index
  #from repeating; that way every per-unit slot is stacked the same way, including ones
  #added later.
  X <- subset_X(X, s[["index"]])

  X[["treat"]] <- .cens_pseudo_treat(s[["n.uncensored"]], s[["n.at.risk"]])

  #The full sample is the target, so it is never reweighted.
  if (is_not_null(X[["weights"]])) {
    X[["weights"]][-seq_len(s[["n.uncensored"]]), ] <- 1
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

  list(X = X, obs = obs)
}

#The row indices that stack a censoring model's two samples: the units still under
#observation, then every at-risk unit. An uncensored unit appearing in both is the point,
#and how many rows each sample contributes comes back with them.
#
#A unit with a missing indicator is not at risk -- in a longitudinal model it was censored
#at an earlier time point -- so it is in neither sample. `.at.risk` says the same thing
#from the other direction: a longitudinal caller has already worked out who was still
#under observation entering this time point, and passes it rather than relying on the
#indicator being blank for everyone else.
.cens_stack_index <- function(C, .at.risk = TRUE) {
  at.risk <- which(.at.risk & !is.na(C))
  uncensored <- which(.at.risk & !is.na(C) & C == 0)

  if (is_null(uncensored)) {
    arg::err("every unit is censored, so there is no sample left to compare against the full one")
  }

  list(index = c(uncensored, at.risk),
       n.uncensored = length(uncensored),
       n.at.risk = length(at.risk))
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

  #Once the two samples exist they are what the balance table's columns are about, so
  #they name those columns: `M.Uncensored` and `M.Full` rather than `M.0` and `M.1`.
  group_labels(treat) <- unname(treat_names(treat))

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
samplesize_cens <- function(C, weights = NULL, s.weights = NULL, subclass = NULL) {
  if (is_null(s.weights)) {
    s.weights <- rep_with(1, C)
  }

  at.risk <- !is.na(C)
  uncensored <- at.risk & C == 0

  #Subclassification is the adjustment, so there is nothing to weight and the interest
  #is in how the two samples divide up. One column per subclass, transposed relative to
  #the unsubclassified table for the same reason the binary one is.
  if (is_not_null(subclass)) {
    subclass <- factor(subclass)
    in.subclass <- !is.na(subclass)

    nn <- make_df(c(levels(subclass), "All"), c("Full", "Uncensored", "Censored"))

    counts <- function(i) {
      c(vapply(levels(subclass),
               function(k) sum(i & in.subclass & subclass == k), numeric(1L)),
        "All" = sum(i))
    }

    nn["Full", ] <- counts(at.risk)
    nn["Uncensored", ] <- counts(uncensored)
    nn["Censored", ] <- counts(at.risk & C == 1)

    small.subclass <- nn["Uncensored", levels(subclass)] <= 1L
    if (any(small.subclass)) {
      arg::wrn("not enough uncensored units in {cli::qty(sum(small.subclass))} subclass{?es} {levels(subclass)[small.subclass]}")
    }

    attr(nn, "tag") <- "Sample sizes by subclass"

    return(nn)
  }

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

  #A `Censored` row that no unit reached says nothing, as with `Discarded`.
  if (nn["Censored", ] == 0) {
    nn <- nn[rownames(nn) != "Censored", , drop = FALSE]
  }

  #No `ss.type`, so no row is starred. That attribute marks which rows of a treatment's
  #table are effective sample sizes, and it earns its keep there because such a table can
  #hold several sets of weights fit by different methods, only some of them weighting.
  #Here the only row that could be an effective sample size is the adjusted one, and the
  #heading says so; the rest are counts of units, which is what their names call them.
  attr(nn, "tag") <- {
    if (is_null(adj.rows)) "Sample sizes"
    else "Effective sample sizes"
  }

  nn
}
