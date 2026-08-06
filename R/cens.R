#' @title Mark a Censoring Indicator
#'
#' @description `.cens()` marks a variable as a censoring indicator rather than a treatment, so that [bal.tab()] assesses the balance of the units still under observation against the full at-risk sample. It is most often used on the left side of a formula, as in `.cens(C) ~ x1 + x2`, but it can also be called directly and supplied to `bal.tab()`'s `treat` argument.
#'
#' @param x a censoring indicator: a numeric variable taking only the values 0 (still under observation) and 1 (censored), a logical variable, or a factor with levels `0`/`1` or `FALSE`/`TRUE`. Missing values are allowed and preserved.
#'
#' @returns `x` coerced to a 0/1 numeric vector of class `treat` (see [`treat-class`]) with a `"treat.type"` attribute of `"censoring"`. Any value other than 0, 1, or `NA` throws an error.
#'
#' Inside a formula the marker is stripped before the formula is processed, so `.cens()` is not actually evaluated there and the treatment name remains that of the indicator itself (e.g., `C` rather than `.cens(C)`).
#'
#' @details
#' Censoring is its own treatment type in \pkg{cobalt}, distinct from binary, multi-category, and continuous treatments, because the question it asks is different. There is no second treatment group to compare against: the target is the full at-risk sample, and what matters is whether the units still under observation resemble it once weighted. See [`class-bal.tab.cens`] for the output this produces and the arguments that control it.
#'
#' This function is deliberately identical in name and contract to `WeightIt::.cens()`, so that the same code works whichever package is attached; \pkg{cobalt} defines its own only to avoid depending on \pkg{WeightIt}. `bal.tab()` recognizes an indicator tagged by either.
#'
#' @section Note on the coding convention:
#' The survival convention is used: `1` means the unit is *censored* (drops out of observation) and `0` means it remains under observation. This is the opposite of an "observed" or "event" indicator.
#'
#' Missing values are permitted because a unit censored at an earlier time point has no later indicator. Units with a missing indicator are not at risk, so they are in neither sample.
#'
#' @seealso
#' * [`class-bal.tab.cens`] for the output of `bal.tab()` with a censoring indicator
#' * [bal.tab()]
#'
#' @examplesIf rlang::is_installed("WeightIt")
#' data("lalonde", package = "cobalt")
#'
#' # A censoring indicator: 1 = lost to follow-up
#' set.seed(1234)
#' lalonde$C <- rbinom(nrow(lalonde), 1,
#'                     prob = plogis(-1.5 + 0.05 * lalonde$age))
#'
#' # Inverse probability of censoring weights
#' W <- WeightIt::weightit(.cens(C) ~ age + educ + race,
#'                         data = lalonde, method = "glm")
#'
#' # Balance of the weighted uncensored units against the
#' # full at-risk sample
#' bal.tab(W, un = TRUE)
#'
#' # The same thing without WeightIt, using the weights directly
#' bal.tab(.cens(C) ~ age + educ + race, data = lalonde,
#'         weights = W$weights, un = TRUE)

#' @rdname cens
#' @export
.cens <- function(x) {
  arg::arg_supplied(x)

  treat.name <- deparse1(substitute(x))

  out <- .make_cens_treat(x)

  attr(out, "treat.type") <- "censoring"
  attr(out, "treat.name") <- treat.name

  #See `?treat-class`: `treat` is the class every processed treatment has, and the one
  #`WeightIt::.cens()` uses, so an indicator tagged by either package is accepted by
  #both.
  set_class(out, c("cobalt.treat", "treat"), .replace = FALSE, .last = FALSE)
}

#Coerce a censoring indicator to a plain 0/1 numeric vector, stripped of attributes
#(1 = censored, 0 = still at risk). `NA`s are preserved: a unit censored at an earlier
#time point has no later indicator.
.make_cens_treat <- function(x) {
  C <- {
    if (is.factor(x) || is.character(x)) suppressWarnings(as.numeric(as.character(x)))
    else if (is.logical(x)) as.numeric(x)
    else suppressWarnings(as.numeric(x))
  }

  #A coercion that introduces new NAs (a "Yes"/"No" factor, say) is an error, not a
  #silently all-missing indicator.
  if (any(is.na(C) & !is.na(x)) || !all(na.rem(C) %in% c(0, 1))) {
    arg::err(c("a censoring indicator must contain only the values {.val {0}} (still at risk) and {.val {1}} (censored).",
               "i" = "logical values are also allowed, as are factors whose levels are {.val {0}} and {.val {1}} or {.val {FALSE}} and {.val {TRUE}}"))
  }

  if (all(is.na(C))) {
    arg::err("the censoring indicator has no non-missing values")
  }

  C
}

#Is this `.cens(C) ~ rhs`? Matches the marker from either package, since a formula
#written for `WeightIt::weightit()` should work here unchanged.
.is_cens_formula <- function(f) {
  if (!rlang::is_formula(f, lhs = TRUE)) {
    return(FALSE)
  }

  lhs <- rlang::f_lhs(f)

  if (!is.call(lhs)) {
    return(FALSE)
  }

  #`WeightIt::.cens(C)` and `cobalt::.cens(C)` are calls to `::`, whose third element
  #is the name.
  fun <- lhs[[1L]]

  if (is.call(fun) && identical(fun[[1L]], quote(`::`))) {
    fun <- fun[[3L]]
  }

  identical(fun, quote(.cens))
}

#`.cens(C) ~ rhs` -> `C ~ rhs`. The marker is stripped rather than evaluated so that the
#treatment keeps the indicator's own name, and so that `terms()` never sees a call it
#would try to resolve as a variable.
.uncens_formula <- function(f) {
  lhs <- rlang::f_lhs(f)

  if (length(lhs) != 2L) {
    arg::err("{.fun .cens} must wrap exactly one censoring indicator, as in {.code .cens(C) ~ x1 + x2}")
  }

  rlang::new_formula(lhs[[2L]], rlang::f_rhs(f), env = rlang::f_env(f))
}

#The censoring marker, removed from a formula if it is there. Returns the formula to
#work with and whether the treatment it names is a censoring indicator. Every entry
#point that reads a treatment out of a formula calls this first.
.strip_cens <- function(f) {
  if (!.is_cens_formula(f)) {
    return(list(f = f, censoring = FALSE))
  }

  list(f = .uncens_formula(f), censoring = TRUE)
}

#The name a censoring model goes by. `WeightIt::weightitMSM()` records the deparsed left
#side of each censoring formula, marker and all (`WeightIt::.cens(C)`); everywhere else in
#cobalt a censoring indicator goes by its own name, which is also what a user would hand
#to `which.time`. Anything that is not a marked left side comes back untouched.
.uncens_name <- function(x) {
  vapply(x, function(nm) {
    f <- tryCatch(rlang::new_formula(str2lang(nm), 1), error = function(e) NULL)

    if (is_null(f) || !.is_cens_formula(f)) {
      return(nm)
    }

    deparse1(rlang::f_lhs(.uncens_formula(f)))
  }, character(1L), USE.NAMES = FALSE)
}

#Is this treatment a censoring indicator? Reads the tag `.cens()` leaves, from either
#package, without requiring the vector to have been processed.
.is_cens <- function(treat) {
  identical(get.treat.type(treat), "censoring")
}
