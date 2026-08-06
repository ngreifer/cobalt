#' @title The `treat` Class
#' @name treat-class
#'
#' @description
#' \pkg{cobalt} represents a processed treatment variable as an object of class `treat`: the treatment values themselves, carrying attributes that record what kind of treatment it is and what its groups are called. Every `bal.tab()` method builds one, and every balance computation reads it rather than re-deriving the same facts from the values.
#'
#' The attributes are
#'
#' \describe{
#'     \item{`treat.type`}{one of `"binary"`, `"multinomial"`, `"continuous"`, or `"censoring"`. This is what decides which balance statistics apply and how the treatment is compared.}
#'     \item{`treat.name`}{the name the treatment was given, if it had one.}
#'     \item{`treat_names`}{for a binary treatment, the display names of the two groups, named `control` and `treated`; for a multi-category treatment, one name per level; for a censoring indicator, `Uncensored` and `Censored`. `NULL` for a continuous treatment.}
#'     \item{`treat_vals`}{the values in the data that those names correspond to.}
#'     \item{`group_labels`}{how the treatment's groups are labelled in the column names of a balance table. `c("0", "1")` for a binary treatment, which is \pkg{cobalt}'s positional convention (`M.0` is the control group's mean); otherwise whatever names the groups have.}
#' }
#'
#' `[` preserves all of them, so a subset of a `treat` is still a `treat`.
#'
#' [.cens()] is the only way to build one directly; it returns a censoring indicator. `WeightIt::.cens()` returns an object with the same class and a compatible `treat.type`, so an indicator tagged by either package is accepted here.
#'
#' The class is shared, so carrying it means only "this is a treatment", not "this carries the attributes above": \pkg{WeightIt} tags its own treatments with it and supplies `treat.type` alone, and even \pkg{cobalt}'s [.cens()] has no group names until the treatment is processed. Anything that reads `treat_names` or `treat_vals` therefore has to check for them rather than for the class, and process the treatment when they are absent.
#'
#' @seealso
#' * [.cens()]
#' * [`class-bal.tab.cens`]
#'
NULL

#Every attribute a `treat` carries. `[.cobalt.treat` copies each one, so a new
#attribute added here survives subsetting without anything else changing.
TREAT_ATTRS <- c("treat.type", "treat.name", "treat_names", "treat_vals",
                 "group_labels")

#The method is registered on `cobalt.treat` rather than on `treat`, even though `treat`
#is the shared class. \pkg{WeightIt} 2.0.0 registers its own `[.treat`, and two packages
#registering the same method means whichever loads second overwrites the other and R
#says so out loud -- on every `library(WeightIt)`, since it imports cobalt. Dispatching
#one class earlier avoids that without giving up the shared contract; when WeightIt
#adopts this class the extra layer can go.
#' @exportS3Method `[` cobalt.treat
`[.cobalt.treat` <- function(x, ...) {
  y <- NextMethod("[")

  for (a in TREAT_ATTRS) {
    attr(y, a) <- .attr(x, a)
  }

  set_class(y, class(x))
}

#What a plot calls each unit's group. An ordinary treatment is labelled by its own
#values; one whose groups are named something more informative than their values -- the
#two samples of a censoring comparison -- is labelled by those names.
treat_labels <- function(treat) {
  labels <- unname(group_labels(treat))

  if (identical(labels, c("0", "1"))) {
    return(treat_vals(treat))
  }

  setNames(labels, names(treat_vals(treat)))
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
#Whether `treat` already carries the group bookkeeping that `process_treat()` adds.
#The `treat` class on its own does not answer this: it is shared with \pkg{WeightIt}
#(see `?treat-class`), whose treatments carry `treat.type` but none of the group
#attributes, and cobalt's own [.cens()] does not add them either until the treatment is
#processed. Anything about to read `treat_names()` or `treat_vals()` must ask this
#rather than `inherits(treat, "treat")`, and process the treatment when it is `FALSE`.
#Continuous treatments have no groups, so they are never "processed" in this sense;
#every caller here is on a branch that has already excluded them.
is_processed_treat <- function(treat) {
  inherits(treat, "treat") && is_not_null(treat_vals(treat))
}
`group_labels<-` <- function(treat, value) {
  `attr<-`(treat, "group_labels", value)
}
#How a treatment's groups are labelled in a balance table's column names. A binary
#treatment uses `c("0", "1")` positionally -- `M.0` has always meant the control
#group's mean -- so that is the default rather than the group names; a treatment that
#wants its groups named says so.
group_labels <- function(treat) {
  .attr(treat, "group_labels") %or% c("0", "1")
}
#The two groups of a binary treatment as a balance table labels them: the label that
#goes into the column name mapped to the value that selects the group.
.bin_groups <- function(treat) {
  setNames(treat_vals(treat)[treat_names(treat)[c("control", "treated")]],
           group_labels(treat))
}
get.treat.type <- function(treat) {
  out <- .attr(treat, "treat.type")

  if (identical(out, "multi-category")) {
    return("multinomial")
  }

  out
}
has.treat.type <- function(treat) {
  is_not_null(get.treat.type(treat))
}
assign.treat.type <- function(treat, use.multi = FALSE, censoring = NULL) {
  #Returns treat with treat.type attribute
  if (is_null(censoring)) {
    censoring <- .is_cens(treat)
  }

  #A censoring indicator is validated rather than classified, and keeps its tag. Unlike
  #a treatment it may take a single value -- a time point at which no unit (or every
  #unit) is censored is degenerate but not malformed -- and it may be missing.
  if (censoring) {
    .make_cens_treat(treat)

    attr(treat, "treat.type") <- "censoring"

    return(treat)
  }

  nunique.treat <- nunique(treat)

  if (nunique.treat < 2L) {
    arg::err("the treatment must have at least two unique values")
  }

  if (!use.multi && nunique.treat == 2L) {
    treat.type <- "binary"
  }
  else if (use.multi || is.character(treat) || is.factor(treat)) {
    treat.type <- "multinomial"

    #Whether it needs coercing, not who owns the class: a processed treatment is
    #already a factor, and re-factoring one would drop the attributes it carries.
    if (!is.factor(treat)) {
      treat <- factor(treat)
    }
  }
  else {
    treat.type <- "continuous"
  }

  attr(treat, "treat.type") <- treat.type

  treat
}

process_treat <- function(treat, ..., keep_values = FALSE, .missing.okay = FALSE) {

  arg::arg_supplied(treat)

  if (inherits(treat, "unprocessed.treat")) {
    attrs <- attributes(treat)
    renamed_original <- setNames(names(treat_vals(treat)), treat_vals(treat))
    treat <- factor(renamed_original[as.character(treat)], levels = renamed_original)

    for (at in c(TREAT_ATTRS, "keep_values", "names")) {
      attr(treat, at) <- attrs[[at]]
    }
  }
  else if (.is_cens(treat)) {
    treat <- .process_cens_treat(treat, ...)
  }
  else {
    # keep_values <- isTRUE(attr(treat, "keep_values")) ||
    #     (has.treat.type(treat) && get.treat.type(treat) == "multinomial")

    treat <- .process_vector(treat, name = "treat",
                             which = "treatment statuses",
                             datalist = list(...), missing.okay = .missing.okay) |>
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

  #The class goes first so that `[.cobalt.treat` is reached before `[.factor`; a binary
  #or multi-category treatment is a factor underneath, and `[.factor` would drop every
  #attribute the class exists to carry. The method calls `NextMethod()`, so the factor
  #method still runs and the levels still survive.
  set_class(treat, c("cobalt.treat", "treat"), .replace = FALSE, .last = FALSE)
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

  #A treatment is undefined for a unit no longer under observation, so once a censoring
  #indicator is anywhere in the list the treatments are allowed to be missing. Which
  #units that covers is not settled here: `base.bal.tab.msm()` accumulates the risk set
  #from the indicators themselves and errors on a unit still at risk whose treatment is
  #unknown. With no indicator in the list a missing treatment stays an error, as for a
  #point treatment.
  .missing.okay <- any(vapply(treat.list, .is_cens, logical(1L)))

  lapply(treat.list, process_treat, ..., .missing.okay = .missing.okay) |>
    setNames(treat.list.names)
}

#A censoring indicator stays a plain 0/1 vector -- `C == 0` has to keep working, and
#nothing about it is ordered or factor-like -- but it is named like the other types so
#that anything reading `treat_names()` gets an answer. The two samples balance is
#actually assessed between are built later, by `.cens_pseudo_treat()`.
.process_cens_treat <- function(treat, ...) {
  treat.name <- .attr(treat, "treat.name")

  treat <- .process_vector(treat, name = "treat",
                           which = "censoring statuses",
                           datalist = list(...), missing.okay = TRUE) |>
    .make_cens_treat()

  attr(treat, "treat.type") <- "censoring"
  attr(treat, "treat.name") <- treat.name

  treat_names(treat) <- setNames(c("Uncensored", "Censored"),
                                 c("control", "treated"))
  treat_vals(treat) <- setNames(c(0, 1), treat_names(treat))
  group_labels(treat) <- unname(treat_names(treat))

  treat
}

subset_treat <- function(x, index) {
  y <- x[index]

  #A censoring indicator has no levels to narrow -- a subset in which nobody is
  #censored is still a censoring indicator -- and reclassifying it would turn it back
  #into a binary treatment.
  if (.is_cens(x)) {
    return(y)
  }

  treat_names(y) <- treat_names(x)[treat_vals(x) %in% unique(y)]
  treat_vals(y) <- treat_vals(x)[treat_vals(x) %in% unique(y)]

  assign.treat.type(y) |>
    set_class(class(x))
}
