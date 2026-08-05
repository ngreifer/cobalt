#Regenerates the on-disk test fixtures in this directory.
#
#Run from the package root with:
#   Rscript tests/testthat/fixtures/make-fixtures.R
#
#Only `cem` gets a saved fixture. Every other supported package is fitted live in
#`helper-fixtures.R`, because cobalt reads other packages' object internals and a
#frozen fixture would keep passing after an upstream layout change. `cem` is
#exempted because it is a frequent install failure (it is not installable on the
#maintainer's machine), and without a fixture `x2base.cem.match()` and
#`get.w.cem.match()` -- roughly 140 lines -- cannot be reached at all.
#
#IMPORTANT: if you have `cem` installed, run this script. It will fit a real
#`cem::cem()` model and overwrite the stand-in below with genuine output. Please
#commit the result. The stand-in exists only so the cem code paths get *some*
#coverage on machines without `cem`; it is not a substitute for the real thing.

library(cobalt)
data("lalonde", package = "cobalt")

out_dir <- file.path("tests", "testthat", "fixtures")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

make_cem_match <- function() {
  vars <- c("age", "educ", "re74")

  if (requireNamespace("cem", quietly = TRUE)) {
    message("cem is installed; fitting a real cem::cem() model.")

    m <- cem::cem(treatment = "treat", data = lalonde[c("treat", vars)],
                  drop = "treat", keep.all = TRUE)

    return(m)
  }

  message("cem is NOT installed; building a documented stand-in instead.")

  #Coarsen each matching variable into bins the way cem does (Sturges' rule for
  #the number of breaks), form strata from the combinations, then weight so that
  #the control units in each stratum reproduce the treated distribution. This is
  #the coarsened exact matching computation, so the resulting weights are real
  #rather than arbitrary.
  coarsen <- function(x) {
    br <- pretty(range(x), n = nclass.Sturges(x))
    cut(x, breaks = br, include.lowest = TRUE, labels = FALSE)
  }

  strata <- interaction(lapply(lalonde[vars], coarsen), drop = TRUE)
  treat <- lalonde$treat

  n_t <- tapply(treat == 1, strata, sum)
  n_c <- tapply(treat == 0, strata, sum)

  #A stratum is matched only if it contains both a treated and a control unit.
  ok <- !is.na(n_t) & !is.na(n_c) & n_t > 0 & n_c > 0
  matched <- ok[as.character(strata)]
  matched[is.na(matched)] <- FALSE

  #Treated units keep a weight of 1; controls are scaled within stratum so that
  #each stratum's control total equals its treated total.
  w <- rep(0, length(treat))
  w[matched & treat == 1] <- 1

  is_c <- matched & treat == 0
  w[is_c] <- (n_t[as.character(strata)][is_c] / n_c[as.character(strata)][is_c])

  #`mstrata` is NA for unmatched units, as in cem's output.
  mstrata <- as.integer(factor(strata))
  is.na(mstrata[!matched]) <- TRUE

  structure(
    list(
      #The components cobalt reads: see x2base.cem.match() and get.w.cem.match().
      w = w,
      groups = factor(treat),
      vars = vars,
      baseline.group = "1",
      mstrata = mstrata,
      matched = matched,
      tab = rbind(All = c(sum(treat == 0), sum(treat == 1)),
                  Matched = c(sum(matched & treat == 0), sum(matched & treat == 1))),
      #Recorded so a reader can tell a stand-in from real cem output.
      .cobalt_standin = TRUE
    ),
    class = "cem.match")
}

cem_match <- make_cem_match()

stopifnot(inherits(cem_match, "cem.match"),
          !all(cobalt:::check_if_zero(cem_match[["w"]])),
          length(cem_match[["w"]]) == nrow(lalonde))

saveRDS(cem_match, file.path(out_dir, "cem_match.rds"), version = 2L)

message(sprintf("wrote %s (%d matched of %d units)",
                file.path(out_dir, "cem_match.rds"),
                sum(cem_match[["w"]] > 0), nrow(lalonde)))
