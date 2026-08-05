#Golden-output harness for the bal.tab() refactor.
#
#The test suite is a strong *shape* gate and a weak *value* gate: test-input-guards.R
#asserts only that each fitted-object class yields a non-empty *.Balance and prints.
#The stages that touch the 21 x2base adapters therefore need something that pins the
#actual numbers and the actual printed characters. That is what this does.
#
#Usage, from the package root:
#
#  Rscript -e 'pkgload::load_all("."); source("_dev/refactor-golden.R"); golden_build()'
#  ... make a change ...
#  Rscript -e 'pkgload::load_all("."); source("_dev/refactor-golden.R"); compare_golden()'
#
#`golden_build()` writes one .rds per cell under _dev/golden/ (gitignored).
#`compare_golden()` reloads and reports the first difference in each cell.
#
#Cells whose fixture cannot be built are recorded as "unavailable" and are compared
#only when both sides agree they are unavailable, so a missing Suggests package
#never looks like a regression.

GOLDEN_DIR <- file.path("_dev", "golden")

# ---- fixture access outside testthat ---------------------------------------

#helper-fixtures.R calls skip_if_not_installed()/skip()/test_path(), which only
#exist under testthat. Stub them so the registry can be sourced as plain code.
.golden_load_fixtures <- function() {
  env <- new.env(parent = globalenv())

  env$skip_if_not_installed <- function(pkg, ...) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("package '", pkg, "' is not installed", call. = FALSE)
    }
  }
  env$skip <- function(message = "") stop(message, call. = FALSE)
  env$test_path <- function(...) file.path("tests", "testthat", ...)

  sys.source(file.path("tests", "testthat", "helper-fixtures.R"), envir = env)

  env
}

.FX <- NULL

#Build a fixture, returning the condition instead of signalling if it fails.
gfx <- function(name) {
  if (is.null(.FX)) {
    .FX <<- .golden_load_fixtures()
  }

  tryCatch(suppressMessages(suppressWarnings(.FX$fx(name))),
           error = function(e) e,
           condition = function(e) e)
}

# ---- capture ----------------------------------------------------------------

#Everything a print() call emits: stdout, messages, warnings, and an error if one
#is thrown. Conditions are captured as text so their wording is pinned too.
.capture_all <- function(expr) {
  msgs <- character()
  wrns <- character()
  err <- NULL

  out <- tryCatch(
    utils::capture.output(
      withCallingHandlers(
        force(expr),
        message = function(m) {
          msgs <<- c(msgs, conditionMessage(m))
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          wrns <<- c(wrns, conditionMessage(w))
          invokeRestart("muffleWarning")
        })),
    error = function(e) {
      err <<- conditionMessage(e)
      character()
    })

  list(output = out,
       messages = .squish(msgs),
       warnings = .squish(wrns),
       error = .squish(err))
}

#Collapse whitespace so cli's console-width wrapping does not register as a diff.
.squish <- function(x) {
  if (length(x) == 0L) {
    return(character())
  }

  gsub("\\s+", " ", trimws(x))
}

# ---- the print-argument variants --------------------------------------------

#Argument sets covering the display logic print() resolves: the `un`/`disp` toggles,
#row filtering, the call block, and a non-default `digits`.
#
#The last four exist to pin column *selection* specifically, which is otherwise only
#lightly exercised: turning a statistic off must drop its column and its threshold
#column, moments must appear in both the unadjusted and adjusted blocks, and naming a
#single aggregating function must drop the other two from a summary table. An argument
#that does not apply to a given shape is ignored by `print()`, and one that is invalid
#for a cell raises a condition that `.capture_all()` records like any other.
.golden_print_variants <- function() {
  list(default = list(),
       un = list(un = TRUE),
       means_sds = list(disp = c("means", "sds")),
       imbalanced = list(imbalanced.only = TRUE),
       call = list(disp.call = TRUE),
       digits5 = list(digits = 5L),
       diff_off = list(disp.diff = FALSE),
       corr_off = list(disp.corr = FALSE),
       un_means = list(un = TRUE, disp = "means"),
       one_agg_fun = list(cluster.fun = "max", imp.fun = "max"))
}

# ---- the grid ---------------------------------------------------------------

#Some objects do not carry their covariates, so `bal.tab()` needs `formula`/`data`.
#Mirrors `.extra_args()` in tests/testthat/test-input-guards.R.
.golden_extra_args <- function(nm) {
  if (identical(nm, "cem_match")) {
    return(list(data = lalonde))
  }

  if (nm %in% c("Match", "optmatch", "ebalance", "designmatch")) {
    return(list(formula = treat ~ age + educ + re74, data = lalonde))
  }

  list()
}

#Each entry is a zero-argument function returning a bal.tab object, so a failure to
#build is caught per cell rather than aborting the sweep. Fixtures are referenced
#through gfx(), which returns a condition when a Suggests package is missing.
golden_grid <- function() {
  covs <- lalonde[c("age", "educ", "race", "married", "nodegree", "re74", "re75")]
  covs_nr <- covs[setdiff(names(covs), "race")]
  treat <- lalonde$treat

  w <- .FX$w_fixed
  sw <- .FX$sw_fixed
  cl <- .FX$cl_idx
  sub <- .FX$sub_idx
  imp <- .FX$imp_idx
  cens <- .FX$cens_idx

  bin <- c("mean.diffs", "variance.ratios", "ks.statistics", "ovl.coefficients")
  cont <- c("correlations", "spearman.correlations", "mean.diffs.target",
            "ks.statistics.target")

  g <- list()

  ## ---- data.frame method, option sweep -------------------------------------

  g$df_default <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled")
  }
  g$df_un <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE)
  }
  g$df_all_stats <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE,
            stats = bin, disp = c("means", "sds"))
  }
  g$df_quick_false <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE,
            quick = FALSE)
  }
  g$df_thresholds <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE,
            stats = bin, thresholds = c(m = .1, v = 2, ks = .05, ovl = .1))
  }
  g$df_abs <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE,
            stats = bin, abs = TRUE)
  }
  g$df_two_weights <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", un = TRUE,
            weights = list(a = w, b = sw), stats = c("mean.diffs", "ks.statistics"),
            thresholds = c(m = .1))
  }
  g$df_two_weights_methods <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", un = TRUE,
            weights = list(a = w, b = sw), method = c("matching", "weighting"))
  }
  g$df_s_weights <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, s.weights = sw,
            un = TRUE, disp = c("means", "sds"))
  }
  g$df_distance <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE,
            distance = data.frame(prop.score = seq_len(nrow(covs)) / nrow(covs)))
  }
  g$df_addl <- function() {
    bal.tab(covs[c("age", "educ")], treat = treat, s.d.denom = "pooled",
            weights = w, un = TRUE, addl = covs["race"])
  }
  g$df_int_poly <- function() {
    bal.tab(covs[c("age", "educ")], treat = treat, s.d.denom = "pooled",
            weights = w, un = TRUE, int = TRUE, poly = 2)
  }
  g$df_subset <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE,
            subset = seq_len(nrow(covs)) %% 3L != 0L)
  }
  g$df_imbalanced_only <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE,
            thresholds = c(m = .05), imbalanced.only = TRUE)
  }
  g$df_binary_std <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE,
            binary = "std", continuous = "raw")
  }

  for (sdd in c("pooled", "treated", "control", "all", "weighted", "hedges")) {
    local({
      .sdd <- sdd
      g[[paste0("df_sdd_", .sdd)]] <<- function() {
        bal.tab(covs, treat = treat, binary = "std", continuous = "std",
                s.d.denom = .sdd, weights = w, s.weights = sw, un = TRUE)
      }
    })
  }

  ## ---- formula and matrix interfaces ---------------------------------------

  g$formula <- function() {
    bal.tab(treat ~ age + educ + race + married + nodegree + re74 + re75,
            data = lalonde, s.d.denom = "pooled", weights = w, un = TRUE)
  }
  g$matrix <- function() {
    bal.tab(as.matrix(splitfactor(covs, drop.first = FALSE)), treat = treat,
            s.d.denom = "pooled", weights = w, un = TRUE)
  }
  #`treat` and `cluster` named as strings rather than passed as vectors. These take
  #a different resolution path, and `treat` in particular is a name that also
  #appears in `...` the whole way down the call chain.
  g$treat_as_string <- function() {
    bal.tab(covs, treat = "treat", data = lalonde, s.d.denom = "pooled",
            weights = w, un = TRUE, stats = c("m", "ks"), thresholds = c(m = .1))
  }
  g$cluster_as_string <- function() {
    bal.tab(reformulate(names(covs_nr), "treat"), data = lalonde,
            s.d.denom = "pooled", weights = w, cluster = "race", un = TRUE)
  }

  ## ---- continuous treatment ------------------------------------------------

  g$cont_default <- function() {
    bal.tab(covs_nr, treat = lalonde$re75, weights = w, un = TRUE)
  }
  g$cont_all_stats <- function() {
    bal.tab(covs_nr[setdiff(names(covs_nr), "re75")], treat = lalonde$re75,
            weights = w, un = TRUE, stats = cont, disp = c("means", "sds"))
  }
  g$cont_thresholds <- function() {
    bal.tab(covs_nr[setdiff(names(covs_nr), "re75")], treat = lalonde$re75,
            weights = w, un = TRUE, stats = c("correlations", "mean.diffs.target"),
            thresholds = c(r = .1, m = .1))
  }
  g$cont_quick_false <- function() {
    bal.tab(covs_nr[setdiff(names(covs_nr), "re75")], treat = lalonde$re75,
            weights = w, un = TRUE, quick = FALSE)
  }

  ## ---- multi-category treatment --------------------------------------------

  g$multi_default <- function() {
    bal.tab(covs_nr, treat = lalonde$race, s.d.denom = "pooled", weights = w,
            un = TRUE)
  }
  g$multi_all_stats <- function() {
    bal.tab(covs_nr, treat = lalonde$race, s.d.denom = "pooled", weights = w,
            un = TRUE, stats = bin, disp = c("means", "sds"))
  }
  g$multi_pairwise_false <- function() {
    bal.tab(covs_nr, treat = lalonde$race, s.d.denom = "pooled", weights = w,
            un = TRUE, pairwise = FALSE)
  }
  g$multi_focal <- function() {
    bal.tab(covs_nr, treat = lalonde$race, s.d.denom = "focal", focal = "white",
            weights = w, un = TRUE)
  }
  g$multi_thresholds <- function() {
    bal.tab(covs_nr, treat = lalonde$race, s.d.denom = "pooled", weights = w,
            un = TRUE, thresholds = c(m = .1), multi.summary = TRUE)
  }
  g$bin_pairwise_false <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, un = TRUE,
            pairwise = FALSE)
  }

  ## ---- cluster --------------------------------------------------------------

  g$cluster_default <- function() {
    bal.tab(covs_nr, treat = treat, s.d.denom = "pooled", weights = w,
            cluster = cl, un = TRUE)
  }
  g$cluster_summary <- function() {
    bal.tab(covs_nr, treat = treat, s.d.denom = "pooled", weights = w,
            cluster = cl, un = TRUE, cluster.summary = TRUE)
  }
  g$cluster_fun_range <- function() {
    bal.tab(covs_nr, treat = treat, s.d.denom = "pooled", weights = w,
            cluster = cl, un = TRUE, cluster.summary = TRUE,
            cluster.fun = c("min", "mean", "max"))
  }
  g$cluster_quick_false <- function() {
    bal.tab(covs_nr, treat = treat, s.d.denom = "pooled", weights = w,
            cluster = cl, un = TRUE, quick = FALSE)
  }
  g$cluster_cont <- function() {
    bal.tab(covs_nr[setdiff(names(covs_nr), "re75")], treat = lalonde$re75,
            weights = w, cluster = cl, un = TRUE, cluster.summary = TRUE)
  }

  ## ---- imputation -----------------------------------------------------------

  g$imp_default <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, imp = imp,
            un = TRUE)
  }
  g$imp_summary <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, imp = imp,
            un = TRUE, imp.summary = TRUE, imp.fun = c("min", "mean", "max"))
  }
  g$imp_quick_false <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", weights = w, imp = imp,
            un = TRUE, quick = FALSE)
  }

  ## ---- subclass -------------------------------------------------------------

  g$subclass_default <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", subclass = sub, un = TRUE)
  }
  g$subclass_summary <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", subclass = sub, un = TRUE,
            subclass.summary = TRUE, disp.subclass = TRUE)
  }
  g$subclass_thresholds <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", subclass = sub, un = TRUE,
            stats = c("mean.diffs", "variance.ratios"),
            thresholds = c(m = .1, v = 2), disp.subclass = TRUE)
  }
  g$subclass_quick_false <- function() {
    bal.tab(covs, treat = treat, s.d.denom = "pooled", subclass = sub, un = TRUE,
            stats = c("mean.diffs", "variance.ratios"), quick = FALSE,
            disp.subclass = TRUE)
  }
  g$subclass_cont <- function() {
    bal.tab(covs, treat = lalonde$re75, subclass = sub, un = TRUE,
            subclass.summary = TRUE, disp.subclass = TRUE,
            disp = c("means", "sds"), thresholds = c(cor = .1))
  }
  g$subclass_cont_quick_false <- function() {
    bal.tab(covs, treat = lalonde$re75, subclass = sub, un = TRUE, quick = FALSE,
            stats = c("correlations", "spearman.correlations"),
            subclass.summary = TRUE, disp.subclass = TRUE)
  }

  ## ---- censoring ------------------------------------------------------------

  #The censoring weights, as a censoring model produces them: zero once censored.
  cens_w <- local({
    p <- fitted(glm(cens ~ age + educ + race, data = lalonde, family = binomial))
    ifelse(cens == 0, 1 / (1 - p), 0)
  })

  g$cens_default <- function() {
    bal.tab(covs, treat = .cens(cens), weights = cens_w, un = TRUE)
  }
  g$cens_all_stats <- function() {
    bal.tab(covs, treat = .cens(cens), weights = cens_w, un = TRUE,
            stats = bin, disp = c("means", "sds"),
            thresholds = c(m = .1, v = 2, ks = .05, ovl = .1))
  }
  g$cens_s_weights <- function() {
    bal.tab(covs, treat = .cens(cens), weights = cens_w, s.weights = sw, un = TRUE,
            s.d.denom = "uncensored", quick = FALSE)
  }
  g$cens_cluster <- function() {
    bal.tab(covs_nr, treat = .cens(cens), weights = cens_w, cluster = cl, un = TRUE,
            cluster.summary = TRUE)
  }
  g$cens_imp <- function() {
    bal.tab(covs_nr, treat = .cens(cens), weights = cens_w, imp = imp, un = TRUE,
            imp.summary = TRUE)
  }
  g$cens_subclass <- function() {
    bal.tab(covs, treat = .cens(cens), subclass = sub, un = TRUE,
            subclass.summary = TRUE, disp.subclass = TRUE,
            disp = c("means", "sds"), thresholds = c(m = .1))
  }
  g$cens_subclass_cluster <- function() {
    bal.tab(covs_nr, treat = .cens(cens), subclass = sub, cluster = cl, un = TRUE)
  }
  #The `weightit_cens` fixture is picked up by the per-object loop below, which covers
  #the object interface and its `quick = FALSE` view.

  ## ---- longitudinal ---------------------------------------------------------

  g$msm_formula_list <- function() {
    bal.tab(list(treat ~ age + educ,
                 nodegree ~ age + educ + treat),
            data = lalonde, s.d.denom = "pooled", weights = w, un = TRUE)
  }
  g$msm_summary <- function() {
    bal.tab(list(treat ~ age + educ,
                 nodegree ~ age + educ + treat),
            data = lalonde, s.d.denom = "pooled", weights = w, un = TRUE,
            msm.summary = TRUE)
  }
  g$msm_matching <- function() {
    bal.tab(list(treat ~ age + educ,
                 nodegree ~ age + educ + treat),
            data = lalonde, s.d.denom = "pooled", weights = w, un = TRUE,
            method = "matching")
  }
  #The cells above pass the value each time point would have defaulted to anyway.
  #These pass a different one, which used to be discarded.
  g$msm_sdd_treated <- function() {
    bal.tab(list(treat ~ age + educ,
                 nodegree ~ age + educ + treat),
            data = lalonde, s.d.denom = "treated", weights = w, un = TRUE,
            msm.summary = TRUE)
  }
  g$msm_cont_sdd <- function() {
    bal.tab(list(re75 ~ age + educ,
                 re78 ~ age + educ + re75),
            data = lalonde, s.d.denom = "weighted", weights = w, un = TRUE,
            msm.summary = TRUE)
  }

  ## ---- nested shapes --------------------------------------------------------

  g$nest_cluster_imp <- function() {
    bal.tab(covs_nr, treat = treat, s.d.denom = "pooled", weights = w,
            cluster = cl, imp = imp, un = TRUE)
  }
  g$nest_cluster_multi <- function() {
    bal.tab(covs_nr, treat = lalonde$race, s.d.denom = "pooled", weights = w,
            cluster = cl, un = TRUE)
  }
  g$nest_cluster_subclass <- function() {
    bal.tab(covs_nr, treat = treat, s.d.denom = "pooled", subclass = sub,
            cluster = cl, un = TRUE)
  }
  g$nest_multi_imp <- function() {
    bal.tab(covs_nr, treat = lalonde$race, s.d.denom = "pooled", weights = w,
            imp = imp, un = TRUE)
  }
  g$nest_msm_cluster <- function() {
    bal.tab(list(treat ~ age + educ,
                 nodegree ~ age + educ + treat),
            data = lalonde, s.d.denom = "pooled", weights = w, cluster = cl,
            un = TRUE)
  }

  ## ---- one cell per fitted-object fixture ----------------------------------

  for (nm in .FX$fx_object_names()) {
    local({
      .nm <- nm
      .extra <- .golden_extra_args(.nm)

      g[[paste0("obj_", .nm)]] <<- function() {
        do.call(bal.tab, c(list(gfx(.nm)), .extra, list(un = TRUE)))
      }
      g[[paste0("obj_", .nm, "_quickF")]] <<- function() {
        do.call(bal.tab, c(list(gfx(.nm)), .extra,
                           list(un = TRUE, quick = FALSE, disp = c("means", "sds"))))
      }
    })
  }

  #A few object cells with extra arguments that reach otherwise-untested branches.
  g$obj_matchit_cluster <- function() {
    bal.tab(gfx("matchit"), cluster = cl, un = TRUE, cluster.summary = TRUE)
  }
  #The `matchit_discard`/`matchit_sub_discard` fixtures are picked up by the loop
  #above; this adds the subclass-display view of the latter, which is where the
  #`Discarded` column of the sample size table shows up.
  g$obj_matchit_sub_discard_disp <- function() {
    bal.tab(gfx("matchit_sub_discard"), un = TRUE, disp.subclass = TRUE)
  }
  g$obj_matchit_addl_int <- function() {
    bal.tab(gfx("matchit"), addl = ~ I(age^2), int = TRUE, un = TRUE)
  }
  g$obj_weightit_thresholds <- function() {
    bal.tab(gfx("weightit"), un = TRUE, thresholds = c(m = .05, v = 2),
            stats = c("mean.diffs", "variance.ratios"))
  }
  g$obj_ps_imp <- function() {
    bal.tab(gfx("ps"), imp = imp, un = TRUE)
  }

  g
}

#Objects for which love.plot()'s data is captured. Restricted to single-stat plots,
#which return a ggplot whose $plot$data is a plain data.frame.
.golden_loveplot_cells <- function() {
  c("df_default", "df_all_stats", "df_two_weights", "cont_default",
    "multi_default", "cluster_summary", "imp_summary", "subclass_summary",
    "msm_summary", "obj_matchit")
}

# ---- build and compare ------------------------------------------------------

golden_build <- function(dir = GOLDEN_DIR, filter = NULL) {
  if (is.null(.FX)) {
    .FX <<- .golden_load_fixtures()
  }

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  grid <- golden_grid()

  if (!is.null(filter)) {
    grid <- grid[grep(filter, names(grid))]
  }

  love_cells <- .golden_loveplot_cells()

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (nm in names(grid)) {
    cell <- .golden_cell(grid[[nm]], capture_love = nm %in% love_cells)
    saveRDS(cell, file.path(dir, paste0(nm, ".rds")), version = 3L)
    cat(sprintf("%-34s %s\n", nm, if (isTRUE(cell$unavailable)) "unavailable" else "ok"))
  }

  cat(sprintf("\n%d cells written to %s\n", length(grid), dir))

  invisible(NULL)
}

.golden_cell <- function(build, capture_love = FALSE) {
  b <- tryCatch(suppressMessages(suppressWarnings(build())),
                error = function(e) e)

  #A cell whose fixture is missing is recorded as unavailable rather than as an
  #error, so a missing Suggests package is not mistaken for a regression.
  if (inherits(b, "condition")) {
    return(list(unavailable = TRUE, reason = .squish(conditionMessage(b))))
  }

  variants <- lapply(.golden_print_variants(), function(a) {
    .capture_all(do.call(print, c(list(b), a)))
  })

  love <- NULL

  if (capture_love) {
    love <- tryCatch({
      p <- suppressMessages(suppressWarnings(love.plot(b)))
      ggplot2::ggplot_build(p)$plot$data
    }, error = function(e) .squish(conditionMessage(e)))
  }

  list(unavailable = FALSE,
       obj = b,
       variants = variants,
       love = love)
}

compare_golden <- function(dir = GOLDEN_DIR, filter = NULL, max_report = 3L) {
  if (is.null(.FX)) {
    .FX <<- .golden_load_fixtures()
  }

  if (!dir.exists(dir)) {
    stop("no golden set at '", dir, "'; run golden_build() first", call. = FALSE)
  }

  grid <- golden_grid()

  if (!is.null(filter)) {
    grid <- grid[grep(filter, names(grid))]
  }

  love_cells <- .golden_loveplot_cells()

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  n_ok <- 0L
  n_skip <- 0L
  bad <- character()

  for (nm in names(grid)) {
    f <- file.path(dir, paste0(nm, ".rds"))

    if (!file.exists(f)) {
      cat(sprintf("%-34s NEW (no golden file)\n", nm))
      next
    }

    old <- readRDS(f)
    new <- .golden_cell(grid[[nm]], capture_love = nm %in% love_cells)

    diffs <- .golden_diff(old, new)

    if (length(diffs) == 0L) {
      if (isTRUE(new$unavailable)) n_skip <- n_skip + 1L
      else n_ok <- n_ok + 1L
      next
    }

    bad <- c(bad, nm)
    cat(sprintf("\n%-34s DIFF\n", nm))
    for (d in utils::head(diffs, max_report)) {
      cat("    ", d, "\n", sep = "")
    }
    if (length(diffs) > max_report) {
      cat(sprintf("     ... and %d more\n", length(diffs) - max_report))
    }
  }

  cat(sprintf("\nidentical: %d | unavailable: %d | differing: %d\n",
              n_ok, n_skip, length(bad)))

  if (length(bad) > 0L) {
    cat("differing cells:", paste(bad, collapse = ", "), "\n")
  }

  invisible(bad)
}

.golden_diff <- function(old, new) {
  out <- character()

  if (!identical(isTRUE(old$unavailable), isTRUE(new$unavailable))) {
    return(sprintf("availability changed: was %s, now %s",
                   if (isTRUE(old$unavailable)) "unavailable" else "ok",
                   if (isTRUE(new$unavailable)) "unavailable" else "ok"))
  }

  if (isTRUE(new$unavailable)) {
    return(out)
  }

  #Objects: all.equal() tolerates floating-point noise but reports structural and
  #attribute changes, which is exactly the sensitivity wanted here.
  cmp <- all.equal(old$obj, new$obj)

  if (!isTRUE(cmp)) {
    out <- c(out, paste0("obj: ", cmp))
  }

  for (v in names(new$variants)) {
    o <- old$variants[[v]]
    n <- new$variants[[v]]

    if (is.null(o)) {
      out <- c(out, sprintf("print(%s): new variant", v))
      next
    }

    out <- c(out, .golden_diff_capture(o, n, sprintf("print(%s)", v)))
  }

  if (!identical(is.null(old$love), is.null(new$love))) {
    out <- c(out, "love.plot: presence changed")
  }
  else if (!is.null(new$love)) {
    cmp <- all.equal(old$love, new$love)

    if (!isTRUE(cmp)) {
      out <- c(out, paste0("love.plot data: ", cmp))
    }
  }

  out
}

.golden_diff_capture <- function(o, n, label) {
  out <- character()

  if (!identical(o$output, n$output)) {
    i <- which(o$output != n$output |
                 seq_along(o$output) > length(n$output))[1L]

    if (is.na(i)) i <- length(n$output)

    out <- c(out, sprintf("%s line %d:\n       was: %s\n       now: %s",
                          label, i,
                          encodeString(o$output[i] %||% "<absent>"),
                          encodeString(n$output[i] %||% "<absent>")))
  }

  for (k in c("messages", "warnings", "error")) {
    if (!identical(o[[k]], n[[k]])) {
      out <- c(out, sprintf("%s %s:\n       was: %s\n       now: %s",
                            label, k,
                            paste(encodeString(o[[k]]), collapse = " | "),
                            paste(encodeString(n[[k]]), collapse = " | ")))
    }
  }

  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# ---- column-spec assertion (Stage 8) ---------------------------------------

#Once .bal_tab_col_spec() exists, this proves it reproduces every balance table's
#column names before anything is allowed to consume it.
check_col_spec <- function(dir = GOLDEN_DIR) {
  if (!exists(".bal_tab_col_spec", envir = asNamespace("cobalt"), inherits = FALSE)) {
    cat("`.bal_tab_col_spec()` does not exist yet; nothing to check\n")
    return(invisible(NULL))
  }

  spec <- get(".bal_tab_col_spec", envir = asNamespace("cobalt"))

  #The three table builders call the spec directly, so comparing their output to it
  #would be tautological. What Stage 9c and `as.data.frame()` actually need is that
  #the spec be recoverable from `print.options` alone, so that is what is checked:
  #re-derive each table's columns from the object's own `p.ops` and compare.
  all_stats <- get("all_STATS", envir = asNamespace("cobalt"))

  .compute <- function(p, type, intersect_stats) {
    stats <- if (intersect_stats) intersect(all_stats(type), p[["stats"]]) else all_stats(type)

    if (isTRUE(p[["quick"]])) c(p[["disp"]], stats)
    else c("means", "sds", stats)
  }

  #One entry per table shape: how to reach it, and what the spec arguments are.
  .expected <- function(x, p) {
    type <- p[["type"]]
    no.adj <- !isTRUE(p[["disp.adj"]])
    wn <- if (no.adj) "Adj" else p[["weight.names"]]
    thr.s <- if (no.adj) "Un" else wn

    out <- list()

    if (is.data.frame(x[["Balance"]])) {
      out[["Balance"]] <- spec(type, p[["compute"]], p[["thresholds"]],
                               samples = c("Un", wn), threshold.samples = thr.s,
                               group.labels = p[["group.labels"]])
    }

    for (nm in names(x)[startsWith(names(x), "Balance.Across.")]) {
      if (!is.data.frame(x[[nm]])) next

      #Only `Balance.Across.Subclass` is an ordinary balance table; the rest are
      #aggregated summaries carrying an aggregating function per statistic.
      if (nm == "Balance.Across.Subclass") {
        out[[nm]] <- spec(type, .compute(p, type, TRUE), p[["thresholds"]],
                          samples = c("Un", "Adj"), threshold.samples = "Adj",
                          group.labels = p[["group.labels"]])
        next
      }

      agg <- switch(nm,
                    "Balance.Across.Clusters" = p[["cluster.fun"]],
                    "Balance.Across.Imputations" = p[["imp.fun"]],
                    "max")
      agg.all <- if (isTRUE(p[["quick"]])) agg else c("min", "mean", "max")

      out[[nm]] <- spec(type, p[["compute"]], p[["thresholds"]],
                        samples = c("Un", wn), quantities = NULL,
                        agg.funs = agg.all,
                        threshold.samples = if (length(agg) != 1L) character(0L) else thr.s,
                        threshold.agg.fun = agg,
                        include.times = nm == "Balance.Across.Times",
                        group.labels = p[["group.labels"]])
    }

    if (is.list(x[["Subclass.Balance"]]) && !is.data.frame(x[["Subclass.Balance"]])) {
      sb <- spec(type, p[["compute"]], p[["thresholds"]], samples = "Adj",
                 group.labels = p[["group.labels"]])

      for (i in names(x[["Subclass.Balance"]])) {
        out[[paste0("Subclass.Balance[[", i, "]]")]] <- sb
      }
    }

    out
  }

  #`bal.tab` objects nest, so walk down to every level that has its own p.ops.
  .walk <- function(x, path, acc = list()) {
    p <- attr(x, "print.options")

    if (is_not_null(p)) {
      expected <- .expected(x, p)

      for (nm in names(expected)) {
        actual <- {
          if (startsWith(nm, "Subclass.Balance"))
            x[["Subclass.Balance"]][[sub("^.*\\[\\[(.*)\\]\\]$", "\\1", nm)]]
          else x[[nm]]
        }

        acc[[paste0(path, "$", nm)]] <- list(expected = expected[[nm]][["name"]],
                                             actual = names(actual))
      }
    }

    for (nm in names(x)) {
      if (!is.list(x[[nm]]) || is.data.frame(x[[nm]])) next

      for (j in seq_along(x[[nm]])) {
        if (!inherits(x[[nm]][[j]], "bal.tab")) next

        acc <- .walk(x[[nm]][[j]], sprintf("%s$%s[[%s]]", path, nm,
                                           names(x[[nm]])[j] %or% j), acc)
      }
    }

    acc
  }

  files <- sort(list.files(dir, pattern = "\\.rds$", full.names = TRUE))
  n_tab <- 0L
  bad <- character(0L)

  for (f in files) {
    cell <- readRDS(f)

    if (isTRUE(cell[["unavailable"]]) || is_null(cell[["obj"]])) next

    tabs <- .walk(cell[["obj"]], basename(f))
    n_tab <- n_tab + length(tabs)

    for (nm in names(tabs)) {
      if (identical(tabs[[nm]][["expected"]], tabs[[nm]][["actual"]])) next

      bad <- c(bad, sprintf("%s\n  spec: %s\n  real: %s", nm,
                            toString(tabs[[nm]][["expected"]]),
                            toString(tabs[[nm]][["actual"]])))
    }
  }

  cat(sprintf("%d tables checked | %d mismatched\n", n_tab, length(bad)))

  if (is_not_null(bad)) {
    cat(paste(utils::head(bad, 5L), collapse = "\n"), "\n")
  }

  invisible(bad)
}
