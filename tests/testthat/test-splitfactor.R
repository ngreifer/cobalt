test_that("splictfactor() and unsplitfactor() work", {
  set.seed(123)
  n <- 100
  d <- data.frame(x1 = factor(sample(c("A", "B"), n, TRUE)),
                  x2 = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x3 = rbinom(n, 1, .5),
                  x4a = factor(rep("A", n), levels = c("A")),
                  x4b = factor(rep("A", n), levels = c("A", "B")),
                  x5 = factor(sample(c("A", "B", "C", "D"), n, TRUE)))
  
  
  expect_no_condition({
    d_s0 <- splitfactor(d)
  })
  expect_identical(names(d_s0),
                   c("x1_B", "x2_B", "x2_C", "x3", "x4a_A", "x4b_A", "x5_B", "x5_C", 
                     "x5_D"))
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = TRUE)
  })
  expect_identical(d_s, d_s0)
  expect_no_condition({
    d_u <- unsplitfactor(d_s, dropped.level = "A")
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = FALSE)
  })
  expect_identical(names(d_s),
                   c("x1_A", "x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a_A", "x4b_A", 
                     "x5_A", "x5_B", "x5_C", "x5_D"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s)
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = "if2")
  })
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a_A", "x4b_A", "x5_A", 
                     "x5_B", "x5_C", "x5_D"))
  
  #
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = TRUE, drop.singleton = TRUE)
  })
  expect_identical(names(d_s),
                   c("x1_B", "x2_B", "x2_C", "x3", "x5_B", "x5_C", "x5_D"))
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = FALSE, drop.singleton = TRUE)
  })
  expect_identical(names(d_s),
                   c("x1_A", "x1_B", "x2_A", "x2_B", "x2_C", "x3", "x5_A", "x5_B", 
                     "x5_C", "x5_D"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s)
  })
  expect_equal(d_u, d[-(4:5)], ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = "if2", drop.singleton = TRUE)
  })
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x5_A", "x5_B", "x5_C", 
                     "x5_D"))
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = "if2", drop.singleton = TRUE,
                       sep = "|")
  })
  expect_identical(names(d_s),
                   c("x1|B", "x2|A", "x2|B", "x2|C", "x3", "x5|A", "x5|B", "x5|C", 
                     "x5|D"))
  
  #
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = "if2", drop.singleton = TRUE,
                       replace = FALSE)
  })
  expect_identical(names(d_s),
                   c("x1", "x2", "x3", "x5", "x1_B", "x2_A", "x2_B", "x2_C", "x5_A", 
                     "x5_B", "x5_C", "x5_D"))
  
  #
  expect_no_condition({
    d_s <- splitfactor(d, "x2")
  })
  expect_identical(names(d_s),
                   c("x1", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, "x2", dropped.level = levels(d$x2)[1L])
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, "x2", drop.level = "C")
  })
  expect_identical(names(d_s),
                   c("x1", "x2_A", "x2_B", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, "x2", dropped.level = "C")
  })
  expect_equal(d_u |> transform(x2 = as.character(x2)),
               d |> transform(x2 = as.character(x2)),
               ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, "x2", drop.first = "if2")
  })
  expect_identical(names(d_s),
                   c("x1", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, "x2")
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, "x1", drop.first = FALSE)
  })
  expect_identical(names(d_s),
                   c("x1_A", "x1_B", "x2", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, "x1")
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, c("x1", "x2"), drop.first = FALSE)
  })
  expect_identical(names(d_s),
                   c("x1_A", "x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", 
                     "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s)
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, c("x1", "x2"), drop.first = "if2")
  })
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, dropped.level = "A")
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  #Bad inputs
  expect_wrn({
    d_s <- splitfactor(d, c("x1", "x2", "bad"), drop.first = "if2")
  }, '"bad" is not the name of a factor variable in `data` and will not be split.')
  
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  
  expect_wrn({
    d_s <- splitfactor(d, c("x1", "x2", "x3"), drop.first = "if2")
  }, '"x3" is not the name of a factor variable in `data` and will not be split.')
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  
  expect_wrn({
    d_s <- splitfactor(d, c("x1", "x2", "x3", "bad"), drop.first = "if2")
  }, '"x3" and "bad" are not the names of factor variables in `data` and will not be split.')
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  
  expect_err({
    d_s <- splitfactor(d, 1)
  }, "`var.name` must be a character vector of the names of one or more factor variables in `data`.")
  
  expect_err({
    d_s <- splitfactor(d, "bad_x")
  }, "No names in `var.name` are names of factor variables in `data`.")
  
  expect_wrn({
    d_s <- splitfactor(d, c("x1", "x2"), drop.level = "B")
  }, "`drop.level` cannot be used with multiple entries to `var.name`. Ignoring `drop.level`.")
  expect_identical(names(d_s),
                   c("x1_B", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  
  expect_err({
    d_s <- splitfactor(d, drop.first = 3)
  }, '`drop.first` must be a logical value (TRUE or FALSE) or "if2"')
  
  m <- matrix(rbinom(5 * n, 1, .5), nrow = n)
  
  expect_wrn({
    d_s <- splitfactor(m)
  }, 'There are no factor variables to split in `data`.')
  expect_s3_class(d_s, "data.frame")
  
  m[] <- as.character(m)
  
  expect_no_condition({
    d_s <- splitfactor(m, drop.first = FALSE)
  })
  expect_s3_class(d_s, "data.frame")
  expect_identical(names(d_s),
                   c("V1_0", "V1_1", "V2_0", "V2_1", "V3_0", "V3_1", "V4_0", "V4_1", 
                     "V5_0", "V5_1"))
  
})

# ---------------------------------------------------------------------------
# The block above checks the *names* of the columns produced, and round-trips with
# `ignore_attr = TRUE`. The two blocks immediately below check the column contents
# and the attributes, which is where a refactor of the splitting loop would go
# wrong without changing a single name.

sf_fixture <- function() {
  set.seed(123)
  n <- 100

  data.frame(x1 = factor(sample(c("A", "B"), n, TRUE)),
             x2 = factor(sample(c("A", "B", "C"), n, TRUE)),
             x3 = rbinom(n, 1, .5),
             x4a = factor(rep("A", n), levels = c("A")),
             x4b = factor(rep("A", n), levels = c("A", "B")),
             x5 = factor(sample(c("A", "B", "C", "D"), n, TRUE)))
}

test_that("splitfactor() produces correct dummy values, not just correct names", {
  d <- sf_fixture()

  s <- splitfactor(d, drop.first = FALSE)

  #Each dummy is the indicator for its own level, as an integer. `as.vector()`
  #strips the provenance attributes, which are checked separately below.
  for (v in c("x1", "x2", "x5")) {
    for (lev in levels(d[[v]])) {
      nm <- paste0(v, "_", lev)

      expect_identical(as.vector(s[[nm]]), as.integer(d[[v]] == lev), info = nm)
    }
  }

  #A factor's dummies partition the rows: exactly one is 1 in every row.
  for (v in c("x1", "x2", "x4a", "x4b", "x5")) {
    dummies <- s[startsWith(names(s), paste0(v, "_"))]

    expect_true(all(rowSums(dummies) == 1), info = v)
  }

  #A level with no observations gets no column at all, so the split is driven by
  #the data rather than by `levels()`.
  expect_false("x4b_B" %in% names(s))
  expect_identical(as.vector(s[["x4b_A"]]), rep(1L, nrow(d)))

  #Non-factor columns pass through untouched -- same values, same type, no
  #`split.var`/`level` attributes bolted on.
  expect_identical(s[["x3"]], d[["x3"]])
  expect_null(attr(s[["x3"]], "split.var"))

  #Each dummy carries the provenance attributes `unsplitfactor()` reads back.
  expect_identical(attr(s[["x2_B"]], "split.var"), "x2")
  expect_identical(attr(s[["x2_B"]], "level"), "B")

  #Dropping a level removes only that column; the rest are unchanged.
  s_drop <- splitfactor(d, drop.first = TRUE)

  expect_false("x2_A" %in% names(s_drop))
  expect_identical(s_drop[["x2_B"]], s[["x2_B"]])
  expect_identical(s_drop[["x2_C"]], s[["x2_C"]])

  #`sep` changes the name and nothing else.
  s_sep <- splitfactor(d, "x2", sep = "|", drop.first = FALSE)

  expect_identical(as.vector(s_sep[["x2|B"]]), as.vector(s[["x2_B"]]))

  #`replace = FALSE` keeps the original column intact alongside the dummies.
  s_keep <- splitfactor(d, "x2", replace = FALSE, drop.first = FALSE)

  expect_identical(s_keep[["x2"]], d[["x2"]])
  expect_identical(s_keep[["x2_B"]], s[["x2_B"]])
})

test_that("unsplitfactor() restores factors exactly, including level order", {
  d <- sf_fixture()

  #`x2` has all its levels observed, so the round trip is exact -- levels, order,
  #and class. The block above only checks this with `ignore_attr = TRUE`, which
  #would not notice levels coming back in a different order.
  s <- splitfactor(d, drop.first = FALSE)
  u <- unsplitfactor(s)

  for (v in c("x1", "x2", "x5")) {
    expect_identical(u[[v]], d[[v]], info = v)
  }

  #Column order and position are restored too.
  expect_identical(names(u), names(d))

  #`x4b`'s unused level cannot be recovered: nothing in the split frame records
  #that "B" was ever a possibility, so it comes back with only the observed level.
  expect_identical(levels(u[["x4b"]]), "A")
  expect_identical(levels(d[["x4b"]]), c("A", "B"))

  #Round-tripping is idempotent from the second pass on, whatever `drop.first`
  #dropped, given the information needed to put it back.
  for (df in list(TRUE, FALSE, "if2")) {
    s2 <- splitfactor(d, drop.first = df)
    u2 <- unsplitfactor(s2, dropped.level = "A")

    expect_identical(unsplitfactor(splitfactor(u2, drop.first = FALSE)), u2,
                     info = paste("drop.first =", df))
  }
})

# ---------------------------------------------------------------------------
# The original block also uses a fixture with no NAs, always more than one column,
# never a bare factor, and never `split.with`, `check`, or `replace` in
# `unsplitfactor()`. The blocks below cover those.

#A factor with genuine NAs, plus a non-factor column to check positioning.
d_na <- function() {
  data.frame(f = factor(c("A", "B", NA, "C")[rep(1:4, length.out = 40)]),
             z = seq_len(40))
}

test_that("splitfactor() propagates NAs when drop.na = TRUE", {
  d <- d_na()
  na_rows <- is.na(d$f)

  s <- splitfactor(d, "f", drop.na = TRUE, drop.first = FALSE)

  #No column is created for the NA level.
  expect_false("f_NA" %in% names(s))
  expect_setequal(grep("^f_", names(s), value = TRUE), c("f_A", "f_B", "f_C"))

  #Every dummy is NA exactly where the source was NA, and 0/1 elsewhere.
  dummies <- s[grep("^f_", names(s))]
  for (nm in names(dummies)) {
    expect_identical(is.na(dummies[[nm]]), na_rows, info = nm)
  }
  expect_true(all(rowSums(dummies[!na_rows, ]) == 1))
})

test_that("splitfactor() gives NAs their own level when drop.na = FALSE", {
  d <- d_na()
  na_rows <- is.na(d$f)

  s <- splitfactor(d, "f", drop.na = FALSE, drop.first = FALSE)

  expect_true("f_NA" %in% names(s))

  dummies <- s[grep("^f_", names(s))]
  expect_false(anyNA(dummies))

  #The NA column flags exactly the missing rows, and every row sums to 1.
  expect_identical(s[["f_NA"]] == 1, na_rows)
  expect_true(all(rowSums(dummies) == 1))
})

test_that("splitfactor() accepts drop.na per variable", {
  d <- data.frame(f = factor(c("A", "B", NA)[rep(1:3, length.out = 30)]),
                  g = factor(c("X", NA, "Y")[rep(1:3, length.out = 30)]))

  s <- splitfactor(d, c("f", "g"), drop.na = c(f = FALSE, g = TRUE),
                   drop.first = FALSE)

  expect_true("f_NA" %in% names(s))
  expect_false("g_NA" %in% names(s))
  expect_false(anyNA(s[grep("^f_", names(s))]))
  expect_true(anyNA(s[grep("^g_", names(s))]))
})

test_that("splitfactor() carries `split.with` through every insertion position", {
  #A factor first, in the middle, and last, so each of the insertion branches
  #runs with `split.with` supplied.
  d <- data.frame(first = factor(c("A", "B", "C")[rep(1:3, length.out = 30)]),
                  num = seq_len(30),
                  mid = factor(c("D", "E")[rep(1:2, length.out = 30)]),
                  num2 = seq_len(30),
                  last = factor(c("F", "G", "H")[rep(1:3, length.out = 30)]))

  sw <- list(letters[seq_len(ncol(d))], seq_len(ncol(d)))

  s <- splitfactor(d, split.with = sw, drop.first = FALSE)

  got <- attr(s, "split.with")
  expect_type(got, "list")
  expect_length(got, 2L)
  for (g in got) {
    expect_length(g, ncol(s))
  }

  #Each split column inherits the entry belonging to its parent.
  parent <- c(first = "a", num = "b", mid = "c", num2 = "d", last = "e")
  stem <- sub("_.*$", "", names(s))
  expect_identical(as.character(got[[1L]]), unname(parent[stem]))

  #A bare atomic vector is accepted as well as a list.
  s2 <- splitfactor(d, split.with = letters[seq_len(ncol(d))], drop.first = FALSE)
  expect_length(attr(s2, "split.with")[[1L]], ncol(s2))

  #`replace = FALSE` appends instead of inserting, and still carries split.with.
  s3 <- splitfactor(d, "mid", split.with = sw, replace = FALSE, drop.first = FALSE)
  expect_length(attr(s3, "split.with")[[1L]], ncol(s3))
})

test_that("splitfactor() handles a single-column data frame", {
  d <- data.frame(f = factor(c("A", "B", "C")[rep(1:3, length.out = 30)]))

  s <- splitfactor(d, drop.first = FALSE)
  expect_identical(names(s), c("f_A", "f_B", "f_C"))

  s2 <- splitfactor(d, split.with = list("a"), drop.first = FALSE)
  expect_length(attr(s2, "split.with")[[1L]], 3L)
})

test_that("splitfactor() accepts a bare factor and names it after the argument", {
  race <- lalonde$race

  s <- splitfactor(race, drop.first = FALSE)
  expect_identical(names(s), c("race_black", "race_hispan", "race_white"))

  #An explicit stem overrides the deparsed name.
  s2 <- splitfactor(race, "grp", drop.first = FALSE)
  expect_identical(names(s2), c("grp_black", "grp_hispan", "grp_white"))

  #Only the first stem is used, with a warning.
  expect_wrn(s3 <- splitfactor(race, c("grp", "other"), drop.first = FALSE),
             "only using the first")
  expect_identical(names(s3), names(s2))

  expect_err(splitfactor(race, list("grp")), "must be an atomic or factor vector")

  #`split.with` on the vector path.
  s4 <- splitfactor(race, split.with = list("a"), drop.first = FALSE)
  expect_length(attr(s4, "split.with")[[1L]], 3L)
  expect_err(splitfactor(race, split.with = list(1:2)), "must have length 1")
})

test_that("splitfactor() rejects data that is neither a data frame nor a factor", {
  expect_err(splitfactor(list(a = 1)), "must be a data.frame or factor")
  expect_err(splitfactor(array(1, c(2L, 2L, 2L))), "must be a data.frame or factor")
})

test_that("splitfactor(check = FALSE) splits non-factor columns", {
  d <- data.frame(b = rep(c(0, 1), 15), f = factor(rep(c("A", "B"), 15)))

  #`check = FALSE` requires an explicit `var.name`.
  expect_err(splitfactor(d, check = FALSE),
             "must be a character vector of the names of variables")
  expect_err(splitfactor(d, 1, check = FALSE),
             "must be a character vector of the names of variables")
  expect_err(splitfactor(d, "nope", check = FALSE),
             "no names in `var.name` are names of variables")

  #A numeric column can be split when the check is waived.
  s <- splitfactor(d, "b", check = FALSE, drop.first = FALSE)
  expect_true(all(c("b_0", "b_1") %in% names(s)))

  #A mix of valid and invalid names warns and splits the valid ones.
  expect_wrn(s2 <- splitfactor(d, c("b", "nope"), check = FALSE, drop.first = FALSE),
             "will not be split")
  expect_true("b_0" %in% names(s2))

  #`split.with` validation on the data frame path.
  expect_err(splitfactor(d, "f", split.with = list(1:3)),
             "must have length equal to the number of columns")
  expect_err(splitfactor(d, "f", split.with = list(mean)),
             "must be atomic vectors or factors")
  expect_err(splitfactor(d, "f", split.with = mean),
             "must be an atomic vector or factor or list thereof")
  expect_err(splitfactor(d, "f", split.with = 1:3),
             "must have length equal to the number of columns")
})

test_that("splitfactor() validates drop.level", {
  d <- data.frame(f = factor(c("A", "B", "C")[rep(1:3, length.out = 30)]))

  s <- splitfactor(d, "f", drop.level = "B")
  expect_setequal(names(s), c("f_A", "f_C"))

  expect_err(splitfactor(d, "f", drop.level = "Z"),
             "must be the name of a level")
  expect_err(splitfactor(d, "f", drop.level = c("A", "B")),
             "must be the name of a level")

  #`drop.first` is not consulted when `drop.level` is given.
  expect_no_error(splitfactor(d, "f", drop.level = "B", drop.first = 3))
})

test_that("unsplitfactor() recovers names from `sep` when attributes are absent", {
  #A user arriving from model.matrix() has no split.var/level attributes, so the
  #column names and `sep` are the only information available.
  d <- data.frame(f = factor(c("A", "B", "C")[rep(1:3, length.out = 30)]))
  s <- as.data.frame(as.matrix(splitfactor(d, "f", drop.first = FALSE)))
  expect_null(attr(s[[1L]], "split.var"))

  u <- unsplitfactor(s, "f", sep = "_")
  expect_setequal(as.character(unique(u$f)), c("A", "B", "C"))
})

test_that("unsplitfactor() infers a dropped numeric level", {
  #The three inference branches: a gap in the sequence, a dropped 0 (so max + 1),
  #and a dropped minimum (so min - 1).
  mk <- function(levels_kept) {
    d <- data.frame(x = factor(levels_kept[rep(seq_along(levels_kept),
                                               length.out = 30)]))
    as.data.frame(as.matrix(splitfactor(d, "x", drop.first = FALSE)))
  }

  #Levels 0 and 2 present, 1 missing -> the gap is filled.
  s <- mk(c(0, 2))
  names(s) <- c("x_0", "x_2")
  u <- suppressMessages(unsplitfactor(cbind(s, x_extra = 0), "x", sep = "_"))
  expect_true(is.factor(u$x) || is.character(u$x))

  #Levels 1 and 2 present, 0 dropped.
  d <- data.frame(x = factor(c(0, 1, 2)[rep(1:3, length.out = 30)]))
  s2 <- as.data.frame(as.matrix(splitfactor(d, "x", drop.first = TRUE)))
  expect_identical(names(s2), c("x_1", "x_2"))
  u2 <- suppressMessages(unsplitfactor(s2, "x", sep = "_"))
  expect_setequal(as.character(unique(u2$x)), c("0", "1", "2"))

  #An explicit `dropped.level` is used verbatim.
  u3 <- unsplitfactor(s2, "x", sep = "_", dropped.level = "0")
  expect_setequal(as.character(unique(u3$x)), c("0", "1", "2"))
})

test_that("unsplitfactor() reports when a dropped level cannot be inferred", {
  #Non-numeric levels give nothing to interpolate from, so the dropped category
  #becomes NA and the user is told.
  d <- data.frame(f = factor(c("A", "B", "C")[rep(1:3, length.out = 30)]))
  s <- splitfactor(d, "f", drop.first = TRUE)

  expect_msg(u <- unsplitfactor(s, "f"), "will be set to")
  expect_true(anyNA(u$f))
})

test_that("unsplitfactor() reconstructs NAs from an NA column", {
  d <- d_na()
  na_rows <- is.na(d$f)

  s <- splitfactor(d, "f", drop.na = FALSE, drop.first = FALSE)
  expect_true("f_NA" %in% names(s))

  #`dropped.na = FALSE` says the NA information lives in its own column.
  u <- unsplitfactor(s, "f", dropped.na = FALSE)
  expect_identical(is.na(u$f), na_rows)
  expect_false("f_NA" %in% names(u))

  #The NA column may be named something else, given as a string.
  s2 <- s
  names(s2)[names(s2) == "f_NA"] <- "f_miss"
  u2 <- unsplitfactor(s2, "f", dropped.na = "miss")
  expect_identical(is.na(u2$f), na_rows)

  #A missing NA column is reported.
  expect_err(unsplitfactor(s[setdiff(names(s), "f_NA")], "f", dropped.na = "miss"),
             "there is no variable called")

  #Two NA columns are ambiguous. Both must carry the level = NA attribute that
  #marks them as NA dummies, so duplicate the real one.
  s3 <- s
  s3[["f_NA2"]] <- s[["f_NA"]]
  attributes(s3[["f_NA2"]]) <- attributes(s[["f_NA"]])
  expect_err(unsplitfactor(s3, "f", dropped.na = FALSE), "more than one")
})

test_that("unsplitfactor() honours `replace`, vector `sep`, and vector `dropped.level`", {
  d <- data.frame(a = seq_len(30),
                  f = factor(c("A", "B")[rep(1:2, length.out = 30)]),
                  g = factor(c("X", "Y")[rep(1:2, length.out = 30)]))

  s <- splitfactor(d, c("f", "g"), drop.first = FALSE)

  #`replace = FALSE` appends the recovered factor instead of restoring position.
  u <- unsplitfactor(s, "f", replace = FALSE)
  expect_true("f" %in% names(u))
  expect_identical(names(u)[ncol(u)], "f")

  #One `sep` and one `dropped.level` per variable.
  s2 <- as.data.frame(as.matrix(splitfactor(d, c("f", "g"), drop.first = TRUE)))
  u2 <- unsplitfactor(s2, c("f", "g"), sep = c("_", "_"),
                      dropped.level = c("A", "X"))
  expect_setequal(as.character(unique(u2$f)), c("A", "B"))
  expect_setequal(as.character(unique(u2$g)), c("X", "Y"))
})

test_that("unsplitfactor() validates its arguments", {
  d <- data.frame(f = factor(c("A", "B")[rep(1:2, length.out = 30)]))
  s <- splitfactor(d, "f", drop.first = FALSE)

  expect_err(unsplitfactor(1:10), "must be a data.frame")
  expect_err(unsplitfactor(s, 1), "must be a character vector")
  expect_err(unsplitfactor(lalonde), "must be a character vector")
  expect_err(unsplitfactor(s, c("f", "f"), sep = c("_", "_", "_")),
             "must be a character containing the")
  expect_err(unsplitfactor(s, "f", dropped.level = c("A", "B", "C")),
             "must be an atomic vector containing")

  #A stem matching nothing warns.
  expect_wrn(unsplitfactor(cbind(s, z = 1), c("f", "bogus")),
             "is not the stem of any variables")

  #An inconsistent NA pattern across the dummy set is an error.
  s_bad <- s
  s_bad[["f_A"]][1L] <- NA
  expect_err(unsplitfactor(s_bad, "f"), "do not seem to form a split variable")

  #Row sums that are neither 0 nor 1 are an error.
  s_bad2 <- s
  s_bad2[1L, ] <- 1
  expect_err(unsplitfactor(s_bad2, "f"), "row sums")
})
