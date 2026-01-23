# Functions for clean argument checking. These are based on {chk} but
# with several modifications to ensure correct and clean printing. These revolve
# around pkg_caller_call(), which tells chk::err() which *user-facing* function
# the error occurred in. These functions require {cli}.

# pkg_caller_call() searches along the call stack to find the function from this
# package that the user called. It gets a list of functions and methods exported
# by this package, moves along the call stack, and returns the call of the
# highest level function that is in this package.
pkg_caller_call <- function() {
  pn <- utils::packageName()
  package.funs <- c(getNamespaceExports(pn),
                    .getNamespaceInfo(asNamespace(pn), "S3methods")[, 3L])
  
  for (i in seq_len(sys.nframe())) {
    e <- sys.call(i)
    
    n <- rlang::call_name(e)
    
    if (is_not_null(n) && n %in% package.funs) {
      return(e)
    }
  }
  
  NULL
}

# .err() is a version of chk::err() that uses pkg_caller_call() to get the correct function
# call since chk::err() has a default that doesn't always work. .wrn() and .msg()
# just call chk::wrn() and chk::msg() but make the syntax consistent.
.err <- function(m, n = NULL, tidy = TRUE, cli = TRUE) {
  if (cli) {
    m <- eval.parent(substitute(cli::format_inline(.m), list(.m = m)))
  }
  
  chk::message_chk(m, n = n, tidy = tidy) |>
    cli::ansi_strwrap() |>
    paste(collapse = "\n") |>
    rlang::abort(call = pkg_caller_call())
}
.wrn <- function(m, n = NULL, tidy = TRUE, immediate = TRUE, cli = TRUE) {
  if (cli) {
    m <- eval.parent(substitute(cli::format_inline(.m), list(.m = m)))
  }
  
  m <- chk::message_chk(m, n = n, tidy = tidy)
  
  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
    rlang::with_options({
      m |>
        cli::ansi_strwrap() |>
        paste(collapse = "\n") |>
        rlang::warn()
    }, warn = 1)
  }
  else {
    m |>
      cli::ansi_strwrap() |>
      paste(collapse = "\n") |>
      rlang::warn()
  }
}
.msg <- function(m, n = NULL, tidy = TRUE, cli = TRUE) {
  if (cli) {
    m <- eval.parent(substitute(cli::format_inline(.m), list(.m = m)))
  }
  
  chk::message_chk(m, n = n, tidy = tidy) |>
    cli::ansi_strwrap() |>
    paste(collapse = "\n") |>
    rlang::inform(tidy = FALSE)
}

# Kind of insane loop to create (at build time) version of all .chk_*
# functions that use .err() instead of chk::abort_chk() internally. All
# .chk_* function now have a version like .chk_*, e.g., .chk_flag(),
# that can be used in package code instead of the chk version.
for (i in getNamespaceExports("chk")) {
  if (!startsWith(i, "chk_")) next
  assign(paste0(".", i), eval(str2expression(sprintf(
    "function(...) {
            tryCatch(chk::%s(...),
                     error = function(e) .err(conditionMessage(e), cli = FALSE))
        }", i
  ))))
}

# Version of .chk_null_or() that isn't bugged.
.chk_null_or <- function(x, chk, ..., x_name = NULL) {
  if (is.null(x_name)) {
    x_name <- deparse1(substitute(x))
  }
  
  x_name <- add_quotes(x_name, "`")
  
  if (is.null(x)) {
    return(invisible(x))
  }
  
  tryCatch(chk(x, ..., x_name = x_name),
           error = function(e) {
             .err(sub("[.]$", " or `NULL`.", conditionMessage(e)), cli = FALSE)
           })
}

.chk_formula <- function(x, sides = NULL, x_name = NULL) {
  if (is.null(sides)) {
    if (rlang::is_formula(x)) {
      return(invisible(x))
    }
    if (is.null(x_name)) {
      x_name <- deparse1(substitute(x))
    }
    .err("{.arg {x_name}} must be a formula")
  }
  
  if (sides == 1) {
    if (rlang::is_formula(x, lhs = FALSE)) {
      return(invisible(x))
    }
    
    if (is.null(x_name)) {
      x_name <- deparse1(substitute(x))
    }
    .err("{.arg {x_name}} must be a formula with no left-hand side")
  }
  
  if (sides == 2) {
    if (rlang::is_formula(x, lhs = TRUE)) {
      return(invisible(x))
    }
    if (is.null(x_name)) {
      x_name <- deparse1(substitute(x))
    }
    .err("{.arg {x_name}} must be a formula with a left-hand side")
  }
  
  stop("`sides` must be NULL, 1, or 2")
}

try_chk <- function(expr, tidy = FALSE, warn = TRUE) {
  .e <- function(e) {
    .err(conditionMessage(e), tidy = tidy, cli = FALSE)
  }
  
  .w <- function(w) {
    .wrn(conditionMessage(w), tidy = tidy, cli = FALSE)
  }
  
  if (warn) {
    tryCatch(expr, error = .e, warning = .w)
  }
  else {
    tryCatch(expr, error = .e)
  }
}
