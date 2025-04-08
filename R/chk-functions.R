# Functions for clean argument checking. These are based on {chk} but
# with several modifications to ensure correct and clean printing. These revolve
# around pkg_caller_call(), which tells chk::err() which *user-facing* function
# the error occurred in.

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
    
    if (is_null(n)) {
      next
    }
    
    if (n %in% package.funs) {
      return(e)
    }
  }
  
  NULL
}

# .err() is a version of chk::err() that uses pkg_caller_call() to get the correct function
# call since chk::err() has a default that doesn't always work. .wrn() and .msg()
# just call chk::wrn() and chk::msg() but make the syntax consistent.
.err <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::abort(paste(strwrap(m), collapse = "\n"),
               call = pkg_caller_call())
}
.wrn <- function(..., n = NULL, tidy = TRUE, immediate = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  
  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
    rlang::with_options({
      rlang::warn(paste(strwrap(m), collapse = "\n"))
    }, warn = 1)
  }
  else {
    rlang::warn(paste(strwrap(m), collapse = "\n"))
  }
}
.msg <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::inform(paste(strwrap(m), collapse = "\n"), tidy = FALSE)
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
                     error = function(e) .err(conditionMessage(e)))
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
             msg <- sub("[.]$", " or `NULL`.",
                        conditionMessage(e))
             .err(msg, .subclass = "chk_error")
           })
}

.chk_formula <- function(x, sides = NULL, x_name = NULL) {
  if (is.null(sides)) {
    if (rlang::is_formula(x)) {
      return(invisible(x))
    }
    if (is.null(x_name)) {
      x_name <- chk::deparse_backtick_chk(substitute(x))
    }
    .err(x_name, " must be a formula",
         x = x)
  }
  else if (sides == 1) {
    if (rlang::is_formula(x, lhs = FALSE)) {
      return(invisible(x))
    }
    if (is.null(x_name)) {
      x_name <- chk::deparse_backtick_chk(substitute(x))
    }
    .err(x_name, " must be a formula with no left-hand side",
         x = x)
  }
  else if (sides == 2) {
    if (rlang::is_formula(x, lhs = TRUE)) {
      return(invisible(x))
    }
    if (is.null(x_name)) {
      x_name <- chk::deparse_backtick_chk(substitute(x))
    }
    .err(x_name, " must be a formula with a left-hand side",
         x = x)
  }
  else stop("`sides` must be NULL, 1, or 2")
}

try_chk <- function(expr, tidy = FALSE, warn = TRUE) {
  .e <- function(e) {
    .err(conditionMessage(e), tidy = tidy)
  }
  
  .w <- function(w) {
    .wrn(conditionMessage(w), tidy = tidy)
  }
  
  if (warn) {
    tryCatch(expr, error = .e, warning = .w)
  }
  else {
    tryCatch(expr, error = .e)
  }
}