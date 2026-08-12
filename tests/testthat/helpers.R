#Test helpers for cobalt.
#
#Messages raised by `arg::err()`/`arg::wrn()`/`arg::msg()` are formatted by cli,
#which capitalizes the first letter, appends a period, converts inline markup
#(e.g., `{.arg x}` to `` `x` ``), and hard-wraps the result to the console width.
#The wrapping means a literal `fixed = TRUE` match against a long message fails,
#and building a regex from the message requires escaping every metacharacter that
#cli may have introduced (`[`, `]`, `{`, `}`, `+`, `?`, `*`, `|`).
#
#`expect_err()`, `expect_wrn()`, and `expect_msg()` sidestep both problems: they
#collapse all whitespace in the *observed* message and then do a literal substring
#match. Never copy a source string from `R/` into these -- use the rendered text.

#Collapse runs of whitespace so matching is invariant to how cli wrapped the message.
squish <- function(x) {
  gsub("\\s+", " ", trimws(paste(x, collapse = " ")))
}

#Shared back end for the three expectations below.
.expect_cnd_text <- function(cnd, text, what) {
  if (is.null(cnd)) {
    testthat::fail(sprintf("Expected %s, but none was signaled.", what))
    return(invisible(NULL))
  }

  if (is.null(text)) {
    testthat::succeed()
    return(invisible(cnd))
  }

  msg <- squish(conditionMessage(cnd))

  #Matching is case-insensitive because cli capitalizes the first letter of every
  #message, which would otherwise make any substring starting at the beginning of
  #the message fail. The point is to identify which message was raised, not to
  #pin its capitalization.
  testthat::expect_true(
    grepl(tolower(text), tolower(msg), fixed = TRUE),
    info = sprintf("%s message did not contain the expected text.\n  expected: %s\n  actual:   %s",
                   what, encodeString(text, quote = "\""), encodeString(msg, quote = "\""))
  )

  invisible(cnd)
}

#Each of the three expectations below asserts that one condition was signaled.
#Anything else the expression prints or signals is incidental to that assertion, and
#letting it through only buries the reporter's own output, so the other two condition
#classes are muffled and printed output is discarded.
#
#That is a real amount of output rather than a stray line or two: `expect_wrn()` and
#`expect_msg()` run the expression to completion, so `expect_wrn(print(b, ...))` --
#the shape most of the `print()` argument checks take -- wrote a whole balance table
#to the console on every call.
#
#`capture.output()` evaluates its argument in the calling frame, so the promise is
#still forced there and assignments inside the expression still take effect there.

#Expect an error whose message contains `text` (literal, whitespace-insensitive).
expect_err <- function(object, text = NULL) {
  cnd <- NULL

  utils::capture.output(
    withCallingHandlers(cnd <- tryCatch({
      force(object)
      NULL
    }, error = function(e) e),
    message = function(m) invokeRestart("muffleMessage"),
    warning = function(w) invokeRestart("muffleWarning")))

  .expect_cnd_text(cnd, text, "error")
}

#Expect a warning whose message contains `text`. The expression still runs to
#completion, so assignments inside it take effect in the calling environment.
expect_wrn <- function(object, text = NULL) {
  cnd <- NULL

  utils::capture.output(
    withCallingHandlers(force(object),
                        message = function(m) invokeRestart("muffleMessage"),
                        warning = function(w) {
                          if (is.null(cnd)) cnd <<- w
                          invokeRestart("muffleWarning")
                        }))

  .expect_cnd_text(cnd, text, "warning")
}

#Expect a message whose message contains `text`. As with `expect_wrn()`, the
#expression runs to completion.
expect_msg <- function(object, text = NULL) {
  cnd <- NULL

  utils::capture.output(
    withCallingHandlers(force(object),
                        warning = function(w) invokeRestart("muffleWarning"),
                        message = function(m) {
                          if (is.null(cnd)) cnd <<- m
                          invokeRestart("muffleMessage")
                        }))

  .expect_cnd_text(cnd, text, "message")
}

#Send plot output to a null device for the duration of the calling test. Needed
#for `print.love.plot()`, `plot.bal.tab()`, and `autoplot.bal.tab()`, which draw
#to the active device; without this they create a stray `Rplots.pdf`.
local_null_device <- function(.env = parent.frame()) {
  grDevices::pdf(NULL)

  #`rlang::defer()` is not exported, so register the cleanup as an `on.exit()`
  #expression in the calling frame instead. This needs no extra dependency.
  do.call(base::on.exit,
          list(quote(grDevices::dev.off()), add = TRUE, after = FALSE),
          envir = .env)

  invisible(NULL)
}

#Clear every `cobalt_*` option for the duration of the calling test, and restore
#the previous state afterwards.
#
#`set.cobalt.options()` writes to global `options()`, so without this a test that
#sets one leaks it into every later test file in the run -- and restoring with
#`options(old)` is not enough, because an option that was previously *unset* cannot
#be un-set that way. `rlang::local_options()` handles the unset case correctly.
local_cobalt_options <- function(.env = parent.frame()) {
  cobalt_opts <- paste0("cobalt_", names(cobalt:::acceptable.options()))
  opts <- setNames(vector("list", length(cobalt_opts)), cobalt_opts)

  do.call(rlang::local_options, c(opts, list(.frame = .env)))

  invisible(NULL)
}

#As above, and additionally pin the base options that affect printed numbers.
local_bal_tab_snapshot <- function(.env = parent.frame()) {
  local_cobalt_options(.env = .env)

  rlang::local_options(digits = 7L, OutDec = ".", useFancyQuotes = FALSE,
                       .frame = .env)

  invisible(NULL)
}

#Gate for tests that are slow or depend on a fragile upstream package. They run
#whenever `COBALT_SLOW_TESTS=true`, and are skipped on CRAN otherwise.
skip_if_slow <- function() {
  if (!identical(Sys.getenv("COBALT_SLOW_TESTS"), "true")) {
    testthat::skip_on_cran()
  }

  invisible(NULL)
}
