#' @title Set and Get Options in `cobalt`
#' 
#' @description Makes it easier to set \pkg{cobalt} options. `set.cobalt.options()` is essentially a wrapper for [options()] but performs several checks, and `get.cobalt.options()` is essentially a wrapper for [getOption()].
#' 
#' @param ... For `set.cobalt.options()`, `bal.tab()` parameters and the values they should take. These should be the name of the parameter in `bal.tab()` without `"cobalt_"` preceding them. See examples. If any values are `NULL`, the corresponding options will be set back to their defaults.
#' 
#' For `get.cobalt.options()`, one or more strings containing the name of a parameter option to be retrieved. See examples. If empty, all available options and their values will be returned.
#' 
#' @param default if `TRUE`, sets all \pkg{cobalt} options not named in `...` to their default values.
#' 
#' @details When an option is set to `NULL`, it is set to its default value. The defaults are not displayed but are listed on the help pages where they appear. Most options correspond to display options, which can be accessed [here][display-options]. Some others (e.g., `continuous` and `binary`) are described on the [bal.tab()] help page.
#' 
#' @seealso 
#' * [options()]
#' * [`display-options`] for some arguments that can be set via options.
#' 
#' @examples 
#' # Set un to be TRUE to always display unadjusted 
#' # balance measures and set binary to "std" to 
#' # produce standardized mean differences for 
#' # binary variables.
#' 
#' set.cobalt.options(un = TRUE, binary = "std")
#' 
#' # Note: the above is equivalent to:
#' # options(cobalt_un = TRUE, cobalt_binary = "std")
#' # but performs some additional checks
#' 
#' get.cobalt.options("un", "binary")
#' 
#' # Note: the above is equivalent to:
#' # getOption("cobalt_un")
#' # getOption("cobalt_binary")
#' 
#' # Return all cobalt options to their defaults
#' 
#' set.cobalt.options(default = TRUE)
#' 
#' # View all available options
#' get.cobalt.options()
#' 

#' @rdname set.cobalt.options
#' @export
set.cobalt.options <- function(..., default = FALSE) {
  
  if (...length() > 0L && (is_null(...names()) || "" %in% ...names())) {
    .err("all arguments must be named")
  }
  # if ("continuous" %in% names(opts)) names(opts)[names(opts) == "continuous"] <- "cont"
  # if ("binary" %in% names(opts)) names(opts)[names(opts) == "binary"] <- "bin"
  
  multiple.allowed <- c("stats", "disp", "cluster.fun", "imp.fun")
  any.string.allowed <- c("int_sep", "factor_sep")
  
  duplicates <- table(...names()) > 1
  
  if (any(duplicates)) {
    .err("{.arg {names(duplicates)[duplicates]}} {?is/are} present more than once in the input to {.fun set.cobalt.options}")
  }
  
  if (!all(...names() %in% names(acceptable.options()))) {
    .wrn("the following are not acceptable options and will be ignored: {.arg {setdiff(...names(), names(acceptable.options()))}}")
  }
  
  opts <- ...mget(names(acceptable.options()))
  
  return.to.default <- NULL
  if (default) {
    return.to.default <- setdiff(names(acceptable.options()), names(opts))
  }
  
  multiple.opts <- NULL
  bad.opts <- NULL
  for (i in names(opts)) {
    if (is_null(opts[[i]])) {
      return.to.default <- c(return.to.default, i)
      opts[[i]] <- NULL
    }
    else {
      if (length(opts[[i]]) > 1L && i %nin% multiple.allowed) {
        multiple.opts <- c(multiple.opts, i)
      }
      
      if (mode(opts[[i]]) != mode(acceptable.options()[[i]]) || 
          (!(is.character(opts[[i]]) && is.character(acceptable.options()[[i]]) &&
             (i %in% any.string.allowed || !anyNA(pmatch(opts[[i]], acceptable.options()[[i]])))) &&
           !all(opts[[i]] %in% acceptable.options()[[i]]))) {
        bad.opts <- c(bad.opts, i)
      }
    }
  }
  
  if (is_not_null(opts)) {
    both.opts <- intersect(multiple.opts, bad.opts)
    multiple.opts <- setdiff(multiple.opts, both.opts)
    bad.opts <- setdiff(bad.opts, both.opts)
    
    problematic.opts <- make_list(c("multiple", "bad", "both"))
    
    problematic.opts[["multiple"]] <- setNames(lapply(multiple.opts, function(i) {
      cli::format_inline("{.arg {i}} must be of length 1")
    }), multiple.opts)
    
    problematic.opts[["bad"]] <- setNames(lapply(bad.opts, function(i) {
      if (i %in% any.string.allowed) cli::format_inline("{.arg {i}} must be a character string")
      else cli::format_inline("{.arg {i}} must be {.or {.val {acceptable.options()[[i]]}}}")
    }), bad.opts)
    
    problematic.opts[["both"]] <- setNames(lapply(both.opts, function(i) {
      if (i %in% any.string.allowed) cli::format_inline("{.arg {i}} must be a character string of length 1")
      else cli::format_inline("{.arg {i}} must be one of {.or {.val {acceptable.options()[[i]]}}}")
    }), both.opts)
    
    problems <- do.call("c", unname(problematic.opts))
    problems <- problems[names(opts)[names(opts) %in% names(problems)]]
    
    if (is_not_null(problems)) {
      .err(do.call("paste", c(list(""), problems, list("\nNo options will be set", sep = "\n"))), cli = FALSE)
    }
    
    names(opts) <- paste0("cobalt_", names(opts))
    options(opts)
  }
  
  if (is_not_null(return.to.default)) {
    options(setNames(replicate(length(return.to.default), NULL),
                     paste0("cobalt_", return.to.default)))
  }
  # if ("continuous" %in% names(opts)) names(acceptable.options)[names(acceptable.options) == "continuous"] <- "cont"
  # if ("binary" %in% names(opts)) names(acceptable.options)[names(acceptable.options) == "binary"] <- "bin"
}

#' @rdname set.cobalt.options
#' @export
get.cobalt.options <- function(...) {
  
  if (...length() == 0L) {
    opts <- names(acceptable.options())
  }
  else {
    opts <- character(...length())
    for (i in seq_len(...length())) {
      if (!is.character(...elt(i))) {
        .err("all arguments must be strings containing the name of an option to return")
      }
      opts[i] <- ...elt(i)
    }
    
    not.in.accept <- opts %nin% names(acceptable.options())
    if (any(not.in.accept)) {
      .err("{.val {opts[not.in.accept]}} {?is/are} not {?an/} acceptible option{?s}")
    }
  }
  
  paste0("cobalt_", opts) |>
    lapply(getOption) |>
    setNames(opts)
}

#set.cobalt.options
acceptable.options <- function() {
  TF <- c(TRUE, FALSE)
  list(stats = c("mean.diffs"),
       un = TF,
       continuous = c("raw", "std"),
       binary = c("raw", "std"),
       imbalanced.only = TF,
       disp = c("means", "sds"),
       disp.means = TF,
       disp.sds = TF,
       disp.v.ratio = TF,
       disp.ks = TF,
       disp.subclass = TF,
       disp.bal.tab = TF,
       cluster.summary = TF,
       cluster.fun = c("min", "mean", "max"),
       imp.summary = TF,
       imp.fun = c("min", "mean", "max"),
       multi.summary = TF,
       msm.summary = TF,
       target.summary = TF,
       subclass.summary = TF,
       int_sep = " * ",
       factor_sep = "_",
       center = TF,
       orth = TF,
       remove_perfect_col = TF,
       disp.call = TF)
}
