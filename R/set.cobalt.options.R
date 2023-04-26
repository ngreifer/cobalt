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
    opts <- list(...)
    if (is_not_null(opts) && (is_null(names(opts)) ||  "" %in% names(opts))) {
        .err("all arguments must be named")
    }
    # if ("continuous" %in% names(opts)) names(opts)[names(opts) == "continuous"] <- "cont"
    # if ("binary" %in% names(opts)) names(opts)[names(opts) == "binary"] <- "bin"
    
    multiple.allowed <- c("stats", "disp", "cluster.fun", "imp.fun")
    any.string.allowed <- c("int_sep", "factor_sep")
    
    if (any(duplicates <- table(names(opts)) > 1)) {
        .err(sprintf("%s present more than once in the input to `set.cobalt.options()`",
                     word_list(names(duplicates)[duplicates], is.are = TRUE, quotes = "`")))
    }
    
    if (any(names(opts) %nin% names(acceptable.options()))) {
        .wrn(sprintf("the following are not acceptable options and will be ignored: %s",
                     word_list(unique(names(opts)[names(opts) %nin% names(acceptable.options())]),
                               quotes = "`")))
        opts <- opts[names(opts) %in% names(acceptable.options())]
    }
    
    if (default) {
        return.to.default <- setdiff(names(acceptable.options()), names(opts))
    }
    else return.to.default <- NULL
    
    multiple.opts <- NULL
    bad.opts <- NULL
    for (i in names(opts)) {
        if (is_null(opts[[i]])) {
            return.to.default <- c(return.to.default, i)
            opts[[i]] <- NULL
        }
        else {
            if (length(opts[[i]]) > 1 && i %nin% multiple.allowed) multiple.opts <- c(multiple.opts, i)
            if (mode(opts[[i]]) != mode(acceptable.options()[[i]]) || 
                (!(is.character(opts[[i]]) && is.character(acceptable.options()[[i]]) && (i %in% any.string.allowed || !anyNA(pmatch(opts[[i]], acceptable.options()[[i]])))) &&
                 !all(opts[[i]] %in% acceptable.options()[[i]]))) bad.opts <- c(bad.opts, i)
        }
    }
    
    if (is_not_null(opts)) {
        both.opts <- intersect(multiple.opts, bad.opts)
        multiple.opts <- multiple.opts[multiple.opts %nin% both.opts]
        bad.opts <- bad.opts[bad.opts %nin% both.opts]
        problematic.opts <- make_list(c("multiple", "bad", "both"))
        problematic.opts[["multiple"]] <- setNames(lapply(multiple.opts, function(i) {
            paste(add_quotes(i, "`"), "must be of length 1.")
        }), multiple.opts)
        problematic.opts[["bad"]] <- setNames(lapply(bad.opts, function(i) {
            if (i %in% any.string.allowed) paste0(add_quotes(i, "`"), " must be a character string.")
            else paste0(add_quotes(i, "`"), " must be ", word_list(acceptable.options()[[i]], quotes = 2*is.character(acceptable.options()[[i]]), and.or = "or"), ".")
        }), bad.opts)
        problematic.opts[["both"]] <- setNames(lapply(both.opts, function(i) {
            if (i %in% any.string.allowed) paste0(add_quotes(i, "`"), " must be a character string of length 1.")
            else paste0(add_quotes(i, "`"), " must be one of ", word_list(acceptable.options()[[i]], quotes = 2*is.character(acceptable.options()[[i]]), and.or = "or"), ".")
        }), both.opts)
        
        problems <- do.call("c", unname(problematic.opts))
        problems <- problems[names(opts)[names(opts) %in% names(problems)]]
        if (is_not_null(problems)) {
            .err(do.call("paste", c(list(""), problems, list("\nNo options will be set", sep = "\n"))))
        }
        
        names(opts) <- paste0("cobalt_", names(opts))
        options(opts)
    }
    
    if (is_not_null(return.to.default)) {
        options(setNames(replicate(length(return.to.default), NULL), paste0("cobalt_", return.to.default)))
    }
    # if ("continuous" %in% names(opts)) names(acceptable.options)[names(acceptable.options) == "continuous"] <- "cont"
    # if ("binary" %in% names(opts)) names(acceptable.options)[names(acceptable.options) == "binary"] <- "bin"
}

#' @rdname set.cobalt.options
#' @export
get.cobalt.options <- function(...) {
    opts <- list(...)
    
    opts <- clear_null(opts)
    if (is_null(opts)) opts <- names(acceptable.options())
    else {
        if (!all(vapply(opts, is.character, logical(1L)))) {
            .err("all arguments must be strings containing the name of an option to return")
        }
        opts <- do.call("c", opts)
        if (any(not.in.accept <- opts %nin% names(acceptable.options()))) {
            .err(sprintf("%s not %s",
                         word_list(opts[not.in.accept], is.are = TRUE, quotes = 2),
                         ngettext(sum(not.in.accept), "an acceptable option", "acceptible options")))
        }
    }
    
    setNames(lapply(paste0("cobalt_", opts), getOption), opts)
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
         remove_perfect_col = TF,
         disp.call = TF)
}
