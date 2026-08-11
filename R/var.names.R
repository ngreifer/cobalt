#' Extract Variable Names from `bal.tab` Objects
#' 
#' @description This function extracts variable names from a `bal.tab` object for use in specifying alternate variable names in [bal.tab()], [print()][print.bal.tab], [format()][extract.bal.tab], [as.data.frame()][extract.bal.tab], and [love.plot()]. Optionally, a file can be written for easy editing of names.
#' 
#' @param b a `bal.tab` object; the output of a call to [bal.tab()].
#' @param type the type of output desired. Can either be `"df"` for a data frame or `"vec"` for a named vector. See "Value". The default is `"vec"` unless `file` is not `NULL`.
#' @param file optional; a file name to save the output if `type = "df"`. See [utils::write.csv()], which `var.name()` calls. Must end in `.csv`.
#' @param minimal whether the output should contain all variable names (i.e., all rows that appear the output of `bal.tab()`) or just the unique base variables. See "Details".
#' 
#' @returns
#' If `type = "vec"`, a character vector with the variable names as the names and the names they are displayed under as the entries.
#'
#' If `type = "df"`, a data frame with two columns called `"old"` and `"new"`, the first with the variable names and the second with the names they are displayed under.
#'
#' When no `var.names` has been applied to the object, every variable is displayed as it is stored and the two agree. When one has, the names it produced are given, so that they can be edited rather than written out again; see Details.
#'
#' If file is not `NULL`, the output will be returned invisibly.
#' 
#' @details
#' The purpose of the function is to make supplying new variable names to the `var.names` argument easier. Rather than manually creating a vector or data frame with all the variable names that one desires to change, one can use `var.names()` to extract variable names from a `bal.tab` object and edit the output. Importantly, the output can be saved to a CSV file, which can be easily edited and read back into R, as demonstrated in the Example. See [`display-options`] for what `var.names` does and for the structures it accepts.
#'
#' When `minimal = TRUE`, only a minimal set of variables will be output. For example, if the variables analyzed in `bal.tab()` are `age`, `race`, and `married`, and `int = TRUE` in `bal.tab()`, many variables will appear in the output, including expansions of the factor variables, the polynomial terms, and the interactions. Rather than renaming all of these variables individually, one can rename just the three base variables, and all variables that arise from them will be accordingly renamed. Setting `minimal = TRUE` requests only these base variables.
#'
#' If a `var.names` was given in the call to `bal.tab()`, the names it produced are what is returned, so that a set of names arrived at once can be edited rather than written out again. Passing the result back to any of the functions that take `var.names` reproduces the names it came from, because the old names it is keyed by are the ones the variables are stored under, which `var.names` never changes.
#'
#' Which names an edit then reaches follows from `minimal`, as it does for any `var.names`: an edit to a base variable in the minimal output reaches every name that variable appears in, while an edit to an entry of the full output changes only the name it is an entry for.
#' 
#' @note
#' Not all programs can properly read the Unicode characters for the polynomial terms when requested. These may appear strange in, e.g., Excel, but R will process the characters correctly.
#'
#' @seealso
#' [`display-options`] for `var.names`, the argument this function supplies
#'
#' @examples 
#' data(lalonde, package = "cobalt")
#' 
#' b1 <- bal.tab(treat ~ age + race + married, data = lalonde,
#'               int = TRUE)
#' v1 <- var.names(b1, type = "vec", minimal = TRUE)
#' v1["age"] <- "Age (Years)"
#' v1["race"] <- "Race/Eth"
#' v1["married"] <- "Married"
#' 
#' love.plot(b1, var.names = v1, factor_sep = ": ")
#'
#' #The same names in the table, set once in `bal.tab()`
#' b1 <- bal.tab(treat ~ age + race + married, data = lalonde,
#'               int = TRUE, var.names = v1)
#' b1
#'
#' #They come back out to be edited rather than rewritten
#' v1 <- var.names(b1, type = "vec", minimal = TRUE)
#' v1
#'
#' v1["married"] <- "Married at baseline"
#' print(b1, var.names = v1)
#' \dontrun{
#' b2 <- bal.tab(treat ~ age + race + married + educ + nodegree +
#'                   re74 + re75 + I(re74==0) + I(re75==0), 
#'               data = lalonde)
#' var.names(b2, file = "varnames.csv")
#' 
#' ##Manually edit the CSV (e.g., in Excel), then save it.
#' v2 <- read.csv("varnames.csv")
#' love.plot(b2, var.names = v2)
#' }
#' 
#' @export 

var.names <- function(b, type, file = NULL, minimal = FALSE) {
  arg::arg_supplied(b)
  arg::arg_is(b, "bal.tab")
  
  co.names <- .bal.tab_co.names(b)

  if (is_null(co.names)) {
    arg::err("no variable names were found in the object. It may be a malformed object")
  }

  arg::arg_flag(minimal)

  old <- {
    if (minimal) unique(unlist(lapply(co.names, function(x) x[["component"]][x[["type"]] == "base"])))
    else names(co.names)
  }

  #The names already on display, so that a set of names given to `bal.tab()` comes back
  #to be edited rather than having to be written out again. With none given, each
  #variable is displayed as it is stored and the two columns agree, as they always did.
  new.labels <- .attr(b, "print.options")[["var.names"]]

  new <- {
    if (is_null(new.labels)) old
    #A base variable's replacement is the one it was given; a whole name's is what the
    #components it is built from come to once each has been replaced.
    else if (minimal) .renamed(old, new.labels)
    else unname(.var_names_map(co.names, new.labels)[old])
  }

  arg::when_not_null(file, arg::arg_string)
  if (is_not_null(file) && !endsWith(file, ".csv")) {
    arg::err("the filename in {.arg file} must end in {.str .csv}")
  }
  
  if (!missing(type)) {
    type <- arg::match_arg(type, c("df", "vec"))
  }
  else if (is_not_null(file)) {
    type <- "df"
  }
  else {
    type <- "vec"
  }
  
  if (type == "df") {
    out <- data.frame(old = old, new = new, row.names = NULL)
    
    if (is_not_null(file)) {
      utils::write.csv(out, file = file, row.names = FALSE)
      return(invisible(out))
    }
  }
  else {
    out <- setNames(new, old)
    
    if (is_not_null(file)) {
      arg::wrn('only {.code type = "df"} is compatible with a file name')
    }
  }

  out
}

#The several structures `var.names` accepts, reduced to the one thing they all say: a
#named character vector pairing an old name with the name to display in its place. Every
#entry point that takes `var.names` -- `bal.tab()`, `print()`, `format()`,
#`as.data.frame()`, `love.plot()` -- reads it through here, so all of them accept the same
#structures and reject the same ones.
.process_var.names <- function(var.names) {
  if (is_null(var.names)) {
    return(NULL)
  }

  if (is.data.frame(var.names)) {
    if (ncol(var.names) == 1L) {
      if (is_null(row.names(var.names))) {
        arg::err("if {.arg var.names} is a data frame with one column, its rows must be named")
      }

      return(setNames(unlist(as.character(var.names[[1L]])),
                      rownames(var.names)))
    }

    if (all(c("old", "new") %in% names(var.names))) {
      return(setNames(unlist(as.character(var.names[["new"]])), var.names[["old"]]))
    }

    if (ncol(var.names) > 2L) {
      arg::wrn("only using first 2 columns of {.arg var.names}")
    }

    return(setNames(unlist(as.character(var.names[[2L]])), var.names[[1L]]))
  }

  if (is.atomic(var.names)) {
    if (is_null(names(var.names))) {
      arg::err("if {.arg var.names} is a vector, its values must be named")
    }

    return(setNames(as.character(var.names), names(var.names)))
  }
  
  if (!is.list(var.names)) {
    arg::err("the argument to {.arg var.names} is not one of the accepted structures. See {.topic cobalt::display-options} for details")
  }
  
  if (!all_apply(var.names, function(z) is.character(z) || is.factor(z))) {
    arg::err("if {.arg var.names} is a list, its values must be the new names of the variables")
  }
  
  if (is_null(names(var.names))) {
    arg::err("if {.arg var.names} is a list, its values must be named")
  }
  
  unlist(var.names)
}

#`var.names` given to a display function adds to what was given to `bal.tab()` and
#replaces any entry it names, so that the two can be built up rather than the later one
#discarding the earlier.
.merge_var.names <- function(stored, new) {
  if (is_null(new)) {
    return(stored)
  }

  if (is_null(stored)) {
    return(new)
  }

  stored[names(new)] <- new

  stored
}

#Every covariate name the object knows about, in the parsed form `co.names` records. The
#top level carries only the first leaf's names, which is all of them except for a
#longitudinal object, whose time points need not share covariates -- so the leaves are
#asked too, the first one to name a covariate winning.
.bal.tab_co.names <- function(x) {
  #Anything else -- including a subclass's balance table, which is a bare data frame --
  #names no covariates of its own.
  if (!inherits(x, "bal.tab")) {
    return(NULL)
  }

  co.names <- .attr(x, "print.options")[["co.names"]]

  out <- co.names

  for (i in which(endsWith(names(x), ".Balance"))) {
    for (z in x[[i]]) {
      z.co.names <- .bal.tab_co.names(z)

      out <- c(out, z.co.names[setdiff(names(z.co.names), names(out))])
    }
  }

  if (is_null(out)) {
    return(NULL)
  }

  #`c()` on a list drops its attributes, and the separators are read back below.
  attr(out, "seps") <- .attr(co.names, "seps")

  out
}

#The separators the covariate names are displayed with: those the object was built with,
#each replaced by one given at display time.
.merge_seps <- function(stored, factor_sep = NULL, int_sep = NULL) {
  if (is_not_null(factor_sep)) {
    arg::arg_string(factor_sep)
    stored["factor"] <- factor_sep
  }

  if (is_not_null(int_sep)) {
    arg::arg_string(int_sep)
    stored["int"] <- int_sep
  }

  stored
}

#What each covariate is displayed as, given a set of replacements and a pair of
#separators: a named character vector from the name a covariate is stored under to the
#name to show in its place.
#
#`co.names` records what each name is made of -- which parts are variables, which are the
#separator between a factor and its level, and which the separator between the two sides
#of an interaction. So a replacement given for a base variable reaches every name it is a
#component of, and not only the name that matches it exactly; and the separators can be
#swapped for others without the names having to be parsed back apart.
#
#Matching is always against the components as stored, so a replacement is keyed by the
#name the covariate goes by in the object however that name is being displayed.
.var_names_map <- function(co.names, new.labels, seps = .attr(co.names, "seps")) {
  out <- vapply(names(co.names), function(i) {
    comp <- co.names[[i]][["component"]]
    type <- co.names[[i]][["type"]]

    if (i %in% names(new.labels) && !is.na(new.labels[i])) {
      return(unname(new.labels[i]))
    }

    if ("isep" %nin% type) {
      return(.display_name(comp, type, new.labels, seps))
    }

    #Each side of the interaction is named on its own, and the sides are rejoined with
    #the separator the display asks for.
    sep.inds <- c(which(type == "isep"), length(comp) + 1L)

    vapply(seq_along(sep.inds), function(k) {
      inds <- {
        if (k == 1L) seq(1L, sep.inds[k] - 1L)
        else seq(sep.inds[k - 1L] + 1L, sep.inds[k] - 1L)
      }

      .display_name(comp[inds], type[inds], new.labels, seps)
    }, character(1L)) |>
      paste(collapse = seps["int"])
  }, character(1L))

  #Two covariates displayed under one name would be indistinguishable in a table and
  #would sit on top of each other in a plot, so it is refused rather than shown.
  if (anyDuplicated(out) > 0L) {
    dup <- unique(out[duplicated(out)])
    arg::err("{.arg var.names} gives the name {.or {.val {dup}}} to more than one variable")
  }

  out
}

#One name, or one side of an interaction, as it is displayed: the whole of it when it was
#given a replacement, and otherwise its components, each base one replaced if it was named
#and the factor separator as the display asks for.
.display_name <- function(comp, type, new.labels, seps) {
  pasted <- paste(comp, collapse = "")

  if (pasted %in% names(new.labels) && !is.na(new.labels[pasted])) {
    return(unname(new.labels[pasted]))
  }

  named <- type == "base" & comp %in% names(new.labels) & !is.na(new.labels[comp])

  comp[named] <- new.labels[comp[named]]
  comp[type == "fsep"] <- seps["factor"]

  paste(comp, collapse = "")
}

#The names in `v` as they are to be displayed. A name with no replacement is kept.
.renamed <- function(v, map) {
  v <- as.character(v)

  new <- unname(map[match(v, names(map))])

  ifelse(is.na(new), v, new)
}

.rename_rows <- function(tab, map) {
  if (is_null(tab) || is_null(rownames(tab))) {
    return(tab)
  }

  rownames(tab) <- .renamed(rownames(tab), map)

  tab
}

#A copy of `x` in which every covariate goes by its display name. Only the names on
#display are rewritten: `co.names` is the record of what the covariates actually are, and
#is what the replacements and the separators are resolved against, so it is left as it is.
#
#The replacements are cleared from the copy's print options once applied. A segmented
#object prints each of its children by printing the child, which would otherwise resolve
#the same replacements a second time against names that have already been replaced. The
#separators need no such clearing: a child asked for the ones it was built with leaves its
#names alone.
.rename_bal.tab <- function(x, new.labels, seps = NULL) {
  co.names <- .bal.tab_co.names(x)

  seps <- seps %or% .attr(co.names, "seps")

  if (is_null(new.labels) && identical(seps, .attr(co.names, "seps"))) {
    return(x)
  }

  .rename_bal.tab_internal(x, .var_names_map(co.names, new.labels, seps))
}

.rename_bal.tab_internal <- function(x, map) {
  #A list of children: `bal.tab` objects for a cluster, imputation, treatment pair, or
  #time point, and plain balance tables for a subclass.
  for (i in which(endsWith(names(x), ".Balance"))) {
    x[[i]] <- lapply(x[[i]], function(z) {
      if (inherits(z, "bal.tab")) .rename_bal.tab_internal(z, map)
      else .rename_rows(z, map)
    })
  }

  for (i in which(names(x) == "Balance" | startsWith(names(x), "Balance.Across."))) {
    x[[i]] <- .rename_rows(x[[i]], map)
  }

  #A balance tally counts variables rather than naming them; the greatest imbalance names
  #one, in a column rather than in the row names.
  for (i in which(startsWith(names(x), "Max.Imbalance"))) {
    if (is_not_null(x[[i]][["Variable"]])) {
      x[[i]][["Variable"]] <- .renamed(x[[i]][["Variable"]], map)
    }
  }

  attr(x, "print.options")[["var.names"]] <- NULL

  x
}
