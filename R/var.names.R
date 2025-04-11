#' @title Extract Variable Names from `bal.tab` Objects
#' 
#' @description This function extracts variable names from a `bal.tab` object for use in specifying alternate variable names in [love.plot()]. Optionally, a file can be written for easy editing of names.
#' 
#' @param b a `bal.tab` object; the output of a call to [bal.tab()].
#' @param type the type of output desired. Can either be `"df"` for a data.frame or `"vec"` for a named vector. See "Value". The default is `"vec"` unless `file` is not `NULL`.
#' @param file optional; a file name to save the output if `type = "df"`. See [utils::write.csv()], which `var.name()` calls. Must end in `.csv`.
#' @param minimal whether the output should contain all variable names (i.e., all rows that appear the output of `bal.tab()`) or just the unique base variables. See "Details".
#' 
#' @returns If `type = "vec"`, a character vector the the variable names both as the names and the entries.
#' 
#' If `type = "df"`, a `data.frame` with two columns called `"old"` and `"new"`, each with the variables as the entries.
#' 
#' If file is not `NULL`, the output will be returned invisibly.
#' 
#' @details The goal of the function is to make supplying new variable names to the `var.names` argument in [love.plot()] easier. Rather than manually creating a vector or `data.frame` with all the variable names that one desires to change, one can use `var.names()` to extract variable names from a `bal.tab` object and edit the output. Importantly, the output can be saved to a CSV file, which can be easily edited and read back into R for use in `love.plot()`, as demonstrated in the Example.
#'
#' When `minimal = TRUE`, only a minimal set of variables will be output. For example, if the variables analyzed in `bal.tab()` are `age`, `race`, and `married`, and `int = TRUE` in `bal.tab()`, many variables will appear in the output, including expansions of the factor variables, the polynomial terms, and the interactions. Rather than renaming all of these variables individually, one can rename just the three base variables, and all variables that arise from them will be accordingly renamed. Setting `minimal = TRUE` requests only these base variables.
#' 
#' @note Not all programs can properly read the Unicode characters for the polynomial terms when requested. These may appear strange in, e.g., Excel, but R will process the characters correctly.
#' 
#' @examples 
#' 
#' data(lalonde, package = "cobalt")
#' 
#' b1 <- bal.tab(treat ~ age + race + married, data = lalonde,
#'               int = TRUE)
#' v1 <- var.names(b1, type = "vec", minimal = TRUE)
#' v1["age"] <- "Age (Years)"
#' v1["race"] <- "Race/Eth"
#' v1["married"] <- "Married"
#' love.plot(b1, var.names = v1)
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
  if (is_null(attr(b, "print.options")[["co.names"]])) {
    .err("no variable names were found in the object. It is probably not a bal.tab object")
  }
  
  .chk_flag(minimal)
  vars <- {
    if (minimal) {
      unique(unlist(lapply(attr(b, "print.options")[["co.names"]],
                           function(x) x[["component"]][x[["type"]] == "base"])))
    }
    else {
      vapply(attr(b, "print.options")[["co.names"]],
             function(x) paste(x[["component"]], collapse = ""),
             character(1L))
    }
  }
  
  .chk_null_or(file, .chk_string)
  if (is_not_null(file) && !endsWith(file, ".csv")) {
    .err("the filename in `file` must end in \".csv\"")
  }
  
  if (!missing(type)) {
    .chk_string(type)
    type <- tolower(type)
    type <- match_arg(type, c("df", "vec"))
  }
  else if (is_not_null(file)) {
    type <- "df"
  }
  else {
    type <- "vec"
  }
  
  out <- switch(type,
                "df" = data.frame(old = vars, new = vars),
                "vec" = setNames(vars, vars))
  
  if (is_null(file)) {
    return(out)
  }
  
  if (type == "df") {
    utils::write.csv(out, file = file, row.names = FALSE)
    return(invisible(out))
  }
  
  .wrn('only `type = "df"` is compatible with a file name')
  out
}
