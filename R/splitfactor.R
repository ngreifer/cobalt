#' Split and Unsplit Factors into Dummy Variables
#' 
#' @description `splitfactor()` splits factor variables into dummy (0/1) variables. This can be useful when functions do not process factor variables well or require numeric matrices to operate. `unsplitfactor()` combines dummy variables into factor variables, undoing the operation of `splitfactor()`.
#' 
#' @param data A `data.frame` containing the variables to be split or unsplit. In `splitfactor()`, can be a factor variable to be split.
#' @param var.name For `splitfactor()`, the names of the factor variables to split. If not specified, will split all factor variables in `data`. If `data` is a factor, the stem for each of the new variables to be created. For `unsplitfactor()`, the name of the previously split factor. If not specified and `data` is the output of a call to `splitfactor()`, all previously split variables will be unsplit.
#' @param drop.level The name of a level of `var.name` for which to drop the dummy variable. Only works if there is only one variable to be split.
#' @param drop.first Whether to drop the first dummy created for each factor. If `"if2"`, will only drop the first category if the factor has exactly two levels. The default is to always drop the first dummy (`TRUE`).
#' @param drop.singleton Whether to drop a factor variable if it only has one level.
#' @param drop.na If `NA`s are present in the variable, how to handle them. If `TRUE`, no new dummy will be created for `NA` values, but all created dummies will have `NA` where the original variable was `NA`. If `FALSE`, `NA` will be treated like any other factor level, given its own column, and the other dummies will have a value of 0 where the original variable is `NA`.
#' @param sep A character separating the the stem from the value of the variable for each dummy. For example, for `"race_black"`, `sep = "_"`.
#' @param replace Whether to replace the original variable(s) with the new variable(s) (`TRUE`) or the append the newly created variable(s) to the end of the data set (`FALSE`).
#' @param split.with A list of vectors or factors with lengths equal to the number of columns of `data` that are to be split in the same way `data` is. See Details.
#' @param check Whether to make sure the variables specified in `var.name` are actually factor (or character) variables. If splitting non-factor (or non-character) variables into dummies, set `check = FALSE`. If `check = FALSE` and `data` is a `data.frame`, an argument to `var.name` must be specified.
#' @param dropped.level The value of each original factor variable whose dummy was dropped when the variable was split. If left empty and a dummy was dropped, the resulting factor will have the value `NA` instead of the dropped value. There should be one entry per variable to unsplit. If no dummy was dropped for a variable, an entry is still required, but it will be ignored.
#' @param dropped.na If `TRUE`, will assume that `NA`s in the variables to be unsplit correspond to `NA` in the unsplit factor (i.e., that `drop.na = TRUE` was specified in `split.factor()`). If `FALSE`, will assume there is a dummy called "var.name_stem_NA" (e.g., "x_NA") that contains 1s where the unsplit factor should be `NA` (i.e., that `drop.na = FALSE` was specified in `split.factor()`. If `NA`s are stored in a different column with the same stem, e.g., "x_miss", that name (e.g., "miss") can be entered instead.
#' 
#' @returns
#' For `splitfactor()`, a `data.frame` containing the original data set with the newly created dummies. For `unsplitfactor()`. a `data.frame` containing the original data set with the newly created factor variables.
#' 
#' @details
#' If there are `NA`s in the variable to be split, the new variables created by `splitfactor()` will have `NA` where the original variable is `NA`.
#' 
#' When using `unsplitfactor()` on a `data.frame` that was generated with `splitfactor()`, the arguments `dropped.na`, and `sep` are unnecessary.
#' 
#' If `split.with` is supplied, the elements will be split in the same way `data` is. For example, if `data` contained a 4-level factor that was to be split, the entries of `split.with` at the same index as the factor and would be duplicated so that resulting entries will have the same length as the number of columns of `data` after being split. The resulting values are stored in the `"split.with"` attribute of the output object. See Examples.
#' 
#' @seealso 
#' [model.matrix()]
#' 
#' @examples 
#' data("lalonde", package = "cobalt")
#' 
#' lalonde.split <- splitfactor(lalonde, "race", 
#'                              replace = TRUE, 
#'                              drop.first = TRUE)
#' # A data set with "race_hispan" and "race_white" instead
#' # of "race".
#' 
#' lalonde.unsplit <- unsplitfactor(lalonde.split, "race", 
#'                                  replace = TRUE,
#'                                  dropped.level = "black")
#' 
#' all.equal(lalonde, lalonde.unsplit) #TRUE
#' 
#' # Demonstrating the use of split.with:
#' to.split <- list(letters[1:ncol(lalonde)],
#'                  1:ncol(lalonde))
#' 
#' lalonde.split <- splitfactor(lalonde, split.with = to.split,
#'                              drop.first = FALSE)
#' attr(lalonde.split, "split.with")
#' 
#' 

#' @export 
splitfactor <- function(data, var.name, drop.level = NULL, drop.first = TRUE,
                        drop.singleton = FALSE, drop.na = TRUE, sep = "_",
                        replace = TRUE, split.with = NULL, check = TRUE) {
  #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
  #Retains all categories unless only 2 levels, in which case only the second level is retained.
  #If variable only has one level, will delete.
  #var.name= the name of the variable to split when data is specified
  #data=data set to be changed
  
  if (is.matrix(data) || length(dim(data)) == 2L) {
    data <- as.data.frame(data)
  }
  
  if (is.data.frame(data)) {
    data <- as.data.frame(data)
    if (check) {
      factor.names <- names(data)[vapply(data, chk::vld_character_or_factor, logical(1L))]
      if (is_null(factor.names)) {
        .wrn("there are no factor variables to split in {.arg data}")
        return(data)
      }
      
      if (missing(var.name)) {
        var.name <- factor.names
      }
      else {
        if (!is.character(var.name)) {
          .err("{.arg var.name} must be a character vector of the names of one or more factor variables in {.arg data}")
        }
        
        if (!any(var.name %in% factor.names)) {
          .err("no names in {.arg var.name} are names of factor variables in {.arg data}")
        }
        
        if (!all(var.name %in% factor.names)) {
          not.in.factor.names <- setdiff(var.name, factor.names)
          .wrn("{.val {not.in.factor.names}} {?is/are} not the name{?s} of {?a/} factor variable{?s} in {.arg data} and will not be split")
        }
        var.name <- intersect(var.name, factor.names)
      }
      
    }
    else {
      if (missing(var.name) || !is.character(var.name)) {
        .err("{.arg var.name} must be a character vector of the names of variables in {.arg data}")
      }
      
      if (!any(var.name %in% names(data))) {
        .err("no names in {.arg var.name} are names of variables in {.arg data}")
      }
      
      if (!all(var.name %in% names(data))) {
        not.in.data.names <- setdiff(var.name, names(data))
        .wrn("{.val {not.in.factor.names}} {?is/are} not the name{?s} of {?a/} factor variable{?s} in {.arg data} and will not be split")
      }
      
      var.name <- intersect(var.name, names(data))
    }
    
    if (is_not_null(split.with)) {
      if (is.list(split.with)) {
        if (!all_apply(split.with, is.atomic)) {
          .err("all entries in {.arg split.with} must must be atomic vectors or factors")
        }
        
        if (!all(lengths(split.with) == ncol(data))) {
          .err("all entries in {.arg split.with} must have length equal to the number of columns of {.arg data}")
        }
      }
      else {
        if (!is.atomic(split.with)) {
          .err("{.arg split.with} must must be an atomic vector or factor or list thereof")
        }
        
        if (length(split.with) != ncol(data)) {
          .err("{.arg split.with} must have length equal to the number of columns of {.arg data}")
        }
        
        split.with <- list(split.with)
      }
    }
  }
  else if (is.atomic(data) && length(dim(data)) <= 1L) {
    dep <- deparse1(substitute(data))
    data <- data.frame(data)
    if (missing(var.name) || is_null(var.name)) {
      names(data) <- dep
    }
    else if (is.atomic(var.name)) {
      if (length(var.name) == 1L) {
        names(data) <- var.name
      }
      else {
        .wrn("only using the first item of {.arg var.name}")
        names(data) <- var.name[1L]
      }
    }
    else {
      .err("{.arg var.name} must be an atomic or factor vector of length 1 with the stem of the new variable")
    }
    var.name <- names(data)
    
    if (is_not_null(split.with)) {
      if (is.list(split.with)) {
        if (!all_apply(split.with, is.atomic)) {
          .err("all entries in {.arg split.with} must must be atomic vectors or factors")
        }
        
        if (!all(lengths(split.with) == ncol(data))) {
          .err("all entries in {.arg split.with} must have length 1")
        }
      }
      else {
        if (!is.atomic(split.with)) {
          .err("{.arg split.with} must must be an atomic vector or factor or list thereof")
        }
        if (length(split.with) != ncol(data)) {
          .err("{.arg split.with} must have length 1")
        }
        split.with <- list(split.with)
      }
    }
  }
  else {
    .err("{.arg data} must a be a data.frame or factor")
  }
  
  if (is_not_null(drop.level) && length(var.name) > 1L) {
    .wrn("{.arg drop.level} cannot be used with multiple entries to {.arg var.name}. Ignoring {.arg drop.level}")
    drop.level <- NULL
  }
  
  .chk_flag(drop.singleton)
  drop.na <- rlang::rep_named(var.name, drop.na)
  
  for (v in var.name) {
    x <- factor(data[v][[1L]])
    anyNAx <- anyNA(x)
    na.level <- if (anyNAx) nlevels(x) + 1L else 0L
    new.levels <- c(levels(x), if (anyNAx) NA_character_)
    
    if (length(new.levels) == 1L && drop.singleton) {
      data[[v]] <- NULL
      next
    }
    
    k <- make_df(paste(v, new.levels, sep = sep), nrow(data),
                 types = "integer")
    
    if (anyNAx) {
      for (i in seq_col(k)) {
        if (i == na.level) k[[i]][is.na(x)] <- 1L
        else k[[i]][!is.na(x) & x == new.levels[i]] <- 1L
      }
    }
    else {
      for (i in seq_col(k)) {
        k[[i]][x == new.levels[i]] <- 1L
      }
    }
    
    if (!anyNAx) {
      drop.na[v] <- FALSE
    }
    else if (drop.na[v]) {
      is.na(k)[is.na(x), -na.level] <- TRUE
    }
    
    dropl <- rlang::rep_named(new.levels, FALSE)
    if (is_not_null(drop.level)) {
      if (!is.character(drop.level) || length(drop.level) != 1L || drop.level %nin% new.levels) {
        .err("{.arg drop} must be the name of a level of {.var {v}} that is to be dropped")
      }
      dropl[drop.level] <- TRUE
    }
    else if (identical(drop.first, "if2") || chk::vld_flag(drop.first)) {
      if ((ncol(k) == 2L && (isTRUE(drop.first) || drop.first == "if2")) ||
          (ncol(k) > 2L && isTRUE(drop.first))) {
        dropl[1L] <- TRUE
      }
    }
    else {
      .err('{.arg drop.first} must be {.val {TRUE}}, {.val {FALSE}}, or {.val {"if2"}}')
    }
    
    if (drop.na[v]) {
      dropl[na.level] <- TRUE
    }
    
    for (i in seq_col(k)) {
      attr(k[[i]], "split.var") <- v
      attr(k[[i]], "level") <- new.levels[i]
    }
    
    k[dropl] <- NULL
    
    for (i in seq_along(split.with)) {
      names(split.with[[i]]) <- names(data)
    }
    
    if (ncol(data) == 1L) {
      data <- k
      for (i in seq_along(split.with)) {
        split.with[[i]] <- rep.int(split.with[[i]][v], ncol(k))
        names(split.with[[i]]) <- names(data)
      }
    }
    else {
      .chk_flag(replace)
      
      if (!replace) {
        data <- setNames(data.frame(data, k, row.names = rownames(data)),
                         c(names(data), names(k)))
        
        for (i in seq_along(split.with)) {
          split.with[[i]] <- c(split.with[[i]], rep(split.with[[i]][v], ncol(k)))
          names(split.with[[i]]) <- names(data)
        }
      }
      else if (names(data)[1L] == v) {
        data <- setNames(data.frame(k, data[names(data) != v], row.names = rownames(data)),
                         c(names(k), names(data)[names(data) != v]))
        
        for (i in seq_along(split.with)) {
          split.with[[i]] <- c(rep(split.with[[i]][v], ncol(k)),
                               split.with[[i]][names(split.with[[i]]) != v])
          names(split.with[[i]]) <- names(data)
        }
      }
      else if (names(data)[ncol(data)] == v) {
        data <- setNames(data.frame(data[names(data) != v], k, row.names = rownames(data)),
                         c(names(data)[names(data) != v], names(k)))
        
        for (i in seq_along(split.with)) {
          split.with[[i]] <- c(split.with[[i]][names(split.with[[i]]) != v],
                               rep(split.with[[i]][v], ncol(k)))
          names(split.with[[i]]) <- names(data)
        }
      }
      else {
        where <- match(v, names(data))
        data <- setNames(data.frame(data[1:(where - 1L)], k, data[(where + 1L):ncol(data)], row.names = rownames(data)),
                         c(names(data)[1:(where - 1L)], names(k), names(data)[(where + 1L):ncol(data)]))
        
        for (i in seq_along(split.with)) {
          split.with[[i]] <- c(split.with[[i]][1:(where - 1L)],
                               rep(split.with[[i]][v], ncol(k)),
                               split.with[[i]][(where + 1L):length(split.with[[i]])])
          names(split.with[[i]]) <- names(data)
        }
      }
    }
  }
  
  attr(data, "split.with") <- split.with
  
  data
}

#' @rdname splitfactor
#' @export
unsplitfactor <- function(data, var.name, dropped.level = NULL, dropped.na = TRUE,
                          sep = "_", replace = TRUE) {
  
  if (!is.data.frame(data)) {
    .err("{.arg data} must be a data.frame containing the variables to unsplit")
  }
  
  if (missing(var.name)) {
    split.dummies <- vapply(data, function(x) {
      is_not_null(.attr(x, "split.var")) && is_not_null(.attr(x, "level")) &&
        all(x %in% c(0L, 1L, NA_integer_))
    }, logical(1L))
    
    if (!any(split.dummies)) {
      .err("{.arg var.name} must be a character vector containing the names of the variables to unsplit")
    }
    
    var.name <- unique(vapply(data[split.dummies], .attr, character(1L), "split.var"))
  }
  else if (!is.character(var.name)) {
    .err("{.arg var.name} must be a character vector containing the names of the variables to unsplit")
  }
  
  if (is_null(sep)) {
    sep <- rlang::rep_named(var.name, "")
  }
  else {
    if (length(sep) %nin% c(1L, length(var.name)) || !is.atomic(sep)) {
      .err("{.arg sep} must be a character containing the seperating character in the names of the split variables. See {.fun unsplitfactor} for details")
    }
    
    if (length(sep) == 1L) {
      sep <- rlang::rep_named(var.name, sep)
    }
    else {
      names(sep) <- var.name
    }
  }
  
  if (is_null(dropped.level)) {
    dropped.level <- NULL
  }
  else {
    if (length(dropped.level) %nin% c(1L, length(var.name)) || !is.atomic(dropped.level)) {
      .err("{.arg dropped.level} must be an atomic vector containing the value of the dropped category of each split variable. See {.fun unsplitfactor} for details")
    }
    
    if (length(dropped.level) == 1L) {
      dropped.level <- rlang::rep_named(var.name, dropped.level)
    }
    else {
      names(dropped.level) <- var.name
    }
  }
  
  not.the.stem <- character()
  
  for (v in var.name) {
    v_sep <- paste0(v, sep[v])
    dropped.level0 <- dropped.level[v]
    split.data <- vapply(data, function(x) identical(.attr(x, "split.var"), v), logical(1L))
    v.is.split <- any(split.data)
    
    var.to.combine <- {
      if (v.is.split) data[split.data]
      else data[startsWith(names(data), v_sep)]
    }
    
    if (is_null(var.to.combine)) {
      not.the.stem <- c(not.the.stem, v_sep)
      next
    }
    
    if (!all(rowSums(apply(var.to.combine, 2L, is.na)) %in% c(0L, ncol(var.to.combine)))) {
      .err("the variables in {.arg data} selected based on {.arg var.name} and {.arg sep} do not seem to form a split variable based on the <NA> pattern")
    }
    
    NA.column <- character()
    
    na.dummy <- vapply(var.to.combine, function(x) {
      s <- .attr(x, "level")
      is_not_null(s) && is.na(s)
    }, logical(1L))
    
    if (v.is.split && any(na.dummy)) {
      dropped.na <- FALSE
    }
    
    .chk_flag(dropped.na)
    if (!dropped.na) {
      NA.column <- {
        if (v.is.split)
          names(var.to.combine)[na.dummy]
        else
          paste0(v, sep[v], {if (isFALSE(dropped.na)) "NA" else dropped.na})
      }
      
      if (length(NA.column) > 1L) {
        .err("there appears to be more than one {.val {NA}} variable for {.var {v}}")
      }
      
      if (!utils::hasName(var.to.combine, NA.column)) {
        .err("there is no variable called {.val {NA.column}} to generate the {.val {NA}} values")
      }
      
      is.na(var.to.combine[var.to.combine[[NA.column]] == 1L, ]) <- TRUE
      var.to.combine[[NA.column]] <- NULL
    }
    
    var.sum <- rowSums(var.to.combine)
    
    if (!all(is.na(var.sum) | check_if_zero(var.sum - 1))) {
      if (!all(is.na(var.sum) | check_if_zero(var.sum - 1) | check_if_zero(var.sum - 0))) {
        .err("the variables in {.arg data} selected based on {.arg var.name} and {.arg sep} do not seem to form a split variable based on the row sums")
      }
      
      #Missing category
      if (is_null(dropped.level[v])) {
        
        k.levels0 <- {
          if (v.is.split) unlist(lapply(var.to.combine, .attr, "level"))
          else substr(names(var.to.combine), nchar(v_sep), nchar(names(var.to.combine)))
        }
        
        if (can_str2num(k.levels0)) {
          k.levels0.num <- as.numeric(k.levels0)
          nums <- seq(min(k.levels0.num), max(k.levels0.num))
          
          dropped.level0 <- {
            if (!all(nums %in% k.levels0.num)) setdiff(nums, k.levels0.num)[1L]
            else if (min(as.numeric(k.levels0)) == 0) max(k.levels0.num) + 1
            else min(k.levels0.num) - 1
          }
          
          dropped.name <- paste0(v_sep, dropped.level0)
        }
        else {
          .msg("the dropped category for {.var {v}} will be set to {.val {NA}}")
          dropped.name <- dropped.level0 <- NA_character_
        }
        
      }
      else {
        dropped.name <- paste0(v_sep, dropped.level[v])
      }
      
      var.to.combine <- setNames(data.frame(1 - var.sum, var.to.combine),
                                 c(dropped.name, names(var.to.combine)))
      
      if (v.is.split) {
        attr(var.to.combine[[1L]], "split.var") <- v
        attr(var.to.combine[[1L]], "level") <- as.character(dropped.level0)
      }
    }
    
    .chk_flag(replace)
    
    k.levels <- {
      if (v.is.split) unlist(lapply(var.to.combine, .attr, "level")) 
      else substr(names(var.to.combine),
                  1L + nchar(v_sep),
                  nchar(names(var.to.combine)))
    }
    
    k <- rep.int(NA_character_, nrow(data))
    for (i in seq_along(k.levels)) {
      k[var.to.combine[[i]] == 1] <- k.levels[i]
    }
    
    k <- factor(k, levels = k.levels)
    
    if (replace) {
      where <- which(names(data) %in% c(names(var.to.combine), NA.column))
      
      data[[where[1L]]] <- k
      remove.cols <- where[-1L]
      if (is_not_null(remove.cols)) data[remove.cols] <- NULL
      names(data)[where[1L]] <- v
    }
    else {
      data[[v]] <- k
    }
  }
  
  if (is_not_null(not.the.stem)) {
    .wrn("{.val {not.the.stem}} {?is/are} not the stem{?s} of any variables in {.arg data} and will be ignored. Ensure {.arg var.name} and {.arg sep} are correct")
  }
  
  data
}
