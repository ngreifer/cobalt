#' @title Split and Unsplit Factors into Dummy Variables
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
#' @returns For `splitfactor()`, a `data.frame` containing the original data set with the newly created dummies. For `unsplitfactor()`. a `data.frame` containing the original data set with the newly created factor variables.
#' 
#' @details If there are `NA`s in the variable to be split, the new variables created by `splitfactor()` will have `NA` where the original variable is `NA`.
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
#' @rdname splitfactor
#' @export 
splitfactor <- function(data, var.name, drop.level = NULL, drop.first = TRUE, drop.singleton = FALSE, drop.na = TRUE, sep = "_", replace = TRUE, split.with = NULL, check = TRUE) {
    #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
    #Retains all categories unless only 2 levels, in which case only the second level is retained.
    #If variable only has one level, will delete.
    #var.name= the name of the variable to split when data is specified
    #data=data set to be changed
    
    if (is.data.frame(data)) {
        data <- as.data.frame(data)
        if (check) {
            factor.names <- names(data)[vapply(data, chk::vld_character_or_factor, logical(1L))]
            if (is_null(factor.names)) {
                .wrn("there are no factor variables to split in `data`")
                return(data)
            }
            
            if (missing(var.name)) {
                var.name <- factor.names
            }
            else if (!is.character(var.name)) {
                .err("`var.name` must be a character vector of the names of one or more factor variables in `data`")
            }
            else if (!any(var.name %in% factor.names)) {
                .err("no names in `var.name` are names of factor variables in `data`")
            }
            else {
                if (any(var.name %nin% factor.names)) {
                    not.in.factor.names <- var.name[var.name %nin% factor.names]
                    .wrn(sprintf("%s not the %s in `data` and will not be split",
                                 word_list(not.in.factor.names, "and", is.are = TRUE, quotes = 2),
                                 ngettext(length(not.in.factor.names),
                                          "name of a factor variable",
                                          "names of factor variables")))
                }
                var.name <- var.name[var.name %in% factor.names]
            }
            
        }
        else {
            if (missing(var.name) || !is.character(var.name)) {
                .err("`var.name` must be a character vector of the names of variables in `data`")
            }
            
            if (!any(var.name %in% names(data))) {
                .err("no names in `var.name` are names of variables in `data`")
            }
            
            if (any(var.name %nin% names(data))) {
                not.in.data.names <- var.name[!var.name %in% names(data)]
                .wrn(sprintf("%s not the %s in `data` and will not be split",
                             word_list(not.in.data.names, "and", is.are = TRUE, quotes = 2),
                             ngettext(length(not.in.data.names),
                                      "name of a variable",
                                      "names of variables")))
            }
            var.name <- var.name[var.name %in% names(data)]
            
        }
        
        if (is_not_null(split.with)) {
            if (is.list(split.with)) {
                if (!all(vapply(split.with, is.atomic, logical(1L)))) 
                    .err("all entries in `split.with` must must be atomic vectors or factors")
                if (!all(lengths(split.with) == ncol(data)))
                    .err("all entries in `split.with` must have length equal to the number of columns of `data`")
            }
            else {
                if (!is.atomic(split.with))
                    .err("`split.with` must must be an atomic vector or factor or list thereof")
                if (length(split.with) != ncol(data))
                    .err("`split.with` must have length equal to the number of columns of `data`")
                split.with <- list(split.with)
            }
        }
    }
    else if (is.atomic(data)) {
        dep <- deparse1(substitute(data))
        data <- data.frame(data)
        if (missing(var.name) || is_null(var.name)) {
            names(data) <- dep
        }
        else if (is.atomic(var.name)) {
            if (length(var.name) == 1) {
                names(data) <- var.name
            }
            else {
                .wrn("only using the first item of `var.name`")
                names(data) <- var.name[1]
            }
        }
        else {
            .err("`var.name` must be an atomic or factor vector of length 1 with the stem of the new variable")
        }
        var.name <- names(data)
        
        if (is_not_null(split.with)) {
            if (is.list(split.with)) {
                if (!all(vapply(split.with, is.atomic, logical(1L)))) 
                    .err("all entries in `split.with` must must be atomic vectors or factors")
                if (!all(lengths(split.with) == ncol(data)))
                    .err("all entries in `split.with` must have length 1")
            }
            else {
                if (!is.atomic(split.with))
                    .err("`split.with` must must be an atomic vector or factor or list thereof")
                if (length(split.with) != ncol(data))
                    .err("`split.with` must have length 1")
                split.with <- list(split.with)
            }
        }
    }
    else {
        .err("`data` must a be a data.frame or factor")
    }
    
    if (is_not_null(drop.level) && length(var.name) > 1) {
        .wrn("`drop.level` cannot be used with multiple entries to `var.name`. Ignoring `drop.level`")
        drop.level <- NULL
    }
    
    .chk_flag(drop.singleton)
    drop.na <- rlang::rep_named(var.name, drop.na)
    
    for (v in var.name) {
        x <- factor(data[v][[1]])
        na.level <- if (anyNA(x)) nlevels(x) + 1 else 0
        new.levels <- c(levels(x), if (na.level) NA_character_)
        
        if (length(new.levels) > 1) {
            k <- make_df(paste(v, new.levels, sep = sep), nrow(data), types = "integer")
            for (i in seq_len(ncol(k))) {
                if (i != na.level) k[!is.na(x) & x == new.levels[i], i] <- 1L
                else k[is.na(x), i] <- 1L
            }
            
            if (na.level == 0) drop.na[v] <- FALSE
            else if (drop.na[v]) is.na(k)[is.na(x), -na.level] <- TRUE
            
        }
        else {
            .chk_flag(drop.singleton)
            if (drop.singleton) {
                data[[v]] <- NULL
                next
            }
            else {
                k <- make_df(paste(v, new.levels[1], sep = sep), length(x), types = "integer")
            }
        }
        
        dropl <- rlang::rep_named(new.levels, FALSE)
        if (is_not_null(drop.level)) {
            if (!is.character(drop.level) || length(drop.level) != 1 || drop.level %nin% new.levels) {
                .err(sprintf("`drop` must be the name of a level of %s that is to be dropped",
                             v))
            }
            dropl[drop.level] <- TRUE
        }
        else {
            if (!identical(drop.first, "if2") && !chk::vld_flag(drop.first)) {
                .err("`drop.first` must be `TRUE`, `FALSE`, or \"if2\"")
            }
            
            if ((ncol(k) == 2 && (drop.first == "if2" || isTRUE(drop.first))) ||
                (ncol(k) > 2 && isTRUE(drop.first))) {
                dropl[1] <- TRUE
            }
        }
        
        if (drop.na[v]) dropl[na.level] <- TRUE
        
        for (i in seq_len(ncol(k))) {
            attr(k[[i]], "split.var") <- v
            attr(k[[i]], "level") <- {if (i == na.level) NA_character_ else new.levels[i]}
        }
        
        k[dropl] <- NULL
        
        if (is_not_null(split.with)) {
            for (i in seq_along(split.with)) names(split.with[[i]]) <- names(data)
        }
        
        if (ncol(data) == 1) {
            data <- k
            if (is_not_null(split.with)) {
                for (i in seq_along(split.with)) {
                    split.with[[i]] <- rep(split.with[[i]][v], ncol(k))
                    names(split.with[[i]]) <- names(data)
                }
            }
        }
        else {
            .chk_flag(replace)
            if (replace) {
                if (match(v, names(data)) == 1){
                    data <- setNames(data.frame(k, data[names(data)!=v], row.names = rownames(data)),
                                     c(names(k), names(data)[names(data)!=v]))
                    if (is_not_null(split.with)) {
                        for (i in seq_along(split.with)) {
                            split.with[[i]] <- c(rep(split.with[[i]][v], ncol(k)), split.with[[i]][names(split.with[[i]]) != v])
                            names(split.with[[i]]) <- names(data)
                        }
                    }
                }
                else if (match(v, names(data)) == ncol(data)) {
                    data <- setNames(data.frame(data[names(data)!=v], k, row.names = rownames(data)),
                                     c(names(data)[names(data)!=v], names(k)))
                    if (is_not_null(split.with)) {
                        for (i in seq_along(split.with)) {
                            split.with[[i]] <- c(split.with[[i]][names(split.with[[i]]) != v], rep(split.with[[i]][v], ncol(k)))
                            names(split.with[[i]]) <- names(data)
                        }
                    }
                }
                else {
                    where <- match(v, names(data))
                    data <- setNames(data.frame(data[1:(where-1)], k, data[(where+1):ncol(data)], row.names = rownames(data)),
                                     c(names(data)[1:(where-1)], names(k), names(data)[(where+1):ncol(data)]))
                    if (is_not_null(split.with)) {
                        for (i in seq_along(split.with)) {
                            split.with[[i]] <- c(split.with[[i]][1:(where-1)], rep(split.with[[i]][v], ncol(k)), split.with[[i]][(where+1):length(split.with[[i]])])
                            names(split.with[[i]]) <- names(data)
                        }
                    }
                }
            }
            else {
                data <- setNames(data.frame(data, k, row.names = rownames(data)),
                                 c(names(data), names(k)))
                if (is_not_null(split.with)) {
                    for (i in seq_along(split.with)) {
                        split.with[[i]] <- c(split.with[[i]], rep(split.with[[i]][v], ncol(k)))
                        names(split.with[[i]]) <- names(data)
                    }
                }
            }
        }
        
    }
    
    attr(data, "split.with") <- split.with
    
    data
}

#' @export
#' @rdname splitfactor
unsplitfactor <- function(data, var.name, dropped.level = NULL, dropped.na = TRUE, sep = "_", replace = TRUE) {
    
    if (!is.data.frame(data)) .err("`data` must be a data.frame containing the variables to unsplit")
    if (missing(var.name)) {
        if (any(split.dummies <- vapply(data, function(x) {
            is_not_null(attr(x, "split.var")) && is_not_null(attr(x, "level")) &&
                all(x %in% c(0L, 1L, NA_integer_))
        }, logical(1L)))) {
            var.name <- unique(vapply(data[split.dummies], attr, character(1L), "split.var"))
        }
        else {
            .err("`var.name` must be a string containing the name of the variables to unsplit")
        }
    }
    else if (!is.character(var.name)) {
        .err("`var.name` must be a string containing the name of the variables to unsplit")
    }
    
    if (is_not_null(sep)) {
        if (!is.atomic(sep)) {
            .err("`sep` must be a character containing the seperating character in the names of the split variables. See `?unsplitfactor` for details")
        }
        
        if (length(sep) == 1L) sep <- rlang::rep_named(var.name, sep)
        else if (length(sep) == length(var.name)) names(sep) <- var.name
        else {
            .err("`sep` must be a character containing the seperating character in the names of the split variables. See `?unsplitfactor` for details")
        }
    }
    else sep <- rlang::rep_named(var.name, "")
    
    if (is_not_null(dropped.level)) {
        if (!is.atomic(dropped.level)) {
            .err("`dropped.level` must be an atomic vector containing the value of the dropped category of each split variable. See ?unsplitfactor for details")
        }
        
        if (length(dropped.level) == 1L) dropped.level <- rlang::rep_named(var.name, dropped.level)
        else if (length(dropped.level) == length(var.name)) names(dropped.level) <- var.name
        else {
            .err("`dropped.level` must be an atomic vector containing the value of the dropped category of each split variable. See ?unsplitfactor for details")
        }
    }
    else dropped.level <- NULL
    
    not.the.stem <- character(0)
    
    for (v in var.name) {
        dropped.level0 <- dropped.level[v]
        v.is.split <- any(split.data <- vapply(data, function(x) isTRUE(attr(x, "split.var", TRUE) == v), logical(1L)))
        
        var.to.combine <- if (v.is.split) data[split.data] else data[startsWith(names(data), paste0(v, sep[v]))]
        
        if (is_null(var.to.combine)) {
            not.the.stem <- c(not.the.stem, paste0(v, sep[v]))
            next
        }
        
        if (any(rowSums(apply(var.to.combine, 2, is.na)) %nin% c(0, ncol(var.to.combine)))) {
            .err("the variables in `data` selected based on `var.name` and `sep` do not seem to form a split variable based on the <NA> pattern")
        }
        NA.column <- character(0)
        if (v.is.split && any(na.dummy <- vapply(var.to.combine, function(x) is_not_null(s <- attr(x, "level", TRUE)) && is.na(s), logical(1L)))) {
            dropped.na <- FALSE
        }
        
        .chk_flag(dropped.na)
        if (!dropped.na) {
            NA.column <- {
                if (v.is.split) {
                    names(var.to.combine)[na.dummy]
                }
                else {
                    paste0(v, sep[v], {if (isFALSE(dropped.na)) "NA" else dropped.na})
                }
            }
            
            if (length(NA.column) > 1) {
                .err(sprintf("there appears to be more than one `NA` variable for %s", v))
            }
            
            if (NA.column %nin% names(var.to.combine)) {
                .err(sprintf("there is no variable called %s to generate the `NA` values",
                             add_quotes(NA.column)))
            }
            
            var.to.combine[var.to.combine[[NA.column]] == 1L,] <- NA_integer_
            var.to.combine[[NA.column]] <- NULL
        }
        
        var.sum <- rowSums(var.to.combine)
        
        if (all(is.na(var.sum) | check_if_zero(var.sum - 1))) {
            #Already unsplit
        }
        else if (all(is.na(var.sum) | check_if_zero(var.sum - 1) | check_if_zero(var.sum - 0))) {
            #Missing category
            
            if (is_null(dropped.level[v])) {
                
                k.levels0 <- if (v.is.split) {
                    unlist(lapply(var.to.combine, function(x) attr(x, "level", TRUE)))
                } else {
                    substr(names(var.to.combine), nchar(paste0(v, sep[v])), nchar(names(var.to.combine)))
                }
                
                if (can_str2num(k.levels0)) {
                    k.levels0.num <- as.numeric(k.levels0)
                    if (any({nums <- seq(min(k.levels0.num), max(k.levels0.num))} %nin% k.levels0.num)) {
                        dropped.level0 <- nums[nums %nin% k.levels0.num][1]
                    }
                    else if (min(as.numeric(k.levels0)) == 0) dropped.level0 <- max(k.levels0.num) + 1
                    else dropped.level0 <- min(k.levels0.num) - 1
                    dropped.name <- paste0(v, sep[v], dropped.level0)
                }
                else {
                    .msg(sprintf("the dropped category for %s will be set to `NA`",
                                 v))
                    dropped.name <- dropped.level0 <- NA_character_
                }
                
            }
            else dropped.name <- paste0(v, sep[v], dropped.level[v])
            
            var.to.combine <- setNames(data.frame(1 - var.sum, var.to.combine),
                                       c(dropped.name, names(var.to.combine)))
            
            if (v.is.split) {
                attr(var.to.combine[[1]], "split.var") <- v
                attr(var.to.combine[[1]], "level") <- as.character(dropped.level0)
            }
            
        }
        else {
            .err("the variables in `data` selected based on `var.name` and `sep` do not seem to form a split variable based on the row sums")
        }
        
        k.levels <- if (v.is.split) {
            unlist(lapply(var.to.combine, function(x) attr(x, "level", TRUE))) 
        }
        else {
            substr(names(var.to.combine), 1 + nchar(paste0(v, sep[v])), nchar(names(var.to.combine)))
        }
        
        k <- rep(NA_character_, nrow(data))
        for (i in seq_along(k.levels)) {
            k[var.to.combine[[i]] == 1] <- k.levels[i]
        }
        
        k <- factor(k, levels = k.levels)
        
        .chk_flag(replace)
        if (replace) {
            where <- which(names(data) %in% c(names(var.to.combine), NA.column))
            
            data[[where[1]]] <- k
            remove.cols <- where[-1]
            if (is_not_null(remove.cols)) data[remove.cols] <- NULL
            names(data)[where[1]] <- v
        }
        else {
            data[[v]] <- k
        }
    }
    
    if (is_not_null(not.the.stem)) {
        .wrn(sprintf("%s not the stem of any variables in `data` and will be ignored. Ensure `var.name` and `sep` are correct",
                     word_list(not.the.stem, is.are = TRUE, quotes = 2)))
    }
    
    data
}
