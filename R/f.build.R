#' @title Convenient Formula Generation
#' 
#' @description `f.build()` returns a [`formula`] of the form `y ~ x1 + x2 + ...` from a data frame input. It can be much quicker to use `f.build()` than to hand-write the precise formula, which may contain errors. It can be used in place of a formula in, for example, [glm()], `matchit()`, or [bal.tab()]. It provides similar functionality to [reformulate()].
#' 
#' @param y the quoted name of the response (left hand side) variable in the formula. Only one variable is supported. If missing, `NULL`, or the empty string (`""`), the formula will have no response variable. If `rhs` is not supplied, `y` will replace `rhs` and `y` will be set to `""`.
#' @param rhs a data frame whose variable names will be the terms on the right hand side of the formula, or a character vector whose values will be the terms on the right hand side of the formula. If missing, the argument to `y` will replace `rhs` and `y` will be set to `""`; in essence, `f.build("x")` is the same as `f.build("", "x")`, both producing `~ x`.
#' 
#' @returns a `formula` object.
#' 
#' @seealso [reformulate()]
#' 
#' @examples 
#' data(lalonde)
#' covs <- subset(lalonde, select = -c(treat, re78))
#' lm(f.build("treat", covs), data = lalonde)

#' @rdname f.build
#' @export 
f.build <- function(y, rhs) {
    if (missing(rhs)) {
        if (missing(y)) {
            .err("the right hand side argument to `f.build()` must be a vector of variable names or a data set with named variables")
        }
        else {
            rhs <- y
            y <- ""
        }
    }
    
    tryCatch(force(y), error = function(e) .err(conditionMessage(e), tidy = FALSE))
    tryCatch(force(rhs), error = function(e) .err(conditionMessage(e), tidy = FALSE))
    
    if (is_mat_like(rhs) && is_not_null(colnames(rhs))) {
        vars <- paste0("`", gsub("`", "", colnames(rhs)), "`")
    }
    else if (is.character(rhs)) {
        vars <- add_quotes(gsub("`", "", rhs), "`")
    }
    else .err("the right hand side argument to `f.build()` must be a vector of variable names or a data set with named variables")
    
    if (missing(y) || is_null(y) || identical(y, "")) y <- NULL
    else if (!is.atomic(y)) {
        .err("the response argument to `f.build()` must be a string containing the response variable")
    }
    
    f <- formula(paste(
        paste(as.character(y), collapse = " + ") , "~", paste(vars, collapse = " + ")
    ))
    # f <- reformulate(vars, y)
    return(f)
}
