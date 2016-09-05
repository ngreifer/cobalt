f.build <- function(y, rhs) {
    if ((is.data.frame(rhs) || is.matrix(rhs)) && length(colnames(rhs)) > 0)
        vars <- colnames(rhs)
    else if (is.character(rhs) && length(rhs) > 0)
        vars <- rhs
    else stop("Right hand side argument to f.build() must be a vector of variable names or a data set with named variables.", call. = FALSE)
    if (!(is.character(y) && length(y) == 1)) stop ("Response argument to f.build() must be the quoted name of the response variable.", call. = FALSE)
  f <- as.formula(paste(y," ~ ",paste(vars,collapse="+")))
  return(f)
}
