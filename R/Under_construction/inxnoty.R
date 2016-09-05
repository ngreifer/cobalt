inxnoty <- function(x, y) {
  #Creates a list or data frame of names in x that are not in y.
  #Useful for subsetting data sets into two groups of variables.
  if (is.character(x)) X <- x else X <- names(x)
  if (is.character(y)) Y <- y else Y <- names(y)
  
  if(is.character(x)) out <- x[-which(X %in% Y)]
  else out <- x[,-which(X %in% Y)]
  return(out)
}