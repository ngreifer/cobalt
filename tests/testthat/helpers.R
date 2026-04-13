#Use regex to make strings invariant to white spaces; also escape perens for perl
.w <- function(x) {
  x <- gsub("(", "\\(", x, fixed = TRUE)
  x <- gsub(")", "\\)", x, fixed = TRUE)
  gsub(" ", "(\\s+)", x, fixed = TRUE)
}