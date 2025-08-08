#Use regex to make strings invariant to white spaces
.w <- function(x) {
  gsub(" ", "(\\s+)", x, fixed = TRUE)
}