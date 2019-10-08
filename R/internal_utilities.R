#Internal Utilities not in SHARED.R
`%+%` <- function(...) {
    a <- list(...)
    if (is.atomic(a[[1]])) do.call(crayon::`%+%`, a)
    else do.call(ggplot2::`%+%`, a)
}
