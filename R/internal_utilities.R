#Internal Utilities not in SHARED.R
`%+%` <- function(...) {
    if (is_(..1, c("atomic", "factor")) && is_(..2, c("atomic", "factor"))) crayon::`%+%`(as.character(..1), as.character(..2))
    else ggplot2::`%+%`(...)
}
