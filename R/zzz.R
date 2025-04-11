.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  b <- utils::packageDate(pkgname)
  
  foo <- paste0(" ", pkgname, " (Version ", v, ", Build Date: ", if (!anyNA(b)) format(b, "%F"), ")")
  packageStartupMessage(foo)
}

globalVariables(c("density", "count"))