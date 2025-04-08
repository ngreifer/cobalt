.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  b <- utils::packageDate(pkgname)
  
  foo <- paste0(" ", pkgname, " (Version ", v, ", Build Date: ", if (!anyNA(b)) format(b, "%F"), ")")
  packageStartupMessage(foo)
}

.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}

globalVariables(c("density", "count"))