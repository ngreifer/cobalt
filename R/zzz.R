.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  b <- utils::packageDate(pkgname)
  
  foo <- sprintf(" %s (Version %s, Build Date: %s)", pkgname, v, if (!anyNA(b)) format(b, "%F") else "unknown")
  packageStartupMessage(foo)
}

globalVariables(c("density", "count"))