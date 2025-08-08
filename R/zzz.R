.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  b <- utils::packageDate(pkgname)
  
  foo <- sprintf(" %s (Version %s, Build Date: %s)",
                 pkgname, v, if (anyNA(b)) "unknown" else format(b, "%F"))
  packageStartupMessage(foo)
}

globalVariables(c("density", "count"))
