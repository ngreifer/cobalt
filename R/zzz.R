.onAttach <- function(libname, pkgname) {
    version <- packageVersion(pkgname)
    BuildDate <- packageDate(pkgname)
    
    foo <- paste0(" ", pkgname, " (Version ", version, ", Build Date: ", if (!anyNA(BuildDate)) format(BuildDate, "%F"), ")")
    packageStartupMessage(foo)
}

.onLoad <- function(libname, pkgname) {
    backports::import(pkgname)
}

globalVariables(c("density", "count"))