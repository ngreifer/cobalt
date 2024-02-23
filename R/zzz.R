.onAttach <- function(libname, pkgname) {
    version <- utils::packageVersion(pkgname)
    BuildDate <- utils::packageDate(pkgname)
    
    foo <- paste0(" ", pkgname, " (Version ", version, ", Build Date: ", if (!anyNA(BuildDate)) format(BuildDate, "%F"), ")")
    packageStartupMessage(foo)
}

.onLoad <- function(libname, pkgname) {
    backports::import(pkgname)
}

globalVariables(c("density", "count"))