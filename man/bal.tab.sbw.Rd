\name{bal.tab.sbwcau}
\alias{bal.tab.sbwcau}
\title{
Balance Statistics for \code{sbw} Objects
}
\description{
Generates balance statistics for \code{sbwcau} objects from \pkg{sbw}.
}
\usage{
\method{bal.tab}{sbwcau}(x, 
        stats,
        int = FALSE,
        poly = 1,
        distance = NULL,
        addl = NULL,
        data = NULL,
        continuous,
        binary,
        s.d.denom,
        thresholds = NULL,
        weights = NULL,
        cluster = NULL,
        imp = NULL,
        pairwise = TRUE,
        s.weights = NULL,
        abs = FALSE,
        subset = NULL,
        quick = TRUE,
        ...)
}

\arguments{
  \item{x}{
an \code{sbwcau} object; the output of a call to \pkgfun{sbw}{sbw}.
}
  \item{stats, int, poly, distance, addl, data, continuous, binary, thresholds, weights, cluster, imp, pairwise, s.weights, abs, subset, quick, \dots}{
see \fun{bal.tab} for details.

See below for a special note on the \code{s.d.denom} argument.
}
The following argument has a special note when used with \code{sbwcau} objects:
  \item{s.d.denom}{
if not specified, \code{bal.tab()} will figure out which one is best based on the \code{par} component of the \code{sbwcau} object: if "att", "treated"; if "atc", "control"; otherwise "pooled".
}
}
\details{
\code{bal.tab.sbwcau()} generates a list of balance summaries for the \code{sbwcau} object given, and functions similarly to \pkgfun{sbw}{summarize}.
}
\value{

If clusters are not specified, an object of class \code{"bal.tab"} containing balance summaries for the \code{sbwcau} object. See \fun{bal.tab} for details.

If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for details.
}

\author{
Noah Greifer 
}

\seealso{
\fun{bal.tab} for details of calculations.
}
\examples{\dontshow{if (requireNamespace("sbw", quietly = TRUE)) force(\{ # examplesIf}
library(sbw); data("lalonde", package = "cobalt")

## Stable balancing weights for the ATT
sbw.out <- sbw(splitfactor(lalonde, drop.first = "if2"),
               ind = "treat",
               bal = list(bal_cov = c("age", "educ", "race_black", 
                                      "race_hispan", "race_white", 
                                      "married", "nodegree", 
                                      "re74", "re75"),
                          bal_alg = FALSE, 
                          bal_tol = .001),
               par = list(par_est = "att"))
               
bal.tab(sbw.out, un = TRUE, poly = 2)
\dontshow{\}) # examplesIf}}
\keyword{tables}
