\name{bal.tab.cem.match}
\alias{bal.tab.cem.match}
\title{
Balance Statistics for \code{cem} Objects
}
\description{
Generates balance statistics for \code{cem.match} objects from \pkg{cem}.
}
\usage{
\method{bal.tab}{cem.match}(x, 
        data,
        stats,
        int = FALSE,
        poly = 1,
        distance = NULL,
        addl = NULL,
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
a \code{cem.match} or \code{cem.match.list} object; the output of a call to \pkgfun{cem}{cem}.
}
  \item{data}{
a data frame containing variables named in other arguments. An argument to \code{data} is \strong{required}. It must be the same data used in the call to \code{cem()} or a \code{mids} object from which the data supplied to \code{datalist} in the {cem()} call originated.
}
  \item{stats, int, poly, distance, addl, continuous, binary, thresholds, weights, cluster, imp, pairwise, s.weights, abs, subset, quick, \dots}{
see \fun{bal.tab} for details.

See below for a special note on the \code{s.d.denom} argument.
}
The following argument has a special note when used with \code{cem.match} or \code{cem.match.list} objects:
  \item{s.d.denom}{
the default is "treated", where the treated group corresponds to the \code{baseline.group} in the call to \code{cem()}.
}
}
\details{
\code{bal.tab.cem.match()} generates a list of balance summaries for the \code{cem.match} object given, and functions similarly to \pkgfun{cem}{imbalance}.
}
\value{
If clusters and imputations are not specified, an object of class \code{"bal.tab"} containing balance summaries for the \code{cem.match} object. See \fun{bal.tab} for details.

If imputations are specified, an object of class \code{"bal.tab.imp"} containing balance summaries for each imputation and a summary of balance across imputations. See \code{\link[=class-bal.tab.imp]{bal.tab.imp}} for details.

If \code{cem()} is used with multi-category treatments, an object of class \code{"bal.tab.multi"} containing balance summaries for each pairwise treatment comparison. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}} for details.

If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for details.
}

\author{
Noah Greifer 
}

\seealso{
\fun{bal.tab} for details of calculations.
}
\examples{\dontshow{if (requireNamespace("cem", quietly = TRUE)) force(\{ # examplesIf}
library(cem); data("lalonde", package = "cobalt")

## Coarsened exact matching
cem.out <- cem("treat", data = lalonde, drop = "re78")

bal.tab(cem.out, data = lalonde, un = TRUE, 
        stats = c("m", "k"))
\dontshow{\}) # examplesIf}}
\keyword{tables}
