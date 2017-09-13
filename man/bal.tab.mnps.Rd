\name{bal.tab.mnps}
\alias{bal.tab.mnps}
\title{
Balance Statistics for Multinomial twang Objects
}
\description{
Will eventually automatically generate balance statistics for \code{mnps} objects from \pkg{twang}. Presently, this page describes how to manually use \pkg{cobalt} to assess balance on \code{mnps} objects.
}
\usage{
\method{bal.tab}{mnps}(mnps, stop.method, ...)
}

\arguments{
  \item{mnps}{
an \code{mnps} object; the output of a call to \code{mnps()} from the \pkg{twang} package.
}
  \item{stop.method}{
a string containing the names of the stopping methods used in the original call to \code{mnps()}. Examples include "es.max" or "ks.mean". \code{bal.tab()} will assess balance for the weights created by those stopping methods. The names can be abbreviated as long as the abbreviations are specific enough. If no stopping methods are provided, \code{bal.tab()} will default to displaying balance for all available stopping methods.
  }
  \item{...}{
further arguments passed to or from other methods.
}
}
\details{
An \code{mnps} object contains a list of \code{ps} objects called "\code{psList}". For estimating the "ATE", the list will contain one \code{ps} object for each treatment group, and each one contains the weights for comparing that treatment group to the others combined. For estmating the "ATT", the list will contain one \code{ps} object for each group that is not the reference group (as specified by the user), and each one contains the weights for comparing that group to the reference group. Thus, to assess balance, it is possible to use \code{bal.tab()} on each item in the list seperately. Otherwise, use \pkg{twang}'s \code{bal.table()} command.

It is possible to visually essess balance on \code{mnps} objects using \code{\link{bal.plot}}. 
}
\value{
Currently, nothing. Throws an error.
}

\author{
Noah Greifer \email{noah@unc.edu}
}

\seealso{
\code{\link{bal.tab}} for details of multinomial balance assessment outside of \code{mnps} objects.
\code{\link{bal.tab.ps}} for \code{ps} objects returned by \pkg{twang}.
}
\examples{
library("twang"); data("lalonde", package = "cobalt")

\dontrun{
## Using mnps() for generalized boosted modeling

lalonde$treat3 <- factor(sample(1:3, nrow(lalonde),
                         replace = TRUE))
mnps.out <- mnps(treat3 ~ age + educ + married + race +
             nodegree + re74 + re75, data = lalonde, 
             stop.method = "es.mean", 
             estimand = "ATE", verbose = FALSE)
             
# Comparing group 1 with all others:
bal.tab(mnps.out$psList[[1]])
}
}
\keyword{design}
