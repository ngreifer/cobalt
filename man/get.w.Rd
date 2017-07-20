\name{get.w}
\alias{get.w}
\alias{get.w.matchit}
\alias{get.w.ps}
\alias{get.w.Match}
\alias{get.w.CBPS}
\alias{get.w.ebalance}
\alias{get.w.optmatch}

\title{
Extract Weights from Preprocessing Objects
}
\description{
Extracts weights from the outputs of preprocessing functions.
}
\usage{
get.w(\dots)

\method{get.w}{matchit}(m, ...)

\method{get.w}{ps}(ps, stop.method, ...)

\method{get.w}{Match}(M, ...)

\method{get.w}{CBPS}(c, use.weights = TRUE, estimand = NULL, ...)

\method{get.w}{ebalance}(e, treat, ...)

\method{get.w}{optmatch}(o, treat, ...)
}

\arguments{
  \item{m, ps, M, c, e, o}{
Output from the corresponding preprocessing packages.
}
  \item{stop.method}{
The name of the stop method used in the original call to \code{ps()} in \pkg{twang}, e.g., \code{"es.mean"}. If empty, will use the first stop method available. Abbreviations allowed.
}
  \item{use.weights}{
\code{logical}; whether to use the weights output by the call to \code{CBPS()} or to calculate weights based on the included propensity scores. \code{CBPS()} standardizes weights by default and normalizes them, so in general the computed weights will not be equal to the weights computed using the IPW or weighting by the odds formulas. In prinicple, weighted computations will not depend on whether the original or calculated weights are used, but if propertes of weights are to compared with weights generated through manual computation, it can be useful to export the calculated weights.
}
  \item{estimand}{
If weights are computed using the propensity score (i.e., if \code{use.weights = TRUE}), which estimand to use to compute the weights. If \code{"ATE"}, weights will compute as \code{1/ps} for the treated group and \code{1/(1-ps)} for the control group. If \code{"ATT"}, weights will compute as \code{1} for the treated group and \code{ps/(1-ps)} for the control group. If \code{NULL}, \code{get.w()} will try to figure out which estimand is desired based on the \code{CBPS} output.
}
  \item{treat}{
A vector of treatment status for each unit. This is required for methods that include \code{treat} as an argument. The treatment variable that was used in the original preprocessing function call should be used.
}
  \item{...}{
further arguments passed to or from other methods.
}
}
\details{
The output of \code{get.w()} can be used in calls to the formula and data frame methods of \code{\link{bal.tab}()} (see example below). In this way, the output of multiple preprocessing packages can be viewed simultaneously and compared. The weights can also be used in \code{weights} statements in regression methods to compute causal effects.

\pkg{twang} has a funciton called \code{get.weights()} that performs the same function on \code{ps} objects but offers finer control.
}
\value{
A vector of weights for each unit. These may be matching weights or balancing weights.
}


\author{
Noah Greifer \email{noah@unc.edu}
}

\seealso{
\code{\link[twang]{get.weights}} in \pkg{twang}.
}
\examples{
data("lalonde", package = "cobalt")
library("MatchIt"); library("CBPS")

m.out <- matchit(treat ~ age + educ + race, data = lalonde) 

c.out <- CBPS(treat ~ age + educ + race, data = lalonde)

bal.tab(treat ~ age + educ + race, data = lalonde,
        weights = data.frame(matchit = get.w(m.out),
                             cbps = get.w(c.out)),
        method = c("matching", "weighting"), 
        estimand = "ATT")
}
