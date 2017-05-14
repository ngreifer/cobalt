\name{lalonde}
\alias{lalonde}
\alias{lalonde_mis}
\docType{data}
\title{Lalonde's National Supported Work Demonstration data}
\description{
One of the datasets used by Dehejia and Wahba in their paper "Causal Effects in Non-Experimental Studies: Reevaluating the Evaluation of Training Programs." Versions of this data set have been used as an example data set in \pkg{MatchIt}, \pkg{twang}, \pkg{Matching}, and \pkg{CBPS}. The data set \code{lalonde_mis} is the same but with some values missing (set to \code{NA}).
}
\usage{
data("lalonde")
data("lalonde_mis")}
\format{
  A data frame with 614 observations on the following 9 variables.
  \describe{
    \item{\code{treat}}{1 if treated in the National Supported Work Demonstration,
                        0 if from the Current Population Survey}
    \item{\code{age}}{age}
    \item{\code{educ}}{years of education}
    \item{\code{race}}{factor; black, Hispanic (hispan), or white}
    \item{\code{married}}{1 if married, 0 otherwise}
    \item{\code{nodegree}}{1 if no degree, 0 otherwise}
    \item{\code{re74}}{earnings in 1974 (pretreatment)}
    \item{\code{re75}}{earnings in 1975 (pretreatment)}
    \item{\code{re78}}{earnings in 1978 (outcome)}
  }
}
\source{
  \url{http://users.nber.org/~rdehejia/data/nswdata2.html}
}
\references{
Lalonde, R. (1986). Evaluating the econometric evaluations of training programs
with experimental data. American Economic Review 76: 604-620.

Dehejia, R.H. and Wahba, S. (1999). Causal Effects in Nonexperimental Studies: 
Re-Evaluating the Evaluation of Training Programs. Journal of the American 
Statistical Association 94: 1053-1062.
}
\keyword{datasets}