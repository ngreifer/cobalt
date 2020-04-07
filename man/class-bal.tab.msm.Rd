\name{class-bal.tab.msm}
\alias{class-bal.tab.msm}
\alias{bal.tab.msm}

\title{Using \code{bal.tab()} with Longitudinal Treatments}
\description{
When using \code{bal.tab()} with longitudinal treatments, the output will be different from the case with point treatments, and there are some options that are common across all \code{bal.tab} methods for dealing with longitudinal data. This page outlines the outputs and options in this case.

There are two main components of the output of \code{bal.tab()} with longitudinal treatments: the time-point-specific balance summary and across-time-points balance summary. The time-point-specific balance summaries are standard point treatment balance summaries at each time point.

The across-time-points balance summary is, for each variable, the greatest imbalance across all time-point-specific balance summaries. If the greatest observed imbalance is tolerable, then all other imbalances for that variable will be tolerable too, so focusing on reducing the greatest imbalance is sufficient for reducing imbalance overall. The balance summary will not be computed if multi-category treatments or multiply imputed data are used.
}
\arguments{
There are two additional arguments for each \code{bal.tab()} method that can handle longitudinal treatments: \code{which.time} and \code{msm.summary}.
\item{which.time}{This is a display option that does not affect computation. If \code{.all} (the default), all time points will be displayed. If \code{.none}, no time points will be displayed. Otherwise, can be a vector of treatment names or indices for which to display balance.}
\item{msm.summary}{This is a display option that does not affect computation. If \code{TRUE}, the balance summary across time points will be displayed. The default is \code{TRUE}, and if \code{which.time} is \code{.none}, it will automatically be set to \code{TRUE}.}
}
\value{
The output is a \code{bal.tab.msm} object, which inherits from \code{bal.tab}. It has the following elements:

\item{Time.Balance}{For each time point, a regular \code{bal.tab} object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified.}
\item{Balance.Across.Times}{The balance summary across time points. This will include the maximum balance statistic(s) for each covariate across all time points.}
\item{Observations}{A table of sample sizes or effective sample sizes for each time point before and after adjustment.}

As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
}

\author{
Noah Greifer 
}

\seealso{
\code{\link{bal.tab}}, \code{\link{bal.tab.list}}, \code{\link{print.bal.tab.msm}}
}
\keyword{tables}