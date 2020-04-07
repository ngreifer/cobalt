\name{class-bal.tab.imp}
\alias{class-bal.tab.imp}
\alias{bal.tab.imp}

\title{Using \code{bal.tab()} with Multiply Imputed Data}

\description{
When using \code{bal.tab()} with multiply imputed data, the output will be different from the case with a single data set. Multiply imputed data can be used with the data.frame, formula, and \code{weightit} \code{bal.tab()} methods, and the \code{mimids} and \code{wimids} methods for \pkg{MatchThem} objects automatically incorporates multiply imputed data. This page outlines the outputs and options available with multiply imputed data.

There are two main components of the output of \code{bal.tab()} with multiply imputed data: the within-imputation balance summaries and the across-imputation balance summary. The within-imputation balance summaries display balance for units within each imputed data set separately. In general, this will not be very useful because interest rarely lies in the qualities of any individual imputed data set.

The across-imputation balance summary pools information across the within-imputation balance summaries to simplify balance assessment. It provides the average, smallest, and largest balance statistic for each covariate across all imputations. This allows you to see how bad the worst imbalance is and what balance looks like on average across the imputations. The summary behaves differently depending on whether \code{abs} is specified as \code{TRUE} or \code{FALSE}. When \code{abs = TRUE}, the across-imputation balance summary will display the mean absolute balance statistics and the maximum absolute balance statistics. When \code{abs = FALSE}, the across-imputation balance summary will display the minimum, mean, and maximum of the balance statistic in its original form.
}
\arguments{
There are four arguments for each \code{bal.tab()} method that can handle multiply imputed data: \code{imp}, \code{which.imp}, \code{imp.summary}, and \code{imp.fun}.
\item{imp}{A vector of imputation membership. This can be factor, character, or numeric vector. This argument is required to let \code{bal.tab()} know that the data is multiply imputed unless \pkg{MatchThem} objects are used. If a \code{data} argument is specified, this can also be the name of a variable in \code{data} that contains imputation membership. If the \code{data} argument is a \code{mids} object, the output of a call to \code{mice()}, \code{imp} does not need to be specified and will automatically be extracted from the \code{mids} object.}
\item{which.imp}{This is a display option that does not affect computation. If \code{.all}, all imputations in \code{imp} will be displayed. If \code{.none} (the default), no imputations will be displayed. Otherwise, can be a vector of imputation indices for which to display balance.}
\item{imp.summary}{This is a display option that does not affect computation. If \code{TRUE}, the balance summary across imputations will be displayed. The default is \code{TRUE}, and if \code{which.imp} is \code{.none}, it will automatically be set to \code{TRUE}.}
\item{imp.fun}{This is a display option that does not affect computation. Can be "min", "mean", or "max" and corresponds to which function is used in the across-imputation summary to combine results across imputations. For example, if \code{imp.fun = "mean"} the mean balance statistic across imputations will be displayed. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display all three. The default when \code{abs = FALSE} in the \code{bal.tab()} call is to display just the mean and max balance statistic.
}
}
\value{
The output is a \code{bal.tab.imp} object, which inherits from \code{bal.tab}. It has the following elements:

\item{Imputation.Balance}{For each imputation, a regular \code{bal.tab} object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified.}
\item{Balance.Across.Imputations}{The balance summary across imputations. This will include the combination of each balance statistic for each covariate across all imputations according to the value of \code{imp.fun}.}
\item{Observations}{A table of sample sizes or effective sample sizes averaged across imputations before and after adjustment.}

As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
}

\author{
Noah Greifer 
}

\seealso{
\code{\link{bal.tab}}, \code{\link{bal.tab.data.frame}}, \code{\link{print.bal.tab.imp}}
}
\keyword{tables}