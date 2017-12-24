\name{bal.tab.subclass}
\alias{bal.tab.subclass}

\title{Using \code{bal.tab()} with Subclassified Data}

\description{
When using \code{bal.tab()} with subclassified data, i.e., data split into subclasses where balance may hold, the output will be different from the standard, non-subclassified case, and there is an additional option for controling display. This page outlines the outputs and options in this case.

There are two main componenets of the output of \code{bal.tab()} with subclassified data: the balance within subclasses and the balance summary across subclasses. The within-subclass balance displays essentially are standard balance displays for each subclass, except that only "adjusted" values are available ,ebcauset he subclassification itself is the adjustment.

The balance summary is, for each variable, a weighted average of the means and mean differences across subclasses. Each subclass statistic is weighed by the proportion of units corresponding to the target population (i.e., treated, control, or total sample) that are in that subclass. Variance ratios and KS statistics do not appear in this summary. Although this summary can be useful as a heuristic, it is important to ensure balance within each subclass as well.

}
\arguments{
There are two arguments for \code{bal.tab()} that relate to subclasses: \code{subclass} and \code{disp.subclass}.
\item{subclass}{For the data.frame and formula methods of \code{bal.tab}, a vector of subclass membership or the name of the variable in \code{data} containing subclass membership. When using subclassification with a function compatible with \pkg{cobalt}, such as \code{matchit()} in \pkg{MatchIt}, this argument can be omitted because the subclass are in the output object.}
\item{disp.subclass}{This is a display option that does not affect computation. If \code{TRUE}, balance for each subclass will be displayed. This can be cumbersome with many subclasses or covariates. The default is \code{TRUE}.}
}
\value{
The output is a \code{bal.tab.subclass} object, which inherits from \code{bal.tab}. It has the following elements:

\item{Subclass.Balance}{A list of data frames containing balance information for each covariate in each subclass.  Each data frame contains the following columns:
\itemize{
\item{\code{Type}: Whether the covariate is binary, continuous, or a measure of distance (e.g., the propensity score).}
\item{\code{M.0.Adj}: The mean of the control group in the subclass.}
\item{\code{M.1.Adj}: The mean of the treated group in the subclass.}
\item{\code{Diff.Adj}: The (standardized) difference in means between the two groups in the subclass.}
\item{\code{M.Threshold}: Whether or not the calculated mean difference exceeds or is within the threshold given by \code{m.threshold}.  If \code{m.threshold} is \code{NULL}, this column will be \code{NA}.}
\item{\code{V.Ratio.Adj}: The ratio of the variances of the two groups in the subclass.  \code{NA} for binary variables.  If less than 1, the reciprocal is reported.}
\item{\code{V.Threshold}: Whether or not the calculated variance ratio exceeds or is within the threshold given by \code{v.threshold} for continuous variables.  If \code{v.threshold} is \code{NULL}, this column will be \code{NA}.}
\item{\code{KS.Adj}: The KS statistic for the two groups in the subclass.  \code{NA} for binary variables.}
\item{\code{KS.Threshold}: Whether or not the calculated KS statistic exceeds or is within the threshold given by \code{ks.threshold} for continuous variables.  If \code{ks.threshold} is \code{NULL}, this column will be \code{NA}.}
}}

\item{Balance.Across.Subclass}{A data frame containing balance statistics for each covariate aggregated across subclasses and for the original sample (i.e., unadjusted).  Variance ratios and KS statistics are not reported here.}
\item{Balanced.Means.Subclass}{If \code{m.threshold} is specified, a table tallying the number of variables in each subclass that exceed or are within the threshold for mean differences.}
\item{Max.Imbalance.Means.Subclass}{If \code{m.threshold} is specified, a table displaying the variable in each subclass with the greatest absolute mean difference.}
\item{Balanced.Variances.Subclass}{If \code{v.threhsold} is specified, a table tallying the number of variables in each subclass that exceed or are within the threshold for variance ratios.}
\item{Max.Imbalance.Variance.Subclass}{If \code{v.threshold} is specified, a table displaying the variable in each subclass with the greatest variance ratio.}
\item{Balanced.KS.Subclass}{If \code{ks.threhsold} is specified, a table tallying the number of variables in each subclass that exceed or are within the threshold for KS statistics.}
\item{Max.Imbalance.KS.Subclass}{If \code{ks.threshold} is specified, a table displaying the variable in each subclass with the greatest KS statistics.}
\item{Subclass.Observations}{A table displaying the sample sizes in each subclass.}
}

\author{
Noah Greifer 
}

\seealso{
\code{link{bal.tab}}, \code{\link{bal.tab.data.frame}}, \code{\link{print.bal.tab.subclass}}
}