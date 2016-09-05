\name{bal.tab.matchit}
\alias{bal.tab.matchit}
\title{
Balance Statistics for Matchit Objects
}
\description{
Generates balance statistics for \code{matchit} objects from \pkg{MatchIt}.
}
\usage{
\method{bal.tab}{matchit}(x, int = FALSE, addl = NULL, continuous = c("std", 
                "raw"), binary = c("raw", "std"), 
                s.d.denom = c("treated", "control", "pooled"), 
                m.threshold = NULL, v.threshold = NULL, un = 
                FALSE, disp.means = FALSE, disp.v.ratio = 
                FALSE, disp.subclass = FALSE, cluster = NULL, 
                which.cluster = NULL, cluster.summary = TRUE, 
                ...)
}

\arguments{
  \item{x}{
a \code{matchit} object; the output of a call to \code{matchit()} from the \pkg{MatchIt} package.
}
  \item{int}{
\code{logical}; whether or not to include 2-way interactions of covariates included in the call to \code{matchit()} and in \code{addl}.
}
  \item{addl}{
a data frame of additional covariates for which to present balance. These may be covariates included in the original dataset but not included in the call to \code{matchit()}. See note in Details.
}
  \item{continuous}{
whether mean differences for continuous variables should be standardized ("std") or raw ("raw"). Default "std". Abbreviations allowed.
}
  \item{binary}{
whether mean differences for binary variables (i.e., difference in proportion) should be standardized ("std") or raw ("raw"). Default "raw". Abbreviations allowed.
}
  \item{s.d.denom}{
whether the denominator for standardized differences (if any are calculated) should be the standard deviation of the treated group ("treated"), the standard deviation of the control group ("control"), or the pooled standard deviation ("pooled"), computed as the square root of the mean of the group variances. Abbreviations allowed. The default is "treated".
}
  \item{m.threshold}{
a numeric value for the threshold for mean differences. .1 is recommended. 
}
  \item{v.threshold}{
a numeric value for the threshold for variance ratios. Will automatically convert to the inverse if less than 1.
}
  \item{un}{
\code{logical}; whether to print statistics for the unadjusted sample as well as for the adjusted sample.
}
  \item{disp.means}{
\code{logical}; whether to print the group means in balance output.
}
  \item{disp.v.ratio}{
\code{logical}; whether to display variance ratios in balance output.
}
  \item{disp.subclass}{
\code{logical}; whether to display balance information for individual subclasses if subclassification is used in conditioning.
}
  \item{cluster}{
a vector of cluster membership.
}
  \item{which.cluster}{
which cluster(s) to display. If \code{NULL}, all clusters in \code{cluster} will be displayed. If \code{NA}, no clusters will be displayed. Otherwise, can be a vector of cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names. 
}
  \item{cluster.summary}{
\code{logical}; whether to display the cluster summary table. If \code{which.cluster} is \code{NULL}, \code{cluster.summary} will be set to \code{TRUE}.
}
  \item{...}{
further arguments passed to or from other methods. They are ignored in this function.
}
}
\details{
\code{bal.tab.matchit()} generates a list of balance summaries for the \code{matchit} object given, and functions similarly to \code{summary.matchit()} in \pkg{MatchIt}. \code{bal.tab()} behaves differently depending on whether subclasses are used in conditioning or not. If they are used, \code{bal.tab()} creates balance statistics for each subclass and for the sample in aggregate. 

The last four arguments of \code{bal.tab()} affect display only; they are passed directly to \code{print.bal.tab} or \code{print.bal.tab.subclass}, and do not affect any calculations or the contents of the \code{bal.tab} object. All balance statistics are calculated whether they are displayed by print or not. The threshold values (\code{m.threshold} and \code{v.threshold}) control whether extra columns should be inserted into the Balance table describing whether the balance statistics in question exceeded or were within the threshold. Including these thresholds also creates summary tables tallying the number of variables that exceeded and were within the threshold and displaying the variables with the greatest imbalance on that balance measure. When subclassification is used, the extra threshold columns are placed within the balance tables for each subclass as well as in the aggregate balance table, and the summary tables display balance for each subclass.

The input to \code{addl} must be a data frame; if more than one variable is included, this is straightforward (i.e., because \code{data[,c("v1", "v2")]} is already a data frame), but if only one variable is used, R will coerce it to a vector, thus making it unfit for input here. To avoid this, simply wrap the input to \code{addl} in \code{data.frame()} or use \code{subset()} if only one variable is to be added. Again, when more than one variable is included, the input is general already a data frame and nothing needs to be done.
}
\value{
If clusters are not specified, an object of class \code{"bal.tab"} containing balance summaries for the \code{matchit} object. If subclassifcation is not used, the following are the elements of \code{bal.tab}:
\item{Balance}{A data frame containing balance information for each covariate.  Balance contains the following columns:
\itemize{
\item{\code{Type}: Whether the covariate is binary, continuous, or a measure of distance (e.g., the propensity score).}
\item{\code{M.C.Un}: The mean of the control group prior to adjusting.}
\item{\code{M.T.Un}: The mean of the treated group prior to adjusting.}
\item{\code{Diff.Un}: The (standardized) difference in means between the two groups prior to adjusting.}
\item{\code{V.Ratio.Un}: The ratio of the variances of the two groups prior to adjusting.  \code{NA} for binary variables.  If less than 1, the reciprocal is reported.}
\item{\code{M.C.Adj}: The mean of the control group after adjusting.}
\item{\code{M.T.Adj}: The mean of the treated group after adjusting.}
\item{\code{Diff.Adj}: The (standardized) difference in means between the two groups after adjusting.}
\item{\code{M.Threshold}: Whether or not the calculated mean difference after adjusting exceeds or is within the threshold given by \code{m.threshold}.  If \code{m.threshold} is \code{NULL}, this column will be \code{NA}.}
\item{\code{V.Ratio.Adj}: The ratio of the variances of the two groups after adjusting.  \code{NA} for binary variables.  If less than 1, the reciprocal is reported.}
\item{\code{V.Threshold}: Whether or not the calculated variance ratio after adjusting exceeds or is within the threshold given by \code{v.threshold} for continuous variables.  If \code{v.threshold} is \code{NULL}, this column will be \code{NA}.}
}}
\item{Balanced.Means}{If \code{m.threshold} is specified, a table tallying the number of variables that exceed or are within the threshold for mean differences.}
\item{Max.Imbalance.Means}{If \code{m.threshold} is specified, a table displaying the variable with the greatest absolute mean difference.}
\item{Balanced.Variances}{If \code{v.threshold} is specified, a table tallying the number of variables that exceed or are within the threshold for variance ratios.}
\item{Max.Imbalance.Variance}{If \code{v.threshold} is specified, a table displaying the variable with the greatest variance ratio.}
\item{Observations}{A table displaying the sample sizes before and after adjusting.}
\item{call}{The original \code{matchit()} call.}
\item{print.options}{A list of print options passed to \code{print.bal.tab}.}

If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. Each balance summary is a balance table as described in \code{Balance} above. The summary of balance across clusters displays the mean, median, and maximum mean difference and variance ratio after adjustment for each covariate across clusters.

If subclassification is used, the following are the elements of \code{bal.tab}:
\item{Subclass.Balance}{A list of data frames containing balance information for each covariate in each subclass.  Each data frame contains the following columns:
\itemize{
\item{\code{Type}: Whether the covariate is binary, continuous, or a measure of distance (e.g., the propensity score).}
\item{\code{M.C.Adj}: The mean of the control group in the subclass.}
\item{\code{M.T.Adj}: The mean of the treated group in the subclass.}
\item{\code{Diff.Adj}: The (standardized) difference in means between the two groups in the subclass.}
\item{\code{M.Threshold}: Whether or not the calculated mean difference exceeds or is within the threshold given by \code{m.threshold}.  If \code{m.threshold} is \code{NULL}, this column will be \code{NA}.}
\item{\code{V.Ratio.Adj}: The ratio of the variances of the two groups in the subclass.  \code{NA} for binary variables.  If less than 1, the reciprocal is reported.}
\item{\code{V.Threshold}: Whether or not the calculated variance ratio exceeds or is within the threshold given by \code{v.threshold} for continuous variables.  If \code{v.threshold} is \code{NULL}, this column will be \code{NA}.}
}}

\item{Balance.Across.Subclass}{A data frame containing balance statistics for each covariate aggregated across subclasses and for the original sample (i.e., unadjusted).  Variance ratios are not reported here.}
\item{Balanced.Means.Subclass}{If \code{m.threshold} is specified, a table tallying the number of variables in each subclass that exceed or are within the threshold for mean differences.}
\item{Max.Imbalance.Means.Subclass}{If \code{m.threshold} is specified, a table displaying the variable in each subclass with the greatest absolute mean difference.}
\item{Balanced.Variances.Subclass}{If \code{v.threhsold} is specified, a table tallying the number of variables in each subclass that exceed or are within the threshold for variance ratios.}
\item{Max.Imbalance.Variance.Subclass}{If \code{v.threshold} is specified, a table displaying the variable in each subclass with the greatest variance ratio.}
\item{Subclass.Observations}{A table displaying the sample sizes in each subclass.}
\item{call}{The original \code{matchit()} call.}
\item{print.options}{A list of print options passed to \code{print.bal.tab}.}
}

\author{
Noah Greifer \email{noah@unc.edu}
}

\seealso{
\code{\link{bal.tab}} for details of calculations.
}
\examples{
library(MatchIt); data("lalonde", package = "cobalt")

## Nearest Neighbor matching
m.out1 <- matchit(treat ~ age + educ + black + hispan + 
                  married + nodegree + re74 + re75, 
                  data = lalonde, method = "nearest")
bal.tab(m.out1, un = TRUE, m.threshold = .1, 
        v.threshold = 2)

## Subclassification
m.out2 <- matchit(treat ~ age + educ + black + hispan + 
                  married + nodegree + re74 + re75, 
                  data = lalonde, method = "subclass")
bal.tab(m.out2, disp.subclass = TRUE)

}
\keyword{design}
