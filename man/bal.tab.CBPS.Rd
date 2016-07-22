\name{bal.tab.CBPS}
\alias{bal.tab.CBPS}

\title{
Balance statistics for CBPS Objects
}
\description{
Generates balance statistics for \code{CBPS} objects from the \pkg{CBPS} package.
}
\usage{
\method{bal.tab}{CBPS}(x, estimand = NULL, int = FALSE, 
           addl = NULL, continuous = c("std", "raw"), 
           binary = c("raw", "std"), s.d.denom, 
           m.threshold = NULL, v.threshold = NULL, 
           un = FALSE, disp.means = FALSE, 
           disp.v.ratio = FALSE, ...)
}

\arguments{
  \item{x}{
a \code{CBPS} object; the output of a call to \code{CBPS()} from the \pkg{CBPS} package.
}
  \item{estimand}{
whether the desired estimand is the "ATT" or "ATE". Used for determining \code{s.d.denom} if not specified. \code{bal.tab()} will determine this automatically if unstandardized weights were requested in the original call to \code{CBPS()}. This argument is only required if standardized weights were requested and if \code{s.d.denom} is left blank.
}
  \item{int}{
\code{logical}; whether or not to include 2-way interactions of covariates included in the call to \code{CBPS()} and in \code{addl}.
}
  \item{addl}{
a data frame of additional covariates for which to present balance. These may be covariates included in the original dataset but not included in the call to \code{CBPS()}. See note in Details.
}
  \item{continuous}{
whether mean differences for continuous variables should be standardized ("std") or raw ("raw"). Default "std". Abbreviations allowed.
}
  \item{binary}{
whether mean differences for binary variables (i.e., difference in proportion) should be standardized ("std") or raw ("raw"). Default "raw". Abbreviations allowed.
}
  \item{s.d.denom}{
whether the denominator for standardized differences (if any are calculated) should be the standard deviation of the treated group ("treated"), the standard deviation of the control group ("control"), or the pooled standard deviation ("pooled"), computed as the square root of the mean of the group variances. Abbreviations allowed. If not specified, \code{bal.tab()} will use "treated" if the estimand of the call to \code{CBPS()} (or that specified by the user to \code{estimand}) is the ATT and "pooled" if the estimand is the ATE.
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
  \item{...}{
further arguments passed to or from other methods. They are ignored in this function.
}
}
\details{
\code{bal.tab.CBPS} generates a list of balance summaries for the \code{CBPS} object given, and functions similarly to \code{balance()} in \pkg{CBPS}.

The last four arguments of \code{bal.tab()} affect display only; they are passed directly to \code{print.bal.tab}, and do not affect any calculations or the contents of the \code{bal.tab} object. All balance statistics are calculated whether they are displayed by \code{print} or not. The threshold values (\code{m.threshold} and \code{v.threshold}) control whether extra columns should be inserted into the Balance table describing whether the balance statistics in question exceeded or were within the threshold. Including these thresholds also creates summary tables tallying the number of variables that exceeded and were within the threshold and displaying the variables with the greatest imbalance on that balance measure.

The input to \code{addl} must be a data frame; if more than one variable is included, this is straightforward (i.e., because \code{data[,c("v1", "v2")]} is already a data frame), but if only one variable is used, R will coerce it to a vector, thus making it unfit for input here. To avoid this, simply wrap the input to \code{addl} in \code{data.frame()} or use \code{subset()} if only one variable is to be added. Again, when more than one variable is included, the input is general already a data frame and nothing needs to be done.

When standardized weights are requested in the call to \code{CBPS()}, variance ratios may not be computed, and if they are computed, they may not be accurate. Mean differences, standardized or not, should be unaffected. For this reason, unstandardized weights are preferred.
}
\value{
An object of class \code{"bal.tab"} containing balance summaries for the \code{CBPS} object. The following are the elements of \code{bal.tab}:

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
\item{call}{The original \code{CBPS()} call.}
\item{print.options}{A list of print options passed to \code{print.bal.tab}.}
}

\author{
Noah Greifer \email{noah@unc.edu}
}


\seealso{
\code{\link{bal.tab}} for details of calculations.
}
\examples{
library(CBPS); data("lalonde", package = "cobalt")

## Using CBPS() for generating covariate balancing propensity score weights
cbps.out <- CBPS(treat ~ age + educ + black + hispan + married + 
             nodegree + re74 + re75, data = lalonde, 
             standardize = FALSE)
bal.tab(cbps.out)

}
\keyword{design}
