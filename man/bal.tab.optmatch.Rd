\name{bal.tab.optmatch}
\alias{bal.tab.optmatch}

\title{
Balance Statistics for optmatch Objects
}
\description{
Generates balance statistics for \code{optmatch} objects from \pkg{optmatch}.
}
\usage{
\method{bal.tab}{optmatch}(optmatch, formula = NULL, data = NULL, treat = NULL, 
    covs = NULL, int = FALSE, distance = NULL, addl = NULL, 
    continuous = c("std", "raw"), binary = c("raw", "std"), 
    s.d.denom = c("treated", "control", "pooled"), 
    m.threshold = NULL, v.threshold = NULL, 
    ks.threshold = NULL, un = FALSE, disp.means = FALSE, 
    disp.v.ratio = FALSE, disp.ks = FALSE, cluster = NULL, 
    which.cluster = NULL, cluster.summary = TRUE, 
    quick = FALSE, ...)
}

\arguments{
  \item{optmatch}{
an \code{optmatch} object; the output of a call to \code{pairmatch()} or \code{fullmatch()} from the \pkg{optmatch} package. This should be a factor vector containing the matching stratum membership for each unit.
}
  \item{formula}{
a \code{formula} with the treatment variable as the response and the covariates for which balanace is to be assessed as the predictors. All named variables must be in \code{data}. See Details.
}
  \item{data}{
a data frame containing all the variables named in \code{formula}. See Details.
}
  \item{treat}{
a vector of treatment statuses. See Details.
}
  \item{covs}{
a data frame of covariate values for which to check balance. See Details.
}
  \item{int}{
\code{logical}; whether or not to include squares and 2-way interactions of covariates included in \code{formula} or \code{covs} and in \code{addl}.
}
  \item{distance}{
Optional; either a vector containing distance values (e.g., propensity scores) for each unit or a string containing the name of the distance variable in \code{data}.
}
  \item{addl}{
a data frame of additional covariates for which to present balance. These may be covariates included in the original dataset but not included in \code{formula} or \code{covs}. In general, it makes more sense to include all desired variables in \code{formula} or \code{covs} than in \code{addl}.
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
  \item{ks.threshold}{
a numeric value for the threshold for Kolmogorov-Smirnov statistics. Must be between 0 and 1. 
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
  \item{disp.ks}{
\code{logical}; whether to display Kolmogorov-Smirnov statistics in balance output.
}
  \item{cluster}{
either a vector containing cluster membserhip for each unit or a string containing the name of the cluster membership variable in \code{data} or the CBPS object. See \code{\link{bal.tab.cluster}} for details.
}
  \item{which.cluster}{
which cluster(s) to display if \code{cluster} is specified. See \code{\link{bal.tab.cluster}} for details.
}
  \item{cluster.summary}{
\code{logical}; whether to display the cluster summary table if \code{cluster} is specified. See \code{\link{bal.tab.cluster}} for details.
}
  \item{quick}{
\code{logical}; if \code{TRUE}, will not compute any values that will not be displayed. Leave \code{FALSE} if computed values not displayed will be used later.
}
  \item{...}{
further arguments passed to or from other methods. They are ignored in this function.
}
}
\details{
\code{bal.tab.optmatch()} generates a list of balance summaries for the optmatch object given.

The input to \code{bal.tab.optmatch()} must include either both \code{formula} and \code{data} or both \code{treat} and \code{covs}. 

All balance statistics are calculated whether they are displayed by print or not, unless \code{quick = TRUE}. The threshold values (\code{m.threshold} and \code{v.threshold}) control whether extra columns should be inserted into the Balance table describing whether the balance statistics in question exceeded or were within the threshold. Including these thresholds also creates summary tables tallying the number of variables that exceeded and were within the threshold and displaying the variables with the greatest imbalance on that balance measure.

The inputs (if any) to \code{covs} must be a data frame; if more than one variable is included, this is straightforward (i.e., because \code{data[,c("v1", "v2")]} is already a data frame), but if only one variable is used (e.g., \code{data[,"v1"]}), R will coerce it to a vector, thus making it unfit for input. To avoid this, simply wrap the input to \code{covs} in \code{data.frame()} or use \code{subset()} if only one variable is to be added. Again, when more than one variable is included, the input is general already a data frame and nothing needs to be done.
}
\value{
If clusters are not specified, an object of class \code{"bal.tab"} containing balance summaries for the \code{optmatch} object. The following are the elements of \code{bal.tab}:
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
\item{call}{\code{NULL}.}
\item{print.options}{A list of print options passed to \code{print.bal.tab}.}

If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. Each balance summary is a balance table as described in \code{Balance} above. The summary of balance across clusters displays the mean, median, and maximum mean difference and variance ratio after adjustment for each covariate across clusters. Minimum statistics are calculated as well, but not displayed. To see these, use the options in \code{\link{print.bal.tab.cluster}}.
}

\author{
Noah Greifer \email{noah@unc.edu}
}

\seealso{
\code{\link{bal.tab}} for details of calculations.
}
\examples{
library("optmatch"); data("lalonde", package = "cobalt")

lalonde$prop.score <- glm(treat ~ age + educ + race + 
            married + nodegree + re74 + re75, 
            data = lalonde, family = binomial)$fitted.values
pm <- pairmatch(treat ~ prop.score, data = lalonde)

## Using formula and data
bal.tab(pm, treat ~ age + educ + race + 
        married + nodegree + re74 + re75, data = lalonde)

## Using treat and covs
covariates <- subset(lalonde, select=-c(treat, re78, prop.score))
bal.tab(pm, treat = lalonde$treat, covs = covariates)

}
\keyword{design}

