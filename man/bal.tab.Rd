\name{bal.tab}
\alias{bal.tab}

\title{
Generates Balance Statistics
}
\description{
Generates balance statistics on covariates in relation to an observed treatment variable. It is a generic function that dispatches to the method corresponding to the class of the first argument. This page links to each method page and documents the calculation and details of aspects all \code{bal.tab} methods use.
}
\usage{
bal.tab(\dots)
}

\arguments{
  \item{\dots}{
arguments passed to other methods. These arguments may be data-related, computation-related, or print-related.
}
}
\details{
\code{bal.tab()} performs various calculations on the the data objects given, and some of these calculations are not transparent on the help pages of the individual methods. This page details the calculations that are used across \code{bal.tab} methods.

\strong{With Binary Point Treatments}

Prior to computation, all variables are checked for variable type, which allows users to differentiate balance statistic calculations based on type using the arguments to \code{continuous} and \code{binary}. First, if a given covariate is numeric and has only 2 levels, it is converted into a binary (0,1) variable. If 0 is a value in the original variable, it retains its value and the other value is converted to 1; otherwise, the lower value is converted to 0 and the other to 1. Next, if the covariate is not numeric or logical (i.e., is a character or factor variable), it will be split into new binary variables, named with the original variable and the value, seperated by an underscore. Otherwise, the covariate will be used as is and treated as a continuous variable.

The default balance statistic for mean differences for continuous variables is the standardized difference, which is the difference in the means divided by a measure of spread (i.e., a d-type effect size measure). This is the default because it puts the mean differences on the same scale for comparison with each other and with a given threshold. It can be helpful to see the raw mean differences if the analyst is familiar with the true significance of the mean difference irrespective of the spread of the variable in the sample. For binary variables, the default balance statistic is the raw difference in proportion. Although standardized differences can be computed, proportion differences for binary variables are already on the same scale, and computing the standardized difference can obscure the true difference in proprotion by dividing the difference in proportion by a number that is itself a function of the observed proportions. For example, if \eqn{XT = .2} and \eqn{XC = .3}, the standardized difference in proportion would be different from that if \eqn{XT = .5} and \eqn{XC = .6}, which seems counterintuitive and not a useful distinction. However, to remain in line with the methodological literature and for comparability with other balance assessment tools (e.g., MatchIt's \code{summary()}), the option to use standardized differences for binary variables remains.

Standardized differences are calculated as follows: the numerator is the mean of the treated group minus the mean of the control group, and the denominator is a measure of spread calculated in accordance with the argument to \code{s.d.denom} or the default of the specific method used. Common approaches in the literature include using the standard deviation of the treated group or using the "pooled" standard deviation (i.e., the square root of the mean of the group variances) in calculating standardized mean differences. The computed spread \code{bal.tab()} uses is always that of the full, unadjusted sample (i.e., before matching, weighting, or subclassification), as recommended by Stuart (2010) and institued in \pkg{MatchIt}, though some analysts and other packages use the spread of the sample in question (i.e., before or after adjustment). One reason to favor the use of the spead of the unadjusted sample is that it prevents the paradoxical situation of adjustment descreasing both the mean difference and the spread of the sample, yielding a larger standardized mean difference than that prior to adjusting even though the adjusted groups are now more similar.

Although the reported mean differences (standardized or otherwise) are presented as they are calculated, when \code{m.threshold} is specified, it is compared to the absolute mean difference for simplicty of reporting. \code{bal.tab()} takes the absolute value of the input to \code{m.threshold} as well before comparing.

Variance ratios are computed within-sample, with the larger of the two variances in the numerator, always yielding values greater than or equal to 1. Variance ratios are not calculated for binary variables since they are only a function of the group proportions and thus provide the same information as differences in proportion. \code{bal.tab()} takes the reciprocal of the input to \code{v.threshold} if it is less than 1.

When matching is used for conditioning, the presented sample sizes are calculated simply by summing the number of observations in each group with matching weights greater than 0. When weighting is used, an "effective sample size" is calculated for each group using the following formula: \eqn{(\sum w)^2 / \sum w^2}, as is used in \pkg{twang}. The effective sample size is "approximately the number of observations from a simple random sample that yields an estimate with sampling variation equal to the sampling variation obtained with the weighted comparison observations" (Ridgeway et al., 2016). The calculated number tends to underestimate the true effective sample size of the weighted samples. The number depends on the variability of the weights, so sometimes trimming units with large weights can actually increase the effective sample size, even though units are being discarded.

When subclassification is used, the balance tables for each subclass stored in \code{$Subclass.Balance} use values calculated as described above. For the aggregate balance table stored in \code{$Balance.Across.Subclass},  the values of each statistic are computed as a weighted average of the statistic across subclasses, weighted by the proportion of units in each subclass. See \code{\link{bal.tab.subclass}} for more details.

\strong{With Continuous Point Treatments}

When continuous treatment variables are considered, the balance statistic calculated is the Pearson correlation between the covariate and treatment. The calculated correlation after adjustment is weighted by the generated weights.

\strong{With Multinomial Point Treatments}

For information on using \code{bal.tab()} with multiple categorical treatments, see \code{\link{bal.tab.multi}}. Essentially, \code{bal.tab()} compares pairs of treatment groups in a standard way.

\strong{With Longitudinal Treatments}

For information on using \code{bal.tab()} with longitudinal treatments, see \code{\link{bal.tab.msm}}. Essentially, \code{bal.tab()} summarizes balance at each time point and summarizes across time points.

\strong{With Clustered or Multiply Imputed Data}

For information on using \code{bal.tab()} with clustered data, see \code{\link{bal.tab.cluster}}. For information on using \code{bal.tab()} with multiply imputed data, see \code{\link{bal.tab.imp}}. 

\strong{Quick}

Calculations can take some time, especially when there are many variables, interactions, or clusters. One reason for this is that \code{bal.tab()} computes all values that it can, even if they are not requested by the user. For example, even if \code{un = FALSE}, which is the default, values for the unadjusted samples are still calculated. This can be useful if the output is to be further examined with \code{print()} or \code{love.plot()} or is to be used in some other way after the original call to \code{bal.tab()}. To avoid these extra calculations if they are not needed, users can set \code{quick = TRUE}, which will often quite dramatically speed up calculation and still display all statistics the user requests. For simple and quick model comparisons, it may be quite useful to do so, but for model reporting and graphical displays, it may be more useful to leave \code{quick = FALSE}, which is the default.
}
\value{
An object of class \code{"bal.tab"}. The use of continuous treatments, subclasses, clusters, and/or imputations will also cause the object to inherit other classes. The class \code{"bal.tab"} has its own \code{\link{print}} method (\code{\link{print.bal.tab}}), which formats the output nicely and in accordance with print-related options given in the call to \code{bal.tab()}, and which can be called with its own options. Each inherited class also has its own \code{print} method.

For scenarios with binary point treatments and no subclasses, imputations, or clusters, the following are the elements of the \code{bal.tab} object:

\item{Balance}{A data frame containing balance information for each covariate.  Balance contains the following columns:
\itemize{
\item{\code{Type}: Whether the covariate is binary, continuous, or a measure of distance (e.g., the propensity score).}
\item{\code{M.C.Un}: The mean of the control group prior to adjusting.}
\item{\code{M.T.Un}: The mean of the treated group prior to adjusting.}
\item{\code{Diff.Un}: The (standardized) difference in means between the two groups prior to adjusting.}
\item{\code{V.Ratio.Un}: The ratio of the variances of the two groups prior to adjusting.  \code{NA} for binary variables.  If less than 1, the reciprocal is reported.}
\item{\code{KS.Un}: The KS statistic of the two groups prior to adjusting.  \code{NA} for binary variables.}
\item{\code{M.C.Adj}: The mean of the control group after adjusting.}
\item{\code{M.T.Adj}: The mean of the treated group after adjusting.}
\item{\code{Diff.Adj}: The (standardized) difference in means between the two groups after adjusting.}
\item{\code{M.Threshold}: Whether or not the calculated mean difference after adjusting exceeds or is within the threshold given by \code{m.threshold}.  If \code{m.threshold} is \code{NULL}, this column will be \code{NA}.}
\item{\code{V.Ratio.Adj}: The ratio of the variances of the two groups after adjusting.  \code{NA} for binary variables.  If less than 1, the reciprocal is reported.}
\item{\code{V.Threshold}: Whether or not the calculated variance ratio after adjusting exceeds or is within the threshold given by \code{v.threshold} for continuous variables.  If \code{v.threshold} is \code{NULL}, this column will be \code{NA}.}
\item{\code{KS.Adj}: The KS statistic of the two groups after adjusting.  \code{NA} for binary variables.}
\item{\code{KS.Threshold}: Whether or not the calculated KS statistic after adjusting exceeds or is within the threshold given by \code{ks.threshold} for continuous variables.  If \code{ks.threshold} is \code{NULL}, this column will be \code{NA}.}
}}
\item{Balanced.Means}{If \code{m.threshold} is specified, a table tallying the number of variables that exceed or are within the threshold for mean differences.}
\item{Max.Imbalance.Means}{If \code{m.threshold} is specified, a table displaying the variable with the greatest absolute mean difference.}
\item{Balanced.Variances}{If \code{v.threshold} is specified, a table tallying the number of variables that exceed or are within the threshold for variance ratios.}
\item{Max.Imbalance.Variance}{If \code{v.threshold} is specified, a table displaying the variable with the greatest variance ratio.}
\item{Balanced.KS}{If \code{ks.threshold} is specified, a table tallying the number of variables that exceed or are within the threshold for KS statistcs.}
\item{Max.Imbalance.KS}{If \code{ks.threshold} is specified, a table displaying the variable with the greatest KS statistic.}
\item{Observations}{A table displaying the sample sizes before and after adjusting.}
\item{call}{The original function call, if adjustment was performed by a function in another package.}
\item{print.options}{A list of print options passed to \code{print.bal.tab}.}

If the treatment is continuous, instead of mean differences, variance ratios, and KS statistics, \code{bal.tab()} will produce correlations between the covariates and the treatment. The corresponding entries in the output will be "Corr.Un", "Corr.Adj", and "R.threshold" (and accordingly for the balance tally and maximum imbalance tables).

If multiple weights are supplied, \code{"Adj"} in \code{Balance} will be replaced by the provided names of the sets of weights, and extra columns will be added for each set of weights. Additional columns and rows for other items in the output will be created as well.

For \code{bal.tab} output with subclassification, see \code{\link{bal.tab.subclass}}.
}

\references{
Ridgeway, G., McCaffrey, D., Morral, A., Burgette, L., & Griffin, B. A. (2016). Toolkit for Weighting and Analysis of Nonequivalent Groups: A tutorial for the twang package. R vignette. RAND.

Stuart, E. A. (2010). Matching Methods for Causal Inference: A Review and a Look Forward. Statistical Science, 25(1), 1-21. http://doi.org/10.1214/09-STS313
}
\author{
Noah Greifer 
}

\seealso{

\code{\link{bal.tab.matchit}} for the method for objects returned by \pkg{MatchIt}.

\code{\link{bal.tab.weightit}} for the method for objects returned by \pkg{WeightIt}.

\code{\link{bal.tab.ps}} for the method for \code{ps} objects returned by \pkg{twang}.

\code{\link{bal.tab.mnps}} for the method for \code{mnps} objects returned by \pkg{twang}.

\code{\link{bal.tab.Match}} for the method for objects returned by \pkg{Matching}.

\code{\link{bal.tab.optmatch}} for the method for objects returned by \pkg{optmatch}.

\code{\link{bal.tab.CBPS}} for the method for objects returned by \pkg{CBPS}.

\code{\link{bal.tab.formula}} and \code{\link{bal.tab.data.frame}} for the methods for \code{formula} and data frame interfaces when the user has covariate values and weights (including matching weights) or subclasses or wants to evaluate balance on an unconditioned data set.

\code{\link{bal.tab.ebalance}} for the method for objects returned by \pkg{ebal}.


}