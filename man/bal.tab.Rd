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

\strong{With Binary Treatments}

Prior to computation, all variables are checked for variable type, which allows users to differentiate balance statistic calculations based on type using the arguments to \code{continuous} and \code{binary}. First, if a given covariate is numeric and has only 2 levels, it is converted into a binary (0,1) variable. If 0 is a value in the original variable, it retains its value and the other value is converted to 1; otherwise, the lower value is converted to 0 and the other to 1. Next, if the covariate is not numeric or logical (i.e., is a character or factor variable), it will be split into new binary variables, named with the original variable and the value, seperated by an underscore. Otherwise, the covariate will be used as is and treated as a continuous variable.

The default balance statistic for mean differences for continuous variables is the standardized difference, which is the difference in the means divided by a measure of spread (i.e., a d-type effect size measure). This is the default because it puts the mean differences on the same scale for comparison with each other and with a given threshold. It can be helpful to see the raw mean differences if the analyst is familiar with the true significance of the mean difference irrespective of the spread of the variable in the sample. For binary variables, the default balance statistic is the raw difference in proportion. Although standardized differences can be computed, proportion differences for binary variables are already on the same scale, and computing the standardized difference can obscure the true difference in proprotion by dividing the difference in proportion by a number that is itself a function of the observed proportions. For example, if \eqn{XT = .2} and \eqn{XC = .3}, the standardized difference in proportion would be different from that if \eqn{XT = .5} and \eqn{XC = .6}, which seems counterintuitive and not a useful distinction. However, to remain in line with the methodological literature and for comparability with other balance assessment tools (e.g., MatchIt's \code{summary()}), the option to use standardized differences for binary variables remains.

Standardized differences are calculated as follows: the numerator is the mean of the treated group minus the mean of the control group, and the denominator is a measure of spread calculated in accordance with the argument to \code{s.d.denom} or the default of the specific method used. Common approaches in the literature include using the standard deviation of the treated group or using the "pooled" standard deviation (i.e., the square root of the mean of the group variances) in calculating standardized mean differences. The computed spread \code{bal.tab()} uses is always that of the full, unadjusted sample (i.e., before matching, weighting, or subclassification), as recommended by Stuart (2010) and institued in \pkg{MatchIt}, though some analysts and other packages use the spread of the sample in question (i.e., before or after adjustment). One reason to favor the use of the spead of the unadjusted sample is that it prevents the paradoxical situation of adjustment descreasing both the mean difference and the spread of the sample, yielding a larger standardized mean difference than that prior to adjusting even though the adjusted groups are now more similar.

Although the reported mean differences (standardized or otherwise) are presented as they are calculated, when \code{m.threshold} is specified, it is compared to the absolute mean difference for simplicty of reporting. \code{bal.tab()} takes the absolute value of the input to \code{m.threshold} as well before comparing.

Variance ratios are computed within-sample, with the larger of the two variances in the numerator, always yielding values greater than or equal to 1. Variance ratios are not calculated for binary variables since they are only a function of the group proportions and thus provide the same information as differences in proportion. \code{bal.tab()} takes the reciprocal of the input to \code{v.threshold} if it is less than 1.

When matching is used for conditioning, the presented sample sizes are calculated simply by summing the number of observations in each group with matching weights greater than 0. When weighting is used, an "effective sample size" is calculated for each group using the following formula: \eqn{(\sum w)^2 / \sum w^2}, as is used in \pkg{twang}. The effective sample size is "approximately the number of observations from a simple random sample that yields an estimate with sampling variation equal to the sampling variation obtained with the weighted comparison observations" (Ridgeway et al., 2016). The calculated number tends to underestimate the true effective sample size of the weighted samples.

When subclassification is used, the balance tables for each subclass stored in \code{$Subclass.Balance} use values calculated as described above. For the aggregate balance table stored in \code{$Balance.Across.Subclass},  the values of each statistic are computed as a weighted average of the statistic across subclasses, weighted by the proportion of units in each subclass, as is done in \pkg{MatchIt}. 

\strong{With Continuous Treatments}

When continuous treatment variables are considered, the balance statistic calculated is the Pearson correlation between the covariate and treatment. Currently only weighting is supported with continuous treatments; the calculated correlation after adjustment is weighted by the generated weights.

\strong{Quick}

Calculations can take some time, especially when there are many variables, interactions, or clusters. One reason for this is that \code{bal.tab()} computes all values that it can, even if they are not requested by the user. For example, even if \code{un = FALSE}, which is the default, values for the unadjusted samples are still calculated. This can be useful if the output is to be further examined with \code{print()} or \code{love.plot()} or is to be used in some other way after the original call to \code{bal.tab()}. To avoid these extra calculations if they are not needed, users can set \code{quick = TRUE}, which will often quite dramatically speed up calculation and still display all statistics the user requests. For simple and quick model comparisons, it may be quite useful to do so, but for model reporting and graphical displays, it may be more useful to leave \code{quick = FALSE}, which is the default.
}
\value{
An object of class \code{"bal.tab"}. The use of cotinuous treatments, subclasses, and/or clusters will also cause the object to inherit other classes. The class \code{"bal.tab"} has its own \code{\link{print}} method (\code{print.bal.tab}), which formats the output nicely and in accordance with print-related options given in the call to \code{bal.tab()}. Each inherited class also has its own \code{print} method.
}
\references{
Ridgeway, G., McCaffrey, D., Morral, A., Burgette, L., & Griffin, B. A. (2016). Toolkit for Weighting and Analysis of Nonequivalent Groups: A tutorial for the twang package. R vignette. RAND.

Stuart, E. A. (2010). Matching Methods for Causal Inference: A Review and a Look Forward. Statistical Science, 25(1), 1-21. http://doi.org/10.1214/09-STS313
}
\author{
Noah Greifer \email{noah@unc.edu}
}

\seealso{

\code{\link{bal.tab.matchit}} for the method for objects returned by \pkg{MatchIt}.

\code{\link{bal.tab.ps}} for the method for objects returned by \pkg{twang}.

\code{\link{bal.tab.Match}} for the method for objects returned by \pkg{Matching}.

\code{\link{bal.tab.CBPS}} for the method for objects returned by \pkg{CBPS}.

\code{\link{bal.tab.formula}} and \code{\link{bal.tab.data.frame}} for the methods for \code{formula} and data frame interfaces when the user has covariate values and weights (including matching weights) or subclasses or wants to evaluate balance on an unconditioned data set.

\code{\link{bal.tab.ebalance}} for the method for objects returned by \pkg{ebal}.

}
\keyword{design}




